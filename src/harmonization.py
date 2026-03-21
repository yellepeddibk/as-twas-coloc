"""
harmonization.py
----------------
GWAS summary-statistics harmonization for the AS TWAS pipeline.

Key responsibilities
--------------------
1. Rename raw columns to the canonical schema (via schemas.py aliases).
2. Validate required columns are present.
3. Optionally lift hg19 -> hg38 (placeholder wrapper).
4. Construct GTEx-style varIDs:  chr_pos_ref_alt_b38.
5. Harmonize alleles (strand flip / allele swap detection).
6. Drop ambiguous palindromic SNPs (AT / CG) when AF disambiguation is not available.
7. Log every harmonization decision clearly.

Pitfalls documented in comments:
- hg19 vs hg38 mismatch -> varID mismatch, 0% model SNPs used
- GTEx uses varID not rsID for model look-up
- Strand flips can masquerade as novel SNPs
- Palindromic SNPs (AT/CG) cannot be reliably strand-flipped without AF
"""

from __future__ import annotations

import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict

import pandas as pd

from src.logging_utils import get_logger
from src.schemas import (
    HARMONIZED_GWAS_REQUIRED,
    apply_column_aliases,
    validate_harmonized_schema,
)

logger = get_logger(__name__)

# Complement bases for strand-flip detection
_COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}
# Ambiguous palindromic pairs (both orderings)
_PALINDROMIC_PAIRS = {frozenset({"A", "T"}), frozenset({"C", "G"})}


# ─────────────────────────────────────────────────────────────
# Public entry point
# ─────────────────────────────────────────────────────────────


def harmonize_gwas(
    df: pd.DataFrame,
    cfg: Dict[str, Any],
    drop_palindromic: bool = True,
    strict_real_mode: bool = False,
) -> pd.DataFrame:
    """Harmonize a raw GWAS DataFrame to the canonical schema.

    Steps performed:
    1. Apply column aliases.
    2. Coerce data types.
    3. Optionally lift coordinates (placeholder).
    4. Construct varID.
    5. Flip alleles where needed.
    6. Drop palindromic SNPs (unless AF disambiguation is available).
    7. Validate final schema.
    8. Log summary statistics.

    Parameters
    ----------
    df:
        Raw GWAS DataFrame (will not be modified in-place).
    cfg:
        Full pipeline configuration dict.
    drop_palindromic:
        Override for palindromic-SNP behaviour (default: use config value).

    Returns
    -------
    pd.DataFrame
        Harmonized GWAS table conforming to the canonical schema.
    """
    df = df.copy()
    gwas_cfg = cfg.get("gwas", {})
    qc: Dict[str, int] = {
        "input_snp_count": len(df),
        "lifted_snp_count": len(df),
        "failed_liftover_count": 0,
        "palindromic_dropped_count": 0,
        "final_retained_snp_count": 0,
    }

    # ── 1. Rename columns ─────────────────────────────────────
    df = apply_column_aliases(df)

    # ── 2. Coerce types ───────────────────────────────────────
    df = _coerce_types(df)
    df = _ensure_beta_from_or(df)
    df = _ensure_total_sample_size(df, gwas_cfg)

    # ── 3. Liftover hg19 -> hg38 ──────────────────────────────
    source_build = gwas_cfg.get("source_build", "hg38")
    if source_build.lower() in ("hg19", "grch37", "b37"):
        logger.warning(
            "Source build is %s.  Liftover to hg38 is required for GTEx varID matching. "
            "Attempting hg19 -> hg38 liftover.",
            source_build,
        )
        df = _liftover_hg19_to_hg38(df, cfg, strict_real_mode=strict_real_mode)
        qc["lifted_snp_count"] = int(df.attrs.get("liftover_mapped_snps", len(df)))
        qc["failed_liftover_count"] = int(df.attrs.get("liftover_failed_snps", 0))
    else:
        logger.info("Source build is %s; no liftover needed.", source_build)

    # ── 4. Construct GTEx-style varID ─────────────────────────
    df = _build_varid(df)

    # ── 5. Strand-flip / allele-swap check ───────────────────
    df = _harmonize_alleles(df)

    # ── 6. Drop palindromic SNPs ──────────────────────────────
    use_drop = gwas_cfg.get("drop_palindromic_snps", True) if drop_palindromic else drop_palindromic
    af_available = gwas_cfg.get("af_disambiguation_available", False)
    if use_drop and not af_available:
        before_drop = len(df)
        df = _drop_palindromic(df)
        qc["palindromic_dropped_count"] = before_drop - len(df)
    elif not use_drop:
        logger.info("Palindromic SNP removal skipped (drop_palindromic=False).")
    else:
        logger.info(
            "AF-based disambiguation is available; retaining palindromic SNPs "
            "(AF-based logic not yet implemented — treating as retained)."
        )

    # ── 7. Remove rows with missing required fields ───────────
    before = len(df)
    required_non_null = ["chrom", "pos", "effect_allele", "non_effect_allele", "beta", "se", "pvalue"]
    df = df.dropna(subset=[c for c in required_non_null if c in df.columns])
    dropped = before - len(df)
    if dropped:
        logger.info("Dropped %d rows with missing values in required columns.", dropped)

    # ── 8. Validate schema ────────────────────────────────────
    validate_harmonized_schema(df)

    # ── 9. Summary ────────────────────────────────────────────
    logger.info(
        "Harmonization complete: %d SNPs retained (source build: %s).",
        len(df),
        source_build,
    )

    qc["final_retained_snp_count"] = len(df)
    df.attrs["harmonization_qc"] = qc
    return df


# ─────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────


def _coerce_types(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce numeric columns and normalise allele strings."""
    for col in ("pos", "N", "N_cases", "N_controls"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in ("beta", "se", "pvalue", "eaf", "maf", "zscore"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Normalise chromosome: strip "chr" prefix, keep as string
    if "chrom" in df.columns:
        df["chrom"] = (
            df["chrom"]
            .astype(str)
            .str.replace(r"^chr", "", regex=True)
            .str.strip()
        )

    # Normalise alleles to uppercase
    for col in ("effect_allele", "non_effect_allele"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.upper().str.strip()

    return df


def _ensure_beta_from_or(df: pd.DataFrame) -> pd.DataFrame:
    """Populate beta from log(OR) when beta is missing and OR is present."""
    if "beta" in df.columns:
        return df
    if "or" not in df.columns:
        return df

    df = df.copy()
    or_vals = pd.to_numeric(df["or"], errors="coerce")
    valid_or = or_vals > 0
    df["beta"] = pd.NA
    df.loc[valid_or, "beta"] = or_vals[valid_or].map(float).map(math.log)
    logger.info(
        "Derived beta from OR for %d rows (invalid OR rows: %d).",
        int(valid_or.sum()),
        int((~valid_or).sum()),
    )
    return df


def _ensure_total_sample_size(df: pd.DataFrame, gwas_cfg: Dict[str, Any]) -> pd.DataFrame:
    """Ensure N exists; fill from config when per-variant N is absent."""
    if "N" in df.columns:
        return df

    n_total = gwas_cfg.get("n_total")
    if n_total is None:
        raise ValueError("Missing N column and gwas.n_total is not set in config.")

    df = df.copy()
    df["N"] = int(n_total)
    logger.info("Filled missing N column from config.gwas.n_total=%s.", n_total)
    return df


def _liftover_hg19_to_hg38(
    df: pd.DataFrame,
    cfg: Dict[str, Any],
    strict_real_mode: bool = False,
) -> pd.DataFrame:
    """Lift coordinates from hg19 to hg38 using CrossMap.

    If CrossMap is not configured and strict_real_mode is False, this falls back
    to passthrough behavior for mock/demo compatibility.
    """
    liftover_cfg = cfg.get("liftover", {})
    enabled = bool(liftover_cfg.get("enabled", False))
    tool = str(liftover_cfg.get("tool", "crossmap")).lower()
    chain_file = liftover_cfg.get("chain_file")
    crossmap_bin = liftover_cfg.get("crossmap_bin", "CrossMap.py")

    if not enabled:
        if strict_real_mode:
            raise RuntimeError(
                "Liftover is required in real mode but liftover.enabled is false in config."
            )
        return _liftover_passthrough(df, "hg19", "hg38")

    if tool not in {"crossmap", "pyliftover"}:
        raise ValueError(f"Unsupported liftover tool '{tool}'. Use 'crossmap' or 'pyliftover'.")

    if not chain_file:
        if strict_real_mode:
            raise RuntimeError("liftover.chain_file is required for CrossMap liftover.")
        return _liftover_passthrough(df, "hg19", "hg38")

    chain_path = Path(chain_file)
    if not chain_path.is_absolute():
        chain_path = Path(cfg.get("_base_dir", ".")) / chain_path
    chain_path = chain_path.resolve()
    if not chain_path.exists():
        if strict_real_mode:
            raise FileNotFoundError(f"Liftover chain file not found: {chain_path}")
        return _liftover_passthrough(df, "hg19", "hg38")

    required = ["chrom", "pos"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Cannot run liftover: missing columns {missing}")

    if tool == "pyliftover":
        return _liftover_with_pyliftover(df, chain_path)

    work_df = df.copy().reset_index(drop=True)
    work_df["_row_id"] = work_df.index.astype(str)

    with tempfile.TemporaryDirectory(prefix="as_liftover_") as tmpdir:
        tmp = Path(tmpdir)
        input_bed = tmp / "input.bed"
        output_bed = tmp / "output.bed"

        bed_df = pd.DataFrame(
            {
                "chrom": "chr" + work_df["chrom"].astype(str).str.replace(r"^chr", "", regex=True),
                "start": (pd.to_numeric(work_df["pos"], errors="coerce") - 1).fillna(-1).astype(int),
                "end": pd.to_numeric(work_df["pos"], errors="coerce").fillna(-1).astype(int),
                "row_id": work_df["_row_id"],
            }
        )
        bed_df = bed_df.loc[(bed_df["start"] >= 0) & (bed_df["end"] > 0)]
        bed_df.to_csv(input_bed, sep="\t", header=False, index=False)

        cmd = [str(crossmap_bin), "bed", str(chain_path), str(input_bed), str(output_bed)]
        logger.info("Running CrossMap liftover: %s", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError as exc:
            fallback_cmd = [sys.executable, "-m", "CrossMap", "bed", str(chain_path), str(input_bed), str(output_bed)]
            logger.warning(
                "CrossMap executable '%s' not found. Retrying with module invocation: %s",
                crossmap_bin,
                " ".join(fallback_cmd),
            )
            try:
                subprocess.run(fallback_cmd, check=True, capture_output=True, text=True)
            except FileNotFoundError as inner_exc:
                if strict_real_mode:
                    raise RuntimeError("CrossMap executable not found. Install CrossMap and configure crossmap_bin.") from inner_exc
                return _liftover_passthrough(df, "hg19", "hg38")
            except subprocess.CalledProcessError as inner_exc:
                inner_msg = (inner_exc.stderr or inner_exc.stdout or "").strip()
                raise RuntimeError(f"CrossMap liftover failed: {inner_msg}") from inner_exc
        except subprocess.CalledProcessError as exc:
            msg = (exc.stderr or exc.stdout or "").strip()
            raise RuntimeError(f"CrossMap liftover failed: {msg}") from exc

        mapped_cols = ["mapped_chrom", "mapped_start", "mapped_end", "row_id"]
        if not output_bed.exists():
            raise RuntimeError("CrossMap completed but did not produce output BED.")

        mapped_df = pd.read_csv(output_bed, sep="\t", header=None, usecols=[0, 1, 2, 3], names=mapped_cols)
        mapped_df["row_id"] = mapped_df["row_id"].astype(str)
        mapped_df["mapped_pos"] = pd.to_numeric(mapped_df["mapped_end"], errors="coerce")
        mapped_df["mapped_chrom"] = mapped_df["mapped_chrom"].astype(str).str.replace(r"^chr", "", regex=True)

    merged = work_df.merge(mapped_df[["row_id", "mapped_chrom", "mapped_pos"]], left_on="_row_id", right_on="row_id", how="left")

    input_count = len(work_df)
    mapped_count = int(merged["mapped_pos"].notna().sum())
    failed_count = input_count - mapped_count

    lifted = merged.loc[merged["mapped_pos"].notna()].copy()
    lifted["chrom"] = lifted["mapped_chrom"]
    lifted["pos"] = lifted["mapped_pos"].astype(int)
    lifted = lifted.drop(columns=["_row_id", "row_id", "mapped_chrom", "mapped_pos"])

    logger.info(
        "Liftover complete: input=%d mapped=%d failed=%d.",
        input_count,
        mapped_count,
        failed_count,
    )

    lifted.attrs["liftover_performed"] = True
    lifted.attrs["liftover_input_snps"] = input_count
    lifted.attrs["liftover_mapped_snps"] = mapped_count
    lifted.attrs["liftover_failed_snps"] = failed_count
    return lifted


def _liftover_with_pyliftover(df: pd.DataFrame, chain_path: Path) -> pd.DataFrame:
    """Lift coordinates with pyliftover using a UCSC chain file."""
    try:
        from pyliftover import LiftOver
    except ImportError as exc:
        raise RuntimeError("pyliftover is not installed. Install it with: pip install pyliftover") from exc

    lo = LiftOver(str(chain_path))
    work_df = df.copy().reset_index(drop=True)

    mapped_chrom = []
    mapped_pos = []
    mapped_ok = []

    for _, row in work_df.iterrows():
        chrom = f"chr{str(row['chrom']).replace('chr', '')}"
        pos = int(row["pos"])
        conv = lo.convert_coordinate(chrom, pos)
        if conv:
            tgt_chr, tgt_pos = conv[0][0], conv[0][1]
            mapped_chrom.append(str(tgt_chr).replace("chr", ""))
            mapped_pos.append(int(round(tgt_pos)))
            mapped_ok.append(True)
        else:
            mapped_chrom.append(None)
            mapped_pos.append(None)
            mapped_ok.append(False)

    work_df["_mapped_ok"] = mapped_ok
    work_df["_mapped_chrom"] = mapped_chrom
    work_df["_mapped_pos"] = mapped_pos

    input_count = len(work_df)
    mapped_count = int(work_df["_mapped_ok"].sum())
    failed_count = input_count - mapped_count

    lifted = work_df.loc[work_df["_mapped_ok"]].copy()
    lifted["chrom"] = lifted["_mapped_chrom"]
    lifted["pos"] = lifted["_mapped_pos"].astype(int)
    lifted = lifted.drop(columns=["_mapped_ok", "_mapped_chrom", "_mapped_pos"])

    logger.info(
        "pyliftover complete: input=%d mapped=%d failed=%d.",
        input_count,
        mapped_count,
        failed_count,
    )

    lifted.attrs["liftover_performed"] = True
    lifted.attrs["liftover_input_snps"] = input_count
    lifted.attrs["liftover_mapped_snps"] = mapped_count
    lifted.attrs["liftover_failed_snps"] = failed_count
    return lifted


def _liftover_passthrough(df: pd.DataFrame, from_build: str, to_build: str) -> pd.DataFrame:
    """Passthrough liftover for non-strict mock/demo flows."""
    logger.warning(
        "[PLACEHOLDER] Liftover %s -> %s not performed. Coordinates are unchanged.",
        from_build,
        to_build,
    )
    passthrough = df.copy()
    passthrough.attrs["liftover_performed"] = False
    passthrough.attrs["liftover_input_snps"] = len(df)
    passthrough.attrs["liftover_mapped_snps"] = len(df)
    passthrough.attrs["liftover_failed_snps"] = 0
    return passthrough


def _build_varid(df: pd.DataFrame) -> pd.DataFrame:
    """Construct a GTEx-style varID column: chr{chrom}_{pos}_{ref}_{alt}_b38.

    The GTEx PredictDB models use varID (not rsID) as the primary key.
    A mismatch here is the most common cause of "0% of model SNPs found".

    Note: ref = non_effect_allele, alt = effect_allele by convention here.
    In practice, ref/alt assignment depends on the reference genome.
    TODO: In real usage, verify ref/alt orientation against the hg38 reference FASTA.
    """
    if "varID" in df.columns:
        logger.info("varID column already present; skipping varID construction.")
        return df

    required = ["chrom", "pos", "non_effect_allele", "effect_allele"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Cannot build varID: missing columns {missing}")

    df["varID"] = (
        "chr"
        + df["chrom"].astype(str)
        + "_"
        + df["pos"].astype(int).astype(str)
        + "_"
        + df["non_effect_allele"]
        + "_"
        + df["effect_allele"]
        + "_b38"
    )
    logger.info("Constructed varID for %d SNPs.", len(df))
    return df


def _is_palindromic(a1: str, a2: str) -> bool:
    """Return True if the allele pair is palindromic (AT or CG)."""
    return frozenset({a1.upper(), a2.upper()}) in _PALINDROMIC_PAIRS


def _strand_flip(allele: str) -> str:
    """Return the complement of a single-character allele."""
    return _COMPLEMENT.get(allele.upper(), allele)


def _harmonize_alleles(df: pd.DataFrame) -> pd.DataFrame:
    """Log palindromic counts and detect obvious multi-character indels.

    Full strand-flip correction requires a reference panel and is left as a
    TODO for real-data runs.  Here we flag the key statistics so the analyst
    is informed.

    Pitfall reminder: Strand ambiguity in palindromic SNPs means a "T" on
    one strand could be an "A" on the reference strand.  Without allele
    frequencies, you cannot know which orientation is correct.
    """
    if "effect_allele" not in df.columns or "non_effect_allele" not in df.columns:
        logger.warning("Allele columns not found; skipping allele harmonization.")
        return df

    # Flag palindromic SNPs
    pal_mask = df.apply(
        lambda r: _is_palindromic(str(r["effect_allele"]), str(r["non_effect_allele"])),
        axis=1,
    )
    n_pal = pal_mask.sum()
    logger.info(
        "Palindromic SNPs (AT/CG): %d / %d (%.1f%%).",
        n_pal,
        len(df),
        100 * n_pal / max(len(df), 1),
    )

    # Flag indels (multi-character alleles)
    indel_mask = (df["effect_allele"].str.len() > 1) | (df["non_effect_allele"].str.len() > 1)
    logger.info("Indels (multi-char alleles): %d.", indel_mask.sum())

    return df


def _drop_palindromic(df: pd.DataFrame) -> pd.DataFrame:
    """Remove palindromic (AT/CG) SNPs that cannot be strand-disambiguated.

    These SNPs are removed by default to prevent misaligned effects.
    Set ``gwas.af_disambiguation_available: true`` in config to retain them
    (you must then implement the AF-based logic).
    """
    if "effect_allele" not in df.columns or "non_effect_allele" not in df.columns:
        return df

    pal_mask = df.apply(
        lambda r: _is_palindromic(str(r["effect_allele"]), str(r["non_effect_allele"])),
        axis=1,
    )
    n_pal = pal_mask.sum()
    df = df.loc[~pal_mask].reset_index(drop=True)
    logger.info(
        "Dropped %d palindromic SNPs (AF disambiguation not available). "
        "%d SNPs remain.",
        n_pal,
        len(df),
    )
    return df


def compute_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """Add a `zscore` column as beta / se if not already present."""
    if "zscore" in df.columns:
        return df
    if "beta" in df.columns and "se" in df.columns:
        df = df.copy()
        df["zscore"] = df["beta"] / df["se"]
        logger.info("Computed zscore = beta / se.")
    return df
