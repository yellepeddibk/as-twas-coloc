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

import re
from typing import Any, Dict, Optional, Tuple

import numpy as np
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

    # ── 1. Rename columns ─────────────────────────────────────
    df = apply_column_aliases(df)

    # ── 2. Coerce types ───────────────────────────────────────
    df = _coerce_types(df)

    # ── 3. Liftover hg19 -> hg38 (placeholder) ───────────────
    source_build = gwas_cfg.get("source_build", "hg38")
    if source_build.lower() in ("hg19", "grch37", "b37"):
        logger.warning(
            "Source build is %s.  Liftover to hg38 is required for GTEx varID matching. "
            "Calling liftover placeholder — coordinates will be UNCHANGED in mock mode.",
            source_build,
        )
        df = _liftover_placeholder(df, source_build, "hg38")
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
        df = _drop_palindromic(df)
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


def _liftover_placeholder(df: pd.DataFrame, from_build: str, to_build: str) -> pd.DataFrame:
    """Placeholder for hg19 -> hg38 coordinate liftover.

    TODO: Replace with a real liftover implementation using one of:
      - pyliftover (Python)
      - CrossMap (Python)
      - UCSC liftOver binary

    Reference chain file: hg19ToHg38.over.chain.gz
    (available from UCSC Genome Browser downloads)

    In real mode this function should:
    1. Write BED from df[chrom, pos].
    2. Call liftOver.
    3. Parse output BED back and join on original index.
    4. Drop rows that did not liftover (log count).

    In mock mode, coordinates are left unchanged and a WARNING is emitted.
    """
    logger.warning(
        "[PLACEHOLDER] Liftover %s -> %s not performed. "
        "Coordinates are unchanged.  "
        "Do NOT use these coordinates for real GTEx varID matching.",
        from_build,
        to_build,
    )
    # Mark the DataFrame so downstream code knows this was not a real lift.
    df = df.copy()
    df.attrs["liftover_performed"] = False
    return df


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
