"""
eqtl_regions.py
---------------
Build gene-level eQTL region files for COLOC from TWAS hits.

Output files follow the pipeline's real-COLOC contract:
  data/interim/eqtl_regions/{tissue}/{gene}.tsv.gz

Each output contains at least:
  varID, beta, se, pvalue
Optional:
  N
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, cast

import pandas as pd

from src.io_utils import ensure_dir
from src.logging_utils import get_logger

logger = get_logger(__name__)

_VARIANT_RE = re.compile(r"^(?:chr)?([^_]+)_([0-9]+)_[^_]+_[^_]+_b(?:37|38)$", re.IGNORECASE)


def normalize_gene_id(gene_id: str) -> str:
    """Normalize ENSG ID by stripping version suffix (e.g., ENSG...1.13 -> ENSG...1)."""
    return str(gene_id).split(".")[0].strip()


def normalize_varid(value: str) -> str:
    """Normalize variant IDs to no-'chr' GTEx-style form: chrom_pos_ref_alt_b38."""
    s = str(value).strip()
    if not s:
        return s
    if s.lower().startswith("chr"):
        s = s[3:]
    return s


def parse_varid_chrom_pos(varid: str) -> Tuple[Optional[str], Optional[int]]:
    """Extract chromosome and position from GTEx-style varID."""
    m = _VARIANT_RE.match(str(varid))
    if not m:
        return None, None
    chrom = m.group(1)
    try:
        pos = int(m.group(2))
    except ValueError:
        return None, None
    return chrom, pos


def load_gene_tss_from_gtf(gtf_path: str | Path) -> pd.DataFrame:
    """Load gene-level TSS coordinates from a GTF/GTF.GZ file."""
    gtf_path = Path(gtf_path)
    if not gtf_path.exists():
        raise FileNotFoundError(f"Gene annotation file not found: {gtf_path}")

    rows: List[Dict[str, object]] = []
    with _open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = parts
            if feature != "gene":
                continue

            gene_id = _extract_gtf_attr(attrs, "gene_id")
            if not gene_id:
                continue

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            chrom = str(chrom).replace("chr", "")
            tss = start_i if strand == "+" else end_i
            rows.append({"gene": normalize_gene_id(gene_id), "chrom": chrom, "tss": tss})

    if not rows:
        raise ValueError(f"No gene rows parsed from GTF: {gtf_path}")

    df = pd.DataFrame(rows).drop_duplicates(subset=["gene"]).reset_index(drop=True)
    logger.info("Loaded %d gene TSS rows from %s", len(df), gtf_path)
    return df


def infer_eqtl_file_for_tissue(eqtl_dir: str | Path, tissue: str) -> Path:
    """Infer the most likely eQTL summary file for a tissue from eqtl_dir."""
    eqtl_dir = Path(eqtl_dir)
    if not eqtl_dir.exists():
        raise FileNotFoundError(f"GTEx eQTL directory not found: {eqtl_dir}")

    candidates = sorted(
        p for p in eqtl_dir.rglob("*")
        if p.is_file() and p.suffix in {".txt", ".tsv", ".gz"}
    )

    tissue_norm = tissue.lower()
    matched = [p for p in candidates if tissue_norm in p.name.lower()]
    if not matched:
        raise FileNotFoundError(
            f"No eQTL summary file found for tissue '{tissue}' under {eqtl_dir}"
        )

    # Prefer files that look like variant-gene pairs tables.
    preferred = [p for p in matched if "variant" in p.name.lower() and "gene" in p.name.lower()]
    return preferred[0] if preferred else matched[0]


def build_eqtl_regions_for_coloc(
    twas_hits_df: pd.DataFrame,
    harmonized_gwas_df: pd.DataFrame,
    gene_tss_df: pd.DataFrame,
    eqtl_dir: str | Path,
    output_dir: str | Path,
    window_bp: int,
    chunk_size: int = 250_000,
) -> pd.DataFrame:
    """Build per-gene per-tissue eQTL region files for real COLOC."""
    required_hits = {"gene", "tissue"}
    missing_hits = required_hits - set(twas_hits_df.columns)
    if missing_hits:
        raise ValueError(f"TWAS hits missing columns: {sorted(missing_hits)}")

    required_gwas = {"varID"}
    missing_gwas = required_gwas - set(harmonized_gwas_df.columns)
    if missing_gwas:
        raise ValueError(f"Harmonized GWAS missing columns: {sorted(missing_gwas)}")

    output_dir = Path(output_dir)
    ensure_dir(output_dir)

    gwas_varids = {normalize_varid(v) for v in harmonized_gwas_df["varID"].astype(str).tolist()}

    hits = twas_hits_df[["gene", "tissue"]].copy()
    hits["gene"] = hits["gene"].map(normalize_gene_id)
    hits = hits.drop_duplicates().reset_index(drop=True)

    gene_tss = gene_tss_df.copy()
    gene_tss["gene"] = gene_tss["gene"].map(normalize_gene_id)

    summary_rows: List[Dict[str, object]] = []

    for tissue_key, tissue_hits in hits.groupby("tissue", sort=True):
        tissue = str(tissue_key)
        tissue_genes = tissue_hits["gene"].tolist()
        gene_coords = gene_tss[gene_tss["gene"].isin(tissue_genes)].copy()

        missing_coord_genes = sorted(set(tissue_genes) - set(gene_coords["gene"].tolist()))
        for gene in missing_coord_genes:
            summary_rows.append({
                "gene": gene,
                "tissue": tissue,
                "status": "missing_gene_coordinate",
                "eqtl_rows_in_window": 0,
                "shared_snps_with_gwas": 0,
                "output_path": "",
            })

        if gene_coords.empty:
            continue

        eqtl_file = infer_eqtl_file_for_tissue(eqtl_dir, tissue)
        logger.info("Using eQTL file for %s: %s", tissue, eqtl_file)

        per_gene_frames: Dict[str, List[pd.DataFrame]] = {g: [] for g in gene_coords["gene"].tolist()}

        for chunk in pd.read_csv(eqtl_file, sep="\t", low_memory=False, chunksize=chunk_size):
            col_map = _resolve_eqtl_columns(chunk.columns)
            varid_col = col_map["varid"]
            beta_col = col_map["beta"]
            se_col = col_map["se"]
            pvalue_col = col_map["pvalue"]
            optional_n = col_map.get("n")
            needed = [varid_col, beta_col, se_col, pvalue_col]
            select_cols = needed + ([optional_n] if optional_n else [])

            work = chunk[select_cols].copy()
            rename = {
                varid_col: "varID",
                beta_col: "beta",
                se_col: "se",
                pvalue_col: "pvalue",
            }
            if optional_n:
                rename[optional_n] = "N"
            work = work.rename(columns=rename)

            work["varID"] = work["varID"].astype(str).map(normalize_varid)
            work = work[work["varID"].isin(gwas_varids)]
            if work.empty:
                continue

            chrom_pos = work["varID"].map(parse_varid_chrom_pos)
            work["chrom"] = chrom_pos.map(lambda x: x[0])
            work["pos"] = chrom_pos.map(lambda x: x[1])
            work = work.dropna(subset=["chrom", "pos", "beta", "se", "pvalue"])
            if work.empty:
                continue

            for _, g in gene_coords.iterrows():
                gene = str(g["gene"])
                chrom = str(g["chrom"])
                tss = int(g["tss"])
                sub = work[
                    (work["chrom"].astype(str) == chrom)
                    & (work["pos"] >= tss - window_bp)
                    & (work["pos"] <= tss + window_bp)
                ]
                if not sub.empty:
                    keep_cols = ["varID", "beta", "se", "pvalue"] + (["N"] if "N" in sub.columns else [])
                    per_gene_frames[gene].append(sub[keep_cols])

        for _, g in gene_coords.iterrows():
            gene = str(g["gene"])
            frames = per_gene_frames.get(gene, [])
            if not frames:
                summary_rows.append({
                    "gene": gene,
                    "tissue": tissue,
                    "status": "no_eqtl_overlap",
                    "eqtl_rows_in_window": 0,
                    "shared_snps_with_gwas": 0,
                    "output_path": "",
                })
                continue

            merged = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["varID"]).reset_index(drop=True)
            out_path = output_dir / tissue / f"{gene}.tsv.gz"
            ensure_dir(out_path.parent)
            merged.to_csv(out_path, sep="\t", index=False, compression="gzip")

            summary_rows.append({
                "gene": gene,
                "tissue": tissue,
                "status": "written",
                "eqtl_rows_in_window": int(len(merged)),
                "shared_snps_with_gwas": int(merged["varID"].nunique()),
                "output_path": str(out_path),
            })

    summary_df = pd.DataFrame(summary_rows).sort_values(["status", "tissue", "gene"]).reset_index(drop=True)
    return summary_df


def _resolve_eqtl_columns(columns: Iterable[str]) -> Dict[str, Optional[str]]:
    names = {str(c).lower(): str(c) for c in columns}

    def pick(candidates: List[str], key: str, required: bool = True) -> Optional[str]:
        for c in candidates:
            if c.lower() in names:
                return names[c.lower()]
        if required:
            raise ValueError(f"Could not resolve required eQTL column '{key}' from columns: {list(columns)}")
        return None

    result: Dict[str, Optional[str]] = {
        "varid": pick(["varID", "variant_id", "variant", "snp"], "varid"),
        "beta": pick(["beta", "slope", "effect_size"], "beta"),
        "se": pick(["se", "slope_se", "beta_se", "effect_size_se"], "se"),
        "pvalue": pick(["pvalue", "pval_nominal", "pval", "p_value"], "pvalue"),
        "n": pick(["N", "n", "ma_samples", "sample_size"], "n", required=False),
    }
    # Required keys are guaranteed non-null by pick(required=True), but keep
    # Optional typing for "n" to satisfy static checkers.
    result["varid"] = cast(str, result["varid"])
    result["beta"] = cast(str, result["beta"])
    result["se"] = cast(str, result["se"])
    result["pvalue"] = cast(str, result["pvalue"])
    return result


def _open_text(path: Path):
    if path.suffix == ".gz":
        import gzip
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _extract_gtf_attr(attr_field: str, key: str) -> Optional[str]:
    m = re.search(rf'{re.escape(key)}\s+"([^"]+)"', attr_field)
    return m.group(1) if m else None
