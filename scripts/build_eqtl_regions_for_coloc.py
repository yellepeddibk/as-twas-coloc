#!/usr/bin/env python3
"""
build_eqtl_regions_for_coloc.py
--------------------------------
Generate gene-level eQTL region files required by real COLOC execution.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add repository root to sys.path so ``src`` is importable
_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from src.config import load_config
from src.eqtl_regions import build_eqtl_regions_for_coloc, load_gene_tss_from_gtf
from src.io_utils import read_tsv_gz


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build per-gene eQTL region files for COLOC")
    p.add_argument("--config", default="config/as.yaml", help="Config path")
    p.add_argument("--base-dir", default=".", help="Repository root")
    p.add_argument(
        "--twas-hits",
        default="data/processed/twas_hits/twas_hits_fdr05.tsv.gz",
        help="TWAS hits file path (tsv.gz)",
    )
    p.add_argument(
        "--harmonized-gwas",
        default="data/interim/gwas_harmonized/as_gwas_hg38_varid.tsv.gz",
        help="Harmonized GWAS file path (tsv.gz)",
    )
    p.add_argument(
        "--gene-annotation",
        required=True,
        help="Path to GENCODE GTF/GTF.GZ annotation for gene TSS lookup",
    )
    p.add_argument(
        "--output-dir",
        default="data/interim/eqtl_regions",
        help="Output directory for per-gene eQTL region files",
    )
    p.add_argument(
        "--summary-out",
        default="data/processed/qc/eqtl_region_build_summary.tsv",
        help="Summary TSV output path",
    )
    p.add_argument(
        "--chunk-size",
        type=int,
        default=250000,
        help="Rows per chunk when scanning eQTL summary files",
    )
    return p.parse_args()


def _resolve(base_dir: Path, path_str: str) -> Path:
    p = Path(path_str)
    return p if p.is_absolute() else base_dir / p


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()

    cfg = load_config(args.config)
    window_bp = int(cfg.get("coloc", {}).get("window_bp", 1_000_000))
    eqtl_dir = _resolve(base_dir, cfg.get("gtex", {}).get("eqtl_dir", "data/reference/gtex_v8_eqtl"))

    twas_hits_path = _resolve(base_dir, args.twas_hits)
    harm_gwas_path = _resolve(base_dir, args.harmonized_gwas)
    gene_ann_path = _resolve(base_dir, args.gene_annotation)
    output_dir = _resolve(base_dir, args.output_dir)
    summary_out = _resolve(base_dir, args.summary_out)

    twas_hits = read_tsv_gz(twas_hits_path)
    harm_gwas = read_tsv_gz(harm_gwas_path)
    gene_tss = load_gene_tss_from_gtf(gene_ann_path)

    summary = build_eqtl_regions_for_coloc(
        twas_hits_df=twas_hits,
        harmonized_gwas_df=harm_gwas,
        gene_tss_df=gene_tss,
        eqtl_dir=eqtl_dir,
        output_dir=output_dir,
        window_bp=window_bp,
        chunk_size=args.chunk_size,
    )

    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_out, sep="\t", index=False)

    written = int((summary["status"] == "written").sum()) if not summary.empty else 0
    print(f"Wrote eQTL region files for {written} gene-tissue pairs")
    print(f"Summary: {summary_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
