"""
manuscript_support.py
---------------------
Utilities that bridge pipeline outputs to the manuscript scaffold.

This module helps:
- generate Table 1 (cohort characteristics) placeholder
- format top-genes tables for LaTeX inclusion
- produce figure captions
- check which result files exist and whether manuscript sections need updating

Do NOT put biological conclusions or numerical claims here until real analyses are done.
All auto-generated content is explicitly labeled as [PLACEHOLDER].
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from src.logging_utils import get_logger

logger = get_logger(__name__)


# ── LaTeX table helpers ───────────────────────────────────────────────────────


def df_to_latex(
    df: pd.DataFrame,
    caption: str = "",
    label: str = "",
    float_fmt: str = "%.3f",
    na_rep: str = "—",
    index: bool = False,
) -> str:
    """Convert a DataFrame to a LaTeX tabular string.

    Parameters
    ----------
    df:
        DataFrame to convert.
    caption:
        Table caption text.
    label:
        LaTeX label (for \\ref{}).
    float_fmt:
        Format string for floating-point values.
    na_rep:
        String for NaN values.
    index:
        Whether to include the DataFrame index.

    Returns
    -------
    str
        LaTeX table environment string.
    """
    latex_str = df.to_latex(
        index=index,
        float_format=float_fmt,
        na_rep=na_rep,
        escape=True,
    )
    if caption or label:
        header = "\\begin{table}[htbp]\n\\centering\n"
        if caption:
            header += f"\\caption{{{caption}}}\n"
        if label:
            header += f"\\label{{{label}}}\n"
        footer = "\\end{table}\n"
        return header + latex_str + footer

    return latex_str


def format_top_genes_for_manuscript(
    top_genes: pd.DataFrame,
    cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Select and rename columns for manuscript Table 2 (top coloc genes).

    Parameters
    ----------
    top_genes:
        Full top-genes DataFrame from ``reporting.build_top_genes_table``.
    cols:
        Subset of columns to include; default selects the most informative.

    Returns
    -------
    pd.DataFrame
        Formatted table for manuscript inclusion.
    """
    default_cols = [
        "gene_id", "gene_name", "tissue",
        "PP4", "PP3_PP4", "PP4_ratio",
        "coloc_hit", "powered_locus",
        "zscore", "pvalue_adj",
    ]
    cols = cols or [c for c in default_cols if c in top_genes.columns]
    df = top_genes[cols].copy()

    rename = {
        "gene_id":      "Gene ID",
        "gene_name":    "Gene",
        "tissue":       "Tissue",
        "PP4":          "PP4",
        "PP3_PP4":      "PP3+PP4",
        "PP4_ratio":    "PP4/(PP3+PP4)",
        "coloc_hit":    "Coloc hit",
        "powered_locus": "Powered",
        "zscore":       "TWAS Z",
        "pvalue_adj":   "TWAS FDR q",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    return df


# ── Artifact checklist ────────────────────────────────────────────────────────


EXPECTED_ARTIFACTS: Dict[str, str] = {
    "data/interim/gwas_harmonized/as_gwas_hg38_varid.tsv.gz": "Harmonized GWAS",
    "data/processed/twas_hits/twas_hits_fdr05.tsv.gz":        "FDR-significant TWAS hits",
    "data/processed/coloc_results/coloc_summary.tsv.gz":       "COLOC summary table",
    "results/tables/top_genes_pp4.tsv":                        "Top coloc genes (PP4 ≥ 0.7)",
    "results/figures/qqplot_twas.png":                         "TWAS QQ-plot",
    "results/figures/manhattan_like_twas.png":                 "TWAS Manhattan-style plot",
    "results/figures/coloc_summary.png":                       "COLOC summary plot",
}


def check_artifact_status(base_dir: str | Path = ".") -> Dict[str, bool]:
    """Return a dict of artifact existence status.

    Parameters
    ----------
    base_dir:
        Repository root directory.

    Returns
    -------
    Dict[str, bool]
        Mapping from artifact path to existence flag.
    """
    base_dir = Path(base_dir)
    status = {}
    for rel_path, label in EXPECTED_ARTIFACTS.items():
        full = base_dir / rel_path
        exists = full.exists()
        status[rel_path] = exists
        symbol = "✓" if exists else "✗"
        logger.info("[%s] %s — %s", symbol, label, rel_path)

    return status


def print_manuscript_readiness(base_dir: str | Path = ".") -> None:
    """Print a readiness report for manuscript figure/table insertion."""
    status = check_artifact_status(base_dir)
    missing = [p for p, exists in status.items() if not exists]
    present = [p for p, exists in status.items() if exists]

    print("\n── Manuscript artifact readiness ───────────────────────────")
    print(f"  Present ({len(present)}):")
    for p in present:
        print(f"    ✓  {p}")
    print(f"  Missing ({len(missing)}) — [PLACEHOLDER] sections need updating:")
    for p in missing:
        print(f"    ✗  {p}")
    print("────────────────────────────────────────────────────────────\n")
