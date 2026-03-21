"""
reporting.py
------------
Generate final summary tables for TWAS + COLOC results.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from src.logging_utils import get_logger

logger = get_logger(__name__)


def build_top_genes_table(
    coloc_summary: pd.DataFrame,
    twas_hits: pd.DataFrame,
    pp4_threshold: float = 0.7,
) -> pd.DataFrame:
    """Join COLOC summary with TWAS hits and filter to coloc-confirmed genes.

    Parameters
    ----------
    coloc_summary:
        Aggregated coloc results (from ``coloc.aggregate_coloc_results``).
    twas_hits:
        FDR-significant TWAS results.
    pp4_threshold:
        Minimum PP4 to include in the top-genes table.

    Returns
    -------
    pd.DataFrame
        Ranked table of coloc-confirmed gene-tissue pairs, sorted by PP4 desc.
    """
    if coloc_summary.empty:
        logger.warning("coloc_summary is empty; returning empty top-genes table.")
        return pd.DataFrame()

    hits = coloc_summary.loc[coloc_summary["PP4"] >= pp4_threshold].copy()
    logger.info(
        "%d gene-tissue pairs pass PP4 >= %.2f threshold.", len(hits), pp4_threshold
    )

    # Merge with TWAS hits for gene name and z-score context
    if not twas_hits.empty and "gene" in twas_hits.columns:
        merge_cols = ["gene", "tissue", "gene_name", "zscore", "pvalue_adj"]
        available = [c for c in merge_cols if c in twas_hits.columns]
        tw = twas_hits[available].rename(columns={"gene": "gene_id"})
        hits = hits.merge(tw, on=["gene_id", "tissue"], how="left")

    hits = hits.sort_values("PP4", ascending=False).reset_index(drop=True)
    return hits


def save_top_genes_table(df: pd.DataFrame, path: str | Path) -> None:
    """Write the top-genes table to a TSV file.

    Parameters
    ----------
    df:
        Top-genes DataFrame.
    path:
        Output path (TSV, uncompressed).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    logger.info("Top-genes table written to %s (%d rows).", path, len(df))


def summarise_run(
    twas_all: pd.DataFrame,
    twas_hits: pd.DataFrame,
    coloc_summary: pd.DataFrame,
    top_genes: pd.DataFrame,
    cfg: Dict[str, Any],
) -> Dict[str, Any]:
    """Return a plain-text summary of the pipeline run for logging / reporting.

    Parameters
    ----------
    twas_all:
        Full aggregated TWAS table (before FDR filter).
    twas_hits:
        FDR-significant TWAS hits.
    coloc_summary:
        All coloc results.
    top_genes:
        Coloc-confirmed top-genes table.
    cfg:
        Full pipeline config.

    Returns
    -------
    Dict[str, Any]
        Summary statistics.
    """
    coloc_cfg = cfg.get("coloc", {})
    correction_cfg = cfg.get("correction", {})

    summary = {
        "n_gene_tissue_tests":    len(twas_all),
        "n_twas_hits_fdr":        len(twas_hits),
        "fdr_alpha":              correction_cfg.get("alpha", 0.05),
        "n_coloc_results":        len(coloc_summary),
        "n_coloc_hits_pp4":       int((coloc_summary["PP4"] >= coloc_cfg.get("pp4_threshold", 0.7)).sum())
                                  if not coloc_summary.empty and "PP4" in coloc_summary.columns else 0,
        "pp4_threshold":          coloc_cfg.get("pp4_threshold", 0.7),
        "n_top_genes":            len(top_genes),
    }

    logger.info("─── Run Summary ───────────────────────────────")
    for k, v in summary.items():
        logger.info("  %-35s %s", k, v)
    logger.info("───────────────────────────────────────────────")

    return summary
