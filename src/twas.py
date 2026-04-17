"""
twas.py
-------
TWAS stage: aggregate S-PrediXcan results and apply multiple-testing correction.

This module does NOT run S-PrediXcan directly; it:
  - provides helpers that post-process S-PrediXcan CSV outputs
  - aggregates results across tissues
  - applies global BH/FDR correction across all gene-tissue pairs
  - filters to TWAS-significant hits

To generate S-PrediXcan inputs, see src/external_tools/spredixcan_wrapper.py.
"""

from __future__ import annotations

import glob as glob_module
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from src.logging_utils import get_logger

logger = get_logger(__name__)

# Expected columns in a S-PrediXcan output CSV
SPREDIXCAN_REQUIRED_COLS = ["gene", "gene_name", "zscore", "pvalue", "pred_perf_pval", "n_snps_used"]
SPREDIXCAN_EFFECT_COLS   = ["effect_size", "effect_size_se"]  # optional


def load_spredixcan_result(path: str | Path, tissue: str) -> pd.DataFrame:
    """Load a single S-PrediXcan output CSV and annotate with tissue name.

    Parameters
    ----------
    path:
        Path to the S-PrediXcan CSV output file.
    tissue:
        Tissue label to add as a column.

    Returns
    -------
    pd.DataFrame
        Results with an additional ``tissue`` column.
    """
    path = Path(path)
    df = pd.read_csv(path, low_memory=False)
    df["tissue"] = tissue

    # Warn if no SNPs were matched — the "0% model SNPs used" failure mode
    if "n_snps_used" in df.columns:
        zero_snp = (df["n_snps_used"] == 0).sum()
        if zero_snp > 0:
            logger.warning(
                "[%s] %d genes have n_snps_used == 0. "
                "This likely indicates a varID / genome-build mismatch. "
                "Check that the GWAS was harmonized to hg38 with GTEx-style varIDs.",
                tissue,
                zero_snp,
            )

    logger.info("Loaded S-PrediXcan results for %s: %d genes.", tissue, len(df))
    return df


def aggregate_twas_results(result_dir: str | Path, tissues: List[str]) -> pd.DataFrame:
    """Aggregate S-PrediXcan CSV outputs from multiple tissues into one table.

    Parameters
    ----------
    result_dir:
        Directory containing per-tissue CSV files named
        ``{tissue}.spredixcan.csv``.
    tissues:
        List of tissue names to load.

    Returns
    -------
    pd.DataFrame
        Combined TWAS table with columns ``gene``, ``tissue``, ``pvalue``,
        ``zscore``, and more.
    """
    result_dir = Path(result_dir)
    frames = []
    for tissue in tissues:
        pattern = result_dir / f"{tissue}.spredixcan.csv"
        matches = list(result_dir.glob(f"{tissue}.spredixcan.csv"))
        if not matches:
            logger.warning("No S-PrediXcan output found for tissue: %s (expected %s)", tissue, pattern)
            continue
        for f in matches:
            frames.append(load_spredixcan_result(f, tissue))

    if not frames:
        logger.warning("No TWAS results loaded. Returning empty DataFrame.")
        return pd.DataFrame(columns=["gene", "gene_name", "tissue", "pvalue", "zscore"])

    combined = pd.concat(frames, ignore_index=True)
    logger.info(
        "Aggregated TWAS results: %d gene-tissue pairs across %d tissues.",
        len(combined),
        combined["tissue"].nunique(),
    )
    return combined


def apply_fdr_correction(df: pd.DataFrame, alpha: float = 0.05, method: str = "fdr_bh") -> pd.DataFrame:
    """Apply multiple-testing correction globally across all gene-tissue pairs.

    Uses statsmodels' ``multipletests`` for BH/FDR correction.

    Parameters
    ----------
    df:
        Aggregated TWAS DataFrame with a ``pvalue`` column.
    alpha:
        FDR significance threshold (default: 0.05).
    method:
        Correction method (default: ``fdr_bh`` = Benjamini-Hochberg).

    Returns
    -------
    pd.DataFrame
        Input DataFrame with added columns ``pvalue_adj`` and ``significant``.
    """
    if "pvalue" not in df.columns:
        raise ValueError("DataFrame must have a 'pvalue' column for FDR correction.")

    # Drop rows where p-value is missing before correction
    valid_mask = df["pvalue"].notna()
    n_missing = (~valid_mask).sum()
    if n_missing:
        logger.warning("FDR correction: %d gene-tissue pairs have missing p-values.", n_missing)

    if int(valid_mask.sum()) == 0:
        logger.warning("FDR correction skipped: no valid p-values were available.")
        out = df.copy()
        out["pvalue_adj"] = np.nan
        out["significant"] = False
        return out

    try:
        from statsmodels.stats.multitest import multipletests
    except ImportError as exc:
        raise ImportError("statsmodels is required for FDR correction. Install it with: pip install statsmodels") from exc

    pvals = df.loc[valid_mask, "pvalue"].values
    reject, pvals_adj, _, _ = multipletests(pvals, alpha=alpha, method=method)

    df = df.copy()
    df["pvalue_adj"] = np.nan
    df.loc[valid_mask, "pvalue_adj"] = pvals_adj
    df["significant"] = False
    df.loc[valid_mask, "significant"] = reject

    n_sig = reject.sum()
    logger.info(
        "FDR correction (%s, alpha=%.2f): %d / %d gene-tissue pairs significant.",
        method,
        alpha,
        n_sig,
        valid_mask.sum(),
    )
    return df


def filter_twas_hits(df: pd.DataFrame) -> pd.DataFrame:
    """Return only significant gene-tissue pairs (``significant == True``).

    Parameters
    ----------
    df:
        TWAS DataFrame with a ``significant`` column (from ``apply_fdr_correction``).

    Returns
    -------
    pd.DataFrame
        Subset of rows where ``significant is True``.
    """
    if "significant" not in df.columns:
        raise ValueError("Call apply_fdr_correction before filter_twas_hits.")
    hits = df.loc[df["significant"]].copy()
    logger.info("TWAS hits (FDR-significant): %d gene-tissue pairs.", len(hits))
    return hits


def make_mock_spredixcan_result(tissue: str, n_genes: int = 20, seed: int = 42) -> pd.DataFrame:
    """Generate a mock S-PrediXcan output DataFrame for testing / demo.

    Parameters
    ----------
    tissue:
        Tissue name.
    n_genes:
        Number of mock genes to generate.
    seed:
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
    """
    rng = np.random.default_rng(seed)
    genes = [f"ENSG{i:011d}" for i in range(1, n_genes + 1)]
    gene_names = [f"GENE{i}" for i in range(1, n_genes + 1)]
    zscores = rng.normal(0, 2, size=n_genes)
    pvalues = 2 * (1 - _normal_cdf(np.abs(zscores)))

    return pd.DataFrame(
        {
            "gene": genes,
            "gene_name": gene_names,
            "zscore": zscores,
            "pvalue": pvalues,
            "effect_size": zscores * 0.05,
            "effect_size_se": rng.uniform(0.01, 0.1, size=n_genes),
            "pred_perf_pval": rng.uniform(0.0001, 0.05, size=n_genes),
            "pred_perf_r2": rng.uniform(0.01, 0.3, size=n_genes),
            "n_snps_used": rng.integers(5, 50, size=n_genes),
            "tissue": tissue,
        }
    )


def _normal_cdf(x: np.ndarray) -> np.ndarray:
    """Standard-normal CDF (elementwise)."""
    from scipy.special import ndtr  # type: ignore[import]
    return ndtr(x)
