"""
visualization.py
----------------
Plotting utilities for TWAS and COLOC results.

All plots use mock/toy data by default and write PNG files to results/figures/.
Real data will flow through the same functions once the pipeline produces outputs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for CI/HPC environments
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.logging_utils import get_logger

logger = get_logger(__name__)


# ── TWAS QQ-plot ─────────────────────────────────────────────────────────────


def plot_twas_qqplot(
    pvalues: Union[pd.Series, np.ndarray],
    output_path: str | Path,
    title: str = "TWAS p-value QQ-plot",
    dpi: int = 150,
) -> None:
    """Generate a quantile-quantile plot of TWAS p-values.

    Parameters
    ----------
    pvalues:
        Array/Series of p-values (should not contain NaN).
    output_path:
        File path for the output PNG.
    title:
        Plot title.
    dpi:
        Figure resolution.
    """
    pvalues = np.asarray(pvalues, dtype=float)
    pvalues = pvalues[~np.isnan(pvalues)]
    pvalues = np.sort(pvalues)
    n = len(pvalues)

    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    observed = -np.log10(pvalues + 1e-300)  # guard against log(0)

    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    ax.scatter(expected, observed, s=8, alpha=0.6, color="steelblue", label="Gene-tissue pairs")
    max_val = max(expected.max(), observed.max()) + 0.5
    ax.plot([0, max_val], [0, max_val], color="red", lw=1, linestyle="--", label="Expected (null)")

    ax.set_xlabel("Expected $-\\log_{10}(p)$")
    ax.set_ylabel("Observed $-\\log_{10}(p)$")
    ax.set_title(title)
    ax.legend(fontsize=8)
    plt.tight_layout()

    _save_figure(fig, output_path)


# ── TWAS Manhattan-like plot ──────────────────────────────────────────────────


def plot_twas_manhattan(
    twas_df: pd.DataFrame,
    output_path: str | Path,
    pvalue_col: str = "pvalue",
    gene_col: str = "gene_name",
    tissue_col: str = "tissue",
    sig_threshold: float = 0.05,
    title: str = "TWAS association (all gene-tissue pairs)",
    dpi: int = 150,
) -> None:
    """Manhattan-style plot of TWAS p-values ordered by gene index.

    Each point is a gene-tissue pair.  Tissues are colour-coded.

    Parameters
    ----------
    twas_df:
        Aggregated TWAS DataFrame.
    output_path:
        Output PNG path.
    pvalue_col / gene_col / tissue_col:
        Column name overrides.
    sig_threshold:
        FDR-adjusted significance threshold line.
    title:
        Plot title.
    dpi:
        Figure resolution.
    """
    df = twas_df.copy()
    df = df.dropna(subset=[pvalue_col])
    df["neg_log10_p"] = -np.log10(df[pvalue_col].clip(lower=1e-300))

    tissues = df[tissue_col].unique() if tissue_col in df.columns else ["all"]
    cmap = plt.get_cmap("tab10")
    tissue_colors = {t: cmap(i % 10) for i, t in enumerate(tissues)}

    fig, ax = plt.subplots(figsize=(10, 4), dpi=dpi)
    for i, (_, row) in enumerate(df.iterrows()):
        tis = row.get(tissue_col, "all") if tissue_col in df.columns else "all"
        ax.scatter(i, row["neg_log10_p"], s=10, color=tissue_colors[tis], alpha=0.6)

    # Significance threshold line
    sig_line = -np.log10(sig_threshold)
    ax.axhline(sig_line, color="red", lw=1, linestyle="--", label=f"FDR α={sig_threshold}")

    ax.set_xlabel("Gene-tissue pair (index)")
    ax.set_ylabel("$-\\log_{10}(p)$")
    ax.set_title(title)
    ax.legend(fontsize=7)
    plt.tight_layout()

    _save_figure(fig, output_path)


# ── COLOC summary plots ───────────────────────────────────────────────────────


def plot_coloc_pp4_barplot(
    coloc_df: pd.DataFrame,
    output_path: str | Path,
    top_n: int = 20,
    pp4_col: str = "PP4",
    gene_col: str = "gene_id",
    tissue_col: str = "tissue",
    title: str = "Top coloc-confirmed genes (PP4)",
    dpi: int = 150,
) -> None:
    """Horizontal bar chart of top gene-tissue pairs ranked by PP4.

    Parameters
    ----------
    coloc_df:
        Aggregated coloc summary DataFrame.
    output_path:
        Output PNG path.
    top_n:
        Number of top genes to show.
    pp4_col / gene_col / tissue_col:
        Column name overrides.
    title:
        Plot title.
    dpi:
        Figure resolution.
    """
    df = coloc_df.dropna(subset=[pp4_col]).sort_values(pp4_col, ascending=False).head(top_n)
    if df.empty:
        logger.warning("No coloc results to plot in pp4_barplot.")
        return

    labels = [
        f"{row[gene_col]} ({row[tissue_col]})"
        if tissue_col in df.columns
        else str(row[gene_col])
        for _, row in df.iterrows()
    ]

    fig, ax = plt.subplots(figsize=(7, max(3, top_n * 0.3)), dpi=dpi)
    bars = ax.barh(range(len(df)), df[pp4_col].values, color="steelblue", edgecolor="white")
    ax.axvline(0.7, color="red", lw=1, linestyle="--", label="PP4 ≥ 0.7")
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("PP4 (posterior probability of shared causal variant)")
    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.legend(fontsize=8)
    plt.tight_layout()

    _save_figure(fig, output_path)


def plot_coloc_pp4_vs_pp3pp4(
    coloc_df: pd.DataFrame,
    output_path: str | Path,
    pp4_col: str = "PP4",
    pp3_pp4_col: str = "PP3_PP4",
    title: str = "PP4 vs PP3+PP4 (coloc results)",
    dpi: int = 150,
) -> None:
    """Scatter plot of PP4 vs PP3+PP4 for all coloc results.

    Useful for visualising the trade-off between powered loci and coloc hits.

    Parameters
    ----------
    coloc_df:
        Aggregated coloc summary DataFrame.
    output_path:
        Output PNG path.
    pp4_col / pp3_pp4_col:
        Column name overrides.
    title:
        Plot title.
    dpi:
        Figure resolution.
    """
    df = coloc_df.dropna(subset=[pp4_col, pp3_pp4_col])
    if df.empty:
        logger.warning("No coloc results to plot in pp4_vs_pp3pp4.")
        return

    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    ax.scatter(df[pp3_pp4_col], df[pp4_col], s=20, alpha=0.6, color="steelblue")
    ax.axhline(0.7, color="red", lw=1, linestyle="--", label="PP4 ≥ 0.7")
    ax.axvline(0.8, color="orange", lw=1, linestyle="--", label="PP3+PP4 ≥ 0.8")
    ax.set_xlabel("PP3 + PP4 (locus power)")
    ax.set_ylabel("PP4 (shared causal signal)")
    ax.set_title(title)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=8)
    plt.tight_layout()

    _save_figure(fig, output_path)


# ── Helpers ───────────────────────────────────────────────────────────────────


def _save_figure(fig: plt.Figure, output_path: str | Path) -> None:
    """Save a matplotlib figure, creating parent directories as needed."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    logger.info("Figure saved to %s.", output_path)


# ── Mock data for demos ───────────────────────────────────────────────────────


def make_mock_twas_df(n: int = 200, seed: int = 0) -> pd.DataFrame:
    """Generate a mock TWAS DataFrame for plotting demos."""
    rng = np.random.default_rng(seed)
    tissues = ["Whole_Blood", "Spleen", "Colon_Sigmoid", "Muscle_Skeletal"]
    n_per = n // len(tissues)
    frames = []
    for tis in tissues:
        zs = rng.normal(0, 2, n_per)
        ps = 2 * (1 - _norm_cdf(np.abs(zs)))
        frames.append(
            pd.DataFrame(
                {
                    "gene_name": [f"G{i}" for i in range(n_per)],
                    "tissue": tis,
                    "zscore": zs,
                    "pvalue": ps,
                }
            )
        )
    return pd.concat(frames, ignore_index=True)


def make_mock_coloc_df(n: int = 50, seed: int = 1) -> pd.DataFrame:
    """Generate a mock coloc summary DataFrame for plotting demos."""
    rng = np.random.default_rng(seed)
    tissues = ["Whole_Blood", "Spleen", "Colon_Sigmoid"]
    rows = []
    for i in range(n):
        pp4 = float(rng.beta(2, 5))
        pp3 = float(rng.beta(1, 5)) * (1 - pp4)
        rows.append(
            {
                "gene_id": f"ENSG{i:011d}",
                "tissue":  tissues[i % len(tissues)],
                "PP4":     pp4,
                "PP3":     pp3,
                "PP3_PP4": pp3 + pp4,
                "PP4_ratio": pp4 / (pp3 + pp4) if (pp3 + pp4) > 0 else 0.0,
            }
        )
    return pd.DataFrame(rows)


def _norm_cdf(x: np.ndarray) -> np.ndarray:
    from scipy.special import ndtr  # type: ignore[import]
    return ndtr(x)
