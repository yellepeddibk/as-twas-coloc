#!/usr/bin/env python3
"""
build_poster_figures.py
-----------------------
Regenerate all four figures for the ISU Honors Poster at 300 DPI
with gene symbols (not ENSG IDs) and poster-appropriate formatting.

Outputs to docs/poster/figures/
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent
TWAS_DIR = ROOT / "results" / "twas"
COLOC_SUMMARY = ROOT / "data" / "processed" / "coloc_results" / "coloc_summary.tsv.gz"
TWAS_HITS = ROOT / "data" / "processed" / "twas_hits" / "twas_hits_fdr05.tsv.gz"
TOP_GENES = ROOT / "results" / "tables" / "top_genes_pp4.tsv"
OUT_DIR = ROOT / "docs" / "poster" / "figures"

DPI = 300

# ── ISU palette ──────────────────────────────────────────────────────────────
ISU_RED = "#C8102E"
ISU_GOLD = "#F1BE48"
ISU_DARK = "#524727"

# Tissue colour map (8 tissues)
TISSUE_COLORS = {
    "Whole_Blood": "#e41a1c",
    "Spleen": "#377eb8",
    "Colon_Sigmoid": "#4daf4a",
    "Colon_Transverse": "#984ea3",
    "Small_Intestine_Terminal_Ileum": "#ff7f00",
    "Lung": "#a65628",
    "Muscle_Skeletal": "#f781bf",
    "Adipose_Subcutaneous": "#999999",
}


def _poster_style():
    """Set matplotlib rcParams for poster-quality figures."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 14,
        "axes.titlesize": 18,
        "axes.labelsize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def _save(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {path}")


# ── 1. TWAS Manhattan ───────────────────────────────────────────────────────

def build_manhattan():
    """Manhattan-style plot of all gene-tissue pairs with gene symbol labels."""
    print("[1/4] TWAS manhattan plot ...")
    frames = []
    for csv_path in sorted(TWAS_DIR.glob("*.spredixcan.csv")):
        tissue = csv_path.stem.replace(".spredixcan", "")
        df = pd.read_csv(csv_path)
        df["tissue"] = tissue
        frames.append(df)
    all_twas = pd.concat(frames, ignore_index=True)
    all_twas = all_twas.dropna(subset=["pvalue"])
    all_twas["neg_log10_p"] = -np.log10(all_twas["pvalue"].clip(lower=1e-300))

    # Sort by tissue for grouping
    all_twas = all_twas.sort_values(["tissue", "gene_name"]).reset_index(drop=True)

    # FDR threshold: use Bonferroni-ish line or BH
    n_tests = len(all_twas)
    bonf_line = -np.log10(0.05 / n_tests)

    fig, ax = plt.subplots(figsize=(14, 5))

    # Plot by tissue
    tissue_list = sorted(all_twas["tissue"].unique())
    x = 0
    tick_positions = []
    tick_labels = []
    for tissue in tissue_list:
        sub = all_twas[all_twas["tissue"] == tissue]
        xs = np.arange(x, x + len(sub))
        color = TISSUE_COLORS.get(tissue, "#888888")
        ax.scatter(xs, sub["neg_log10_p"].values, s=12, alpha=0.6,
                   color=color, label=tissue.replace("_", " "), rasterized=True)
        mid = x + len(sub) // 2
        tick_positions.append(mid)
        tick_labels.append(tissue.replace("_", "\n"))
        x += len(sub) + 5  # gap between tissues

    # Significance line
    # Use FDR 0.05 on hits file to get the actual threshold
    twas_hits = pd.read_csv(TWAS_HITS, sep="\t")
    if len(twas_hits) > 0:
        max_hit_p = twas_hits["pvalue"].max()
        fdr_line = -np.log10(max_hit_p)
        ax.axhline(fdr_line, color=ISU_RED, lw=1.5, linestyle="--",
                    label=f"FDR < 0.05 ({len(twas_hits)} hits)", zorder=1)

    # Label top hits
    top_genes = all_twas.nlargest(8, "neg_log10_p")
    for _, row in top_genes.iterrows():
        idx = row.name  # original index
        # Find the x position for this point
        tissue = row["tissue"]
        tissue_start = 0
        for t in tissue_list:
            if t == tissue:
                break
            tissue_start += len(all_twas[all_twas["tissue"] == t]) + 5
        within_tissue = all_twas[all_twas["tissue"] == tissue].index.get_loc(idx)
        x_pos = tissue_start + within_tissue
        gene_label = row["gene_name"] if pd.notna(row.get("gene_name")) else row["gene"]
        ax.annotate(
            gene_label,
            (x_pos, row["neg_log10_p"]),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center", fontsize=9, fontstyle="italic",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
        )

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=9)
    ax.set_ylabel("$-\\log_{10}(p)$")
    ax.set_title("Transcriptome-Wide Association Study: All Gene-Tissue Pairs")
    ax.legend(loc="upper right", ncol=2, fontsize=8, framealpha=0.9)
    fig.tight_layout()
    _save(fig, "manhattan_twas.png")


# ── 2. QQ-plot ───────────────────────────────────────────────────────────────

def build_qqplot():
    """QQ-plot of all TWAS p-values."""
    print("[2/4] TWAS QQ-plot ...")
    frames = []
    for csv_path in sorted(TWAS_DIR.glob("*.spredixcan.csv")):
        df = pd.read_csv(csv_path)
        frames.append(df)
    all_twas = pd.concat(frames, ignore_index=True)
    pvals = all_twas["pvalue"].dropna().values
    pvals = np.sort(pvals)
    n = len(pvals)

    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    observed = -np.log10(pvals + 1e-300)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=14, alpha=0.5, color="steelblue",
               edgecolors="none", label=f"Gene-tissue pairs (n={n:,})")
    max_val = max(expected.max(), observed.max()) + 0.5
    ax.plot([0, max_val], [0, max_val], color=ISU_RED, lw=1.5,
            linestyle="--", label="Expected (null)")

    # Confidence band (95%)
    ci_upper = []
    ci_lower = []
    for i in range(1, n + 1):
        from scipy.stats import beta as beta_dist
        ci_lower.append(-np.log10(beta_dist.ppf(0.975, i, n - i + 1)))
        ci_upper.append(-np.log10(beta_dist.ppf(0.025, i, n - i + 1)))
    ax.fill_between(expected, ci_lower, ci_upper, alpha=0.15, color="grey",
                    label="95% CI")

    ax.set_xlabel("Expected $-\\log_{10}(p)$")
    ax.set_ylabel("Observed $-\\log_{10}(p)$")
    ax.set_title("TWAS P-value QQ-Plot")
    ax.legend(fontsize=10)
    fig.tight_layout()
    _save(fig, "qqplot_twas.png")


# ── 3. COLOC PP4 bar chart ──────────────────────────────────────────────────

def build_coloc_barplot():
    """Horizontal bar chart of top coloc gene-tissue pairs by PP4."""
    print("[3/4] COLOC PP4 bar plot ...")
    coloc = pd.read_csv(COLOC_SUMMARY, sep="\t")
    top = pd.read_csv(TOP_GENES, sep="\t")

    # Merge gene_name into coloc
    if "gene_name" not in coloc.columns:
        gene_map = top[["gene_id", "gene_name"]].drop_duplicates()
        coloc = coloc.merge(gene_map, on="gene_id", how="left")
        coloc["gene_name"] = coloc["gene_name"].fillna(coloc["gene_id"])

    hits = coloc[coloc["PP4"] >= 0.70].sort_values("PP4", ascending=True).copy()
    hits["label"] = hits["gene_name"] + " (" + hits["tissue"].str.replace("_", " ") + ")"

    fig, ax = plt.subplots(figsize=(8, max(4, len(hits) * 0.35)))
    colors = [TISSUE_COLORS.get(t, "#888888") for t in hits["tissue"]]
    ax.barh(range(len(hits)), hits["PP4"].values, color=colors,
            edgecolor="white", linewidth=0.5)
    ax.axvline(0.70, color=ISU_RED, lw=1.5, linestyle="--", label="PP4 = 0.70")
    ax.set_yticks(range(len(hits)))
    ax.set_yticklabels(hits["label"].values, fontsize=10)
    ax.set_xlabel("PP4 (Posterior Probability of Shared Causal Variant)")
    ax.set_title(f"Colocalization-Confirmed Gene-Tissue Pairs (n={len(hits)})")
    ax.set_xlim(0.5, 1.02)
    ax.legend(fontsize=10, loc="lower right")
    fig.tight_layout()
    _save(fig, "coloc_pp4_barplot.png")


# ── 4. PP4 vs PP3+PP4 scatter ───────────────────────────────────────────────

def build_pp4_scatter():
    """Scatter of PP4 vs PP3+PP4 with top genes labelled."""
    print("[4/4] PP4 vs PP3+PP4 scatter ...")
    coloc = pd.read_csv(COLOC_SUMMARY, sep="\t")
    top = pd.read_csv(TOP_GENES, sep="\t")

    # Merge gene_name
    if "gene_name" not in coloc.columns:
        gene_map = top[["gene_id", "gene_name"]].drop_duplicates()
        coloc = coloc.merge(gene_map, on="gene_id", how="left")
        coloc["gene_name"] = coloc["gene_name"].fillna(coloc["gene_id"])

    fig, ax = plt.subplots(figsize=(7, 7))

    # Background points
    non_hits = coloc[coloc["PP4"] < 0.70]
    ax.scatter(non_hits["PP3_PP4"], non_hits["PP4"], s=30, alpha=0.35,
               color="#AAAAAA", edgecolors="none", label="Non-hits", zorder=2)

    # Hits
    hits = coloc[coloc["PP4"] >= 0.70]
    colors = [TISSUE_COLORS.get(t, "#888888") for t in hits["tissue"]]
    ax.scatter(hits["PP3_PP4"], hits["PP4"], s=60, alpha=0.85,
               c=colors, edgecolors="black", linewidth=0.5,
               label=f"Coloc hits (n={len(hits)})", zorder=3)

    # Threshold lines
    ax.axhline(0.70, color=ISU_RED, lw=1.5, linestyle="--", label="PP4 = 0.70", zorder=1)
    ax.axvline(0.80, color=ISU_GOLD, lw=1.5, linestyle="--", label="PP3+PP4 = 0.80", zorder=1)

    # Quadrant shading
    ax.fill_between([0.80, 1.05], 0.70, 1.05, alpha=0.06, color=ISU_RED, zorder=0)

    # Label top 10 hits
    top_hits = hits.nlargest(10, "PP4")
    for _, row in top_hits.iterrows():
        ax.annotate(
            row["gene_name"],
            (row["PP3_PP4"], row["PP4"]),
            textcoords="offset points",
            xytext=(8, 4),
            fontsize=9, fontstyle="italic",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
        )

    ax.set_xlabel("PP3 + PP4 (Locus Power)")
    ax.set_ylabel("PP4 (Shared Causal Signal)")
    ax.set_title("Colocalization Evidence: PP4 vs Locus Power")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=10, loc="upper left")
    fig.tight_layout()
    _save(fig, "coloc_pp4_vs_pp3pp4.png")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    _poster_style()
    print(f"Output directory: {OUT_DIR}")
    build_manhattan()
    build_qqplot()
    build_coloc_barplot()
    build_pp4_scatter()
    print("\nAll poster figures generated.")


if __name__ == "__main__":
    main()
