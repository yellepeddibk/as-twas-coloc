"""
coloc.py
--------
COLOC stage: build coloc-ready datasets and interpret posterior probabilities.

Key design decisions documented here
-------------------------------------
* Use beta + varbeta (= se**2), NOT z-scores, for coloc.abf.
  coloc's ABF formula requires varbeta explicitly.
* Do NOT pre-filter locus SNPs by p-value before running coloc.
  Removing SNPs biases the posterior probabilities toward H4.
* Operate on one genomic region at a time (1 gene × 1 tissue at a time).
* Store priors (p1, p2, p12) in config and record them in every output JSON.
* Support coloc.abf as the primary scaffold path.
* Placeholder support for coloc.susie for multi-causal regions (not enabled by default).

Coloc interpretation metrics used in this pipeline
----------------------------------------------------
PP4      >= 0.7   -> coloc hit (shared causal signal)
PP3+PP4  >= 0.8   -> powered locus (enough signal to distinguish H3/H4)
PP4/(PP3+PP4) >= 0.9 -> strong preference for shared signal over distinct signals

References
----------
Giambartolomei et al. (2014) PLOS Genetics — coloc.abf
Wallace (2021) PLOS Genetics — coloc.susie
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from src.logging_utils import get_logger

logger = get_logger(__name__)


# ── Varbeta computation ───────────────────────────────────────────────────────


def compute_varbeta(df: pd.DataFrame, se_col: str = "se") -> pd.DataFrame:
    """Add a ``varbeta`` column as ``se ** 2``.

    This is required by coloc.abf.  Never pass se directly to coloc.

    Parameters
    ----------
    df:
        DataFrame with a standard-error column.
    se_col:
        Name of the SE column (default: ``se``).

    Returns
    -------
    pd.DataFrame
        Copy of *df* with a ``varbeta`` column added.
    """
    if se_col not in df.columns:
        raise ValueError(f"SE column '{se_col}' not found in DataFrame.")
    df = df.copy()
    df["varbeta"] = df[se_col] ** 2
    logger.debug("Computed varbeta = %s ** 2 for %d rows.", se_col, len(df))
    return df


# ── Locus extraction ──────────────────────────────────────────────────────────


def extract_gwas_locus(
    gwas_df: pd.DataFrame,
    chrom: str,
    center_pos: int,
    window_bp: int = 1_000_000,
) -> pd.DataFrame:
    """Extract GWAS summary statistics for a ±window_bp locus.

    Parameters
    ----------
    gwas_df:
        Harmonized GWAS DataFrame.
    chrom:
        Chromosome (no "chr" prefix).
    center_pos:
        Center position (e.g., gene TSS) in hg38 bp.
    window_bp:
        Half-window size in bp (default: 1 Mb each side).

    Returns
    -------
    pd.DataFrame
        Subset of rows within the locus.
    """
    chrom = str(chrom).replace("chr", "")
    mask = (
        (gwas_df["chrom"].astype(str) == chrom)
        & (gwas_df["pos"] >= center_pos - window_bp)
        & (gwas_df["pos"] <= center_pos + window_bp)
    )
    locus = gwas_df.loc[mask].copy()
    logger.debug(
        "GWAS locus chr%s ±%d bp around %d: %d SNPs.",
        chrom,
        window_bp,
        center_pos,
        len(locus),
    )
    return locus


def build_coloc_dataset(
    gwas_locus: pd.DataFrame,
    eqtl_locus: pd.DataFrame,
    n_gwas: int,
    n_cases: Optional[int],
    n_controls: Optional[int],
    s_fraction: Optional[float] = None,
    dataset_type: str = "cc",
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Build coloc dataset1 (GWAS) and dataset2 (eQTL) dictionaries.

    These dicts are passed to the R coloc.abf wrapper.

    Parameters
    ----------
    gwas_locus:
        GWAS summary stats for the locus.  Must contain
        ``varID``, ``beta``, ``se``, ``pvalue``.
    eqtl_locus:
        eQTL summary stats for the locus.  Must contain
        ``varID``, ``beta``, ``se``, ``pvalue``, ``N``.
    n_gwas:
        Total GWAS sample size.
    n_cases / n_controls:
        Case/control counts (for case-control GWAS).  Used to compute
        the ``s`` parameter (case fraction) required by coloc for cc studies.
    s_fraction:
        Override for case fraction.  If None, computed from n_cases/n_gwas.
    dataset_type:
        ``"cc"`` (case-control) or ``"quant"`` (quantitative trait).

    Returns
    -------
    Tuple[Dict, Dict]
        (dataset1, dataset2) ready for R coloc.abf.
    """
    # ── Align on shared varIDs ─────────────────────────────────
    shared_ids = set(gwas_locus["varID"]) & set(eqtl_locus["varID"])
    if len(shared_ids) == 0:
        logger.warning(
            "No shared varIDs between GWAS locus and eQTL locus. "
            "Check genome build and varID construction."
        )
        return {}, {}

    g = gwas_locus.set_index("varID").loc[list(shared_ids)].copy()
    e = eqtl_locus.set_index("varID").loc[list(shared_ids)].copy()

    # Ensure consistent SNP ordering
    g = g.sort_index()
    e = e.sort_index()

    # ── Compute varbeta ────────────────────────────────────────
    g["varbeta"] = g["se"] ** 2
    e["varbeta"] = e["se"] ** 2

    # ── Case fraction ──────────────────────────────────────────
    if s_fraction is None and n_cases is not None and n_gwas > 0:
        s_fraction = n_cases / n_gwas
    elif s_fraction is None:
        # TODO: update s_fraction once real case/control counts are confirmed
        logger.warning(
            "Case fraction (s) not computable; n_cases or n_gwas missing.  "
            "Using placeholder s=0.5.  Update config with real case/control counts."
        )
        s_fraction = 0.5

    # ── Dataset dicts ──────────────────────────────────────────
    snp_ids = g.index.tolist()

    dataset1: Dict[str, Any] = {
        "snp":     snp_ids,
        "beta":    g["beta"].tolist(),
        "varbeta": g["varbeta"].tolist(),
        "N":       n_gwas,
        "type":    dataset_type,
    }
    if dataset_type == "cc":
        dataset1["s"] = s_fraction

    dataset2: Dict[str, Any] = {
        "snp":     snp_ids,
        "beta":    e["beta"].tolist(),
        "varbeta": e["varbeta"].tolist(),
        "N":       int(e["N"].median()) if "N" in e.columns else 200,  # GTEx typical N
        "type":    "quant",  # eQTL is always quantitative
    }

    logger.info(
        "Built coloc datasets: %d shared SNPs, dataset1 type=%s, s=%.3f.",
        len(snp_ids),
        dataset_type,
        s_fraction if dataset_type == "cc" else float("nan"),
    )
    return dataset1, dataset2


# ── Coloc result interpretation ───────────────────────────────────────────────


def interpret_coloc_result(result: Dict[str, Any], cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Annotate a raw coloc.abf result dict with interpretation flags.

    Parameters
    ----------
    result:
        Dict with keys ``PP.H0.abf``, ``PP.H1.abf``, ``PP.H2.abf``,
        ``PP.H3.abf``, ``PP.H4.abf`` (as returned by coloc.abf in R).
    cfg:
        Full pipeline config (reads thresholds from ``coloc`` section).

    Returns
    -------
    Dict[str, Any]
        Annotated result dict with derived metrics and interpretation flags.
    """
    coloc_cfg = cfg.get("coloc", {})
    pp4_thresh  = coloc_cfg.get("pp4_threshold", 0.7)
    pp34_thresh = coloc_cfg.get("pp3_pp4_threshold", 0.8)
    ratio_thresh = coloc_cfg.get("pp4_ratio_threshold", 0.9)

    pp3 = result.get("PP.H3.abf", 0.0)
    pp4 = result.get("PP.H4.abf", 0.0)
    pp3_pp4 = pp3 + pp4
    ratio = pp4 / pp3_pp4 if pp3_pp4 > 0 else 0.0

    annotated = dict(result)
    annotated["PP3_PP4"]         = pp3_pp4
    annotated["PP4_ratio"]       = ratio
    annotated["coloc_hit"]       = pp4 >= pp4_thresh
    annotated["powered_locus"]   = pp3_pp4 >= pp34_thresh
    annotated["strong_shared"]   = ratio >= ratio_thresh
    annotated["priors_used"]     = {
        "p1":  coloc_cfg.get("p1",  1e-4),
        "p2":  coloc_cfg.get("p2",  1e-4),
        "p12": coloc_cfg.get("p12", 1e-5),
    }

    logger.debug(
        "Coloc interpretation: PP4=%.3f, PP3+PP4=%.3f, ratio=%.3f, "
        "coloc_hit=%s, powered=%s, strong=%s.",
        pp4, pp3_pp4, ratio,
        annotated["coloc_hit"],
        annotated["powered_locus"],
        annotated["strong_shared"],
    )
    return annotated


def save_coloc_result(result: Dict[str, Any], path: str | Path) -> None:
    """Serialize a coloc result dict to JSON.

    Parameters
    ----------
    result:
        Coloc result dict (from interpret_coloc_result).
    path:
        Output file path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        json.dump(result, fh, indent=2, default=_json_default)
    logger.debug("Saved coloc result to %s.", path)


def _json_default(obj: Any) -> Any:
    """JSON serialization fallback for numpy types."""
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


# ── Aggregation ───────────────────────────────────────────────────────────────


def aggregate_coloc_results(coloc_dir: str | Path) -> pd.DataFrame:
    """Collect all per-gene coloc JSON files and build a summary table.

    Parameters
    ----------
    coloc_dir:
        Top-level directory containing ``{tissue}/{gene_id}.json`` files.

    Returns
    -------
    pd.DataFrame
        Summary table with one row per gene-tissue pair, sorted by PP4 descending.
    """
    coloc_dir = Path(coloc_dir)
    rows = []
    for json_path in sorted(coloc_dir.rglob("*.json")):
        try:
            with json_path.open() as fh:
                res = json.load(fh)
        except Exception as exc:
            logger.warning("Could not load coloc JSON %s: %s", json_path, exc)
            continue

        # Infer gene and tissue from path structure: {tissue}/{gene_id}.json
        gene_id = json_path.stem
        tissue  = json_path.parent.name
        row = {
            "gene_id": gene_id,
            "tissue":  tissue,
            "PP0": res.get("PP.H0.abf"),
            "PP1": res.get("PP.H1.abf"),
            "PP2": res.get("PP.H2.abf"),
            "PP3": res.get("PP.H3.abf"),
            "PP4": res.get("PP.H4.abf"),
            "PP3_PP4":    res.get("PP3_PP4"),
            "PP4_ratio":  res.get("PP4_ratio"),
            "coloc_hit":  res.get("coloc_hit"),
            "powered_locus": res.get("powered_locus"),
            "strong_shared": res.get("strong_shared"),
        }
        rows.append(row)

    if not rows:
        logger.warning("No coloc JSON files found in %s.", coloc_dir)
        return pd.DataFrame()

    df = pd.DataFrame(rows).sort_values("PP4", ascending=False).reset_index(drop=True)
    logger.info("Aggregated %d coloc results.", len(df))
    return df


# ── Mock coloc result (for testing/demo) ─────────────────────────────────────


def make_mock_coloc_result(pp4: float = 0.85, pp3: float = 0.10) -> Dict[str, Any]:
    """Return a mock coloc.abf result dict for testing.

    Parameters
    ----------
    pp4:
        Posterior probability of hypothesis H4 (shared causal variant).
    pp3:
        Posterior probability of H3 (distinct causal variants).

    Returns
    -------
    Dict[str, Any]
    """
    total_remaining = max(0.0, 1.0 - pp4 - pp3)
    return {
        "PP.H0.abf": round(total_remaining * 0.01, 6),
        "PP.H1.abf": round(total_remaining * 0.33, 6),
        "PP.H2.abf": round(total_remaining * 0.66, 6),
        "PP.H3.abf": pp3,
        "PP.H4.abf": pp4,
        "nsnps":     500,
    }
