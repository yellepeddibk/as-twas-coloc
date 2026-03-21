"""
pipeline.py
-----------
Top-level pipeline orchestrator for the AS TWAS + COLOC pipeline.

Pipeline stages
---------------
1.  Load raw GWAS summary statistics
2.  Harmonize GWAS
3.  Save harmonized GWAS to data/interim/
4.  Run S-PrediXcan (command-builder / dry-run in mock mode)
5.  Aggregate TWAS results across tissues
6.  Apply BH/FDR correction globally
7.  Filter TWAS-significant hits
8.  (Placeholder) Extract locus-level eQTL regional data for COLOC
9.  (Placeholder) Build COLOC-ready datasets
10. Run COLOC (mock by default; real when R + coloc installed)
11. Rank genes by PP4, PP3+PP4, PP4/(PP3+PP4)
12. Output ranked coloc-prioritized genes
13. Generate summary plots
14. Save provenance manifest

In mock mode (``mock=True``), external tools are not invoked:
- S-PrediXcan commands are built but not executed
- TWAS results are generated synthetically for demonstration
- COLOC is run on synthetic data via the mock interface
"""

from __future__ import annotations

import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from src.config import get_coloc_priors, get_tissues, load_config
from src.io_utils import ensure_dir, read_gwas, read_tsv_gz, write_tsv_gz
from src.harmonization import harmonize_gwas
from src.twas import (
    aggregate_twas_results,
    apply_fdr_correction,
    filter_twas_hits,
    make_mock_spredixcan_result,
)
from src.coloc import (
    aggregate_coloc_results,
    interpret_coloc_result,
    make_mock_coloc_result,
    save_coloc_result,
)
from src.reporting import build_top_genes_table, save_top_genes_table, summarise_run
from src.visualization import (
    make_mock_coloc_df,
    make_mock_twas_df,
    plot_coloc_pp4_barplot,
    plot_coloc_pp4_vs_pp3pp4,
    plot_twas_manhattan,
    plot_twas_qqplot,
)
from src.provenance import build_manifest, save_manifest
from src.logging_utils import get_logger
from src.external_tools.spredixcan_wrapper import build_spredixcan_command, run_spredixcan
from src.external_tools.coloc_r_wrapper import run_coloc_abf
from src.validation import check_gtex_model_varid_overlap

logger = get_logger(__name__)

# Run-time identifier (used in manifest filename and output naming)
_RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def run_pipeline(
    config_path: Optional[str | Path] = None,
    mock: bool = True,
    base_dir: str | Path = ".",
) -> Dict[str, Any]:
    """Execute the full TWAS + COLOC pipeline.

    Parameters
    ----------
    config_path:
        Path to YAML config.  If None, uses ``config/as.yaml``.
    mock:
        Run in mock/demo mode (no external tools, synthetic data).
    base_dir:
        Repository root (used for relative path resolution).

    Returns
    -------
    Dict[str, Any]
        Summary dict with counts and output paths.
    """
    base_dir = Path(base_dir).resolve()
    cfg = load_config(config_path)
    cfg["_base_dir"] = str(base_dir)
    paths = cfg.get("paths", {})

    # ── Create output directories ─────────────────────────────
    for key in (
        "data_interim", "data_processed",
        "results_twas", "results_coloc",
        "results_tables", "results_figures", "results_manifests",
    ):
        ensure_dir(base_dir / paths.get(key, key))

    logger.info("=" * 60)
    logger.info("AS TWAS + COLOC Pipeline")
    logger.info("Mock mode: %s", mock)
    logger.info("Run ID:    %s", _RUN_ID)
    logger.info("=" * 60)

    # ── Stage 1: Load raw GWAS ────────────────────────────────
    logger.info("Stage 1: Load raw GWAS")
    gwas_input = _resolve_gwas_input_path(base_dir, cfg)

    if mock:
        if not gwas_input.exists():
            logger.info("[MOCK] GWAS input not found; generating synthetic GWAS data.")
            raw_gwas = _make_mock_gwas_df(cfg)
        else:
            raw_gwas = read_gwas(gwas_input)
    else:
        if not gwas_input.exists():
            raise FileNotFoundError(
                f"Real-mode GWAS input not found: {gwas_input}. "
                "Provide a real summary statistics file or run with --mock."
            )
        raw_gwas = read_gwas(gwas_input)

    # ── Stage 2: Harmonize GWAS ───────────────────────────────
    logger.info("Stage 2: Harmonize GWAS")
    harmonized_gwas = harmonize_gwas(raw_gwas, cfg, strict_real_mode=not mock)

    # ── Stage 3: Save harmonized GWAS ─────────────────────────
    logger.info("Stage 3: Save harmonized GWAS")
    harm_dir = base_dir / paths.get("data_interim", "data/interim") / "gwas_harmonized"
    harm_path = harm_dir / "as_gwas_hg38_varid.tsv.gz"
    write_tsv_gz(harmonized_gwas, harm_path)

    # ── Stage 3b: Save harmonization QC summary ──────────────
    qc_dir = base_dir / paths.get("data_processed", "data/processed") / "qc"
    ensure_dir(qc_dir)
    qc_summary = _build_harmonization_qc(harmonized_gwas)
    qc_path = qc_dir / "gwas_harmonization_qc.tsv"
    pd.DataFrame([qc_summary]).to_csv(qc_path, sep="\t", index=False)
    logger.info("Saved GWAS harmonization QC summary to %s", qc_path)

    # ── Stage 3c: GTEx overlap sanity check ──────────────────
    overlap = _run_gtex_overlap_check(harmonized_gwas, cfg, mock=mock)
    if overlap is not None:
        overlap_tissue = overlap.get("tissue", "unknown_tissue")
        overlap_out = qc_dir / f"gwas_gtex_overlap_{overlap_tissue}.tsv"
        pd.DataFrame([overlap]).to_csv(overlap_out, sep="\t", index=False)
        logger.info("Saved GTEx overlap sanity check to %s", overlap_out)

    # ── Stage 4: Run S-PrediXcan ──────────────────────────────
    logger.info("Stage 4: S-PrediXcan (mock=%s)", mock)
    tissues = get_tissues(cfg)
    twas_dir = base_dir / paths.get("results_twas", "results/twas")

    if mock:
        # Write synthetic S-PrediXcan outputs so aggregation works
        for tissue in tissues:
            mock_df = make_mock_spredixcan_result(tissue)
            mock_out = twas_dir / f"{tissue}.spredixcan.csv"
            mock_df.to_csv(mock_out, index=False)
            logger.info("[MOCK] Wrote synthetic TWAS for %s.", tissue)
    else:
        # Build and execute real S-PrediXcan commands
        from src.external_tools.spredixcan_wrapper import run_all_tissues
        run_all_tissues(cfg, tissues=tissues, dry_run=False)

    # ── Stage 5: Aggregate TWAS ───────────────────────────────
    logger.info("Stage 5: Aggregate TWAS results")
    twas_all = aggregate_twas_results(twas_dir, tissues)

    # ── Stage 6: FDR correction ───────────────────────────────
    logger.info("Stage 6: FDR correction")
    correction_cfg = cfg.get("correction", {})
    twas_all = apply_fdr_correction(
        twas_all,
        alpha=correction_cfg.get("alpha", 0.05),
        method=correction_cfg.get("method", "fdr_bh"),
    )

    # ── Stage 7: Filter TWAS hits ─────────────────────────────
    logger.info("Stage 7: Filter TWAS-significant hits")
    twas_hits = filter_twas_hits(twas_all)
    hits_path = base_dir / paths.get("data_processed", "data/processed") / "twas_hits" / "twas_hits_fdr05.tsv.gz"
    write_tsv_gz(twas_hits, hits_path)

    # ── Stages 8-10: COLOC ────────────────────────────────────
    logger.info("Stages 8-10: COLOC (mock=%s)", mock)
    coloc_priors = get_coloc_priors(cfg)
    coloc_out_dir = base_dir / paths.get("results_coloc", "results/coloc")

    if mock:
        _run_mock_coloc(twas_hits, cfg, coloc_out_dir)
    else:
        raise RuntimeError(
            "Real-mode COLOC path is not fully wired yet. "
            "No mock fallback is allowed in real mode. "
            "Complete eQTL regional extraction + coloc wiring before running --no-mock end-to-end."
        )

    # ── Stage 11-12: Rank and output genes ────────────────────
    logger.info("Stages 11-12: Aggregate and rank coloc results")
    coloc_summary = aggregate_coloc_results(coloc_out_dir)

    if not coloc_summary.empty:
        coloc_summary_path = (
            base_dir / paths.get("data_processed", "data/processed")
            / "coloc_results" / "coloc_summary.tsv.gz"
        )
        write_tsv_gz(coloc_summary, coloc_summary_path)

    top_genes = build_top_genes_table(
        coloc_summary,
        twas_hits,
        pp4_threshold=cfg.get("coloc", {}).get("pp4_threshold", 0.7),
    )
    top_genes_path = base_dir / paths.get("results_tables", "results/tables") / "top_genes_pp4.tsv"
    save_top_genes_table(top_genes, top_genes_path)

    # ── Stage 13: Plots ───────────────────────────────────────
    logger.info("Stage 13: Generate plots")
    fig_dir = base_dir / paths.get("results_figures", "results/figures")
    _generate_plots(twas_all, coloc_summary, fig_dir, cfg, mock=mock)

    # ── Stage 14: Provenance manifest ─────────────────────────
    logger.info("Stage 14: Save provenance manifest")
    manifest = build_manifest(cfg, input_files=[harm_path])
    save_manifest(
        manifest,
        manifest_dir=base_dir / paths.get("results_manifests", "results/manifests"),
        run_id=_RUN_ID,
    )

    summary = summarise_run(twas_all, twas_hits, coloc_summary, top_genes, cfg)
    logger.info("Pipeline complete.")
    return summary


# ── Internal helpers ──────────────────────────────────────────────────────────


def _make_mock_gwas_df(cfg: Dict[str, Any], n_snps: int = 200) -> pd.DataFrame:
    """Generate a minimal mock GWAS DataFrame for demonstration."""
    import numpy as np

    gwas_cfg = cfg.get("gwas", {})
    rng = np.random.default_rng(42)

    chroms = rng.choice([str(c) for c in range(1, 23)], size=n_snps)
    positions = rng.integers(1_000_000, 250_000_000, size=n_snps)
    betas = rng.normal(0, 0.05, size=n_snps)
    ses = rng.uniform(0.02, 0.1, size=n_snps)
    pvals = 2 * (1 - _norm_cdf_arr(abs(betas / ses)))
    alleles_a = rng.choice(["A", "C", "G", "T"], size=n_snps)
    # Avoid palindromic by ensuring alt != complement
    _comp = {"A": "C", "C": "A", "G": "T", "T": "G"}
    alleles_b = [_comp[a] for a in alleles_a]

    return pd.DataFrame(
        {
            "CHR":  chroms,
            "BP":   positions,
            "A1":   alleles_a,
            "A2":   alleles_b,
            "BETA": betas,
            "SE":   ses,
            "P":    pvals,
            "N":    gwas_cfg.get("n_total", 10619),
        }
    )


def _run_mock_coloc(
    twas_hits: pd.DataFrame,
    cfg: Dict[str, Any],
    coloc_out_dir: Path,
) -> None:
    """Run mock COLOC for each TWAS hit and save JSON results."""
    import numpy as np

    if twas_hits.empty:
        logger.info("No TWAS hits; skipping COLOC stage.")
        return

    rng = np.random.default_rng(0)
    for _, row in twas_hits.head(10).iterrows():  # limit to 10 in mock mode
        gene_id = str(row.get("gene", "GENE_UNKNOWN"))
        tissue  = str(row.get("tissue", "tissue_unknown"))

        # Generate mock posteriors (biased toward H4 for a subset)
        pp4 = float(rng.beta(5, 2)) if rng.random() > 0.5 else float(rng.beta(1, 5))
        pp3 = float(rng.beta(1, 5)) * (1.0 - pp4)
        raw_result = make_mock_coloc_result(pp4=pp4, pp3=pp3)
        annotated = interpret_coloc_result(raw_result, cfg)

        out_path = coloc_out_dir / tissue / f"{gene_id}.json"
        save_coloc_result(annotated, out_path)

    logger.info("[MOCK] Saved COLOC results for %d gene-tissue pairs.", min(len(twas_hits), 10))


def _generate_plots(
    twas_df: pd.DataFrame,
    coloc_df: pd.DataFrame,
    fig_dir: Path,
    cfg: Dict[str, Any],
    mock: bool = True,
) -> None:
    """Generate all summary figures."""
    viz_cfg = cfg.get("visualization", {})
    dpi = viz_cfg.get("dpi", 150)

    # In real mode we fail fast if required results are missing instead of plotting mock data.
    if not mock and (twas_df.empty or coloc_df.empty):
        raise RuntimeError(
            "Real-mode plotting requires non-empty TWAS and COLOC tables. "
            "Mock fallback plotting is disabled in real mode."
        )

    # Use mock data for plots in demo mode if outputs are empty
    plot_twas = twas_df if not twas_df.empty else make_mock_twas_df()
    plot_coloc = coloc_df if not coloc_df.empty else make_mock_coloc_df()

    plot_twas_qqplot(
        pvalues=plot_twas["pvalue"].dropna(),
        output_path=fig_dir / "qqplot_twas.png",
        dpi=dpi,
    )
    plot_twas_manhattan(
        twas_df=plot_twas,
        output_path=fig_dir / "manhattan_like_twas.png",
        dpi=dpi,
    )
    plot_coloc_pp4_barplot(
        coloc_df=plot_coloc,
        output_path=fig_dir / "coloc_summary.png",
        dpi=dpi,
    )
    plot_coloc_pp4_vs_pp3pp4(
        coloc_df=plot_coloc,
        output_path=fig_dir / "coloc_pp4_vs_pp3pp4.png",
        dpi=dpi,
    )


def _norm_cdf_arr(x):
    from scipy.special import ndtr
    return ndtr(x)


def _resolve_gwas_input_path(base_dir: Path, cfg: Dict[str, Any]) -> Path:
    """Resolve GWAS input file path from config, supporting absolute or relative paths."""
    raw_path = Path(cfg.get("gwas", {}).get("input_file", "data/raw/as_gwas.tsv.gz"))
    if raw_path.is_absolute():
        return raw_path
    return (base_dir / raw_path).resolve()


def _build_harmonization_qc(harmonized_gwas: pd.DataFrame) -> Dict[str, Any]:
    """Return harmonization QC metrics expected for first real-data milestone."""
    attrs_qc = harmonized_gwas.attrs.get("harmonization_qc", {})
    return {
        "input_snp_count": int(attrs_qc.get("input_snp_count", len(harmonized_gwas))),
        "lifted_snp_count": int(attrs_qc.get("lifted_snp_count", len(harmonized_gwas))),
        "failed_liftover_count": int(attrs_qc.get("failed_liftover_count", 0)),
        "palindromic_dropped_count": int(attrs_qc.get("palindromic_dropped_count", 0)),
        "final_retained_snp_count": int(attrs_qc.get("final_retained_snp_count", len(harmonized_gwas))),
    }


def _run_gtex_overlap_check(
    harmonized_gwas: pd.DataFrame,
    cfg: Dict[str, Any],
    mock: bool,
) -> Optional[Dict[str, Any]]:
    """Run overlap sanity check against one GTEx model namespace."""
    gtex_cfg = cfg.get("gtex", {})
    sp_cfg = cfg.get("spredixcan", {})
    validation_cfg = cfg.get("validation", {})

    tissue = validation_cfg.get("gtex_overlap_tissue", "Whole_Blood")
    min_fraction = float(validation_cfg.get("gtex_overlap_min_fraction", 0.0))

    model_dir = Path(gtex_cfg.get("model_dir", "data/reference/gtex_v8_models"))
    if not model_dir.is_absolute():
        model_dir = Path(cfg.get("_base_dir", ".")) / model_dir

    pattern = sp_cfg.get("model_db_pattern", "{model_dir}/gtex_v8_mashr_{tissue}.db")
    model_path = Path(pattern.format(model_dir=model_dir, tissue=tissue))

    if not model_path.exists():
        if mock:
            logger.warning("Skipping GTEx overlap check in mock mode; model DB missing: %s", model_path)
            return None
        raise FileNotFoundError(
            f"GTEx overlap check requires model DB for {tissue}: {model_path}"
        )

    overlap = check_gtex_model_varid_overlap(harmonized_gwas, model_path, tissue=tissue)
    if not mock and overlap["overlap_fraction"] < min_fraction:
        raise RuntimeError(
            f"GTEx overlap fraction {overlap['overlap_fraction']:.6f} is below configured minimum {min_fraction:.6f}"
        )
    return overlap
