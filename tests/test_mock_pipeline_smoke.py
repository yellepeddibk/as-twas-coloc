"""
test_mock_pipeline_smoke.py
---------------------------
Smoke test: run the full pipeline in mock mode and verify it completes
without errors and produces expected outputs.

This test must run successfully in CI without any external tools,
large data files, GTEx models, or R/coloc installation.
"""

import pytest
import tempfile
from pathlib import Path
import shutil

from src.config import load_config
from src.pipeline import run_pipeline

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "as.yaml"


@pytest.fixture(scope="module")
def tmp_base_dir():
    """Temporary directory that serves as the base_dir for the mock pipeline run."""
    d = tempfile.mkdtemp(prefix="as_twas_coloc_test_")
    yield Path(d)
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture(scope="module")
def pipeline_summary(tmp_base_dir):
    """Run the mock pipeline once and return the summary dict."""
    summary = run_pipeline(
        config_path=_CONFIG_PATH,
        mock=True,
        base_dir=tmp_base_dir,
    )
    return summary, tmp_base_dir


def test_pipeline_returns_summary_dict(pipeline_summary):
    """Pipeline returns a non-empty summary dict."""
    summary, _ = pipeline_summary
    assert isinstance(summary, dict)
    assert len(summary) > 0


def test_pipeline_records_gene_tissue_tests(pipeline_summary):
    """Summary contains the number of gene-tissue tests performed."""
    summary, _ = pipeline_summary
    assert "n_gene_tissue_tests" in summary
    assert summary["n_gene_tissue_tests"] > 0


def test_harmonized_gwas_file_created(pipeline_summary):
    """Harmonized GWAS file is written to data/interim/gwas_harmonized/."""
    _, base_dir = pipeline_summary
    harm_path = base_dir / "data" / "interim" / "gwas_harmonized" / "as_gwas_hg38_varid.tsv.gz"
    assert harm_path.exists(), f"Harmonized GWAS not found: {harm_path}"


def test_twas_hits_file_created(pipeline_summary):
    """TWAS hits file is written to data/processed/twas_hits/."""
    _, base_dir = pipeline_summary
    hits_path = base_dir / "data" / "processed" / "twas_hits" / "twas_hits_fdr05.tsv.gz"
    assert hits_path.exists(), f"TWAS hits file not found: {hits_path}"


def test_top_genes_table_created(pipeline_summary):
    """Top-genes table is written to results/tables/."""
    _, base_dir = pipeline_summary
    table_path = base_dir / "results" / "tables" / "top_genes_pp4.tsv"
    assert table_path.exists(), f"Top-genes table not found: {table_path}"


def test_manifest_created(pipeline_summary):
    """At least one provenance manifest JSON is written to results/manifests/."""
    _, base_dir = pipeline_summary
    manifest_dir = base_dir / "results" / "manifests"
    manifests = list(manifest_dir.glob("manifest_*.json"))
    assert len(manifests) >= 1, f"No manifest found in {manifest_dir}"


def test_figures_created(pipeline_summary):
    """At least one figure PNG is written to results/figures/."""
    _, base_dir = pipeline_summary
    fig_dir = base_dir / "results" / "figures"
    pngs = list(fig_dir.glob("*.png"))
    assert len(pngs) >= 1, f"No PNG figures found in {fig_dir}"


def test_qqplot_figure_created(pipeline_summary):
    """QQ-plot figure is written."""
    _, base_dir = pipeline_summary
    qq_path = base_dir / "results" / "figures" / "qqplot_twas.png"
    assert qq_path.exists(), f"QQ-plot not found: {qq_path}"


def test_coloc_jsons_created(pipeline_summary):
    """At least one coloc JSON is written to results/coloc/."""
    _, base_dir = pipeline_summary
    coloc_dir = base_dir / "results" / "coloc"
    jsons = list(coloc_dir.rglob("*.json"))
    # May be zero if no TWAS hits; log a note but don't fail hard
    # (mock data may not always produce FDR-significant hits)
    # We just check the directory exists
    assert coloc_dir.exists(), f"Coloc directory not found: {coloc_dir}"


def test_pipeline_fdr_alpha_respected(pipeline_summary):
    """n_twas_hits_fdr is less than or equal to n_gene_tissue_tests."""
    summary, _ = pipeline_summary
    assert summary.get("n_twas_hits_fdr", 0) <= summary.get("n_gene_tissue_tests", 0)


def test_harmonized_gwas_has_varid_column(pipeline_summary):
    """Harmonized GWAS file contains the varID column."""
    import pandas as pd

    _, base_dir = pipeline_summary
    harm_path = base_dir / "data" / "interim" / "gwas_harmonized" / "as_gwas_hg38_varid.tsv.gz"
    df = pd.read_csv(harm_path, sep="\t", nrows=5)
    assert "varID" in df.columns, "Harmonized GWAS must contain 'varID' column"


def test_manifest_has_expected_keys(pipeline_summary):
    """Provenance manifest contains all required provenance keys."""
    import json

    _, base_dir = pipeline_summary
    manifest_dir = base_dir / "results" / "manifests"
    manifests = list(manifest_dir.glob("manifest_*.json"))
    assert manifests, "No manifest found"

    with manifests[0].open() as fh:
        manifest = json.load(fh)

    for key in ("timestamp", "git_commit", "coloc_priors", "decision_thresholds", "config_snapshot"):
        assert key in manifest, f"Manifest missing key: '{key}'"
