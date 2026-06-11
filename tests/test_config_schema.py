"""Schema tests for the pipeline YAML config.

Guards against silent structural drift in ``config/as.yaml`` so downstream
loaders, the harmonization/TWAS/COLOC stages, and the prereq checker can rely
on the same shape across runs. These complement ``test_config_loading.py`` by
asserting per-section invariants and cross-field consistency.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.config import load_config

_REPO_ROOT = Path(__file__).resolve().parent.parent
_CONFIG_PATH = _REPO_ROOT / "config" / "as.yaml"


@pytest.fixture(scope="module")
def cfg() -> dict:
    assert _CONFIG_PATH.exists(), f"Missing config file: {_CONFIG_PATH}"
    return load_config(_CONFIG_PATH)


def test_project_section_well_formed(cfg):
    project = cfg["project"]
    for key in ("name", "disease", "gwas_trait", "genome_build"):
        assert key in project, f"Missing project key: {key}"
        assert isinstance(project[key], str) and project[key], f"{key} must be non-empty string"
    assert project["genome_build"] == "hg38", \
        "Pipeline targets hg38; update harmonization if this changes"


def test_paths_section_has_required_keys(cfg):
    paths = cfg["paths"]
    required_dirs = (
        "data_raw",
        "data_interim",
        "data_processed",
        "data_reference",
        "results_twas",
        "results_coloc",
        "results_tables",
        "results_figures",
        "results_manifests",
    )
    for key in required_dirs:
        assert key in paths, f"Missing paths key: {key}"
        assert isinstance(paths[key], str) and paths[key], f"{key} must be non-empty string"

    required_outputs = ("gwas_harmonized", "twas_hits", "coloc_summary")
    for key in required_outputs:
        assert key in paths, f"Missing paths output key: {key}"
        assert isinstance(paths[key], str) and paths[key], f"{key} must be non-empty string"


def test_gwas_section_has_required_keys(cfg):
    gwas = cfg["gwas"]
    for key in ("input_file", "source_build", "n_total", "n_cases", "n_controls", "case_fraction"):
        assert key in gwas, f"Missing GWAS key: {key}"

    assert isinstance(gwas["input_file"], str) and gwas["input_file"]
    assert gwas["source_build"] in ("hg18", "hg19", "hg38"), \
        f"Unsupported source_build: {gwas['source_build']}"
    for key in ("n_total", "n_cases", "n_controls"):
        assert isinstance(gwas[key], int) and gwas[key] > 0, f"{key} must be positive int"


def test_gwas_sample_size_is_consistent(cfg):
    gwas = cfg["gwas"]
    assert gwas["n_cases"] + gwas["n_controls"] == gwas["n_total"], \
        "n_cases + n_controls must equal n_total"
    expected_fraction = gwas["n_cases"] / gwas["n_total"]
    actual_fraction = float(gwas["case_fraction"])
    assert abs(actual_fraction - expected_fraction) < 1e-3, (
        f"case_fraction ({actual_fraction}) does not match "
        f"n_cases / n_total ({expected_fraction:.4f})"
    )
    assert 0.0 < actual_fraction < 1.0, "case_fraction must be in (0, 1)"


def test_gwas_harmonization_flags_are_bool(cfg):
    gwas = cfg["gwas"]
    for key in ("drop_palindromic_snps", "af_disambiguation_available"):
        assert key in gwas, f"Missing harmonization flag: {key}"
        assert isinstance(gwas[key], bool), f"{key} must be bool"


def test_liftover_section_well_formed(cfg):
    liftover = cfg["liftover"]
    assert isinstance(liftover.get("enabled"), bool), "liftover.enabled must be bool"
    if liftover["enabled"]:
        chain = liftover.get("chain_file")
        assert isinstance(chain, str) and chain, \
            "liftover.chain_file required when liftover.enabled is true"
        tool = liftover.get("tool")
        assert tool in ("pyliftover", "crossmap"), \
            f"Unsupported liftover tool: {tool}"


def test_validation_section_well_formed(cfg):
    validation = cfg["validation"]
    assert isinstance(validation["gtex_overlap_tissue"], str) and validation["gtex_overlap_tissue"]
    fraction = float(validation["gtex_overlap_min_fraction"])
    assert 0.0 < fraction <= 1.0, "gtex_overlap_min_fraction must be in (0, 1]"


def test_gtex_section_well_formed(cfg):
    gtex = cfg["gtex"]
    for key in ("model_dir", "eqtl_dir"):
        assert key in gtex, f"Missing GTEx key: {key}"
        assert isinstance(gtex[key], str) and gtex[key]

    tissues = gtex["tissues"]
    assert isinstance(tissues, list) and len(tissues) > 0, "GTEx tissues must be a non-empty list"
    assert all(isinstance(t, str) and t for t in tissues), \
        "Each tissue must be a non-empty string"
    assert len(tissues) == len(set(tissues)), "GTEx tissues must be unique"


def test_spredixcan_section_has_column_mappings(cfg):
    spredixcan = cfg["spredixcan"]
    for key in ("script", "model_db_pattern", "covariance_pattern", "gwas_file", "output_dir"):
        assert key in spredixcan, f"Missing S-PrediXcan key: {key}"
        assert isinstance(spredixcan[key], str) and spredixcan[key]

    column_keys = (
        "snp_column",
        "model_db_snp_key",
        "effect_allele_column",
        "non_effect_allele_column",
        "beta_column",
        "se_column",
        "pvalue_column",
    )
    for key in column_keys:
        assert key in spredixcan, f"Missing S-PrediXcan column mapping: {key}"
        assert isinstance(spredixcan[key], str) and spredixcan[key]

    for pattern_key in ("model_db_pattern", "covariance_pattern"):
        pattern = spredixcan[pattern_key]
        assert "{tissue}" in pattern, f"{pattern_key} must contain '{{tissue}}' placeholder"


def test_correction_section_well_formed(cfg):
    correction = cfg["correction"]
    assert correction["method"] in ("fdr_bh", "fdr_by", "bonferroni", "holm"), \
        f"Unsupported correction method: {correction['method']}"
    alpha = float(correction["alpha"])
    assert 0.0 < alpha < 1.0, "correction.alpha must be in (0, 1)"


def test_coloc_priors_are_floats(cfg):
    coloc = cfg["coloc"]
    for key in ("p1", "p2", "p12"):
        assert key in coloc, f"Missing coloc prior: {key}"
        value = float(coloc[key])
        assert 0.0 < value < 1.0, f"coloc.{key} must be in (0, 1)"


def test_coloc_window_and_thresholds(cfg):
    coloc = cfg["coloc"]
    window = coloc["window_bp"]
    assert isinstance(window, int) and window > 0, "coloc.window_bp must be positive int"

    for key in ("pp4_threshold", "pp3_pp4_threshold", "pp4_ratio_threshold"):
        value = float(coloc[key])
        assert 0.0 < value <= 1.0, f"coloc.{key} must be in (0, 1]"

    assert isinstance(coloc["coloc_susie_available"], bool), \
        "coloc.coloc_susie_available must be bool"
    assert isinstance(coloc["r_script"], str) and coloc["r_script"], \
        "coloc.r_script must be a non-empty string"
    assert isinstance(coloc["output_dir"], str) and coloc["output_dir"], \
        "coloc.output_dir must be a non-empty string"


def test_visualization_section_well_formed(cfg):
    viz = cfg["visualization"]
    dpi = viz["dpi"]
    assert isinstance(dpi, int) and dpi > 0, "visualization.dpi must be positive int"
    assert viz["figure_format"] in ("png", "pdf", "svg"), \
        f"Unsupported figure_format: {viz['figure_format']}"


def test_provenance_section_well_formed(cfg):
    provenance = cfg["provenance"]
    assert isinstance(provenance["manifest_dir"], str) and provenance["manifest_dir"], \
        "provenance.manifest_dir must be a non-empty string"
    assert isinstance(provenance["record_checksums"], bool), \
        "provenance.record_checksums must be bool"
