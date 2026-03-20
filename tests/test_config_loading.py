"""
test_config_loading.py
----------------------
Tests for YAML configuration loading and access utilities.
"""

import pytest
from pathlib import Path

from src.config import load_config, get_coloc_priors, get_tissues


_CONFIG_PATH = Path(__file__).parent.parent / "config" / "as.yaml"


def test_load_config_returns_dict():
    """load_config returns a non-empty dictionary."""
    cfg = load_config(_CONFIG_PATH)
    assert isinstance(cfg, dict)
    assert len(cfg) > 0


def test_load_config_has_required_top_level_keys():
    """Config contains expected top-level sections."""
    cfg = load_config(_CONFIG_PATH)
    for key in ("project", "paths", "gwas", "gtex", "coloc", "correction"):
        assert key in cfg, f"Missing top-level config key: '{key}'"


def test_load_config_missing_file_raises():
    """load_config raises FileNotFoundError for a non-existent path."""
    with pytest.raises(FileNotFoundError):
        load_config("/nonexistent/path/config.yaml")


def test_get_coloc_priors_returns_expected_keys():
    """get_coloc_priors returns a dict with p1, p2, p12."""
    cfg = load_config(_CONFIG_PATH)
    priors = get_coloc_priors(cfg)
    assert set(priors.keys()) == {"p1", "p2", "p12"}


def test_get_coloc_priors_values_are_positive():
    """Coloc priors must be positive probabilities."""
    cfg = load_config(_CONFIG_PATH)
    priors = get_coloc_priors(cfg)
    for key, val in priors.items():
        assert val > 0, f"Prior {key} must be > 0"
        assert val < 1, f"Prior {key} must be < 1"


def test_get_coloc_priors_p12_less_than_p1_p2():
    """p12 should be smaller than p1 and p2 (standard coloc assumption)."""
    cfg = load_config(_CONFIG_PATH)
    priors = get_coloc_priors(cfg)
    assert priors["p12"] < priors["p1"], "p12 should be < p1"
    assert priors["p12"] < priors["p2"], "p12 should be < p2"


def test_get_tissues_returns_list():
    """get_tissues returns a non-empty list."""
    cfg = load_config(_CONFIG_PATH)
    tissues = get_tissues(cfg)
    assert isinstance(tissues, list)
    assert len(tissues) > 0


def test_get_tissues_contains_whole_blood():
    """Whole_Blood should be in the default tissue list."""
    cfg = load_config(_CONFIG_PATH)
    tissues = get_tissues(cfg)
    assert "Whole_Blood" in tissues


def test_gwas_config_has_n_total():
    """GWAS config must specify n_total for coloc case fraction computation."""
    cfg = load_config(_CONFIG_PATH)
    assert "n_total" in cfg["gwas"], "n_total required for case fraction"
    assert cfg["gwas"]["n_total"] > 0


def test_correction_config_has_fdr_method():
    """Correction section must specify FDR method."""
    cfg = load_config(_CONFIG_PATH)
    assert "method" in cfg["correction"]
    assert cfg["correction"]["method"] in ("fdr_bh", "fdr_by", "bonferroni", "holm")


def test_coloc_thresholds_reasonable():
    """Coloc thresholds must be in [0, 1] range."""
    cfg = load_config(_CONFIG_PATH)
    coloc = cfg["coloc"]
    for key in ("pp4_threshold", "pp3_pp4_threshold", "pp4_ratio_threshold"):
        val = coloc[key]
        assert 0.0 < val <= 1.0, f"{key} = {val} is out of [0,1] range"
