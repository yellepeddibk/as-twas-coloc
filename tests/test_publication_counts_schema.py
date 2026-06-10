"""Schema tests for the frozen publication counts JSON.

Guards against silent format drift in ``config/publication_counts.json`` so the
validator script can always parse it the same way.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_COUNTS_PATH = _REPO_ROOT / "config" / "publication_counts.json"


@pytest.fixture(scope="module")
def counts() -> dict:
    assert _COUNTS_PATH.exists(), f"Missing publication counts file: {_COUNTS_PATH}"
    with _COUNTS_PATH.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def test_top_level_sections_present(counts):
    for key in ("gwas_harmonization", "twas", "coloc"):
        assert key in counts, f"Missing top-level section: {key}"
        assert isinstance(counts[key], dict), f"Section must be a dict: {key}"


def test_gwas_harmonization_counts_are_non_negative_ints(counts):
    section = counts["gwas_harmonization"]
    for key in (
        "input_snp_count",
        "lifted_snp_count",
        "palindromic_dropped_count",
        "final_retained_snp_count",
    ):
        assert key in section, f"Missing harmonization key: {key}"
        value = section[key]
        assert isinstance(value, int), f"{key} must be int, got {type(value).__name__}"
        assert value >= 0, f"{key} must be non-negative, got {value}"


def test_gwas_funnel_is_monotonic(counts):
    section = counts["gwas_harmonization"]
    assert section["lifted_snp_count"] <= section["input_snp_count"], \
        "lifted_snp_count cannot exceed input_snp_count"
    assert section["final_retained_snp_count"] <= section["lifted_snp_count"], \
        "final_retained_snp_count cannot exceed lifted_snp_count"


def test_twas_counts_have_expected_types(counts):
    section = counts["twas"]
    for key in ("raw_gene_tissue_rows", "valid_gene_tissue_tests", "fdr_significant_hits"):
        assert key in section, f"Missing TWAS key: {key}"
        assert isinstance(section[key], int) and section[key] >= 0, \
            f"{key} must be non-negative int"

    assert "tissues" in section, "Missing TWAS tissues list"
    assert isinstance(section["tissues"], list)
    assert all(isinstance(t, str) and t for t in section["tissues"]), \
        "tissues must be a list of non-empty strings"
    assert len(section["tissues"]) == len(set(section["tissues"])), \
        "tissues list must be unique"


def test_twas_counts_are_internally_consistent(counts):
    section = counts["twas"]
    assert section["valid_gene_tissue_tests"] <= section["raw_gene_tissue_rows"], \
        "valid tests cannot exceed raw rows"
    assert section["fdr_significant_hits"] <= section["valid_gene_tissue_tests"], \
        "FDR hits cannot exceed valid tests"


def test_coloc_section_well_formed(counts):
    section = counts["coloc"]
    assert "pp4_threshold" in section, "Missing coloc pp4_threshold"
    threshold = section["pp4_threshold"]
    assert isinstance(threshold, (int, float)) and 0.0 <= float(threshold) <= 1.0, \
        "pp4_threshold must be in [0, 1]"

    for key in ("coloc_confirmed_pairs", "unique_genes", "whole_blood_pairs"):
        assert key in section, f"Missing coloc key: {key}"
        assert isinstance(section[key], int) and section[key] >= 0, \
            f"{key} must be non-negative int"

    assert section["unique_genes"] <= section["coloc_confirmed_pairs"], \
        "unique_genes cannot exceed coloc_confirmed_pairs"
    assert section["whole_blood_pairs"] <= section["coloc_confirmed_pairs"], \
        "whole_blood_pairs cannot exceed coloc_confirmed_pairs"


def test_validator_can_parse_counts_file():
    """The validator script must be importable and able to load the JSON."""
    from scripts.validate_publication_counts import _load_counts

    parsed = _load_counts(_COUNTS_PATH)
    assert isinstance(parsed, dict)
    assert {"gwas_harmonization", "twas", "coloc"}.issubset(parsed.keys())
