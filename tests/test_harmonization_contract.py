"""
test_harmonization_contract.py
--------------------------------
Tests for GWAS harmonization logic and schema validation.
"""

import pytest
import pandas as pd
import numpy as np

from src.schemas import (
    HARMONIZED_GWAS_REQUIRED,
    apply_column_aliases,
    validate_harmonized_schema,
)
from src.harmonization import (
    harmonize_gwas,
    _is_palindromic,
    _strand_flip,
    _drop_palindromic,
    _build_varid,
    _coerce_types,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────


def make_minimal_raw_gwas(n: int = 10) -> pd.DataFrame:
    """Minimal raw GWAS DataFrame with standard column names."""
    rng = np.random.default_rng(42)
    betas = rng.normal(0, 0.05, n)
    ses = rng.uniform(0.01, 0.05, n)
    pvals = 2 * (1.0 - 1.0 / (1.0 + np.exp(-np.abs(betas / ses))))
    # Cycle through non-palindromic allele pairs to fill any n
    _a1_pool = ["A", "C", "G", "T"]
    _a2_pool = ["C", "A", "T", "C"]  # none of these form AT/CG pairs
    a1 = [_a1_pool[i % len(_a1_pool)] for i in range(n)]
    a2 = [_a2_pool[i % len(_a2_pool)] for i in range(n)]
    return pd.DataFrame(
        {
            "CHR":  [str(i % 22 + 1) for i in range(n)],
            "BP":   rng.integers(1_000_000, 200_000_000, n),
            "A1":   a1,
            "A2":   a2,
            "BETA": betas,
            "SE":   ses,
            "P":    pvals,
            "N":    [10619] * n,
        }
    )


def make_minimal_config() -> dict:
    return {
        "gwas": {
            "source_build": "hg38",
            "drop_palindromic_snps": True,
            "af_disambiguation_available": False,
            "n_total": 10619,
            "n_cases": 4925,
            "n_controls": 5694,
        },
        "paths": {},
        "coloc": {},
    }


# ── Column alias tests ────────────────────────────────────────────────────────


def test_apply_column_aliases_basic():
    """Standard column names are mapped to canonical names."""
    df = pd.DataFrame({"CHR": ["1"], "BP": [100], "A1": ["A"], "A2": ["C"], "BETA": [0.1], "SE": [0.01], "P": [0.05]})
    result = apply_column_aliases(df)
    assert "chrom" in result.columns
    assert "pos" in result.columns
    assert "effect_allele" in result.columns
    assert "non_effect_allele" in result.columns
    assert "beta" in result.columns
    assert "se" in result.columns
    assert "pvalue" in result.columns


def test_apply_column_aliases_no_overwrite():
    """Alias mapping does not overwrite an already-canonical column."""
    df = pd.DataFrame({"chrom": ["1"], "CHR": ["2"]})
    result = apply_column_aliases(df)
    # canonical column should be unchanged
    assert "chrom" in result.columns
    assert result["chrom"].iloc[0] == "1"


# ── Palindromic SNP tests ─────────────────────────────────────────────────────


@pytest.mark.parametrize("a1,a2,expected", [
    ("A", "T", True),
    ("T", "A", True),
    ("C", "G", True),
    ("G", "C", True),
    ("A", "C", False),
    ("G", "T", False),
    ("A", "G", False),
])
def test_is_palindromic(a1, a2, expected):
    assert _is_palindromic(a1, a2) == expected


def test_drop_palindromic_removes_at_cg_snps():
    """_drop_palindromic removes AT and CG SNPs."""
    df = pd.DataFrame(
        {
            "effect_allele":     ["A", "C", "G", "A"],
            "non_effect_allele": ["T", "G", "C", "C"],  # first 3 palindromic, last not
        }
    )
    result = _drop_palindromic(df)
    assert len(result) == 1
    assert result.iloc[0]["effect_allele"] == "A"
    assert result.iloc[0]["non_effect_allele"] == "C"


def test_drop_palindromic_preserves_non_palindromic():
    """Non-palindromic SNPs are all kept."""
    df = pd.DataFrame(
        {
            "effect_allele":     ["A", "C", "G", "T"],
            "non_effect_allele": ["C", "T", "A", "G"],
        }
    )
    result = _drop_palindromic(df)
    assert len(result) == 4


# ── Strand flip ───────────────────────────────────────────────────────────────


@pytest.mark.parametrize("allele,complement", [
    ("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"),
])
def test_strand_flip(allele, complement):
    assert _strand_flip(allele) == complement


# ── varID construction ────────────────────────────────────────────────────────


def test_build_varid_format():
    """varID must follow chr{chrom}_{pos}_{ref}_{alt}_b38 format."""
    df = pd.DataFrame(
        {
            "chrom": ["1", "6"],
            "pos": [1000000, 31000000],
            "non_effect_allele": ["A", "C"],
            "effect_allele": ["G", "T"],
        }
    )
    result = _build_varid(df)
    assert result.loc[0, "varID"] == "chr1_1000000_A_G_b38"
    assert result.loc[1, "varID"] == "chr6_31000000_C_T_b38"


def test_build_varid_idempotent():
    """_build_varid does not overwrite an existing varID column."""
    df = pd.DataFrame({"varID": ["existing_id"], "chrom": ["1"], "pos": [100],
                       "non_effect_allele": ["A"], "effect_allele": ["C"]})
    result = _build_varid(df)
    assert result.loc[0, "varID"] == "existing_id"


# ── Schema validation ─────────────────────────────────────────────────────────


def test_validate_harmonized_schema_passes():
    """A fully populated harmonized DataFrame passes schema validation."""
    df = pd.DataFrame(
        {c: ["x"] for c in HARMONIZED_GWAS_REQUIRED}
    )
    # Should not raise
    validate_harmonized_schema(df)


def test_validate_harmonized_schema_missing_column():
    """Missing required columns raise ValueError."""
    df = pd.DataFrame({"varID": ["x"], "chrom": ["1"]})
    with pytest.raises(ValueError, match="missing required columns"):
        validate_harmonized_schema(df)


# ── End-to-end harmonize_gwas ─────────────────────────────────────────────────


def test_harmonize_gwas_e2e():
    """harmonize_gwas on minimal toy data produces a valid harmonized DataFrame."""
    raw = make_minimal_raw_gwas(20)
    cfg = make_minimal_config()
    harmonized = harmonize_gwas(raw, cfg)

    # Check schema
    validate_harmonized_schema(harmonized)

    # varIDs should be non-empty and correctly formatted
    assert harmonized["varID"].str.startswith("chr").all()
    assert harmonized["varID"].str.endswith("_b38").all()

    # No palindromic SNPs should remain (drop_palindromic_snps=True)
    from src.harmonization import _is_palindromic
    palindromic_mask = harmonized.apply(
        lambda r: _is_palindromic(r["effect_allele"], r["non_effect_allele"]), axis=1
    )
    assert palindromic_mask.sum() == 0, "Palindromic SNPs were not removed"


def test_harmonize_gwas_hg19_warning(caplog):
    """harmonize_gwas logs a warning when source_build is hg19."""
    import logging
    raw = make_minimal_raw_gwas(5)
    cfg = make_minimal_config()
    cfg["gwas"]["source_build"] = "hg19"

    with caplog.at_level(logging.WARNING, logger="src.harmonization"):
        harmonize_gwas(raw, cfg)

    assert any("liftover" in msg.lower() or "hg19" in msg.lower() for msg in caplog.messages)
