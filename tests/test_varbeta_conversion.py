"""
test_varbeta_conversion.py
--------------------------
Tests for the varbeta = se**2 conversion required by coloc.abf.

This is a critical correctness test: passing se instead of varbeta to coloc
produces incorrect posterior probabilities and is a very common mistake.
"""

import numpy as np
import pandas as pd
import pytest

from src.coloc import compute_varbeta


def test_varbeta_equals_se_squared():
    """varbeta must equal se**2 for every row."""
    se_values = np.array([0.05, 0.10, 0.20, 0.03, 0.15])
    df = pd.DataFrame({"se": se_values, "beta": np.ones(len(se_values))})
    result = compute_varbeta(df)

    np.testing.assert_allclose(
        result["varbeta"].values,
        se_values ** 2,
        rtol=1e-10,
        err_msg="varbeta must equal se**2",
    )


def test_varbeta_is_strictly_positive():
    """varbeta must be positive for all valid se values."""
    df = pd.DataFrame({"se": [0.01, 0.05, 0.5, 2.0]})
    result = compute_varbeta(df)
    assert (result["varbeta"] > 0).all()


def test_varbeta_column_is_new():
    """compute_varbeta adds the varbeta column (does not overwrite input df)."""
    df = pd.DataFrame({"se": [0.1, 0.2], "beta": [0.5, -0.3]})
    result = compute_varbeta(df)
    # Original should not have varbeta
    assert "varbeta" not in df.columns
    assert "varbeta" in result.columns


def test_varbeta_missing_se_raises():
    """compute_varbeta raises ValueError if SE column is missing."""
    df = pd.DataFrame({"beta": [0.1, 0.2]})
    with pytest.raises(ValueError, match="not found in DataFrame"):
        compute_varbeta(df)


def test_varbeta_custom_se_col():
    """compute_varbeta works with a non-default SE column name."""
    df = pd.DataFrame({"standard_error": [0.05, 0.10]})
    result = compute_varbeta(df, se_col="standard_error")
    np.testing.assert_allclose(result["varbeta"].values, [0.0025, 0.01], rtol=1e-10)


def test_varbeta_not_equal_to_se():
    """Sanity check: varbeta != se (this would be the common mistake)."""
    se_values = np.array([0.05, 0.10, 0.20])
    df = pd.DataFrame({"se": se_values})
    result = compute_varbeta(df)
    # varbeta should NOT equal se (unless se happens to be 0 or 1, which it won't be here)
    assert not np.allclose(result["varbeta"].values, se_values), \
        "varbeta should not equal se — common coloc pitfall!"


def test_build_coloc_dataset_includes_varbeta():
    """build_coloc_dataset output dicts contain varbeta, not se."""
    from src.coloc import build_coloc_dataset

    n = 5
    varids = [f"chr1_{i}_A_C_b38" for i in range(n)]
    se = 0.05

    gwas_locus = pd.DataFrame(
        {
            "varID": varids,
            "beta": np.random.default_rng(0).normal(0, 0.05, n),
            "se": [se] * n,
            "pvalue": [0.05] * n,
        }
    )
    eqtl_locus = pd.DataFrame(
        {
            "varID": varids,
            "beta": np.random.default_rng(1).normal(0, 0.1, n),
            "se": [se] * n,
            "pvalue": [0.1] * n,
            "N": [200] * n,
        }
    )

    ds1, ds2 = build_coloc_dataset(
        gwas_locus=gwas_locus,
        eqtl_locus=eqtl_locus,
        n_gwas=10619,
        n_cases=4925,
        n_controls=5694,
    )

    assert "varbeta" in ds1, "dataset1 must contain varbeta key"
    assert "varbeta" in ds2, "dataset2 must contain varbeta key"

    # Confirm varbeta ≈ se**2
    np.testing.assert_allclose(
        np.array(ds1["varbeta"]),
        np.full(n, se ** 2),
        rtol=1e-8,
    )
