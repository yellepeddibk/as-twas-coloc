"""
test_gtex_overlap_validation.py
--------------------------------
Tests for GTEx model namespace overlap sanity checks.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd

from src.validation import check_gtex_model_varid_overlap


def test_check_gtex_model_varid_overlap_counts(tmp_path: Path):
    db_path = tmp_path / "model.db"
    with sqlite3.connect(db_path) as conn:
        conn.execute("CREATE TABLE weights (gene TEXT, rsid TEXT)")
        conn.executemany(
            "INSERT INTO weights (gene, rsid) VALUES (?, ?)",
            [
                ("G1", "chr1_100_A_G_b38"),
                ("G2", "chr1_200_C_T_b38"),
                ("G3", "chr1_300_G_A_b38"),
            ],
        )
        conn.commit()

    gwas = pd.DataFrame(
        {
            "varID": [
                "chr1_100_A_G_b38",
                "chr1_200_C_T_b38",
                "chr2_400_A_C_b38",
            ]
        }
    )

    result = check_gtex_model_varid_overlap(gwas, db_path, tissue="Whole_Blood")
    assert result["model_variant_column"] == "rsid"
    assert result["gwas_unique_varids"] == 3
    assert result["model_unique_varids"] == 3
    assert result["overlap_varids"] == 2
    assert abs(result["overlap_fraction"] - (2 / 3)) < 1e-9
