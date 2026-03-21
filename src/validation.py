"""
validation.py
-------------
Validation helpers for early real-data checks.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any, Dict

import pandas as pd

from src.logging_utils import get_logger

logger = get_logger(__name__)


def check_gtex_model_varid_overlap(
    harmonized_gwas: pd.DataFrame,
    model_db_path: str | Path,
    tissue: str,
) -> Dict[str, Any]:
    """Compute varID overlap between harmonized GWAS and one GTEx model DB."""
    model_db_path = Path(model_db_path)
    if not model_db_path.exists():
        raise FileNotFoundError(f"GTEx model DB not found: {model_db_path}")

    if "varID" not in harmonized_gwas.columns:
        raise ValueError("Harmonized GWAS is missing varID column.")

    with sqlite3.connect(model_db_path) as conn:
        table_info = conn.execute("PRAGMA table_info(weights)").fetchall()
        if not table_info:
            raise RuntimeError(f"weights table not found in model DB: {model_db_path}")

        available_cols = {row[1] for row in table_info}
        candidate_cols = ["varID", "varid", "rsid"]
        variant_col = next((c for c in candidate_cols if c in available_cols), None)
        if variant_col is None:
            raise RuntimeError(
                f"No recognized variant column found in weights table. Available columns: {sorted(available_cols)}"
            )

        rows = conn.execute(f"SELECT DISTINCT {variant_col} FROM weights WHERE {variant_col} IS NOT NULL").fetchall()

    model_varids = {str(r[0]) for r in rows if r and r[0] is not None}
    gwas_varids = set(harmonized_gwas["varID"].astype(str).dropna())
    overlap = gwas_varids & model_varids

    overlap_fraction = (len(overlap) / len(gwas_varids)) if gwas_varids else 0.0
    result = {
        "tissue": tissue,
        "model_db": str(model_db_path),
        "model_variant_column": variant_col,
        "gwas_unique_varids": len(gwas_varids),
        "model_unique_varids": len(model_varids),
        "overlap_varids": len(overlap),
        "overlap_fraction": overlap_fraction,
    }

    logger.info(
        "GTEx overlap (%s): overlap=%d gwas=%d fraction=%.6f",
        tissue,
        result["overlap_varids"],
        result["gwas_unique_varids"],
        result["overlap_fraction"],
    )
    return result
