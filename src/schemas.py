"""
schemas.py
----------
Canonical GWAS data contract and schema validation for the pipeline.

The harmonized GWAS table must contain the columns defined in
``HARMONIZED_GWAS_REQUIRED`` before any downstream analysis begins.
"""

from __future__ import annotations

from typing import List

import pandas as pd

from src.logging_utils import get_logger

logger = get_logger(__name__)

# ── Required columns in the harmonized GWAS table ────────────────────────────

HARMONIZED_GWAS_REQUIRED: List[str] = [
    "varID",          # GTEx-style variant ID: chr_pos_ref_alt_b38
    "chrom",          # chromosome (no "chr" prefix, e.g. "1")
    "pos",            # hg38 position (integer)
    "effect_allele",  # allele whose effect is reported in `beta`
    "non_effect_allele",
    "beta",           # effect size (log-OR or continuous)
    "se",             # standard error of beta
    "pvalue",         # two-sided p-value
    "N",              # total sample size
]

# Optional columns that are used when present
HARMONIZED_GWAS_OPTIONAL: List[str] = [
    "rsid",
    "zscore",
    "eaf",   # effect-allele frequency
    "maf",
    "N_cases",
    "N_controls",
]

# ── Raw GWAS column aliases ───────────────────────────────────────────────────
# Map common column names found in published summary stats to the canonical schema.
RAW_COLUMN_ALIASES: dict[str, str] = {
    # rsID synonyms
    "SNP": "rsid",
    "snp": "rsid",
    "MarkerName": "rsid",
    "markername": "rsid",
    "ID": "rsid",
    # chromosome
    "CHR": "chrom",
    "chr": "chrom",
    "chromosome": "chrom",
    # position
    "BP": "pos",
    "bp": "pos",
    "POS": "pos",
    "position": "pos",
    # alleles
    "A1": "effect_allele",
    "a1": "effect_allele",
    "ALT": "effect_allele",
    "A2": "non_effect_allele",
    "a2": "non_effect_allele",
    "REF": "non_effect_allele",
    # effect
    "BETA": "beta",
    "b": "beta",
    "Effect": "beta",
    "OR": "or",         # will be log-transformed later
    # se
    "SE": "se",
    "StdErr": "se",
    # p-value
    "P": "pvalue",
    "p": "pvalue",
    "P_BOLT_LMM_INF": "pvalue",
    "P-value": "pvalue",
    "Pvalue": "pvalue",
    # sample size
    "Neff": "N",
    "n": "N",
    # frequencies
    "FRQ_A": "eaf",
    "EAF": "eaf",
    "MAF": "maf",
}


def validate_harmonized_schema(df: pd.DataFrame) -> None:
    """Raise ValueError if any required column is missing from *df*.

    Parameters
    ----------
    df:
        DataFrame representing the harmonized GWAS table.

    Raises
    ------
    ValueError
        Lists all missing required columns.
    """
    missing = [c for c in HARMONIZED_GWAS_REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(
            f"Harmonized GWAS table is missing required columns: {missing}. "
            f"Present columns: {list(df.columns)}"
        )
    logger.debug("Schema validation passed; %d rows, %d columns", len(df), len(df.columns))


def apply_column_aliases(df: pd.DataFrame) -> pd.DataFrame:
    """Rename raw GWAS columns to canonical names using ``RAW_COLUMN_ALIASES``.

    Only renames columns that are present and not already canonical.

    Parameters
    ----------
    df:
        Raw GWAS DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with renamed columns.
    """
    rename_map = {
        raw: canonical
        for raw, canonical in RAW_COLUMN_ALIASES.items()
        if raw in df.columns and canonical not in df.columns
    }
    if rename_map:
        logger.info("Renaming columns: %s", rename_map)
        df = df.rename(columns=rename_map)
    return df
