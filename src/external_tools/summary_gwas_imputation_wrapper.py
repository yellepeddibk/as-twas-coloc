"""
summary_gwas_imputation_wrapper.py
----------------------------------
Placeholder wrapper for summary-gwas-imputation (MetaXcan preprocessing tool).

summary-gwas-imputation (https://github.com/hakyimlab/summary-gwas-imputation)
imputes GWAS summary statistics to GTEx variant positions and constructs
GTEx-style varIDs.  It is a required preprocessing step when the GWAS panel
has incomplete overlap with the GTEx eQTL models.

Real execution requires:
  1. Clone: https://github.com/hakyimlab/summary-gwas-imputation
  2. Reference panel: GTEx v8 1000G reference (available from the MetaXcan team)
  3. GWAS summary stats in the expected input format

This wrapper:
  - Provides a command-builder skeleton
  - Returns a clearly-labeled placeholder result in mock mode
  - Does NOT perform any real computation

TODO: Implement real summary-gwas-imputation support once reference panel
      data are acquired and the imputation step is validated.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from src.logging_utils import get_logger

logger = get_logger(__name__)


def build_gwas_imputation_command(
    gwas_file: str | Path,
    output_dir: str | Path,
    parquet_genotype: str | Path,
    parquet_genotype_metadata: str | Path,
    chromosome: int,
    script_dir: str | Path = "external/summary-gwas-imputation/src",
) -> List[str]:
    """Build the command for summary-gwas-imputation.

    Parameters
    ----------
    gwas_file:
        Input GWAS summary statistics (harmonized, hg38).
    output_dir:
        Directory for imputed output files.
    parquet_genotype:
        Path to the reference genotype in Parquet format.
    parquet_genotype_metadata:
        Path to variant metadata Parquet file.
    chromosome:
        Chromosome number (1-22) to impute.
    script_dir:
        Directory containing the summary-gwas-imputation Python scripts.

    Returns
    -------
    List[str]
        Command list for subprocess.
    """
    script = Path(script_dir) / "gwas_summary_imputation.py"

    logger.warning(
        "[PLACEHOLDER] summary-gwas-imputation command builder.  "
        "Script path: %s  (may not exist without real installation).",
        script,
    )

    cmd = [
        "python", str(script),
        "--gwas_file",                    str(gwas_file),
        "--output_folder",                str(output_dir),
        "--parquet_genotype",             str(parquet_genotype),
        "--parquet_genotype_metadata",    str(parquet_genotype_metadata),
        "--chromosome",                   str(chromosome),
        "--regularization",               "0.1",
        "--frequency_filter",             "0.01",
        "--sub_batches",                  "10",
        "--sub_batch",                    "0",
    ]
    return cmd


def run_gwas_imputation_placeholder(
    gwas_file: str | Path,
    output_dir: str | Path,
) -> Dict[str, Any]:
    """Return a clearly-labeled placeholder result for GWAS imputation.

    Parameters
    ----------
    gwas_file:
        Input GWAS file path.
    output_dir:
        Intended output directory.

    Returns
    -------
    Dict[str, Any]
        Placeholder result dict.
    """
    logger.warning(
        "[PLACEHOLDER] summary-gwas-imputation not run.  "
        "Output would go to %s.",
        output_dir,
    )
    return {
        "status":     "placeholder",
        "method":     "summary-gwas-imputation",
        "gwas_file":  str(gwas_file),
        "output_dir": str(output_dir),
        "note": (
            "summary-gwas-imputation requires the GTEx v8 reference panel "
            "(not included in this repository).  "
            "See README for data acquisition instructions."
        ),
    }
