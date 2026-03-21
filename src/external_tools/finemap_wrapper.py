"""
finemap_wrapper.py
------------------
Placeholder wrapper for FINEMAP fine-mapping tool.

FINEMAP (Benner et al. 2016) performs Bayesian fine-mapping of GWAS loci
to identify credible sets of causal variants.  It can be used to escalate
analysis for top coloc-confirmed loci.

Real execution requires:
  1. FINEMAP binary (http://www.christianbenner.com)
  2. In-sample LD matrices for each locus
  3. GWAS summary statistics in FINEMAP format
  4. Zfile, bfile, and ld files for each locus

This placeholder:
  - Documents the expected interface
  - Provides a command-builder skeleton
  - Returns clearly-labeled placeholder output

TODO: Implement real FINEMAP support after coloc identifies top loci.
      Consider using SuSiE as an alternative (susie_wrapper.py).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from src.logging_utils import get_logger

logger = get_logger(__name__)


def build_finemap_command(
    z_file: str | Path,
    ld_file: str | Path,
    output_dir: str | Path,
    n_samples: int,
    n_causal_snps: int = 5,
    finemap_binary: str = "finemap",
) -> List[str]:
    """Build a FINEMAP command-line call.

    Parameters
    ----------
    z_file:
        Path to the z-score input file (FINEMAP format).
    ld_file:
        Path to the LD matrix file.
    output_dir:
        Directory for FINEMAP output files.
    n_samples:
        GWAS sample size.
    n_causal_snps:
        Maximum number of causal SNPs to consider (--n-causal-snps).
    finemap_binary:
        Name or path of the FINEMAP binary.

    Returns
    -------
    List[str]
        Command list for subprocess.

    Notes
    -----
    TODO: This is a skeleton builder.  Verify FINEMAP input format requirements
          and add the --log, --sss (stochastic shotgun search) flags as needed.
    """
    logger.warning(
        "[PLACEHOLDER] FINEMAP command builder is a scaffold only.  "
        "Do not use without verifying input file formats."
    )
    cmd = [
        str(finemap_binary),
        "--sss",
        "--in-files",        str(z_file),
        "--ld-file",         str(ld_file),
        "--out-prefix",      str(Path(output_dir) / "finemap_out"),
        "--n-samples",       str(n_samples),
        "--n-causal-snps",   str(n_causal_snps),
    ]
    return cmd


def run_finemap_placeholder(locus_id: str) -> Dict[str, Any]:
    """Return a clearly-labeled placeholder FINEMAP result.

    Parameters
    ----------
    locus_id:
        Identifier for the locus (e.g., ``chr6_31000000``).

    Returns
    -------
    Dict[str, Any]
        Placeholder result.
    """
    logger.warning("[PLACEHOLDER] FINEMAP not run for locus %s.", locus_id)
    return {
        "locus_id": locus_id,
        "status":   "placeholder",
        "method":   "FINEMAP",
        "note":     (
            "FINEMAP requires binary installation and in-sample LD matrices.  "
            "This is a placeholder result only.  "
            "Run FINEMAP after coloc identifies top priority loci."
        ),
        "credible_sets": None,
    }
