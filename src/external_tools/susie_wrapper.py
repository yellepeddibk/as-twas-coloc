"""
susie_wrapper.py
----------------
Placeholder wrapper for coloc.susie (SuSiE-based multi-signal colocalization).

coloc.susie is appropriate when a locus may have multiple independent causal
variants.  It requires SuSiE fine-mapping results as input.

Real execution requires:
  1. R + coloc package (>= 5.2, which includes coloc.susie)
  2. susieR package: install.packages("susieR")
  3. Fine-mapped LD matrices for the locus (e.g., from UK Biobank or 1000G)

This placeholder:
  - Documents the interface
  - Returns a clearly-labeled placeholder result
  - Does NOT perform any actual computation

TODO: Implement real coloc.susie support after initial coloc.abf analyses
      identify loci that may be multi-causal (multiple credible sets).
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from src.logging_utils import get_logger

logger = get_logger(__name__)


def run_coloc_susie(
    dataset1: Dict[str, Any],
    dataset2: Dict[str, Any],
    ld_matrix: Optional[Any] = None,
    priors: Optional[Dict[str, float]] = None,
    mock: bool = True,
) -> Dict[str, Any]:
    """Placeholder for coloc.susie for multi-causal loci.

    Parameters
    ----------
    dataset1:
        GWAS coloc dataset dict.
    dataset2:
        eQTL coloc dataset dict.
    ld_matrix:
        LD correlation matrix for the locus (required for real execution).
    priors:
        Coloc prior probabilities dict.
    mock:
        If True (default), return a placeholder result and log a warning.

    Returns
    -------
    Dict[str, Any]
        Placeholder result dict with a clear warning label.

    Raises
    ------
    NotImplementedError
        If mock=False (real implementation not yet available).
    """
    if mock:
        logger.warning(
            "[PLACEHOLDER] coloc.susie not yet implemented. "
            "Returning placeholder result.  "
            "Implement real coloc.susie support after initial coloc.abf analyses."
        )
        return {
            "status": "placeholder",
            "method": "coloc.susie",
            "note": (
                "coloc.susie requires SuSiE fine-mapping and LD matrices. "
                "Not yet implemented in this scaffold.  "
                "Run coloc.abf first; escalate to coloc.susie for multi-causal loci."
            ),
            "PP.H4.abf": None,
        }

    raise NotImplementedError(
        "Real coloc.susie execution is not yet implemented.  "
        "Install susieR and implement the R interface in run_coloc_susie.R."
    )
