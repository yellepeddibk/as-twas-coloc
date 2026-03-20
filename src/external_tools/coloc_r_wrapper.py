"""
coloc_r_wrapper.py
------------------
Python interface to the R coloc package via subprocess.

Real execution requires:
  1. R installed (>= 4.0)
  2. coloc R package: install.packages("coloc")
  3. Optional: coloc.susie for multi-causal regions

This wrapper:
  - Accepts Python dataset dicts (from coloc.build_coloc_dataset)
  - Serializes them to temporary JSON
  - Invokes an R helper script (src/external_tools/run_coloc.R)
  - Parses the JSON result back into Python

In mock mode (mock=True), returns a mock posterior dict without calling R.

Pitfall reminders (in code comments):
- Always pass varbeta (= se**2), not se, to coloc in R.
- Do not pre-filter region SNPs by p-value before calling coloc.
- coloc.abf is the primary method; coloc.susie is an optional extension.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from src.logging_utils import get_logger

logger = get_logger(__name__)

# Path to the companion R script (relative to repo root)
_DEFAULT_R_SCRIPT = Path(__file__).parent / "run_coloc.R"


def run_coloc_abf(
    dataset1: Dict[str, Any],
    dataset2: Dict[str, Any],
    priors: Dict[str, float],
    r_script: Optional[str | Path] = None,
    mock: bool = False,
) -> Dict[str, Any]:
    """Run coloc.abf for a single gene-locus pair.

    Parameters
    ----------
    dataset1:
        GWAS coloc dataset dict (keys: snp, beta, varbeta, N, type, [s]).
    dataset2:
        eQTL coloc dataset dict (keys: snp, beta, varbeta, N, type).
    priors:
        Dict with keys p1, p2, p12.
    r_script:
        Path to the R helper script.  Defaults to ``run_coloc.R`` in this dir.
    mock:
        If True, return a mock result without calling R.

    Returns
    -------
    Dict[str, Any]
        Coloc posterior probabilities (PP.H0 .. PP.H4 and nsnps).

    Raises
    ------
    RuntimeError
        If R is not available and mock=False.
    """
    if mock or (not dataset1 or not dataset2):
        logger.info("[MOCK] Returning mock coloc.abf result (R not called).")
        from src.coloc import make_mock_coloc_result
        return make_mock_coloc_result()

    r_script = Path(r_script) if r_script else _DEFAULT_R_SCRIPT
    if not r_script.exists():
        raise RuntimeError(
            f"R coloc script not found at {r_script}. "
            "Ensure run_coloc.R exists and R + coloc package are installed."
        )

    # Validate varbeta is present (common pitfall: passing se instead)
    for ds_name, ds in (("dataset1", dataset1), ("dataset2", dataset2)):
        if "varbeta" not in ds:
            raise ValueError(
                f"{ds_name} is missing 'varbeta'.  "
                "Compute varbeta = se**2 before calling coloc.  "
                "Do NOT pass se directly to coloc."
            )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        input_path  = tmpdir / "coloc_input.json"
        output_path = tmpdir / "coloc_output.json"

        payload = {
            "dataset1": dataset1,
            "dataset2": dataset2,
            "priors":   priors,
        }
        with input_path.open("w") as fh:
            json.dump(payload, fh)

        cmd = ["Rscript", str(r_script), str(input_path), str(output_path)]
        logger.info("Running coloc via R: %s", " ".join(cmd))

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=300,
            )
        except FileNotFoundError:
            raise RuntimeError(
                "Rscript not found. Install R and the coloc package before running coloc."
            )
        except subprocess.CalledProcessError as exc:
            logger.error("R coloc failed:\n%s\n%s", exc.stdout, exc.stderr)
            raise

        with output_path.open() as fh:
            coloc_result = json.load(fh)

    logger.info("coloc.abf completed: PP4=%.3f", coloc_result.get("PP.H4.abf", float("nan")))
    return coloc_result
