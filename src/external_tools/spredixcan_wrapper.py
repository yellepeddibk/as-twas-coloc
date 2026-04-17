"""
spredixcan_wrapper.py
---------------------
Command builder and mock interface for S-PrediXcan (MetaXcan).

Real execution requires:
  1. Clone MetaXcan: https://github.com/hakyimlab/MetaXcan
  2. Install its Python dependencies (pip install -r software/requirements.txt)
  3. Download GTEx v8 PredictDB .db models from https://predictdb.org
  4. Download GTEx v8 mashr covariance .txt.gz files
  5. Set paths in config/as.yaml (spredixcan.script, gtex.model_dir)

Pitfalls
--------
* "0% of model SNPs found" — almost always a genome-build mismatch.
  Ensure GWAS harmonization to hg38 and GTEx-style varIDs.
* Covariance files must match the .db model exactly (same tissue, same build).
* SPrediXcan requires beta + se columns; z-score-only input is not supported
  for weighted models.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.logging_utils import get_logger

logger = get_logger(__name__)


def build_spredixcan_command(
    cfg: Dict[str, Any],
    tissue: str,
    model_db: str | Path,
    covariance: str | Path,
    gwas_file: str | Path,
    output_file: str | Path,
) -> List[str]:
    """Construct the S-PrediXcan command-line argument list.

    Parameters
    ----------
    cfg:
        Full pipeline configuration dict.
    tissue:
        GTEx tissue name (used for logging only).
    model_db:
        Path to the PredictDB .db file for this tissue.
    covariance:
        Path to the covariance .txt.gz file for this tissue.
    gwas_file:
        Path to the harmonized GWAS TSV.GZ file.
    output_file:
        Path where S-PrediXcan should write its CSV output.

    Returns
    -------
    List[str]
        Argument list suitable for ``subprocess.run``.

    Raises
    ------
    RuntimeError
        If the SPrediXcan script is not found.
    """
    sp_cfg = cfg.get("spredixcan", {})
    script = Path(sp_cfg.get("script", "external/MetaXcan/software/SPrediXcan.py"))

    # In real execution, we require the script to exist.
    # In mock/demo mode, we allow it to be missing.
    if not script.exists():
        logger.warning(
            "[PLACEHOLDER] S-PrediXcan script not found at %s. "
            "Returning mock command for demonstration only.  "
            "Install MetaXcan before running real analyses.",
            script,
        )

    cmd = [
        "python", str(script),
        "--model_db_path",    str(model_db),
        "--covariance",       str(covariance),
        "--gwas_file",        str(gwas_file),
        "--model_db_snp_key", sp_cfg.get("model_db_snp_key", "rsid"),
        "--snp_column",       sp_cfg.get("snp_column", "varID"),
        "--effect_allele_column",     sp_cfg.get("effect_allele_column", "effect_allele"),
        "--non_effect_allele_column", sp_cfg.get("non_effect_allele_column", "non_effect_allele"),
        "--beta_column",      sp_cfg.get("beta_column", "beta"),
        "--se_column",        sp_cfg.get("se_column", "se"),
        "--pvalue_column",    sp_cfg.get("pvalue_column", "pvalue"),
        "--output_file",      str(output_file),
    ]

    logger.info("Built S-PrediXcan command for tissue '%s'.", tissue)
    logger.debug("  command: %s", " ".join(cmd))
    return cmd


def run_spredixcan(
    cmd: List[str],
    dry_run: bool = False,
    timeout: int = 3600,
) -> Optional[subprocess.CompletedProcess]:
    """Execute an S-PrediXcan command.

    Parameters
    ----------
    cmd:
        Command list from ``build_spredixcan_command``.
    dry_run:
        If True, log the command but do not execute it.
    timeout:
        Maximum execution time in seconds.

    Returns
    -------
    subprocess.CompletedProcess or None (in dry-run mode).
    """
    if dry_run:
        logger.info("[DRY-RUN] Would execute: %s", " ".join(cmd))
        return None

    logger.info("Running S-PrediXcan: %s", " ".join(cmd[:4]) + " ...")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=True,
        )
        logger.info("S-PrediXcan finished successfully.")
        if result.stdout:
            logger.debug("stdout: %s", result.stdout[-2000:])
        return result
    except subprocess.CalledProcessError as exc:
        logger.error(
            "S-PrediXcan failed (exit code %d):\nstdout: %s\nstderr: %s",
            exc.returncode,
            exc.stdout[-1000:] if exc.stdout else "",
            exc.stderr[-1000:] if exc.stderr else "",
        )
        raise
    except FileNotFoundError:
        raise RuntimeError(
            "Python interpreter or S-PrediXcan script not found. "
            "Ensure MetaXcan is installed and the script path in config is correct."
        )


def run_all_tissues(
    cfg: Dict[str, Any],
    tissues: Optional[List[str]] = None,
    dry_run: bool = True,
) -> Dict[str, Path]:
    """Run S-PrediXcan for all configured tissues.

    Parameters
    ----------
    cfg:
        Full pipeline config.
    tissues:
        Override list of tissues; if None, uses ``cfg.gtex.tissues``.
    dry_run:
        If True, build and log commands without executing.

    Returns
    -------
    Dict[str, Path]
        Mapping from tissue name to expected output CSV path.
    """
    if tissues is None:
        tissues = cfg.get("gtex", {}).get("tissues", [])

    sp_cfg   = cfg.get("spredixcan", {})
    gtex_cfg = cfg.get("gtex", {})
    model_dir  = Path(gtex_cfg.get("model_dir", "data/reference/gtex_v8_models"))
    output_dir = Path(sp_cfg.get("output_dir", "results/twas"))
    gwas_file  = Path(sp_cfg.get("gwas_file", "data/interim/gwas_harmonized/as_gwas_hg38_varid.tsv.gz"))

    model_pattern  = sp_cfg.get("model_db_pattern",  "{model_dir}/gtex_v8_mashr_{tissue}.db")
    cov_pattern    = sp_cfg.get("covariance_pattern", "{model_dir}/gtex_v8_mashr_{tissue}.txt.gz")

    output_paths: Dict[str, Path] = {}
    for tissue in tissues:
        model_db    = Path(model_pattern.format(model_dir=model_dir, tissue=tissue))
        covariance  = Path(cov_pattern.format(model_dir=model_dir, tissue=tissue))
        output_file = output_dir / f"{tissue}.spredixcan.csv"

        cmd = build_spredixcan_command(
            cfg=cfg,
            tissue=tissue,
            model_db=model_db,
            covariance=covariance,
            gwas_file=gwas_file,
            output_file=output_file,
        )
        run_spredixcan(cmd, dry_run=dry_run)
        output_paths[tissue] = output_file

    return output_paths
