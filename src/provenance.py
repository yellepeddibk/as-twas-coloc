"""
provenance.py
-------------
Record a run manifest (JSON) capturing metadata for reproducibility.

Each pipeline run saves a manifest to results/manifests/ containing:
- timestamp
- git commit hash (if available)
- config snapshot
- tool version placeholders
- input file paths and checksum placeholders
- decision thresholds
- coloc priors
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.logging_utils import get_logger

logger = get_logger(__name__)


def _git_commit() -> Optional[str]:
    """Return the current git commit hash, or None if not in a git repo."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except Exception:
        return None


def _file_md5(path: str | Path) -> Optional[str]:
    """Compute MD5 checksum of a file; return None if file not found."""
    path = Path(path)
    if not path.exists():
        return None
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def build_manifest(
    cfg: Dict[str, Any],
    input_files: Optional[List[str | Path]] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build a provenance manifest dict for the current run.

    Parameters
    ----------
    cfg:
        Full pipeline configuration dict.
    input_files:
        List of input file paths whose checksums should be recorded.
    extra:
        Additional key-value pairs to merge into the manifest.

    Returns
    -------
    Dict[str, Any]
        Manifest dictionary.
    """
    coloc_cfg = cfg.get("coloc", {})
    correction_cfg = cfg.get("correction", {})
    gwas_cfg = cfg.get("gwas", {})

    checksums: Dict[str, Optional[str]] = {}
    for f in (input_files or []):
        f = Path(f)
        checksums[str(f)] = _file_md5(f)

    manifest: Dict[str, Any] = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(),
        "python_version": sys.version,
        "config_snapshot": cfg,
        "tool_versions": {
            "numpy": _try_version("numpy"),
            "pandas": _try_version("pandas"),
            "scipy": _try_version("scipy"),
            "statsmodels": _try_version("statsmodels"),
            # External tools are not installed in CI; record as None
            "MetaXcan/SPrediXcan": None,  # TODO: record after installation
            "R": None,                     # TODO: record R version
            "coloc_R_package": None,       # TODO: record coloc package version
        },
        "input_files": {str(f): str(f) for f in (input_files or [])},
        "input_checksums": checksums,
        "decision_thresholds": {
            "fdr_alpha": correction_cfg.get("alpha", 0.05),
            "fdr_method": correction_cfg.get("method", "fdr_bh"),
            "pp4_threshold": coloc_cfg.get("pp4_threshold", 0.7),
            "pp3_pp4_threshold": coloc_cfg.get("pp3_pp4_threshold", 0.8),
            "pp4_ratio_threshold": coloc_cfg.get("pp4_ratio_threshold", 0.9),
        },
        "coloc_priors": {
            "p1":  coloc_cfg.get("p1",  1e-4),
            "p2":  coloc_cfg.get("p2",  1e-4),
            "p12": coloc_cfg.get("p12", 1e-5),
        },
        "gwas_metadata": {
            "source_build": gwas_cfg.get("source_build"),
            "n_total":      gwas_cfg.get("n_total"),
            "n_cases":      gwas_cfg.get("n_cases"),
            "n_controls":   gwas_cfg.get("n_controls"),
        },
    }

    if extra:
        manifest.update(extra)

    return manifest


def save_manifest(manifest: Dict[str, Any], manifest_dir: str | Path, run_id: str) -> Path:
    """Write the manifest to a JSON file.

    Parameters
    ----------
    manifest:
        Manifest dict from ``build_manifest``.
    manifest_dir:
        Directory to write the manifest file into.
    run_id:
        Unique run identifier (e.g., timestamp string) used in the filename.

    Returns
    -------
    Path
        Path to the written manifest file.
    """
    manifest_dir = Path(manifest_dir)
    manifest_dir.mkdir(parents=True, exist_ok=True)

    out_path = manifest_dir / f"manifest_{run_id}.json"
    with out_path.open("w") as fh:
        json.dump(manifest, fh, indent=2, default=_json_default)

    logger.info("Provenance manifest saved to %s.", out_path)
    return out_path


def _try_version(package: str) -> Optional[str]:
    """Return the version of an installed package, or None."""
    try:
        import importlib.metadata
        return importlib.metadata.version(package)
    except Exception:
        return None


def _json_default(obj: Any) -> Any:
    """JSON serialization fallback."""
    if hasattr(obj, "isoformat"):
        return obj.isoformat()
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)
