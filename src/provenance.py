"""
provenance.py
-------------
Record a run manifest (JSON) capturing metadata for reproducibility.

Each pipeline run saves a manifest to results/manifests/ containing:
- timestamp
- git commit hash (if available)
- config snapshot
- Python and external tool versions
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


def _git_commit_for_repo(repo_path: str | Path) -> Optional[str]:
    """Return the git commit hash for a specific repository path."""
    repo_path = Path(repo_path)
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_path), "rev-parse", "HEAD"],
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


def _command_output(args: List[str]) -> Optional[str]:
    """Run a command and return stripped stdout or stderr, or None on failure."""
    try:
        result = subprocess.run(
            args,
            capture_output=True,
            text=True,
            check=True,
        )
    except Exception:
        return None

    stdout = result.stdout.strip()
    stderr = result.stderr.strip()
    return stdout or stderr or None


def _resolve_path(path: str | Path, cfg: Dict[str, Any]) -> Path:
    """Resolve a repo-relative path using cfg['_base_dir'] when available."""
    path = Path(path)
    if path.is_absolute():
        return path

    base_dir = Path(cfg.get("_base_dir", ".")).resolve()
    return base_dir / path


def _find_git_root(start_path: str | Path) -> Optional[Path]:
    """Walk upward from a path until a git root is found."""
    start_path = Path(start_path).resolve()
    candidates = [start_path] + list(start_path.parents)
    for candidate in candidates:
        if (candidate / ".git").exists():
            return candidate
    return None


def _metaxcan_version(cfg: Dict[str, Any]) -> Optional[str]:
    """Return a MetaXcan/SPrediXcan repository commit if available."""
    sp_cfg = cfg.get("spredixcan", {})
    script_path = _resolve_path(
        sp_cfg.get("script", "external/MetaXcan/software/SPrediXcan.py"),
        cfg,
    )
    git_root = _find_git_root(script_path.parent)
    if git_root is None:
        return None

    commit = _git_commit_for_repo(git_root)
    if commit is None:
        return None

    return f"git:{commit[:7]}"


def _rscript_version() -> Optional[str]:
    """Return the installed Rscript version string, or None."""
    output = _command_output(["Rscript", "--version"])
    if output is None:
        return None
    return output.splitlines()[0].strip()


def _r_package_version(package: str) -> Optional[str]:
    """Return the installed R package version, or None if unavailable."""
    output = _command_output(
        ["Rscript", "-e", f"cat(as.character(packageVersion('{package}')))"],
    )
    return output.strip() if output else None


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
            "MetaXcan/SPrediXcan": _metaxcan_version(cfg),
            "R": _rscript_version(),
            "coloc_R_package": _r_package_version("coloc"),
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
