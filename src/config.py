"""
config.py
---------
Load and validate the pipeline YAML configuration.
"""

import os
from pathlib import Path
from typing import Any, Dict

import yaml

from src.logging_utils import get_logger

logger = get_logger(__name__)

DEFAULT_CONFIG_PATH = Path(__file__).parent.parent / "config" / "as.yaml"


def load_config(config_path: str | Path | None = None) -> Dict[str, Any]:
    """Load YAML configuration and return as a nested dictionary.

    Parameters
    ----------
    config_path:
        Path to the YAML config file.  If *None*, uses
        ``config/as.yaml`` relative to the repository root.

    Returns
    -------
    Dict[str, Any]
        Parsed configuration dictionary.

    Raises
    ------
    FileNotFoundError
        If the config file does not exist.
    """
    if config_path is None:
        config_path = DEFAULT_CONFIG_PATH

    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with config_path.open() as fh:
        cfg = yaml.safe_load(fh)

    logger.info("Loaded config from %s", config_path)
    return cfg


def resolve_path(cfg: Dict[str, Any], key: str, base_dir: str | Path | None = None) -> Path:
    """Resolve a path from the config, optionally prepending a base directory.

    Parameters
    ----------
    cfg:
        Full config dict.
    key:
        Dot-separated key path into cfg, e.g. ``"paths.data_raw"``.
    base_dir:
        If provided, the resolved path is relative to this directory.

    Returns
    -------
    Path
    """
    parts = key.split(".")
    value = cfg
    for part in parts:
        value = value[part]

    p = Path(value)  # type: ignore[arg-type]
    if base_dir is not None and not p.is_absolute():
        p = Path(base_dir) / p

    return p


def get_coloc_priors(cfg: Dict[str, Any]) -> Dict[str, float]:
    """Return the coloc prior probabilities from config."""
    coloc = cfg.get("coloc", {})
    return {
        "p1": float(coloc.get("p1", 1e-4)),
        "p2": float(coloc.get("p2", 1e-4)),
        "p12": float(coloc.get("p12", 1e-5)),
    }


def get_tissues(cfg: Dict[str, Any]) -> list:
    """Return the list of GTEx tissues from config."""
    return cfg.get("gtex", {}).get("tissues", [])
