"""
io_utils.py
-----------
I/O helpers for reading and writing pipeline data files.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import Optional

import pandas as pd

from src.logging_utils import get_logger

logger = get_logger(__name__)


def read_gwas(path: str | Path, sep: str = "\t", **kwargs) -> pd.DataFrame:
    """Read a GWAS summary statistics file (TSV or TSV.GZ).

    Parameters
    ----------
    path:
        File path.  Gzipped files are detected automatically.
    sep:
        Field separator (default: tab).
    **kwargs:
        Additional arguments forwarded to ``pd.read_csv``.

    Returns
    -------
    pd.DataFrame
    """
    path = Path(path)
    logger.info("Reading GWAS file: %s", path)
    if not path.exists():
        raise FileNotFoundError(f"GWAS file not found: {path}")

    df = pd.read_csv(path, sep=sep, low_memory=False, **kwargs)
    logger.info("  Loaded %d rows × %d columns", len(df), len(df.columns))
    return df


def write_tsv_gz(df: pd.DataFrame, path: str | Path, **kwargs) -> None:
    """Write a DataFrame to a gzipped TSV file, creating parent directories.

    Parameters
    ----------
    df:
        DataFrame to write.
    path:
        Output file path (will be gzipped regardless of extension).
    **kwargs:
        Additional arguments forwarded to ``DataFrame.to_csv``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip", **kwargs)
    logger.info("Wrote %d rows to %s", len(df), path)


def read_tsv_gz(path: str | Path, **kwargs) -> pd.DataFrame:
    """Read a gzipped TSV file produced by this pipeline.

    Parameters
    ----------
    path:
        File path.
    **kwargs:
        Additional arguments forwarded to ``pd.read_csv``.

    Returns
    -------
    pd.DataFrame
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False, **kwargs)
    logger.info("Read %d rows from %s", len(df), path)
    return df


def ensure_dir(path: str | Path) -> Path:
    """Create directory (and parents) if it does not exist.

    Parameters
    ----------
    path:
        Directory path.

    Returns
    -------
    Path
    """
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def safe_output_path(
    directory: str | Path,
    filename: str,
    overwrite: bool = True,
) -> Path:
    """Construct an output path, optionally checking for existing files.

    Parameters
    ----------
    directory:
        Target directory.
    filename:
        File name.
    overwrite:
        If False and the file already exists, raises FileExistsError.

    Returns
    -------
    Path
    """
    directory = ensure_dir(directory)
    p = directory / filename
    if not overwrite and p.exists():
        raise FileExistsError(f"Output file already exists: {p}")
    return p
