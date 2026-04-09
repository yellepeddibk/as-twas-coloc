#!/usr/bin/env python3
"""
check_real_analysis_prereqs.py
------------------------------
Validate prerequisites for real TWAS + COLOC poster analysis.

Checks:
- MetaXcan script present
- Rscript available
- GTEx model DB + covariance files for poster tissues
- GTEx eQTL directory exists and is non-empty
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
import sys

_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from src.config import load_config

POSTER_TISSUES = [
    "Whole_Blood",
    "Spleen",
    "Small_Intestine_Terminal_Ileum",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Check prerequisites for real TWAS + COLOC runs")
    p.add_argument("--config", default="config/as.yaml", help="Config path")
    p.add_argument("--base-dir", default=".", help="Repository root")
    return p.parse_args()


def _ok(flag: bool) -> str:
    return "OK" if flag else "MISSING"


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()
    cfg = load_config(args.config)

    gtex_cfg = cfg.get("gtex", {})
    sp_cfg = cfg.get("spredixcan", {})

    model_dir = Path(gtex_cfg.get("model_dir", "data/reference/gtex_v8_models"))
    eqtl_dir = Path(gtex_cfg.get("eqtl_dir", "data/reference/gtex_v8_eqtl"))
    if not model_dir.is_absolute():
        model_dir = base_dir / model_dir
    if not eqtl_dir.is_absolute():
        eqtl_dir = base_dir / eqtl_dir

    model_pattern = sp_cfg.get("model_db_pattern", "{model_dir}/gtex_v8_mashr_{tissue}.db")
    cov_pattern = sp_cfg.get("covariance_pattern", "{model_dir}/gtex_v8_mashr_{tissue}.txt.gz")

    script_path = Path(sp_cfg.get("script", "external/MetaXcan/software/SPrediXcan.py"))
    if not script_path.is_absolute():
        script_path = base_dir / script_path

    rows = []

    rows.append(("MetaXcan script", str(script_path), script_path.exists()))
    rows.append(("Rscript executable", "Rscript", shutil.which("Rscript") is not None))
    rows.append(("GTEx model dir", str(model_dir), model_dir.exists()))
    rows.append(("GTEx eQTL dir", str(eqtl_dir), eqtl_dir.exists() and any(eqtl_dir.iterdir()) if eqtl_dir.exists() else False))

    for tissue in POSTER_TISSUES:
        db_path = Path(model_pattern.format(model_dir=model_dir, tissue=tissue))
        cov_path = Path(cov_pattern.format(model_dir=model_dir, tissue=tissue))
        rows.append((f"Model DB ({tissue})", str(db_path), db_path.exists()))
        rows.append((f"Covariance ({tissue})", str(cov_path), cov_path.exists()))

    missing = [r for r in rows if not r[2]]

    print("Real analysis prerequisite check")
    print("=" * 40)
    for name, target, ok in rows:
        print(f"[{_ok(ok):7}] {name}: {target}")

    print("=" * 40)
    if missing:
        print(f"Missing items: {len(missing)}")
        return 1

    print("All required prerequisites are present.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
