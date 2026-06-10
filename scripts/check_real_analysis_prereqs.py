#!/usr/bin/env python3
"""
check_real_analysis_prereqs.py
------------------------------
Validate prerequisites for configured real TWAS + COLOC analyses.

Checks:
- AS GWAS input present
- liftover chain present when configured
- MetaXcan script present
- Rscript available
- R coloc/jsonlite packages available
- GTEx model DB + covariance files for configured tissues
- GTEx eQTL summary file present for configured tissues
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import shutil
import subprocess
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
    p.add_argument(
        "--report-tsv",
        help="Optional path for a machine-readable prerequisite report TSV",
    )
    p.add_argument(
        "--skip-r-package-check",
        action="store_true",
        help="Only check that Rscript exists; skip R package import checks",
    )
    tissue_group = p.add_mutually_exclusive_group()
    tissue_group.add_argument(
        "--tissues",
        help="Comma-separated tissue override (default: use gtex.tissues from config)",
    )
    tissue_group.add_argument(
        "--poster-tissues",
        action="store_true",
        help="Check only the poster tissue subset",
    )
    return p.parse_args()


@dataclass(frozen=True)
class PrereqRow:
    category: str
    name: str
    target: str
    ok: bool
    required: bool = True
    detail: str = ""


def _ok(flag: bool) -> str:
    return "OK" if flag else "MISSING"


def _selected_tissues(cfg: dict, args: argparse.Namespace) -> list[str]:
    """Resolve which tissues should be checked for required GTEx assets."""
    if getattr(args, "tissues", None):
        tissues = [t.strip() for t in args.tissues.split(",") if t.strip()]
        if tissues:
            return tissues

    if getattr(args, "poster_tissues", False):
        return POSTER_TISSUES

    tissues = cfg.get("gtex", {}).get("tissues", [])
    if not tissues:
        raise ValueError("No tissues configured under gtex.tissues")
    return tissues


def _resolve(base_dir: Path, path_value: str | Path) -> Path:
    path = Path(path_value)
    return path if path.is_absolute() else base_dir / path


def _nonempty_dir(path: Path) -> bool:
    return path.exists() and path.is_dir() and any(path.iterdir())


def _r_package_available(package: str) -> bool:
    if shutil.which("Rscript") is None:
        return False

    cmd = [
        "Rscript",
        "-e",
        f"quit(status = ifelse(requireNamespace('{package}', quietly=TRUE), 0, 1))",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except OSError:
        return False
    return result.returncode == 0


def _find_eqtl_file_for_tissue(eqtl_dir: Path, tissue: str) -> Path | None:
    if not eqtl_dir.exists() or not eqtl_dir.is_dir():
        return None

    tissue_norm = tissue.lower()
    valid_suffixes = (".gz", ".txt", ".tsv")
    candidates = sorted(
        path
        for path in eqtl_dir.rglob("*")
        if path.is_file()
        and tissue_norm in path.name.lower()
        and path.name.lower().endswith(valid_suffixes)
    )
    if not candidates:
        return None

    preferred = [path for path in candidates if "allpairs" in path.name.lower()]
    return preferred[0] if preferred else candidates[0]


def _build_prereq_rows(
    cfg: dict,
    base_dir: Path,
    tissues: list[str],
    *,
    check_r_packages: bool = True,
) -> list[PrereqRow]:
    project_cfg = cfg.get("project", {})
    gwas_cfg = cfg.get("gwas", {})
    gtex_cfg = cfg.get("gtex", {})
    liftover_cfg = cfg.get("liftover", {})
    sp_cfg = cfg.get("spredixcan", {})

    rows: list[PrereqRow] = []

    gwas_input = _resolve(base_dir, gwas_cfg.get("input_file", "data/raw/as_gwas.tsv.gz"))
    rows.append(PrereqRow("input", "AS GWAS summary statistics", str(gwas_input), gwas_input.exists()))

    source_build = str(gwas_cfg.get("source_build", "")).lower()
    target_build = str(project_cfg.get("genome_build", "")).lower()
    liftover_required = bool(liftover_cfg.get("enabled", False)) and source_build and target_build and source_build != target_build
    if liftover_required:
        chain_path = _resolve(base_dir, liftover_cfg.get("chain_file", "data/reference/chains/hg19ToHg38.over.chain.gz"))
        rows.append(PrereqRow("liftover", "LiftOver chain", str(chain_path), chain_path.exists()))

    script_path = _resolve(base_dir, sp_cfg.get("script", "external/MetaXcan/software/SPrediXcan.py"))
    rows.append(PrereqRow("tool", "MetaXcan script", str(script_path), script_path.exists()))

    rscript_path = shutil.which("Rscript")
    rows.append(PrereqRow("tool", "Rscript executable", rscript_path or "Rscript", rscript_path is not None))
    if check_r_packages:
        for package in ("coloc", "jsonlite"):
            rows.append(PrereqRow("r-package", package, f"R package: {package}", _r_package_available(package)))

    model_dir = _resolve(base_dir, gtex_cfg.get("model_dir", "data/reference/gtex_v8_models"))
    eqtl_dir = _resolve(base_dir, gtex_cfg.get("eqtl_dir", "data/reference/gtex_v8_eqtl"))
    rows.append(PrereqRow("reference", "GTEx model directory", str(model_dir), model_dir.exists() and model_dir.is_dir()))
    rows.append(PrereqRow("reference", "GTEx eQTL directory", str(eqtl_dir), _nonempty_dir(eqtl_dir)))

    model_pattern = sp_cfg.get("model_db_pattern", "{model_dir}/gtex_v8_mashr_{tissue}.db")
    cov_pattern = sp_cfg.get("covariance_pattern", "{model_dir}/gtex_v8_mashr_{tissue}.txt.gz")

    for tissue in tissues:
        db_path = Path(model_pattern.format(model_dir=model_dir, tissue=tissue))
        cov_path = Path(cov_pattern.format(model_dir=model_dir, tissue=tissue))
        eqtl_path = _find_eqtl_file_for_tissue(eqtl_dir, tissue)

        rows.append(PrereqRow("model", f"Model DB ({tissue})", str(db_path), db_path.exists()))
        rows.append(PrereqRow("model", f"Covariance ({tissue})", str(cov_path), cov_path.exists()))
        rows.append(
            PrereqRow(
                "eqtl",
                f"eQTL summary ({tissue})",
                str(eqtl_path or eqtl_dir / f"{tissue}.allpairs.txt.gz"),
                eqtl_path is not None,
            )
        )

    return rows


def _write_report_tsv(rows: list[PrereqRow], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["category", "name", "target", "ok", "required", "detail"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "category": row.category,
                "name": row.name,
                "target": row.target,
                "ok": row.ok,
                "required": row.required,
                "detail": row.detail,
            })


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = base_dir / config_path
    cfg = load_config(config_path)

    tissues = _selected_tissues(cfg, args)
    rows = _build_prereq_rows(
        cfg,
        base_dir,
        tissues,
        check_r_packages=not args.skip_r_package_check,
    )
    missing = [row for row in rows if row.required and not row.ok]

    if args.report_tsv:
        report_path = _resolve(base_dir, args.report_tsv)
        _write_report_tsv(rows, report_path)

    print("Real analysis prerequisite check")
    print("=" * 40)
    print(f"Tissues checked: {', '.join(tissues)}")
    for row in rows:
        print(f"[{_ok(row.ok):7}] {row.category} / {row.name}: {row.target}")

    print("=" * 40)
    if args.report_tsv:
        print(f"Report: {report_path}")
    if missing:
        print(f"Missing items: {len(missing)}")
        return 1

    print("All required prerequisites are present.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
