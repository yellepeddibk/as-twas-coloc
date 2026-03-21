#!/usr/bin/env python3
"""
run_gwas_intake_validation.py
-----------------------------
Run the first real-data milestone only:
1. Load real AS GWAS summary statistics
2. Harmonize to hg38 with GTEx-style varIDs
3. Save harmonized file
4. Save harmonization QC summary
5. Run one GTEx model varID overlap sanity check
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from src.config import load_config
from src.harmonization import harmonize_gwas
from src.io_utils import ensure_dir, read_gwas, write_tsv_gz
from src.logging_utils import get_logger
from src.validation import check_gtex_model_varid_overlap

logger = get_logger("run_gwas_intake_validation")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run real GWAS intake + harmonization validation milestone")
    parser.add_argument("--config", default="config/as.yaml", help="Path to YAML configuration")
    parser.add_argument("--base-dir", default=".", help="Repository root")
    return parser.parse_args()


def _resolve_gwas_input(base_dir: Path, cfg: dict) -> Path:
    p = Path(cfg.get("gwas", {}).get("input_file", "data/raw/as_gwas.tsv.gz"))
    if p.is_absolute():
        return p
    return (base_dir / p).resolve()


def _build_harmonization_qc(harmonized_gwas: pd.DataFrame) -> dict:
    qc = harmonized_gwas.attrs.get("harmonization_qc", {})
    return {
        "input_snp_count": int(qc.get("input_snp_count", len(harmonized_gwas))),
        "lifted_snp_count": int(qc.get("lifted_snp_count", len(harmonized_gwas))),
        "failed_liftover_count": int(qc.get("failed_liftover_count", 0)),
        "palindromic_dropped_count": int(qc.get("palindromic_dropped_count", 0)),
        "final_retained_snp_count": int(qc.get("final_retained_snp_count", len(harmonized_gwas))),
    }


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()
    cfg = load_config(args.config)
    cfg["_base_dir"] = str(base_dir)

    gwas_input = _resolve_gwas_input(base_dir, cfg)
    if not gwas_input.exists():
        raise FileNotFoundError(f"GWAS input file not found: {gwas_input}")

    raw = read_gwas(gwas_input)
    harmonized = harmonize_gwas(raw, cfg, strict_real_mode=True)

    harm_out = base_dir / "data" / "interim" / "gwas_harmonized" / "as_gwas_hg38_varid.tsv.gz"
    write_tsv_gz(harmonized, harm_out)
    logger.info("Harmonized GWAS written: %s", harm_out)

    qc_dir = base_dir / "data" / "processed" / "qc"
    ensure_dir(qc_dir)

    qc = _build_harmonization_qc(harmonized)
    qc_path = qc_dir / "gwas_harmonization_qc.tsv"
    pd.DataFrame([qc]).to_csv(qc_path, sep="\t", index=False)
    logger.info("Harmonization QC written: %s", qc_path)

    tissue = cfg.get("validation", {}).get("gtex_overlap_tissue", "Whole_Blood")
    model_dir = Path(cfg.get("gtex", {}).get("model_dir", "data/reference/gtex_v8_models"))
    if not model_dir.is_absolute():
        model_dir = base_dir / model_dir

    model_pattern = cfg.get("spredixcan", {}).get("model_db_pattern", "{model_dir}/gtex_v8_mashr_{tissue}.db")
    model_path = Path(model_pattern.format(model_dir=model_dir, tissue=tissue))
    overlap = check_gtex_model_varid_overlap(harmonized, model_path, tissue=tissue)

    overlap_path = qc_dir / f"gwas_gtex_overlap_{tissue}.tsv"
    pd.DataFrame([overlap]).to_csv(overlap_path, sep="\t", index=False)
    logger.info("GTEx overlap QC written: %s", overlap_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
