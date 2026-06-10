"""Focused tests for real-analysis prerequisite tissue selection."""

from argparse import Namespace
from pathlib import Path

import pytest

from scripts.check_real_analysis_prereqs import (
    POSTER_TISSUES,
    _build_prereq_rows,
    _find_eqtl_file_for_tissue,
    _selected_tissues,
    _write_report_tsv,
)


def test_selected_tissues_defaults_to_configured_tissues():
    cfg = {"gtex": {"tissues": ["Whole_Blood", "Spleen"]}}
    args = Namespace(tissues=None, poster_tissues=False)

    assert _selected_tissues(cfg, args) == ["Whole_Blood", "Spleen"]


def test_selected_tissues_uses_cli_override():
    cfg = {"gtex": {"tissues": ["Whole_Blood", "Spleen"]}}
    args = Namespace(tissues="Lung, Colon_Sigmoid", poster_tissues=False)

    assert _selected_tissues(cfg, args) == ["Lung", "Colon_Sigmoid"]


def test_selected_tissues_can_use_poster_subset():
    cfg = {"gtex": {"tissues": ["Whole_Blood", "Spleen"]}}
    args = Namespace(tissues=None, poster_tissues=True)

    assert _selected_tissues(cfg, args) == POSTER_TISSUES


def test_selected_tissues_requires_configured_tissues_when_not_overridden():
    args = Namespace(tissues=None, poster_tissues=False)

    with pytest.raises(ValueError, match="gtex.tissues"):
        _selected_tissues({}, args)


def test_find_eqtl_file_prefers_allpairs_file(tmp_path):
    eqtl_dir = tmp_path / "eqtl"
    eqtl_dir.mkdir()
    fallback = eqtl_dir / "Whole_Blood.small.tsv.gz"
    preferred = eqtl_dir / "Whole_Blood.allpairs.txt.gz"
    fallback.write_text("fallback\n")
    preferred.write_text("preferred\n")

    assert _find_eqtl_file_for_tissue(eqtl_dir, "Whole_Blood") == preferred


def test_build_prereq_rows_checks_real_run_assets(tmp_path, monkeypatch):
    base_dir = tmp_path
    gwas_input = base_dir / "data" / "raw" / "as_gwas.tsv.gz"
    chain_file = base_dir / "data" / "reference" / "chains" / "hg19ToHg38.over.chain.gz"
    metaxcan_script = base_dir / "external" / "MetaXcan" / "software" / "SPrediXcan.py"
    model_dir = base_dir / "data" / "reference" / "gtex_v8_models"
    eqtl_dir = base_dir / "data" / "reference" / "gtex_v8_eqtl"

    for path in (
        gwas_input,
        chain_file,
        metaxcan_script,
        model_dir / "gtex_v8_mashr_Whole_Blood.db",
        model_dir / "gtex_v8_mashr_Whole_Blood.txt.gz",
        eqtl_dir / "Whole_Blood.allpairs.txt.gz",
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("ok\n")

    monkeypatch.setattr("scripts.check_real_analysis_prereqs.shutil.which", lambda cmd: "/usr/bin/Rscript" if cmd == "Rscript" else None)

    cfg = {
        "project": {"genome_build": "hg38"},
        "gwas": {"input_file": "data/raw/as_gwas.tsv.gz", "source_build": "hg19"},
        "liftover": {"enabled": True, "chain_file": "data/reference/chains/hg19ToHg38.over.chain.gz"},
        "gtex": {
            "model_dir": "data/reference/gtex_v8_models",
            "eqtl_dir": "data/reference/gtex_v8_eqtl",
        },
        "spredixcan": {
            "script": "external/MetaXcan/software/SPrediXcan.py",
            "model_db_pattern": "{model_dir}/gtex_v8_mashr_{tissue}.db",
            "covariance_pattern": "{model_dir}/gtex_v8_mashr_{tissue}.txt.gz",
        },
    }

    rows = _build_prereq_rows(cfg, base_dir, ["Whole_Blood"], check_r_packages=False)
    status = {(row.category, row.name): row.ok for row in rows}

    assert status[("input", "AS GWAS summary statistics")]
    assert status[("liftover", "LiftOver chain")]
    assert status[("tool", "MetaXcan script")]
    assert status[("tool", "Rscript executable")]
    assert status[("model", "Model DB (Whole_Blood)")]
    assert status[("model", "Covariance (Whole_Blood)")]
    assert status[("eqtl", "eQTL summary (Whole_Blood)")]


def test_write_report_tsv_writes_machine_readable_rows(tmp_path):
    from scripts.check_real_analysis_prereqs import PrereqRow

    out_path = tmp_path / "reports" / "prereqs.tsv"
    _write_report_tsv(
        [PrereqRow("tool", "Rscript executable", "/usr/bin/Rscript", True)],
        out_path,
    )

    text = out_path.read_text()
    assert "category\tname\ttarget\tok\trequired\tdetail" in text
    assert "tool\tRscript executable\t/usr/bin/Rscript\tTrue\tTrue" in text