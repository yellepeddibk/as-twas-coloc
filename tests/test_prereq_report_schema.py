"""Schema tests for the real-analysis prerequisite report TSV.

Guarantees the TSV format produced by ``scripts/check_real_analysis_prereqs.py``
stays stable so downstream consumers can rely on the columns.
"""

from __future__ import annotations

import csv
from pathlib import Path

from scripts.check_real_analysis_prereqs import PrereqRow, _write_report_tsv

_EXPECTED_FIELDS = ["category", "name", "target", "ok", "required", "detail"]


def test_report_tsv_has_required_columns(tmp_path):
    out_path = tmp_path / "report.tsv"
    rows = [
        PrereqRow("input", "AS GWAS summary statistics", "/data/gwas.tsv.gz", True),
        PrereqRow("tool", "Rscript executable", "/usr/bin/Rscript", False, detail="not found"),
    ]
    _write_report_tsv(rows, out_path)

    with out_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        assert reader.fieldnames == _EXPECTED_FIELDS, (
            f"Unexpected columns: {reader.fieldnames}"
        )
        body = list(reader)

    assert len(body) == 2
    assert body[0]["category"] == "input"
    assert body[0]["ok"] == "True"
    assert body[1]["detail"] == "not found"


def test_report_tsv_creates_parent_directory(tmp_path):
    nested = tmp_path / "deeply" / "nested" / "report.tsv"
    _write_report_tsv(
        [PrereqRow("tool", "Rscript executable", "/usr/bin/Rscript", True)],
        nested,
    )
    assert nested.exists()


def test_report_tsv_distinguishes_required_and_optional_rows(tmp_path):
    out_path = tmp_path / "report.tsv"
    _write_report_tsv(
        [
            PrereqRow("r-package", "coloc", "R package: coloc", True, required=True),
            PrereqRow("r-package", "susieR", "R package: susieR", False, required=False),
        ],
        out_path,
    )

    with out_path.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    required_flags = {row["name"]: row["required"] for row in rows}
    assert required_flags["coloc"] == "True"
    assert required_flags["susieR"] == "False"
