"""Tests for validate_publication_counts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from scripts.validate_publication_counts import (
    CountCheck,
    collect_count_checks,
    write_report_tsv,
)


def _write_twas_csv(out_dir: Path, tissue: str, p_values: list[float | None]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        {
            "gene": [f"ENSG{i:09d}" for i in range(len(p_values))],
            "gene_name": [f"G{i}" for i in range(len(p_values))],
            "tissue": tissue,
            "zscore": [0.5] * len(p_values),
            "pvalue": p_values,
        }
    )
    df.to_csv(out_dir / f"{tissue}.spredixcan.csv", index=False)


def _write_harm_qc(path: Path, **kwargs: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([kwargs]).to_csv(path, sep="\t", index=False)


def _write_twas_hits(path: Path, n_hits: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        {
            "gene": [f"ENSG{i:09d}" for i in range(n_hits)],
            "tissue": ["Whole_Blood"] * n_hits,
            "pvalue": [1e-6] * n_hits,
        }
    )
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def _write_coloc_summary(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression="gzip")


@pytest.fixture
def synthetic_artifacts(tmp_path: Path) -> dict:
    twas_dir = tmp_path / "results" / "twas"
    _write_twas_csv(twas_dir, "Whole_Blood", [1e-9, 1e-3, None])
    _write_twas_csv(twas_dir, "Spleen", [1e-2, 0.5])

    twas_hits = tmp_path / "data" / "processed" / "twas_hits" / "twas_hits_fdr05.tsv.gz"
    _write_twas_hits(twas_hits, n_hits=2)

    coloc_summary = tmp_path / "data" / "processed" / "coloc_results" / "coloc_summary.tsv.gz"
    _write_coloc_summary(
        coloc_summary,
        [
            {"gene_id": "ENSG000001", "tissue": "Whole_Blood", "PP4": 0.95},
            {"gene_id": "ENSG000002", "tissue": "Spleen",      "PP4": 0.80},
            {"gene_id": "ENSG000001", "tissue": "Spleen",      "PP4": 0.10},
        ],
    )

    harm_qc = tmp_path / "data" / "processed" / "qc" / "gwas_harmonization_qc.tsv"
    _write_harm_qc(
        harm_qc,
        input_snp_count=100,
        lifted_snp_count=99,
        failed_liftover_count=1,
        palindromic_dropped_count=10,
        final_retained_snp_count=89,
    )

    return {
        "twas_dir": twas_dir,
        "twas_hits": twas_hits,
        "coloc_summary": coloc_summary,
        "harmonization_qc": harm_qc,
    }


def test_all_counts_match(synthetic_artifacts):
    counts = {
        "gwas_harmonization": {
            "input_snp_count": 100,
            "final_retained_snp_count": 89,
        },
        "twas": {
            "raw_gene_tissue_rows": 5,
            "valid_gene_tissue_tests": 4,
            "fdr_significant_hits": 2,
            "tissues": ["Whole_Blood", "Spleen"],
        },
        "coloc": {
            "pp4_threshold": 0.7,
            "coloc_confirmed_pairs": 2,
            "unique_genes": 2,
            "whole_blood_pairs": 1,
        },
    }

    rows = collect_count_checks(
        counts,
        twas_dir=synthetic_artifacts["twas_dir"],
        twas_hits_path=synthetic_artifacts["twas_hits"],
        coloc_summary_path=synthetic_artifacts["coloc_summary"],
        harmonization_qc_path=synthetic_artifacts["harmonization_qc"],
    )

    assert all(row.ok for row in rows), [
        (row.category, row.name, row.expected, row.actual) for row in rows if not row.ok
    ]


def test_drift_is_detected(synthetic_artifacts):
    counts = {
        "twas": {"fdr_significant_hits": 999},
        "coloc": {"pp4_threshold": 0.7, "coloc_confirmed_pairs": 999},
    }

    rows = collect_count_checks(
        counts,
        twas_dir=synthetic_artifacts["twas_dir"],
        twas_hits_path=synthetic_artifacts["twas_hits"],
        coloc_summary_path=synthetic_artifacts["coloc_summary"],
        harmonization_qc_path=synthetic_artifacts["harmonization_qc"],
    )

    drift = [row for row in rows if not row.ok]
    drift_names = {row.name for row in drift}
    assert "fdr_significant_hits" in drift_names
    assert "coloc_confirmed_pairs" in drift_names


def test_report_tsv_writes_machine_readable_rows(tmp_path):
    out_path = tmp_path / "reports" / "publication_counts_report.tsv"
    write_report_tsv(
        [CountCheck("twas", "fdr_significant_hits", 248, 248, True)],
        out_path,
    )

    text = out_path.read_text()
    assert "category\tname\texpected\tactual\tok" in text
    assert "twas\tfdr_significant_hits\t248\t248\tTrue" in text


def test_collects_only_categories_present_in_counts(synthetic_artifacts):
    counts = {"twas": {"fdr_significant_hits": 2}}

    rows = collect_count_checks(
        counts,
        twas_dir=synthetic_artifacts["twas_dir"],
        twas_hits_path=synthetic_artifacts["twas_hits"],
        coloc_summary_path=synthetic_artifacts["coloc_summary"],
        harmonization_qc_path=synthetic_artifacts["harmonization_qc"],
    )

    categories = {row.category for row in rows}
    assert categories == {"twas"}
