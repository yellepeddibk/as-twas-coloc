#!/usr/bin/env python3
"""
validate_publication_counts.py
-------------------------------
Validate that the live result artifacts still match the frozen publication
counts recorded in ``config/publication_counts.json``.

This is a publication-readiness gate: it prevents silent drift between the
manuscript narrative and the actual files in the repository.

Usage
-----
  python scripts/validate_publication_counts.py \
      --config config/as.yaml \
      --base-dir . \
      --counts config/publication_counts.json \
      --report-tsv results/manifests/publication_counts_report.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd

_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_REPO_ROOT))


@dataclass(frozen=True)
class CountCheck:
    category: str
    name: str
    expected: Any
    actual: Any
    ok: bool


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Validate frozen publication counts against live artifacts")
    p.add_argument("--base-dir", default=".", help="Repository root")
    p.add_argument(
        "--counts",
        default="config/publication_counts.json",
        help="Path to frozen publication counts JSON",
    )
    p.add_argument(
        "--twas-dir",
        default="results/twas",
        help="Directory containing per-tissue S-PrediXcan CSV outputs",
    )
    p.add_argument(
        "--twas-hits",
        default="data/processed/twas_hits/twas_hits_fdr05.tsv.gz",
        help="FDR-significant TWAS hits TSV path",
    )
    p.add_argument(
        "--coloc-summary",
        default="data/processed/coloc_results/coloc_summary.tsv.gz",
        help="Aggregated coloc summary TSV path",
    )
    p.add_argument(
        "--harmonization-qc",
        default="data/processed/qc/gwas_harmonization_qc.tsv",
        help="GWAS harmonization QC TSV path",
    )
    p.add_argument(
        "--report-tsv",
        help="Optional path for a machine-readable report TSV",
    )
    return p.parse_args()


def _resolve(base_dir: Path, value: str | Path) -> Path:
    p = Path(value)
    return p if p.is_absolute() else base_dir / p


def _load_counts(counts_path: Path) -> dict:
    with counts_path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _aggregate_twas(twas_dir: Path) -> pd.DataFrame:
    csvs = sorted(twas_dir.glob("*.spredixcan.csv"))
    if not csvs:
        raise FileNotFoundError(f"No S-PrediXcan CSVs found in {twas_dir}")
    frames = []
    for csv_path in csvs:
        df = pd.read_csv(csv_path)
        df["_tissue"] = csv_path.stem.replace(".spredixcan", "")
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def _harmonization_check(qc_path: Path, expected: dict) -> Iterable[CountCheck]:
    qc = pd.read_csv(qc_path, sep="\t")
    if qc.empty:
        raise ValueError(f"Harmonization QC file is empty: {qc_path}")
    row = qc.iloc[0]
    for key in (
        "input_snp_count",
        "lifted_snp_count",
        "palindromic_dropped_count",
        "final_retained_snp_count",
    ):
        if key in expected:
            actual = int(row[key])
            yield CountCheck("gwas_harmonization", key, expected[key], actual, actual == expected[key])


def _twas_checks(twas_dir: Path, twas_hits_path: Path, expected: dict) -> Iterable[CountCheck]:
    twas = _aggregate_twas(twas_dir)
    twas_hits = pd.read_csv(twas_hits_path, sep="\t")

    if "raw_gene_tissue_rows" in expected:
        actual = int(len(twas))
        yield CountCheck("twas", "raw_gene_tissue_rows", expected["raw_gene_tissue_rows"], actual, actual == expected["raw_gene_tissue_rows"])

    if "valid_gene_tissue_tests" in expected:
        actual = int(twas["pvalue"].notna().sum())
        yield CountCheck("twas", "valid_gene_tissue_tests", expected["valid_gene_tissue_tests"], actual, actual == expected["valid_gene_tissue_tests"])

    if "fdr_significant_hits" in expected:
        actual = int(len(twas_hits))
        yield CountCheck("twas", "fdr_significant_hits", expected["fdr_significant_hits"], actual, actual == expected["fdr_significant_hits"])

    if "tissues" in expected:
        actual_tissues = sorted(twas["_tissue"].unique().tolist())
        expected_tissues = sorted(list(expected["tissues"]))
        yield CountCheck("twas", "tissues", expected_tissues, actual_tissues, actual_tissues == expected_tissues)


def _coloc_checks(coloc_summary_path: Path, expected: dict) -> Iterable[CountCheck]:
    coloc = pd.read_csv(coloc_summary_path, sep="\t")
    threshold = float(expected.get("pp4_threshold", 0.7))
    pp4_hits = coloc[coloc["PP4"] >= threshold]

    if "coloc_confirmed_pairs" in expected:
        actual = int(len(pp4_hits))
        yield CountCheck("coloc", "coloc_confirmed_pairs", expected["coloc_confirmed_pairs"], actual, actual == expected["coloc_confirmed_pairs"])

    if "unique_genes" in expected:
        actual = int(pp4_hits["gene_id"].nunique())
        yield CountCheck("coloc", "unique_genes", expected["unique_genes"], actual, actual == expected["unique_genes"])

    if "whole_blood_pairs" in expected:
        actual = int((pp4_hits["tissue"] == "Whole_Blood").sum())
        yield CountCheck("coloc", "whole_blood_pairs", expected["whole_blood_pairs"], actual, actual == expected["whole_blood_pairs"])


def collect_count_checks(
    counts: dict,
    twas_dir: Path,
    twas_hits_path: Path,
    coloc_summary_path: Path,
    harmonization_qc_path: Path,
) -> list[CountCheck]:
    """Build the full list of count checks. Importable for tests."""
    rows: list[CountCheck] = []

    harm_expected = counts.get("gwas_harmonization") or {}
    if harm_expected and harmonization_qc_path.exists():
        rows.extend(_harmonization_check(harmonization_qc_path, harm_expected))

    twas_expected = counts.get("twas") or {}
    if twas_expected:
        rows.extend(_twas_checks(twas_dir, twas_hits_path, twas_expected))

    coloc_expected = counts.get("coloc") or {}
    if coloc_expected:
        rows.extend(_coloc_checks(coloc_summary_path, coloc_expected))

    return rows


def write_report_tsv(rows: list[CountCheck], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["category", "name", "expected", "actual", "ok"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "category": row.category,
                "name": row.name,
                "expected": row.expected,
                "actual": row.actual,
                "ok": row.ok,
            })


def _ok(flag: bool) -> str:
    return "OK" if flag else "DRIFT"


def main() -> int:
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()

    counts_path = _resolve(base_dir, args.counts)
    twas_dir = _resolve(base_dir, args.twas_dir)
    twas_hits_path = _resolve(base_dir, args.twas_hits)
    coloc_summary_path = _resolve(base_dir, args.coloc_summary)
    harmonization_qc_path = _resolve(base_dir, args.harmonization_qc)

    counts = _load_counts(counts_path)
    rows = collect_count_checks(
        counts,
        twas_dir,
        twas_hits_path,
        coloc_summary_path,
        harmonization_qc_path,
    )

    if args.report_tsv:
        report_path = _resolve(base_dir, args.report_tsv)
        write_report_tsv(rows, report_path)

    print("Publication count validation")
    print("=" * 40)
    drift = [r for r in rows if not r.ok]
    for row in rows:
        print(f"[{_ok(row.ok):5}] {row.category} / {row.name}: expected={row.expected} actual={row.actual}")
    print("=" * 40)

    if args.report_tsv:
        print(f"Report: {report_path}")

    if drift:
        print(f"Drift detected in {len(drift)} count(s)")
        return 1

    print("All publication counts match frozen artifacts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
