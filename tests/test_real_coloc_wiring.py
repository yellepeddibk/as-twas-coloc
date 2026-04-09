"""
Test real-mode COLOC wiring in pipeline orchestrator.
"""

from pathlib import Path

import pandas as pd
import pytest

from src.pipeline import _run_real_coloc


def _minimal_cfg(eqtl_pattern: str) -> dict:
    return {
        "gwas": {
            "n_total": 1000,
            "n_cases": 450,
            "n_controls": 550,
            "case_fraction": 0.45,
        },
        "coloc": {
            "eqtl_region_pattern": eqtl_pattern,
            "r_script": "src/external_tools/run_coloc.R",
            "pp4_threshold": 0.7,
            "pp3_pp4_threshold": 0.8,
            "pp4_ratio_threshold": 0.9,
        },
    }


def test_run_real_coloc_writes_json(tmp_path, monkeypatch):
    harmonized_gwas = pd.DataFrame(
        {
            "varID": ["1_100_A_G_b38", "1_200_C_T_b38"],
            "beta": [0.12, -0.05],
            "se": [0.04, 0.03],
            "pvalue": [1e-4, 2e-3],
            "chrom": ["1", "1"],
            "pos": [100, 200],
            "effect_allele": ["A", "C"],
            "non_effect_allele": ["G", "T"],
            "N": [1000, 1000],
        }
    )
    twas_hits = pd.DataFrame(
        {
            "gene": ["ENSG000001"],
            "tissue": ["Whole_Blood"],
        }
    )

    eqtl_dir = tmp_path / "data" / "interim" / "eqtl_regions" / "Whole_Blood"
    eqtl_dir.mkdir(parents=True, exist_ok=True)
    eqtl_path = eqtl_dir / "ENSG000001.tsv.gz"
    pd.DataFrame(
        {
            "varID": ["1_100_A_G_b38", "1_200_C_T_b38"],
            "beta": [0.2, -0.1],
            "se": [0.05, 0.04],
            "pvalue": [5e-6, 7e-4],
            "N": [350, 350],
        }
    ).to_csv(eqtl_path, sep="\t", index=False, compression="gzip")

    def _fake_run_coloc_abf(dataset1, dataset2, priors, r_script=None, mock=False):
        assert len(dataset1["snp"]) == 2
        assert len(dataset2["snp"]) == 2
        return {
            "PP.H0.abf": 0.01,
            "PP.H1.abf": 0.02,
            "PP.H2.abf": 0.03,
            "PP.H3.abf": 0.14,
            "PP.H4.abf": 0.80,
            "nsnps": 2,
        }

    monkeypatch.setattr("src.pipeline.run_coloc_abf", _fake_run_coloc_abf)

    coloc_out_dir = tmp_path / "results" / "coloc"
    cfg = _minimal_cfg("data/interim/eqtl_regions/{tissue}/{gene}.tsv.gz")

    _run_real_coloc(
        twas_hits=twas_hits,
        harmonized_gwas=harmonized_gwas,
        cfg=cfg,
        coloc_priors={"p1": 1e-4, "p2": 1e-4, "p12": 1e-5},
        coloc_out_dir=coloc_out_dir,
        base_dir=tmp_path,
    )

    out_json = coloc_out_dir / "Whole_Blood" / "ENSG000001.json"
    assert out_json.exists(), f"Expected coloc JSON output: {out_json}"


def test_run_real_coloc_raises_when_all_eqtl_missing(tmp_path):
    harmonized_gwas = pd.DataFrame(
        {
            "varID": ["1_100_A_G_b38"],
            "beta": [0.12],
            "se": [0.04],
            "pvalue": [1e-4],
            "chrom": ["1"],
            "pos": [100],
            "effect_allele": ["A"],
            "non_effect_allele": ["G"],
            "N": [1000],
        }
    )
    twas_hits = pd.DataFrame(
        {
            "gene": ["ENSG_MISSING"],
            "tissue": ["Whole_Blood"],
        }
    )

    coloc_out_dir = tmp_path / "results" / "coloc"
    cfg = _minimal_cfg("data/interim/eqtl_regions/{tissue}/{gene}.tsv.gz")

    with pytest.raises(RuntimeError, match="Real COLOC could not run"):
        _run_real_coloc(
            twas_hits=twas_hits,
            harmonized_gwas=harmonized_gwas,
            cfg=cfg,
            coloc_priors={"p1": 1e-4, "p2": 1e-4, "p12": 1e-5},
            coloc_out_dir=coloc_out_dir,
            base_dir=tmp_path,
        )
