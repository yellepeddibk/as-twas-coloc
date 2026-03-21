"""
test_spredixcan_command_builder.py
-----------------------------------
Tests for the S-PrediXcan command builder.
"""

import pytest
from pathlib import Path

from src.external_tools.spredixcan_wrapper import build_spredixcan_command, run_all_tissues


def make_minimal_config() -> dict:
    return {
        "spredixcan": {
            "script": "external/MetaXcan/software/SPrediXcan.py",
            "snp_column": "varID",
            "effect_allele_column": "effect_allele",
            "non_effect_allele_column": "non_effect_allele",
            "beta_column": "beta",
            "se_column": "se",
            "pvalue_column": "pvalue",
            "output_dir": "results/twas",
            "gwas_file": "data/interim/gwas_harmonized/as_gwas_hg38_varid.tsv.gz",
            "model_db_pattern": "{model_dir}/gtex_v8_mashr_{tissue}.db",
            "covariance_pattern": "{model_dir}/gtex_v8_mashr_{tissue}.txt.gz",
        },
        "gtex": {
            "model_dir": "data/reference/gtex_v8_models",
            "tissues": ["Whole_Blood", "Spleen"],
        },
    }


def test_build_spredixcan_command_returns_list():
    """build_spredixcan_command returns a non-empty list."""
    cfg = make_minimal_config()
    cmd = build_spredixcan_command(
        cfg=cfg,
        tissue="Whole_Blood",
        model_db="data/reference/gtex_v8_models/gtex_v8_mashr_Whole_Blood.db",
        covariance="data/reference/gtex_v8_models/gtex_v8_mashr_Whole_Blood.txt.gz",
        gwas_file="data/interim/gwas_harmonized/as_gwas_hg38_varid.tsv.gz",
        output_file="results/twas/Whole_Blood.spredixcan.csv",
    )
    assert isinstance(cmd, list)
    assert len(cmd) > 0


def test_build_spredixcan_command_contains_required_flags():
    """Command contains the key S-PrediXcan flags."""
    cfg = make_minimal_config()
    cmd = build_spredixcan_command(
        cfg=cfg,
        tissue="Whole_Blood",
        model_db="model.db",
        covariance="cov.txt.gz",
        gwas_file="gwas.tsv.gz",
        output_file="output.csv",
    )
    cmd_str = " ".join(cmd)
    assert "--model_db_path" in cmd_str
    assert "--covariance" in cmd_str
    assert "--gwas_file" in cmd_str
    assert "--snp_column" in cmd_str
    assert "--output_file" in cmd_str
    assert "varID" in cmd_str, "snp_column should be varID (GTEx-style)"


def test_build_spredixcan_command_uses_varid():
    """snp_column must be varID (not rsID) to match GTEx PredictDB models."""
    cfg = make_minimal_config()
    cmd = build_spredixcan_command(
        cfg=cfg,
        tissue="Spleen",
        model_db="model.db",
        covariance="cov.txt.gz",
        gwas_file="gwas.tsv.gz",
        output_file="output.csv",
    )
    # Find --snp_column value
    idx = cmd.index("--snp_column")
    snp_col_value = cmd[idx + 1]
    assert snp_col_value == "varID", (
        "snp_column must be 'varID' to match GTEx PredictDB models. "
        "Using rsID will result in 0% of model SNPs being found."
    )


def test_build_spredixcan_command_effect_allele_columns():
    """Command includes both effect_allele and non_effect_allele columns."""
    cfg = make_minimal_config()
    cmd = build_spredixcan_command(
        cfg=cfg,
        tissue="Whole_Blood",
        model_db="model.db",
        covariance="cov.txt.gz",
        gwas_file="gwas.tsv.gz",
        output_file="output.csv",
    )
    cmd_str = " ".join(cmd)
    assert "--effect_allele_column" in cmd_str
    assert "--non_effect_allele_column" in cmd_str


def test_run_all_tissues_dry_run_returns_paths():
    """run_all_tissues in dry-run mode returns correct output paths without executing."""
    cfg = make_minimal_config()
    result = run_all_tissues(cfg, dry_run=True)
    assert isinstance(result, dict)
    assert "Whole_Blood" in result
    assert "Spleen" in result
    # Output paths should end in .spredixcan.csv
    for tissue, path in result.items():
        assert str(path).endswith(".spredixcan.csv"), (
            f"Expected .spredixcan.csv for {tissue}, got {path}"
        )


def test_build_spredixcan_command_output_path_matches_tissue():
    """Output file path includes the tissue name."""
    cfg = make_minimal_config()
    tissue = "Colon_Sigmoid"
    output_file = f"results/twas/{tissue}.spredixcan.csv"
    cmd = build_spredixcan_command(
        cfg=cfg,
        tissue=tissue,
        model_db="model.db",
        covariance="cov.txt.gz",
        gwas_file="gwas.tsv.gz",
        output_file=output_file,
    )
    cmd_str = " ".join(cmd)
    assert tissue in cmd_str
