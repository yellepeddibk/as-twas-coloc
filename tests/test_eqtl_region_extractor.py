"""Tests for gene-level eQTL region extractor."""

from pathlib import Path

import pandas as pd

from src.eqtl_regions import (
    build_eqtl_regions_for_coloc,
    load_gene_tss_from_gtf,
    normalize_varid,
)


def test_normalize_varid_removes_chr_prefix():
    assert normalize_varid("chr1_123_A_G_b38") == "1_123_A_G_b38"
    assert normalize_varid("1_123_A_G_b38") == "1_123_A_G_b38"


def test_build_eqtl_regions_for_coloc_writes_expected_file(tmp_path: Path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG000001"; gene_name "GENE1";',
                '',
            ]
        ),
        encoding="utf-8",
    )

    gene_tss = load_gene_tss_from_gtf(gtf)

    twas_hits = pd.DataFrame(
        {
            "gene": ["ENSG000001"],
            "tissue": ["Whole_Blood"],
        }
    )
    gwas = pd.DataFrame(
        {
            "varID": ["1_150_A_G_b38", "1_300_C_T_b38"],
        }
    )

    eqtl_dir = tmp_path / "eqtl"
    eqtl_dir.mkdir(parents=True, exist_ok=True)
    eqtl_file = eqtl_dir / "Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"
    pd.DataFrame(
        {
            "variant_id": ["chr1_150_A_G_b38", "chr1_700_A_G_b38"],
            "slope": [0.2, 0.5],
            "slope_se": [0.1, 0.2],
            "pval_nominal": [1e-5, 0.2],
            "ma_samples": [300, 300],
        }
    ).to_csv(eqtl_file, sep="\t", index=False, compression="gzip")

    out_dir = tmp_path / "out"
    summary = build_eqtl_regions_for_coloc(
        twas_hits_df=twas_hits,
        harmonized_gwas_df=gwas,
        gene_tss_df=gene_tss,
        eqtl_dir=eqtl_dir,
        output_dir=out_dir,
        window_bp=1000,
        chunk_size=100,
    )

    written = summary[summary["status"] == "written"]
    assert len(written) == 1

    out_file = out_dir / "Whole_Blood" / "ENSG000001.tsv.gz"
    assert out_file.exists()

    out_df = pd.read_csv(out_file, sep="\t")
    assert set(["varID", "beta", "se", "pvalue"]).issubset(out_df.columns)
    assert "1_150_A_G_b38" in set(out_df["varID"].tolist())
