# as-twas-coloc

> **Towards Coloc-Confirmed TWAS of Ankylosing Spondylitis to Prioritize Tissue-Specific Drug Targets**
>
> A reproducible bioinformatics research pipeline for honors thesis and manuscript preparation.

---

## Project Overview

This repository implements a modular, reproducible pipeline that:

1. Harmonizes ankylosing spondylitis (AS) GWAS summary statistics to hg38 with GTEx-style varIDs
2. Runs S-PrediXcan (MetaXcan) across multiple GTEx v8 tissues to identify disease-associated genes (TWAS)
3. Applies global Benjamini-Hochberg FDR correction across all gene-tissue pairs
4. Prioritizes TWAS hits using Bayesian colocalization (COLOC) with PP4 ≥ 0.7
5. Reports additional coloc metrics: PP3+PP4 (locus power) and PP4/(PP3+PP4) (signal specificity)
6. Saves provenance metadata, plots, and summary tables for manuscript integration

**Important**: This repository runs in lightweight **mock/demo mode** by default. Real scientific execution requires downloading GTEx v8 reference data and installing external genomics tools (see [Prerequisites](#prerequisites)).

---

## Scientific Motivation

Ankylosing spondylitis (AS) is a chronic inflammatory arthropathy with strong HLA associations and complex non-HLA genetic architecture. While GWAS has identified >100 risk loci, the causal genes and effector tissues remain poorly characterized. Transcriptome-wide association studies (TWAS) using eQTL prediction models (S-PrediXcan) can link GWAS signals to gene expression, but TWAS alone does not distinguish co-localized signals from coincidental proximity. Bayesian colocalization (coloc.abf) tests whether the GWAS and eQTL signals at a locus share a single causal variant, providing stronger evidence for target gene prioritization.

This pipeline integrates TWAS and COLOC to produce a ranked list of candidate genes with tissue-specific biological support — directly informing drug target selection for AS.

---

## Repository Structure

```
as-twas-coloc/
├── README.md
├── requirements.txt
├── .gitignore
├── config/
│   └── as.yaml                   # Pipeline configuration
├── data/
│   ├── raw/                      # Original GWAS inputs (not committed)
│   ├── interim/                  # Intermediate outputs (harmonized GWAS, eQTL regions)
│   ├── processed/                # Final cleaned datasets (TWAS hits, coloc summary)
│   └── reference/                # GTEx v8 models + eQTL data (not committed, large)
├── results/
│   ├── twas/                     # Per-tissue S-PrediXcan CSV outputs
│   ├── coloc/                    # Per-gene coloc JSON outputs
│   ├── tables/                   # Final summary tables (top_genes_pp4.tsv)
│   ├── figures/                  # Plots (QQ, Manhattan-style, COLOC barplot)
│   └── manifests/                # Provenance JSON manifests
├── src/
│   ├── pipeline.py               # Top-level pipeline orchestrator
│   ├── config.py                 # Config loading and access
│   ├── io_utils.py               # File I/O helpers
│   ├── logging_utils.py          # Structured logging
│   ├── schemas.py                # Canonical GWAS schema + column aliases
│   ├── harmonization.py          # GWAS harmonization (liftover, varID, alleles)
│   ├── twas.py                   # TWAS aggregation + FDR correction
│   ├── coloc.py                  # COLOC data preparation + result interpretation
│   ├── reporting.py              # Summary tables
│   ├── provenance.py             # Run manifest generation
│   ├── visualization.py          # Plotting utilities
│   ├── manuscript_support.py     # LaTeX table helpers + artifact checklist
│   └── external_tools/
│       ├── spredixcan_wrapper.py           # S-PrediXcan command builder
│       ├── coloc_r_wrapper.py              # Python → R coloc interface
│       ├── run_coloc.R                     # R coloc.abf script
│       ├── susie_wrapper.py               # coloc.susie placeholder
│       ├── finemap_wrapper.py             # FINEMAP placeholder
│       └── summary_gwas_imputation_wrapper.py  # Imputation placeholder
├── scripts/
│   └── run_pipeline.py           # CLI entry point
├── tests/
│   ├── test_harmonization_contract.py
│   ├── test_varbeta_conversion.py
│   ├── test_spredixcan_command_builder.py
│   ├── test_mock_pipeline_smoke.py
│   └── test_config_loading.py
├── notebooks/
│   └── 01_mock_pipeline_walkthrough.ipynb
├── docs/
│   ├── manuscript/
│   │   ├── manuscript.tex        # LaTeX manuscript scaffold
│   │   └── references.bib        # BibTeX reference database
│   └── poster/
│       └── poster_outline.md     # Poster scaffold for conferences
└── .github/
    └── workflows/
        └── ci.yml                # GitHub Actions CI
```

---

## Setup Instructions

### 1. Clone and install Python dependencies

```bash
git clone https://github.com/yellepeddibk/as-twas-coloc.git
cd as-twas-coloc
pip install -r requirements.txt
```

Python ≥ 3.10 is required. A virtual environment is recommended.

### 2. Mock / Demo Run (no external data required)

```bash
python scripts/run_pipeline.py --mock
```

This generates synthetic GWAS data, runs all pipeline stages using mock data, and writes outputs to:
- `results/figures/` — TWAS QQ-plot, Manhattan-style plot, COLOC barplots
- `results/tables/top_genes_pp4.tsv` — top coloc-confirmed genes
- `results/manifests/` — provenance JSON

### 3. Run Tests

```bash
python -m pytest tests/ -v
```

All tests run in mock mode and do not require external data.

---

## Prerequisites for Real Execution

The following external resources are **NOT bundled** in this repository and must be acquired separately.

### AS GWAS Summary Statistics

Download from the GWAS Catalog or the original publication:
- Recommended: [Ellinghaus et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5657534/) (IGAS meta-analysis)
- Or: [International Genetics of AS Consortium (IGAS)](http://www.igas.org/)
- GWAS Catalog accession: `GCST003452` or newer AS studies
- Place in: `data/raw/as_gwas.tsv.gz`

### GTEx v8 PredictDB Models

Download from [predictdb.org](https://predictdb.org):
- GTEx v8 mashr-based prediction models: `*.db` files (one per tissue)
- GTEx v8 mashr covariance files: `*.txt.gz` files
- Place in: `data/reference/gtex_v8_models/`

### GTEx v8 eQTL Summary Statistics

Required for COLOC regional data extraction.
- Download from [GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex/qtl)
- All tissues single-tissue eQTL summary stats (`GTEx_Analysis_v8_eQTL.tar.gz`)
- Place in: `data/reference/gtex_v8_eqtl/`

### MetaXcan / S-PrediXcan

```bash
git clone https://github.com/hakyimlab/MetaXcan.git external/MetaXcan
pip install -r external/MetaXcan/software/requirements.txt
```

Update `config/as.yaml` → `spredixcan.script` to point to the installed script.

### R + coloc

```r
install.packages(c("coloc", "jsonlite"))
# Optional for multi-signal loci:
install.packages("susieR")
```

Update `config/as.yaml` → `coloc.r_script` as needed.

---

## Mock vs Real Execution

| Feature | Mock mode | Real mode |
|---------|-----------|-----------|
| GWAS data | Synthetic (200 SNPs) | Real AS GWAS summary stats |
| S-PrediXcan | Synthetic outputs | Requires MetaXcan + GTEx models |
| COLOC | Synthetic posteriors | Requires R + coloc + GTEx eQTL |
| Liftover | Placeholder warning | Requires pyliftover / CrossMap |
| varID | Constructed from mock positions | Matched to real GTEx varIDs |
| Figures | Generated from mock data | Generated from real results |
| CI | ✓ Always passes | Not run in CI |

---

## Harmonization Pitfalls

1. **hg19 vs hg38 mismatch**: Most published AS GWAS use hg19. GTEx v8 uses hg38. A mismatch causes "0% of model SNPs found" in S-PrediXcan. The pipeline warns explicitly and calls a liftover placeholder.

2. **GTEx varID vs rsID**: GTEx PredictDB models use `chr{chrom}_{pos}_{ref}_{alt}_b38` as the primary key, not rsID. The pipeline constructs varIDs from the harmonized positions.

3. **Allele alignment and strand flips**: The harmonization module flags strand-ambiguous situations. Full correction requires a reference panel.

4. **Palindromic SNPs (AT/CG)**: Cannot be reliably strand-flipped without allele frequencies. Dropped by default (`gwas.drop_palindromic_snps: true`). Set `af_disambiguation_available: true` to retain them (AF-based logic must be implemented).

5. **Pre-filtering SNPs before COLOC**: Do NOT filter SNPs by p-value before running coloc.abf. This inflates PP4. The pipeline preserves full locus coverage.

6. **SE vs varbeta**: coloc.abf requires `varbeta = se**2`, not `se` directly. The pipeline computes this explicitly.

---

## COLOC Interpretation Metrics

| Metric | Formula | Threshold | Interpretation |
|--------|---------|-----------|----------------|
| PP4 | — | ≥ 0.7 | Primary coloc hit: shared causal variant |
| PP3+PP4 | PP3 + PP4 | ≥ 0.8 | Powered locus: enough signal to distinguish H3/H4 |
| PP4 ratio | PP4 / (PP3+PP4) | ≥ 0.9 | Strong preference for shared signal over distinct signals |

**Coloc hypotheses:**
- H0: No association with either trait
- H1: Association with GWAS only
- H2: Association with eQTL only
- H3: Association with both, distinct causal variants
- H4: Association with both, **shared** causal variant (coloc)

---

## Pipeline Outputs → Manuscript Sections

| Output file | Manuscript use |
|------------|----------------|
| `results/tables/top_genes_pp4.tsv` | Table 2: Coloc-prioritized genes |
| `results/figures/qqplot_twas.png` | Figure 2: TWAS QQ-plot |
| `results/figures/manhattan_like_twas.png` | Figure 3: TWAS association overview |
| `results/figures/coloc_summary.png` | Figure 4: COLOC PP4 barplot |
| `data/processed/twas_hits/twas_hits_fdr05.tsv.gz` | Supplementary Table S1 |
| `data/processed/coloc_results/coloc_summary.tsv.gz` | Supplementary Table S2 |
| `results/manifests/manifest_*.json` | Methods: reproducibility statement |

The manuscript scaffold lives at `docs/manuscript/manuscript.tex`. Update the Results section placeholders after completing real analyses.

---

## Reproducibility Notes

- All random seeds are fixed for mock data generation
- Each run saves a provenance manifest with timestamp, git commit hash, config snapshot, coloc priors, and input file checksums
- Pipeline stages are independently rerunnable from their input files
- Output file naming is deterministic

---

## Limitations and TODOs

- [ ] Implement real hg19→hg38 liftover (pyliftover or CrossMap)
- [ ] Implement real eQTL regional data extraction from GTEx eQTL files
- [ ] Implement real coloc.abf R execution with GTEx eQTL data
- [ ] Implement coloc.susie for multi-causal loci (after initial coloc.abf screen)
- [ ] Implement FINEMAP for top priority loci
- [ ] Validate varID construction against hg38 reference FASTA
- [ ] Implement AF-based palindromic SNP disambiguation
- [ ] Add sensitivity analysis plots for top coloc loci
- [ ] Run summary-gwas-imputation for GWAS panels with incomplete GTEx overlap

---

## Citation

If this pipeline scaffold contributes to your research, please cite the underlying methods:

- S-PrediXcan: Barbeira et al. (2018) *Nature Communications*
- COLOC: Giambartolomei et al. (2014) *PLOS Genetics*; Wallace (2021) *PLOS Genetics*
- GTEx v8: GTEx Consortium (2020) *Science*
- AS GWAS: Ellinghaus et al. (2016) *Science* (or the specific GWAS used)

