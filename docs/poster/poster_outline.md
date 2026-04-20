# Poster Outline: Colocalization-Confirmed TWAS of Ankylosing Spondylitis

**Event**: ISU Honors Poster Symposium, May 6, 2026

**Status**: COMPLETE - all content finalized with real results

---

## Header

**Title**: Colocalization-Confirmed Transcriptome-Wide Association Study of Ankylosing Spondylitis to Prioritize Tissue-Specific Drug Targets

**Author**: Bhargav Yellepeddi

**Advisor**: Dr. Ratul Chowdhury

**Department**: Chemical and Biological Engineering, Iowa State University

---

## 1. Introduction (Column 1)

- AS: chronic inflammatory arthritis of spine/sacroiliac joints, 0.1-0.5% prevalence
- Strong genetic component (h^2 ~70-90%), HLA-B*27 as major risk factor
- GWAS identified >100 risk loci, but most are non-coding
- TWAS integrates GWAS with gene expression models to identify causal genes
- Colocalization distinguishes shared causal signals from LD artifacts

**Key question**: Which genes have genetically regulated expression that causally colocalizes with AS risk variants, and in which tissues?

---

## 2. Methods (Column 1)

Pipeline (5 stages):
1. GWAS QC and Harmonization: AS GWAS GCST005529 (N=25,764), liftOver hg38, palindromic SNP removal
2. S-PrediXcan TWAS: GTEx v8 Mashr models, 8 tissues, 2,316 gene-tissue tests
3. FDR Correction: Benjamini-Hochberg, alpha=0.05
4. Bayesian Colocalization: coloc.abf, PP4 >= 0.70 threshold
5. Prioritized gene list output

**GWAS QC numbers**: 123,941 input SNPs -> 123,925 lifted -> 14,032 palindromic removed -> 109,893 harmonized

**8 Tissues**: Whole Blood, Spleen, Lung, Colon Sigmoid, Colon Transverse, Small Intestine, Muscle Skeletal, Adipose Subcutaneous

---

## 3. TWAS Results (Column 2)

- 2,316 gene-tissue pairs tested across 8 GTEx v8 tissues
- 248 significant associations after BH-FDR correction (alpha=0.05)
- Strong departure from null in QQ-plot, consistent with polygenic architecture

**Figures**: Manhattan-style plot (`docs/poster/figures/manhattan_twas.png`), QQ-plot (`docs/poster/figures/qqplot_twas.png`)

---

## 4. Colocalization Results (Column 2-3)

- 248 TWAS hits tested with coloc.abf
- **35 gene-tissue pairs (18 unique genes)** with PP4 >= 0.70
- Confirms shared causal variant between GWAS and eQTL

**Top colocalized genes** (by PP4):

| Gene | Top Tissue | PP4 | z-score | FDR p |
|------|-----------|-----|---------|-------|
| ERAP1 | Whole Blood | 0.999 | 13.28 | 2.5e-38 |
| GPR25 | Whole Blood | 0.997 | -7.03 | 6.0e-11 |
| NPIPB7 | Whole Blood | 0.993 | -5.16 | 5.9e-06 |
| CCDC116 | Spleen | 0.988 | -5.16 | 5.9e-06 |
| UBE2L3 | Whole Blood | 0.983 | 5.32 | 2.5e-06 |
| CARD9 | Whole Blood | 0.974 | 4.94 | 1.7e-05 |
| FAS | Colon Sigmoid | 0.973 | 4.42 | 1.9e-04 |
| MYH7B | Colon Transverse | 0.973 | 4.42 | 1.9e-04 |
| RNF123 | Adipose Subcut. | 0.963 | 4.47 | 1.5e-04 |
| PRSS53 | Whole Blood | 0.976 | -3.68 | 3.3e-03 |

**Figures**: PP4 bar plot (`docs/poster/figures/coloc_pp4_barplot.png`), PP4 vs PP3+PP4 scatter (`docs/poster/figures/coloc_pp4_vs_pp3pp4.png`)

---

## 5. Discussion (Column 3)

- ERAP1 (aminopeptidase, MHC-I antigen processing) validates pipeline as established AS gene
- ERAP1 colocalizes in both Whole Blood and Small Intestine, highlighting tissue-specific regulation
- Novel candidates: GPR25, CCDC116, NPIPB7 warrant experimental follow-up
- CARD9 and FAS implicate innate immune and apoptotic pathways beyond HLA-B27
- Colocalization reduced 248 TWAS hits to 35 confirmed pairs, filtering LD-driven false positives

---

## 6. Conclusion (Column 3)

Integrating S-PrediXcan TWAS with Bayesian colocalization across 8 tissues identified 18 unique genes with strong evidence of shared causal variants with AS risk loci. This two-stage approach provides a principled filter against LD contamination and prioritizes tissue-specific regulatory targets for therapeutic development.

---

## 7. Future Directions (Column 3)

- Fine-mapping with SuSiE to identify credible causal variants
- Validation in independent AS cohorts
- Drug repurposing analysis for top gene targets
- Extension to additional immune-mediated diseases

---

## 8. Acknowledgements (Column 3)

- Dr. Ratul Chowdhury (advisor and mentorship)
- GTEx Consortium (v8 expression models)
- NHGRI-EBI GWAS Catalog (GCST005529)
- S-PrediXcan and coloc.abf software authors
- QR code links to: https://github.com/yellepeddibk/as-twas-coloc

---

## Build Instructions

```bash
# Generate poster figures (300 DPI, gene symbols)
python scripts/build_poster_figures.py

# Build poster PDF
cd docs/poster && make poster
```

**Format**: 48" x 36" landscape, LaTeX beamerposter with ISU brand theme
**Files**: `poster.tex`, `beamerthemeISU.sty`, `Makefile`, `figures/`
