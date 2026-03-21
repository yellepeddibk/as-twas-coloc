# Poster Outline: Towards Coloc-Confirmed TWAS of Ankylosing Spondylitis

**Conference target**: [PLACEHOLDER — e.g., ASHG Annual Meeting, EULAR, ACR Convergence]

---

## Header

**Title**: Towards Colocalisation-Confirmed Transcriptome-Wide Association Study of Ankylosing Spondylitis to Prioritize Tissue-Specific Drug Targets

**Authors**: [PLACEHOLDER: Author Name, Institution]

**Contact**: [PLACEHOLDER: email / ORCID]

---

## 1. Background / Motivation (Top Left)

- AS: chronic inflammatory arthropathy, ~0.1–0.5% prevalence
- Strong genetic component (h² ~70–90%), HLA-B*27 as major risk factor
- GWAS has identified >100 risk loci — but **which genes?** in **which tissues?**
- TWAS links GWAS signals to gene expression using eQTL prediction models
- Coloc distinguishes shared signals from coincidental proximity

**Key question**: Which genes are causally dysregulated in AS, and in which tissues?

---

## 2. Methods (Top Center / Top Right)

```
AS GWAS summary stats
        │
        ▼
   Harmonization
   (hg38, GTEx varIDs,
    palindromic SNP filter)
        │
        ▼
 S-PrediXcan × [N tissues]
 (GTEx v8 mashr models)
        │
        ▼
  Global BH/FDR correction
  across all gene-tissue pairs
        │
        ▼
  FDR-significant hits
  (q < 0.05)
        │
        ▼
  coloc.abf per locus
  (beta + varbeta, no SNP pre-filter)
        │
        ▼
  PP4 ≥ 0.7 → Coloc hits
  PP3+PP4 ≥ 0.8 → Powered locus
  PP4/(PP3+PP4) ≥ 0.9 → Strong shared signal
```

**Tools**: S-PrediXcan (MetaXcan), coloc R package, GTEx v8, Python pipeline

---

## 3. Results (Center)

> **[PLACEHOLDER — Do NOT fill until real analyses are complete]**

- N gene-tissue pairs tested across N tissues
- N TWAS-significant hits (FDR q < 0.05)
- N coloc-confirmed hits (PP4 ≥ 0.7)

**[INSERT: TWAS QQ-plot]** — `results/figures/qqplot_twas.png`

**[INSERT: PP4 barplot of top genes]** — `results/figures/coloc_summary.png`

---

## 4. Top Coloc-Confirmed Genes (Center Right)

> **[PLACEHOLDER — Insert table after real analyses]**

| Gene | Tissue | PP4 | PP3+PP4 |
|------|--------|-----|---------|
| — | — | — | — |

**[INSERT: top_genes_pp4.tsv content after analyses]**

---

## 5. Conclusions (Bottom Left)

> **[PLACEHOLDER — Do NOT write scientific conclusions before analyses]**

- This pipeline integrates TWAS and COLOC in a reproducible, modular framework
- Enables tissue-specific candidate gene prioritization for AS
- [Real conclusions to follow after analysis completion]

---

## 6. Next Steps (Bottom Center)

- [ ] Download and process real AS GWAS summary statistics
- [ ] Run S-PrediXcan across GTEx v8 tissues
- [ ] Run coloc.abf for TWAS hits
- [ ] Functional annotation of top genes
- [ ] Drug target mapping (OpenTargets, DGIdb)
- [ ] coloc.susie for multi-causal loci

---

## 7. Acknowledgements + QR Code (Bottom Right)

- [PLACEHOLDER: Funding source]
- [PLACEHOLDER: Advisor / PI acknowledgement]
- Data: GTEx Consortium, GWAS Catalog
- Code: https://github.com/yellepeddibk/as-twas-coloc

**[INSERT: QR code linking to GitHub repository]**

---

## Poster Design Notes

- Dimensions: 48" × 36" (landscape) or as required by conference
- Colour scheme: Professional; consider blues/greys for scientific poster
- Keep text concise — use bullet points and diagrams
- Methods flowchart (Section 2) should be a visual box diagram, not text
- All figures should be high-resolution (≥ 300 DPI)
- Poster template: consider using [Inkscape](https://inkscape.org) or conference LaTeX template
