# Supplementary Materials for Zhu et al.
**Manuscript ID:** EJSOBI-D-26-00507

This repository contains the supplementary data and code for:  
**DNA-based nematode community profiles reflect composition by taxon mass rather than abundance**

## Repository structure

| Folder | Description |
|--------|-------------|
| `1 Raw_Data/` | Raw ASV tables, morphological measurements, and count data |
| `2 Processed_Data/` | Processed community data for statistical analyses |
| `3 Code/` | R scripts for data processing and statistical analyses |

Detailed descriptions of each file are provided in the `README.txt` file within each folder.

---

## Software requirements

All R scripts were run using R version 4.5.1 with the following packages:

- `vegan` (v2.6-8)
- `tidyverse` (v2.0.0)
- `glmmTMB` (v1.1.10)
- `DHARMa` (v0.4.7)

Linear regressions were performed in SPSS Statistics 22. Figures were finalized in OriginPro 2024.

---

## Running the code

Scripts are numbered in the order of execution:

1. Rarefaction.R` – normalize sequencing depth across samples
2. Relative composition.R` – calculate relative read abundance
3. adjust_p.R` – FDR correction for multiple comparisons
4. LMMs.R` – beta mixed‑effects models
5. Bland-Altman analysis.R` – robust Bland‑Altman analysis

All scripts expect input files to be located in the `../2 Processed_Data/` directory.

---

## Data availability

Raw sequencing reads have been deposited in the NCBI Sequence Read Archive under accession number **PRJNA1394764**.
