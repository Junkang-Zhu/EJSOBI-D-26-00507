Supplementary Materials for:

Mass proportion as a new perspective for interpreting DNA‑based nematode community profiles

Zhu et al.

================================================================================
Quick reference
================================================================================

| File / Folder | Description |
|---------------|-------------|
| 1 Raw_Data/ | Raw sequencing and morphological data |
| 2 Processed_Data/ | Processed community data for statistical analyses |
| 3 Code/ | R scripts for all analyses |
| | |
| 1 ASV PR2 and SILVA Raw.xlsx | ASV tables (PR2 + SILVA 138 annotations) |
| 2 ASV PR2 AllEukaryotes.xlsx | Raw ASV tables (PR2, all eukaryotes) |
| 3 ASV PR2 Nematoda.xlsx | Nematode‑only ASV tables (PR2) |
| 4 Nematode body size raw data.xlsx | Individual body length, width, mass (18 samples) |
| 5 Nematode count abundance raw data.xlsx | Morphological counts + moisture + dry soil abundance |
| Processed data.xlsx | Morphology counts, biomass, average mass, eDNA, nematode‑DNA |
| 1_Rarefaction.R | Rarefaction script |
| 2_Relative composition.R | Relative abundance calculation |
| 3_adjust_p.R | FDR correction |
| 4_LMMs.R | Beta mixed‑effects models |
| 5_Bland‑Altman.R | Bland‑Altman analysis |

================================================================================
Folder: 1 Raw_Data
================================================================================

This folder contains all raw data files used in this study, including ASV tables
from two reference databases, morphological measurements, and nematode counts.

-------------------------------------------------------------------------------
File list and descriptions
-------------------------------------------------------------------------------

1. 1 ASV PR2 and SILVA Raw.xlsx

   ASV tables annotated against the PR2 and SILVA 138 reference databases.
   This file contains two sheets:
     - Sheet 1: PR2 annotation results
     - Sheet 2: SILVA 138 annotation results
   Each sheet includes ASV IDs, taxonomic assignments (from kingdom to genus/
   species where available), and read counts per sample.

2. 2 ASV PR2 AllEukaryotes.xlsx

   Raw ASV tables annotated against the PR2 database, separately for the two
   DNA sources. This file contains two sheets:
     - Sheet 1: eDNA samples
     - Sheet 2: Nematode‑extracted DNA samples
   Each sheet includes all ASVs assigned to all eukaryotic taxa, with read
   counts per sample.

3. 3 ASV PR2 Nematoda.xlsx

   Subset of the PR2‑annotated ASV tables retaining only sequences assigned to
   Nematoda. This file contains two sheets:
     - Sheet 1: eDNA samples (nematode ASVs only)
     - Sheet 2: Nematode‑extracted DNA samples (nematode ASVs only)
   These filtered tables were used as the starting point for downstream
   community analyses after removal of non‑target taxa.

4. 4 Nematode body size raw data.xlsx

   Individual‑level morphological measurements for all identified nematodes.
   This file contains 18 sheets, one for each sample (Z1–Z18). Each sheet
   includes:
     - Sample ID
     - Genus/Family assignment
     - Body length (L, in µm)
     - Maximum body diameter (D, in µm)
     - Calculated individual fresh mass (W, in µg), based on Andrássy (1956):
       W = L × D² / (1.6 × 10⁶)

5. 5 Nematode count abundance raw data.xlsx

   Morphological count data and soil moisture content. This file contains:
     - Nematode counts per taxon per sample
     - Soil moisture content (%) determined by the aluminum box method
     - Converted nematode abundance expressed as individuals per 100 g dry soil

-------------------------------------------------------------------------------
Notes
-------------------------------------------------------------------------------

- All raw data are provided as originally recorded or directly exported from
  sequencing and morphological analyses, without any preprocessing or
  aggregation, to ensure full transparency and reproducibility.
- Family‑level aggregation and calculation of relative abundances were
  performed separately in R and Excel; the processed data files are provided in
  the 02_Processed_Data folder.
- Taxonomic assignments were based on the PR2 database (version 5.2) and
  SILVA 138. The PR2 annotations were used for all downstream analyses due to
  its higher resolution for nematode taxa.


================================================================================
Folder: 2 Processed Data
================================================================================

This folder contains processed community data used as input for statistical
analyses and figure generation. All data have been aggregated at the genus
level and organized into a single Excel file with multiple sheets.

-------------------------------------------------------------------------------
File list and descriptions
-------------------------------------------------------------------------------

1. Processed data.xlsx

   Consolidated dataset containing processed morphological and DNA-based
   community data. This file contains five sheets:

   ------------------------------------------------------------
   Sheet 1: Morphological identification
   ------------------------------------------------------------
   Genus‑level nematode counts per sample based on morphological
   identification. Values represent the number of individuals
   of each genus identified in the subsample (≥250 individuals
   per sample, or all individuals when fewer than 250 were
   present).

   ------------------------------------------------------------
   Sheet 2: Total biomass (μg)
   ------------------------------------------------------------
   Sum of individual fresh mass (μg) for each genus per sample,
   calculated from body length and width measurements using the
   Andrássy (1956) formula.

   ------------------------------------------------------------
   Sheet 3: Average mass (μg)
   ------------------------------------------------------------
   Mean individual fresh mass (μg) for each genus per sample,
   calculated as total biomass divided by the number of
   individuals identified for that genus.

   ------------------------------------------------------------
   Sheet 4: eDNA metabarcoding
   ------------------------------------------------------------
   Genus‑level relative read abundance from soil eDNA
   metabarcoding for each sample. Values represent the
   proportion of reads assigned to each genus relative to the
   total nematode reads per sample.

   ------------------------------------------------------------
   Sheet 5: Nematode-extracted DNA
   ------------------------------------------------------------
   Genus‑level relative read abundance from nematode‑extracted
   DNA metabarcoding for each sample. Values represent the
   proportion of reads assigned to each genus relative to the
   total nematode reads per sample.

-------------------------------------------------------------------------------
Notes
-------------------------------------------------------------------------------

- Taxonomic classification of nematode genera was cross‑validated against the
  Nemaplex online database (http://nemaplex.ucdavis.edu) for both morphological
  and molecular datasets.
- All processed data are derived from the raw files provided in the
  01_Raw_Data folder.
- These processed data files serve as the direct input for the R analysis
  scripts provided in the 03_Code folder.


================================================================================
Folder: 3 Code
================================================================================

This folder contains all R scripts used for data processing, statistical
analyses, and figure generation in this study. Scripts are numbered in the
order of execution.

-------------------------------------------------------------------------------
File list and descriptions
-------------------------------------------------------------------------------

1. 1_Rarefaction.R

   Rarefaction of ASV tables to normalize sequencing depth across samples.
   Input: Raw ASV tables (from 01_Raw_Data or processed tables from
          02_Processed_Data)
   Output: Rarefied ASV tables for diversity analyses
   Key functions: rrarefy() in the vegan package

2. 2_Relative composition.R

   Calculation of relative read abundance for each taxon per sample.
   Input: Filtered nematode ASV tables
   Output: Relative abundance matrices at genus and family levels
   Key functions: tidyverse package (dplyr, tidyr)

3. 3_adjust_p.R

   Benjamini-Hochberg false discovery rate (FDR) correction for all p-values
   derived from multiple linear regressions (Figures 2, 4-6).
   Input: Raw p-values from SPSS linear regressions
   Output: Adjusted p-values
   Key functions: p.adjust() in base R with method = "BH"

4. 4_LMMs.R

   Beta mixed-effects models to separate the independent effects of abundance
   and mean individual mass on DNA relative read abundance.
   Input: Processed community data (from 02_Processed_Data)
   Output: Model summaries, fixed and random effects estimates, diagnostic plots
   Key functions: glmmTMB() in the glmmTMB package; simulateResiduals() in the
                 DHARMa package

5. 5_Bland-Altman analysis.R

   Robust Bland-Altman analysis comparing DNA-based relative read abundance
   against two reference metrics: morphology-based relative abundance and
   mass proportion.
   Input: Processed community data (from 02_Processed_Data)
   Output: LOA summaries (median differences, 2.5% and 97.5% percentiles),
           bootstrap tests for LOA width differences, Mann-Whitney U tests
           for absolute residuals, and exported data for OriginPro figures
   Key functions: Base R functions for quantiles, medians, and bootstrapping

-------------------------------------------------------------------------------
Software requirements
-------------------------------------------------------------------------------

All scripts were run using R version 4.5.1 with the following packages:

- vegan (v2.6-8)       : rarefaction
- tidyverse (v2.0.0)   : data manipulation and relative abundance calculation
- glmmTMB (v1.1.10)    : beta mixed-effects models
- DHARMa (v0.4.7)      : model diagnostics

-------------------------------------------------------------------------------
Notes
-------------------------------------------------------------------------------

- All scripts assume that the input data files are located in the
  ../02_Processed_Data/ directory.
- Output files are saved in the same directory as the scripts unless otherwise
  specified in the code.
- Figures were finalized in OriginPro 2024; R scripts export intermediate data
  for figure generation.
- Linear regressions (R², p-values, and regression equations) were performed
  in SPSS Statistics 22 and are therefore not included in these R scripts.


================================================================================
Corresponding author
================================================================================

Yongxin Lin
Fujian Provincial Key Laboratory for Subtropical Resources and Environment
Fujian Normal University, Fuzhou 350117, China
E-mail: yxlin@fjnu.edu.cn