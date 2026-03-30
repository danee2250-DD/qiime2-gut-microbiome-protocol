# R Analysis Protocol
## Rat Gut Microbiome Study — 16S rRNA Analysis

**Author:** Daneesha  
**R Version:** 4.5.2  
**Last Updated:** March 2026

---

## Overview

This protocol documents the complete downstream R analysis pipeline for the rat gut microbiome study. Analysis covers alpha diversity, beta diversity, taxonomic composition, differential abundance, longitudinal trends, and functional prediction (PICRUSt2).

---

## Study Design

| Parameter | Details |
|---|---|
| Organism | Rat gut microbiome |
| Model | Spontaneously Hypertensive Rat (SHR) |
| Samples | 108 total |
| Subjects | 6 rats |
| Timepoints | Day 0, Day 28, Day 56 |
| Treatment groups | MSC, EV, Combo, Positive control, Negative control, Normal |
| Reference group | Positive control (gold standard drug treatment) |
| ASVs | 5,817 (post Deblur denoising) |
| Rarefaction depth | 4,000 reads per sample |

---

## R Packages Required
```r
install.packages(c(
  "tidyverse", "vegan", "pheatmap", "RColorBrewer",
  "ggplot2", "scales", "rstatix", "nlme", "ggpubr",
  "patchwork", "dunn.test", "ggrepel"
))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("phyloseq")

install.packages("ggpicrust2")
```

---

## File Structure
```
QIIME_Analysis/
├── R_exports/
│   ├── Metadata.tsv
│   ├── feature-table/
│   │   └── feature-table-clean.tsv
│   ├── taxonomy/
│   │   └── taxonomy.tsv
│   ├── tree/
│   │   └── tree.nwk
│   ├── alpha/
│   │   ├── faith_pd/alpha-diversity.tsv
│   │   ├── shannon/alpha-diversity.tsv
│   │   ├── evenness/alpha-diversity.tsv
│   │   └── observed_features/alpha-diversity.tsv
│   ├── beta/
│   │   ├── unweighted_unifrac/distance-matrix.tsv
│   │   ├── weighted_unifrac/distance-matrix.tsv
│   │   └── bray_curtis/distance-matrix.tsv
│   ├── ancombc/
│   │   ├── Day_0/
│   │   ├── Day_28/
│   │   ├── Day_56/
│   │   ├── l6_Day_0/
│   │   ├── l6_Day_28/
│   │   └── l6_Day_56/
│   └── picrust2/
│       ├── pathway_abundance.tsv
│       ├── EC_metagenome.tsv
│       ├── KO_metagenome.tsv
│       └── nsti_scores.tsv
├── R_scripts/
│   ├── 00_setup.R
│   ├── 01_alpha_diversity.R
│   ├── 02_beta_diversity.R
│   ├── 03_taxonomy.R
│   ├── 04_differential_abundance.R
│   ├── 05_longitudinal.R
│   └── 06_picrust2.R
└── Plots/
    ├── alpha_diversity/
    ├── beta_diversity/
    ├── taxonomy/
    ├── differential_abundance/
    ├── longitudinal/
    └── picrust2/
```

---

## Running Order

Always run scripts in order. Each script depends on objects created by `00_setup.R`.
```r
# Step 1 — Run once per R session (required)
source("/home/daneesha/QIIME_Analysis/R_scripts/00_setup.R")

# Step 2 — Run each analysis script
source("/home/daneesha/QIIME_Analysis/R_scripts/01_alpha_diversity.R")
source("/home/daneesha/QIIME_Analysis/R_scripts/02_beta_diversity.R")
source("/home/daneesha/QIIME_Analysis/R_scripts/03_taxonomy.R")
source("/home/daneesha/QIIME_Analysis/R_scripts/04_differential_abundance.R")
source("/home/daneesha/QIIME_Analysis/R_scripts/05_longitudinal.R")
source("/home/daneesha/QIIME_Analysis/R_scripts/06_picrust2.R")
```

---

## Script Descriptions

### 00_setup.R — Setup and Data Loading
- Loads all required R packages
- Sets base paths and creates plot directories
- Loads metadata, feature table, taxonomy, phylogenetic tree
- Builds phyloseq object (5,817 ASVs × 108 samples)
- Defines group colour and shape palettes

### 01_alpha_diversity.R — Within-Sample Diversity
**Metrics:** Faith PD, Shannon, Pielou Evenness, Observed ASVs, Chao1, Simpson, InvSimpson

**Statistical methods:**
- Linear Mixed Effects Model (LME) — accounts for repeated measures per rat
- Kruskal-Wallis test — between-group comparison per timepoint
- Dunn post-hoc test with Benjamini-Hochberg (BH) FDR correction

**Outputs:** Boxplots per metric (with and without stats), combined thesis figure

### 02_beta_diversity.R — Between-Sample Diversity
**Distance metrics:** Unweighted UniFrac, Weighted UniFrac, Bray-Curtis

**Statistical methods:**
- PERMANOVA (adonis2, 999 permutations) — group and timepoint effects
- ANOSIM — community structure differences
- Betadisper — homogeneity of dispersion

**Outputs:** PCoA plots by group and timepoint, statistics CSV

### 03_taxonomy.R — Taxonomic Composition
**Levels:** Phylum, Family, Genus

**Outputs:**
- Stacked bar plots by group
- Stacked bar plots faceted by group + timepoint
- Genus heatmap (top 15)
- Firmicutes:Bacteroidota ratio plot

### 04_differential_abundance.R — Differential Taxa
**Method:** ANCOM-BC (run in QIIME2)
**Reference group:** Positive control per timepoint
**FDR correction:** Benjamini-Hochberg

**Outputs:**
- Volcano plots per timepoint (genus level)
- DA heatmap Day 56
- Bubble plot all timepoints

### 05_longitudinal.R — Temporal Trends
**Outputs:**
- Line plots for all 7 alpha metrics (individual rat lines + mean ± SE)
- Combined main thesis figure (Shannon + Faith PD + Observed)
- Firmicutes:Bacteroidota ratio longitudinal
- Top 6 genera longitudinal trends

### 06_picrust2.R — Functional Prediction
**Input:** PICRUSt2 outputs (pathway_abundance.tsv, EC_metagenome.tsv, KO_metagenome.tsv)

**Outputs:**
- NSTI quality histogram
- MetaCyc pathway bar plots (by group, by group + timepoint)
- MetaCyc pathway heatmap
- Functional PCoA (Bray-Curtis, PERMANOVA)
- KEGG Level 2 category bar plot
- KEGG Level 2 hierarchical clustering heatmap
- LEfSe-style significant pathway bar plot
- Spearman correlation heatmap (taxa vs pathways)

---

## Statistical Methods Summary

| Analysis | Method | Correction |
|---|---|---|
| Alpha diversity (within-group) | Wilcoxon rank-sum | BH FDR |
| Alpha diversity (between-group) | Kruskal-Wallis + Dunn | BH FDR |
| Alpha diversity (longitudinal) | Linear Mixed Effects | — |
| Beta diversity | PERMANOVA + ANOSIM | — |
| Differential abundance | ANCOM-BC | BH FDR |
| Functional prediction DAA | LinDA | BH FDR |
| Functional community | PERMANOVA Bray-Curtis | — |
| Taxa-function correlation | Spearman | BH FDR |

---

## Colour Palette

| Group | Colour | Shape |
|---|---|---|
| MSC | #2196F3 (Blue) | Circle |
| EV | #4CAF50 (Green) | Triangle up |
| Combo | #FF9800 (Orange) | Square |
| Positive | #E91E63 (Pink) | Diamond |
| Negative | #9C27B0 (Purple) | Star |
| Normal | #4E342E (Brown) | Triangle down |

---

## Key Results Summary

| Analysis | Key Finding |
|---|---|
| Alpha diversity | All groups equal at Day 0; MSC/EV significant richness increase by Day 56 |
| Beta diversity | PERMANOVA significant (p<0.01); Unweighted UniFrac R²=0.125 |
| Taxonomy | Firmicutes + Bacteroidota dominant (~90%); F:B ratio converges to ~1.0 by Day 56 |
| Differential abundance | Combo most similar to Positive control at Day 56 (1 difference) |
| Functional prediction | PERMANOVA R²=0.096, p=0.004; EV+Combo cluster with Positive control |

---

## References

- McMurdie PJ & Holmes S (2013). phyloseq: An R Package for Reproducible Interactive Analysis. *PLoS ONE*, 8(4).
- Oksanen J et al. (2022). vegan: Community Ecology Package. R package version 2.6-4.
- Lin H & Peddada SD (2020). Analysis of compositions of microbiomes with bias correction. *Nature Communications*, 11, 3514.
- Douglas GM et al. (2020). PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology*, 38, 685–688.
- Yang C et al. (2023). ggpicrust2: an R package for PICRUSt2 analysis. *Bioinformatics*, 39(8).
