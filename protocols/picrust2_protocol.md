# PICRUSt2 Functional Prediction Protocol
## Rat Gut Microbiome Study — 16S rRNA Analysis

**Author:** Daneesha  
**PICRUSt2 Version:** 2.5.2  
**Placement Tool:** epa-ng  
**Last Updated:** March 2026

---

## Overview

PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States) predicts the functional potential of microbial communities from 16S rRNA amplicon sequencing data. This protocol documents the complete functional prediction pipeline applied to rat gut microbiome data.

---

## System Requirements

- Ubuntu 24.04 (WSL2)
- RAM: minimum 16 GB allocated to WSL2
- Conda: Miniconda or Anaconda
- QIIME2 pipeline must be completed first (see README.md)

### WSL2 Memory Configuration
Edit `C:\Users\<username>\.wslconfig`:
```
[wsl2]
memory=16GB
processors=8
swap=8GB
```

---

## Step 1 — Install PICRUSt2
```bash
conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.5.2
conda activate picrust2
picrust2_pipeline.py --version
```

---

## Step 2 — Export Required Files from QIIME2
```bash
conda activate qiime2-amplicon-2025.7

# Export representative sequences
qiime tools export \
  --input-path /home/daneesha/QIIME_Analysis/results/rep-seqs-deblur.qza \
  --output-path /home/daneesha/QIIME_Analysis/results/picrust2_input

# Export feature table
qiime tools export \
  --input-path /home/daneesha/QIIME_Analysis/results/table-deblur.qza \
  --output-path /home/daneesha/QIIME_Analysis/results/picrust2_input

# Convert feature table to TSV
biom convert \
  -i /home/daneesha/QIIME_Analysis/results/picrust2_input/feature-table.biom \
  -o /home/daneesha/QIIME_Analysis/results/picrust2_input/feature-table.tsv \
  --to-tsv

echo "Export complete!"
ls /home/daneesha/QIIME_Analysis/results/picrust2_input/
```

---

## Step 3 — Run PICRUSt2 Pipeline
```bash
conda activate picrust2

picrust2_pipeline.py \
  -s /home/daneesha/QIIME_Analysis/results/picrust2_input/dna-sequences.fasta \
  -i /home/daneesha/QIIME_Analysis/results/picrust2_input/feature-table.biom \
  -o /home/daneesha/QIIME_Analysis/results/picrust2_results \
  -p 8 \
  --placement_tool epa-ng \
  --verbose

echo "PICRUSt2 pipeline complete!"
```

**Parameters:**
| Parameter | Value | Description |
|---|---|---|
| `-p` | 8 | Number of CPU threads |
| `--placement_tool` | epa-ng | Sequence placement algorithm |
| `--verbose` | - | Show progress messages |

**Note:** epa-ng is used instead of the default pplacer due to compatibility issues with Ubuntu 24.04.

---

## Step 4 — Quality Assessment (NSTI)
```bash
# Decompress NSTI scores
gunzip -c /home/daneesha/QIIME_Analysis/results/picrust2_results/marker_predicted_and_nsti.tsv.gz \
  > /home/daneesha/QIIME_Analysis/R_exports/picrust2/nsti_scores.tsv

# Check NSTI summary
awk 'NR>1 {sum+=$3; count++} END {print "Mean NSTI:", sum/count}' \
  /home/daneesha/QIIME_Analysis/R_exports/picrust2/nsti_scores.tsv
```

**NSTI Thresholds:**
| Range | Quality |
|---|---|
| < 0.06 | Excellent |
| 0.06 – 0.15 | Good |
| 0.15 – 2.0 | Acceptable |
| > 2.0 | Poor — exclude |

**Study Results:** Mean NSTI = 0.31 | Median = 0.225 | All 5,813 ASVs below threshold ✅

---

## Step 5 — Export PICRUSt2 Outputs
```bash
mkdir -p /home/daneesha/QIIME_Analysis/R_exports/picrust2

# MetaCyc pathway abundance
gunzip -c /home/daneesha/QIIME_Analysis/results/picrust2_results/pathways_out/path_abun_unstrat.tsv.gz \
  > /home/daneesha/QIIME_Analysis/R_exports/picrust2/pathway_abundance.tsv

# EC number abundance
gunzip -c /home/daneesha/QIIME_Analysis/results/picrust2_results/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
  > /home/daneesha/QIIME_Analysis/R_exports/picrust2/EC_metagenome.tsv

# KO abundance
gunzip -c /home/daneesha/QIIME_Analysis/results/picrust2_results/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
  > /home/daneesha/QIIME_Analysis/R_exports/picrust2/KO_metagenome.tsv

# NSTI scores
gunzip -c /home/daneesha/QIIME_Analysis/results/picrust2_results/marker_predicted_and_nsti.tsv.gz \
  > /home/daneesha/QIIME_Analysis/R_exports/picrust2/nsti_scores.tsv

echo "All PICRUSt2 files exported!"
ls /home/daneesha/QIIME_Analysis/R_exports/picrust2/
```

---

## Step 6 — R Analysis

See `R_scripts/06_picrust2.R` for complete downstream analysis including:
- NSTI quality assessment
- MetaCyc pathway bar plots
- MetaCyc pathway heatmap
- Functional PCoA (Bray-Curtis)
- KEGG Level 2 category bar plot
- KEGG Level 2 hierarchical clustering heatmap
- LEfSe-style significant pathway bar plot
- Spearman correlation — taxa vs pathways

---

## Output Files

| File | Description |
|---|---|
| `pathway_abundance.tsv` | MetaCyc pathway abundances per sample |
| `EC_metagenome.tsv` | Enzyme Commission abundances per sample |
| `KO_metagenome.tsv` | KEGG Orthology abundances per sample |
| `nsti_scores.tsv` | NSTI quality scores per ASV |

---

## Key Results

| Metric | Value |
|---|---|
| Total ASVs | 5,813 |
| Mean NSTI | 0.31 |
| Median NSTI | 0.225 |
| ASVs below threshold | 100% |
| MetaCyc pathways predicted | 420+ |
| PERMANOVA R² (treatment) | 0.096 |
| PERMANOVA p-value | 0.004 |

---

## References

- Douglas GM et al. (2020). PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology*, 38, 685–688.
- Yang C et al. (2023). ggpicrust2: an R package for PICRUSt2 predicted functional profile analysis and visualization. *Bioinformatics*, 39(8).
- Wright RJ & Langille MGI (2025). PICRUSt2-SC: an update to the reference database. *Bioinformatics*, 41(5).
