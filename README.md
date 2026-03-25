# QIIME 2 Gut Microbiome Analysis Protocol
## Paired-End 16S rRNA Amplicon Sequencing — Deblur Pipeline

**Author:** Daneesha  
**QIIME 2 Version:** qiime2-amplicon-2025.7  
**Last Updated:** March 2026  
**Reference Tutorial:** [QIIME 2 Moving Pictures Tutorial](https://docs.qiime2.org/2024.10/tutorials/moving-pictures/)

---

## Study Design

| Parameter | Details |
|---|---|
| Organism | Rat gut microbiome |
| Samples | 108 total |
| Subjects | 6 rats |
| Timepoints | Day 0, Day 28, Day 56 |
| Treatment groups | MSC, EV, Combo, Positive control, Negative control, Normal |
| Sequencing | Paired-end Illumina, 251 bp reads |
| Target region | 16S rRNA V4 (515F/806R) |
| Denoising method | Deblur |
| Taxonomic database | SILVA 138 |

---

## Required Files

| File | Description |
|---|---|
| `manifest_paired.tsv` | Tab-separated manifest with forward and reverse read paths per sample |
| `Metadata.tsv` | Tab-separated sample metadata with `#q2:types` row |
| Raw `.fq` files | Paired-end reads named `{sample-id}-R1.fq` and `{sample-id}-R2.fq` |

### Manifest format (manifest_paired.tsv)
```
sample-id	forward-absolute-filepath	reverse-absolute-filepath
M1	/home/daneesha/QIIME_Analysis/Raw_Data/M1-R1.fq	/home/daneesha/QIIME_Analysis/Raw_Data/M1-R2.fq
M2	/home/daneesha/QIIME_Analysis/Raw_Data/M2-R1.fq	/home/daneesha/QIIME_Analysis/Raw_Data/M2-R2.fq
...
```

### Metadata format (Metadata.tsv)
```
sample-id	subject	treatment	timepoints
#q2:types	categorical	categorical	categorical
M1	Rat 1	MSC-0	Day 0
M2	Rat 1	MSC-0	Day 0
...
```

### Generate manifest and convert metadata using Python
```bash
# Generate manifest_paired.tsv directly on server
python3 -c "
import pandas as pd
samples = ['M' + str(i) for i in range(1, 109)]
rows = []
for s in samples:
    rows.append({
        'sample-id': s,
        'forward-absolute-filepath': f'/home/daneesha/QIIME_Analysis/Raw_Data/{s}-R1.fq',
        'reverse-absolute-filepath': f'/home/daneesha/QIIME_Analysis/Raw_Data/{s}-R2.fq'
    })
df = pd.DataFrame(rows)
df.to_csv('/home/daneesha/QIIME_Analysis/Raw_Data/manifest_paired.tsv', sep='\t', index=False)
print('manifest_paired.tsv saved')
"

# Install openpyxl if needed
pip install openpyxl

# Convert Metadata.xlsx to TSV
python3 -c "
import pandas as pd
df = pd.read_excel('/home/daneesha/QIIME_Analysis/Raw_Data/Metadata.xlsx')
df.columns = ['sample-id', 'subject', 'treatment', 'timepoints']
df.iloc[0] = ['#q2:types', 'categorical', 'categorical', 'categorical']
df.to_csv('/home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv', sep='\t', index=False)
print('Metadata.tsv saved')
"
```

---

## Setup

```bash
# Activate QIIME 2 environment
conda activate qiime2-amplicon-2025.7

# Create results directory
mkdir -p /home/daneesha/QIIME_Analysis/results
cd /home/daneesha/QIIME_Analysis
```

---

## Step 1 — Import Raw Reads

Imports pre-demultiplexed paired-end reads into a QIIME 2 artifact using the manifest file.

```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /home/daneesha/QIIME_Analysis/Raw_Data/manifest_paired.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path /home/daneesha/QIIME_Analysis/results/demux-paired-end.qza
```

**Output:** `demux-paired-end.qza`

---

## Step 2 — Visualize Read Quality

Generates an interactive quality plot to assess per-base quality scores before filtering.

```bash
qiime demux summarize \
  --i-data /home/daneesha/QIIME_Analysis/results/demux-paired-end.qza \
  --o-visualization /home/daneesha/QIIME_Analysis/results/demux-summary.qzv
```

**Output:** `demux-summary.qzv`

> **Review at [view.qiime2.org](https://view.qiime2.org):** Open the Interactive Quality Plot tab. Note the bp position where median quality drops below Q25 — this informs your trim length for Deblur. In this study, reads were 251 bp with high quality throughout.

---

## Step 3a — Quality Filter

Pre-filters reads based on Phred quality scores. Required before Deblur denoising.

```bash
qiime quality-filter q-score \
  --i-demux /home/daneesha/QIIME_Analysis/results/demux-paired-end.qza \
  --o-filtered-sequences /home/daneesha/QIIME_Analysis/results/demux-filtered.qza \
  --o-filter-stats /home/daneesha/QIIME_Analysis/results/demux-filter-stats.qza

qiime metadata tabulate \
  --m-input-file /home/daneesha/QIIME_Analysis/results/demux-filter-stats.qza \
  --o-visualization /home/daneesha/QIIME_Analysis/results/demux-filter-stats.qzv
```

**Outputs:** `demux-filtered.qza`, `demux-filter-stats.qza`, `demux-filter-stats.qzv`

> **Review `demux-filter-stats.qzv`:** Confirm the majority of reads per sample passed filtering. Flag any sample that lost more than 50% of reads.

---

## Step 3b — Join Paired-End Reads

Merges forward and reverse reads into single joined sequences. Required before Deblur as Deblur operates on single-end or joined reads only.

```bash
qiime vsearch merge-pairs \
  --i-demultiplexed-seqs /home/daneesha/QIIME_Analysis/results/demux-filtered.qza \
  --o-merged-sequences /home/daneesha/QIIME_Analysis/results/demux-merged.qza \
  --o-unmerged-sequences /home/daneesha/QIIME_Analysis/results/demux-unmerged.qza

qiime demux summarize \
  --i-data /home/daneesha/QIIME_Analysis/results/demux-merged.qza \
  --o-visualization /home/daneesha/QIIME_Analysis/results/demux-merged-summary.qzv
```

**Outputs:** `demux-merged.qza`, `demux-unmerged.qza`, `demux-merged-summary.qzv`

> **Review `demux-merged-summary.qzv`:** Check the Interactive Quality Plot tab on merged reads. Find the bp position where quality drops below Q25 — this is your `--p-trim-length` for Deblur. In this study, merged reads were 250 bp with consistently high quality.

---

## Step 4 — Deblur Denoising

Denoises merged reads using pre-trained 16S error profiles to generate high-quality Amplicon Sequence Variants (ASVs). Trim length was set to 250 bp based on merged read quality plots.

```bash
qiime deblur denoise-16S \
  --i-demultiplexed-seqs /home/daneesha/QIIME_Analysis/results/demux-merged.qza \
  --p-trim-length 250 \
  --o-representative-sequences /home/daneesha/QIIME_Analysis/results/rep-seqs.qza \
  --o-table /home/daneesha/QIIME_Analysis/results/table.qza \
  --p-sample-stats \
  --o-stats /home/daneesha/QIIME_Analysis/results/deblur-stats.qza \
  --p-jobs-to-start 1

qiime deblur visualize-stats \
  --i-deblur-stats /home/daneesha/QIIME_Analysis/results/deblur-stats.qza \
  --o-visualization /home/daneesha/QIIME_Analysis/results/deblur-stats.qzv
```

**Outputs:** `rep-seqs.qza`, `table.qza`, `deblur-stats.qza`, `deblur-stats.qzv`

> **Review `deblur-stats.qzv`:** Confirm most samples retain a good number of reads after denoising (ideally >1,000 reads per sample). If many samples drop to near zero, reduce `--p-trim-length` by 20–30 bp and re-run.

**Results in this study:**
- Total samples: 108 (all retained)
- Unique ASVs: 5,817
- Minimum reads per sample: 4,342
- Median reads per sample: 11,157
- Maximum reads per sample: 25,861

---

## Step 5 — Summarize Feature Table

Generates summaries of the feature table and representative sequences for quality assessment and rarefaction depth selection.

```bash
qiime feature-table summarize \
  --i-table /home/daneesha/QIIME_Analysis/results/table.qza \
  --o-visualization /home/daneesha/QIIME_Analysis/results/table.qzv \
  --m-sample-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data /home/daneesha/QIIME_Analysis/results/rep-seqs.qza \
  --o-visualization /home/daneesha/QIIME_Analysis/results/rep-seqs.qzv
```

**Outputs:** `table.qzv`, `rep-seqs.qzv`

> **Review `table.qzv` — Interactive Sample Detail tab:**
> - **Rarefaction depth:** Highest depth that retains all samples → used in Step 7 (`--p-sampling-depth`)
> - **Max depth:** Maximum frequency per sample → used in Step 8 (`--p-max-depth`)
>
> In this study: rarefaction depth = **4,000**, max depth = **25,000**

---

## Step 6 — Build Phylogenetic Tree

Constructs a rooted phylogenetic tree from representative sequences using MAFFT alignment and FastTree. Required for phylogenetic diversity metrics (Faith's PD, UniFrac).

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /home/daneesha/QIIME_Analysis/results/rep-seqs.qza \
  --o-alignment /home/daneesha/QIIME_Analysis/results/aligned-rep-seqs.qza \
  --o-masked-alignment /home/daneesha/QIIME_Analysis/results/masked-aligned-rep-seqs.qza \
  --o-tree /home/daneesha/QIIME_Analysis/results/unrooted-tree.qza \
  --o-rooted-tree /home/daneesha/QIIME_Analysis/results/rooted-tree.qza
```

**Outputs:** `aligned-rep-seqs.qza`, `masked-aligned-rep-seqs.qza`, `unrooted-tree.qza`, `rooted-tree.qza`

Pipeline steps:
1. MAFFT multiple sequence alignment
2. Masking of hypervariable positions
3. FastTree maximum-likelihood phylogenetic tree
4. Midpoint rooting

---

## Step 7 — Alpha and Beta Diversity

Computes core diversity metrics including alpha diversity (within-sample richness) and beta diversity (between-sample composition). Rarefies all samples to a uniform depth of 4,000 reads.

```bash
# Core diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny /home/daneesha/QIIME_Analysis/results/rooted-tree.qza \
  --i-table /home/daneesha/QIIME_Analysis/results/table.qza \
  --p-sampling-depth 4000 \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --output-dir /home/daneesha/QIIME_Analysis/results/core-metrics-results

# Alpha diversity significance — by treatment
qiime diversity alpha-group-significance \
  --i-alpha-diversity /home/daneesha/QIIME_Analysis/results/core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --o-visualization /home/daneesha/QIIME_Analysis/results/core-metrics-results/faith-pd-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity /home/daneesha/QIIME_Analysis/results/core-metrics-results/shannon_vector.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --o-visualization /home/daneesha/QIIME_Analysis/results/core-metrics-results/shannon-significance.qzv

# Beta diversity significance — by treatment
qiime diversity beta-group-significance \
  --i-distance-matrix /home/daneesha/QIIME_Analysis/results/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization /home/daneesha/QIIME_Analysis/results/core-metrics-results/unifrac-treatment-significance.qzv \
  --p-pairwise

# Beta diversity significance — by timepoint
qiime diversity beta-group-significance \
  --i-distance-matrix /home/daneesha/QIIME_Analysis/results/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --m-metadata-column timepoints \
  --o-visualization /home/daneesha/QIIME_Analysis/results/core-metrics-results/unifrac-timepoints-significance.qzv \
  --p-pairwise
```

**Outputs directory:** `core-metrics-results/`

Alpha diversity metrics computed: Shannon diversity, Faith's Phylogenetic Diversity, Observed Features, Pielou's Evenness

Beta diversity metrics computed: Jaccard, Bray-Curtis, Unweighted UniFrac, Weighted UniFrac

---

## Step 8 — Alpha Rarefaction

Plots species richness as a function of sequencing depth to assess whether sequencing depth was sufficient to capture community diversity.

```bash
qiime diversity alpha-rarefaction \
  --i-table /home/daneesha/QIIME_Analysis/results/table.qza \
  --i-phylogeny /home/daneesha/QIIME_Analysis/results/rooted-tree.qza \
  --p-max-depth 25000 \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --o-visualization /home/daneesha/QIIME_Analysis/results/alpha-rarefaction.qzv
```

**Output:** `alpha-rarefaction.qzv`

> **Review `alpha-rarefaction.qzv`:** Group by `treatment` and `timepoints`. Curves that plateau indicate sufficient sequencing depth. Curves still rising suggest undersequencing.

---

## Step 9 — Taxonomic Classification

Assigns taxonomy to each ASV using a pre-trained Naive Bayes classifier against the SILVA 138 database (515F/806R region, 99% OTUs).

```bash
# Download SILVA 138 classifier
wget -c --tries=5 \
  -O /home/daneesha/QIIME_Analysis/results/silva-138-classifier.qza \
  "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"

# Classify ASVs
qiime feature-classifier classify-sklearn \
  --i-classifier /home/daneesha/QIIME_Analysis/results/silva-138-classifier.qza \
  --i-reads /home/daneesha/QIIME_Analysis/results/rep-seqs.qza \
  --o-classification /home/daneesha/QIIME_Analysis/results/taxonomy.qza \
  --p-n-jobs 1

# Visualize taxonomy table
qiime metadata tabulate \
  --m-input-file /home/daneesha/QIIME_Analysis/results/taxonomy.qza \
  --o-visualization /home/daneesha/QIIME_Analysis/results/taxonomy.qzv

# Interactive taxa bar plots
qiime taxa barplot \
  --i-table /home/daneesha/QIIME_Analysis/results/table.qza \
  --i-taxonomy /home/daneesha/QIIME_Analysis/results/taxonomy.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --o-visualization /home/daneesha/QIIME_Analysis/results/taxa-bar-plots.qzv
```

**Outputs:** `taxonomy.qza`, `taxonomy.qzv`, `taxa-bar-plots.qzv`

> **Review `taxa-bar-plots.qzv`:** Set level to Level 2 (Phylum) and sort by `treatment` then `timepoints`. For rat gut microbiome, expect dominance of Firmicutes and Bacteroidota. Look for shifts in the Firmicutes:Bacteroidota ratio across treatment groups.

---

## Step 10 — Differential Abundance Testing (ANCOM-BC)

Identifies microbial features significantly enriched or depleted relative to the Positive control (gold standard drug treatment) at each timepoint. Analysis performed at both ASV level and genus level (taxonomy level 6).

All comparisons are made against the **Positive control** as the reference group, representing the known effective treatment outcome.

```bash
# Day 0 — ASV level
qiime composition ancombc \
  --i-table /home/daneesha/QIIME_Analysis/results/Day_0-table.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --p-formula 'treatment' \
  --p-reference-levels 'treatment::Positive-0' \
  --o-differentials /home/daneesha/QIIME_Analysis/results/ancombc-Day_0.qza

# Day 28 — ASV level
qiime composition ancombc \
  --i-table /home/daneesha/QIIME_Analysis/results/Day_28-table.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --p-formula 'treatment' \
  --p-reference-levels 'treatment::Positive-28' \
  --o-differentials /home/daneesha/QIIME_Analysis/results/ancombc-Day_28.qza

# Day 56 — ASV level
qiime composition ancombc \
  --i-table /home/daneesha/QIIME_Analysis/results/Day_56-table.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --p-formula 'treatment' \
  --p-reference-levels 'treatment::Positive-56' \
  --o-differentials /home/daneesha/QIIME_Analysis/results/ancombc-Day_56.qza

# ASV-level barplots
for SAFE in Day_0 Day_28 Day_56; do
  qiime composition da-barplot \
    --i-data /home/daneesha/QIIME_Analysis/results/ancombc-${SAFE}.qza \
    --p-significance-threshold 0.05 \
    --o-visualization /home/daneesha/QIIME_Analysis/results/da-barplot-${SAFE}.qzv
  echo "Done ASV barplot: ${SAFE}"
done

# Collapse to genus level (taxonomy level 6)
for TIMEPOINT in "Day 0" "Day 28" "Day 56"; do
  SAFE=$(echo $TIMEPOINT | tr ' ' '_')
  qiime feature-table filter-samples \
    --i-table /home/daneesha/QIIME_Analysis/results/table.qza \
    --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
    --p-where "[timepoints]='${TIMEPOINT}'" \
    --o-filtered-table /home/daneesha/QIIME_Analysis/results/${SAFE}-table.qza

  qiime taxa collapse \
    --i-table /home/daneesha/QIIME_Analysis/results/${SAFE}-table.qza \
    --i-taxonomy /home/daneesha/QIIME_Analysis/results/taxonomy.qza \
    --p-level 6 \
    --o-collapsed-table /home/daneesha/QIIME_Analysis/results/${SAFE}-table-l6.qza
done

# Day 0 — genus level
qiime composition ancombc \
  --i-table /home/daneesha/QIIME_Analysis/results/Day_0-table-l6.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --p-formula 'treatment' \
  --p-reference-levels 'treatment::Positive-0' \
  --o-differentials /home/daneesha/QIIME_Analysis/results/l6-ancombc-Day_0.qza

# Day 28 — genus level
qiime composition ancombc \
  --i-table /home/daneesha/QIIME_Analysis/results/Day_28-table-l6.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --p-formula 'treatment' \
  --p-reference-levels 'treatment::Positive-28' \
  --o-differentials /home/daneesha/QIIME_Analysis/results/l6-ancombc-Day_28.qza

# Day 56 — genus level
qiime composition ancombc \
  --i-table /home/daneesha/QIIME_Analysis/results/Day_56-table-l6.qza \
  --m-metadata-file /home/daneesha/QIIME_Analysis/Raw_Data/Metadata.tsv \
  --p-formula 'treatment' \
  --p-reference-levels 'treatment::Positive-56' \
  --o-differentials /home/daneesha/QIIME_Analysis/results/l6-ancombc-Day_56.qza

# Genus-level barplots
for SAFE in Day_0 Day_28 Day_56; do
  qiime composition da-barplot \
    --i-data /home/daneesha/QIIME_Analysis/results/l6-ancombc-${SAFE}.qza \
    --p-significance-threshold 0.05 \
    --p-level-delimiter ';' \
    --o-visualization /home/daneesha/QIIME_Analysis/results/l6-da-barplot-${SAFE}.qzv
  echo "Done genus barplot: ${SAFE}"
done
```

**Outputs per timepoint:**
- `ancombc-{timepoint}.qza` — ASV-level differential abundance
- `da-barplot-{timepoint}.qzv` — ASV-level visualization
- `l6-ancombc-{timepoint}.qza` — Genus-level differential abundance
- `l6-da-barplot-{timepoint}.qzv` — Genus-level visualization

> **Interpreting results:** Positive log-fold change (LFC) = enriched relative to Positive control. Negative LFC = depleted. Only features with q-value < 0.05 are shown. If no features appear, the group was not significantly different from the Positive control at that timepoint — this is a biological result, not an error.

---

## Output Files Summary

| File | Step | Description |
|---|---|---|
| `demux-paired-end.qza` | 1 | Imported raw reads |
| `demux-summary.qzv` | 2 | Read quality visualization |
| `demux-filtered.qza` | 3a | Quality-filtered reads |
| `demux-merged.qza` | 3b | Merged paired-end reads |
| `rep-seqs.qza` | 4 | Deblur representative sequences (ASVs) |
| `table.qza` | 4 | Deblur feature table |
| `deblur-stats.qzv` | 4 | Deblur denoising statistics |
| `table.qzv` | 5 | Feature table summary |
| `rooted-tree.qza` | 6 | Rooted phylogenetic tree |
| `core-metrics-results/` | 7 | Alpha and beta diversity outputs |
| `alpha-rarefaction.qzv` | 8 | Rarefaction curves |
| `taxonomy.qza` | 9 | Taxonomic classifications |
| `taxa-bar-plots.qzv` | 9 | Interactive taxonomy bar plots |
| `ancombc-{timepoint}.qza` | 10 | ASV-level differential abundance |
| `l6-ancombc-{timepoint}.qza` | 10 | Genus-level differential abundance |
| `da-barplot-{timepoint}.qzv` | 10 | ASV-level ANCOM-BC barplots |
| `l6-da-barplot-{timepoint}.qzv` | 10 | Genus-level ANCOM-BC barplots |

---

## Decision Points Reference

| Step | Visualization to review | Value to extract | Used in |
|---|---|---|---|
| After Step 2 | `demux-summary.qzv` | Quality profile overview | Context only |
| After Step 3b | `demux-merged-summary.qzv` | Trim length (bp where Q < 25) | `--p-trim-length` Step 4 |
| After Step 4 | `deblur-stats.qzv` | Reads passing per sample | Diagnostic |
| After Step 5 | `table.qzv` | Rarefaction depth (min reads retaining all samples) | `--p-sampling-depth` Step 7 |
| After Step 5 | `table.qzv` | Max depth (highest sample frequency) | `--p-max-depth` Step 8 |

---

## Notes

- The `UserWarning: pkg_resources is deprecated` message appears on every command and can be safely ignored throughout the analysis
- DADA2 was considered but requires >8 GB RAM — Deblur was used as a memory-efficient alternative that produces equivalent quality results
- All comparisons in ANCOM-BC are relative to the **Positive control (gold standard drug treatment)** group
- The `--p-significance-threshold` can be relaxed to 0.1 for exploratory analysis if few features appear at 0.05

---

## Citations

- QIIME 2: Bolyen et al. (2019) Nature Biotechnology. https://doi.org/10.1038/s41587-019-0209-9
- Deblur: Amir et al. (2017) mSystems. https://doi.org/10.1128/mSystems.00191-16
- SILVA: Quast et al. (2013) Nucleic Acids Research. https://doi.org/10.1093/nar/gks1219
- ANCOM-BC: Lin & Peddada (2020) Nature Communications. https://doi.org/10.1038/s41467-020-17041-7
- FastTree: Price et al. (2010) PLoS ONE. https://doi.org/10.1371/journal.pone.0009490
- MAFFT: Katoh & Standley (2013) Molecular Biology and Evolution. https://doi.org/10.1093/molbev/mst010
