# ============================================
# 00_setup.R
# Gut Microbiome Analysis — Setup & Data Load
# Author: Daneesha
# Date: March 2026
# ============================================

# --- Install packages if needed ---
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "DESeq2", "microbiome"), ask = FALSE)

install.packages(c(
  "tidyverse", "vegan", "ggplot2", "ggpubr",
  "pheatmap", "RColorBrewer", "patchwork",
  "scales", "reshape2", "ape", "nlme",
  "rstatix", "data.table"
))

# --- Load all libraries ---
library(phyloseq)
library(ggplot2)
library(vegan)
library(ggpubr)
library(pheatmap)
library(tidyverse)
library(patchwork)
library(ape)
library(nlme)
library(rstatix)
library(data.table)
library(RColorBrewer)

cat("All packages loaded!\n")

# --- Set paths ---
base_dir    <- "/home/daneesha/QIIME_Analysis"
data_dir    <- file.path(base_dir, "R_exports")
plots_dir   <- file.path(base_dir, "Plots")
scripts_dir <- file.path(base_dir, "R_scripts")

# Create plot subfolders
dir.create(file.path(plots_dir, "alpha_diversity/without_stats"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_dir, "alpha_diversity/with_stats"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_dir, "beta_diversity"),                recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_dir, "taxonomy"),                      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_dir, "differential_abundance"),        recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_dir, "longitudinal"),                  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_dir, "picrust2"),                      recursive = TRUE, showWarnings = FALSE)

cat("Folder structure ready!\n")

# --- Load feature table ---
lines <- readLines(file.path(data_dir, "feature-table/feature-table-clean.tsv"))
header_line <- which(grepl("^OTU_ID", lines))
lines[header_line] <- gsub("^OTU_ID", "OTU_ID", lines[header_line])
writeLines(lines[header_line:length(lines)],
           file.path(data_dir, "feature-table/temp.tsv"))

otu_table <- fread(file.path(data_dir, "feature-table/temp.tsv"),
                   header = TRUE, sep = "\t", check.names = FALSE)
otu_table <- as.data.frame(otu_table)
rownames(otu_table) <- otu_table$OTU_ID
otu_table$OTU_ID <- NULL
if("taxonomy" %in% colnames(otu_table)){
  otu_table <- otu_table[, !colnames(otu_table) %in% "taxonomy"]
}
cat("Feature table:", nrow(otu_table), "ASVs x", ncol(otu_table), "samples\n")

# --- Load taxonomy ---
taxonomy <- read.table(
  file.path(data_dir, "taxonomy/taxonomy.tsv"),
  header = TRUE, sep = "\t", row.names = 1
)
tax_matrix <- taxonomy %>%
  select(Taxon) %>%
  separate(Taxon,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep = "; ", fill = "right") %>%
  as.matrix()
rownames(tax_matrix) <- rownames(taxonomy)
cat("Taxonomy:", nrow(tax_matrix), "ASVs\n")

# --- Load metadata ---
metadata <- read.table(
  file.path(data_dir, "Metadata.tsv"),
  header = TRUE, sep = "\t", row.names = 1
)
metadata <- metadata[rownames(metadata) != "#q2:types", ]
metadata$timepoints <- factor(metadata$timepoints,
                              levels = c("Day 0", "Day 28", "Day 56"))
metadata$treatment <- factor(metadata$treatment)
metadata$group <- factor(gsub("-.*", "", as.character(metadata$treatment)),
                         levels = c("MSC","EV","Combo","Positive","Negative","Normal"))
cat("Metadata:", nrow(metadata), "samples\n")

# --- Load phylogenetic tree ---
tree <- read_tree(file.path(data_dir, "tree/tree.nwk"))
cat("Tree:", Ntip(tree), "tips\n")

# --- Build phyloseq object ---
ps <- phyloseq(
  otu_table(as.matrix(otu_table), taxa_are_rows = TRUE),
  tax_table(tax_matrix),
  sample_data(metadata),
  phy_tree(tree)
)
print(ps)
cat("\nPhyloseq object created successfully!\n")

# --- Define consistent color palette ---
group_colors <- c(
  "MSC"      = "#2196F3",
  "EV"       = "#4CAF50",
  "Combo"    = "#FF9800",
  "Positive" = "#E91E63",
  "Negative" = "#9C27B0",
  "Normal"   = "#795548"
)

timepoint_colors <- c(
  "Day 0"  = "#FFA726",
  "Day 28" = "#42A5F5",
  "Day 56" = "#66BB6A"
)

cat("\nSetup complete! Ready for analysis.\n")
cat("=== Summary ===\n")
cat("Samples:", nsamples(ps), "\n")
cat("ASVs:", ntaxa(ps), "\n")
cat("Treatment groups:", levels(metadata$group), "\n")
cat("Timepoints:", levels(metadata$timepoints), "\n")

# Shape palettes (added for beta diversity PCoA)
group_shapes <- c(
  MSC = 16, EV = 17, Combo = 15,
  Positive = 18, Negative = 8, Normal = 25
)
timepoint_shapes <- c(
  "Day 0" = 16, "Day 28" = 17, "Day 56" = 15
)
