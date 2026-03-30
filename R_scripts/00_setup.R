# ============================================
# 00_setup.R
# Setup, Data Loading and Phyloseq Construction
# Rat Gut Microbiome Study — 16S rRNA Analysis
# Author: Daneesha
# QIIME2 Version: qiime2-amplicon-2025.7
# ============================================

# ============================================
# LOAD LIBRARIES
# ============================================

library(phyloseq)
library(tidyverse)
library(vegan)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(rstatix)
library(nlme)
library(ggpubr)
library(patchwork)
library(dunn.test)

cat("All libraries loaded successfully\n")

# ============================================
# SET BASE PATH
# ============================================

base_path <- "/home/daneesha/QIIME_Analysis"
data_path <- file.path(base_path, "R_exports")
plots_path <- file.path(base_path, "Plots")

# Create plot subdirectories
dir.create(file.path(plots_path, "alpha_diversity/without_stats"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_path, "alpha_diversity/with_stats"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_path, "beta_diversity"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_path, "taxonomy"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_path, "differential_abundance"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_path, "longitudinal"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(plots_path, "picrust2"),
           recursive = TRUE, showWarnings = FALSE)

cat("Directory structure created\n")

# ============================================
# LOAD METADATA
# ============================================

metadata <- read.table(
  file.path(data_path, "Metadata.tsv"),
  header = TRUE, sep = "\t", row.names = 1
)
metadata <- metadata[rownames(metadata) != "#q2:types", ]
metadata$timepoints <- factor(metadata$timepoints,
                               levels = c("Day 0", "Day 28", "Day 56"))
metadata$treatment  <- as.factor(metadata$treatment)
metadata$group <- factor(
  gsub("-.*", "", as.character(metadata$treatment)),
  levels = c("MSC", "EV", "Combo", "Positive", "Negative", "Normal")
)

cat("Metadata loaded:", nrow(metadata), "samples\n")

# ============================================
# LOAD FEATURE TABLE
# ============================================

feature_table_raw <- read.table(
  file.path(data_path, "feature-table/feature-table-clean.tsv"),
  header = TRUE, sep = "\t", row.names = 1,
  comment.char = "", check.names = FALSE
)

otu_mat <- as.matrix(feature_table_raw)
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

cat("Feature table loaded:", nrow(OTU), "ASVs x",
    ncol(OTU), "samples\n")

# ============================================
# LOAD TAXONOMY
# ============================================

taxonomy_raw <- read.table(
  file.path(data_path, "taxonomy/taxonomy.tsv"),
  header = TRUE, sep = "\t", row.names = 1,
  stringsAsFactors = FALSE
)

# Parse taxonomy string into individual ranks
tax_levels <- c("Kingdom","Phylum","Class","Order",
                 "Family","Genus","Species")

tax_parsed <- str_extract_all(
  taxonomy_raw$Taxon,
  "(?<=__)[^;]+"
)

tax_mat <- t(sapply(tax_parsed, function(x) {
  length(x) <- 7
  x
}))
colnames(tax_mat) <- tax_levels
rownames(tax_mat) <- rownames(taxonomy_raw)

tax_mat[is.na(tax_mat)] <- "Unclassified"
tax_mat <- gsub("^\\s+|\\s+$", "", tax_mat)

TAX <- tax_table(tax_mat)

cat("Taxonomy loaded:", nrow(TAX), "ASVs\n")

# ============================================
# LOAD PHYLOGENETIC TREE
# ============================================

tree <- read_tree(
  file.path(data_path, "tree/tree.nwk")
)

cat("Phylogenetic tree loaded\n")

# ============================================
# BUILD PHYLOSEQ OBJECT
# ============================================

sample_data_df <- sample_data(metadata)

ps <- phyloseq(OTU, TAX, sample_data_df, tree)

cat("Phyloseq object created\n")
cat("  ASVs:", ntaxa(ps), "\n")
cat("  Samples:", nsamples(ps), "\n")
cat("  Sample variables:", sample_variables(ps), "\n")

# ============================================
# COLOUR AND SHAPE PALETTES
# ============================================

group_colors <- c(
  MSC      = "#2196F3",
  EV       = "#4CAF50",
  Combo    = "#FF9800",
  Positive = "#E91E63",
  Negative = "#9C27B0",
  Normal   = "#4E342E"
)

group_shapes <- c(
  MSC      = 16,
  EV       = 17,
  Combo    = 15,
  Positive = 18,
  Negative = 8,
  Normal   = 25
)

timepoint_colors <- c(
  "Day 0"  = "#FFA726",
  "Day 28" = "#42A5F5",
  "Day 56" = "#66BB6A"
)

cat("\n=== SETUP COMPLETE ===\n")
cat("Phyloseq object: ps\n")
cat("Groups:", levels(metadata$group), "\n")
cat("Timepoints:", levels(metadata$timepoints), "\n")

