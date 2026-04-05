# ============================================
# 00_setup.R
# Setup, Data Loading and Phyloseq Construction
# Rat Gut Microbiome Study — 16S rRNA Analysis
# Author: Daneesha
# QIIME2 Version: qiime2-amplicon-2025.7
# Last Updated: April 2026
# ============================================

# ============================================
# LOAD LIBRARIES
# ============================================

required_packages <- c(
  "phyloseq", "tidyverse", "vegan", "pheatmap",
  "RColorBrewer", "ggplot2", "scales", "rstatix",
  "nlme", "ggpubr", "patchwork", "dunn.test",
  "ggrepel", "ComplexHeatmap", "circlize", "grid"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("phyloseq", "ComplexHeatmap")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

cat("All libraries loaded successfully\n")

# ============================================
# SET BASE PATHS
# ============================================

base_path  <- "/home/daneesha/QIIME_Analysis"
data_path  <- file.path(base_path, "R_exports")
plots_path <- file.path(base_path, "Plots")

# Create all plot subdirectories
dirs_to_create <- c(
  "alpha_diversity/without_stats",
  "alpha_diversity/with_stats",
  "beta_diversity",
  "taxonomy",
  "differential_abundance",
  "longitudinal",
  "picrust2"
)

for (d in dirs_to_create) {
  dir.create(file.path(plots_path, d),
             recursive = TRUE, showWarnings = FALSE)
}

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

# Read feature table — skip comment lines starting with #
ft_path <- file.path(data_path, "feature-table/feature-table-clean.tsv")

# Read all lines first to handle QIIME2 header format
ft_lines <- readLines(ft_path)

# Remove lines starting with # (QIIME2 comment lines)
ft_lines <- ft_lines[!grepl("^#", ft_lines)]

# Write cleaned lines to temp file and read as table
tmp_file <- tempfile(fileext = ".tsv")
writeLines(ft_lines, tmp_file)

feature_table_raw <- read.table(
  tmp_file,
  header = TRUE, sep = "\t",
  row.names = 1, check.names = FALSE
)

file.remove(tmp_file)

otu_mat <- as.matrix(feature_table_raw)
OTU     <- otu_table(otu_mat, taxa_are_rows = TRUE)

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

# Align samples between feature table and metadata
common_samples <- intersect(colnames(OTU), rownames(metadata))
cat("Common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common samples between feature table and metadata!",
       " Check sample names.")
}

OTU_aligned  <- prune_samples(common_samples, OTU)
meta_aligned <- metadata[common_samples, , drop = FALSE]
sample_data_df <- sample_data(meta_aligned)

ps <- phyloseq(OTU_aligned, TAX, sample_data_df, tree)

# Add group variable to phyloseq sample data
sample_data(ps)$group <- factor(
  gsub("-.*", "", as.character(sample_data(ps)$treatment)),
  levels = c("MSC","EV","Combo","Positive","Negative","Normal")
)

cat("Phyloseq object created\n")
cat("  ASVs:", ntaxa(ps), "\n")
cat("  Samples:", nsamples(ps), "\n")
cat("  Sample variables:", paste(sample_variables(ps), collapse=", "), "\n")

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
cat("Groups:", paste(levels(metadata$group), collapse=", "), "\n")
cat("Timepoints:", paste(levels(metadata$timepoints), collapse=", "), "\n")
cat("\nTo run full analysis:\n")
cat("  source('03_taxonomy.R')\n")
cat("  source('04_differential_abundance.R')\n")
cat("  source('05_longitudinal.R')\n")
cat("  source('06_picrust2.R')\n")
