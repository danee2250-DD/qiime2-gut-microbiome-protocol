# ============================================
# 99_session_info.R
# Reproducibility: record R and package versions
# Rat Gut Microbiome Study - 16S rRNA Analysis
# Author: Daneesha
# Last Updated: July 2026
# ============================================
# Run this AFTER 00_setup.R (so all analysis packages are loaded).
# It writes session_info.txt to the repo root, capturing the exact
# R version and package versions used, so results can be reproduced.
#
#   source("00_setup.R")
#   source("99_session_info.R")
# ============================================

out_file <- file.path(base_path, "session_info.txt")

sink(out_file)
cat("=================================================\n")
cat(" SESSION INFO - Rat Gut Microbiome 16S Analysis\n")
cat(" Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=================================================\n\n")

cat("--- Key analysis package versions ---\n")
key_pkgs <- c("phyloseq", "vegan", "rstatix", "nlme", "ggpubr",
              "tidyverse", "ggplot2", "patchwork", "DESeq2",
              "ComplexHeatmap", "microbiome")
for (p in key_pkgs) {
  v <- tryCatch(as.character(packageVersion(p)),
                error = function(e) "NOT INSTALLED")
  cat(sprintf("  %-16s %s\n", p, v))
}

cat("\n--- Random seed used ---\n")
cat("  set.seed(123) in 01_alpha_diversity.R and 02_beta_diversity.R\n")

cat("\n--- Full sessionInfo() ---\n\n")
print(sessionInfo())
sink()

cat("Session info written to:", out_file, "\n")
