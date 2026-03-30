# ============================================
# 02_beta_diversity.R
# Beta Diversity Analysis
# Rat Gut Microbiome Study — 16S rRNA Analysis
# Author: Daneesha
# ============================================

# Run 00_setup.R first before this script
# source("00_setup.R")

# ============================================
# SET PATHS
# ============================================

beta_path <- file.path(data_path, "beta")
plot_path  <- file.path(plots_path, "beta_diversity")

# ============================================
# LOAD DISTANCE MATRICES
# ============================================

unweighted_unifrac <- read.table(
  file.path(beta_path, "unweighted_unifrac/distance-matrix.tsv"),
  header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE
)

weighted_unifrac <- read.table(
  file.path(beta_path, "weighted_unifrac/distance-matrix.tsv"),
  header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE
)

bray_curtis <- read.table(
  file.path(beta_path, "bray_curtis/distance-matrix.tsv"),
  header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE
)

cat("Distance matrices loaded\n")

# Align sample order with metadata
sample_order <- rownames(metadata)

unweighted_mat <- as.dist(unweighted_unifrac[sample_order, sample_order])
weighted_mat   <- as.dist(weighted_unifrac[sample_order, sample_order])
bray_mat       <- as.dist(bray_curtis[sample_order, sample_order])

meta_beta <- metadata %>%
  rownames_to_column("SampleID") %>%
  mutate(
    group = factor(group,
                   levels = c("MSC","EV","Combo",
                               "Positive","Negative","Normal")),
    timepoints = factor(timepoints,
                        levels = c("Day 0","Day 28","Day 56"))
  )

# ============================================
# PERMANOVA ANALYSIS
# ============================================

cat("Running PERMANOVA tests...\n")

set.seed(123)

perm_uw_group <- adonis2(unweighted_mat ~ group,
                          data = meta_beta, permutations = 999)
perm_uw_time  <- adonis2(unweighted_mat ~ timepoints,
                          data = meta_beta, permutations = 999)

perm_w_group  <- adonis2(weighted_mat ~ group,
                          data = meta_beta, permutations = 999)
perm_w_time   <- adonis2(weighted_mat ~ timepoints,
                          data = meta_beta, permutations = 999)

perm_bc_group <- adonis2(bray_mat ~ group,
                          data = meta_beta, permutations = 999)
perm_bc_time  <- adonis2(bray_mat ~ timepoints,
                          data = meta_beta, permutations = 999)

# ANOSIM
anosim_uw <- anosim(unweighted_mat, meta_beta$group,
                     permutations = 999)
anosim_w  <- anosim(weighted_mat,   meta_beta$group,
                     permutations = 999)
anosim_bc <- anosim(bray_mat,       meta_beta$group,
                     permutations = 999)

# Betadisper (homogeneity of dispersion)
bd_uw <- betadisper(unweighted_mat, meta_beta$group)
bd_w  <- betadisper(weighted_mat,   meta_beta$group)
bd_bc <- betadisper(bray_mat,       meta_beta$group)

# Save statistics
beta_stats <- data.frame(
  Metric = c("Unweighted UniFrac","Unweighted UniFrac",
             "Weighted UniFrac","Weighted UniFrac",
             "Bray-Curtis","Bray-Curtis"),
  Comparison = rep(c("Treatment Group","Timepoint"), 3),
  R2      = c(round(perm_uw_group$R2[1], 3),
               round(perm_uw_time$R2[1], 3),
               round(perm_w_group$R2[1], 3),
               round(perm_w_time$R2[1], 3),
               round(perm_bc_group$R2[1], 3),
               round(perm_bc_time$R2[1], 3)),
  P_value = c(perm_uw_group$`Pr(>F)`[1],
               perm_uw_time$`Pr(>F)`[1],
               perm_w_group$`Pr(>F)`[1],
               perm_w_time$`Pr(>F)`[1],
               perm_bc_group$`Pr(>F)`[1],
               perm_bc_time$`Pr(>F)`[1]),
  ANOSIM_R = c(round(anosim_uw$statistic, 3), NA,
                round(anosim_w$statistic, 3),  NA,
                round(anosim_bc$statistic, 3), NA),
  ANOSIM_p = c(anosim_uw$signif, NA,
                anosim_w$signif,  NA,
                anosim_bc$signif, NA)
)

write.csv(beta_stats,
          file.path(plot_path, "beta_diversity_statistics.csv"),
          row.names = FALSE)

cat("PERMANOVA complete\n")
print(beta_stats)

# ============================================
# PCoA PLOT FUNCTION
# ============================================

plot_pcoa <- function(dist_mat, title_label,
                       color_by = "group",
                       permanova_r2, permanova_p) {

  pcoa_result <- cmdscale(dist_mat, k = 2, eig = TRUE)

  eig  <- pcoa_result$eig
  var1 <- round(eig[1] / sum(eig[eig > 0]) * 100, 1)
  var2 <- round(eig[2] / sum(eig[eig > 0]) * 100, 1)

  pcoa_df <- data.frame(
    PC1      = pcoa_result$points[, 1],
    PC2      = pcoa_result$points[, 2],
    SampleID = rownames(pcoa_result$points)
  ) %>%
    left_join(meta_beta, by = "SampleID")

  if (color_by == "group") {
    p <- ggplot(pcoa_df,
                aes(x = PC1, y = PC2,
                    color = group, shape = group)) +
      scale_color_manual(values = group_colors,
                         name = "Treatment") +
      scale_shape_manual(values = group_shapes,
                         name = "Treatment") +
      stat_ellipse(aes(group = group), type = "t",
                   level = 0.95, linetype = "dashed",
                   linewidth = 0.7, fill = NA)
  } else {
    p <- ggplot(pcoa_df,
                aes(x = PC1, y = PC2,
                    color = timepoints, shape = timepoints)) +
      scale_color_manual(values = timepoint_colors,
                         name = "Timepoint") +
      stat_ellipse(aes(group = timepoints), type = "t",
                   level = 0.95, linetype = "dashed",
                   linewidth = 0.7, fill = NA)
  }

  p <- p +
    geom_point(size = 3.5, alpha = 0.85) +
    labs(
      title    = title_label,
      subtitle = paste0("PERMANOVA: R²=", permanova_r2,
                         ", p=", permanova_p),
      x = paste0("PC1 (", var1, "%)"),
      y = paste0("PC2 (", var2, "%)")
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )

  return(p)
}

# ============================================
# GENERATE AND SAVE PCoA PLOTS
# ============================================

cat("Generating PCoA plots...\n")

# Unweighted UniFrac — by group
p_uw_group <- plot_pcoa(
  unweighted_mat,
  "Unweighted UniFrac PCoA — by Treatment Group",
  color_by     = "group",
  permanova_r2 = round(perm_uw_group$R2[1], 3),
  permanova_p  = perm_uw_group$`Pr(>F)`[1]
)
ggsave(file.path(plot_path, "unweighted_unifrac_pcoa_group.png"),
       p_uw_group, width = 10, height = 7, dpi = 300)
cat("1/6 Unweighted UniFrac by group saved\n")

# Unweighted UniFrac — by timepoint
p_uw_time <- plot_pcoa(
  unweighted_mat,
  "Unweighted UniFrac PCoA — by Timepoint",
  color_by     = "timepoint",
  permanova_r2 = round(perm_uw_time$R2[1], 3),
  permanova_p  = perm_uw_time$`Pr(>F)`[1]
)
ggsave(file.path(plot_path, "unweighted_unifrac_pcoa_timepoint.png"),
       p_uw_time, width = 10, height = 7, dpi = 300)
cat("2/6 Unweighted UniFrac by timepoint saved\n")

# Weighted UniFrac — by group
p_w_group <- plot_pcoa(
  weighted_mat,
  "Weighted UniFrac PCoA — by Treatment Group",
  color_by     = "group",
  permanova_r2 = round(perm_w_group$R2[1], 3),
  permanova_p  = perm_w_group$`Pr(>F)`[1]
)
ggsave(file.path(plot_path, "weighted_unifrac_pcoa_group.png"),
       p_w_group, width = 10, height = 7, dpi = 300)
cat("3/6 Weighted UniFrac by group saved\n")

# Weighted UniFrac — by timepoint
p_w_time <- plot_pcoa(
  weighted_mat,
  "Weighted UniFrac PCoA — by Timepoint",
  color_by     = "timepoint",
  permanova_r2 = round(perm_w_time$R2[1], 3),
  permanova_p  = perm_w_time$`Pr(>F)`[1]
)
ggsave(file.path(plot_path, "weighted_unifrac_pcoa_timepoint.png"),
       p_w_time, width = 10, height = 7, dpi = 300)
cat("4/6 Weighted UniFrac by timepoint saved\n")

# Bray-Curtis — by group
p_bc_group <- plot_pcoa(
  bray_mat,
  "Bray-Curtis PCoA — by Treatment Group",
  color_by     = "group",
  permanova_r2 = round(perm_bc_group$R2[1], 3),
  permanova_p  = perm_bc_group$`Pr(>F)`[1]
)
ggsave(file.path(plot_path, "bray_curtis_pcoa_group.png"),
       p_bc_group, width = 10, height = 7, dpi = 300)
cat("5/6 Bray-Curtis by group saved\n")

# Bray-Curtis — by timepoint
p_bc_time <- plot_pcoa(
  bray_mat,
  "Bray-Curtis PCoA — by Timepoint",
  color_by     = "timepoint",
  permanova_r2 = round(perm_bc_time$R2[1], 3),
  permanova_p  = perm_bc_time$`Pr(>F)`[1]
)
ggsave(file.path(plot_path, "bray_curtis_pcoa_timepoint.png"),
       p_bc_time, width = 10, height = 7, dpi = 300)
cat("6/6 Bray-Curtis by timepoint saved\n")

cat("\n=== BETA DIVERSITY COMPLETE ===\n")
cat("Plots saved to:", plot_path, "\n")
cat("Files:", length(list.files(plot_path)), "\n")

