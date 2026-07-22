# ============================================
# 02_beta_diversity.R
# Beta Diversity Analysis (combined group + timepoint plots)
# Rat Gut Microbiome Study - 16S rRNA Analysis
# Author: Daneesha
# Last Updated: July 2026
# ============================================
# This version:
#  * Uses the variable names from 00_setup.R: data_dir, plots_dir.
#  * Reproducibility: set.seed(123) before the PERMANOVA block so the
#    999-permutation PERMANOVA / ANOSIM p-values are identical each run.
#  * ONE PCoA per distance metric, showing BOTH dimensions at once:
#      colour = treatment group (6), shape = timepoint (3).
#  * Both ellipse sets: solid coloured ellipses per group, plus grey
#    dashed ellipses per timepoint, so they are visually separable.
#  * Subtitle shows PERMANOVA for BOTH group and timepoint.
#  * Saves 3 separate figures (one per metric) AND 1 combined stacked
#    figure of all three.
#  * Requires group_colors, group_shapes, timepoint_colors,
#    timepoint_shapes from 00_setup.R.
# ============================================

# Run 00_setup.R first before this script
# source("00_setup.R")

# ============================================
# SET PATHS
# ============================================

beta_path <- file.path(data_dir, "beta")
plot_path <- file.path(plots_dir, "beta_diversity")
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

# ============================================
# LOAD DISTANCE MATRICES
# ============================================

unweighted_unifrac <- read.table(
  file.path(beta_path, "unweighted_unifrac/distance-matrix.tsv"),
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

weighted_unifrac <- read.table(
  file.path(beta_path, "weighted_unifrac/distance-matrix.tsv"),
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

bray_curtis <- read.table(
  file.path(beta_path, "bray_curtis/distance-matrix.tsv"),
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

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
# PERMANOVA ANALYSIS  (computed once, seeded for reproducibility)
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

# ANOSIM (group only)
anosim_uw <- anosim(unweighted_mat, meta_beta$group, permutations = 999)
anosim_w  <- anosim(weighted_mat,   meta_beta$group, permutations = 999)
anosim_bc <- anosim(bray_mat,       meta_beta$group, permutations = 999)

# ------------------------------------------------------------
# Save the statistics table (group + timepoint, all three metrics)
# ------------------------------------------------------------
beta_stats <- data.frame(
  key = c("uw_group","uw_time","w_group","w_time","bc_group","bc_time"),
  Metric = c("Unweighted UniFrac","Unweighted UniFrac",
             "Weighted UniFrac","Weighted UniFrac",
             "Bray-Curtis","Bray-Curtis"),
  Comparison = rep(c("Treatment Group","Timepoint"), 3),
  R2 = c(round(perm_uw_group$R2[1], 3), round(perm_uw_time$R2[1], 3),
         round(perm_w_group$R2[1], 3),  round(perm_w_time$R2[1], 3),
         round(perm_bc_group$R2[1], 3), round(perm_bc_time$R2[1], 3)),
  P_value = c(perm_uw_group$`Pr(>F)`[1], perm_uw_time$`Pr(>F)`[1],
              perm_w_group$`Pr(>F)`[1],  perm_w_time$`Pr(>F)`[1],
              perm_bc_group$`Pr(>F)`[1], perm_bc_time$`Pr(>F)`[1]),
  ANOSIM_R = c(round(anosim_uw$statistic, 3), NA,
               round(anosim_w$statistic, 3),  NA,
               round(anosim_bc$statistic, 3), NA),
  ANOSIM_p = c(anosim_uw$signif, NA,
               anosim_w$signif,  NA,
               anosim_bc$signif, NA),
  stringsAsFactors = FALSE
)

write.csv(beta_stats,
          file.path(plot_path, "beta_diversity_statistics.csv"),
          row.names = FALSE)

cat("PERMANOVA complete\n")
print(beta_stats)

# Read stats back so subtitles match the CSV exactly
beta_authoritative <- read.csv(
  file.path(plot_path, "beta_diversity_statistics.csv"),
  stringsAsFactors = FALSE)

fmt_p <- function(p) formatC(as.numeric(p), format = "g", digits = 3)
fmt_r <- function(r) formatC(as.numeric(r), format = "g", digits = 3)

# Build the two-line PERMANOVA subtitle for a metric from the CSV
perm_subtitle <- function(group_key, time_key) {
  g <- beta_authoritative[beta_authoritative$key == group_key, ]
  t <- beta_authoritative[beta_authoritative$key == time_key, ]
  paste0(
    "Group: PERMANOVA R2=", fmt_r(g$R2), ", p=", fmt_p(g$P_value),
    ifelse(as.numeric(g$P_value) < 0.05, "*", ""),
    "   |   Timepoint: PERMANOVA R2=", fmt_r(t$R2), ", p=", fmt_p(t$P_value),
    ifelse(as.numeric(t$P_value) < 0.05, "*", "")
  )
}

# ============================================
# COMBINED PCoA PLOT FUNCTION
# colour = group, shape = timepoint, both ellipse sets
# ============================================

plot_pcoa_combined <- function(dist_mat, title_label,
                                group_key, time_key) {

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

  ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
    # grey dashed ellipses per TIMEPOINT (drawn first, underneath)
    stat_ellipse(aes(group = timepoints),
                 type = "t", level = 0.95,
                 linetype = "dashed", linewidth = 0.5,
                 colour = "grey55") +
    # coloured solid ellipses per GROUP
    stat_ellipse(aes(group = group, colour = group),
                 type = "t", level = 0.95,
                 linetype = "solid", linewidth = 0.6) +
    # points: colour = group, shape = timepoint
    geom_point(aes(colour = group, shape = timepoints),
               size = 2.8, alpha = 0.9) +
    scale_color_manual(values = group_colors, name = "Treatment group") +
    scale_shape_manual(values = timepoint_shapes, name = "Timepoint") +
    labs(
      title    = title_label,
      subtitle = perm_subtitle(group_key, time_key),
      caption  = paste0("Each point = one animal. Colour = treatment group, ",
                        "shape = timepoint.\nSolid ellipses = 95% CI per group, ",
                        "grey dashed ellipses = 95% CI per timepoint."),
      x = paste0("PC1 (", var1, "%)"),
      y = paste0("PC2 (", var2, "%)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, colour = "grey30"),
      plot.caption  = element_text(hjust = 0, size = 7.5, colour = "grey50"),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

# ============================================
# GENERATE PLOTS
# ============================================

cat("Generating combined PCoA plots...\n")

p_bray <- plot_pcoa_combined(bray_mat,
  "Bray-Curtis PCoA", "bc_group", "bc_time")
ggsave(file.path(plot_path, "pcoa_bray_curtis_combined.png"),
       p_bray, width = 10, height = 7, dpi = 300, bg = "white")
cat("1/3 Bray-Curtis saved\n")

p_wuf <- plot_pcoa_combined(weighted_mat,
  "Weighted UniFrac PCoA", "w_group", "w_time")
ggsave(file.path(plot_path, "pcoa_weighted_unifrac_combined.png"),
       p_wuf, width = 10, height = 7, dpi = 300, bg = "white")
cat("2/3 Weighted UniFrac saved\n")

p_uwuf <- plot_pcoa_combined(unweighted_mat,
  "Unweighted UniFrac PCoA", "uw_group", "uw_time")
ggsave(file.path(plot_path, "pcoa_unweighted_unifrac_combined.png"),
       p_uwuf, width = 10, height = 7, dpi = 300, bg = "white")
cat("3/3 Unweighted UniFrac saved\n")

# Combined stacked figure (all three, shared legend)
combined <- (p_bray / p_wuf / p_uwuf) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Beta Diversity: All 108 Samples",
    theme = theme(plot.title = element_text(hjust = 0.5,
                                            face = "bold", size = 15))
  )

ggsave(file.path(plot_path, "beta_diversity_combined_all.png"),
       combined, width = 10, height = 20, dpi = 300, bg = "white")
cat("Combined stacked figure saved\n")

cat("\n=== BETA DIVERSITY COMPLETE ===\n")
cat("Plots saved to:", plot_path, "\n")
cat("Files:", length(list.files(plot_path)), "\n")
