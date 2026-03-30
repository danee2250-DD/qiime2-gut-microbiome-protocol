# ============================================
# 01_alpha_diversity.R
# Alpha Diversity Analysis
# Rat Gut Microbiome Study — 16S rRNA Analysis
# Author: Daneesha
# ============================================

# Run 00_setup.R first before this script
# source("00_setup.R")

# ============================================
# LOAD ALPHA DIVERSITY DATA
# ============================================

alpha_path <- file.path(data_path, "alpha")
plot_path  <- file.path(plots_path, "alpha_diversity")

# Load all alpha diversity metrics
faith_pd <- read.table(
  file.path(alpha_path, "faith_pd/alpha-diversity.tsv"),
  header = TRUE, sep = "\t", row.names = 1
)
shannon <- read.table(
  file.path(alpha_path, "shannon/alpha-diversity.tsv"),
  header = TRUE, sep = "\t", row.names = 1
)
evenness <- read.table(
  file.path(alpha_path, "evenness/alpha-diversity.tsv"),
  header = TRUE, sep = "\t", row.names = 1
)
observed <- read.table(
  file.path(alpha_path, "observed_features/alpha-diversity.tsv"),
  header = TRUE, sep = "\t", row.names = 1
)

# Combine all metrics with metadata
alpha_df <- metadata %>%
  rownames_to_column("SampleID") %>%
  left_join(faith_pd    %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  left_join(shannon     %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  left_join(evenness    %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  left_join(observed    %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  mutate(
    group      = factor(group,
                        levels = c("MSC","EV","Combo",
                                   "Positive","Negative","Normal")),
    timepoints = factor(timepoints,
                        levels = c("Day 0","Day 28","Day 56"))
  )

# Add richness estimates from phyloseq
chao1_df <- estimate_richness(ps, measures = c("Chao1","Simpson","InvSimpson")) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = gsub("\\.", "-", SampleID))

alpha_df <- alpha_df %>%
  left_join(chao1_df, by = "SampleID")

cat("Alpha diversity data loaded:", nrow(alpha_df), "samples\n")
cat("Metrics available:", colnames(alpha_df), "\n")

# ============================================
# DEFINE METRICS FOR PLOTTING
# ============================================

metrics <- list(
  list(col = "faith_pd",          label = "Faith Phylogenetic Diversity"),
  list(col = "shannon_entropy",   label = "Shannon Diversity Index"),
  list(col = "pielou_evenness",   label = "Pielou Evenness"),
  list(col = "observed_features", label = "Observed ASVs (Richness)"),
  list(col = "Chao1",             label = "Chao1 Richness"),
  list(col = "Simpson",           label = "Simpson Index"),
  list(col = "InvSimpson",        label = "Inverse Simpson Index")
)

# ============================================
# STATISTICAL TESTS
# ============================================

cat("Running statistical tests...\n")

stats_results <- list()

for (m in metrics) {
  metric_col <- m$col
  metric_label <- m$label

  # Linear Mixed Effects Model (repeated measures)
  lme_formula <- as.formula(
    paste0(metric_col, " ~ group * timepoints")
  )
  lme_model <- tryCatch(
    lme(lme_formula,
        random = ~ 1 | subject,
        data   = alpha_df,
        na.action = na.omit),
    error = function(e) NULL
  )

  # Kruskal-Wallis by timepoint
  kw_results <- alpha_df %>%
    group_by(timepoints) %>%
    kruskal_test(as.formula(paste0(metric_col, " ~ group"))) %>%
    mutate(metric = metric_label)

  # Dunn post-hoc test
  dunn_results <- alpha_df %>%
    group_by(timepoints) %>%
    dunn_test(as.formula(paste0(metric_col, " ~ group")),
              p.adjust.method = "BH") %>%
    mutate(metric = metric_label)

  stats_results[[metric_col]] <- list(
    lme   = lme_model,
    kw    = kw_results,
    dunn  = dunn_results
  )
}

cat("Statistical tests complete\n")

# ============================================
# PLOT FUNCTION
# ============================================

plot_alpha <- function(metric_col, metric_label,
                        show_stats = FALSE) {

  plot_data <- alpha_df %>%
    filter(!is.na(.data[[metric_col]])) %>%
    mutate(group = factor(group,
                           levels = c("MSC","EV","Combo",
                                      "Positive","Negative","Normal")))

  p <- ggplot(plot_data,
              aes(x = group,
                  y = .data[[metric_col]],
                  fill = group,
                  color = group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA,
                 linewidth = 0.6) +
    geom_jitter(width = 0.15, size = 1.8,
                alpha = 0.7, shape = 16) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~ timepoints, scales = "free_y", nrow = 1) +
    labs(
      title = metric_label,
      x     = "Treatment Group",
      y     = metric_label
    ) +
    theme_bw() +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
      legend.position  = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )

  if (show_stats) {
    dunn <- stats_results[[metric_col]]$dunn
    stat_sig <- dunn %>%
      filter(p.adj < 0.05) %>%
      mutate(y.position = max(plot_data[[metric_col]],
                               na.rm = TRUE) * 1.05)

    if (nrow(stat_sig) > 0) {
      p <- p + stat_pvalue_manual(
        stat_sig,
        label          = "p.adj.signif",
        tip.length     = 0.01,
        step.increase  = 0.08,
        hide.ns        = TRUE
      )
    }
  }
  return(p)
}

# ============================================
# GENERATE AND SAVE PLOTS
# ============================================

cat("Generating alpha diversity plots...\n")

for (i in seq_along(metrics)) {
  m <- metrics[[i]]
  cat(i, "/", length(metrics), m$label, "\n")

  # Without stats
  p_no_stats <- plot_alpha(m$col, m$label, show_stats = FALSE)
  ggsave(
    file.path(plot_path, "without_stats",
              paste0(m$col, "_boxplot.png")),
    p_no_stats, width = 12, height = 6, dpi = 300
  )

  # With stats
  p_with_stats <- plot_alpha(m$col, m$label, show_stats = TRUE)
  ggsave(
    file.path(plot_path, "with_stats",
              paste0(m$col, "_boxplot.png")),
    p_with_stats, width = 12, height = 6, dpi = 300
  )
}

# ============================================
# COMBINED FIGURE — MAIN THESIS PLOT
# ============================================

p_shannon <- plot_alpha("shannon_entropy",
                         "Shannon Diversity Index", FALSE)
p_faith   <- plot_alpha("faith_pd",
                         "Faith Phylogenetic Diversity", FALSE)
p_obs     <- plot_alpha("observed_features",
                         "Observed ASVs (Richness)", FALSE)

p_combined <- p_shannon / p_faith / p_obs +
  plot_annotation(
    title = "Alpha Diversity — All Treatment Groups Across Timepoints",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  )

ggsave(
  file.path(plot_path, "alpha_diversity_combined.png"),
  p_combined, width = 14, height = 14, dpi = 300
)

cat("\n=== ALPHA DIVERSITY COMPLETE ===\n")
cat("Plots saved to:", plot_path, "\n")
cat("Without stats:", length(list.files(
  file.path(plot_path, "without_stats"))), "files\n")
cat("With stats:", length(list.files(
  file.path(plot_path, "with_stats"))), "files\n")

