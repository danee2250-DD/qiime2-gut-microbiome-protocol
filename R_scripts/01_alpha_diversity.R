# ============================================
# 01_alpha_diversity.R
# Alpha Diversity Analysis
# Rat Gut Microbiome Study — 16S rRNA Analysis
# Author: Daneesha
# Last Updated: April 2026
# ============================================

# Run 00_setup.R first
# source("00_setup.R")

alpha_path <- file.path(data_path, "alpha")
plot_path  <- file.path(plots_path, "alpha_diversity")

faith_pd <- read.table(
  file.path(alpha_path, "faith_pd/alpha-diversity.tsv"),
  header = TRUE, sep = "\t", row.names = 1)
shannon <- read.table(
  file.path(alpha_path, "shannon/alpha-diversity.tsv"),
  header = TRUE, sep = "\t", row.names = 1)
evenness <- read.table(
  file.path(alpha_path, "evenness/alpha-diversity.tsv"),
  header = TRUE, sep = "\t", row.names = 1)
observed <- read.table(
  file.path(alpha_path, "observed_features/alpha-diversity.tsv"),
  header = TRUE, sep = "\t", row.names = 1)

alpha_df <- metadata %>%
  rownames_to_column("SampleID") %>%
  left_join(faith_pd  %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  left_join(shannon   %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  left_join(evenness  %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  left_join(observed  %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  mutate(
    group      = factor(group,
                        levels = c("MSC","EV","Combo",
                                   "Positive","Negative","Normal")),
    timepoints = factor(timepoints,
                        levels = c("Day 0","Day 28","Day 56"))
  )

chao1_df <- estimate_richness(ps,
             measures = c("Chao1","Simpson","InvSimpson")) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = gsub("\\.", "-", SampleID))

alpha_df <- alpha_df %>% left_join(chao1_df, by = "SampleID")
cat("Alpha diversity data loaded:", nrow(alpha_df), "samples\n")

metrics <- list(
  list(col = "faith_pd",          label = "Faith Phylogenetic Diversity"),
  list(col = "shannon_entropy",   label = "Shannon Diversity Index"),
  list(col = "pielou_evenness",   label = "Pielou Evenness"),
  list(col = "observed_features", label = "Observed ASVs (Richness)"),
  list(col = "Chao1",             label = "Chao1 Richness Estimator"),
  list(col = "Simpson",           label = "Simpson Index"),
  list(col = "InvSimpson",        label = "Inverse Simpson Index")
)

cat("Running statistical tests...\n")
stats_results <- list()

for (m in metrics) {
  metric_col   <- m$col
  metric_label <- m$label

  # LME model
  lme_model <- tryCatch(
    lme(as.formula(paste0(metric_col, " ~ group * timepoints")),
        random = ~ 1 | subject,
        data = alpha_df, na.action = na.omit),
    error = function(e) NULL
  )

  # Kruskal-Wallis between groups per timepoint
  kw_results <- alpha_df %>%
    group_by(timepoints) %>%
    kruskal_test(as.formula(paste0(metric_col, " ~ group"))) %>%
    mutate(metric = metric_label)

  # Dunn post-hoc between groups per timepoint
  dunn_results <- alpha_df %>%
    group_by(timepoints) %>%
    dunn_test(
      as.formula(paste0(metric_col, " ~ group")),
      p.adjust.method = "BH"
    ) %>%
    mutate(metric = metric_label)

  # Wilcoxon within each group across timepoints
  # ALL pairs including ns
  wilcox_results <- alpha_df %>%
    group_by(group) %>%
    wilcox_test(
      as.formula(paste0(metric_col, " ~ timepoints")),
      p.adjust.method = "BH",
      paired = FALSE
    ) %>%
    add_xy_position(x = "timepoints") %>%
    mutate(metric = metric_label)

  stats_results[[metric_col]] <- list(
    lme    = lme_model,
    kw     = kw_results,
    dunn   = dunn_results,
    wilcox = wilcox_results
  )
}

cat("Statistical tests complete\n")

# Save stats CSVs
kw_all <- bind_rows(lapply(stats_results, function(x) x$kw))
write.csv(kw_all,
          file.path(plot_path, "kruskal_wallis_results.csv"),
          row.names = FALSE)

dunn_all <- bind_rows(lapply(stats_results, function(x) x$dunn))
write.csv(dunn_all,
          file.path(plot_path, "dunn_posthoc_results.csv"),
          row.names = FALSE)

# ============================================
# PLOT FUNCTION
# Faceted by GROUP — x-axis = timepoints
# Brackets = Wilcoxon within group (ALL shown)
# Subtitle = KW p-values between groups
# ============================================

plot_alpha <- function(metric_col, metric_label,
                        show_stats = FALSE) {

  plot_data <- alpha_df %>%
    filter(!is.na(.data[[metric_col]])) %>%
    mutate(
      group      = factor(group,
                          levels = c("MSC","EV","Combo",
                                     "Positive","Negative","Normal")),
      timepoints = factor(timepoints,
                          levels = c("Day 0","Day 28","Day 56"))
    )

  # KW p-values for subtitle
  kw <- stats_results[[metric_col]]$kw
  make_kw_label <- function(tp) {
    p   <- round(kw$p[kw$timepoints == tp], 3)
    sig <- ifelse(p < 0.05, "*", "")
    paste0(tp, ": p=", p, sig)
  }
  kw_subtitle <- paste(
    make_kw_label("Day 0"),
    make_kw_label("Day 28"),
    make_kw_label("Day 56"),
    sep = "  |  "
  )

  p <- ggplot(plot_data,
              aes(x     = timepoints,
                  y     = .data[[metric_col]],
                  fill  = group,
                  color = group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA,
                 linewidth = 0.6) +
    geom_jitter(width = 0.15, size = 1.8,
                alpha = 0.7, shape = 16) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~ group, nrow = 2, scales = "free_y") +
    labs(
      title    = metric_label,
      subtitle = kw_subtitle,
      x        = "Timepoint",
      y        = metric_label,
      caption  = paste0(
        "Pairwise Wilcoxon test with BH correction within each group\n",
        "Kruskal-Wallis p-values between all groups per timepoint\n",
        "ns: p>0.05  *: p<=0.05  **: p<=0.01  ***: p<=0.001"
      )
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9,
                                    colour = "grey40"),
      plot.caption  = element_text(size = 7.5, hjust = 0,
                                    colour = "grey50"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
      legend.position  = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )

  if (show_stats) {
    wilcox   <- stats_results[[metric_col]]$wilcox
    # Show ALL brackets including ns
    if (nrow(wilcox) > 0) {
      p <- p + stat_pvalue_manual(
        wilcox,
        label         = "p.adj.signif",
        tip.length    = 0.01,
        step.increase = 0.08,
        hide.ns       = FALSE
      )
    }
  }
  return(p)
}

cat("Generating alpha diversity plots...\n")

for (i in seq_along(metrics)) {
  m <- metrics[[i]]
  cat(i, "/", length(metrics), m$label, "\n")

  p_no_stats <- plot_alpha(m$col, m$label, show_stats = FALSE)
  ggsave(
    file.path(plot_path, "without_stats",
              paste0(m$col, "_boxplot.png")),
    p_no_stats, width = 14, height = 8, dpi = 300, bg = "white"
  )

  p_with_stats <- plot_alpha(m$col, m$label, show_stats = TRUE)
  ggsave(
    file.path(plot_path, "with_stats",
              paste0(m$col, "_boxplot.png")),
    p_with_stats, width = 14, height = 8, dpi = 300, bg = "white"
  )
}

# Combined figure
p_shannon <- plot_alpha("shannon_entropy",
                         "Shannon Diversity Index", TRUE)
p_faith   <- plot_alpha("faith_pd",
                         "Faith Phylogenetic Diversity", TRUE)
p_obs     <- plot_alpha("observed_features",
                         "Observed ASVs (Richness)", TRUE)

p_combined <- p_shannon / p_faith / p_obs +
  plot_annotation(
    title = "Alpha Diversity — All Treatment Groups Across Timepoints",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  )

ggsave(
  file.path(plot_path, "alpha_diversity_combined.png"),
  p_combined, width = 14, height = 18, dpi = 300, bg = "white"
)

cat("\n=== ALPHA DIVERSITY COMPLETE ===\n")
cat("Plots saved to:", plot_path, "\n")
cat("Without stats:", length(list.files(
  file.path(plot_path, "without_stats"))), "files\n")
cat("With stats:", length(list.files(
  file.path(plot_path, "with_stats"))), "files\n")
