# ============================================
# 05_longitudinal.R
# Longitudinal Analysis
# Rat Gut Microbiome Study — 16S rRNA Analysis
# Author: Daneesha
# ============================================

# Run 00_setup.R first before this script
# source("00_setup.R")

# ============================================
# SET PATHS
# ============================================

plot_path  <- file.path(plots_path, "longitudinal")
alpha_path <- file.path(data_path, "alpha")

# ============================================
# LOAD AND PREPARE ALPHA DIVERSITY DATA
# ============================================

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

chao1_df <- estimate_richness(ps,
             measures = c("Chao1","Simpson","InvSimpson")) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = gsub("\\.", "-", SampleID))

alpha_long <- metadata %>%
  rownames_to_column("SampleID") %>%
  left_join(faith_pd  %>% rownames_to_column("SampleID"),
            by = "SampleID") %>%
  left_join(shannon   %>% rownames_to_column("SampleID"),
            by = "SampleID") %>%
  left_join(evenness  %>% rownames_to_column("SampleID"),
            by = "SampleID") %>%
  left_join(observed  %>% rownames_to_column("SampleID"),
            by = "SampleID") %>%
  left_join(chao1_df, by = "SampleID") %>%
  mutate(
    group = factor(group,
                   levels = c("MSC","EV","Combo",
                               "Positive","Negative","Normal")),
    timepoints = factor(timepoints,
                        levels = c("Day 0","Day 28","Day 56"))
  )

cat("Longitudinal data prepared:", nrow(alpha_long), "samples\n")

# ============================================
# LONGITUDINAL LINE PLOT FUNCTION
# ============================================

plot_longitudinal <- function(metric_col, metric_label) {

  # Summary statistics per group per timepoint
  summary_df <- alpha_long %>%
    filter(!is.na(.data[[metric_col]])) %>%
    group_by(group, timepoints) %>%
    summarise(
      mean_val = mean(.data[[metric_col]], na.rm = TRUE),
      se_val   = sd(.data[[metric_col]], na.rm = TRUE) /
                   sqrt(n()),
      .groups  = "drop"
    )

  # Individual rat trajectories
  individual_df <- alpha_long %>%
    filter(!is.na(.data[[metric_col]]),
           !is.na(subject))

  p <- ggplot(summary_df,
              aes(x = timepoints, y = mean_val,
                  color = group, group = group)) +

    # Individual rat lines (faint background)
    geom_line(
      data = individual_df,
      aes(x = timepoints,
          y = .data[[metric_col]],
          group = interaction(subject, group)),
      alpha = 0.15, linewidth = 0.4
    ) +

    # Mean line
    geom_line(linewidth = 1.3) +

    # Mean point with shape
    geom_point(
      aes(shape = group, fill = group),
      size = 4, stroke = 0.4
    ) +

    # Error bars
    geom_errorbar(
      aes(ymin = mean_val - se_val,
          ymax = mean_val + se_val),
      width = 0.12, linewidth = 0.7
    ) +

    scale_color_manual(values = group_colors,
                        name = "Treatment") +
    scale_fill_manual(values  = group_colors,
                       name = "Treatment") +
    scale_shape_manual(values = group_shapes,
                        name = "Treatment") +

    labs(
      title = metric_label,
      x     = "Timepoint",
      y     = metric_label
    ) +
    theme_bw() +
    theme(
      plot.title       = element_text(hjust = 0.5,
                                       face = "bold", size = 12),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )

  return(p)
}

# ============================================
# GENERATE LONGITUDINAL PLOTS — ALL METRICS
# ============================================

cat("Generating longitudinal plots...\n")

metrics <- list(
  list(col = "faith_pd",          label = "Faith Phylogenetic Diversity"),
  list(col = "shannon_entropy",   label = "Shannon Diversity Index"),
  list(col = "pielou_evenness",   label = "Pielou Evenness"),
  list(col = "observed_features", label = "Observed ASVs (Richness)"),
  list(col = "Chao1",             label = "Chao1 Richness"),
  list(col = "Simpson",           label = "Simpson Index"),
  list(col = "InvSimpson",        label = "Inverse Simpson Index")
)

for (i in seq_along(metrics)) {
  m <- metrics[[i]]
  p <- plot_longitudinal(m$col, m$label)
  ggsave(
    file.path(plot_path,
              paste0("longitudinal_", m$col, ".png")),
    p, width = 10, height = 6, dpi = 300
  )
  cat(i, "/", length(metrics), m$label, "saved\n")
}

# ============================================
# COMBINED MAIN THESIS FIGURE
# ============================================

cat("Generating combined thesis figure...\n")

p_shannon <- plot_longitudinal("shannon_entropy",
                                "Shannon Diversity Index")
p_faith   <- plot_longitudinal("faith_pd",
                                "Faith Phylogenetic Diversity")
p_obs     <- plot_longitudinal("observed_features",
                                "Observed ASVs (Richness)")

p_combined <- (p_shannon + p_faith + p_obs) +
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(
    title = "Alpha Diversity Longitudinal Trends — All Treatment Groups",
    theme = theme(
      plot.title = element_text(hjust = 0.5,
                                 face = "bold", size = 14)
    )
  )

ggsave(
  file.path(plot_path, "longitudinal_combined_main.png"),
  p_combined, width = 12, height = 14, dpi = 300
)
cat("Combined thesis figure saved\n")

# ============================================
# FIRMICUTES TO BACTEROIDOTA RATIO
# ============================================

cat("Generating F:B ratio longitudinal plot...\n")

ps_phylum <- tax_glom(ps, "Phylum")
ps_phylum_rel <- transform_sample_counts(
  ps_phylum, function(x) x / sum(x)
)

fb_long <- psmelt(ps_phylum_rel) %>%
  filter(Phylum %in% c("Firmicutes","Bacteroidota")) %>%
  mutate(
    group = factor(
      gsub("-.*", "", as.character(treatment)),
      levels = c("MSC","EV","Combo",
                  "Positive","Negative","Normal")
    ),
    timepoints = factor(timepoints,
                         levels = c("Day 0","Day 28","Day 56"))
  ) %>%
  group_by(Sample, group, timepoints, subject, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from  = Phylum,
              values_from = Abundance,
              values_fill = 0) %>%
  mutate(FB_ratio = (Firmicutes + 1e-10) /
                     (Bacteroidota + 1e-10))

fb_summary <- fb_long %>%
  group_by(group, timepoints) %>%
  summarise(
    mean_ratio = mean(FB_ratio, na.rm = TRUE),
    se_ratio   = sd(FB_ratio, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

p_fb <- ggplot(fb_summary,
                aes(x = timepoints, y = mean_ratio,
                    color = group, group = group)) +

  # Individual rat lines
  geom_line(
    data = fb_long,
    aes(x = timepoints, y = FB_ratio,
        group = interaction(subject, group)),
    alpha = 0.12, linewidth = 0.4
  ) +

  geom_line(linewidth = 1.3) +
  geom_point(aes(shape = group, fill = group),
             size = 4, stroke = 0.4) +
  geom_errorbar(
    aes(ymin = mean_ratio - se_ratio,
        ymax = mean_ratio + se_ratio),
    width = 0.12, linewidth = 0.7
  ) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "grey40", linewidth = 0.8) +

  scale_color_manual(values = group_colors, name = "Treatment") +
  scale_fill_manual(values  = group_colors, name = "Treatment") +
  scale_shape_manual(values = group_shapes, name = "Treatment") +

  labs(
    title    = "Firmicutes to Bacteroidota (F:B) Ratio — Across Timepoints",
    subtitle = "Dashed line indicates F:B ratio = 1 (balanced community)",
    x        = "Timepoint",
    y        = "F:B Ratio (mean ± SE)"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5,
                                  face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 9,
                                  color = "grey40"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_path, "fb_ratio_longitudinal.png"),
       p_fb, width = 10, height = 6, dpi = 300)
cat("F:B ratio longitudinal plot saved\n")

# ============================================
# ============================================
# GENUS LONGITUDINAL — TOP 5, 6, 10 GENERA
# ============================================

cat("Generating genus longitudinal plots...\n")

ps_genus <- tax_glom(ps, "Genus")
ps_genus_rel <- transform_sample_counts(
  ps_genus, function(x) x / sum(x)
)

genus_long <- psmelt(ps_genus_rel) %>%
  filter(!grepl("Unclassified|uncultured",
                Genus, ignore.case = TRUE),
         !is.na(Genus)) %>%
  mutate(
    group = factor(
      gsub("-.*", "", as.character(treatment)),
      levels = c("MSC","EV","Combo",
                  "Positive","Negative","Normal")
    ),
    timepoints = factor(timepoints,
                         levels = c("Day 0","Day 28","Day 56"))
  )

# ── Plot 1: Faceted by GROUP, lines = genera ─────────────
plot_genus_by_group <- function(n_genera, nrow_facet) {

  top_genera <- genus_long %>%
    group_by(Genus) %>%
    summarise(mean_abund = mean(Abundance)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = n_genera) %>%
    pull(Genus)

  genus_sum <- genus_long %>%
    filter(Genus %in% top_genera) %>%
    group_by(Genus, group, timepoints) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(
      Genus = factor(Genus, levels = top_genera),
      group = factor(group,
                     levels = c("MSC","EV","Combo",
                                "Positive","Negative","Normal"))
    )

  genus_line_colors <- setNames(
    colorRampPalette(c(
      "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
      "#A65628","#F781BF","#66C2A5","#FC8D62","#8DA0CB"
    ))(n_genera),
    top_genera
  )

  ggplot(genus_sum,
         aes(x = timepoints, y = mean_abund,
             color = Genus, group = Genus)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3, shape = 16) +
    scale_color_manual(values = genus_line_colors, name = "Genus") +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, NA)) +
    facet_wrap(~ group, nrow = nrow_facet, scales = "fixed") +
    labs(
      title = paste0("Top ", n_genera,
                     " Genera — Longitudinal Change by Treatment Group"),
      x = "Timepoint",
      y = "Mean Relative Abundance"
    ) +
    theme_bw() +
    theme(
      plot.title       = element_text(hjust = 0.5,
                                      face = "bold", size = 13),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold", size = 10),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(face = "italic", size = 8.5),
      panel.grid.minor = element_blank()
    )
}

# ── Plot 2: Faceted by GENUS, lines = treatment groups ───
plot_genus_by_genus <- function(n_genera, nrow_facet) {

  top_genera <- genus_long %>%
    group_by(Genus) %>%
    summarise(mean_abund = mean(Abundance)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = n_genera) %>%
    pull(Genus)

  genus_sum <- genus_long %>%
    filter(Genus %in% top_genera) %>%
    group_by(Genus, group, timepoints) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(
      Genus = factor(Genus, levels = top_genera),
      group = factor(group,
                     levels = c("MSC","EV","Combo",
                                "Positive","Negative","Normal"))
    )

  ggplot(genus_sum,
         aes(x = timepoints, y = mean_abund,
             color = group, group = group)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3, shape = 16) +
    scale_color_manual(values = group_colors, name = "Treatment") +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(~ Genus, nrow = nrow_facet,
               scales = "free_y",
               labeller = label_wrap_gen(width = 20)) +
    labs(
      title = paste0("Top ", n_genera,
                     " Genera — Longitudinal Abundance Trends"),
      x = "Timepoint",
      y = "Relative Abundance"
    ) +
    theme_bw() +
    theme(
      plot.title       = element_text(hjust = 0.5,
                                      face = "bold", size = 13),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold.italic", size = 9),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

# ── Generate all 6 plots ─────────────────────────────────

# Style 1: Faceted by GROUP
p_grp_top5 <- plot_genus_by_group(n_genera = 5, nrow_facet = 1)
ggsave(file.path(plot_path, "genus_longitudinal_top5_bygroup.png"),
       p_grp_top5, width = 18, height = 5, dpi = 300, bg = "white")
cat("By-group Top 5 saved\n")

p_grp_top6 <- plot_genus_by_group(n_genera = 6, nrow_facet = 2)
ggsave(file.path(plot_path, "genus_longitudinal_top6_bygroup.png"),
       p_grp_top6, width = 14, height = 8, dpi = 300, bg = "white")
cat("By-group Top 6 saved\n")

p_grp_top10 <- plot_genus_by_group(n_genera = 10, nrow_facet = 2)
ggsave(file.path(plot_path, "genus_longitudinal_top10_bygroup.png"),
       p_grp_top10, width = 18, height = 9, dpi = 300, bg = "white")
cat("By-group Top 10 saved\n")

# Style 2: Faceted by GENUS
p_gen_top5 <- plot_genus_by_genus(n_genera = 5, nrow_facet = 1)
ggsave(file.path(plot_path, "genus_longitudinal_top5_bygenus.png"),
       p_gen_top5, width = 18, height = 5, dpi = 300, bg = "white")
cat("By-genus Top 5 saved\n")

p_gen_top6 <- plot_genus_by_genus(n_genera = 6, nrow_facet = 2)
ggsave(file.path(plot_path, "genus_longitudinal_top6_bygenus.png"),
       p_gen_top6, width = 14, height = 8, dpi = 300, bg = "white")
cat("By-genus Top 6 saved\n")

p_gen_top10 <- plot_genus_by_genus(n_genera = 10, nrow_facet = 2)
ggsave(file.path(plot_path, "genus_longitudinal_top10_bygenus.png"),
       p_gen_top10, width = 20, height = 9, dpi = 300, bg = "white")
cat("By-genus Top 10 saved\n")

# ============================================
# SAVE DATA CSVs
# ============================================

write.csv(fb_long,
          file.path(plot_path, "fb_ratio_data.csv"),
          row.names = FALSE)

genus_summary <- genus_long %>%
  filter(Genus %in% (genus_long %>%
    group_by(Genus) %>%
    summarise(mean_abund = mean(Abundance)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = 10) %>%
    pull(Genus))) %>%
  group_by(Genus, group, timepoints) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE),
            .groups = "drop")

write.csv(genus_summary,
          file.path(plot_path, "genus_abundance_longitudinal.csv"),
          row.names = FALSE)

cat("\n=== LONGITUDINAL ANALYSIS COMPLETE ===\n")
cat("Plots saved to:", plot_path, "\n")
cat("Files:", length(list.files(plot_path)), "\n")
