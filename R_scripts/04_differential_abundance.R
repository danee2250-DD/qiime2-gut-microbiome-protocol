# ============================================
# 04_differential_abundance.R
# Differential Abundance Analysis (ANCOM-BC)
# Rat Gut Microbiome Study — 16S rRNA Analysis
# Author: Daneesha
# ============================================

# Run 00_setup.R first before this script
# source("00_setup.R")

# ============================================
# SET PATHS
# ============================================

ancombc_path <- file.path(data_path, "ancombc")
plot_path    <- file.path(plots_path, "differential_abundance")

# ============================================
# LOAD ANCOM-BC GENUS LEVEL RESULTS
# ============================================

timepoints_list <- c("Day_0", "Day_28", "Day_56")

load_ancombc_genus <- function(timepoint) {
  folder <- file.path(ancombc_path,
                       paste0("l6_", timepoint))
  files  <- list.files(folder, pattern = "\.tsv$",
                        full.names = TRUE)
  if (length(files) == 0) return(NULL)

  results <- lapply(files, function(f) {
    df <- read.table(f, header = TRUE,
                     sep = "	", check.names = FALSE)
    df$comparison <- gsub("_vs_.*|\.tsv", "",
                           basename(f))
    df$timepoint  <- gsub("_", " ", timepoint)
    df
  })
  bind_rows(results)
}

genus_results <- bind_rows(
  lapply(timepoints_list, load_ancombc_genus)
)

cat("ANCOM-BC genus results loaded\n")
cat("  Rows:", nrow(genus_results), "\n")
cat("  Columns:", colnames(genus_results), "\n")

# ============================================
# PARSE AND CLEAN GENUS NAMES
# ============================================

# Identify taxon column
taxon_col <- intersect(c("taxon","id","feature"),
                        colnames(genus_results))[1]

genus_results <- genus_results %>%
  mutate(
    genus_name = str_extract(.data[[taxon_col]],
                              "(?<=g__)[^;]+"),
    genus_name = ifelse(is.na(genus_name) | genus_name == "",
                         .data[[taxon_col]], genus_name),
    genus_name = gsub("_", " ", trimws(genus_name))
  )

# Identify LFC and q-value columns
lfc_col <- grep("^lfc", colnames(genus_results),
                 value = TRUE, ignore.case = TRUE)[1]
q_col   <- grep("^q|^p_adj|^padj",
                 colnames(genus_results),
                 value = TRUE, ignore.case = TRUE)[1]
p_col   <- grep("^p_val|^pval|^p$",
                 colnames(genus_results),
                 value = TRUE, ignore.case = TRUE)[1]

cat("LFC column:", lfc_col, "\n")
cat("Q column:", q_col, "\n")

genus_results <- genus_results %>%
  rename(lfc = all_of(lfc_col),
         q   = all_of(q_col)) %>%
  mutate(
    timepoint   = factor(timepoint,
                          levels = c("Day 0","Day 28","Day 56")),
    significant = !is.na(q) & q < 0.05,
    direction   = case_when(
      significant & lfc > 0 ~ "Enriched vs Positive",
      significant & lfc < 0 ~ "Depleted vs Positive",
      TRUE ~ "Not significant"
    )
  )

cat("Significant taxa:", sum(genus_results$significant), "\n")

# ============================================
# VOLCANO PLOT FUNCTION
# ============================================

plot_volcano <- function(df, timepoint_label) {

  plot_df <- df %>%
    filter(timepoint == timepoint_label) %>%
    mutate(
      neg_log_q  = -log10(q + 1e-10),
      label_name = ifelse(significant &
                             abs(lfc) > 1,
                           genus_name, "")
    )

  n_enriched <- sum(plot_df$direction ==
                      "Enriched vs Positive", na.rm = TRUE)
  n_depleted <- sum(plot_df$direction ==
                      "Depleted vs Positive", na.rm = TRUE)

  ggplot(plot_df,
         aes(x = lfc, y = neg_log_q,
             color = direction)) +
    geom_point(size = 2, alpha = 0.75) +
    ggrepel::geom_text_repel(
      aes(label = label_name),
      size = 2.8, max.overlaps = 15,
      show.legend = FALSE
    ) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dotted", color = "grey60") +
    scale_color_manual(
      values = c(
        "Enriched vs Positive" = "#E41A1C",
        "Depleted vs Positive" = "#377EB8",
        "Not significant"      = "grey70"
      ),
      name = "Direction"
    ) +
    facet_wrap(~ comparison, nrow = 2) +
    labs(
      title    = paste0("Differential Abundance — ",
                         timepoint_label,
                         " vs Positive Control"),
      subtitle = paste0("ANCOM-BC | BH FDR | Enriched: ",
                         n_enriched, " | Depleted: ",
                         n_depleted),
      x        = "Log Fold Change (LFC)",
      y        = "-log10(q-value)",
      caption  = "Dashed: q=0.05 | Dotted: |LFC|=0.5"
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(hjust = 0.5,
                                    face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9,
                                    color = "grey40"),
      plot.caption  = element_text(size = 7.5, hjust = 0),
      strip.background = element_rect(fill = "grey90"),
      strip.text    = element_text(face = "bold", size = 10),
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# ============================================
# GENERATE VOLCANO PLOTS
# ============================================

cat("Generating volcano plots...\n")

for (tp in c("Day 0","Day 28","Day 56")) {
  p_vol <- plot_volcano(genus_results, tp)
  ggsave(
    file.path(plot_path,
              paste0("volcano_genus_",
                      gsub(" ", "_", tp), ".png")),
    p_vol, width = 14, height = 10, dpi = 300
  )
  cat("Volcano saved:", tp, "\n")
}

# ============================================
# DA HEATMAP — DAY 56
# ============================================

cat("Generating DA heatmap Day 56...\n")

heatmap_df <- genus_results %>%
  filter(timepoint == "Day 56",
         significant == TRUE,
         !is.na(genus_name)) %>%
  dplyr::select(genus_name, comparison, lfc) %>%
  group_by(genus_name, comparison) %>%
  summarise(lfc = mean(lfc), .groups = "drop") %>%
  pivot_wider(names_from  = comparison,
              values_from = lfc,
              values_fill = 0) %>%
  column_to_rownames("genus_name")

if (nrow(heatmap_df) > 0) {
  p_heatmap <- pheatmap(
    as.matrix(heatmap_df),
    cluster_rows  = TRUE,
    cluster_cols  = FALSE,
    color = colorRampPalette(
      c("#2166AC","white","#B2182B"))(100),
    breaks       = seq(-3, 3, length.out = 101),
    main         = "Differential Abundance — Day 56 vs Positive Control\n(ANCOM-BC, q<0.05)",
    fontsize      = 9,
    fontsize_row  = 7.5,
    fontsize_col  = 9,
    border_color  = "white",
    angle_col     = 45,
    silent        = TRUE
  )

  png(file.path(plot_path, "da_heatmap_day56.png"),
      width = 12, height = 10, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(p_heatmap$gtable)
  dev.off()
  cat("DA heatmap Day 56 saved\n")
}

# ============================================
# BUBBLE PLOT — ALL TIMEPOINTS
# ============================================

cat("Generating bubble plot...\n")

bubble_df <- genus_results %>%
  filter(significant == TRUE,
         !is.na(genus_name)) %>%
  mutate(abs_lfc = abs(lfc),
         neg_log_q = -log10(q + 1e-10),
         comparison = factor(comparison,
                              levels = c("MSC","EV","Combo",
                                          "Negative","Normal")))

if (nrow(bubble_df) > 0) {
  p_bubble <- ggplot(bubble_df,
                      aes(x = comparison,
                          y = reorder(genus_name, lfc),
                          size = abs_lfc,
                          color = lfc,
                          shape = timepoint)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      name     = "LFC"
    ) +
    scale_size_continuous(
      range = c(1.5, 8),
      name  = "|LFC|"
    ) +
    scale_shape_manual(
      values = c("Day 0" = 16,
                  "Day 28" = 17,
                  "Day 56" = 15),
      name = "Timepoint"
    ) +
    labs(
      title    = "Significant Differential Taxa vs Positive Control",
      subtitle = "Size = effect size | Color = direction | Shape = timepoint",
      x        = "Treatment Group",
      y        = "Genus"
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(hjust = 0.5,
                                    face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9,
                                    color = "grey40"),
      axis.text.y   = element_text(size = 7.5),
      axis.text.x   = element_text(size = 10),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(plot_path,
                    "da_bubble_all_timepoints.png"),
         p_bubble, width = 12, height = 12, dpi = 300)
  cat("Bubble plot saved\n")
}

# ============================================
# SAVE SIGNIFICANT TAXA SUMMARY
# ============================================

sig_summary <- genus_results %>%
  filter(significant == TRUE) %>%
  dplyr::select(genus_name, comparison,
                timepoint, lfc, q, direction) %>%
  arrange(timepoint, comparison, q)

write.csv(sig_summary,
          file.path(plot_path,
                     "significant_taxa_all_timepoints.csv"),
          row.names = FALSE)

sig_day56 <- sig_summary %>%
  filter(timepoint == "Day 56")

write.csv(sig_day56,
          file.path(plot_path,
                     "significant_taxa_day56.csv"),
          row.names = FALSE)

cat("Significant taxa CSVs saved\n")
cat("  All timepoints:", nrow(sig_summary), "\n")
cat("  Day 56 only:", nrow(sig_day56), "\n")

cat("\n=== DIFFERENTIAL ABUNDANCE COMPLETE ===\n")
cat("Plots saved to:", plot_path, "\n")
cat("Files:", length(list.files(plot_path)), "\n")

