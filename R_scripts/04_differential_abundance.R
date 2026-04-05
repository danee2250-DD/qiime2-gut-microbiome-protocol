# ============================================
# 04_differential_abundance.R
# Differential Abundance Analysis (ANCOM-BC)
# Rat Gut Microbiome Study — 16S rRNA Analysis
# Author: Daneesha
# Last Updated: 2026
# ============================================

# Run 00_setup.R first
# source("00_setup.R")

# ============================================
# PACKAGES
# ============================================

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")

library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================
# SET PATHS
# ============================================

ancombc_path <- file.path(data_path, "ancombc")
plot_path    <- file.path(plots_path, "differential_abundance")
dir.create(plot_path, showWarnings = FALSE, recursive = TRUE)

# ============================================
# LOAD ANCOM-BC GENUS LEVEL RESULTS
# ============================================

timepoints_list <- c("Day_0", "Day_28", "Day_56")

load_ancombc_genus <- function(timepoint) {
  folder  <- file.path(ancombc_path, paste0("l6_", timepoint))
  lfc_f   <- file.path(folder, "lfc_slice.csv")
  q_f     <- file.path(folder, "q_val_slice.csv")
  if (!file.exists(lfc_f) | !file.exists(q_f)) return(NULL)

  lfc_wide <- read.csv(lfc_f, check.names = FALSE)
  q_wide   <- read.csv(q_f,   check.names = FALSE)

  lfc_long <- lfc_wide %>%
    dplyr::select(id, starts_with("treatment")) %>%
    pivot_longer(-id, names_to = "comparison", values_to = "lfc")

  q_long <- q_wide %>%
    dplyr::select(id, starts_with("treatment")) %>%
    pivot_longer(-id, names_to = "comparison", values_to = "q")

  combined <- left_join(lfc_long, q_long, by = c("id", "comparison")) %>%
    mutate(
      comparison = gsub("^treatment", "", comparison),
      comparison = gsub("-\\d+$", "", comparison),
      timepoint  = gsub("_", " ", timepoint)
    )
  combined
}

genus_results <- bind_rows(lapply(timepoints_list, load_ancombc_genus))

cat("ANCOM-BC genus results loaded\n")
cat("  Rows:", nrow(genus_results), "\n")
cat("  Comparisons:", paste(unique(genus_results$comparison), collapse=", "), "\n")
cat("  Timepoints:", paste(unique(genus_results$timepoint), collapse=", "), "\n")

# ============================================
# CLEAN GENUS NAMES — STRICT GENUS ONLY
# ============================================

clean_genus_name <- function(tax_string) {
  # Extract genus level (g__)
  genus <- str_extract(tax_string, "(?<=g__)[^;]+")
  genus <- trimws(gsub("_", " ", genus))
  # If genus empty or whitespace after g__
  if (is.na(genus) || genus == "" || genus == " ") {
    return(NA_character_)
  }
  # Remove square brackets
  genus <- gsub("\\[|\\]", "", genus)
  # Hide uninformative names
  hide_names <- c("uncultured", "Incertae Sedis",
                  "uncultured bacterium", "metagenome", "", " ")
  if (genus %in% hide_names) return(NA_character_)
  return(genus)
}
genus_results <- genus_results %>%
  mutate(genus_name = sapply(id, clean_genus_name))

genus_results <- genus_results %>%
  mutate(
    timepoint   = factor(timepoint,
                         levels = c("Day 0", "Day 28", "Day 56")),
    significant = !is.na(q) & q < 0.05,
    direction   = case_when(
      significant & lfc >  0 ~ "Enriched vs Positive",
      significant & lfc <  0 ~ "Depleted vs Positive",
      TRUE                   ~ "Not significant"
    )
  )

genus_results_clean <- genus_results %>%
  filter(
    !is.na(genus_name),
    !grepl("\\(f\\)$", genus_name),
    !grepl(";", genus_name),
    !genus_name %in% c("uncultured", "Incertae Sedis",
                       "uncultured bacterium", "metagenome",
                       "d  Bacteria", "d Bacteria")
  )

cat("Significant taxa:", sum(genus_results$significant), "\n")
cat("Clean genus results:", nrow(genus_results_clean), "rows\n")

# ============================================
# VOLCANO PLOT FUNCTION
# ============================================

plot_volcano <- function(df, timepoint_label) {
  plot_df <- df %>%
    filter(timepoint == timepoint_label) %>%
    mutate(
      neg_log_q  = -log10(q + 1e-10),
      label_name = ifelse(significant & !is.na(genus_name),
                          genus_name, NA_character_),
      direction  = factor(direction,
                          levels = c("Not significant",
                                     "Depleted vs Positive",
                                     "Enriched vs Positive"))
    )

  counts_df <- plot_df %>%
    group_by(comparison) %>%
    summarise(
      n_enriched = sum(direction == "Enriched vs Positive", na.rm = TRUE),
      n_depleted = sum(direction == "Depleted vs Positive", na.rm = TRUE),
      x_pos      = max(abs(lfc), na.rm = TRUE) * 0.95,
      y_pos      = max(neg_log_q, na.rm = TRUE) * 0.97,
      .groups    = "drop"
    ) %>%
    mutate(count_label = paste0("Enriched: ", n_enriched,
                                "  Depleted: ", n_depleted))

  ggplot(plot_df, aes(x = lfc, y = neg_log_q)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey45", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted",
               colour = "grey55", linewidth = 0.4) +
    geom_point(data = ~ filter(.x, direction == "Not significant"),
               colour = "grey78", fill = "grey88",
               shape = 21, size = 1.8, alpha = 0.55, stroke = 0.15) +
    geom_point(data = ~ filter(.x, direction != "Not significant"),
               aes(fill = direction),
               colour = "white", shape = 21,
               size = 3.2, alpha = 0.92, stroke = 0.45) +
    geom_text_repel(
      aes(label = label_name, colour = direction),
      size = 2.7, fontface = "italic",
      box.padding = 0.5, point.padding = 0.35,
      segment.colour = "grey50", segment.size = 0.3,
      segment.alpha = 0.65, min.segment.length = 0.15,
      max.overlaps = Inf, seed = 42, na.rm = TRUE,
      force = 3, force_pull = 0.5, show.legend = FALSE
    ) +
    geom_text(data = counts_df,
              aes(x = x_pos, y = y_pos, label = count_label),
              inherit.aes = FALSE, hjust = 1, vjust = 1,
              size = 3, colour = "grey20", fontface = "bold") +
    scale_fill_manual(
      values = c("Enriched vs Positive" = "#C62828",
                 "Depleted vs Positive" = "#1565C0"),
      name  = "Direction",
      guide = guide_legend(
        override.aes = list(size = 4, stroke = 0.5, colour = "white")
      )
    ) +
    scale_colour_manual(
      values = c("Enriched vs Positive" = "#C62828",
                 "Depleted vs Positive" = "#1565C0",
                 "Not significant"      = "grey60"),
      guide = "none"
    ) +
    facet_wrap(~ comparison, nrow = 2) +
    labs(
      title   = paste0("Differential Abundance — ", timepoint_label,
                       " vs Positive Control"),
      x       = "Log Fold Change (LFC)",
      y       = expression(-log[10](italic(q)-value)),
      caption = "Dashed: q = 0.05  |  Dotted: |LFC| = 0.5  |  Unresolved taxonomy hidden"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.caption     = element_text(size = 7.5, hjust = 0, colour = "grey50"),
      strip.background = element_rect(fill = "grey92", colour = "grey70"),
      strip.text       = element_text(face = "bold", size = 10.5),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 9),
      legend.text      = element_text(size = 8.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey94")
    )
}

# ============================================
# GENERATE VOLCANO PLOTS
# ============================================

cat("Generating volcano plots...\n")
for (tp in c("Day 0", "Day 28", "Day 56")) {
  p_vol <- plot_volcano(genus_results, tp)
  ggsave(
    filename = file.path(plot_path,
                         paste0("volcano_genus_", gsub(" ", "_", tp), ".png")),
    plot = p_vol, width = 16, height = 11, dpi = 300, bg = "white"
  )
  cat("Volcano saved:", tp, "\n")
}

# ============================================
# DA HEATMAP — DAY 56 ONLY
# ============================================

cat("Generating DA heatmap Day 56...\n")

all_groups <- c("MSC", "EV", "Combo", "Negative", "Normal")

heatmap_df <- genus_results_clean %>%
  filter(timepoint == "Day 56", significant == TRUE) %>%
  dplyr::select(genus_name, comparison, lfc) %>%
  group_by(genus_name, comparison) %>%
  summarise(lfc = mean(lfc), .groups = "drop") %>%
  complete(genus_name, comparison = all_groups,
           fill = list(lfc = 0)) %>%
  pivot_wider(names_from  = comparison,
              values_from = lfc,
              values_fill = 0) %>%
  column_to_rownames("genus_name")

col_order  <- intersect(all_groups, colnames(heatmap_df))
heatmap_df <- heatmap_df[, col_order, drop = FALSE]

sig_df <- genus_results_clean %>%
  filter(timepoint == "Day 56") %>%
  dplyr::select(genus_name, comparison, q) %>%
  group_by(genus_name, comparison) %>%
  slice(1) %>% ungroup() %>%
  complete(genus_name, comparison = all_groups,
           fill = list(q = 1)) %>%
  pivot_wider(names_from  = comparison,
              values_from = q,
              values_fill = 1) %>%
  column_to_rownames("genus_name")

common_rows <- intersect(rownames(heatmap_df), rownames(sig_df))
sig_df      <- sig_df[common_rows, col_order, drop = FALSE]
heatmap_df  <- heatmap_df[common_rows, , drop = FALSE]

if (nrow(heatmap_df) > 0) {
  lfc_mat   <- as.matrix(heatmap_df)
  sig_mat   <- as.matrix(sig_df)
  lfc_range <- ceiling(max(abs(lfc_mat), na.rm = TRUE))

  col_fun <- colorRamp2(
    c(-lfc_range, -lfc_range/2, 0, lfc_range/2, lfc_range),
    c("#1565C0", "#90CAF9", "white", "#EF9A9A", "#C62828")
  )

  group_cols <- c(MSC="#2196F3", EV="#4CAF50", Combo="#FF9800",
                  Negative="#9C27B0", Normal="#4E342E")

  top_anno <- HeatmapAnnotation(
    Group = col_order,
    col   = list(Group = group_cols[col_order]),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 8, fontface = "bold"),
    simple_anno_size     = unit(0.5, "cm"),
    show_legend          = FALSE
  )

  cell_fn <- function(j, i, x, y, width, height, fill) {
    q_val  <- sig_mat[i, j]
    lfc_val <- lfc_mat[i, j]
    if (is.na(q_val) || q_val >= 0.05 || abs(lfc_val) < 0.01) {
      grid.rect(x, y, width, height,
                gp = gpar(fill = "grey95", col = "white", lwd = 0.5))
      return()
    }
    txt_col <- ifelse(abs(lfc_val) > 2, "white", "grey20")
    if (q_val < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 7, col = txt_col))
    } else if (q_val < 0.01) {
      grid.text("**",  x, y, gp = gpar(fontsize = 7, col = txt_col))
    } else {
      grid.points(x, y, pch = 16, size = unit(2.5, "pt"),
                  gp = gpar(col = txt_col))
    }
  }

  ht <- Heatmap(
    lfc_mat, name = "LFC", col = col_fun, cell_fun = cell_fn,
    top_annotation = top_anno,
    cluster_rows = TRUE, clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    row_dend_side = "left", row_dend_width = unit(2, "cm"),
    show_row_names = TRUE, row_names_side = "right",
    row_names_gp = gpar(fontsize = 9, fontface = "italic"),
    row_names_max_width = unit(7, "cm"),
    cluster_columns = FALSE, show_column_names = TRUE,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 11, fontface = "bold",
                            col = unname(group_cols[col_order])),
    column_names_rot = 0,
    width  = unit(9, "cm"),
    height = unit(max(0.55 * nrow(lfc_mat), 8), "cm"),
    rect_gp = gpar(col = "white", lwd = 1), border = TRUE,
    heatmap_legend_param = list(
      title = "Log Fold\nChange",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(4, "cm"), border = "grey40"
    )
  )

  sig_lgd <- Legend(
    labels = c("q < 0.001", "q < 0.01", "q < 0.05", "Not significant"),
    title = "Significance",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    graphics = list(
      function(x,y,w,h) grid.text("***", x, y, gp=gpar(fontsize=9, col="grey20")),
      function(x,y,w,h) grid.text("**",  x, y, gp=gpar(fontsize=9, col="grey20")),
      function(x,y,w,h) grid.points(x, y, pch=16, size=unit(4,"pt"),
                                     gp=gpar(col="grey20")),
      function(x,y,w,h) grid.rect(x, y, w*0.8, h*0.8,
                                   gp=gpar(fill="grey95", col="grey70"))
    )
  )

  png(file.path(plot_path, "da_heatmap_day56.png"),
      width = 3600, height = 2800, res = 300, bg = "white")
  draw(ht,
       annotation_legend_list = list(sig_lgd),
       padding = unit(c(12, 5, 5, 5), "mm"),
       column_title = "Differential Abundance — Day 56 vs Positive Control\n(ANCOM-BC2 | Genus level | Reference: Positive control)",
       column_title_gp = gpar(fontsize = 13, fontface = "bold"))
  dev.off()
  cat("DA heatmap Day 56 saved\n")
} else {
  cat("No significant taxa for heatmap\n")
}

# ============================================
# DA HEATMAP — DAY 0 VS DAY 56 SIDE BY SIDE
# ============================================

cat("Generating Day 0 vs Day 56 heatmap...\n")

sig_taxa_day56 <- genus_results_clean %>%
  filter(timepoint == "Day 56", significant == TRUE) %>%
  pull(genus_name) %>% unique()

build_lfc_matrix <- function(tp, taxa) {
  genus_results_clean %>%
    filter(timepoint == tp, genus_name %in% taxa) %>%
    dplyr::select(genus_name, comparison, lfc) %>%
    group_by(genus_name, comparison) %>%
    summarise(lfc = mean(lfc), .groups = "drop") %>%
    complete(genus_name, comparison = all_groups, fill = list(lfc = 0)) %>%
    pivot_wider(names_from = comparison, values_from = lfc,
                values_fill = 0) %>%
    column_to_rownames("genus_name") %>%
    .[, intersect(all_groups, colnames(.)), drop = FALSE]
}

build_sig_matrix <- function(tp, taxa) {
  genus_results_clean %>%
    filter(timepoint == tp, genus_name %in% taxa) %>%
    dplyr::select(genus_name, comparison, q) %>%
    group_by(genus_name, comparison) %>%
    slice(1) %>% ungroup() %>%
    complete(genus_name, comparison = all_groups, fill = list(q = 1)) %>%
    pivot_wider(names_from = comparison, values_from = q,
                values_fill = 1) %>%
    column_to_rownames("genus_name") %>%
    .[, intersect(all_groups, colnames(.)), drop = FALSE]
}

lfc_d0  <- as.matrix(build_lfc_matrix("Day 0",  sig_taxa_day56))
lfc_d56 <- as.matrix(build_lfc_matrix("Day 56", sig_taxa_day56))
sig_d0  <- as.matrix(build_sig_matrix("Day 0",  sig_taxa_day56))
sig_d56 <- as.matrix(build_sig_matrix("Day 56", sig_taxa_day56))

common_rows <- intersect(rownames(lfc_d0), rownames(lfc_d56))
lfc_d0  <- lfc_d0[common_rows, ]
lfc_d56 <- lfc_d56[common_rows, ]
sig_d0  <- sig_d0[common_rows, ]
sig_d56 <- sig_d56[common_rows, ]

lfc_range <- ceiling(max(abs(c(lfc_d0, lfc_d56)), na.rm = TRUE))
col_fun2  <- colorRamp2(
  c(-lfc_range, -lfc_range/2, 0, lfc_range/2, lfc_range),
  c("#1565C0", "#90CAF9", "white", "#EF9A9A", "#C62828")
)

group_cols2 <- c(MSC="#2196F3", EV="#4CAF50", Combo="#FF9800",
                 Negative="#9C27B0", Normal="#4E342E")

make_top_anno <- function() {
  HeatmapAnnotation(
    Group = all_groups,
    col   = list(Group = group_cols2[all_groups]),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 8, fontface = "bold"),
    simple_anno_size     = unit(0.5, "cm"),
    show_legend          = FALSE
  )
}

make_cell_fn <- function(sig_mat, lfc_mat) {
  function(j, i, x, y, width, height, fill) {
    q_val   <- sig_mat[i, j]
    lfc_val <- lfc_mat[i, j]
    if (is.na(q_val) || q_val >= 0.05 || abs(lfc_val) < 0.01) {
      grid.rect(x, y, width, height,
                gp = gpar(fill = "grey95", col = "white", lwd = 0.5))
      return()
    }
    txt_col <- ifelse(abs(lfc_val) > 2, "white", "grey20")
    if (q_val < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 7, col = txt_col))
    } else if (q_val < 0.01) {
      grid.text("**",  x, y, gp = gpar(fontsize = 7, col = txt_col))
    } else {
      grid.points(x, y, pch = 16, size = unit(2.5, "pt"),
                  gp = gpar(col = txt_col))
    }
  }
}

ht_d0 <- Heatmap(
  lfc_d0, name = "LFC ",
  col = col_fun2, cell_fun = make_cell_fn(sig_d0, lfc_d0),
  top_annotation = make_top_anno(),
  column_title = "Day 0 (Baseline)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold", col = "grey30"),
  cluster_rows = TRUE, clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_side = "left", row_dend_width = unit(2, "cm"),
  show_row_names = FALSE,
  cluster_columns = FALSE, show_column_names = FALSE,
  width = unit(7, "cm"),
  height = unit(max(0.55 * nrow(lfc_d0), 8), "cm"),
  rect_gp = gpar(col = "white", lwd = 0.8), border = TRUE,
  show_heatmap_legend = FALSE
)

ht_d56 <- Heatmap(
  lfc_d56, name = "Log Fold Change",
  col = col_fun2, cell_fun = make_cell_fn(sig_d56, lfc_d56),
  top_annotation = make_top_anno(),
  column_title = "Day 56 (Endpoint)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold", col = "grey30"),
  cluster_rows = FALSE,
  show_row_names = TRUE, row_names_side = "right",
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  row_names_max_width = unit(7, "cm"),
  cluster_columns = FALSE, show_column_names = FALSE,
  width = unit(7, "cm"),
  height = unit(max(0.55 * nrow(lfc_d56), 8), "cm"),
  rect_gp = gpar(col = "white", lwd = 0.8), border = TRUE,
  heatmap_legend_param = list(
    title = "Log Fold\nChange",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(4, "cm"), border = "grey40"
  )
)

sig_lgd2 <- Legend(
  labels = c("q < 0.001", "q < 0.01", "q < 0.05", "Not significant"),
  title = "Significance",
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  graphics = list(
    function(x,y,w,h) grid.text("***", x, y, gp=gpar(fontsize=9, col="grey20")),
    function(x,y,w,h) grid.text("**",  x, y, gp=gpar(fontsize=9, col="grey20")),
    function(x,y,w,h) grid.points(x, y, pch=16, size=unit(4,"pt"),
                                   gp=gpar(col="grey20")),
    function(x,y,w,h) grid.rect(x, y, w*0.8, h*0.8,
                                 gp=gpar(fill="grey95", col="grey70"))
  )
)

group_lgd <- Legend(
  labels = c("MSC", "EV", "Combo", "Negative", "Normal"),
  title = "Treatment Group",
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 9),
  legend_gp = gpar(fill = c("#2196F3","#4CAF50","#FF9800",
                             "#9C27B0","#4E342E"))
)

png(file.path(plot_path, "da_heatmap_day0_vs_day56.png"),
    width = 5000, height = 2800, res = 300, bg = "white")
draw(ht_d0 + ht_d56,
     annotation_legend_list = list(sig_lgd2, group_lgd),
     padding = unit(c(12, 5, 5, 5), "mm"),
     column_title = "Differential Abundance vs Positive Control — Day 0 (Baseline) vs Day 56 (Endpoint)\n(ANCOM-BC2 | Genus level | Taxa anchored on Day 56 significance)",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     gap = unit(1, "cm"), merge_legends = FALSE)
dev.off()
cat("Day 0 vs Day 56 heatmap saved\n")

# ============================================
# LOLLIPOP PLOTS — One per group, Day 56
# ============================================

cat("Generating lollipop plots...\n")

lollipop_day56 <- genus_results_clean %>%
  filter(significant == TRUE, timepoint == "Day 56") %>%
  mutate(
    comparison = factor(comparison,
                        levels = c("MSC","EV","Combo","Negative","Normal")),
    direction  = factor(direction,
                        levels = c("Depleted vs Positive",
                                   "Enriched vs Positive"))
  )

plot_lollipop_group <- function(df, group_name, group_colour) {
  group_df <- df %>%
    filter(comparison == group_name) %>%
    slice_min(q, n = 10, with_ties = FALSE) %>%
    mutate(genus_name = factor(genus_name,
                               levels = genus_name[order(lfc)]))
  if (nrow(group_df) == 0) {
    cat("  No significant taxa for", group_name, "— skipping\n")
    return(NULL)
  }
  n_enriched <- sum(group_df$direction == "Enriched vs Positive")
  n_depleted <- sum(group_df$direction == "Depleted vs Positive")

  ggplot(group_df, aes(x = lfc, y = genus_name, color = direction)) +
    geom_segment(aes(x = 0, xend = lfc, y = genus_name, yend = genus_name),
                 linewidth = 1, alpha = 0.75) +
    geom_point(aes(size = -log10(q + 1e-10)), shape = 16, alpha = 0.95) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.7) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted",
               colour = "grey60", linewidth = 0.4) +
    scale_color_manual(
      values = c("Enriched vs Positive" = "#C62828",
                 "Depleted vs Positive" = "#1565C0"),
      name = "Direction", drop = FALSE
    ) +
    scale_size_continuous(range = c(3, 8),
                          name = expression(-log[10](q))) +
    labs(
      title    = paste0(group_name, " vs Positive Control — Day 56"),
      subtitle = paste0("Enriched: ", n_enriched,
                        "  |  Depleted: ", n_depleted,
                        "  |  Top 10 by significance"),
      x = "Log Fold Change (LFC)", y = NULL,
      caption = "Dotted: |LFC| = 0.5  |  q < 0.05  |  Genus level only"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold",
                                   size = 13, colour = group_colour),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
      plot.caption  = element_text(size = 7.5, hjust = 0, colour = "grey50"),
      axis.text.y   = element_text(size = 10, face = "italic"),
      axis.text.x   = element_text(size = 9),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 9),
      legend.text      = element_text(size = 8.5),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey92")
    )
}

grp_colours <- c(MSC="#2196F3", EV="#4CAF50", Combo="#FF9800",
                 Negative="#9C27B0", Normal="#4E342E")

for (grp in c("MSC","EV","Combo","Negative","Normal")) {
  p <- plot_lollipop_group(lollipop_day56, grp, grp_colours[grp])
  if (!is.null(p)) {
    ggsave(
      filename = file.path(plot_path,
                           paste0("da_lollipop_", grp, "_day56.png")),
      plot = p, width = 10, height = 7, dpi = 300, bg = "white"
    )
    cat("  Lollipop saved:", grp, "\n")
  }
}

# ============================================
# BUBBLE PLOT — TAXA SIGNIFICANT AT 2+ TIMEPOINTS
# ============================================

cat("Generating bubble plot...\n")

taxa_2tp <- genus_results_clean %>%
  filter(significant == TRUE) %>%
  group_by(genus_name) %>%
  summarise(n_timepoints = n_distinct(timepoint), .groups = "drop") %>%
  filter(n_timepoints >= 2) %>%
  pull(genus_name)

bubble_df <- genus_results_clean %>%
  filter(significant == TRUE, genus_name %in% taxa_2tp) %>%
  mutate(
    abs_lfc    = abs(lfc),
    comparison = factor(comparison,
                        levels = c("MSC","EV","Combo","Negative","Normal"))
  )

if (nrow(bubble_df) > 0) {
  p_bubble <- ggplot(bubble_df,
                     aes(x = comparison,
                         y = reorder(genus_name, lfc),
                         size = abs_lfc, color = lfc,
                         shape = timepoint)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(low="#1565C0", mid="white", high="#C62828",
                          midpoint = 0, name = "LFC") +
    scale_size_continuous(range = c(2, 9), name = "|LFC|") +
    scale_shape_manual(
      values = c("Day 0"=16, "Day 28"=17, "Day 56"=15),
      name = "Timepoint"
    ) +
    labs(
      title    = "Taxa Significant at 2+ Timepoints vs Positive Control",
      subtitle = "Size = effect size  |  Colour = direction  |  Shape = timepoint",
      x = "Treatment Group", y = "Genus",
      caption = "Only taxa with q < 0.05 at 2 or more timepoints shown"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
      plot.caption  = element_text(size = 7.5, hjust = 0, colour = "grey50"),
      axis.text.y   = element_text(size = 9, face = "italic"),
      axis.text.x   = element_text(size = 10, face = "bold"),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey94")
    )
  ggsave(file.path(plot_path, "da_bubble_all_timepoints.png"),
         p_bubble, width = 12, height = 10, dpi = 300, bg = "white")
  cat("Bubble plot saved\n")
}

# ============================================
# SAVE SIGNIFICANT TAXA CSVs
# ============================================

sig_summary <- genus_results_clean %>%
  filter(significant == TRUE) %>%
  dplyr::select(genus_name, comparison, timepoint, lfc, q, direction) %>%
  arrange(timepoint, comparison, q)

write.csv(sig_summary,
          file.path(plot_path, "significant_taxa_all_timepoints.csv"),
          row.names = FALSE)

sig_day56 <- sig_summary %>% filter(timepoint == "Day 56")
write.csv(sig_day56,
          file.path(plot_path, "significant_taxa_day56.csv"),
          row.names = FALSE)

cat("Significant taxa CSVs saved\n")
cat("  All timepoints:", nrow(sig_summary), "\n")
cat("  Day 56 only:",    nrow(sig_day56),   "\n")
cat("\n=== DIFFERENTIAL ABUNDANCE COMPLETE ===\n")
cat("Plots saved to:", plot_path, "\n")
cat("Files:", length(list.files(plot_path)), "\n")
