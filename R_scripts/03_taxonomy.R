# ============================================
# 03_taxonomy.R
# Taxonomic Composition Analysis
# Rat Gut Microbiome Study — 16S rRNA Analysis
# Author: Daneesha
# ============================================

# Run 00_setup.R first before this script
# source("00_setup.R")

# ============================================
# SET PATHS
# ============================================

plot_path <- file.path(plots_path, "taxonomy")

# ============================================
# PREPARE TAXONOMY DATA
# ============================================

# Add metadata to phyloseq sample data
sample_data(ps)$group <- factor(
  gsub("-.*", "", as.character(sample_data(ps)$treatment)),
  levels = c("MSC","EV","Combo","Positive","Negative","Normal")
)
sample_data(ps)$timepoints <- factor(
  sample_data(ps)$timepoints,
  levels = c("Day 0","Day 28","Day 56")
)

cat("Phyloseq object ready\n")
cat("  ASVs:", ntaxa(ps), "\n")
cat("  Samples:", nsamples(ps), "\n")

# ============================================
# PHYLUM LEVEL ANALYSIS
# ============================================

# Agglomerate to Phylum level
ps_phylum <- tax_glom(ps, "Phylum")
ps_phylum_rel <- transform_sample_counts(
  ps_phylum, function(x) x / sum(x)
)

# Melt to long format
phylum_long <- psmelt(ps_phylum_rel) %>%
  filter(!grepl("Unclassified", Phylum, ignore.case = TRUE)) %>%
  mutate(
    group = factor(group,
                   levels = c("MSC","EV","Combo",
                               "Positive","Negative","Normal")),
    timepoints = factor(timepoints,
                        levels = c("Day 0","Day 28","Day 56"))
  )

# Top 10 phyla
top_phyla <- phylum_long %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

phylum_long_top <- phylum_long %>%
  mutate(Phylum = ifelse(Phylum %in% top_phyla,
                          Phylum, "Other"))

# Phylum color palette
phylum_colors <- setNames(
  c(brewer.pal(10, "Paired"), "grey70"),
  c(top_phyla, "Other")
)

# Phylum bar plot by group
phylum_group_summary <- phylum_long_top %>%
  group_by(group, Phylum) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

p_phylum_group <- ggplot(phylum_group_summary,
                          aes(x = group, y = mean_abund,
                              fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill",
           width = 0.8, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = phylum_colors, name = "Phylum") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    title = "Phylum-Level Taxonomic Composition — by Treatment Group",
    x = "Treatment Group", y = "Relative Abundance (%)"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 9),
    legend.text     = element_text(size = 8),
    panel.grid      = element_blank()
  )

ggsave(file.path(plot_path, "phylum_by_group.png"),
       p_phylum_group, width = 12, height = 7, dpi = 300)
cat("1/7 Phylum by group saved\n")

# ============================================
# FAMILY LEVEL ANALYSIS
# ============================================

ps_family <- tax_glom(ps, "Family")
ps_family_rel <- transform_sample_counts(
  ps_family, function(x) x / sum(x)
)

family_long <- psmelt(ps_family_rel) %>%
  filter(!grepl("Unclassified", Family, ignore.case = TRUE)) %>%
  mutate(
    group = factor(group,
                   levels = c("MSC","EV","Combo",
                               "Positive","Negative","Normal")),
    timepoints = factor(timepoints,
                        levels = c("Day 0","Day 28","Day 56"))
  )

top_families <- family_long %>%
  group_by(Family) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 15) %>%
  pull(Family)

family_long_top <- family_long %>%
  mutate(Family = ifelse(Family %in% top_families,
                          Family, "Other"))

family_colors <- setNames(
  c(colorRampPalette(brewer.pal(12, "Set3"))(15), "grey70"),
  c(top_families, "Other")
)

family_group_summary <- family_long_top %>%
  group_by(group, Family) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

p_family_group <- ggplot(family_group_summary,
                          aes(x = group, y = mean_abund,
                              fill = Family)) +
  geom_bar(stat = "identity", position = "fill",
           width = 0.8, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = family_colors, name = "Family") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    title = "Family-Level Taxonomic Composition — by Treatment Group",
    x = "Treatment Group", y = "Relative Abundance (%)"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 9),
    legend.text     = element_text(size = 7.5),
    legend.key.size = unit(0.4, "cm"),
    panel.grid      = element_blank()
  )

ggsave(file.path(plot_path, "family_by_group.png"),
       p_family_group, width = 12, height = 7, dpi = 300)
cat("2/7 Family by group saved\n")

# ============================================
# GENUS LEVEL ANALYSIS
# ============================================

ps_genus <- tax_glom(ps, "Genus")
ps_genus_rel <- transform_sample_counts(
  ps_genus, function(x) x / sum(x)
)

genus_long <- psmelt(ps_genus_rel) %>%
  filter(!grepl("Unclassified|uncultured", Genus,
                ignore.case = TRUE),
         !is.na(Genus)) %>%
  mutate(
    group = factor(group,
                   levels = c("MSC","EV","Combo",
                               "Positive","Negative","Normal")),
    timepoints = factor(timepoints,
                        levels = c("Day 0","Day 28","Day 56"))
  )

top_genera <- genus_long %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

genus_long_top <- genus_long %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other"))

genus_colors <- setNames(
  c(colorRampPalette(c(
    "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
    "#A65628","#F781BF","#66C2A5","#FC8D62","#8DA0CB",
    "#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"
  ))(15), "grey70"),
  c(top_genera, "Other")
)

genus_group_summary <- genus_long_top %>%
  group_by(group, Genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

p_genus_group <- ggplot(genus_group_summary,
                         aes(x = group, y = mean_abund,
                             fill = Genus)) +
  geom_bar(stat = "identity", position = "fill",
           width = 0.8, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = genus_colors, name = "Genus") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    title = "Genus-Level Taxonomic Composition — by Treatment Group",
    x = "Treatment Group", y = "Relative Abundance (%)"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 9),
    legend.text     = element_text(size = 7.5),
    legend.key.size = unit(0.4, "cm"),
    panel.grid      = element_blank()
  )

ggsave(file.path(plot_path, "genus_by_group.png"),
       p_genus_group, width = 12, height = 7, dpi = 300)
cat("3/7 Genus by group saved\n")

# Genus facet by group and timepoint
genus_timepoint_summary <- genus_long_top %>%
  group_by(group, timepoints, Genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

p_genus_timepoint <- ggplot(genus_timepoint_summary,
                             aes(x = timepoints, y = mean_abund,
                                 fill = Genus)) +
  geom_bar(stat = "identity", position = "fill",
           width = 0.8, color = "white", linewidth = 0.08) +
  scale_fill_manual(values = genus_colors, name = "Genus") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  facet_wrap(~ group, nrow = 2) +
  labs(
    title = "Genus-Level Taxonomic Composition — by Group and Timepoint",
    x = "Timepoint", y = "Relative Abundance (%)"
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 8),
    legend.text      = element_text(size = 6.5),
    legend.key.size  = unit(0.35, "cm"),
    panel.grid       = element_blank()
  )

ggsave(file.path(plot_path, "genus_by_group_timepoint.png"),
       p_genus_timepoint, width = 14, height = 8, dpi = 300)
cat("4/7 Genus by group+timepoint saved\n")

# ============================================
# GENUS HEATMAP — TOP 15 GENERA
# ============================================

genus_heatmap_mat <- genus_long %>%
  filter(Genus %in% top_genera) %>%
  group_by(group, Genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = group,
              values_from = mean_abund,
              values_fill = 0) %>%
  column_to_rownames("Genus")

p_genus_heatmap <- pheatmap(
  as.matrix(genus_heatmap_mat),
  scale        = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  main         = "Top 15 Genera — Relative Abundance Heatmap",
  fontsize      = 10,
  fontsize_row  = 8,
  fontsize_col  = 10,
  border_color  = "white",
  angle_col     = "45"
  
)

png(file.path(plot_path, "genus_heatmap.png"),
    width = 12, height = 9, units = "in", res = 300)

draw(p_genus_heatmap)
dev.off()
cat("5/7 Genus heatmap saved\n")

# ============================================
# FIRMICUTES TO BACTEROIDOTA RATIO
# ============================================

fb_ratio <- phylum_long %>%
  filter(Phylum %in% c("Firmicutes","Bacteroidota")) %>%
  group_by(Sample, group, timepoints, subject, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Phylum,
              values_from = Abundance,
              values_fill = 0) %>%
  mutate(
    FB_ratio = (Firmicutes + 1e-10) / (Bacteroidota + 1e-10)
  )

fb_summary <- fb_ratio %>%
  group_by(group, timepoints) %>%
  summarise(
    mean_ratio = mean(FB_ratio, na.rm = TRUE),
    se_ratio   = sd(FB_ratio, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

p_fb_ratio <- ggplot(fb_summary,
                      aes(x = timepoints, y = mean_ratio,
                          color = group, group = group)) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(shape = group), size = 3.5) +
  geom_errorbar(aes(ymin = mean_ratio - se_ratio,
                     ymax = mean_ratio + se_ratio),
                width = 0.15, linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = group_colors, name = "Treatment") +
  scale_shape_manual(values = group_shapes, name = "Treatment") +
  labs(
    title    = "Firmicutes to Bacteroidota (F:B) Ratio — Across Timepoints",
    subtitle = "Dashed line indicates F:B ratio = 1",
    x        = "Timepoint",
    y        = "F:B Ratio (mean ± SE)"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_path, "fb_ratio_longitudinal.png"),
       p_fb_ratio, width = 10, height = 6, dpi = 300)
cat("6/7 F:B ratio plot saved\n")

# Save F:B ratio data
write.csv(fb_ratio,
          file.path(plot_path, "fb_ratio_data.csv"),
          row.names = FALSE)

# ============================================
# SAVE GENUS ABUNDANCE DATA
# ============================================

write.csv(genus_long_top,
          file.path(plot_path, "genus_abundance_data.csv"),
          row.names = FALSE)
cat("7/7 Data CSVs saved\n")

cat("\n=== TAXONOMY ANALYSIS COMPLETE ===\n")
cat("Plots saved to:", plot_path, "\n")
cat("Files:", length(list.files(plot_path)), "\n")

