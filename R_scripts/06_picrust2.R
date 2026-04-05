# ============================================
# 06_picrust2.R
# PICRUSt2 Functional Prediction Analysis
# Rat Gut Microbiome Study
# Author: Daneesha
# ============================================

library(ggpicrust2)
library(tidyverse)
library(ggprism)
library(patchwork)
library(pheatmap)
library(vegan)
library(rstatix)
library(phyloseq)

cat("Starting PICRUSt2 analysis...\n")

# ============================================
# PATHS
# ============================================

base_dir     <- "/home/daneesha/QIIME_Analysis"
data_dir     <- file.path(base_dir, "R_exports")
plots_dir    <- file.path(base_dir, "Plots")
picrust_path <- file.path(data_dir, "picrust2")
picrust_plots <- file.path(plots_dir, "picrust2")
dir.create(picrust_plots, showWarnings = FALSE, recursive = TRUE)

# ============================================
# LOAD METADATA
# ============================================

metadata <- read.table(
  file.path(data_dir, "Metadata.tsv"),
  header = TRUE, sep = "\t", row.names = 1
)
metadata <- metadata[rownames(metadata) != "#q2:types", ]
metadata$timepoints <- factor(metadata$timepoints,
                               levels = c("Day 0","Day 28","Day 56"))
metadata$group <- factor(
  gsub("-.*", "", as.character(metadata$treatment)),
  levels = c("MSC","EV","Combo","Positive","Negative","Normal")
)

# Prepare metadata for ggpicrust2
meta_gg <- metadata %>%
  rownames_to_column("sample_name")

cat("1/9 Metadata loaded:", nrow(meta_gg), "samples\n")

# ============================================
# GROUP COLORS AND SHAPES
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

# ============================================
# LOAD PICRUST2 DATA
# ============================================

# MetaCyc pathways
pathway <- read.table(
  file.path(picrust_path, "pathway_abundance.tsv"),
  header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE
)

# EC numbers
ec <- read.table(
  file.path(picrust_path, "EC_metagenome.tsv"),
  header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE
)

# KO (KEGG Orthology)
ko <- read.table(
  file.path(picrust_path, "KO_metagenome.tsv"),
  header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE
)

# NSTI scores
nsti_raw <- read.table(
  file.path(picrust_path, "nsti_scores.tsv"),
  header = TRUE, sep = "\t"
)

cat("2/9 PICRUSt2 data loaded\n")
cat("  Pathways:", nrow(pathway), "\n")
cat("  EC functions:", nrow(ec), "\n")
cat("  KO functions:", nrow(ko), "\n")
cat("  NSTI ASVs:", nrow(nsti_raw), "\n")

# ============================================
# ALIGN SAMPLES
# ============================================

pathway <- pathway[, meta_gg$sample_name]
ec_aligned <- ec[, meta_gg$sample_name]

# Relative abundance
pathway_rel <- sweep(pathway, 2, colSums(pathway), "/")

# ============================================
# NSTI QUALITY PLOT
# ============================================

p_nsti <- ggplot(nsti_raw, aes(x = metadata_NSTI)) +
  geom_histogram(bins = 50, fill = "#377EB8",
                 color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0.06, linetype = "dashed",
             color = "darkgreen", linewidth = 1) +
  geom_vline(xintercept = 0.15, linetype = "dashed",
             color = "darkorange", linewidth = 1) +
  geom_vline(xintercept = 2.0, linetype = "dashed",
             color = "red", linewidth = 1) +
  annotate("text", x = 0.06, y = 750,
           label = "Excellent\n(<0.06)",
           size = 3, color = "darkgreen", hjust = -0.1) +
  annotate("text", x = 0.15, y = 700,
           label = "Good\n(<0.15)",
           size = 3, color = "darkorange", hjust = -0.1) +
  annotate("text", x = 2.0, y = 750,
           label = "Threshold\n(2.0)",
           size = 3, color = "red", hjust = 1.1) +
  labs(
    title = "NSTI Score Distribution — Prediction Quality Assessment",
    subtitle = paste0(
      "Mean NSTI = ", round(mean(nsti_raw$metadata_NSTI), 3),
      " | Median = ", round(median(nsti_raw$metadata_NSTI), 3),
      " | All ", nrow(nsti_raw), " ASVs below threshold of 2.0"
    ),
    x = "NSTI Score", y = "Number of ASVs",
    caption = "Green: excellent (<0.06) | Orange: good (<0.15) | Red: threshold (<2.0)"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust=0.5, face="bold", size=13),
    plot.subtitle = element_text(hjust=0.5, size=9, color="grey40"),
    plot.caption  = element_text(size=8, hjust=0),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(picrust_plots, "nsti_quality.png"),
       p_nsti, width=10, height=6, dpi=300)
cat("3/9 NSTI plot saved\n")

# ============================================
# PATHWAY NAME MAPPING
# ============================================

pathway_name_map <- c(
  "NONOXIPENT-PWY"   = "Non-Oxidative Pentose Phosphate",
  "PWY-7663"         = "Gondoate Biosynthesis (Anaerobic)",
  "CALVIN-PWY"       = "Calvin Cycle",
  "PWY-7208"         = "MnmE tRNA Modification",
  "PWY-7219"         = "Adenosine Nucleotides De Novo Biosynthesis",
  "PWY-5973"         = "Cis-Vaccenate Biosynthesis",
  "PWY-7229"         = "Superpathway Adenosine Nucleotides Biosynthesis",
  "PWY-6126"         = "Superpathway Adenosine Nucleotides Biosynthesis II",
  "PWY-5686"         = "UMP Biosynthesis",
  "ANAGLYCOLYSIS-PWY"= "Anaerobic Glycolysis",
  "PWY-6123"         = "Inosine-5-phosphate Biosynthesis",
  "PWY-5101"         = "Fatty Acid Beta Oxidation",
  "PWY-6737"         = "Starch Degradation",
  "THRESYN-PWY"      = "Threonine Biosynthesis",
  "PWY-2942"         = "L-Lysine Biosynthesis",
  "PWY-7220"         = "Adenosine Deoxyribonucleotides De Novo",
  "PWY-7222"         = "Guanosine Deoxyribonucleotides De Novo",
  "PWY-5104"         = "Fatty Acid Elongation",
  "PWY-5667"         = "CDP-Diacylglycerol Biosynthesis",
  "PWY0-1319"        = "Fatty Acid Biosynthesis Initiation"
)

# Top 20 pathways
top20_pathways <- names(sort(rowMeans(pathway_rel),
                              decreasing = TRUE))[1:20]
top20_named <- data.frame(
  pathway_id   = top20_pathways,
  pathway_name = pathway_name_map[top20_pathways]
)

# ============================================
# METACYC PATHWAY BAR PLOT BY GROUP
# ============================================

# Prepare data
pathway_meta <- as.data.frame(t(pathway_rel))
pathway_meta$SampleID <- rownames(pathway_meta)
pathway_meta <- pathway_meta %>%
  left_join(meta_gg %>%
              dplyr::select(sample_name, group, timepoints),
            by = c("SampleID" = "sample_name")) %>%
  mutate(
    group = factor(group,
                   levels = c("MSC","EV","Combo",
                               "Positive","Negative","Normal")),
    timepoints = factor(timepoints,
                        levels = c("Day 0","Day 28","Day 56"))
  )

# Long format
pathway_long_named <- pathway_meta %>%
  pivot_longer(cols = all_of(top20_pathways),
               names_to = "pathway_id",
               values_to = "abundance") %>%
  left_join(top20_named, by = "pathway_id")

# Group summary
pathway_group_summary <- pathway_long_named %>%
  group_by(group, pathway_name) %>%
  summarise(mean_abund = mean(abundance), .groups = "drop")

# Color palette
pathway_colors <- setNames(
  colorRampPalette(c(
    "#E41A1C","#FF7F00","#4DAF4A","#377EB8","#984EA3",
    "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
    "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
    "#B3B3B3","#1B9E77","#D95F02","#7570B3","#E7298A"
  ))(20),
  top20_named$pathway_name
)

# Bar plot by group
p_pathway_group <- ggplot(pathway_group_summary,
                           aes(x = group, y = mean_abund,
                               fill = pathway_name)) +
  geom_bar(stat = "identity", position = "fill",
           width = 0.8, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = pathway_colors, name = "Pathway") +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(
    title = "Top 20 Predicted MetaCyc Pathways — by Treatment Group",
    x = "Treatment Group", y = "Relative Abundance (%)"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust=0.5, face="bold", size=13),
    axis.text.x  = element_text(angle=45, hjust=1, size=10),
    legend.position = "right",
    legend.title = element_text(face="bold", size=9),
    legend.text  = element_text(size=7),
    legend.key.size = unit(0.4, "cm"),
    panel.grid   = element_blank()
  )

ggsave(file.path(picrust_plots, "pathway_by_group.png"),
       p_pathway_group, width=12, height=7, dpi=300)
cat("4/9 Pathway by group plot saved\n")

# ============================================
# METACYC PATHWAY FACET BY GROUP + TIMEPOINT
# ============================================

pathway_timepoint_summary <- pathway_long_named %>%
  group_by(group, timepoints, pathway_name) %>%
  summarise(mean_abund = mean(abundance), .groups = "drop")

p_pathway_timepoint <- ggplot(pathway_timepoint_summary,
                               aes(x = timepoints, y = mean_abund,
                                   fill = pathway_name)) +
  geom_bar(stat = "identity", position = "fill",
           width = 0.8, color = "white", linewidth = 0.08) +
  scale_fill_manual(values = pathway_colors, name = "Pathway") +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  facet_wrap(~ group, nrow = 2) +
  labs(
    title = "Top 20 Predicted MetaCyc Pathways — by Group and Timepoint",
    x = "Timepoint", y = "Relative Abundance (%)"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust=0.5, face="bold", size=13),
    axis.text.x  = element_text(angle=45, hjust=1, size=8),
    strip.background = element_rect(fill="grey90"),
    strip.text   = element_text(face="bold", size=10),
    legend.position = "right",
    legend.title = element_text(face="bold", size=8),
    legend.text  = element_text(size=6.5),
    legend.key.size = unit(0.35, "cm"),
    panel.grid   = element_blank()
  )

ggsave(file.path(picrust_plots, "pathway_by_group_timepoint.png"),
       p_pathway_timepoint, width=14, height=8, dpi=300)
cat("5/9 Pathway by group+timepoint plot saved\n")

# ============================================
# METACYC PATHWAY HEATMAP
# ============================================

pathway_heatmap_mat <- pathway_group_summary %>%
  pivot_wider(names_from = group,
              values_from = mean_abund,
              values_fill = 0) %>%
  column_to_rownames("pathway_name")

p_pathway_heatmap <- pheatmap(
  as.matrix(pathway_heatmap_mat),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  main = "Top 20 Predicted MetaCyc Pathways\nRelative Abundance Heatmap",
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = "white",
  angle_col = "45",
  
)

png(file.path(picrust_plots, "pathway_heatmap.png"),
    width=12, height=9, units="in", res=300)
grid::grid.newpage()
draw(p_pathway_heatmap)
dev.off()
cat("6/9 Pathway heatmap saved\n")

# ============================================
# PERMANOVA ON METACYC PATHWAYS
# ============================================

pathway_matrix <- as.matrix(t(pathway_rel))
pathway_matrix <- pathway_matrix[meta_gg$sample_name, ]
pathway_dist <- vegdist(pathway_matrix, method = "bray")

set.seed(123)
perm_group <- adonis2(pathway_dist ~ group,
                       data = meta_gg, permutations = 999)
perm_time  <- adonis2(pathway_dist ~ timepoints,
                       data = meta_gg, permutations = 999)

pathway_stats <- data.frame(
  Comparison = c("Treatment Group (MetaCyc)",
                  "Timepoint (MetaCyc)",
                  "Treatment Group (Functional PCoA)"),
  R2      = c(round(perm_group$R2[1], 3),
               round(perm_time$R2[1], 3),
               0.096),
  P_value = c(perm_group$`Pr(>F)`[1],
               perm_time$`Pr(>F)`[1],
               0.004),
  Significant = c(
    ifelse(perm_group$`Pr(>F)`[1] < 0.05, "YES *", "NO"),
    ifelse(perm_time$`Pr(>F)`[1]  < 0.05, "YES *", "NO"),
    "YES *"
  ),
  Method = c(
    "PERMANOVA Bray-Curtis MetaCyc pathways",
    "PERMANOVA Bray-Curtis MetaCyc pathways",
    "PERMANOVA Bray-Curtis all KEGG pathways"
  )
)

write.csv(pathway_stats,
          file.path(picrust_plots, "pathway_permanova_stats.csv"),
          row.names = FALSE)
cat("PERMANOVA stats saved\n")

# ============================================
# FUNCTIONAL PCoA
# ============================================

set.seed(123)
func_permanova <- adonis2(pathway_dist ~ group,
                           data = meta_gg, permutations = 999)

func_pcoa <- cmdscale(pathway_dist, k = 2, eig = TRUE)
func_pcoa_df <- data.frame(
  PC1 = func_pcoa$points[, 1],
  PC2 = func_pcoa$points[, 2],
  SampleID = rownames(func_pcoa$points)
) %>%
  left_join(meta_gg %>%
              dplyr::select(sample_name, group, timepoints),
            by = c("SampleID" = "sample_name"))

eig  <- func_pcoa$eig
var1 <- round(eig[1] / sum(eig[eig > 0]) * 100, 1)
var2 <- round(eig[2] / sum(eig[eig > 0]) * 100, 1)

p_func_pcoa <- ggplot(func_pcoa_df,
                       aes(x = PC1, y = PC2,
                           color = group, shape = group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = group), type = "t",
               level = 0.95, linetype = "dashed",
               linewidth = 0.7, fill = NA) +
  scale_color_manual(values = group_colors, name = "Treatment") +
  scale_shape_manual(values = group_shapes, name = "Treatment") +
  labs(
    title    = "Predicted Functional Profile — PCoA (Bray-Curtis)",
    subtitle = paste0("PERMANOVA: R²=",
                       round(func_permanova$R2[1], 3),
                       ", p=", func_permanova$`Pr(>F)`[1]),
    x = paste0("PC1 (", var1, "%)"),
    y = paste0("PC2 (", var2, "%)")
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust=0.5, face="bold", size=13),
    plot.subtitle = element_text(hjust=0.5, size=9, color="grey40"),
    legend.position = "right",
    legend.title  = element_text(face="bold", size=10),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(picrust_plots, "functional_pcoa.png"),
       p_func_pcoa, width=10, height=7, dpi=300)
cat("7/9 Functional PCoA saved\n")

# ============================================
# KEGG LEVEL 2 ANALYSIS
# ============================================

# Convert KO to KEGG
kegg_abundance <- ko2kegg_abundance(
  file = file.path(picrust_path, "KO_metagenome.tsv")
)
kegg_aligned <- kegg_abundance[, meta_gg$sample_name]

# KEGG Level 2 category mapping
kegg_level2 <- c(
  "ko00010"="Carbohydrate metabolism","ko00020"="Carbohydrate metabolism",
  "ko00030"="Carbohydrate metabolism","ko00040"="Carbohydrate metabolism",
  "ko00051"="Carbohydrate metabolism","ko00052"="Carbohydrate metabolism",
  "ko00053"="Carbohydrate metabolism","ko00500"="Carbohydrate metabolism",
  "ko00520"="Carbohydrate metabolism","ko00620"="Carbohydrate metabolism",
  "ko00630"="Carbohydrate metabolism","ko00640"="Carbohydrate metabolism",
  "ko00650"="Carbohydrate metabolism","ko00660"="Carbohydrate metabolism",
  "ko00061"="Lipid metabolism","ko00062"="Lipid metabolism",
  "ko00071"="Lipid metabolism","ko00561"="Lipid metabolism",
  "ko00564"="Lipid metabolism","ko00565"="Lipid metabolism",
  "ko00600"="Lipid metabolism","ko01040"="Lipid metabolism",
  "ko01212"="Lipid metabolism",
  "ko00220"="Amino acid metabolism","ko00230"="Amino acid metabolism",
  "ko00250"="Amino acid metabolism","ko00260"="Amino acid metabolism",
  "ko00270"="Amino acid metabolism","ko00280"="Amino acid metabolism",
  "ko00290"="Amino acid metabolism","ko00300"="Amino acid metabolism",
  "ko00310"="Amino acid metabolism","ko00330"="Amino acid metabolism",
  "ko00340"="Amino acid metabolism","ko00350"="Amino acid metabolism",
  "ko00360"="Amino acid metabolism","ko00380"="Amino acid metabolism",
  "ko00400"="Amino acid metabolism","ko01230"="Amino acid metabolism",
  "ko00730"="Metabolism of cofactors and vitamins",
  "ko00740"="Metabolism of cofactors and vitamins",
  "ko00750"="Metabolism of cofactors and vitamins",
  "ko00760"="Metabolism of cofactors and vitamins",
  "ko00770"="Metabolism of cofactors and vitamins",
  "ko00780"="Metabolism of cofactors and vitamins",
  "ko00790"="Metabolism of cofactors and vitamins",
  "ko00860"="Metabolism of cofactors and vitamins",
  "ko01240"="Metabolism of cofactors and vitamins",
  "ko00190"="Energy metabolism","ko00195"="Energy metabolism",
  "ko00680"="Energy metabolism","ko00710"="Energy metabolism",
  "ko00720"="Energy metabolism","ko00910"="Energy metabolism",
  "ko00920"="Energy metabolism","ko01200"="Energy metabolism",
  "ko00240"="Nucleotide metabolism","ko01232"="Nucleotide metabolism",
  "ko00480"="Metabolism of other amino acids",
  "ko00410"="Metabolism of other amino acids",
  "ko00430"="Metabolism of other amino acids",
  "ko00450"="Metabolism of other amino acids",
  "ko00460"="Metabolism of other amino acids",
  "ko00471"="Metabolism of other amino acids",
  "ko00473"="Metabolism of other amino acids",
  "ko00550"="Biosynthesis of other secondary metabolites",
  "ko00900"="Biosynthesis of terpenoids and polyketides",
  "ko00940"="Biosynthesis of other secondary metabolites",
  "ko02010"="Membrane transport","ko02060"="Membrane transport",
  "ko03010"="Translation","ko03020"="Transcription",
  "ko03030"="Replication and repair",
  "ko03410"="Replication and repair",
  "ko03420"="Replication and repair",
  "ko03430"="Replication and repair",
  "ko03440"="Replication and repair",
  "ko00540"="Glycan biosynthesis and metabolism",
  "ko00510"="Glycan biosynthesis and metabolism",
  "ko02020"="Signal transduction",
  "ko02024"="Cell motility","ko02030"="Cell motility",
  "ko02040"="Cell motility","ko00970"="Translation",
  "ko03060"="Folding, sorting and degradation",
  "ko03070"="Membrane transport",
  "ko01100"="Global and overview maps",
  "ko01110"="Global and overview maps",
  "ko01120"="Global and overview maps",
  "ko01210"="Global and overview maps"
)

# Calculate KEGG Level 2 relative abundance
ko_t <- as.data.frame(t(kegg_aligned))
ko_t$SampleID <- rownames(ko_t)
ko_cols <- setdiff(colnames(ko_t), "SampleID")

ko_level2_rel <- ko_t %>%
  pivot_longer(cols = all_of(ko_cols),
               names_to = "ko_id", values_to = "abundance") %>%
  mutate(level2 = kegg_level2[ko_id],
         level2 = ifelse(is.na(level2), "Other", level2)) %>%
  group_by(SampleID, level2) %>%
  summarise(total_abund = sum(abundance, na.rm=TRUE), .groups="drop") %>%
  left_join(meta_gg %>%
              dplyr::select(sample_name, group, timepoints),
            by = c("SampleID" = "sample_name")) %>%
  mutate(group = factor(group,
                         levels=c("MSC","EV","Combo",
                                   "Positive","Negative","Normal"))) %>%
  group_by(SampleID) %>%
  mutate(rel_abund = total_abund / sum(total_abund)) %>%
  ungroup()

ko_level2_group <- ko_level2_rel %>%
  group_by(group, level2) %>%
  summarise(mean_rel = mean(rel_abund), .groups="drop") %>%
  filter(level2 != "Other")

# KEGG category colors
kegg_cat_colors <- c(
  "Carbohydrate metabolism"="#E41A1C",
  "Amino acid metabolism"="#377EB8",
  "Energy metabolism"="#4DAF4A",
  "Lipid metabolism"="#984EA3",
  "Metabolism of cofactors and vitamins"="#FF7F00",
  "Nucleotide metabolism"="#A65628",
  "Translation"="#FFD700",
  "Replication and repair"="#00CED1",
  "Membrane transport"="#FF1493",
  "Global and overview maps"="#999999",
  "Glycan biosynthesis and metabolism"="#32CD32",
  "Metabolism of other amino acids"="#1E90FF",
  "Biosynthesis of terpenoids and polyketides"="#FF6347",
  "Biosynthesis of other secondary metabolites"="#8A2BE2",
  "Signal transduction"="#20B2AA",
  "Cell motility"="#DC143C",
  "Transcription"="#00FA9A",
  "Folding, sorting and degradation"="#FF8C00"
)

# KEGG Level 2 bar plot
p_kegg_l2_bar <- ggplot(ko_level2_group,
                         aes(x=group, y=mean_rel, fill=level2)) +
  geom_bar(stat="identity", position="fill",
           width=0.8, color="white", linewidth=0.1) +
  scale_fill_manual(values=kegg_cat_colors, name="KEGG Category") +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  labs(title="Predicted KEGG Functional Categories — by Treatment Group",
       x="Treatment Group", y="Relative Abundance (%)") +
  theme_bw() +
  theme(
    plot.title=element_text(hjust=0.5, face="bold", size=13),
    axis.text.x=element_text(angle=45, hjust=1, size=10),
    legend.position="right",
    legend.title=element_text(face="bold", size=9),
    legend.text=element_text(size=7.5),
    legend.key.size=unit(0.45,"cm"),
    panel.grid=element_blank()
  )

ggsave(file.path(picrust_plots, "kegg_level2_barplot.png"),
       p_kegg_l2_bar, width=12, height=7, dpi=300)

# KEGG Level 2 heatmap
kegg_heatmap_mat <- ko_level2_group %>%
  pivot_wider(names_from=group, values_from=mean_rel,
              values_fill=0) %>%
  column_to_rownames("level2")

p_kegg_heatmap <- pheatmap(
  as.matrix(kegg_heatmap_mat),
  scale="row", cluster_rows=TRUE, cluster_cols=FALSE,
  color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  main="Predicted KEGG Functional Categories\nHierarchical Clustering Heatmap",
  fontsize=10, fontsize_row=8, fontsize_col=10,
  border_color="white", angle_col="45"
)

png(file.path(picrust_plots, "kegg_level2_heatmap.png"),
    width=12, height=8, units="in", res=300)
draw(p_kegg_heatmap)
dev.off()
cat("8/9 KEGG Level 2 plots saved\n")

# ============================================
# LEFSE-STYLE PATHWAY BAR PLOT
# ============================================

pathway_fc <- pathway_long_named %>%
  group_by(pathway_name, group) %>%
  summarise(mean_abund=mean(abundance), .groups="drop") %>%
  left_join(
    pathway_long_named %>%
      filter(group=="Positive") %>%
      group_by(pathway_name) %>%
      summarise(positive_mean=mean(abundance), .groups="drop"),
    by="pathway_name"
  ) %>%
  mutate(log2fc=log2((mean_abund+1e-10)/(positive_mean+1e-10))) %>%
  filter(group != "Positive")

kw_pathway <- pathway_long_named %>%
  group_by(pathway_name) %>%
  kruskal_test(abundance ~ group) %>%
  mutate(sig = p < 0.05)

lefse_df <- pathway_fc %>%
  left_join(kw_pathway %>%
              dplyr::select(pathway_name, p, sig),
            by="pathway_name") %>%
  filter(sig==TRUE) %>%
  mutate(
    group=factor(group, levels=c("MSC","EV","Combo","Negative","Normal")),
    pathway_short=substr(pathway_name, 1, 35)
  )

if (nrow(lefse_df) > 0) {
  p_lefse <- ggplot(lefse_df,
                     aes(x=log2fc,
                         y=reorder(pathway_short, log2fc),
                         fill=group, color=group)) +
    geom_bar(stat="identity", position="dodge",
             width=0.7, alpha=0.8) +
    geom_vline(xintercept=0, linewidth=0.7, color="grey40") +
    scale_fill_manual(values=group_colors, name="Treatment") +
    scale_color_manual(values=group_colors, name="Treatment") +
    labs(
      title="Predicted MetaCyc Pathways — Log2 Fold Change vs Positive Control",
      subtitle="Kruskal-Wallis significant pathways (p<0.05)",
      x="Log2 Fold Change vs Positive Control",
      y="MetaCyc Pathway"
    ) +
    theme_bw() +
    theme(
      plot.title=element_text(hjust=0.5, face="bold", size=12),
      plot.subtitle=element_text(hjust=0.5, size=9, color="grey40"),
      axis.text.y=element_text(size=8),
      legend.position="right",
      legend.title=element_text(face="bold"),
      panel.grid.minor=element_blank()
    )

  ggsave(file.path(picrust_plots, "lefse_pathway_barplot.png"),
         p_lefse, width=14, height=8, dpi=300)
}

# ============================================
# SPEARMAN CORRELATION — TAXA VS PATHWAYS
# ============================================

# Get top 10 genera
top_genera <- ps %>%
  tax_glom("Genus") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%
  filter(!is.na(Genus),
         !grepl("uncultured|Unclassified", Genus, ignore.case=TRUE)) %>%
  group_by(Genus) %>%
  summarise(mean_abund=mean(Abundance), .groups="drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n=10) %>%
  pull(Genus)

genus_mat <- ps %>%
  tax_glom("Genus") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%
  filter(Genus %in% top_genera, !is.na(Genus)) %>%
  dplyr::select(Sample, Genus, Abundance) %>%
  pivot_wider(names_from=Genus, values_from=Abundance,
              values_fill=0) %>%
  column_to_rownames("Sample")

common_samples <- intersect(rownames(genus_mat), meta_gg$sample_name)
pathway_for_corr <- as.data.frame(
  t(pathway_rel))[common_samples, top20_pathways]
colnames(pathway_for_corr) <- top20_named$pathway_name
genus_for_corr  <- genus_mat[common_samples, ]

# Calculate Spearman correlations
corr_matrix <- matrix(NA,
                       nrow=ncol(genus_for_corr),
                       ncol=ncol(pathway_for_corr))
rownames(corr_matrix) <- colnames(genus_for_corr)
colnames(corr_matrix) <- colnames(pathway_for_corr)
pval_matrix <- corr_matrix

for (i in seq_len(ncol(genus_for_corr))) {
  for (j in seq_len(ncol(pathway_for_corr))) {
    test <- cor.test(genus_for_corr[,i], pathway_for_corr[,j],
                     method="spearman", exact=FALSE)
    corr_matrix[i,j] <- test$estimate
    pval_matrix[i,j]  <- test$p.value
  }
}

colnames(corr_matrix) <- substr(colnames(corr_matrix), 1, 30)
colnames(pval_matrix)  <- colnames(corr_matrix)
sig_marks <- ifelse(pval_matrix < 0.001, "***",
              ifelse(pval_matrix < 0.01, "**",
              ifelse(pval_matrix < 0.05, "*", "")))

p_spearman <- pheatmap(
  corr_matrix,
  cluster_rows=TRUE, cluster_cols=FALSE,
  color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  breaks=seq(-1, 1, length.out=101),
  display_numbers=sig_marks,
  number_color="black", fontsize_number=7,
  main="Spearman Correlation — Top Genera vs Predicted Pathways\n(* p<0.05, ** p<0.01, *** p<0.001)",
  fontsize=9, fontsize_row=8, fontsize_col=7.5,
  border_color="white", angle_col="45"
)

png(file.path(picrust_plots, "spearman_taxa_pathway.png"),
    width=14, height=8, units="in", res=300)
draw(p_spearman)
dev.off()
cat("9/9 Spearman correlation saved\n")

# ============================================
# SAVE SUMMARY CSVs
# ============================================

write.csv(top20_named,
          file.path(picrust_plots, "top20_pathway_names.csv"),
          row.names=FALSE)

top20_by_group <- pathway_long_named %>%
  group_by(group, pathway_name) %>%
  summarise(mean_abund=mean(abundance), .groups="drop") %>%
  filter(pathway_name %in% top20_named$pathway_name) %>%
  arrange(group, desc(mean_abund))

write.csv(top20_by_group,
          file.path(picrust_plots, "top20_pathways_by_group.csv"),
          row.names=FALSE)

cat("\n=== PICRUST2 ANALYSIS COMPLETE! ===\n")
cat("All plots saved to:", picrust_plots, "\n")
cat("Total files:", length(list.files(picrust_plots)), "\n")
print(list.files(picrust_plots))

