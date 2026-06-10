# =========================================================================== #
# FigureABC: UMAP of cluster functional annotations, elbow plot,                 #
#         clustree, and marker gene heatmap                                   #
# =========================================================================== #

# ====================== Paths ======================
seurat_rds        <- "path/to/01merge.rds"
marker_txt        <- "path/to/04marker_cluster_m01_l00_r05.txt"
output_dir        <- "path/to/F2"
output_umap       <- file.path(output_dir, "A.pdf")
output_elbow      <- file.path(output_dir, "A_elbow_plot.pdf")
output_clustree   <- file.path(output_dir, "A_clustree.pdf")
output_heatmap    <- file.path(output_dir, "B_marker_genes_heatmap.pdf")

# ====================== Packages ======================
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggrepel)
library(clustree)
library(viridis)
library(grid)

# ====================== Load data ======================
scRNA_merge <- readRDS(seurat_rds)

# ======================================================
# Part A1: UMAP with functional annotation labels
# ======================================================

# ---------- Create functional annotation table ----------
cluster_data <- data.frame(
  Cluster = 0:49,
  product = c(
    "Non-Specialized / Quiet Community",
    "Core Metabolism & Redox Homeostasis",
    "Mixed-Acid Fermenters & Biosynthetic Hubs",
    "Low-Activity Primary Carbohydrate Consumers",
    "Integrated Stress Response & Cellular Maintenance",
    "Gluconeogenic Specialists",
    "Anaerobic Amino Acid Fermenters (Deaminating)",
    "Canonical Heat Shock Response",
    "Lysine-Derived Butyrate & Glutarate Producers",
    "Broad-Spectrum Sugar Transporters & Utilizers",
    "Starch Degradation Specialists",
    "N-Acetylglucosamine (GlcNAc) Utilizers",
    "β-Lactam Resistant Stress Survivors",
    "Glycogen-Storing Fermenters",
    "Propionate Producers (Methylmalonyl-CoA Pathway)",
    "Anabolic NADPH & Redox Powerhouse",
    "Motile Butyrogenic Precursors",
    "Aspartate-Coupled Gluconeogenesis",
    "Stickland Reaction Fermenters (Glycine Reducers)",
    "Cold Shock & Metabolic Adaptation",
    "Sulfidogenic H₂S Producers (Taurine/Choline)",
    "Generalist Hyper-Metabolizers",
    "Oligopeptide Scavengers",
    "Tryptophan & Indole Metabolism Specialists",
    "Protein Quality Control & Refolding",
    "Histidine Degradation & One-Carbon Metabolism",
    "Glycogen-Mobilizing Anaerobic Respirers",
    "Actively Dividing & Heat-Stressed Population",
    "Dissimilatory Sulfate/Thiosulfate Reducers",
    "Specialized Peroxide Scavengers",
    "Osmo-Protective Amino Acid Fermenters",
    "Glycerol Catabolism Specialists",
    "Motile Symbiont Butyrate Producers",
    "Lactose-Specific PTS Utilizers",
    "Amino Acid-Dependent pH Regulators",
    "High-Flux Propionate & Energy Producers",
    "Sporulating & Dormant Cells",
    "Alternative Pathway Propionate Producers",
    "Nitrogenous Nutrient Storage & Reservoirs",
    "Core Protein Stress Responsome",
    "Fatty Acid β-Oxidation Specialists",
    "Pro-Atherogenic TMA Producers (Choline Utilizers)",
    "Lactate-Crossfeeding & Energy Optimizers",
    "Dual-Substrate Choline & Organosulfur Metabolizers",
    "Acid-Resistant GABAergic Neuronal Modulators",
    "Stickland Reaction Fermenters (Proline Oxidizers)",
    "Alternative Electron Acceptor Sulfur Reducers",
    "Motile Mucosal Glycerol Consumers",
    "Cobalt/Nickel Cofactor Scavengers",
    "Polyamine Biosynthesis & Transporters"
  ),
  Rationale = c(
    "Indicates a lack of defining functional markers, suggesting a generalized or inactive state.",
    "Highlights the central carbon metabolic pathways (PPP, Glycolysis) combined with the critical role of reactive oxygen species (ROS) detoxification.",
    "Emphasizes the specific fermentation type and the concurrent high activity of protein synthesis machinery.",
    "Distinguishes this group from Cluster 2 by its lower metabolic activity while maintaining a similar core function.",
    "A comprehensive name covering the response to protein, DNA, and membrane stress.",
    "Directly names the key anabolic pathway for sugar synthesis from non-carbohydrate precursors.",
    "Specifies the energy-generating process and the primary mechanism (deamination) linked to ammonia production.",
    "Uses the standard term for the conserved protein-folding stress response mediated by major chaperones.",
    "Precisely defines the substrate (Lysine) and the two key metabolic outputs.",
    "Describes the capability to recognize and import a wide range of simple sugars.",
    "Highlights the specialization for breaking down a complex polymer, a key step in dietary fiber metabolism.",
    "Uses the specific biochemical name for the amino sugar, indicating a niche role in utilizing host-derived or environmental glycans.",
    "Identifies the specific antibiotic resistance mechanism and the general stress-tolerant phenotype.",
    "Combines the two key functions: energy storage (glycogen synthesis) and energy generation (fermentation).",
    "Names the product and specifies the primary bacterial biosynthetic pathway.",
    "Focuses on the outcome of a robust PPP: generating reducing power (NADPH) for biosynthesis and antioxidant defense.",
    "Indicates motility and the potential for butyrate production, though the pathway may not be complete.",
    "Describes the unique link between the degradation of a specific amino acid and the initiation of sugar synthesis.",
    "Uses the specific biochemical term for this type of amino acid pair fermentation, specifying glycine as the electron acceptor.",
    "Names the primary stressor and indicates a broader metabolic adjustment beyond just the shock response.",
    "Uses the term 'sulfidogenic' to describe H₂S production and specifies the primary substrates.",
    "Conveys a broad metabolic capacity operating at a very high rate across multiple core processes.",
    "Evokes the role of efficiently foraging for and importing short-chain peptides.",
    "Clearly identifies the substrate and the production of a critical host-active metabolite.",
    "A more precise name for the disaggregation and refolding function of ClpB et al.",
    "Highlights the specific amino acid catabolism and its connection to the vital folate cycle.",
    "Emphasizes the breakdown of stored glycogen coupled with anaerobic energy generation via the POR system.",
    "Describes the physiological state (division) and the concomitant stress response often seen in fast-growing cells.",
    "Uses the correct term for anaerobic respiration using sulfur compounds as terminal electron acceptors.",
    "Highlights the niche role in detoxifying hydrogen peroxide, a key host-derived antimicrobial.",
    "Links the fermentation of glycine/betaine to its role in maintaining osmotic balance.",
    "A direct name for the specific function of importing and initiating the metabolism of glycerol.",
    "Positions this cluster as a key mutualist, combining motility with a complete and important beneficial metabolic pathway.",
    "Specifies the sugar and the dedicated phosphotransferase system used for its uptake and phosphorylation.",
    "Focuses on the functional outcome of the transamination and ion transport activities.",
    "Conveys the exceptionally high metabolic throughput and the primary energetic and metabolic outputs.",
    "Describes the developmental state and the resulting metabolically inactive, resistant form.",
    "Suggests a distinct pathway from Cluster 14/35, potentially using different precursors.",
    "Emphasizes the role in acquiring and storing nitrogen (BCAAs, glutathione) and carbon (glycogen).",
    "'Responsome' indicates the complete set of components required for a maximal stress response to proteotoxic damage.",
    "Uses the standard biochemical term for the fatty acid catabolism pathway.",
    "Connects the bacterial metabolism of choline to the host disease pathway (atherogenesis) via TMAO production.",
    "Highlights the role in consuming a microbial metabolite (lactate) and optimizing energy yield via respiration.",
    "Indicates the ability to utilize two distinct classes of substrates (choline and taurine-like compounds).",
    "Combines the survival mechanism (acid resistance) with the production of a neuroactive molecule (GABA), implicating a direct role in the gut-brain axis.",
    "Specifies proline as the electron donor in the Stickland reaction.",
    "Suggests a distinct sulfidogenic pathway from Cluster 28, potentially using different electron acceptors.",
    "Implies a niche within the mucus layer, combining motility with glycerol utilization.",
    "Describes the specific micronutrients targeted and their ultimate purpose as enzyme cofactors.",
    "A direct name for the synthesis and export of polyamines, which are crucial for host cellular functions."
  ),
  stringsAsFactors = FALSE
)

cluster_data$lable_PFG <- paste(cluster_data$Cluster, cluster_data$product, sep = ": ")

# ---------- Merge functional annotations into Seurat object ----------
scRNA_merge@meta.data <- left_join(
  scRNA_merge@meta.data, cluster_data,
  by = c("seurat_clusters" = "Cluster"),
  relationship = "many-to-many"
)
rownames(scRNA_merge@meta.data) <- scRNA_merge@meta.data$barcode

# ---------- Extract UMAP coordinates and metadata ----------
umap_data <- as.data.frame(scRNA_merge@reductions$umap@cell.embeddings)
colnames(umap_data) <- c("UMAP_1", "UMAP_2")

metadata <- scRNA_merge@meta.data
metadata$seurat_clusters <- as.factor(metadata$seurat_clusters)
metadata$lable_PFG <- as.character(metadata$lable_PFG)

# ---------- Compute cluster centers (median) ----------
cluster_centers <- metadata %>%
  bind_cols(umap_data) %>%
  group_by(Cluster = seurat_clusters) %>%
  summarise(
    UMAP1 = median(UMAP_1, na.rm = TRUE),
    UMAP2 = median(UMAP_2, na.rm = TRUE),
    .groups = 'drop'
  )

# ---------- Prepare label mapping in custom order ----------
desired_order <- c(1, 2, 3, 5, 13, 15, 17, 21, 26, 35,
                   4, 7, 19, 24, 29, 39, 9, 10, 11, 22,
                   31, 33, 6, 18, 25, 30, 34, 38, 45, 20,
                   28, 40, 41, 43, 46, 8, 14, 23, 37, 42,
                   44, 49, 16, 27, 32, 36, 47, 12, 48, 0)

label_mapping <- metadata %>%
  distinct(lable_PFG, seurat_clusters) %>%
  group_by(lable_PFG) %>%
  summarise(clusters = paste(sort(unique(seurat_clusters)), collapse = ","),
            .groups = 'drop') %>%
  ungroup() %>%
  mutate(clusters_num = as.numeric(str_extract(clusters, "^\\d+"))) %>%
  mutate(clusters_num = factor(clusters_num, levels = desired_order)) %>%
  arrange(clusters_num) %>%
  mutate(display_label = str_replace(lable_PFG, "^[^_]+_", "")) %>%
  select(-clusters_num)

metadata$lable_PFG <- factor(metadata$lable_PFG,
                             levels = unique(label_mapping$lable_PFG))

# ---------- Generate colors ----------
color_list <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF",
                "#FFA500", "#800080", "#008000", "#000080", "#800000", "#008080",
                "#FF6347", "#40E0D0", "#EE82EE", "#F5DEB3", "#000000")
num_groups <- n_distinct(metadata$lable_PFG)
set.seed(1234)
colors <- sample(colorRampPalette(color_list)(num_groups))

# ---------- Axes limits ----------
x_range <- range(umap_data$UMAP_1)
y_range <- range(umap_data$UMAP_2)
max_span <- max(diff(x_range), diff(y_range)) / 2
x_mid <- mean(x_range)
y_mid <- mean(y_range)

# ---------- UMAP raster plot ----------
base_plot <- ggplot() +
  ggrastr::geom_point_rast(
    data = cbind(umap_data, metadata),
    aes(x = UMAP_1, y = UMAP_2, color = lable_PFG),
    size = 0.8, alpha = 0.7, shape = 16,
    raster.dpi = 300, dev = "cairo"
  ) +
  scale_color_manual(values = colors) +
  coord_fixed(
    xlim = c(x_mid - max_span, x_mid + max_span),
    ylim = c(y_mid - max_span, y_mid + max_span)
  ) +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 7, alpha = 1, shape = 16)
  )) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 22, margin = margin(t = 4, b = 4)),
    legend.title = element_blank(),
    legend.key = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  )

final_plot <- base_plot +
  ggrepel::geom_label_repel(
    data = cluster_centers,
    aes(x = UMAP1, y = UMAP2, label = Cluster),
    size = 5, fontface = "bold", family = "Arial",
    color = "black", fill = alpha("white", 0.8),
    box.padding = 0.5, segment.color = "grey50",
    segment.size = 0.3, min.segment.length = 0.2,
    max.overlaps = Inf, seed = 123
  )

ggsave(output_umap, plot = final_plot, width = 26, height = 12,
       bg = "transparent", limitsize = FALSE, device = cairo_pdf)

# ======================================================
# Part A2: Elbow plot and clustering tree
# ======================================================

# ---------- Elbow plot ----------
p_elbow <- ElbowPlot(scRNA_merge, ndims = 100) +
  geom_vline(xintercept = 40, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(title = "PCA Elbow Plot (dims 1:40 used)") +
  theme_minimal()

ggsave(output_elbow, plot = p_elbow, device = "pdf", width = 8, height = 5)

# ---------- Clustering tree ----------
# Backup original clusters
original_cluster <- scRNA_merge$seurat_clusters

scRNA_merge <- FindClusters(
  scRNA_merge,
  resolution = seq(0.2, 1.6, by = 0.2),
  verbose = FALSE
)

p_tree <- clustree(scRNA_merge, prefix = "SCT_snn.res.") +
  ggtitle("Clustering tree across resolutions") +
  theme(legend.position = "bottom")

ggsave(output_clustree, plot = p_tree, device = "pdf", width = 12, height = 10)

# Restore original clusters
scRNA_merge$seurat_clusters <- original_cluster

# ======================================================
# Part B: Marker gene heatmap per cluster
# ======================================================

# ---------- Load marker table ----------
marker <- read.table(marker_txt, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
marker$cluster <- as.numeric(as.character(marker$cluster))
marker <- marker %>% arrange(cluster)

# ---------- Select top 10 unique genes per cluster ----------
top_markers <- data.frame()
for (clust in 0:49) {
  clust_genes <- marker %>%
    filter(cluster == clust) %>%
    slice_head(n = 10)
  if (nrow(clust_genes) > 0) {
    top_markers <- rbind(top_markers, clust_genes)
  }
}

# Remove duplicate genes across clusters
gene_counts <- table(top_markers$gene)
duplicate_genes <- names(gene_counts[gene_counts > 1])
cluster_gene_count <- rep(0, 50)
names(cluster_gene_count) <- as.character(0:49)

final_markers <- data.frame()
seen_genes <- c()

# Non-duplicated genes first
non_duplicate_genes <- top_markers %>% filter(!gene %in% duplicate_genes)
for (i in 1:nrow(non_duplicate_genes)) {
  gene <- non_duplicate_genes$gene[i]
  cluster <- as.character(non_duplicate_genes$cluster[i])
  if (!gene %in% seen_genes) {
    final_markers <- rbind(final_markers, non_duplicate_genes[i, ])
    seen_genes <- c(seen_genes, gene)
    cluster_gene_count[cluster] <- cluster_gene_count[cluster] + 1
  }
}

# Handle duplicated genes by assigning to cluster with fewest genes
for (gene_name in duplicate_genes) {
  gene_occurrences <- top_markers %>% filter(gene == gene_name)
  cluster_counts <- sapply(gene_occurrences$cluster, function(clust) {
    cluster_gene_count[as.character(clust)]
  })
  min_index <- which.min(cluster_counts)
  selected_cluster <- gene_occurrences$cluster[min_index]
  selected_row <- gene_occurrences[min_index, ]
  if (!gene_name %in% seen_genes) {
    final_markers <- rbind(final_markers, selected_row)
    seen_genes <- c(seen_genes, gene_name)
    cluster_gene_count[as.character(selected_cluster)] <- cluster_gene_count[as.character(selected_cluster)] + 1
  }
}

final_markers <- final_markers %>% arrange(cluster)
unique_genes <- unique(final_markers$gene)

# ---------- Average expression per cluster ----------
expr_matrix <- GetAssayData(scRNA_merge, slot = "data")
genes_in_matrix <- intersect(unique_genes, rownames(expr_matrix))
final_markers <- final_markers %>% filter(gene %in% genes_in_matrix)
unique_genes <- unique(final_markers$gene)

cell_clusters <- Idents(scRNA_merge)
avg_expr <- data.frame()
for (clust in 0:49) {
  clust_cells <- names(cell_clusters[cell_clusters == clust])
  if (length(clust_cells) > 0) {
    clust_expr <- rowMeans(expr_matrix[unique_genes, clust_cells, drop = FALSE])
    avg_expr <- rbind(avg_expr, data.frame(
      gene = names(clust_expr),
      cluster = as.character(clust),
      avg_expression = as.numeric(clust_expr)
    ))
  } else {
    avg_expr <- rbind(avg_expr, data.frame(
      gene = unique_genes,
      cluster = as.character(clust),
      avg_expression = 0
    ))
  }
}

heatmap_data <- avg_expr %>%
  pivot_wider(names_from = cluster, values_from = avg_expression, values_fill = 0) %>%
  column_to_rownames("gene")

# Order genes as per final_markers
final_markers <- final_markers[!duplicated(final_markers$gene), ]
gene_order <- rev(final_markers$gene)
heatmap_data <- heatmap_data[gene_order, ]

# Row-wise z-score
heatmap_data_scaled <- t(scale(t(heatmap_data)))
heatmap_data_scaled[is.na(heatmap_data_scaled)] <- 0

plot_data <- heatmap_data_scaled %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "cluster", values_to = "expression")

gene_cluster_map <- setNames(final_markers$cluster, final_markers$gene)
plot_data$gene_cluster <- gene_cluster_map[plot_data$gene]

# ---------- Custom x-axis order ----------
x_axis_order <- c(0, 1, 2, 3, 5, 13, 15, 17, 21, 26, 35,
                  4, 7, 19, 24, 29, 39, 9, 10, 11, 22,
                  31, 33, 6, 18, 25, 30, 34, 38, 45, 20,
                  28, 40, 41, 43, 46, 8, 14, 23, 37, 42,
                  44, 49, 16, 27, 32, 36, 47, 12, 48)

all_clusters <- 0:49
missing_clusters <- setdiff(all_clusters, x_axis_order)
if (length(missing_clusters) > 0) {
  x_axis_order <- c(x_axis_order, missing_clusters)
}
x_axis_order <- as.character(x_axis_order)

final_markers$cluster <- factor(final_markers$cluster, levels = x_axis_order)
final_markers <- final_markers %>% arrange(cluster)
gene_order <- rev(final_markers$gene)

# Re-align heatmap rows
heatmap_data <- heatmap_data[gene_order, ]

plot_data <- heatmap_data_scaled %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "cluster", values_to = "expression")

plot_data$cluster <- factor(plot_data$cluster, levels = x_axis_order)
plot_data$gene_cluster <- gene_cluster_map[plot_data$gene]
plot_data$gene <- factor(plot_data$gene, levels = gene_order)

# ---------- Y-axis labels (top 2 genes per cluster) ----------
cluster_labels <- sapply(0:49, function(clust) {
  genes <- cluster_gene_lists[[as.character(clust)]]
  if (is.null(genes) || length(genes) == 0) {
    return(as.character(clust))
  } else {
    return(paste(head(genes, 2), collapse = ", "))
  }
})

# Positions for y-axis breaks
y_label_positions <- c()
y_labels <- c()
for (clust in x_axis_order) {
  clust_genes <- names(gene_cluster_map[gene_cluster_map == clust])
  if (length(clust_genes) > 0) {
    mid_index <- ifelse(length(clust_genes) %% 2 == 1,
                        (length(clust_genes) + 1) / 2,
                        length(clust_genes) / 2 + 1)
    y_label_positions <- c(y_label_positions, clust_genes[mid_index])
    y_labels <- c(y_labels, cluster_labels[as.numeric(clust) + 1])
  }
}

# ---------- Heatmap plot ----------
p_heatmap <- ggplot(plot_data, aes(x = cluster, y = gene, fill = expression)) +
  geom_tile(height = 1, width = 1, color = NA) +
  scale_fill_gradientn(
    colors = c("#4C9AC9","#C8E0EF","#E7F2F6","white","#FEE8DD","#F5B99E","#E2745E"),
    values = scales::rescale(c(-1, -0.67, -0.33, 0, 0.33, 0.67, 1)),
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Scaled expression"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(
    expand = c(0, 0),
    breaks = y_label_positions,
    labels = y_labels
  ) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 9, face = "italic", color = "black"),
    axis.title = element_text(size = 12),
    legend.title = element_text(family = "Arial", size = 10),
    legend.text = element_text(family = "Arial", size = 8),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

# ---------- Save heatmap ----------
n_genes <- length(unique(plot_data$gene))
pdf_height <- max(6, n_genes * 0.04)
pdf_width <- 12

cairo_pdf(output_heatmap, width = pdf_width, height = pdf_height, bg = "transparent")
dev.off()

# =========================================================================== #
# Figure DE: Phylum UMAP and stacked bar plots of Genus / Phylum composition   #
#           across clusters; Species / Phylum composition                     #
# =========================================================================== #

# ====================== Paths ======================
seurat_rds       <- "path/to/01merge.rds"
output_dir       <- "path/to/F2"
output_phylum_umap_png <- file.path(output_dir, "E_Phylum_umap.png")
output_genus_phylum_pdf <- file.path(output_dir, "E_Combined_Genus_Phylum_cluster.pdf")
output_species_phylum_pdf <- file.path(output_dir, "E_Combined_Species_Phylum_cluster.pdf")

# ====================== Packages ======================
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(patchwork)
library(ggpointdensity)

# ====================== Load data ======================
scRNA_merge <- readRDS(seurat_rds)

# ====================== Harmonize phylum names ======================
scRNA_merge@meta.data <- scRNA_merge@meta.data %>%
  mutate(Phylum = case_when(
    Phylum %in% c("Firmicutes_A", "Firmicutes_B", "Firmicutes_C", "Firmicutes_G") ~ "Firmicutes",
    TRUE ~ Phylum
  ))

# Mark phyla with <50 cells as "Others"
serial_cellnum_list <- scRNA_merge@meta.data %>%
  group_by(Phylum) %>%
  summarise(cellnum = n(), .groups = 'drop') %>%
  arrange(Phylum) %>%
  mutate(serial = as.integer(factor(Phylum)) - 1)

scRNA_merge@meta.data <- left_join(scRNA_merge@meta.data, serial_cellnum_list,
                                   by = "Phylum", relationship = "many-to-many")
scRNA_merge@meta.data$cellnum_Phylum <- paste(scRNA_merge@meta.data$Phylum,
                                              "(", scRNA_merge@meta.data$cellnum, ")", sep = "")
scRNA_merge@meta.data$Phylum[scRNA_merge@meta.data$cellnum < 50] <- "Others"
scRNA_merge@meta.data <- scRNA_merge@meta.data %>%
  select(-cellnum, -serial, -cellnum_Phylum)

# Recalculate final cell counts per phylum after grouping
serial_cellnum_list <- scRNA_merge@meta.data %>%
  group_by(Phylum) %>%
  summarise(cellnum = n(), .groups = 'drop') %>%
  arrange(Phylum) %>%
  mutate(serial = as.integer(factor(Phylum)) - 1)

scRNA_merge@meta.data <- left_join(scRNA_merge@meta.data, serial_cellnum_list,
                                   by = "Phylum", relationship = "many-to-many")
scRNA_merge@meta.data$cellnum_Phylum <- paste(scRNA_merge@meta.data$Phylum,
                                              "(", scRNA_merge@meta.data$cellnum, ")", sep = "")
scRNA_merge@meta.data <- scRNA_merge@meta.data %>% select(-cellnum, -serial)
rownames(scRNA_merge@meta.data) <- scRNA_merge@meta.data$barcode

# ======================================================
# Part 1: Phylum UMAP
# ======================================================

# Generate colors for phyla
phyla_all <- sort(unique(scRNA_merge$Phylum))
num_phyla <- length(phyla_all)
if (num_phyla <= 12) {
  colors_phylum <- brewer.pal(num_phyla, "Set2")
} else {
  colors_phylum <- colorRampPalette(brewer.pal(12, "Set2"))(num_phyla)
}
set.seed(1)
colors_phylum <- sample(colors_phylum)

# UMAP without labels
pic5 <- DimPlot(scRNA_merge, reduction = "umap", group.by = "Phylum",
                repel = TRUE, pt.size = 0.5, raster = FALSE,
                label = FALSE) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(family = "Times New Roman")) +
  theme_dr(xlength = 0.1, ylength = 0.1) +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        legend.box.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "gray40", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank()) +
  scale_color_manual(values = colors_phylum) +
  coord_fixed() +
  guides(x = "none", y = "none")

ggsave(output_phylum_umap_png, plot = pic5, width = 8, height = 8, dpi = 300,
       bg = "white", limitsize = FALSE)

# ======================================================
# Part 2: Stacked bar plots (Genus/Phylum and Species/Phylum)
# ======================================================

# Define custom cluster order for x-axis
x_axis_order <- c(0, 1, 2, 3, 5, 13, 15, 17, 21, 26, 35,
                  4, 7, 19, 24, 29, 39, 9, 10, 11, 22,
                  31, 33, 6, 18, 25, 30, 34, 38, 45, 20,
                  28, 40, 41, 43, 46, 8, 14, 23, 37, 42,
                  44, 49, 16, 27, 32, 36, 47, 12, 48)

sc_df <- scRNA_merge@meta.data

# ---------- Genus vs Phylum ----------
genus_counts <- sc_df %>%
  count(Genus) %>%
  arrange(desc(n))
top_20_genus <- head(genus_counts$Genus, 20)

sc_df <- sc_df %>%
  mutate(Genus_Grouped = if_else(Genus %in% top_20_genus,
                                 as.character(Genus),
                                 "Other Genus"))
genus_levels <- c(setdiff(unique(sc_df$Genus_Grouped), "Other Genus"), "Other Genus")
sc_df$Genus_Grouped <- factor(sc_df$Genus_Grouped, levels = genus_levels)

cluster_genus_data <- sc_df %>%
  count(seurat_clusters, Genus_Grouped) %>%
  group_by(seurat_clusters) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

num_genus <- length(levels(sc_df$Genus_Grouped)) - 1
genus_colors <- colorRampPalette(brewer.pal(12, "Paired"))(num_genus)
genus_colors <- c(genus_colors, "white")
names(genus_colors) <- levels(sc_df$Genus_Grouped)

Genus_plot <- ggplot(cluster_genus_data,
                     aes(x = factor(seurat_clusters, levels = x_axis_order),
                         y = percent, fill = Genus_Grouped)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = percent_format()) +
  scale_fill_manual(values = genus_colors, guide = guide_legend(ncol = 1)) +
  labs(x = NULL, y = "Percentage", fill = "Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  )

# Phylum bar plot (same phylum colors as before)
cluster_phylum_data <- sc_df %>%
  count(seurat_clusters, Phylum) %>%
  group_by(seurat_clusters) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

phylum_totals <- cluster_phylum_data %>%
  group_by(Phylum) %>%
  summarise(total_n = sum(n)) %>%
  arrange(desc(total_n))

phylum_levels <- phylum_totals$Phylum
cluster_phylum_data$Phylum <- factor(cluster_phylum_data$Phylum, levels = phylum_levels)

# Use same color mapping for phyla as in UMAP
color_mapping <- setNames(colors_phylum, phyla_all)

Phylum_plot <- ggplot(cluster_phylum_data,
                      aes(x = factor(seurat_clusters, levels = x_axis_order),
                          y = percent, fill = Phylum)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = percent_format()) +
  scale_fill_manual(values = color_mapping, guide = guide_legend(ncol = 1)) +
  labs(x = "Cluster", y = "Percentage", fill = "Phylum") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(),
    axis.ticks.x = element_line(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

combined_genus_phylum <- Phylum_plot / Genus_plot + plot_layout(heights = c(1, 1))
ggsave(output_genus_phylum_pdf, plot = combined_genus_phylum,
       width = 15, height = 12, bg = "transparent", limitsize = FALSE, device = cairo_pdf)

# ---------- Species vs Phylum ----------
species_counts <- sc_df %>%
  count(Species) %>%
  arrange(desc(n))
top_20_species <- head(species_counts$Species, 20)

sc_df <- sc_df %>%
  mutate(Species_Grouped = if_else(Species %in% top_20_species,
                                   as.character(Species),
                                   "Other Species"))
species_levels <- c(setdiff(unique(sc_df$Species_Grouped), "Other Species"), "Other Species")
sc_df$Species_Grouped <- factor(sc_df$Species_Grouped, levels = species_levels)

cluster_species_data <- sc_df %>%
  count(seurat_clusters, Species_Grouped) %>%
  group_by(seurat_clusters) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

num_species <- length(levels(sc_df$Species_Grouped)) - 1
species_colors <- colorRampPalette(brewer.pal(12, "Paired"))(num_species)
species_colors <- c(species_colors, "white")
names(species_colors) <- levels(sc_df$Species_Grouped)

Species_plot <- ggplot(cluster_species_data,
                       aes(x = factor(seurat_clusters, levels = x_axis_order),
                           y = percent, fill = Species_Grouped)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = percent_format()) +
  scale_fill_manual(values = species_colors, guide = guide_legend(ncol = 1)) +
  labs(x = NULL, y = "Percentage", fill = "Species") +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  )

# Phylum bar (reuse same data but with Times New Roman)
Phylum_plot_species <- ggplot(cluster_phylum_data,
                              aes(x = factor(seurat_clusters, levels = x_axis_order),
                                  y = percent, fill = Phylum)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = percent_format()) +
  scale_fill_manual(values = color_mapping, guide = guide_legend(ncol = 1)) +
  labs(x = "Cluster", y = "Percentage", fill = "Phylum") +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    axis.text.x = element_text(family = "Times New Roman"),
    axis.ticks.x = element_line(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

combined_species_phylum <- Species_plot / Phylum_plot_species + plot_layout(heights = c(1, 1))
ggsave(output_species_phylum_pdf, plot = combined_species_phylum,
       width = 15, height = 12, bg = "transparent", limitsize = FALSE, device = cairo_pdf)


# ============================================================
# Figure FG
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggpubr)       # provides geom_signif
library(extrafont)    # for loading system fonts

# ============================================================
# Set file paths (replace "path/to/" with actual paths)
# ============================================================
metadata_path <- "path/to/metadata.rds"
output_dir    <- "path/to/output/"

# ============================================================
# Read data
# ============================================================
metadata_raw <- readRDS(metadata_path)

# Define the 8 functional categories and their cluster assignments
cluster_categories <- list(
  "Energy & Carbon"          = c(1, 2, 3, 5, 13, 15, 17, 21, 26, 35),
  "Stress & Maintenance"     = c(4, 7, 19, 24, 29, 39),
  "Substrate & Fermentation" = c(6, 9, 10, 11, 18, 20, 22, 25, 28, 30, 31, 33, 34, 38, 40, 41, 43, 45, 46),
  "Beneficial Metabolite"    = c(8, 14, 23, 37, 42, 44, 49),
  "Cell Growth"              = c(16, 27, 32, 36, 47),
  "Defense & Resistance"     = c(12),
  "Cofactor & Micronutrient" = c(48),
  "Unspecified"              = c(0)
)

# Order of functional categories for plots
cluster_order <- c(
  "Cofactor & Micronutrient",
  "Cell Growth",
  "Beneficial Metabolite",
  "Substrate & Fermentation",
  "Unspecified",
  "Energy & Carbon",
  "Stress & Maintenance",
  "Defense & Resistance"
)

# Helper function to map seurat_clusters to categories
map_to_category <- function(cluster_vec) {
  res <- rep(NA_character_, length(cluster_vec))
  for (cat in names(cluster_categories)) {
    res[cluster_vec %in% cluster_categories[[cat]]] <- cat
  }
  # Any remaining NA should be "Unspecified" (catches 0 and others)
  res[is.na(res)] <- "Unspecified"
  return(res)
}

# ============================================================
# Part 1: Phylum-level stacked area chart
# ============================================================
metadata_phylum <- metadata_raw

# Rename phyla with < 50 cells to "Others"
phyla_counts <- metadata_phylum %>%
  group_by(Phylum) %>%
  summarise(cellnum = n(), .groups = 'drop')

low_abundant <- phyla_counts$Phylum[phyla_counts$cellnum < 50]
metadata_phylum$Phylum[metadata_phylum$Phylum %in% low_abundant] <- "Others"
metadata_phylum$Phylum <- factor(metadata_phylum$Phylum)

# Add functional category
metadata_phylum$cluster_category <- factor(
  map_to_category(as.numeric(as.character(metadata_phylum$seurat_clusters))),
  levels = cluster_order
)

# Build complete grid and compute percentages
all_combinations_phylum <- expand.grid(
  cluster_category = levels(metadata_phylum$cluster_category),
  Phylum           = levels(metadata_phylum$Phylum)
)

plot_data1 <- metadata_phylum %>%
  group_by(cluster_category, Phylum) %>%
  summarise(count = n(), .groups = 'drop') %>%
  right_join(all_combinations_phylum, by = c("cluster_category", "Phylum")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(cluster_category) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Color mapping for phyla
phyla_ordered <- levels(metadata_phylum$Phylum)
num_phyla <- length(phyla_ordered)
if (num_phyla <= 12) {
  phylum_colors <- brewer.pal(num_phyla, "Set2")
} else {
  set2_colors <- brewer.pal(12, "Set2")
  phylum_colors <- colorRampPalette(set2_colors)(num_phyla)
}
set.seed(1)
phylum_colors <- sample(phylum_colors)
color_mapping1 <- setNames(phylum_colors, phyla_ordered)

# Prepare x-axis numeric positions
plot_data1 <- plot_data1 %>%
  mutate(
    cluster_category = factor(cluster_category, levels = cluster_order),
    category_num = as.numeric(cluster_category)
  )
x_breaks1 <- unique(plot_data1$category_num)

p1 <- ggplot(plot_data1, aes(x = category_num, y = percentage, fill = Phylum)) +
  geom_area(position = "stack", alpha = 0.3) +
  geom_col(position = "stack", width = 0.5, alpha = 1) +
  labs(x = "", y = "% among all cells") +
  theme_minimal() +
  theme(
    panel.grid        = element_blank(),
    axis.line.x       = element_blank(),
    axis.line.y       = element_line(color = "black"),
    axis.ticks.y      = element_line(color = "black"),
    plot.margin       = margin(5, 5, 5, 5),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    text              = element_text(family = "Times New Roman"),
    axis.text         = element_text(family = "Times New Roman"),
    axis.title        = element_text(family = "Times New Roman"),
    legend.text       = element_text(family = "Times New Roman"),
    legend.title      = element_text(family = "Times New Roman"),
    plot.title        = element_text(hjust = 0.5, size = 12)
  ) +
  scale_x_continuous(breaks = x_breaks1, labels = NULL) +
  scale_fill_manual(values = color_mapping1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off")

# ============================================================
# Part 2: Genus-level stacked area chart
# ============================================================
metadata_genus <- metadata_raw

# Map functional categories
metadata_genus$cluster_category <- factor(
  map_to_category(as.numeric(as.character(metadata_genus$seurat_clusters))),
  levels = cluster_order
)

# Identify top 20 genera, group the rest as "Other Genus"
genus_counts <- metadata_genus %>%
  count(Genus) %>%
  arrange(desc(n))
top20_genera <- head(genus_counts$Genus, 20)
metadata_genus$Species_Grouped <- if_else(
  metadata_genus$Genus %in% top20_genera,
  as.character(metadata_genus$Genus),
  "Other Genus"
)
species_levels <- c(setdiff(unique(metadata_genus$Species_Grouped), "Other Genus"), "Other Genus")
metadata_genus$Species_Grouped <- factor(metadata_genus$Species_Grouped, levels = species_levels)

# Complete grid and percentages
all_combinations_genus <- expand.grid(
  cluster_category = levels(metadata_genus$cluster_category),
  Species_Grouped  = levels(metadata_genus$Species_Grouped)
)

plot_data2 <- metadata_genus %>%
  group_by(cluster_category, Species_Grouped) %>%
  summarise(count = n(), .groups = 'drop') %>%
  right_join(all_combinations_genus, by = c("cluster_category", "Species_Grouped")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(cluster_category) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Color mapping for genera
num_species <- length(species_levels) - 1  # exclude "Other Genus"
species_colors <- colorRampPalette(brewer.pal(12, "Paired"))(num_species)
species_colors <- c(species_colors, "#F0F0F0")  # light gray for Other Genus
color_mapping2 <- setNames(species_colors, species_levels)

# Reorder so that "Other Genus" is plotted first (at bottom of stack)
plot_data2 <- plot_data2 %>%
  mutate(
    cluster_category = factor(cluster_category, levels = cluster_order),
    category_num     = as.numeric(cluster_category),
    Species_Grouped_plot = factor(Species_Grouped, levels = rev(species_levels))
  )

x_breaks2 <- unique(plot_data2$category_num)

p2 <- ggplot(plot_data2, aes(x = category_num, y = percentage, fill = Species_Grouped_plot)) +
  geom_area(position = "stack", alpha = 0.3) +
  geom_col(position = "stack", width = 0.5, alpha = 1) +
  labs(x = "", y = "% among all cells") +
  theme_minimal() +
  theme(
    panel.grid        = element_blank(),
    axis.line.x       = element_blank(),
    axis.line.y       = element_line(color = "black"),
    axis.ticks.y      = element_line(color = "black"),
    plot.margin       = margin(5, 5, 60, 5),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1, family = "Times New Roman"),
    text              = element_text(family = "Times New Roman"),
    axis.text         = element_text(family = "Times New Roman"),
    axis.title        = element_text(family = "Times New Roman"),
    legend.text       = element_text(family = "Times New Roman"),
    legend.title      = element_text(family = "Times New Roman"),
    plot.title        = element_text(hjust = 0.5, size = 12)
  ) +
  scale_x_continuous(breaks = x_breaks2, labels = cluster_order) +
  scale_fill_manual(
    values = color_mapping2,
    name   = "Genus",
    guide  = guide_legend(ncol = 1),
    breaks = species_levels,
    labels = species_levels
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off")

# ============================================================
# Combine Phylum and Genus charts (G_combined)
# ============================================================
combined_plot <- p1 / p2 + plot_layout(heights = c(1, 1))

ggsave(
  filename = paste0(output_dir, "F_combined_Genus.pdf"),
  plot     = combined_plot,
  width    = 8,
  height   = 8,
  dpi      = 500,
  bg       = "white",
  limitsize = FALSE,
  device   = cairo_pdf
)

# ============================================================
# Part 3: Boxplots for Bacteroidota, Firmicutes, Proteobacteria (G)
# ============================================================
# Use metadata_genus (contains original Phylum and cluster_category)
metadata_box <- metadata_genus

# Manual colors for the 8 functional categories
manual_colors_8 <- c(
  "#3FA9F5",  # Energy & Carbon
  "#BDCCD4",  # Stress & Maintenance
  "#FF7BAC",  # Substrate & Fermentation
  "#7AC943",  # Beneficial Metabolite
  "#FF931E",  # Cell Growth
  "#00FFFF",  # Defense & Resistance
  "#FCEE21",  # Cofactor & Micronutrient
  "#C69C6D"   # Unspecified
)
category_color_mapping <- setNames(manual_colors_8, cluster_order)

# Function to create boxplot for a single phylum
create_phylum_boxplot <- function(phylum_name, data, cat_order, color_map) {
  # Compute abundance per cluster
  cluster_data <- data %>%
    group_by(seurat_clusters) %>%
    summarise(
      total_cells  = n(),
      phylum_cells = sum(Phylum == phylum_name, na.rm = TRUE),
      .groups      = 'drop'
    ) %>%
    mutate(abundance = (phylum_cells / total_cells) * 100) %>%
    left_join(
      data %>% select(seurat_clusters, cluster_category) %>% distinct(),
      by = "seurat_clusters"
    ) %>%
    mutate(cluster_category = factor(cluster_category, levels = cat_order))
  
  # Only perform pairwise t-test if at least 2 categories with >=2 clusters
  cat_counts <- cluster_data %>%
    group_by(cluster_category) %>%
    summarise(n_clusters = n(), .groups = 'drop')
  eligible <- cat_counts$cluster_category[cat_counts$n_clusters >= 2]
  
  sig_data <- NULL
  if (length(eligible) >= 2) {
    eligible_data <- cluster_data %>% filter(cluster_category %in% eligible)
    combos <- combn(as.character(eligible), 2, simplify = FALSE)
    test_results <- lapply(combos, function(groups) {
      g1 <- groups[1]; g2 <- groups[2]
      v1 <- eligible_data$abundance[eligible_data$cluster_category == g1]
      v2 <- eligible_data$abundance[eligible_data$cluster_category == g2]
      if (length(v1) >= 2 && length(v2) >= 2) {
        tt <- t.test(v1, v2, var.equal = FALSE)
        data.frame(group1 = g1, group2 = g2, p.value = tt$p.value, stringsAsFactors = FALSE)
      } else {
        NULL
      }
    })
    all_tests <- bind_rows(test_results)
    if (nrow(all_tests) > 0) {
      sig <- all_tests %>% filter(p.value < 0.05)
      if (nrow(sig) > 0) {
        sig_data <- sig %>%
          mutate(
            p.signif = case_when(
              p.value < 0.001 ~ "***",
              p.value < 0.01  ~ "**",
              p.value < 0.05  ~ "*",
              TRUE            ~ "ns"
            )
          ) %>%
          filter(p.signif != "ns")
      }
    }
  }
  
  max_y <- max(cluster_data$abundance, na.rm = TRUE)
  
  p <- ggplot(cluster_data, aes(x = cluster_category, y = abundance, fill = cluster_category)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.7, color = "black") +
    labs(x = "", y = "Abundance (%)", title = phylum_name) +
    theme_classic() +
    theme(
      text          = element_text(family = "Times New Roman"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 12),
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y  = element_text(size = 12),
      legend.position = "none"
    ) +
    scale_fill_manual(values = color_map, drop = FALSE)
  
  if (!is.null(sig_data) && nrow(sig_data) > 0) {
    for (i in seq_len(nrow(sig_data))) {
      p <- p + geom_signif(
        comparisons = list(c(sig_data$group1[i], sig_data$group2[i])),
        annotations = sig_data$p.signif[i],
        y_position  = max_y * (1 + 0.05 * i),
        tip_length  = 0.01,
        vjust       = 0.5,
        textsize    = 4,
        family      = "Times New Roman"
      )
    }
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15 + 0.05 * nrow(sig_data))))
  }
  return(p)
}

# Generate plots for the three phyla
phyla_list <- c("Bacteroidota", "Firmicutes", "Proteobacteria")
plot_list_G <- lapply(phyla_list, function(ph) {
  create_phylum_boxplot(ph, metadata_box, cluster_order, category_color_mapping)
})

final_plot_G <- plot_list_G[[1]] | plot_list_G[[2]] | plot_list_G[[3]]
final_plot_G <- final_plot_G +
  plot_layout(ncol = 3, widths = c(1, 1, 1)) &
  theme(plot.margin = margin(10, 10, 10, 10))

ggsave(
  filename = paste0(output_dir, "G.pdf"),
  plot     = final_plot_G,
  width    = 14,
  height   = 6,
  bg       = "transparent",
  limitsize = FALSE,
  device   = cairo_pdf
)

# ============================================================
# Part 4: Boxplots for groups A, B, C (Bacteroidota & Firmicutes) (G_down)
# ============================================================
# Define groups A, B, C (Substrate & Fermentation sub-clusters)
group_ABC <- list(
  "A" = c(9, 10, 11, 22, 31, 33),
  "B" = c(6, 18, 25, 30, 34, 38, 45),
  "C" = c(20, 28, 40, 41, 43, 46)
)

# Filter metadata_genus to only clusters in A, B, C
metadata_abc <- metadata_genus %>%
  filter(seurat_clusters %in% unlist(group_ABC)) %>%
  mutate(
    cluster_category = case_when(
      seurat_clusters %in% group_ABC[["A"]] ~ "A",
      seurat_clusters %in% group_ABC[["B"]] ~ "B",
      seurat_clusters %in% group_ABC[["C"]] ~ "C",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cluster_category)) %>%
  mutate(cluster_category = factor(cluster_category, levels = c("A", "B", "C")))

# Manual colors for A, B, C
abc_colors <- c("#66C1A4", "#8C9FCA", "#FB8C62")
abc_color_map <- setNames(abc_colors, c("A", "B", "C"))

# Legend labels
abc_labels <- c(
  "A" = "A: Complex & Simple Carbohydrate Utilizers",
  "B" = "B: Amino Acid & Nitrogenous Compound Fermenters",
  "C" = "C: Lipid & Organosulfur Compound Metabolizers"
)

# Boxplot function adapted for A, B, C groups
create_abc_boxplot <- function(phylum_name, data, cat_levels, color_map, leg_labels) {
  cluster_data <- data %>%
    group_by(seurat_clusters) %>%
    summarise(
      total_cells  = n(),
      phylum_cells = sum(Phylum == phylum_name, na.rm = TRUE),
      .groups      = 'drop'
    ) %>%
    mutate(abundance = (phylum_cells / total_cells) * 100) %>%
    left_join(
      data %>% select(seurat_clusters, cluster_category) %>% distinct(),
      by = "seurat_clusters"
    ) %>%
    mutate(cluster_category = factor(cluster_category, levels = cat_levels))
  
  # t-test between groups
  cat_counts <- cluster_data %>%
    group_by(cluster_category) %>%
    summarise(n_clusters = n(), .groups = 'drop')
  eligible <- cat_counts$cluster_category[cat_counts$n_clusters >= 2]
  
  sig_data <- NULL
  if (length(eligible) >= 2) {
    eligible_data <- cluster_data %>% filter(cluster_category %in% eligible)
    combos <- combn(as.character(eligible), 2, simplify = FALSE)
    test_results <- lapply(combos, function(groups) {
      g1 <- groups[1]; g2 <- groups[2]
      v1 <- eligible_data$abundance[eligible_data$cluster_category == g1]
      v2 <- eligible_data$abundance[eligible_data$cluster_category == g2]
      if (length(v1) >= 2 && length(v2) >= 2) {
        tt <- t.test(v1, v2, var.equal = FALSE)
        data.frame(group1 = g1, group2 = g2, p.value = tt$p.value, stringsAsFactors = FALSE)
      } else {
        NULL
      }
    })
    all_tests <- bind_rows(test_results)
    if (nrow(all_tests) > 0) {
      sig <- all_tests %>% filter(p.value < 0.05)
      if (nrow(sig) > 0) {
        sig_data <- sig %>%
          mutate(
            p.signif = case_when(
              p.value < 0.001 ~ "***",
              p.value < 0.01  ~ "**",
              p.value < 0.05  ~ "*",
              TRUE            ~ "ns"
            )
          ) %>%
          filter(p.signif != "ns")
      }
    }
  }
  
  max_y <- max(cluster_data$abundance, na.rm = TRUE)
  
  p <- ggplot(cluster_data, aes(x = cluster_category, y = abundance, fill = cluster_category)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.7, color = "black") +
    labs(x = "", y = "Abundance (%)", title = phylum_name, fill = "") +
    theme_classic() +
    theme(
      text           = element_text(family = "Times New Roman"),
      axis.text.x    = element_text(angle = 0, hjust = 1, size = 12),
      plot.title     = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.y   = element_text(size = 12),
      legend.position = "bottom",
      legend.title   = element_text(size = 12, face = "bold"),
      legend.text    = element_text(size = 10)
    ) +
    scale_fill_manual(values = color_map, labels = leg_labels, drop = FALSE)
  
  if (!is.null(sig_data) && nrow(sig_data) > 0) {
    for (i in seq_len(nrow(sig_data))) {
      p <- p + geom_signif(
        comparisons = list(c(sig_data$group1[i], sig_data$group2[i])),
        annotations = sig_data$p.signif[i],
        y_position  = max_y * (1 + 0.05 * i),
        tip_length  = 0.01,
        vjust       = 0.5,
        textsize    = 4,
        family      = "Times New Roman"
      )
    }
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15 + 0.05 * nrow(sig_data))))
  }
  return(p)
}

# Create plots for Bacteroidota and Firmicutes
plot_abc_bac <- create_abc_boxplot("Bacteroidota", metadata_abc, c("A","B","C"), abc_color_map, abc_labels)
plot_abc_firm <- create_abc_boxplot("Firmicutes", metadata_abc, c("A","B","C"), abc_color_map, abc_labels)

final_plot_abc <- plot_abc_bac | plot_abc_firm
final_plot_abc <- final_plot_abc +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  theme(
    plot.margin        = margin(10, 10, 10, 10),
    text               = element_text(family = "Times New Roman"),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 12, face = "bold", margin = margin(b = 5)),
    legend.text        = element_text(size = 10, margin = margin(r = 10)),
    legend.spacing.x   = unit(0.5, "cm"),
    legend.box.margin  = margin(10, 0, 0, 0)
  )

ggsave(
  filename = paste0(output_dir, "G_down.pdf"),
  plot     = final_plot_abc,
  width    = 10,
  height   = 6,
  bg       = "transparent",
  limitsize = FALSE,
  device   = cairo_pdf
)



# ===========================================================================
# Figure H: Circular phylogenetic tree of bacterial species, with bar chart
#         showing normalized Gini coefficient (cell count) colored by
#         dominant cluster proportion
# ===========================================================================

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ineq)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(readr)
library(scales)
library(tibble)

# Set file paths (replace "path/to/" with actual paths)
metadata_path   <- "path/to/metadata.rds"
genome_map_path <- "path/to/01genome_map_v2.1.tsv"
tree_path       <- "path/to/reference_genomes.tre.treefile"
abundance_out   <- "path/to/H_species_cluster_abundance.csv"
figure_out      <- "path/to/H_gini_circular_tree.pdf"

# Read data
metadata_raw <- readRDS(metadata_path)

# ============================================================
# 1. Save raw species-cluster abundance matrix (all species)
# ============================================================
cluster_species_counts_raw <- as.data.frame(
  table(Cluster = metadata_raw$seurat_clusters, Species = metadata_raw$Species)
)
species_wide_raw <- dcast(cluster_species_counts_raw, Cluster ~ Species,
                          value.var = "Freq", fill = 0)
species_matrix_raw <- t(as.matrix(species_wide_raw[, -1]))
colnames(species_matrix_raw) <- species_wide_raw$Cluster
species_abundance <- sweep(species_matrix_raw, 1, rowSums(species_matrix_raw), FUN = "/")
species_abundance[is.na(species_abundance)] <- 0
write.csv(species_abundance, file = abundance_out, quote = FALSE)

# ============================================================
# 2. Gini coefficient computation (species with ≥100 cells)
# ============================================================
meta_filtered <- metadata_raw %>%
  group_by(Species) %>%
  filter(n() >= 100) %>%
  ungroup()

cluster_species_counts <- as.data.frame(
  table(Cluster = meta_filtered$seurat_clusters, Species = meta_filtered$Species)
)
species_wide <- dcast(cluster_species_counts, Cluster ~ Species,
                      value.var = "Freq", fill = 0)
species_matrix <- t(as.matrix(species_wide[, -1]))
colnames(species_matrix) <- species_wide$Cluster
cluster_total_cells <- colSums(species_matrix)

# Gini (cell count)
species_gini_cell <- apply(species_matrix, 1, function(x) {
  x_pos <- x[x > 0]
  if (length(x_pos) <= 1) return(ifelse(sum(x_pos) == 0, NA, 1))
  Gini(x_pos)
})

# Gini (nFeature_SCT, if available)
if ("nFeature_SCT" %in% colnames(meta_filtered)) {
  weighted_data <- meta_filtered %>%
    group_by(seurat_clusters, Species) %>%
    summarise(total_nFeature_SCT = sum(nFeature_SCT, na.rm = TRUE), .groups = 'drop')
  weighted_wide <- weighted_data %>%
    pivot_wider(names_from = Species, values_from = total_nFeature_SCT, values_fill = 0)
  species_matrix_weighted <- t(as.matrix(weighted_wide[, -1]))
  colnames(species_matrix_weighted) <- weighted_wide$seurat_clusters
  species_gini_weighted <- apply(species_matrix_weighted, 1, function(x) {
    x_pos <- x[x > 0]
    if (length(x_pos) <= 1) return(ifelse(sum(x_pos) == 0, NA, 1))
    Gini(x_pos)
  })
} else {
  species_gini_weighted <- rep(NA, length(species_gini_cell))
  names(species_gini_weighted) <- names(species_gini_cell)
}

# Results data frame
results <- data.frame(
  Species = names(species_gini_cell),
  Gini_CellCount = species_gini_cell,
  Gini_nFeature_SCT = species_gini_weighted[match(names(species_gini_cell), names(species_gini_weighted))],
  Total_Cells = rowSums(species_matrix),
  Num_Clusters = apply(species_matrix > 0, 1, sum),
  stringsAsFactors = FALSE
)

# Dominant cluster information
for (i in seq_len(nrow(results))) {
  sp <- results$Species[i]
  cnt <- species_matrix[sp, ]
  tot <- sum(cnt)
  if (tot > 0) {
    max_idx <- which.max(cnt)
    dom_cluster <- colnames(species_matrix)[max_idx]
    dom_prop <- cnt[max_idx] / tot
    abund_in_dom <- cnt[max_idx] / cluster_total_cells[dom_cluster]
    results$Dominant_Cluster[i] <- dom_cluster
    results$Dominant_Proportion[i] <- round(dom_prop, 3)
    results$Abundance_in_Dominant_Cluster[i] <- round(abund_in_dom, 5)
  } else {
    results$Dominant_Cluster[i] <- NA
    results$Dominant_Proportion[i] <- NA
    results$Abundance_in_Dominant_Cluster[i] <- NA
  }
}

# Filter and normalize Gini (0-1 scaling)
results <- results %>% filter(Total_Cells > 1000)
gini_min <- min(results$Gini_CellCount, na.rm = TRUE)
gini_max <- max(results$Gini_CellCount, na.rm = TRUE)
results$Gini_CellCount <- (results$Gini_CellCount - gini_min) / (gini_max - gini_min)

# ============================================================
# 3. Prepare tree annotation data
# ============================================================
# Genome-to-species mapping
df_map <- read_tsv(genome_map_path, col_names = FALSE)
folder_path <- "path/to/reference_genomes"
file_names <- list.files(folder_path, pattern = "\\.fna\\.gz$")
match_strings <- gsub("\\.fna\\.gz$", "", file_names)
df_map <- df_map %>% filter(X1 %in% match_strings)

metadata_known <- metadata_raw %>% filter(Species != "Unknown")

species_info <- metadata_known %>%
  group_by(Phylum, Family, Species) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(abundance = (count / sum(count)) * 100) %>%
  select(-count) %>%
  arrange(desc(abundance)) %>%
  left_join(df_map, by = c("Species" = "X2")) %>%
  rename(Genome = X1)

merged_df <- species_info %>%
  left_join(results, by = "Species") %>%
  select(Species, Genome, Gini_CellCount, Dominant_Proportion) %>%
  filter(!is.na(Gini_CellCount), !is.na(Dominant_Proportion), !is.na(Genome)) %>%
  column_to_rownames("Genome")

# ============================================================
# 4. Circular phylogenetic tree with bar layer
# ============================================================
tree <- read.tree(tree_path)
tree$tip.label <- trimws(tree$tip.label)

# Retain only matched tips
merged_df <- merged_df[rownames(merged_df) %in% tree$tip.label, , drop = FALSE]

# Base tree
p_base <- ggtree(tree, layout = "circular")
tree_data <- p_base$data
tree_data$label <- trimws(tree_data$label)

# Merge annotation data
merged_df <- merged_df %>% mutate(label = rownames(merged_df))
tree_data_merged <- tree_data %>% left_join(merged_df, by = "label")
p_base$data <- tree_data_merged
p_base$data$Species <- sub("^([^ ])[^ ]+", "\\1.", p_base$data$Species)

# Add bar chart (Gini) colored by dominant proportion
p <- p_base +
  geom_fruit(
    geom = geom_col,
    mapping = aes(x = Gini_CellCount, y = label, fill = Dominant_Proportion),
    width = 0.5,
    orientation = "y",
    offset = 0.1
  ) +
  scale_fill_gradientn(
    colors = c("#7B92C7", "#FFD47F", "#F7C1CF"),
    limits = c(0.2, 1),
    oob = squish,
    na.value = "gray90",
    name = "Dominant\nProportion"
  ) +
  geom_tiplab(
    aes(label = Species),
    offset = 0.25,
    align = TRUE,
    linetype = NA,
    size = 3
  )

# Save
ggsave(figure_out, p, width = 10, height = 10, limitsize = FALSE)

# =========================================================================== #
# Figure I: Ridge plot of Gini coefficient by phylum (all phyla with ≥3 species)
# Figure J: Violin plot of Gini coefficient by genus in selected phyla
#           (Bacteroidota, Firmicutes, Proteobacteria, Fusobacteriota,
#            Desulfobacterota), with corresponding bubble plot of functional
#            category composition
# Figure K: Lollipop plot of Gini coefficient for Bacteroides species,
#           colored by normalized abundance
# =========================================================================== #

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(forcats)
library(rstatix)
library(ggpubr)
library(RColorBrewer)
library(tibble)

# Set file paths (replace "path/to/" with actual paths)
metadata_path <- "path/to/metadata.rds"
gini_path     <- "path/to/H_gini.txt"
output_dir    <- "path/to/output/"

# Read data (shared across figures)
meta_data <- readRDS(metadata_path)
results <- read.table(
  file = gini_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE
)

# Remove unknown species
results <- subset(results, Species != "Unknown")

# ============================================================
# Figure I: Ridge plot – Gini by phylum (all phyla)
# ============================================================

# Compute phylum cell count and order by abundance (ascending)
phylum_abundance <- meta_data %>%
  filter(Species != "Unknown") %>%
  group_by(Phylum) %>%
  summarise(CellCount = n(), .groups = "drop")
phylum_ordered <- phylum_abundance %>%
  arrange(CellCount) %>%
  pull(Phylum)

# Map phylum to each species
species_phylum <- meta_data %>%
  distinct(Species, Phylum) %>%
  filter(Species != "Unknown", !is.na(Phylum))

# Prepare data for ridge plot
ridge_data <- results %>%
  inner_join(species_phylum, by = "Species") %>%
  drop_na(Gini_CellCount) %>%
  mutate(Phylum = factor(Phylum, levels = phylum_ordered)) %>%
  group_by(Phylum) %>%
  filter(n_distinct(Species) >= 3) %>%
  ungroup() %>%
  mutate(Phylum = fct_drop(Phylum))

# Perform Dunn's test
dunn_res <- ridge_data %>%
  dunn_test(Gini_CellCount ~ Phylum, p.adjust.method = "BH") %>%
  filter(p.adj < 0.05) %>%
  mutate(signif = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE           ~ "ns"
  ))

# Map factor levels to numeric y-positions
phyla_levels <- levels(ridge_data$Phylum)
y_map <- data.frame(Phylum = phyla_levels, y_num = seq_along(phyla_levels))
ridge_data <- ridge_data %>% left_join(y_map, by = "Phylum")

# Prepare segment data for significant comparisons
seg_data <- NULL
if (nrow(dunn_res) > 0) {
  seg_data <- dunn_res %>%
    left_join(y_map, by = c("group1" = "Phylum")) %>% rename(y1 = y_num) %>%
    left_join(y_map, by = c("group2" = "Phylum")) %>% rename(y2 = y_num) %>%
    mutate(
      y_start = pmin(y1, y2),
      y_end   = pmax(y1, y2),
      y_mid   = (y_start + y_end) / 2
    ) %>%
    arrange(y_mid) %>%
    mutate(offset = 1.05 + 0.03 * (row_number() - 1))
  xlim_max <- max(1.05 + 0.03 * (nrow(seg_data) - 1), 1.2) + 0.05
} else {
  xlim_max <- 1.2
}

# Ridge plot
p_ridge <- ggplot(ridge_data, aes(x = Gini_CellCount, y = y_num, fill = Phylum)) +
  geom_density_ridges(alpha = 0.7, scale = 0.9, rel_min_height = 0.01,
                      quantile_lines = TRUE, quantile_fun = median) +
  scale_fill_viridis_d(option = "D") +
  labs(x = "Gini Coefficient", y = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15, face = "bold")) +
  scale_x_continuous(limits = c(0.4, xlim_max), breaks = seq(0.4, 1, 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(phyla_levels), labels = phyla_levels,
                     expand = expansion(mult = c(0.05, 0.15)))

if (!is.null(seg_data) && nrow(seg_data) > 0) {
  p_ridge <- p_ridge +
    geom_segment(data = seg_data,
                 aes(x = offset, xend = offset, y = y_start, yend = y_end),
                 inherit.aes = FALSE, color = "black", size = 0.5) +
    geom_text(data = seg_data,
              aes(x = offset + 0.02, y = y_mid, label = signif),
              inherit.aes = FALSE, size = 5, hjust = 0, vjust = 0.5)
}

# Save Figure I
ggsave(file.path(output_dir, "I_Gini_Ridgeplot_AllPhyla.pdf"),
       p_ridge, width = 12, height = 8, device = cairo_pdf)

# ============================================================
# Figure J: Violin + Bubble for selected phyla
# ============================================================

# Select five phyla
selected_phyla <- c("Bacteroidota", "Firmicutes", "Proteobacteria",
                    "Fusobacteriota", "Desulfobacterota")

# Prepare data: merge Gini with taxonomy, keep only selected phyla
violin_data <- results %>%
  left_join(meta_data %>% distinct(Species, Phylum, Genus), by = "Species") %>%
  drop_na(Phylum, Genus, Gini_CellCount) %>%
  filter(Phylum %in% selected_phyla) %>%
  group_by(Genus) %>%
  filter(n_distinct(Species) >= 2) %>%
  ungroup()

# Order Phylum factor
violin_data$Phylum <- factor(violin_data$Phylum, levels = selected_phyla)

# Genus color palette (consistent across plots)
all_genera <- unique(violin_data$Genus)
num_genera <- length(all_genera)
genus_colors <- hcl(h = seq(15, 375, length = num_genera + 1),
                    l = 65, c = 100)[1:num_genera]
names(genus_colors) <- all_genera

# Create Genus factor levels in the original order (for x-axis)
violin_data <- violin_data %>%
  mutate(Genus = factor(Genus, levels = all_genera))

# Compute phylum means for reference lines
phylum_means <- violin_data %>%
  group_by(Phylum) %>%
  summarise(mean_gini = mean(Gini_CellCount, na.rm = TRUE), .groups = "drop")

# Violin plot
p_violin <- ggplot(violin_data, aes(x = Genus, y = Gini_CellCount)) +
  geom_violin(aes(fill = Genus), trim = FALSE, width = 0.8) +
  geom_boxplot(fill = "white", width = 0.15, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Species), position = position_jitter(width = 0.15, height = 0),
              size = 1, alpha = 0.6) +
  geom_hline(data = phylum_means, aes(yintercept = mean_gini),
             color = "gray80", linetype = "dashed", size = 0.8, alpha = 0.8) +
  facet_grid(. ~ Phylum, scales = "free_x", space = "free_x") +
  labs(x = "", y = "Gini Coefficient") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = "black")
  ) +
  scale_fill_manual(values = genus_colors) +
  scale_color_viridis_d(option = "D") +
  scale_y_continuous(limits = c(0.2, max(violin_data$Gini_CellCount, na.rm = TRUE)),
                     breaks = seq(0.2, 1.0, by = 0.2),
                     labels = sprintf("%.1f", seq(0.2, 1.0, by = 0.2)),
                     expand = expansion(mult = c(0.05, 0.15)))

ggsave(file.path(output_dir, "J_Gini_Violin_Selected_Phyla.pdf"),
       p_violin, width = 20, height = 5, device = cairo_pdf)

# Bubble plot: functional category composition per genus
# Ensure meta_data has cluster_category (must have been added earlier)
if (!"cluster_category" %in% colnames(meta_data)) {
  stop("meta_data must contain 'cluster_category' column. Run functional annotation first.")
}

selected_genera <- unique(violin_data$Genus)
genus_total <- meta_data %>%
  filter(Genus %in% selected_genera) %>%
  group_by(Genus) %>%
  summarise(total_cells = n(), .groups = "drop")

bubble_data <- meta_data %>%
  filter(Genus %in% selected_genera) %>%
  group_by(Genus, cluster_category) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  left_join(genus_total, by = "Genus") %>%
  mutate(proportion = cell_count / total_cells * 100) %>%
  left_join(meta_data %>% distinct(Genus, Phylum), by = "Genus")

bubble_data$Phylum <- factor(bubble_data$Phylum, levels = selected_phyla)
bubble_data$cluster_category <- factor(bubble_data$cluster_category)
bubble_data$Genus <- factor(bubble_data$Genus, levels = all_genera)

p_bubble <- ggplot(bubble_data, aes(x = Genus, y = cluster_category)) +
  geom_point(aes(size = proportion, color = Genus), alpha = 0.8) +
  facet_grid(. ~ Phylum, scales = "free_x", space = "free_x") +
  scale_color_manual(values = genus_colors, guide = "none") +
  scale_size_continuous(
    name = "Abundance",
    range = c(1, 12),
    breaks = pretty(range(bubble_data$proportion, na.rm = TRUE), n = 5),
    labels = function(x) format(x / 100, digits = 2)
  ) +
  labs(x = "", y = "") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 10, unit = "pt")
  ) +
  guides(size = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5))

ggsave(file.path(output_dir, "J_Bubble_Plot.pdf"),
       p_bubble, width = 20, height = 5, device = cairo_pdf)

# ============================================================
# Figure K: Lollipop plot of Gini for Bacteroides species
# ============================================================

# Species abundance
species_abund <- meta_data %>%
  filter(Species != "Unknown") %>%
  group_by(Species) %>%
  summarise(Abundance = n(), .groups = "drop")

# Bacteroides species list
bacteroides_spp <- meta_data %>%
  filter(Genus == "Bacteroides") %>%
  distinct(Species) %>%
  pull(Species)

# Prepare data
bact_gini <- results %>%
  filter(Species %in% bacteroides_spp, Species != "Unknown") %>%
  select(Species, Gini_CellCount) %>%
  drop_na(Gini_CellCount) %>%
  left_join(species_abund, by = "Species") %>%
  arrange(desc(Gini_CellCount)) %>%
  mutate(Species = factor(Species, levels = Species))

# Normalize abundance (min-max within Bacteroides)
min_ab <- min(bact_gini$Abundance)
max_ab <- max(bact_gini$Abundance)
bact_gini <- bact_gini %>%
  mutate(Abundance_norm = if (max_ab - min_ab > 0) {
    (Abundance - min_ab) / (max_ab - min_ab)
  } else 0.5)

# Lollipop plot
p_lolli <- ggplot(bact_gini, aes(x = Gini_CellCount, y = Species)) +
  geom_segment(aes(x = 0, xend = Gini_CellCount, y = Species, yend = Species),
               color = "black", size = 15) +
  geom_point(aes(color = Abundance_norm), size = 7) +
  labs(x = "Gini Coefficient", y = "", color = "Normalized Abundance") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 15, face = "italic", color = "black"),
    panel.grid.major.y = element_blank()
  ) +
  scale_x_continuous(limits = c(0.6, 0.9), breaks = seq(0.6, 0.9, 0.1)) +
  scale_color_gradientn(colors = brewer.pal(11, "RdGy")[1:5])

ggsave(file.path(output_dir, "K_Bacteroides_Gini.pdf"),
       p_lolli, width = 8, height = max(4, nrow(bact_gini) * 0.3), device = cairo_pdf)