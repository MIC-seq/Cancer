# =========================================================================== #
# Figure A: Dendrogram of clusters (Spearman/Ward.D2) and phylum pie charts 
# per super-cluster.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(Seurat)
library(dendextend)
library(cluster)
library(dplyr)
library(ggplot2)
library(patchwork)

# ------------------------------ File paths ----------------------------------
seurat_obj_path <- "path/to/01merge.rds"
metadata_path  <- "path/to/metadata.rds"
output_dir     <- "path/to/F3"

# -------------------------- Load and prepare data ---------------------------
seurat_obj <- readRDS(seurat_obj_path)
metadata   <- readRDS(metadata_path)
seurat_obj@meta.data <- metadata
rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$barcode

# ------------------- Cluster average expression & dendrogram ----------------
cluster_averages <- AverageExpression(
  seurat_obj, assays = "RNA", slot = "data",
  group.by = "seurat_clusters", verbose = FALSE
)
avg_matrix <- cluster_averages$RNA

# Subset to highly variable genes if available
if (length(VariableFeatures(seurat_obj)) > 0) {
  hv_genes <- VariableFeatures(seurat_obj)
  avg_matrix <- avg_matrix[rownames(avg_matrix) %in% hv_genes, ]
}

# Hierarchical clustering
cor_matrix <- cor(avg_matrix, method = "spearman")
dist_matrix <- as.dist(1 - cor_matrix)
hc <- hclust(dist_matrix, method = "ward.D2")
dend <- as.dendrogram(hc)

# Determine optimal number of clusters using silhouette width
silhouette_scores <- sapply(2:15, function(k) {
  clusters <- cutree(hc, k)
  if (length(unique(clusters)) > 1) {
    sil <- silhouette(clusters, dist_matrix)
    return(mean(sil[, 3]))
  }
  return(0)
})
optimal_k <- which.max(silhouette_scores) + 1
cluster_groups <- cutree(hc, k = optimal_k)

# -------------- Phylum grouping and super-cluster assignment ----------------
# Merge Firmicutes sub-phyla
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(Phylum = case_when(
    Phylum %in% c("Firmicutes_A", "Firmicutes_B", "Firmicutes_C", "Firmicutes_G") ~ "Firmicutes",
    TRUE ~ Phylum
  ))

# Label phyla with fewer than 50 cells as "Others"
phylum_counts_meta <- seurat_obj@meta.data %>%
  group_by(Phylum) %>%
  summarise(cellnum = n(), .groups = "drop")
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  left_join(phylum_counts_meta, by = "Phylum") %>%
  mutate(Phylum = ifelse(cellnum < 50, "Others", Phylum)) %>%
  select(-cellnum)

# Add super-cluster labels
group_map <- data.frame(
  seurat_clusters = names(cluster_groups),
  super_cluster = as.character(cluster_groups)
)
seurat_obj@meta.data$seurat_clusters_char <- as.character(seurat_obj@meta.data$seurat_clusters)
seurat_obj@meta.data$super_cluster <- group_map$super_cluster[
  match(seurat_obj@meta.data$seurat_clusters_char, group_map$seurat_clusters)
]

# ------------ Compute phylum proportions per super-cluster ------------------
meta_data <- seurat_obj@meta.data

phylum_counts <- meta_data %>%
  filter(!is.na(Phylum) & !is.na(super_cluster)) %>%
  group_by(super_cluster, Phylum) %>%
  summarise(cell_count = n(), .groups = "drop")

total_cells <- phylum_counts %>%
  group_by(super_cluster) %>%
  summarise(total = sum(cell_count), .groups = "drop")

phylum_prop <- phylum_counts %>%
  left_join(total_cells, by = "super_cluster") %>%
  mutate(proportion = cell_count / total * 100)

# Phylum colour palette
phylum_colors <- c(
  'Bacteroidota'      = '#66C2A5',
  'Desulfobacterota'  = '#E78AC3',
  'Firmicutes'        = '#B3B3B3',
  'Fusobacteriota'    = '#FC8D62',
  'Others'            = '#FFD92F',
  'Proteobacteria'    = '#8DA0CB',
  'Synergistota'      = '#E5C494',
  'Verrucomicrobiota' = '#A6D854'
)

# ------------------------ Pie chart function --------------------------------
plot_pie <- function(super_id, data, color_map) {
  sub_data <- data %>% filter(super_cluster == super_id)
  sub_data <- sub_data %>%
    arrange(desc(proportion)) %>%
    mutate(
      label_legend = paste0(Phylum, " (", round(proportion, 1), "%)"),
      color = color_map[Phylum]
    )
  
  ggplot(sub_data, aes(x = "", y = proportion, fill = Phylum)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.2) +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(
      values = setNames(sub_data$color, sub_data$Phylum),
      labels = setNames(sub_data$label_legend, sub_data$Phylum),
      name = "Phylum"
    ) +
    labs(title = paste("Super Cluster", super_id)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
      legend.title = element_text(size = 10, color = "black"),
      legend.text = element_text(size = 8, color = "black"),
      legend.position = "right",
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    ) +
    guides(fill = guide_legend(title = "Phylum", ncol = 1))
}

# ---------------------- Generate and save pie charts ------------------------
super_list <- sort(unique(phylum_prop$super_cluster))
if (length(super_list) > 3) {
  super_list <- super_list[1:3]   # Plot only the first three if more exist
}

pie_plots <- lapply(super_list, function(s) plot_pie(s, phylum_prop, phylum_colors))

combined_plot <- wrap_plots(pie_plots, nrow = 1) +
  plot_annotation(
    title = "Phylum composition per super cluster",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
  )

ggsave(file.path(output_dir, "A_phylum_pie.pdf"),
       plot = combined_plot, width = 12, height = 5, dpi = 300,
       limitsize = FALSE, bg = "transparent")

# ------------------ Coloured dendrogram ------------------------------------
# Colour mapping: functional category -> cluster numbers
color_mapping <- list(
  "#3FA9F5" = c(1, 2, 3, 5, 13, 15, 17, 21, 26, 35),
  "#BDCCD4" = c(4, 7, 19, 24, 29, 39),
  "#FF931E" = c(6, 9, 10, 11, 18, 20, 22, 25, 28, 30, 31, 33, 34, 38, 40, 41, 43, 45, 46),
  "#FF7BAC" = c(8, 14, 23, 37, 42, 44, 49),
  "#00FFFF" = c(16, 27, 32, 36, 47),
  "#7AC943" = c(12),
  "#FCEE21" = c(48),
  "#C69C6D" = c(0)
)

# Assign a colour to each leaf label
sample_labels <- as.numeric(labels(dend))
color_vector <- rep("black", length(sample_labels))
for (color_name in names(color_mapping)) {
  indices <- which(sample_labels %in% color_mapping[[color_name]])
  color_vector[indices] <- color_name
}
label_to_color <- setNames(color_vector, sample_labels)

# Cell count per cluster (for internal node colour resolution)
cluster_counts <- table(meta_data$seurat_clusters)

# Recursive colour assignment
assign_colors_to_nodes <- function(node, label_to_color_map, cluster_counts) {
  if (is.leaf(node)) {
    node_label <- as.numeric(attr(node, "label"))
    color <- label_to_color_map[as.character(node_label)]
    count <- cluster_counts[as.character(node_label)]
    if (is.na(count)) count <- 0
    return(list(node = node, color = color, count = count))
  } else {
    left_result  <- assign_colors_to_nodes(node[[1]], label_to_color_map, cluster_counts)
    right_result <- assign_colors_to_nodes(node[[2]], label_to_color_map, cluster_counts)
    
    left_color  <- left_result$color
    right_color <- right_result$color
    left_count  <- left_result$count
    right_count <- right_result$count
    
    if (left_color != right_color) {
      node_color <- if (left_count > right_count) left_color else right_color
    } else {
      node_color <- left_color
    }
    node_count <- left_count + right_count
    
    node[[1]] <- left_result$node
    node[[2]] <- right_result$node
    
    attr(node, "edgePar") <- list(col = node_color, lwd = 2)
    attr(node[[1]], "edgePar") <- list(col = left_color, lwd = 2)
    attr(node[[2]], "edgePar") <- list(col = right_color, lwd = 2)
    
    return(list(node = node, color = node_color, count = node_count))
  }
}

result <- assign_colors_to_nodes(dend, label_to_color, cluster_counts)
dend_colored <- result$node
dend_colored <- set(dend_colored, "labels_cex", 0.8)
dend_colored <- set(dend_colored, "labels_col", "black")

# Save dendrogram
pdf(file.path(output_dir, "A_dendrogram.pdf"),
    width = 3000/300, height = 1600/300,
    bg = "transparent", family = "serif")
par(mar = c(4, 8, 4, 10), las = 1)
height_range <- range(attr(dend_colored, "height"))
y_ticks <- seq(0, ceiling(height_range[2]/0.1)*0.1, by = 0.1)
plot(dend_colored, main = "", xlab = "", ylab = "", sub = "",
     cex.main = 1.2, cex.lab = 1.0, axes = FALSE)
axis(side = 4, at = y_ticks, labels = sprintf("%.1f", y_ticks),
     cex.axis = 0.9, las = 1)
dev.off()

# ---------------------------- Colour legend ---------------------------------
color_names <- c(
  "Core Energy Metabolism & Central Carbon Processing",
  "Stress Response & Cellular Maintenance",
  "Substrate Utilization & Specialized Fermentation",
  "Beneficial Metabolite Production (Host-Microbe Interactions)",
  "Cell Growth, Division & Motility",
  "Defense & Resistance",
  "Cofactor & Micronutrient Procurement",
  "Unspecified & Generalist Functions"
)
color_values <- c("#3FA9F5", "#BDCCD4", "#FF931E", "#FF7BAC",
                  "#00FFFF", "#7AC943", "#FCEE21", "#C69C6D")

pdf(file.path(output_dir, "A_legend.pdf"),
    width = 2000/300, height = 1000/300,
    bg = "transparent")
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(0, 1), ylim = c(0, 1))
legend("center",
       legend = color_names,
       fill = color_values,
       border = NA,
       bty = "n",
       cex = 1.2,
       xpd = TRUE,
       title = "Functional Categories",
       title.cex = 1.3,
       title.font = 2,
       text.width = 0.7,
       y.intersp = 1.3)
dev.off()




# =========================================================================== #
# Figure B: Dot plot of clinical feature associations with cluster proportions.
# For each clinical feature (coded 0/1), log2(ratio of mean proportions) and
# Wilcoxon test significance are shown per cluster.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------ File paths ----------------------------------
seurat_obj_path <- "path/to/01merge.rds"
metadata_path  <- "path/to/metadata.rds"
output_dir     <- "path/to/F3"
wide_prop_path <- file.path(output_dir, "B_cell_proportion_wide.csv")
results_path   <- file.path(output_dir, "B_logistic_regression_all_results_baseline.csv")
figure_path    <- file.path(output_dir, "B_clinical_features_logistic_dotplot_with_grid_baseline.pdf")

# -------------------------- Load and prepare data ---------------------------
seurat_obj <- readRDS(seurat_obj_path)
metadata   <- readRDS(metadata_path)
seurat_obj@meta.data <- metadata
rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$barcode

# Subset to Pre-treatment samples
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  filter(treatment_condition == "Pre")

# -------------------------- 1. Parameters -----------------------------------
clinical_features <- c("PRPDSD", "MMR_situation", "hepatic_metastases",
                       "Gender", "PPI", "Antibiotics")

coding_map <- list(
  PRPDSD = list(keep = c("PR", "PD"),
                mapping = c("PD" = 0, "PR" = 1)),
  MMR_situation = list(keep = c("dMMR", "pMMR"),
                       mapping = c("pMMR" = 0, "dMMR" = 1)),
  hepatic_metastases = list(keep = c("Metastatic", "Non-metastatic"),
                            mapping = c("Non-metastatic" = 0, "Metastatic" = 1)),
  Gender = list(keep = c("Male", "Female"),
                mapping = c("Female" = 0, "Male" = 1)),
  PPI = list(keep = c("PPI", "Non-PPI"),
             mapping = c("Non-PPI" = 0, "PPI" = 1)),
  Antibiotics = list(keep = c("Antibiotic", "Non-Antibiotic"),
                     mapping = c("Non-Antibiotic" = 0, "Antibiotic" = 1))
)

clinical_feature_order <- c("Gender", "Antibiotics", "PPI", "MMR_situation",
                            "hepatic_metastases", "PRPDSD")
alpha <- 0.05

# -------------------- 2. Extract sample-level clinical info -----------------
meta <- seurat_obj@meta.data

sample_clinical <- meta %>%
  dplyr::select(orig.ident, all_of(clinical_features)) %>%
  distinct() %>%
  group_by(orig.ident) %>%
  summarise(across(everything(), ~ {
    uniq_vals <- unique(na.omit(.x))
    if (length(uniq_vals) > 1) {
      warning(paste("Sample", cur_group()$orig.ident, "has multiple values for", cur_column()))
      return(NA)
    } else if (length(uniq_vals) == 0) {
      return(NA)
    } else {
      return(uniq_vals[1])
    }
  })) %>%
  ungroup()

# ----------------------- 3. Encode clinical features ------------------------
regression_data_list <- list()
for (feature in clinical_features) {
  cfg <- coding_map[[feature]]
  if (is.null(cfg)) next
  df <- sample_clinical %>%
    dplyr::select(orig.ident, !!feature := all_of(feature)) %>%
    filter(!is.na(!!sym(feature))) %>%
    filter(!!sym(feature) %in% cfg$keep)
  df$binary <- as.integer(cfg$mapping[as.character(df[[feature]])])
  regression_data_list[[feature]] <- df
}

# --------------- 4. Compute cluster proportions per sample ------------------
sample_total <- meta %>%
  group_by(orig.ident) %>%
  summarise(total_cells = n())

cluster_counts <- meta %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(cluster_cells = n(), .groups = "drop")

cluster_prop <- cluster_counts %>%
  left_join(sample_total, by = "orig.ident") %>%
  mutate(proportion = cluster_cells / total_cells) %>%
  dplyr::select(orig.ident, cluster = seurat_clusters, proportion)

all_samples  <- unique(meta$orig.ident)
all_clusters <- sort(unique(meta$seurat_clusters))

complete_prop <- expand.grid(orig.ident = all_samples,
                             cluster = all_clusters,
                             stringsAsFactors = FALSE) %>%
  left_join(cluster_prop, by = c("orig.ident", "cluster")) %>%
  mutate(proportion = ifelse(is.na(proportion), 0, proportion))

wide_prop <- complete_prop %>%
  pivot_wider(id_cols = orig.ident,
              names_from = cluster,
              values_from = proportion,
              values_fill = 0)

write.csv(wide_prop, file = wide_prop_path, row.names = FALSE)

# -------- 5. Compute log2(Ratio) and Wilcoxon test significance -------------
results <- list()

for (feat in names(regression_data_list)) {
  clin_df <- regression_data_list[[feat]]
  
  for (cl in all_clusters) {
    prop_df <- complete_prop %>%
      filter(cluster == cl, orig.ident %in% clin_df$orig.ident) %>%
      dplyr::select(orig.ident, proportion)
    
    model_df <- clin_df %>%
      left_join(prop_df, by = "orig.ident")
    
    # Mean proportions per binary group
    group_means <- model_df %>%
      group_by(binary) %>%
      summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop")
    
    mean1 <- group_means$mean_prop[group_means$binary == 1]
    mean0 <- group_means$mean_prop[group_means$binary == 0]
    
    if (length(mean1) == 0) mean1 <- NA
    if (length(mean0) == 0) mean0 <- NA
    
    ratio_val <- ifelse(!is.na(mean1) & !is.na(mean0) & mean0 > 0,
                        mean1 / mean0,
                        NA)
    log2_ratio <- log2(ratio_val)
    # Cap infinite values at ±4
    log2_ratio <- ifelse(is.infinite(log2_ratio),
                         sign(log2_ratio) * 4,
                         log2_ratio)
    
    # Wilcoxon test
    pr_props <- model_df$proportion[model_df$binary == 1]
    pd_props <- model_df$proportion[model_df$binary == 0]
    
    if (length(pr_props) >= 3 & length(pd_props) >= 3) {
      wt <- wilcox.test(pr_props, pd_props)
      p_val <- wt$p.value
    } else {
      p_val <- NA
    }
    
    sig_symbol <- case_when(
      is.na(p_val) ~ "",
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE ~ ""
    )
    
    results[[paste(feat, cl, sep = "_")]] <- data.frame(
      Clinical_Feature = feat,
      Cluster = as.character(cl),
      Coefficient = log2_ratio,
      P_value = p_val,
      Significance = sig_symbol,
      stringsAsFactors = FALSE
    )
  }
}

all_results <- bind_rows(results)
write.csv(all_results, file = results_path, row.names = FALSE)

# ------------------------ 6. Prepare heatmap data ---------------------------
# Order clusters using the dendrogram from Figure A (hc must exist)
if (exists("hc") && !is.null(hc)) {
  dend_order <- labels(as.dendrogram(hc))
  cluster_order <- intersect(dend_order, all_clusters)
  cluster_order <- c(cluster_order, setdiff(all_clusters, cluster_order))
} else {
  cluster_order <- as.character(sort(all_clusters))
}

all_results$Cluster <- factor(all_results$Cluster, levels = cluster_order)
all_results$Clinical_Feature <- factor(all_results$Clinical_Feature,
                                       levels = clinical_feature_order)

# Y-axis labels with category information
y_labels <- sapply(clinical_feature_order, function(feat) {
  cfg <- coding_map[[feat]]
  if (!is.null(cfg)) {
    cat_1 <- names(cfg$mapping)[cfg$mapping == 1]
    cat_0 <- names(cfg$mapping)[cfg$mapping == 0]
    if (length(cat_1) == 1 && length(cat_0) == 1) {
      categories <- paste(cat_1, cat_0, sep = " / ")
    } else {
      categories <- paste(cfg$keep, collapse = " / ")
    }
    label <- paste0(feat, " (", categories, ")")
  } else {
    label <- feat
  }
  return(label)
})
names(y_labels) <- clinical_feature_order

# Size scaling (log2 ratio absolute value capped at 4)
custom_limits <- c(-4, 4)
min_abs <- 1
max_abs <- 4
min_size <- 2
max_size <- 10

all_results <- all_results %>%
  mutate(abs_coef = abs(Coefficient),
         size_value = case_when(
           abs_coef <= min_abs ~ min_size,
           abs_coef >= max_abs ~ max_size,
           TRUE ~ scales::rescale(abs_coef, to = c(min_size, max_size),
                                  from = c(min_abs, max_abs))
         ))

# ------------------------------ 7. Draw dot plot ----------------------------
p <- ggplot(all_results, aes(x = Cluster, y = Clinical_Feature)) +
  geom_tile(color = "gray80", fill = NA, size = 0.5) +
  geom_point(aes(color = Coefficient, size = size_value), alpha = 0.8) +
  geom_text(aes(label = Significance), size = 5, color = "black",
            vjust = 0.5, hjust = 0.5) +
  scale_color_gradientn(
    colours = c("#5a68ae", "white", "#f1656d"),
    values = scales::rescale(c(-4, 0, 4)),
    limits = custom_limits,
    oob = scales::squish,
    na.value = "grey80",
    name = "log2(Ratio)"
  ) +
  scale_size_identity(
    name = "|log2(Ratio)|",
    breaks = c(min_size, (min_size + max_size) / 2, max_size),
    labels = c(paste0("≤ ", min_abs),
               paste0(round((min_abs + max_abs) / 2, 1)),
               paste0("≥ ", max_abs)),
    guide = guide_legend(override.aes = list(shape = 16))
  ) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0), position = "right",
                   labels = y_labels) +
  labs(title = "", subtitle = "") +
  theme_minimal() +
  theme(
    text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
                               colour = "black"),
    axis.text.x.top = element_text(margin = margin(b = 5)),
    axis.text.y = element_text(colour = "black", hjust = 0, size = 10),
    axis.text.y.right = element_text(margin = margin(l = 5), size = 10),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(colour = "black"),
    legend.title = element_text(colour = "black")
  )

ggsave(figure_path, plot = p, width = 23, height = 3, dpi = 300,
       limitsize = FALSE)



# =========================================================================== #
# Figure C: Heatmap of tissue distribution (R/oe) for PRPDSD × treatment 
# combinations across clusters.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(Startrac)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(scales)

# ------------------------------ File paths ----------------------------------
output_dir  <- "path/to/F3"
output_path <- file.path(output_dir, "C.pdf")

# -------------------------- Prepare metadata --------------------------------
# Create combined PRPDSD_treatment variable
seurat_obj@meta.data$PRPDSDtreatment <- paste(
  seurat_obj@meta.data$PRPDSD,
  seurat_obj@meta.data$treatment_condition,
  sep = "_"
)

C_metadata <- seurat_obj@meta.data
C_metadata <- subset(C_metadata, PRPDSD != "Unknow")
C_metadata <- subset(C_metadata, treatment_condition != "Unknow")
C_metadata <- subset(C_metadata, treatment_condition != "Post2")

# ------------ Calculate tissue distribution (R/oe) with chi-square ----------
R_oe <- calTissueDist(C_metadata,
                      byPatient = FALSE,
                      colname.cluster = "PRPDSDtreatment",
                      colname.patient = "orig.ident",
                      colname.tissue = "seurat_clusters",
                      method = "chisq",
                      min.rowSum = 0)

# Melt matrix for ggplot
R_oe_matrix <- as.matrix(R_oe)
R_oe_melted <- melt(R_oe_matrix)
colnames(R_oe_melted) <- c("Feature_Value", "Cluster", "Value")

# Assign text labels based on R/oe thresholds
R_oe_melted$Label <- ifelse(R_oe_melted$Value > 1, "+++",
                     ifelse(R_oe_melted$Value > 0.8 & R_oe_melted$Value <= 1, "++",
                     ifelse(R_oe_melted$Value >= 0.2 & R_oe_melted$Value <= 0.8, "+",
                     ifelse(R_oe_melted$Value > 0 & R_oe_melted$Value < 0.2, "+/−",
                     ifelse(R_oe_melted$Value == 0, "-", NA)))))

# Colour palette (BuPu)
my_palette <- colorRampPalette(brewer.pal(n = 9, name = "BuPu")[1:5])(100)
color_limits <- c(0, 2)

# Order clusters by dendrogram from Figure A (hc must exist)
dend_order <- labels(as.dendrogram(hc))
R_oe_melted$Cluster <- factor(R_oe_melted$Cluster, levels = dend_order)

# Order rows (feature values) as specified
feature_order <- c("PD_Post", "SD_Post", "PR_Post", "PD_Pre", "SD_Pre", "PR_Pre")
R_oe_melted$Feature_Value <- factor(R_oe_melted$Feature_Value, levels = feature_order)

# ------------------------------ Draw heatmap ---------------------------------
p <- ggplot(R_oe_melted, aes(x = Cluster, y = Feature_Value, fill = Value)) +
  geom_tile(color = "gray70") +
  geom_text(aes(label = Label), color = "black", size = 3) +
  scale_fill_gradientn(
    colours = my_palette,
    name = "R/oe",
    limits = color_limits,
    oob = scales::squish
  ) +
  scale_y_discrete(
    expand = c(0, 0),
    position = "right",
    labels = function(x) gsub("_", " ", x)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 10),
    axis.text.y.right = element_text(colour = "black", hjust = 0, size = 12),
    axis.text.y.left = element_blank(),
    axis.title = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    panel.grid = element_blank(),
    legend.text = element_text(colour = "black"),
    legend.title = element_text(colour = "black")
  )

ggsave(output_path, plot = p, width = 16, height = 6 * 0.3, dpi = 300,
       limitsize = FALSE)



# =========================================================================== #
# Figure D: Heatmap of cell-abundance differences between PRPDSD groups and
# treatment conditions (uncorrected Wilcoxon p-values).
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ------------------------------ File paths ----------------------------------
seurat_obj_path <- "path/to/01merge.rds"
metadata_path  <- "path/to/metadata.rds"
output_dir     <- "path/to/F3"
csv_path       <- file.path(output_dir, "D_cell_abundance_stats.csv")
figure_path    <- file.path(output_dir, "D_cellassample_abundance.pdf")

# -------------------------- Load and prepare data ---------------------------
seurat_obj <- readRDS(seurat_obj_path)
metadata   <- readRDS(metadata_path)
seurat_obj@meta.data <- metadata
rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$barcode

# Filter out "Unknow" categories
seurat_known <- subset(seurat_obj,
                       subset = PRPDSD != "Unknow" & treatment_condition != "Unknow")

# ----------- Function: cell-abundance contrast matrix (uncorrected) ---------
create_gene_expression_contrast_matrix <- function(seurat_obj) {
  meta_data <- seurat_obj@meta.data
  clusters  <- sort(unique(meta_data$seurat_clusters))

  # Initialize result matrices
  results_matrix <- matrix(0, nrow = 6, ncol = length(clusters))
  pvalue_matrix  <- matrix(1, nrow = 6, ncol = length(clusters))
  rownames(results_matrix) <- rownames(pvalue_matrix) <-
    c("PR_vs_PD", "PR_vs_SD", "SD_vs_PD",
      "PR_Pre_vs_Post", "SD_Pre_vs_Post", "PD_Pre_vs_Post")
  colnames(results_matrix) <- colnames(pvalue_matrix) <- clusters

  # Total cells per sample
  sample_total <- meta_data %>%
    group_by(orig.ident) %>%
    summarise(total_n = n(), .groups = "drop")

  # Process each cluster
  for (cluster in clusters) {
    cluster_counts <- meta_data %>%
      filter(seurat_clusters == cluster) %>%
      group_by(orig.ident) %>%
      summarise(n = n(), .groups = "drop")

    all_samples <- unique(meta_data$orig.ident)
    full_counts <- expand.grid(orig.ident = all_samples) %>%
      left_join(cluster_counts, by = "orig.ident") %>%
      mutate(n = ifelse(is.na(n), 0, n)) %>%
      left_join(sample_total, by = "orig.ident") %>%
      mutate(prop = n / total_n)   # raw proportions, no pseudocount

    sample_info <- meta_data %>%
      select(orig.ident, PRPDSD, treatment_condition) %>%
      distinct()

    prop_df <- full_counts %>%
      left_join(sample_info, by = "orig.ident")

    # Extract group vectors
    PR_samps  <- prop_df %>% filter(PRPDSD == "PR") %>% pull(prop)
    PD_samps  <- prop_df %>% filter(PRPDSD == "PD") %>% pull(prop)
    SD_samps  <- prop_df %>% filter(PRPDSD == "SD") %>% pull(prop)
    PR_Pre    <- prop_df %>% filter(PRPDSD == "PR", treatment_condition == "Pre")  %>% pull(prop)
    PR_Post   <- prop_df %>% filter(PRPDSD == "PR", treatment_condition == "Post") %>% pull(prop)
    SD_Pre    <- prop_df %>% filter(PRPDSD == "SD", treatment_condition == "Pre")  %>% pull(prop)
    SD_Post   <- prop_df %>% filter(PRPDSD == "SD", treatment_condition == "Post") %>% pull(prop)
    PD_Pre    <- prop_df %>% filter(PRPDSD == "PD", treatment_condition == "Pre")  %>% pull(prop)
    PD_Post   <- prop_df %>% filter(PRPDSD == "PD", treatment_condition == "Post") %>% pull(prop)

    # Safe Wilcoxon test
    safe_wilcox <- function(x, y) {
      if (length(x) < 2 | length(y) < 2) return(list(diff = 0, p = 1))
      if (sd(x) == 0 & sd(y) == 0) return(list(diff = 0, p = 1))
      tryCatch({
        test <- wilcox.test(x, y)
        list(diff = median(x) - median(y), p = test$p.value)
      }, error = function(e) list(diff = 0, p = 1))
    }

    res <- safe_wilcox(PR_samps, PD_samps)
    results_matrix["PR_vs_PD", as.character(cluster)] <- res$diff
    pvalue_matrix["PR_vs_PD", as.character(cluster)]  <- res$p

    res <- safe_wilcox(PR_samps, SD_samps)
    results_matrix["PR_vs_SD", as.character(cluster)] <- res$diff
    pvalue_matrix["PR_vs_SD", as.character(cluster)]  <- res$p

    res <- safe_wilcox(SD_samps, PD_samps)
    results_matrix["SD_vs_PD", as.character(cluster)] <- res$diff
    pvalue_matrix["SD_vs_PD", as.character(cluster)]  <- res$p

    res <- safe_wilcox(PR_Pre, PR_Post)
    results_matrix["PR_Pre_vs_Post", as.character(cluster)] <- res$diff
    pvalue_matrix["PR_Pre_vs_Post", as.character(cluster)]  <- res$p

    res <- safe_wilcox(SD_Pre, SD_Post)
    results_matrix["SD_Pre_vs_Post", as.character(cluster)] <- res$diff
    pvalue_matrix["SD_Pre_vs_Post", as.character(cluster)]  <- res$p

    res <- safe_wilcox(PD_Pre, PD_Post)
    results_matrix["PD_Pre_vs_Post", as.character(cluster)] <- res$diff
    pvalue_matrix["PD_Pre_vs_Post", as.character(cluster)]  <- res$p
  }

  return(list(diff_matrix = results_matrix, pvalue_matrix = pvalue_matrix))
}

# ------------------ Run abundance analysis (uncorrected) ---------------------
analysis_results <- create_gene_expression_contrast_matrix(seurat_known)
diff_matrix   <- analysis_results$diff_matrix
pvalue_matrix <- analysis_results$pvalue_matrix

# Prepare long-format data for plotting
diff_df <- as.data.frame(diff_matrix) %>%
  rownames_to_column(var = "Comparison") %>%
  pivot_longer(cols = -Comparison, names_to = "Cluster", values_to = "Difference") %>%
  mutate(Cluster    = factor(Cluster, levels = colnames(diff_matrix)),
         Comparison = factor(Comparison, levels = rownames(diff_matrix)))

pvalue_df <- as.data.frame(pvalue_matrix) %>%
  rownames_to_column(var = "Comparison") %>%
  pivot_longer(cols = -Comparison, names_to = "Cluster", values_to = "p_value") %>%
  mutate(Cluster    = factor(Cluster, levels = colnames(pvalue_matrix)),
         Comparison = factor(Comparison, levels = rownames(pvalue_matrix)))

plot_data <- diff_df %>%
  left_join(pvalue_df, by = c("Comparison", "Cluster"))

# Add significance labels (uncorrected)
plot_data <- plot_data %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    sig_color  = ifelse(significance == "ns", "gray50", "black"),
    sig_factor = factor(significance, levels = c("***", "**", "*", "ns", "NA"))
  )

write.csv(plot_data, file = csv_path, row.names = FALSE)

# --------------- Set cluster order from dendrogram (Figure A) ---------------
dend_order <- labels(as.dendrogram(hc))
plot_data$Cluster <- factor(plot_data$Cluster, levels = dend_order)

# Colour limits and row order
color_limits <- c(-0.01, 0.01)
Comparison_order <- rev(c("PR_vs_PD", "PR_vs_SD", "SD_vs_PD",
                          "PR_Pre_vs_Post", "SD_Pre_vs_Post", "PD_Pre_vs_Post"))
plot_data$Comparison <- factor(plot_data$Comparison, levels = Comparison_order)

# ------------------------------ Draw heatmap ---------------------------------
p <- ggplot(plot_data, aes(x = Cluster, y = Comparison, fill = Difference)) +
  geom_tile(color = "gray70") +
  geom_text(aes(label = ifelse(significance != "ns", significance, "")),
            color = "black", size = 3) +
  scale_fill_gradient2(
    low = "#313695",
    mid = "white",
    high = "#ff7f00",
    midpoint = 0,
    name = "Abundance Difference",
    limits = color_limits,
    oob = scales::squish
  ) +
  scale_y_discrete(expand = c(0, 0), position = "right",
                   labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Cell Abundance Differences Between Groups (uncorrected p-values)") +
  theme_minimal() +
  theme(
    text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 10),
    axis.text.y.right = element_text(colour = "black", hjust = 0, size = 12),
    axis.text.y.left = element_blank(),
    axis.title = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    panel.grid = element_blank(),
    legend.text = element_text(colour = "black"),
    legend.title = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave(figure_path, plot = p, width = 16, height = 6 * 0.3, dpi = 300,
       limitsize = FALSE)



# =========================================================================== #
# Figure E: Heatmap of Spearman correlations between cluster proportions and
# targeted serum metabolites (left), and log2 fold-change (PR vs. PD) of the
# same metabolites in pre-treatment samples (right).
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(reshape2)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(readxl)
library(sva)
library(scales)

# ------------------------------ File paths ----------------------------------
output_dir     <- "path/to/F3"
cor_r_path     <- file.path(output_dir, "E_01cluster_metabolome_correlation_R0103.csv")
cor_p_path     <- file.path(output_dir, "E_01cluster_metabolome_correlation_P0103.csv")
metab1_path    <- "path/to/metabolome/01.xlsx"
metab2_path    <- "path/to/metabolome/03.xlsx"
clinical_path  <- "path/to/Clinical_information.rds"
out_heatmap    <- file.path(output_dir, "E_metabolome.pdf")
out_log2fc     <- file.path(output_dir, "E_metabolome_log2FC_column.pdf")

# ================= Part 1: Correlation heatmap ==============================

# Read correlation and p-value matrices
df0103R <- read.csv(cor_r_path, header = TRUE, row.names = 1, check.names = FALSE)
df0103P <- read.csv(cor_p_path, header = TRUE, row.names = 1, check.names = FALSE)

# Melt to long format
cor_df <- melt(as.matrix(df0103R), varnames = c("Cluster", "Metabolite"),
               value.name = "Correlation")
p_df   <- melt(as.matrix(df0103P), varnames = c("Cluster", "Metabolite"),
               value.name = "Pvalue")

# Merge and assign significance stars
plot_df <- merge(cor_df, p_df, by = c("Cluster", "Metabolite"))
plot_df$Significance <- ifelse(plot_df$Pvalue <= 0.005, "***",
                        ifelse(plot_df$Pvalue <= 0.01, "**",
                        ifelse(plot_df$Pvalue <= 0.05, "*", "")))

# Remove metabolites without any significant correlation
plot_df <- plot_df %>%
  group_by(Metabolite) %>%
  filter(any(Significance != "")) %>%
  ungroup()

# Filter to a curated set of target metabolites
target_metabolites <- c(
  "L-Tryptophan", "Butyric acid", "Succinic acid", "L-Arginine", "serotonin",
  "Kynurenine", "Indole", "Indole-3-acetic acid", "3-Hydroxybutyric acid",
  "Lactic acid", "Phosphoenolpyruvic acid", "2-Oxobutanoic acid",
  "L-Citrulline", "L-Ornithine", "5-Hydroxytryptophan", "Tryptamine",
  "Kynurenic acid", "Anthranilic acid", "Indoxylsulfate", "Melatonin",
  "3-Indolepropionic acid", "2-Hydroxybutyric acid", "Creatine", "Creatinine",
  "Urea", "Fumaric acid", "Quinolinic acid", "Xanthurenic acid",
  "Tryptophol", "3-Guanidinopropionic acid", "Itaconic acid"
)

plot_df <- plot_df %>% filter(Metabolite %in% target_metabolites)

# Clean cluster names (remove "cluster" prefix)
plot_df$Cluster <- str_remove_all(plot_df$Cluster, "cluster")

# Set cluster order (must match the dendrogram from Figure A)
dend_order <- c("22", "18", "31", "45", "38", "13", "23", "33", "11", "9",
                "16", "48", "30", "8", "32", "47", "42", "34", "40", "39",
                "15", "6", "10", "3", "44", "29", "2", "1", "7", "21", "41",
                "37", "14", "35", "17", "5", "27", "20", "43", "28", "46",
                "49", "12", "4", "24", "19", "26", "25", "0", "36")
plot_df$Cluster <- factor(plot_df$Cluster, levels = dend_order)

# Metabolite order (reverse of target list for y-axis)
plot_df$Metabolite <- factor(plot_df$Metabolite, levels = rev(target_metabolites))

# Custom colour palette for correlation heatmap
my_colors <- c("#8c0152", "#dc77ac", "#fadeed", "#fdf5f2", "white",
               "#f5faf3", "#e5f3ce", "#b7df84", "#276419")

# Draw correlation heatmap
p_cor <- ggplot(plot_df, aes(x = Cluster, y = Metabolite, fill = Correlation)) +
  geom_tile(color = "gray70") +
  geom_text(aes(label = ifelse(Significance != "", Significance, "")),
            color = "black", size = 3) +
  scale_fill_gradientn(
    colours = rev(my_colors),
    limits = c(-0.5, 0.5),
    oob = scales::squish,
    name = "Correlation coefficient"
  ) +
  scale_y_discrete(expand = c(0, 0), position = "right",
                   labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 10),
    axis.text.y.right = element_text(colour = "black", hjust = 0, size = 12),
    axis.text.y.left = element_blank(),
    axis.title = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    panel.grid = element_blank(),
    legend.text = element_text(colour = "black"),
    legend.title = element_text(colour = "black")
  )

ggsave(out_heatmap, plot = p_cor, width = 16, height = 25 * 0.3, dpi = 300,
       limitsize = FALSE)


# ================= Part 2: log2 fold-change (PR/PD) =========================

# Load and harmonize two metabolite datasets
metab01 <- read_excel(metab1_path)
metab02 <- read_excel(metab2_path)

all_cols <- union(colnames(metab01), colnames(metab02))
metab01 <- metab01 %>%
  mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>%
  select(all_of(all_cols))
metab02 <- metab02 %>%
  mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>%
  select(all_of(all_cols))

# Identify common metabolites and combine with batch information
common_metabs <- intersect(names(metab01)[-1], names(metab02)[-1])
combined <- bind_rows(
  metab01 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 1),
  metab02 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 2)
) %>%
  mutate(Batch = as.integer(Batch),
         BioID = ifelse(Samples %in% intersect(metab01$Samples, metab02$Samples),
                        Samples, paste0("uniq_", row_number())))

# ComBat batch correction (log10-transformed internally)
combat_matrix <- combined %>%
  select(all_of(common_metabs)) %>%
  as.matrix() %>%
  {log10(. + 1)} %>%
  t() %>%
  ComBat(batch = combined$Batch,
         mod = model.matrix(~1, data = combined),
         par.prior = TRUE) %>%
  t() %>%
  {10^. - 1} %>%
  {.[. < 0] <- 0; .}

# Consolidate duplicate samples and create final matrix
final_result <- combat_matrix %>%
  as.data.frame() %>%
  mutate(Sample = combined$Samples,
         Batch = combined$Batch,
         BioID = combined$BioID) %>%
  group_by(BioID) %>%
  summarise(
    Sample = first(Sample),
    across(all_of(common_metabs), mean),
    .groups = "drop"
  ) %>%
  arrange(factor(Sample, levels = c(metab01$Samples,
                                    setdiff(metab02$Samples, metab01$Samples)))) %>%
  select(Sample, all_of(common_metabs)) %>%
  column_to_rownames("Sample")

# Map metabolite names to codes
mapping_list <- data.frame(
  original_name = colnames(final_result),
  new_name = paste0("metabolite_", seq_len(ncol(final_result)))
)
colnames(final_result) <- mapping_list$new_name

# Add pseudocount and filter for targeted metabolites
final_result <- final_result %>% mutate(across(everything(), ~ .x + 1))
target_codes <- mapping_list$new_name[mapping_list$original_name %in% target_metabolites]
available_codes <- intersect(target_codes, colnames(final_result))
if (length(available_codes) == 0) stop("No target metabolites found in metabolomics data.")

# Load clinical data and subset to PR/PD pre-treatment
group <- readRDS(clinical_path) %>%
  filter(PRPDSD %in% c("PR", "PD"), treatment_condition == "Pre")

common_samples <- intersect(rownames(final_result), rownames(group))
if (length(common_samples) == 0) stop("No common samples between metabolites and clinical data.")

meta_sub <- final_result[common_samples, available_codes, drop = FALSE]
group_factor <- factor(group[common_samples, "PRPDSD"], levels = c("PD", "PR"))

# Compute mean per group and log2 fold-change
mean_PR <- colMeans(meta_sub[group_factor == "PR", , drop = FALSE], na.rm = TRUE)
mean_PD <- colMeans(meta_sub[group_factor == "PD", , drop = FALSE], na.rm = TRUE)

code_to_name <- mapping_list[mapping_list$new_name %in% names(mean_PR), ]
log2fc_df <- data.frame(
  code = names(mean_PR),
  PR_mean = mean_PR,
  PD_mean = mean_PD,
  stringsAsFactors = FALSE
) %>%
  left_join(code_to_name, by = c("code" = "new_name")) %>%
  filter(original_name %in% target_metabolites) %>%
  mutate(log2FC = log2(PR_mean / PD_mean),
         Metabolite = factor(original_name, levels = rev(target_metabolites)),
         Comparison = "log2FC (PR/PD)") %>%
  na.omit()

# Draw log2FC column heatmap
p_log2fc <- ggplot(log2fc_df, aes(x = Comparison, y = Metabolite, fill = log2FC)) +
  geom_tile(color = "gray70") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = expression(log[2]~Fold~Change)
  ) +
  scale_y_discrete(expand = c(0, 0), position = "right",
                   labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 12),
    axis.text.y.right = element_text(colour = "black", hjust = 0, size = 12),
    axis.text.y.left = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(colour = "black"),
    legend.title = element_text(colour = "black")
  )

ggsave(out_log2fc, plot = p_log2fc, width = 5, height = 25 * 0.3, dpi = 300,
       limitsize = FALSE)



# =========================================================================== #
# Figure F: Upper-triangle heatmap of Spearman correlations among cluster
# proportions. Rows and columns are ordered by hierarchical clustering.
# Point colour = correlation coefficient (r), size = BH-adjusted p-value.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(WGCNA)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

# ------------------------------ File paths ----------------------------------
prop_csv_path  <- "path/to/prop_sample1_clusters.csv"
output_pdf     <- "path/to/F3/F_cluster_correlation_percentage2.pdf"

# --------------------------- Load proportion matrix -------------------------
prop_matrix <- read.csv(prop_csv_path, row.names = 1, check.names = FALSE)
percent_mat <- as.matrix(prop_matrix)

# ----------- Spearman correlation with WGCNA (fast pairwise) ----------------
cor_res <- WGCNA::corAndPvalue(percent_mat,
                               method = "spearman",
                               use = "pairwise.complete.obs",
                               alternative = "two.sided")
cor_mat <- cor_res$cor
p_mat_raw <- cor_res$p

# Benjamini-Hochberg correction
p_mat_adj <- matrix(p.adjust(as.vector(p_mat_raw), method = "BH"),
                    nrow = nrow(p_mat_raw), ncol = ncol(p_mat_raw))
rownames(p_mat_adj) <- rownames(p_mat_raw)
colnames(p_mat_adj) <- colnames(p_mat_raw)

# Melt to long format
cor_melt <- melt(cor_mat, varnames = c("Cluster1", "Cluster2"))
p_melt   <- melt(p_mat_adj, varnames = c("Cluster1", "Cluster2"))
plot_data <- cbind(cor_melt, P = p_melt$value)

# Set diagonal to 0 (self-correlation)
plot_data$value[plot_data$Cluster1 == plot_data$Cluster2] <- 0

# P-value categories for point size
plot_data <- plot_data %>%
  mutate(P_category = case_when(
    P <= 0.0001 ~ "≤0.0001",
    P <= 0.001  ~ "≤0.001",
    P <= 0.05   ~ "≤0.05",
    TRUE        ~ ">0.05"
  ),
  P_category = factor(P_category,
                      levels = c("≤0.0001", "≤0.001", "≤0.05", ">0.05"),
                      ordered = TRUE))

# ------------------ Hierarchical clustering of clusters ---------------------
# Build full symmetric matrix (clusters 0-49)
full_mat <- matrix(NA, 50, 50,
                   dimnames = list(as.character(0:49), as.character(0:49)))
for (i in seq_len(nrow(plot_data))) {
  full_mat[as.character(plot_data$Cluster1[i]),
           as.character(plot_data$Cluster2[i])] <- plot_data$value[i]
}
full_mat[lower.tri(full_mat)] <- t(full_mat)[lower.tri(full_mat)]
diag(full_mat) <- 0

dist_mat <- as.dist(1 - abs(full_mat))
hc <- hclust(dist_mat, method = "complete")
cluster_order <- hc$order
cluster_labels <- as.character(0:49)[cluster_order]

# Mapping original label -> new position
cluster_mapping <- setNames(seq_along(cluster_labels), cluster_labels)

# --------------- Prepare data in the new clustered order --------------------
all_pairs <- expand.grid(Cluster1 = cluster_labels,
                         Cluster2 = cluster_labels,
                         stringsAsFactors = FALSE)
all_pairs$value <- apply(all_pairs, 1, function(r) full_mat[r[1], r[2]])

# Build lookup for P and P_category (both directions)
lookup <- plot_data %>%
  select(Cluster1, Cluster2, P, P_category)
lookup_both <- bind_rows(
  lookup,
  lookup %>% rename(Cluster1 = Cluster2, Cluster2 = Cluster1)
) %>% distinct(Cluster1, Cluster2, .keep_all = TRUE)

all_pairs <- all_pairs %>%
  left_join(lookup_both, by = c("Cluster1", "Cluster2")) %>%
  mutate(Cluster1_idx = cluster_mapping[Cluster1],
         Cluster2_idx = cluster_mapping[Cluster2])

# Keep only upper triangle (including diagonal)
plot_upper <- all_pairs %>%
  filter(Cluster1_idx <= Cluster2_idx) %>%
  mutate(Cluster1 = factor(Cluster1, levels = cluster_labels),
         Cluster2 = factor(Cluster2, levels = cluster_labels))

# ------------ Grid lines for upper‑triangle cells ---------------------------
n <- length(cluster_labels)
hline_data <- data.frame()
for (y in seq_len(n)) {
  hline_data <- rbind(hline_data,
                      data.frame(y = y - 0.5,
                                 x_start = 0.5,
                                 x_end = y - 0.5))
}
hline_data <- rbind(hline_data, data.frame(y = n + 0.5,
                                           x_start = 0.5,
                                           x_end = n + 0.5))

vline_data <- data.frame()
for (x in seq_len(n)) {
  vline_data <- rbind(vline_data,
                      data.frame(x = x - 0.5,
                                 y_start = x - 0.5,
                                 y_end = n + 0.5))
}

diag_top <- data.frame(
  y = (0:49) + 0.5,
  x_start = (0:49) + 0.5,
  x_end = (0:49) + 1.5
)
diag_right <- data.frame(
  x = (0:49) + 1.5,
  y_start = (0:49) + 0.5,
  y_end = (0:49) + 1.5
)

# ----------------------- Functional category colours ------------------------
color_mapping <- list(
  "#3FA9F5" = c(1, 2, 3, 5, 13, 15, 17, 21, 26, 35),
  "#BDCCD4" = c(4, 7, 19, 24, 29, 39),
  "#FF931E" = c(6, 9, 10, 11, 18, 20, 22, 25, 28, 30, 31, 33, 34, 38, 40, 41, 43, 45, 46),
  "#FF7BAC" = c(8, 14, 23, 37, 42, 44, 49),
  "#00FFFF" = c(16, 27, 32, 36, 47),
  "#7AC943" = c(12),
  "#FCEE21" = c(48),
  "#C69C6D" = c(0)
)

color_df <- data.frame()
for (col in names(color_mapping)) {
  for (num in color_mapping[[col]]) {
    color_df <- rbind(color_df, data.frame(label = as.character(num), color = col))
  }
}
label_colors <- sapply(cluster_labels, function(lab) {
  r <- color_df[color_df$label == lab, "color"]
  if (length(r)) as.character(r) else "black"
})

# ----------------------------- Draw heatmap ---------------------------------
p <- ggplot(plot_upper, aes(x = Cluster1, y = Cluster2)) +
  # Upper-triangle grid lines
  geom_segment(data = hline_data,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               color = "gray60", linewidth = 0.3, alpha = 0.7) +
  geom_segment(data = vline_data,
               aes(x = x, xend = x, y = y_start, yend = y_end),
               color = "gray60", linewidth = 0.3, alpha = 0.7) +
  geom_segment(data = diag_top,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               color = "gray60", linewidth = 0.3, alpha = 0.7) +
  geom_segment(data = diag_right,
               aes(x = x, xend = x, y = y_start, yend = y_end),
               color = "gray60", linewidth = 0.3, alpha = 0.7) +
  # Points
  geom_point(aes(size = P_category, fill = value),
             color = "black", alpha = 0.8, shape = 21, stroke = 0.7) +
  # Fill scale: RdYlBu reversed
  scale_fill_gradientn(
    colors = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
    limits = c(-0.4, 0.6),
    name = "Correlation coefficient (r)",
    breaks = c(-0.3, 0, 0.3, 0.6),
    guide = guide_colorbar(barwidth = 0.8, barheight = 10,
                           title.position = "top"),
    oob = squish
  ) +
  # Size scale for adjusted p-value
  scale_size_manual(
    name = "Adjusted P value",
    values = c("≤0.0001" = 4, "≤0.001" = 3.5, "≤0.05" = 3, ">0.05" = 2.5),
    drop = FALSE,
    guide = guide_legend(override.aes = list(shape = 21, color = "gray60",
                                             fill = "white", stroke = 0.7),
                         title.position = "top")
  ) +
  # Axes: top x, left y, with functional colours
  scale_x_discrete(breaks = cluster_labels, drop = FALSE,
                   expand = c(0, 0), position = "top") +
  scale_y_discrete(breaks = cluster_labels, drop = FALSE,
                   expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.x.top = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 10, color = label_colors[cluster_labels]),
    axis.text.x.bottom = element_blank(),
    axis.text.y = element_text(hjust = 1, size = 10,
                               color = label_colors[cluster_labels]),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20),
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    text = element_text(family = "sans")
  ) +
  coord_fixed() +
  ggtitle("")

ggsave(output_pdf, p, width = 10, height = 10, device = cairo_pdf)