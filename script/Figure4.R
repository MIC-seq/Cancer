# =========================================================================== #
# Figure A: NMF-based subtyping of baseline samples using species-level 
# relative abundance. Consensus matrix heatmap shows three distinct clusters.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(readxl)
library(NMF)

# ------------------------------ File paths ----------------------------------
input_xlsx  <- "path/to/A_prop_sample1_clusters.xlsx"
output_dir  <- "path/to/F4"
rank_pdf1   <- file.path(output_dir, "A_1_rank_estimation_plot.pdf")
rank_pdf2   <- file.path(output_dir, "A_2_consensus_map_rank_estimation.pdf")
final_pdf   <- file.path(output_dir, "A_3_final_consensus_heatmap.pdf")
csv_cons    <- file.path(output_dir, "A_4_consensus_matrix.csv")
csv_group   <- file.path(output_dir, "A_5_sample_cluster_assignment.csv")
csv_ordered <- file.path(output_dir, "A_6_ordered_feature_matrix.csv")

# --------------------------- 1. Load abundance matrix -----------------------
species_cluster <- read_excel(input_xlsx, sheet = 1)
species_cluster <- as.data.frame(species_cluster)
rownames(species_cluster) <- species_cluster[, 1]
species_cluster <- species_cluster[, -1, drop = FALSE]

# -------------------- 2. Define clinical response groups --------------------
baseline <- c("A1733","A1516","A1532","A1536","A1562","A1555","A1541","A1641",
              "A1746","A1764","A1821","A1823","A2067","A1593","A1741","A1697",
              "A1822","A1869","A1550","A1582","A1668","A1575","A1710","A1756",
              "A2039","A1787","A2038","A1685")

PR_B    <- c("A1733","A1536","A1562","A1593","A1697","A1741","A1822","A1869","A1787","A2038")
PD_B    <- c("A1516","A1541","A1555","A1641")
SD_PD_B <- c("A1532","A1764","A1746","A1823")
PR_T1   <- c("A1955","A1686","A2199","A1709","A1843","A2162","A2148")
PR_T2   <- c("A1955","A1997","A2307","A2291","A2163")
PD_T1   <- c("A1747","A1748","A1763","A1607")
SD_PD_T1 <- c("A2291","A1645","A2044","A2037","A2161")

# Build sample annotation data frame
clinical_response <- data.frame(sample = rownames(species_cluster),
                                response = "NA")
clinical_response$response[clinical_response$sample %in% PR_B]    <- "PR_B"
clinical_response$response[clinical_response$sample %in% PD_B]    <- "PD_B"
clinical_response$response[clinical_response$sample %in% SD_PD_B] <- "SD_PD_B"
clinical_response$response[clinical_response$sample %in% PR_T1]   <- "PR_T1"
clinical_response$response[clinical_response$sample %in% PR_T2]   <- "PR_T2"
clinical_response$response[clinical_response$sample %in% PD_T1]   <- "PD_T1"
clinical_response$response[clinical_response$sample %in% SD_PD_T1] <- "SD_PD_T1"

clinical_response_new <- data.frame(response = clinical_response$response)
rownames(clinical_response_new) <- clinical_response$sample

# ---------------- 3. Transpose baseline matrix for NMF ---------------------
species_cluster_t <- t(species_cluster[baseline, ])
clinical_response_new1 <- clinical_response_new[colnames(species_cluster_t), , drop = FALSE]
colnames(clinical_response_new1) <- "response"

# -------------------- 4. Estimate optimal rank (2-10) -----------------------
ranks <- 2:10
estim <- nmf(species_cluster_t, ranks, nrun = 50)

pdf(rank_pdf1, width = 8, height = 6)
plot(estim)
dev.off()

pdf(rank_pdf2, width = 10, height = 8)
consensusmap(estim)
dev.off()

# ------------------ 5. Final NMF with rank = 3 ------------------------------
seed <- 2025820
nmf_fit <- nmf(species_cluster_t, rank = 3, nrun = 50,
               seed = seed, method = "brunet")

# Extract signature features
index <- extractFeatures(nmf_fit, "max")
sig.order <- unlist(index)
NMF_exp <- species_cluster_t[sig.order, ]
NMF_exp <- na.omit(NMF_exp)

# Subtype assignment
group <- predict(nmf_fit)

# Colour definitions
jco <- c("#2874C5", "#EABF00", "#C6524A", "#868686")
ann_colors <- list(
  cluster  = c("1" = jco[1], "2" = jco[2], "3" = jco[3], "4" = jco[4]),
  response = c(
    "NA"        = "yellow",
    "PR_B"      = "green",
    "PR_T1"     = "yellowgreen",
    "PR_T2"     = "darkgreen",
    "PD_B"      = "red",
    "PD_T1"     = "darkred",
    "SD_PD_B"   = "lightblue",
    "SD_PD_T1"  = "darkblue"
  )
)

# --------------- 6. Final consensus heatmap ---------------------------------
pdf(final_pdf, width = 12, height = 10)
consensusmap(nmf_fit,
             annRow = clinical_response_new1,
             annCol = data.frame("cluster" = group[colnames(NMF_exp)]),
             annColors = ann_colors)
dev.off()

# ---------------------- 7. Export result tables -----------------------------
consensus_mat <- consensus(nmf_fit)
write.csv(consensus_mat, csv_cons, quote = FALSE)

group_df <- data.frame(Sample = names(group), Cluster = group)
write.csv(group_df, csv_group, row.names = FALSE, quote = FALSE)

write.csv(NMF_exp, csv_ordered, quote = FALSE)



# =========================================================================== #
# Figure B: Heatmap of z-score-normalized cluster proportions across selected
# baseline samples. Clusters are ordered by hierarchical clustering (Euclidean
# distance, Ward.D2).
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(dplyr)
library(tidyr)
library(pheatmap)

# ------------------------------ File paths ----------------------------------
metadata_path <- "path/to/metadata.rds"
output_dir    <- "path/to/F4"
output_pdf    <- file.path(output_dir, "B.pdf")

# -------------------------- Load and prepare data ---------------------------
sc_df <- readRDS(metadata_path)

selected_samples <- c("A2067", "A1575", "A1550", "A1685", "A1746", "A1823",
                      "A1787", "A1536", "A1710", "A1593", "A1741", "A1582",
                      "A1641", "A1869", "A2039", "A1541", "A1516", "A1555",
                      "A2038", "A1562", "A1697", "A1764", "A1532", "A1668",
                      "A1822", "A1756", "A1733", "A1821")

sc_selected <- sc_df %>%
  filter(orig.ident %in% selected_samples) %>%
  mutate(orig.ident = factor(orig.ident, levels = selected_samples))

# Compute cluster proportions per sample
cluster_prop <- sc_selected %>%
  group_by(orig.ident) %>%
  mutate(total_cells = n()) %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(proportion = n() / first(total_cells), .groups = "drop") %>%
  mutate(seurat_clusters = as.character(seurat_clusters))

# Pivot to wide matrix (clusters x samples)
prop_wide <- cluster_prop %>%
  pivot_wider(names_from = orig.ident,
              values_from = proportion,
              values_fill = 0)

prop_matrix <- as.matrix(prop_wide %>% select(-seurat_clusters))
rownames(prop_matrix) <- prop_wide$seurat_clusters

# Z-score normalise across rows (per cluster)
prop_scaled <- t(scale(t(prop_matrix)))
prop_scaled[is.na(prop_scaled)] <- 0

# Hierarchical clustering of clusters
hclust_result <- hclust(dist(prop_scaled, method = "euclidean"),
                        method = "ward.D2")

# Colour palette (blue – white – orange – red)
my_colors <- colorRampPalette(c("#4575B4", "#FFFFFF", "#FDAE61", "#D73027"))(100)
legend_breaks <- seq(-1, 2, by = 0.5)

# ------------------------------ Draw heatmap ---------------------------------
p <- pheatmap(
  prop_scaled,
  color = my_colors,
  breaks = seq(-1, 2, length.out = 100),
  cluster_rows = hclust_result,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  scale = "none",
  angle_col = 45,
  fontsize_row = 8,
  fontsize_col = 8,
  treeheight_row = 50,
  legend_breaks = legend_breaks,
  legend_labels = legend_breaks,
  main = ""
)

# Save as PDF
pdf(output_pdf, width = 10, height = 8)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()




# =========================================================================== #
# Figure C: Lollipop charts of Cohen's d comparing cluster proportions in each
# NMF subtype versus all other subtypes, using pre-treatment samples.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(dplyr)
library(ggplot2)
library(effsize)
library(ggrepel)
library(scales)

# ------------------------------ File paths ----------------------------------
sc_df_path   <- "path/to/01metadata.rds"
group_path   <- "path/to/Clinical_information.rds"
output_dir   <- "path/to/F4"
out_pdf1     <- file.path(output_dir, "C1.pdf")
out_pdf2     <- file.path(output_dir, "C2.pdf")
out_pdf3     <- file.path(output_dir, "C3.pdf")

# -------------------------- Load and prepare data ---------------------------
sc_df <- readRDS(sc_df_path)

group <- readRDS(group_path) %>%
  filter(treatment_condition == "Pre") %>%
  rename(model = subtype) %>%
  filter(model != "Unknow")

# Keep only cells from common samples
common_samples <- intersect(unique(sc_df$orig.ident), group$sample)
sc_df_filtered <- sc_df[sc_df$orig.ident %in% common_samples, ]

# Merge subtype information
sc_df_merged <- left_join(sc_df_filtered, group,
                          by = c("orig.ident" = "sample"))

# Compute cluster proportions per sample
prop_df <- sc_df_merged %>%
  group_by(orig.ident, model, seurat_clusters) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(total_cells = sum(n_cells),
         prop = n_cells / total_cells) %>%
  ungroup()

# --------------- Function: compute Cohen's d for one subtype -----------------
compute_cohen_d <- function(prop_df, target_model) {
  prop_target <- prop_df %>%
    filter(model == target_model) %>%
    select(seurat_clusters, orig.ident, prop) %>%
    rename(target_prop = prop)

  prop_control <- prop_df %>%
    filter(model != target_model) %>%
    select(seurat_clusters, orig.ident, prop) %>%
    rename(control_prop = prop)

  clusters <- unique(prop_df$seurat_clusters)
  cluster_d_values <- data.frame()

  for (cl in clusters) {
    target  <- prop_target %>% filter(seurat_clusters == cl) %>% pull(target_prop)
    control <- prop_control %>% filter(seurat_clusters == cl) %>% pull(control_prop)

    mean_target  <- if (length(target)  > 0) mean(target,  na.rm = TRUE) else 0
    mean_control <- if (length(control) > 0) mean(control, na.rm = TRUE) else 0

    n_target  <- length(target)
    n_control <- length(control)

    d <- NA; ci_lower <- NA; ci_upper <- NA
    significant <- "No"
    direction   <- "Non-positive"
    cluster_label <- NA

    if (n_target > 1 && n_control > 1) {
      suppressMessages(
        tryCatch({
          cd <- cohen.d(target, control)
          d        <- cd$estimate
          ci_lower <- cd$conf.int[1]
          ci_upper <- cd$conf.int[2]
          if (!is.na(ci_lower) && !is.na(ci_upper))
            significant <- ifelse(ci_lower > 0 | ci_upper < 0, "Yes", "No")
          if (!is.na(d))
            direction <- ifelse(d > 0, "Positive", "Non-positive")
          if (!is.na(d) && d >= 0.5)
            cluster_label <- paste("Cluster", cl)
        }, error = function(e) NULL)
      )
    }

    cluster_d_values <- rbind(cluster_d_values, data.frame(
      Cluster         = cl,
      Mean_Target     = mean_target,
      Mean_Control    = mean_control,
      Cohens_d        = d,
      CI_lower        = ci_lower,
      CI_upper        = ci_upper,
      N_cluster       = n_target,
      N_other         = n_control,
      Significant     = significant,
      Effect_Direction = direction,
      Cluster_Label   = cluster_label,
      stringsAsFactors = FALSE
    ))
  }
  return(cluster_d_values)
}

# ---------- Function: draw lollipop chart for one subtype -------------------
draw_lollipop <- function(cluster_d_values, subtype_title, positive_color) {
  ggplot(cluster_d_values, aes(y = Cohens_d, x = Cluster, color = Effect_Direction)) +
    geom_segment(aes(y = 0, yend = Cohens_d, x = Cluster, xend = Cluster),
                 linewidth = 0.8, alpha = 0.7) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
    geom_text_repel(aes(y = Cohens_d, label = Cluster_Label),
                    color = "black", size = 3.5,
                    direction = "y", nudge_y = 0.1,
                    segment.size = 0.2, max.overlaps = Inf,
                    min.segment.length = 0, box.padding = 0.5) +
    scale_color_manual(values = c("Positive" = positive_color,
                                  "Non-positive" = "gray")) +
    labs(y = "Cohen's d", x = NULL, title = subtype_title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_y_continuous(breaks = function(y) {
      brk <- breaks_extended()(y)
      unique(sort(c(brk, 0.5)))
    }, expand = expansion(mult = c(0.05, 0.15)))
}

# -------------------- Compute and plot for each subtype ---------------------
model_list <- unique(prop_df$model)

for (m in model_list) {
  result <- compute_cohen_d(prop_df, m)

  subtype_title <- paste("Subtype", m)
  pos_color <- switch(m,
                      "1" = "yellowgreen",
                      "2" = "red",
                      "3" = "royalblue",
                      "gray50")

  p <- draw_lollipop(result, subtype_title, pos_color)

  out_file <- switch(m,
                     "1" = out_pdf1,
                     "2" = out_pdf2,
                     "3" = out_pdf3,
                     file.path(output_dir, paste0("C", m, ".pdf")))
  ggsave(out_file, plot = p, width = 15, height = 5, dpi = 300, bg = "white")
}



# =========================================================================== #
# Figure D: Alluvial plot showing the association between NMF subtypes and
# PRPDSD response categories in pre-treatment samples. Flows are colored by
# standardized residuals from a chi‑square test; asterisks mark significant
# residuals (*** p<0.001, ** p<0.01, * p<0.05).
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(dplyr)
library(ggplot2)
library(ggalluvial)

# ------------------------------ File paths ----------------------------------
metadata_path <- "path/to/metadata.rds"
output_dir    <- "path/to/F4"
output_pdf    <- file.path(output_dir, "D_subtype_PRPDSD_association_plot.pdf")

# -------------------------- Load and prepare data ---------------------------
sc_df <- readRDS(metadata_path)

sc_df <- sc_df %>%
  filter(treatment_condition == "Pre") %>%
  filter(subtype != "Unknow") %>%
  filter(PRPDSD != "Unknow") %>%
  mutate(PRPDSD = factor(PRPDSD, levels = c("PR", "SD", "PD")))

sc_df <- sc_df %>%
  select(orig.ident, subtype, PRPDSD) %>%
  drop_na() %>%
  distinct() %>%
  mutate(subtype = paste0("Subtype ", subtype))

# ------------------- Chi-square test and residuals -------------------------
contingency_table <- table(sc_df$subtype, sc_df$PRPDSD)
chi_test <- chisq.test(contingency_table)
chi_residuals <- chi_test$stdres

# --------------------- Prepare alluvial plot data ---------------------------
freq_table <- as.data.frame(as.table(contingency_table))
colnames(freq_table) <- c("subtype", "PRPDSD", "count")

sankey_data <- as.data.frame(as.table(chi_residuals))
colnames(sankey_data) <- c("subtype", "PRPDSD", "residual")

sankey_data <- sankey_data %>%
  left_join(freq_table, by = c("subtype", "PRPDSD")) %>%
  mutate(
    sig_stars = case_when(
      abs(residual) > 3.29 ~ "***",
      abs(residual) > 2.58 ~ "**",
      abs(residual) > 1.96 ~ "*",
      TRUE ~ ""
    ),
    flow = count,
    direction = ifelse(residual > 0, "Positive", "Negative")
  )

# Compute vertical positions for significance labels
sankey_data <- sankey_data %>%
  group_by(subtype) %>%
  arrange(PRPDSD) %>%
  mutate(
    cumulative_flow = cumsum(flow),
    start_position = cumulative_flow - flow,
    mid_position = start_position + flow / 2
  ) %>%
  ungroup()

# ----------------------- Custom colour gradient -----------------------------
custom_colors <- c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7",
                   "#fddbc7", "#ef8a62", "#b2182b")

# --------------------------- Draw alluvial plot -----------------------------
p <- ggplot(sankey_data,
            aes(axis1 = subtype, axis2 = PRPDSD, y = flow)) +
  geom_alluvium(
    aes(fill = residual),
    width = 1/8,
    alpha = 0.8,
    curve_type = "sigmoid"
  ) +
  geom_stratum(width = 1/8, fill = "grey95", color = "grey60") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  geom_text(data = sankey_data %>% filter(sig_stars != ""),
            aes(x = 1.2, y = mid_position, label = sig_stars),
            size = 5, fontface = "bold", color = "black",
            hjust = 1.2, check_overlap = TRUE) +
  scale_x_discrete(limits = c("subtype", "PRPDSD"),
                   expand = c(0.05, 0.05)) +
  scale_fill_gradientn(
    name = "Standardized\nResidual",
    colors = custom_colors,
    limits = c(-1, 1),
    oob = scales::oob_squish,
    breaks = seq(-1, 1, by = 0.5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    title = "",
    subtitle = paste0("Chi-square test: χ² = ", round(chi_test$statistic, 2),
                      ", df = ", chi_test$parameter,
                      ", p = ", format.pval(chi_test$p.value, digits = 3)),
    y = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  )

ggsave(output_pdf, plot = p, width = 8, height = 8, device = "pdf")



# =========================================================================== #
# Figure E: Kaplan-Meier curves for progression-free survival (PFS) stratified
# by NMF subtypes in pre-treatment samples.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(survival)
library(survminer)
library(dplyr)

# ------------------------------ File paths ----------------------------------
clinical_path <- "path/to/Clinical_information.rds"
output_dir    <- "path/to/F4"
output_pdf    <- file.path(output_dir, "E_PFS.pdf")

# -------------------------- Load and prepare data ---------------------------
group <- readRDS(clinical_path) %>%
  filter(treatment_condition == "Pre") %>%
  filter(PFS_situation != "Unknow") %>%
  filter(subtype != "Unknow") %>%
  filter(time != "Unknow") %>%
  mutate(
    time = as.numeric(time),
    PFS_situation = as.factor(PFS_situation),
    cluster = as.factor(subtype)
  )

# Define event: PFS_situation == 1 as event, 0 as censored
group$event <- ifelse(group$PFS_situation == 1, 1, 0)

# Fit Kaplan-Meier curves
fit_simple <- survfit(Surv(time, event) ~ cluster, data = group)

# Plot
pic <- ggsurvplot(
  fit_simple,
  surv.median.line = "hv",
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Time (days)",
  ylab = "PFS Probability",
  legend.title = "",
  legend.labs = c("Subtype 1", "Subtype 2", "Subtype 3"),
  ggtheme = theme_minimal(),
  break.x.by = 100,
  palette = c("#bc5148", "#3090a1", "#4daf4a"),
  title = "Kaplan-Meier Curve for PFS"
)

# Save as PDF
pdf(output_pdf, width = 8, height = 8, bg = "transparent")
print(pic$plot)
dev.off()




# =========================================================================== #
# Figure F: NMF-based functional module decomposition of microbial clusters.
# Consensus heatmaps, basis weight plots, and module networks are shown.
# Figure G: Radar charts of functional module (CM) distribution across
# NMF subtypes in pre-treatment samples.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(CoVarNet)
library(dplyr)
library(stringr)
library(igraph)
library(RColorBrewer)
library(circlize)
library(ggsci)
library(reshape2)
library(ggplot2)
library(tidyr)

# ------------------------------ File paths ----------------------------------
clinical_path   <- "path/to/Clinical_information.rds"
meta_path       <- "path/to/01metadata.rds"
prop_csv_path   <- "path/to/prop_sample1_clusters.csv"
output_dir      <- "path/to/F4"
f_kvalue_pdf    <- file.path(output_dir, "F_Kvalue.pdf")
f_heatmap_pdf   <- file.path(output_dir, "F_01heatmap.pdf")
f_column_pdf    <- file.path(output_dir, "F_02column.pdf")
f_network_pdf   <- file.path(output_dir, "F_04net.pdf")
g_dist_pdf      <- file.path(output_dir, "G_03model.pdf")
g_radar_pdf     <- file.path(output_dir, "G_03model_radar_filled.pdf")

# -------------------------- Load clinical data ------------------------------
group <- readRDS(clinical_path) %>%
  filter(treatment_condition == "Pre")

# =========================================================================== #
#                          Figure F: NMF analysis
# =========================================================================== #

# --------------------- Load and prepare cluster matrix ----------------------
cluster_sample <- read.csv(prop_csv_path, row.names = 1)
colnames(cluster_sample) <- gsub("^X", "cluster", colnames(cluster_sample))
cluster_sample <- t(cluster_sample)

common_samples <- intersect(colnames(cluster_sample), group$sample)
cluster_sample <- cluster_sample[, common_samples, drop = FALSE]
mat_fq_raw <- as.data.frame(cluster_sample)

# Annotate major functional categories
mat_fq_raw <- mat_fq_raw %>%
  mutate(majorCluster = case_when(
    str_extract(rownames(.), "\\d+") %in% c("1","2","3","5","13","15","17","21","26","35") ~
      "Core Energy Metabolism & Central Carbon Processing",
    str_extract(rownames(.), "\\d+") %in% c("4","7","19","24","29","39") ~
      "Stress Response & Cellular Maintenance",
    str_extract(rownames(.), "\\d+") %in% c("9","10","11","22","31","33","6","18","25","30","34","38","45","20","28","40","41","43","46") ~
      "Substrate Utilization & Specialized Fermentation",
    str_extract(rownames(.), "\\d+") %in% c("8","14","23","37","42","44","49") ~
      "Beneficial Metabolite Production (Host-Microbe Interactions)",
    str_extract(rownames(.), "\\d+") %in% c("16","27","32","36","47") ~
      "Cell Growth, Division & Motility",
    str_extract(rownames(.), "\\d+") %in% c("12") ~ "Defense & Resistance",
    str_extract(rownames(.), "\\d+") %in% c("48") ~ "Cofactor & Micronutrient Procurement",
    str_extract(rownames(.), "\\d+") %in% c("0")  ~ "Unspecified & Generalist Functions",
    TRUE ~ "Unknown"
  )) %>%
  relocate(majorCluster, .before = 1)

# Normalize (min-max)
mat_fq_norm <- freq_normalize(mat_fq_raw, normalize = "minmax")

# ----------- Rank estimation (2–20) ----------------------------------------
res <- nmf(mat_fq_norm, rank = 2:20, method = "nsNMF",
           seed = rep(123456, 6), .options = "vp")
pdf(f_kvalue_pdf, width = 4000/300, height = 2000/300, bg = "transparent")
plot(res)
dev.off()

# ----------- Final NMF (K = 4) ---------------------------------------------
K <- 4
NMF_K <- nmf(mat_fq_norm, K, method = "nsNMF", seed = rep(77, 6), nrun = 30)

# Reorder modules by total contribution
module_importance <- rowSums(coef(NMF_K))
sorted_order <- order(module_importance, decreasing = TRUE)
if (!is.null(sorted_order)) {
  n <- length(sorted_order)
  name_base <- if (n <= 9) paste0("CM0", 1:n) else paste0("CM", sprintf("%02d", 1:n))
  reordered_names <- name_base[sorted_order]
  colnames(basis(NMF_K)) <- reordered_names
  rownames(coef(NMF_K)) <- reordered_names
}

# Weight heatmap and column plot
pdf(f_heatmap_pdf, width = 3000/300, height = 2000/300, bg = "transparent")
gr.weight_all(NMF_K)
dev.off()

pdf(f_column_pdf, width = 3000/300, height = 2000/300, bg = "transparent")
gr.weight_top(NMF_K)
dev.off()

# ----------------------- Prepare meta for networks --------------------------
meta <- readRDS(meta_path)
meta <- meta[meta$orig.ident %in% group$sample, ]
meta <- meta %>%
  left_join(group %>% select(sample, subtype), by = c("orig.ident" = "sample")) %>%
  mutate(
    majorCluster = case_when(
      seurat_clusters %in% c("1","2","3","5","13","15","17","21","26","35") ~ "Core Energy Metabolism & Central Carbon Processing",
      seurat_clusters %in% c("4","7","19","24","29","39") ~ "Stress Response & Cellular Maintenance",
      seurat_clusters %in% c("9","10","11","22","31","33","6","18","25","30","34","38","45","20","28","40","41","43","46") ~ "Substrate Utilization & Specialized Fermentation",
      seurat_clusters %in% c("8","14","23","37","42","44","49") ~ "Beneficial Metabolite Production (Host-Microbe Interactions)",
      seurat_clusters %in% c("16","27","32","36","47") ~ "Cell Growth, Division & Motility",
      seurat_clusters %in% c("12") ~ "Defense & Resistance",
      seurat_clusters %in% c("48") ~ "Cofactor & Micronutrient Procurement",
      seurat_clusters %in% c("0")  ~ "Unspecified & Generalist Functions",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(seurat_clusters, orig.ident, majorCluster, subtype) %>%
  rename(sampleID = orig.ident, subCluster = seurat_clusters,
         majorCluster = majorCluster, subtype = subtype)

meta <- meta %>%
  rownames_to_column("cell_id") %>%
  add_count(majorCluster, subCluster) %>%
  group_by(subCluster) %>%
  filter(n == max(n)) %>%
  select(-n) %>%
  ungroup() %>%
  column_to_rownames("cell_id") %>%
  mutate(subCluster = paste0("cluster", subCluster))

# --------------------- Pairwise correlation and network ---------------------
cor_pair <- pair_correlation(mat_fq_raw, method = "pearson")

# Custom network function (silent, no cat/print)
sm_cm_network <- function(NMFres, Corres, ...) {
  var_args <- list(...)
  top_n <- if (!is.null(var_args[["top_n"]])) var_args[["top_n"]] else 10
  corr  <- if (!is.null(var_args[["corr"]]))  var_args[["corr"]]  else 0.2
  fdr   <- if (!is.null(var_args[["pval_fdr"]])) var_args[["pval_fdr"]] else 0.05

  apply_top_n       <- if (!is.null(var_args[["apply_top_n"]]))       var_args[["apply_top_n"]]       else TRUE
  apply_corr        <- if (!is.null(var_args[["apply_corr"]]))        var_args[["apply_corr"]]        else TRUE
  apply_fdr         <- if (!is.null(var_args[["apply_fdr"]]))         var_args[["apply_fdr"]]         else TRUE
  apply_spe         <- if (!is.null(var_args[["apply_spe"]]))         var_args[["apply_spe"]]         else TRUE
  remove_isolated   <- if (!is.null(var_args[["remove_isolated"]]))   var_args[["remove_isolated"]]   else TRUE

  # Build annotation
  listA <- Corres[, c("subCluster1", "majorCluster1")]
  listB <- Corres[, c("subCluster2", "majorCluster2")]
  colnames(listA) <- c("subCluster", "majorCluster")
  colnames(listB) <- c("subCluster", "majorCluster")
  ann <- rbind(listA, listB)
  ann <- ann[!duplicated(ann$subCluster), ]
  rownames(ann) <- ann$subCluster

  w <- basis(NMFres)
  if (is.null(colnames(w))) {
    colnames(w) <- sprintf("CM%02d", 1:ncol(w))
    rownames(coef(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
    colnames(basis(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
  }
  w <- w[, order(colnames(w))]

  weight <- reshape2::melt(w)
  colnames(weight) <- c("subCluster", "cm", "weight")
  weight <- weight %>% group_by(cm) %>% arrange(cm, desc(weight)) %>% ungroup()

  if (apply_top_n) {
    weight_filtered <- weight %>% group_by(cm) %>% arrange(desc(weight)) %>%
      slice(1:top_n) %>% ungroup()
  } else {
    weight_filtered <- weight
  }

  meta_sc <- data.frame(subCluster = c(Corres$subCluster1, Corres$subCluster2),
                        majorCluster = c(Corres$majorCluster1, Corres$majorCluster2))
  meta_sc <- meta_sc[!duplicated(meta_sc$subCluster), ]
  rownames(meta_sc) <- meta_sc$subCluster
  meta_sc <- meta_sc[order(meta_sc$subCluster), ]

  pl_df <- Corres
  if (apply_corr) pl_df <- pl_df[pl_df$correlation > corr, ]
  if (apply_fdr)  pl_df <- pl_df[pl_df$pval_fdr <= fdr, ]
  if (apply_spe) {
    n_all <- nrow(meta_sc)
    spe_cutoff <- 1 - ((top_n - 1) * 2 - 1) / ((n_all - 1) * 2 - 1)
    pl_df <- pl_df[pl_df$spe >= spe_cutoff, ]
  }

  graph <- graph_from_data_frame(pl_df, directed = FALSE)
  node_global <- data.frame(subCluster = V(graph)$name,
                            majorCluster = meta_sc$majorCluster[match(V(graph)$name, meta_sc$subCluster)])
  node_global <- node_global[order(node_global$subCluster), ]

  node_each <- data.frame()
  edge_each <- data.frame()

  for (cm in unique(weight_filtered$cm)) {
    cm_nodes <- weight_filtered$subCluster[weight_filtered$cm == cm]
    sub_graph <- subgraph(graph, vids = intersect(V(graph)$name, cm_nodes))
    if (remove_isolated && length(V(sub_graph)) > 0) {
      sub_graph <- delete.vertices(sub_graph, V(sub_graph)[degree(sub_graph) == 0])
    }
    if (length(V(sub_graph)) != 0) {
      tmp1 <- data.frame(cm = cm, subCluster = V(sub_graph)$name)
      node_each <- rbind(node_each, tmp1)
      tmp2 <- pl_df[(pl_df$subCluster1 %in% tmp1$subCluster) &
                    (pl_df$subCluster2 %in% tmp1$subCluster), ]
      edge_each <- rbind(edge_each, data.frame(cm = cm, tmp2))
    }
  }

  if (nrow(node_each) > 0) {
    node_each$majorCluster <- meta_sc$majorCluster[match(node_each$subCluster, meta_sc$subCluster)]
  }
  if (nrow(node_each) > 0) {
    weight_final <- merge(node_each, weight_filtered)
    weight_final <- weight_final[order(weight_final$cm, -weight_final$weight), ]
    rownames(weight_final) <- 1:nrow(weight_final)
  } else {
    weight_final <- data.frame()
  }

  list(global = list(node = node_global, edge = pl_df),
       each = list(node = node_each, edge = edge_each),
       raw = NMFres, filter = weight_final, ann = ann,
       filter_params = list(top_n = top_n, corr = corr, fdr = fdr,
                            apply_top_n = apply_top_n, apply_corr = apply_corr,
                            apply_fdr = apply_fdr, apply_spe = apply_spe,
                            remove_isolated = remove_isolated))
}

topnum <- floor(nrow(basis(NMF_K)) / ncol(basis(NMF_K)))
network <- sm_cm_network(
  NMF_K, cor_pair,
  corr = 0.2, top_n = topnum,
  apply_fdr = FALSE, apply_corr = TRUE,
  apply_top_n = TRUE, apply_spe = TRUE, remove_isolated = TRUE
)

# Filter network nodes by weight threshold > 0.05 (multi-module assignment)
basis_matrix <- basis(NMF_K)
weight_matrix <- as.data.frame(basis_matrix[, c("CM01","CM02","CM03","CM04")])
weight_threshold <- 0.05

important_clusters <- lapply(c("CM01","CM02","CM03","CM04"), function(mod) {
  rownames(weight_matrix)[weight_matrix[[mod]] > weight_threshold]
})
names(important_clusters) <- c("CM01","CM02","CM03","CM04")

filter_each_by_assignment <- function(each_obj, important_clusters) {
  cluster_to_module <- setNames(
    rep(names(important_clusters), lengths(important_clusters)),
    unlist(important_clusters)
  )
  # Filter nodes
  filtered_nodes <- each_obj$node[
    mapply(function(cl, cm) cl %in% important_clusters[[cm]],
           each_obj$node$subCluster, each_obj$node$cm), ]
  # Filter edges
  if (nrow(each_obj$edge) > 0) {
    keep <- mapply(function(c1, c2, cm) {
      c1 %in% important_clusters[[cm]] && c2 %in% important_clusters[[cm]]
    }, each_obj$edge$subCluster1, each_obj$edge$subCluster2, each_obj$edge$cm)
    filtered_edges <- each_obj$edge[keep, ]
  } else {
    filtered_edges <- each_obj$edge
  }
  list(node = filtered_nodes, edge = filtered_edges)
}

each <- filter_each_by_assignment(network$each, important_clusters)

# Custom network plot with legend
gr.igraph_each_with_legend <- function(each, ...) {
  var_args <- list(...)
  Layout <- if (!is.null(var_args[["Layout"]])) var_args[["Layout"]] else layout_in_circle
  node <- each$node; edge <- each$edge
  edge <- edge[, c("subCluster1","subCluster2","cm","correlation",
                   "pval","pval_fdr","spe","majorCluster1","majorCluster2")]
  graph <- graph_from_data_frame(edge, directed = FALSE)
  meta_sc <- data.frame(subCluster = c(edge$subCluster1, edge$subCluster2),
                        majorCluster = c(edge$majorCluster1, edge$majorCluster2))
  meta_sc <- meta_sc[!duplicated(meta_sc$subCluster), ]
  rownames(meta_sc) <- meta_sc$subCluster
  V(graph)$majorCluster <- meta_sc[V(graph)$name, "majorCluster"]

  colors <- RColorBrewer::brewer.pal(length(unique(V(graph)$majorCluster)), "Set3")
  names(colors) <- unique(V(graph)$majorCluster)
  V(graph)$frame.color <- V(graph)$color <- colors[V(graph)$majorCluster]

  col_fun <- circlize::colorRamp2(
    breaks = quantile(E(graph)$spe, probs = seq(0, 1, length = 6)),
    colors = (ggsci::pal_material("grey", n = 10))(10)[3:8]
  )
  E(graph)$color <- col_fun(E(graph)$spe)
  E(graph)$width <- 1

  par(mfrow = c(2, 6), mar = c(0, 0, 0, 0) + 0.5)
  for (cm in unique(node$cm)) {
    sub_graph <- subgraph(graph, vids = intersect(V(graph)$name,
                          node$subCluster[node$cm == cm]))
    sub_graph <- delete.vertices(sub_graph, V(sub_graph)[degree(sub_graph) == 0])
    plot.igraph(sub_graph, layout = Layout, xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
                vertex.size = 50, vertex.label.cex = 5/8,
                vertex.label.color = "black", edge.curved = FALSE)
    title(cm, cex.main = 7/8, line = -0.5)
  }

  plot.new()
  par(fig = c(0, 1, 0, 1), new = TRUE)
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  legend("center", legend = names(colors), pch = 21, pt.bg = colors, col = colors,
         pt.cex = 1.5, cex = 0.7, y.intersp = 1.5, bty = "n", title = "",
         ncol = min(3, length(colors)))

  par(fig = c(0, 1, 0, 0.33), new = TRUE)
  spe_values <- quantile(E(graph)$spe, probs = seq(0, 1, length = 5))
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  legend("center", legend = round(spe_values, 2), lwd = 3,
         col = col_fun(spe_values), cex = 0.7, bty = "n",
         title = "Edge (spe)", horiz = TRUE)
}

pdf(f_network_pdf, width = 3000/300, height = 1000/300, bg = "transparent")
set.seed(1234)
gr.igraph_each_with_legend(each, Layout = layout_in_circle)
dev.off()

# =========================================================================== #
#                          Figure G: Radar plots
# =========================================================================== #

# ------------------ CM distribution by subtype (bar plot) -------------------
pdf(g_dist_pdf, width = 3000/300, height = 2000/300, bg = "transparent")
set.seed(123)
gr.distribution(NMF_K, meta = meta, group = "subtype")
dev.off()

# ------------------------- Radar chart (mean CMs) ---------------------------
h <- coef(NMF_K)
h_sub <- h[c("CM01","CM02","CM03","CM04"), ]
h_sub_t <- t(h_sub)

common_samples <- intersect(rownames(h_sub_t), meta$sampleID)
h_sub_t <- h_sub_t[common_samples, ]
meta_sub <- meta[match(common_samples, meta$sampleID), ]
subtype_vec <- meta_sub$subtype

radar_mean <- h_sub_t %>%
  as.data.frame() %>%
  mutate(subtype = subtype_vec) %>%
  group_by(subtype) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

radar_long <- radar_mean %>%
  pivot_longer(cols = -subtype, names_to = "variable", values_to = "value")

# Polar to Cartesian
polar_to_cartesian <- function(r, theta_deg) {
  theta_rad <- theta_deg * pi / 180
  data.frame(x = r * cos(theta_rad), y = r * sin(theta_rad))
}

variables <- c("CM01","CM02","CM03","CM04")
n_vars <- length(variables)
angles_deg <- seq(0, 360 - 360/n_vars, length.out = n_vars)

radar_points <- radar_long %>%
  rowwise() %>%
  mutate(coord = list(polar_to_cartesian(value, angles_deg[match(variable, variables)]))) %>%
  unnest(coord) %>%
  select(subtype, variable, x, y)

curve_data <- radar_points %>%
  group_by(subtype) %>%
  arrange(match(variable, variables)) %>%
  summarise(x = c(x, first(x)), y = c(y, first(y)), .groups = "drop") %>%
  group_by(subtype) %>%
  mutate(xend = lead(x), yend = lead(y), .after = last_col()) %>%
  ungroup() %>%
  filter(!is.na(xend))

max_val <- max(radar_mean[, variables])
max_plot <- max_val * 1.05
custom_breaks <- seq(0, 6, length.out = 5)
circle_radii <- custom_breaks[custom_breaks <= max_plot]

circles <- lapply(circle_radii, function(r) {
  theta <- seq(0, 360, length.out = 100)
  polar_to_cartesian(r, theta) %>% mutate(radius = r)
}) %>% bind_rows()

rays <- data.frame()
for (angle in angles_deg) {
  end <- polar_to_cartesian(max_plot, angle)
  rays <- rbind(rays, data.frame(x = c(0, end$x), y = c(0, end$y), angle = angle))
}
rays_all <- bind_rows(lapply(unique(radar_mean$subtype), function(st) rays %>% mutate(subtype = st)))

label_r <- max_plot * 1.08
label_coords <- data.frame(
  variable = variables,
  angle = angles_deg,
  x = label_r * cos(angles_deg * pi/180),
  y = label_r * sin(angles_deg * pi/180)
)
label_all <- bind_rows(lapply(unique(radar_mean$subtype), function(st) label_coords %>% mutate(subtype = st)))

text_radial <- data.frame(
  r = circle_radii, angle = 0,
  x = circle_radii * cos(0), y = circle_radii * sin(0),
  label = as.character(circle_radii)
)
text_radial_all <- bind_rows(lapply(unique(radar_mean$subtype), function(st) text_radial %>% mutate(subtype = st)))

polygon_data <- radar_points %>%
  group_by(subtype) %>%
  arrange(match(variable, variables)) %>%
  summarise(x = c(x, first(x)), y = c(y, first(y)), .groups = "drop")

p_radar <- ggplot() +
  geom_path(data = circles, aes(x = x, y = y, group = radius),
            color = "gray50", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = rays_all, aes(x = x, y = y, group = interaction(angle, subtype)),
            color = "gray50", linetype = "dashed", linewidth = 0.5) +
  geom_polygon(data = polygon_data, aes(x = x, y = y, fill = subtype),
               alpha = 0.4, show.legend = FALSE) +
  geom_curve(data = curve_data,
             aes(x = x, y = y, xend = xend, yend = yend, color = subtype),
             curvature = -0.1, linewidth = 1, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = radar_points, aes(x = x, y = y, color = subtype),
             size = 2, show.legend = FALSE) +
  geom_text(data = label_all, aes(x = x, y = y, label = variable),
            size = 4, fontface = "bold") +
  geom_text(data = text_radial_all, aes(x = x, y = y, label = label),
            size = 3, hjust = -0.2, vjust = -0.5) +
  coord_fixed(ratio = 1) +
  facet_wrap(~ subtype, ncol = 3) +
  theme_void() +
  theme(strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  labs(title = "CM Distribution by Subtype (Filled Curves)")

pdf(g_radar_pdf, width = 12, height = 4, bg = "white")
print(p_radar)
dev.off()



# =========================================================================== #
# Figure H: Spearman correlation dot plot between functional module (CM)
# abundance and selected metabolite levels. Left columns show SuperClass and
# Class annotations.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(tidyverse)
library(readxl)
library(sva)
library(patchwork)
library(RColorBrewer)

# ------------------------------ File paths ----------------------------------
mapping1_path   <- "path/to/BQ_JWQ20240911_CJ300_RFB.xlsx"
mapping2_path   <- "path/to/BQ-JWQ20250312-LC-RFB-2-定量结果.xlsx"
metab1_path     <- "path/to/metabolome/01.xlsx"
metab2_path     <- "path/to/metabolome/03.xlsx"
prop_csv_path   <- "path/to/prop_sample1_clusters.csv"
clinical_path   <- "path/to/Clinical_information.rds"
hmdb_csv_path   <- "path/to/python_hmdb_metabolites_all.csv"
output_pdf      <- "path/to/F4/H.pdf"

# ------------------------- Load metabolite annotations ----------------------
map1 <- read_excel(mapping1_path) %>% select(1, 9)
colnames(map1) <- c("Metabolome", "HMDBID")
map2 <- read_excel(mapping2_path) %>% select(1, 9)
colnames(map2) <- c("Metabolome", "HMDBID")
mapping_all <- bind_rows(map1, map2) %>% distinct()
# Replace "/" in HMDBID with Unknown_1, Unknown_2, ...
idx <- which(mapping_all$HMDBID == "/")
mapping_all$HMDBID[idx] <- paste0("Unknown_", seq_along(idx))

name_map <- setNames(mapping_all$HMDBID, mapping_all$Metabolome)

# --------------------- Load and merge metabolite datasets -------------------
met01 <- read_excel(metab1_path)
met02 <- read_excel(metab2_path)

# Rename columns to HMDBID using the map
rename_to_hmdb <- function(df, nm) {
  old <- names(df)
  new <- nm[old]
  new[is.na(new)] <- old[is.na(new)]
  names(df) <- new
  df
}
met01 <- rename_to_hmdb(met01, name_map)
met02 <- rename_to_hmdb(met02, name_map)

# Harmonize columns
all_cols <- union(colnames(met01), colnames(met02))
met01 <- met01 %>% mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>% select(all_of(all_cols))
met02 <- met02 %>% mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>% select(all_of(all_cols))

# Prepare for batch correction
dup_samples <- intersect(met01$Samples, met02$Samples)
common_metabs <- intersect(names(met01)[-1], names(met02)[-1])

combined <- bind_rows(
  met01 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 1L),
  met02 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 2L)
) %>%
  mutate(BioID = ifelse(Samples %in% dup_samples, Samples, paste0("uniq_", row_number())))

# ComBat correction
combat_mat <- combined %>%
  select(all_of(common_metabs)) %>%
  as.matrix() %>%
  {log10(. + 1)} %>%
  t() %>%
  ComBat(batch = combined$Batch, mod = model.matrix(~1, data = combined), par.prior = TRUE) %>%
  t() %>%
  {10^. - 1}
combat_mat[combat_mat < 0] <- 0

# Consolidate duplicates and finalize matrix
metab_final <- combat_mat %>%
  as.data.frame() %>%
  mutate(Sample = combined$Samples, BioID = combined$BioID) %>%
  group_by(BioID) %>%
  summarise(Sample = first(Sample), across(all_of(common_metabs), mean), .groups = "drop") %>%
  arrange(factor(Sample, levels = c(met01$Samples, setdiff(met02$Samples, dup_samples)))) %>%
  select(Sample, all_of(common_metabs)) %>%
  column_to_rownames("Sample")

# Filter metabolites by prevalence (>0.5 in at least 20% of samples)
metab_final <- metab_final[, colSums(metab_final > 0.5) >= max(5, 0.2 * nrow(metab_final))]

# ------------------------ Build CM abundance matrix -------------------------
prop_mat <- read.csv(prop_csv_path, row.names = 1, check.names = FALSE)

cm_clusters <- list(
  CM01 = c(1, 15, 18, 20, 25, 29, 30, 41, 44),
  CM02 = c(10, 2, 3, 39, 46, 49, 6, 7, 8),
  CM03 = c(11, 13, 30, 32, 40, 47, 8, 9),
  CM04 = c(0, 12, 22, 30, 31, 36, 38, 42)
)

sample_cm <- data.frame(row.names = rownames(prop_mat))
for (cm_name in names(cm_clusters)) {
  cols <- as.character(cm_clusters[[cm_name]])
  valid_cols <- intersect(cols, colnames(prop_mat))
  sample_cm[[cm_name]] <- if (length(valid_cols)) rowSums(prop_mat[, valid_cols, drop = FALSE]) else 0
}
sample_cm <- sample_cm[, c("CM01", "CM02", "CM03", "CM04")]

# --------------------- Load clinical data, subset Pre -----------------------
group <- readRDS(clinical_path) %>% filter(treatment_condition == "Pre")

# ------------------ Common samples and Spearman correlation -----------------
common <- intersect(rownames(metab_final), intersect(rownames(sample_cm), rownames(group)))
metab_sub <- metab_final[common, ]
cm_sub    <- sample_cm[common, ]

grid_df <- expand.grid(CM = colnames(cm_sub), HMDBID = colnames(metab_sub), stringsAsFactors = FALSE)
grid_df$R <- NA; grid_df$P <- NA

for (i in seq_len(nrow(grid_df))) {
  res <- cor.test(cm_sub[[grid_df$CM[i]]], metab_sub[[grid_df$HMDBID[i]]],
                  method = "spearman", exact = FALSE)
  grid_df$R[i] <- res$estimate
  grid_df$P[i] <- res$p.value
}

# Keep only metabolites with at least one significant correlation
sig_hmdb <- unique(grid_df$HMDBID[grid_df$P < 0.05])
grid_df <- grid_df %>% filter(HMDBID %in% sig_hmdb)

grid_df$Signif <- ifelse(grid_df$P < 0.001, "***",
                         ifelse(grid_df$P < 0.01, "**",
                                ifelse(grid_df$P < 0.05, "*", "")))

# Map back to metabolite names
grid_df <- grid_df %>%
  left_join(mapping_all %>% select(HMDBID, Metabolome), by = "HMDBID")

# Filter to a specific set of metabolites
target_metab <- c("Acetic acid", "Propionic acid", "Butyric acid", "L-Proline",
                  "L-Leucine", "L-Arginine", "Indole", "Tryptamine", "Oleic acid",
                  "Palmitic acid", "Deoxycholic Acid-3-Sulfate", "Ursodeoxycholic acid",
                  "Lithocholic acid", "Chenodeoxycholic acid", "Suberic acid",
                  "Phenylacetic acid", "Homovanillic acid", "L-Dopa", "Tyramine",
                  "Sucrose", "Sorbitol", "Taurine", "Acetylcholine", "NAD+",
                  "Flavin adenine dinucleotide", "Flavin mononucleotide",
                  "Adenosine monophosphate", "Nicotinamide mononucleotide",
                  "Dephospho coenzyme A", "Thiamine diphosphate")
grid_df <- grid_df %>% filter(Metabolome %in% target_metab)

# ---------------------- Add HMDB classification info ------------------------
hmdb_class <- read.csv(hmdb_csv_path, stringsAsFactors = FALSE)
grid_df <- grid_df %>% left_join(hmdb_class, by = c("HMDBID" = "HMDB_ID"))
grid_df[is.na(grid_df)] <- "Undefined"

# -------------------------- Generate the plot -------------------------------
metabolome_order <- unique(grid_df$Metabolome)
grid_df$Metabolome <- factor(grid_df$Metabolome, levels = metabolome_order)

# SuperClass colour block
n_super <- length(unique(grid_df$SuperClass))
super_colors <- if (n_super <= 12) {
  brewer.pal(12, "Set3")[1:n_super]
} else {
  colorRampPalette(brewer.pal(12, "Set3"))(n_super)
}

p_super <- ggplot(grid_df, aes(y = Metabolome, x = "SuperClass", fill = SuperClass)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_manual(values = super_colors, name = "SuperClass") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(limits = metabolome_order, expand = c(0, 0)) +
  labs(x = "SuperClass") +
  theme_void() +
  theme(axis.text.x.top = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 20, hjust = 1),
        axis.title.x = element_text(size = 18, color = "black"),
        legend.position = "right",
        plot.margin = margin(0, 0, 0, 0))

# Class colour block
n_class <- length(unique(grid_df$Class))
class_colors <- if (n_class <= 12) {
  brewer.pal(12, "Set3")[1:n_class]
} else {
  colorRampPalette(brewer.pal(12, "Set3"))(n_class)
}

p_class <- ggplot(grid_df, aes(y = Metabolome, x = "Class", fill = Class)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_manual(values = class_colors, name = "Class") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(limits = metabolome_order, expand = c(0, 0)) +
  labs(x = "Class") +
  theme_void() +
  theme(axis.text.x.top = element_text(size = 16, color = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 18, color = "black"),
        legend.position = "right",
        plot.margin = margin(0, 0, 0, 0))

# Correlation dot plot
p_dot <- ggplot(grid_df, aes(x = CM, y = Metabolome)) +
  geom_hline(yintercept = seq(0.5, length(unique(grid_df$Metabolome)) + 0.5, by = 1),
             color = "black", linewidth = 0.5) +
  geom_vline(xintercept = seq(0.5, length(unique(grid_df$CM)) + 0.5, by = 1),
             color = "black", linewidth = 0.5) +
  geom_point(aes(fill = R), shape = 21, size = 10, color = "black", stroke = 0) +
  geom_text(aes(label = Signif), size = 5, color = "black") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1), name = "R") +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(position = "top", expand = c(0, 0.5)) +
  scale_y_discrete(limits = metabolome_order, expand = c(0, 0.5)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.text.x.top = element_text(angle = 45, hjust = 0, size = 20),
        axis.text.y = element_text(size = 20),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.key.height = unit(1.5, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(t = 30, r = 5, b = 5, l = 5, unit = "mm")) +
  coord_fixed(ratio = 0.8)

# Combine
combined_plot <- (p_super + p_class + p_dot) +
  plot_layout(nrow = 1, widths = c(0.25, 0.25, 1), guides = 'collect') +
  plot_annotation(theme = theme(plot.margin = margin(t = 10, r = 5, b = 5, l = 5, unit = "mm")))

ggsave(output_pdf, combined_plot, bg = "white", width = 14, height = 12, dpi = 300)



# =========================================================================== #
# Figure I: Heatmap of top 10 species per functional module (CM) with 
#           z-score normalized cell counts (left), and Gini coefficient 
#           lollipop chart (right).
# Figure J: Boxplots of CBEA enrichment scores for each CM across NMF subtypes
#           using metagenomic and scRNA-seq species abundance data.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(patchwork)
library(CBEA)
library(BiocSet)
library(TreeSummarizedExperiment)
library(rstatix)
library(ggpubr)
library(tibble)
library(stringr)

# ------------------------------ File paths ----------------------------------
gini_file       <- "path/to/H_gini.txt"
meta_path       <- "path/to/01metadata.rds"
clinical_path   <- "path/to/Clinical_information.rds"
species_csv     <- "path/to/prop_sample1_species.csv"
metagenome_file <- "path/to/DNA_merged_abundance_table.txt"
output_dir      <- "path/to/F4"
i_top_sp_file   <- file.path(output_dir, "I_topSP.txt")
i_heatmap_pdf   <- file.path(output_dir, "I_heatmap_species_abundance.pdf")
j_stocm_file    <- file.path(output_dir, "J_StoCM.txt")
j_metagenome_pdf <- file.path(output_dir, "J_metagenome_CEBA_baseline_03.pdf")
j_scRNA_pdf     <- file.path(output_dir, "J_scRNA_CEBA_baseline_03.pdf")

# =========================================================================== #
#                           Figure I: Species heatmap & Gini
# =========================================================================== #

# Read Gini coefficients
results <- read.table(gini_file, header = TRUE, sep = "\t",
                      quote = "", stringsAsFactors = FALSE)

# --------------- Species–CM mapping (for dominant cluster) -------------------
cm_clusters <- list(CM01 = c(1, 15, 18, 20, 25, 29, 30, 41, 44),
                    CM02 = c(10, 2, 3, 39, 46, 49, 6, 7, 8),
                    CM03 = c(11, 13, 30, 32, 40, 47, 8, 9),
                    CM04 = c(0, 12, 22, 30, 31, 36, 38, 42))

# Remove duplicated cluster assignments (keep only first occurrence)
all_numbers <- unlist(cm_clusters)
duplicated_numbers <- all_numbers[duplicated(all_numbers)]
cm_clusters <- lapply(cm_clusters, function(x) setdiff(x, duplicated_numbers))

# Convert to data frame and merge with results
cm_df <- cm_clusters %>%
  enframe(name = "CM", value = "Cluster") %>%
  unnest(Cluster)

results <- results %>%
  left_join(cm_df, by = c("Dominant_Cluster" = "Cluster")) %>%
  mutate(CM = ifelse(is.na(CM), "Unclassified", CM))

# --------------- Select top 10 species per CM by Gini -----------------------
results_top <- results %>%
  group_by(CM) %>%
  arrange(desc(Gini_CellCount), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write.table(results_top, file = i_top_sp_file, row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")

# ---------------- Prepare metadata and cell counts per CM -------------------
meta <- readRDS(meta_path)
group <- readRDS(clinical_path) %>%
  filter(treatment_condition == "Pre")

meta <- meta[meta$orig.ident %in% group$sample, ]
meta <- meta %>%
  left_join(group %>% select(sample, subtype), by = c("orig.ident" = "sample"))

# CM classification (non-overlapping clusters)
cm_clusters_clean <- list(
  CM01 = c(1, 18, 29),
  CM02 = c(10, 2, 3, 39, 46, 49, 6, 7),
  CM03 = c(11, 13, 32, 40, 47, 9),
  CM04 = c(0, 12, 22, 31, 36, 38, 42)
)

mapping <- stack(cm_clusters_clean)
colnames(mapping) <- c("seurat_clusters", "CM_class")
meta <- merge(meta, mapping, by = "seurat_clusters", all.x = TRUE)

# Count cells per CM and species
count_data <- meta %>%
  group_by(CM_class, Species) %>%
  summarise(count = n(), .groups = "drop")

count_wide <- count_data %>%
  pivot_wider(names_from = CM_class, values_from = count, values_fill = 0) %>%
  as.data.frame()
rownames(count_wide) <- count_wide$Species
count_wide <- count_wide[, -1, drop = FALSE]

# -------------- Filter and order target species -----------------------------
target_species <- rev(c(
  "Bilophila wadsworthia", "Bacteroides fragilis",
  "Enterocloster bolteae", "Bacteroides salyersiae",
  "Bacteroides sp002491635", "Enterocloster clostridioformis",
  "Bacteroides ovatus", "Bacteroides intestinalis",
  "Phocaeicola vulgatus", "Phocaeicola massiliensis",
  "Phocaeicola sp000432735", "Phocaeicola sp000436795",
  "Bacteroides sp003545565", "Parabacteroides distasonis",
  "Desulfovibrio sp900319575", "Paraprevotella clara",
  "Bacteroides sp900765785", "Sutterella wadsworthensis_A",
  "Bacteroides sp900547205", "Faecalibacterium prausnitzii_G",
  "Faecalibacterium sp900551435", "Enterocloster sp900541315",
  "Caecibacter sp003467125", "Enterocloster sp900543885",
  "Faecalibacterium prausnitzii_C", "Fusicatenibacter saccharivorans",
  "Faecalibacterium prausnitzii", "Faecalibacterium prausnitzii_D",
  "UBA11774 sp003507655", "Clostridium_Q sp900547735",
  "Klebsiella pneumoniae", "Hungatella sp005845265",
  "UBA9502 sp003506385", "UBA3402 sp003478355",
  "Enterocloster sp000155435", "Anaerotignum faecicola",
  "Citrobacter europaeus", "Citrobacter portucalensis",
  "Klebsiella variicola", "Enterocloster sp900555045"
))

present_species <- intersect(target_species, rownames(count_wide))
mat_selected <- count_wide[present_species, , drop = FALSE]
row_order <- target_species[target_species %in% present_species]
mat_selected <- mat_selected[row_order, , drop = FALSE]

desired_cm_order <- c("CM01", "CM02", "CM03", "CM04")
mat_selected <- mat_selected[, desired_cm_order, drop = FALSE]

# Z-score normalisation (per species)
mat_scaled <- t(scale(t(mat_selected)))
mat_scaled <- as.matrix(mat_scaled)

# Prepare plotting data
plot_data <- mat_scaled %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  pivot_longer(cols = starts_with("CM"), names_to = "CM_class",
               values_to = "scaled_value") %>%
  mutate(Species = factor(Species, levels = row_order),
         CM_class = factor(CM_class, levels = desired_cm_order))

# Gini data (aligned with species order)
gini_raw <- unique(results[, c("Species", "Gini_CellCount")])
gini_data <- data.frame(Species = levels(plot_data$Species))
gini_data <- merge(gini_data, gini_raw, by = "Species", all.x = TRUE, sort = FALSE)
gini_data$Species <- factor(gini_data$Species, levels = levels(plot_data$Species))

# --------------------------- Heatmap (left) ---------------------------------
my_colors <- colorRampPalette(brewer.pal(9, "RdYlGn"))(100)

p_heatmap <- ggplot(plot_data, aes(x = CM_class, y = Species, fill = scaled_value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = rev(my_colors), limits = c(-1, 1),
                       oob = scales::squish, name = "Scaled value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(r = 0, l = 5, t = 5, b = 5)) +
  labs(x = "", y = "")

# -------------------------- Lollipop chart (right) --------------------------
p_bar <- ggplot(gini_data, aes(x = Gini_CellCount, y = Species)) +
  geom_segment(aes(x = 0, xend = Gini_CellCount, y = Species, yend = Species),
               color = "steelblue", linewidth = 1, na.rm = TRUE) +
  geom_point(color = "steelblue", size = 8, na.rm = TRUE) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red", linewidth = 0.6) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(l = 0, r = 5, t = 5, b = 5)) +
  labs(x = "Gini coefficient", y = NULL)

# Combine and save
combined_i <- p_heatmap + p_bar +
  plot_layout(widths = c(3, 1), guides = "collect") &
  theme(legend.position = "right")

ggsave(i_heatmap_pdf, plot = combined_i, width = 12, height = 12, device = "pdf")


# =========================================================================== #
#           Figure J: CBEA enrichment analysis (metagenome & scRNA)
# =========================================================================== #

# ---- J.1 Prepare Species-to-CM mapping from metadata -----------------------
meta_data <- readRDS(meta_path) %>%
  group_by(Species) %>%
  filter(n() >= 100) %>%
  ungroup()

meta_data$module <- "Unclassified"
cm_clusters_j <- list(CM01 = c(1, 18, 29),
                       CM02 = c(10, 2, 3, 39, 46, 49, 6, 7),
                       CM03 = c(11, 13, 32, 40, 47, 9),
                       CM04 = c(0, 12, 22, 31, 36, 38, 42))

for (mod_name in names(cm_clusters_j)) {
  clusters_in_mod <- cm_clusters_j[[mod_name]]
  meta_data$module[meta_data$seurat_clusters %in% clusters_in_mod] <- mod_name
}

# Count per species–module
count_table <- table(meta_data$Species, meta_data$module)
prop_table <- prop.table(count_table, margin = 1)

result_df <- data.frame(
  Species = rownames(prop_table),
  CM = colnames(prop_table)[apply(prop_table, 1, which.max)],
  Dominant_CMProp = apply(prop_table, 1, max),
  stringsAsFactors = FALSE
)
rownames(result_df) <- NULL

# Merge with Gini results and filter
results_j <- read.table(gini_file, header = TRUE, sep = "\t",
                        quote = "", stringsAsFactors = FALSE)
cm_species_df <- inner_join(result_df, results_j, by = "Species")
cm_species_df <- cm_species_df %>%
  select(Species, CM, Dominant_CMProp, Gini_CellCount, Total_Cells)

write.table(cm_species_df, file = j_stocm_file, row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")

main_group <- subset(cm_species_df, Dominant_CMProp > 0.3) %>%
  select(Species, MainGroup = CM) %>%
  column_to_rownames(var = "Species")

# BiocSet object
elementset_df <- data.frame(element = rownames(main_group),
                            set = main_group$MainGroup,
                            stringsAsFactors = FALSE)
set <- BiocSet_from_elementset(elementset_df)

# ---- J.2 Metagenome CBEA ---------------------------------------------------
meta_df <- read.delim(metagenome_file, check.names = FALSE, row.names = 1)
meta_df <- meta_df %>%
  rownames_to_column(var = "Taxonomy") %>%
  filter(str_count(Taxonomy, "\\|") == 6, str_detect(Taxonomy, "\\|s__")) %>%
  mutate(Species = str_replace(Taxonomy, ".*\\|s__", "")) %>%
  select(-Taxonomy) %>%
  column_to_rownames("Species") %>%
  t() %>%
  as.data.frame()

rownames(meta_df) <- gsub("_metaphlan", "", rownames(meta_df))
specie_sample <- t(meta_df)

# Subset to pre-treatment common samples
group_j <- readRDS(clinical_path) %>%
  filter(treatment_condition == "Pre", PRPDSD != "Unknow")

common_samples <- intersect(colnames(specie_sample), group_j$sample)
specie_sample <- specie_sample[, common_samples, drop = FALSE]

rownames(specie_sample) <- gsub("_", " ", rownames(specie_sample))
rownames(specie_sample) <- gsub(" sp_", " sp", rownames(specie_sample))
rownames(specie_sample) <- sub("_.*", "", rownames(specie_sample))

seq_meta <- TreeSummarizedExperiment(
  assays = list(`16SrRNA` = as.matrix(specie_sample)),
  colData = DataFrame(sample_id = colnames(specie_sample)),
  rowData = DataFrame(SPECIES = rownames(specie_sample))
)

set.seed(123)
mod_meta <- cbea(obj = seq_meta, set = set, output = "zscore",
                 abund_values = "16SrRNA",
                 distr = "norm", parametric = TRUE,
                 adj = TRUE, thresh = 0.05, n_perm = 1000)

enrichment_meta <- do.call(cbind, mod_meta$R) %>%
  as.data.frame() %>%
  t()
colnames(enrichment_meta) <- names(mod_meta$R)

# ---- J.3 scRNA-seq CBEA -----------------------------------------------------
specie_sc <- read.csv(species_csv)
colnames(specie_sc)[1] <- "patient_id"
colnames(specie_sc) <- gsub("^X", "cluster", colnames(specie_sc))
specie_sc <- specie_sc %>%
  column_to_rownames(var = "patient_id") %>%
  t()

common_samples <- intersect(colnames(specie_sc), group_j$sample)
specie_sc <- specie_sc[, common_samples, drop = FALSE]

clean_species_names <- function(species_names) {
  cleaned <- gsub("\\.", " ", species_names)
  cleaned <- gsub("_[A-Z][a-zA-Z]*(?=\\s|$)", "", cleaned, perl = TRUE)
  cleaned
}
rownames(specie_sc) <- clean_species_names(rownames(specie_sc))

seq_sc <- TreeSummarizedExperiment(
  assays = list(`16SrRNA` = as.matrix(specie_sc)),
  colData = DataFrame(sample_id = colnames(specie_sc)),
  rowData = DataFrame(SPECIES = rownames(specie_sc))
)

set.seed(123)
mod_sc <- cbea(obj = seq_sc, set = set, output = "zscore",
               abund_values = "16SrRNA",
               distr = "norm", parametric = TRUE,
               adj = TRUE, thresh = 0.05, n_perm = 1000)

enrichment_sc <- do.call(cbind, mod_sc$R) %>%
  as.data.frame() %>%
  t()
colnames(enrichment_sc) <- names(mod_sc$R)

# ---- J.4 Plot boxplots for both data types ---------------------------------
plot_cbea <- function(enrichment_df, group_df, title_label) {
  common <- intersect(colnames(enrichment_df), group_df$sample)
  enrichment_sub <- enrichment_df[, common, drop = FALSE]

  enrichment_long <- as.data.frame(t(enrichment_sub)) %>%
    mutate(sample = rownames(.)) %>%
    pivot_longer(cols = -sample, names_to = "CM", values_to = "enrichment_score") %>%
    filter(CM != "Unclassified")

  enrichment_plot <- enrichment_long %>%
    inner_join(group_df %>% select(sample, subtype), by = "sample") %>%
    mutate(subtype = as.factor(subtype),
           CM = factor(CM, levels = c("CM01", "CM02", "CM03", "CM04")))

  stat_test <- enrichment_plot %>%
    group_by(CM) %>%
    wilcox_test(enrichment_score ~ subtype, p.adjust.method = "BH") %>%
    add_xy_position(x = "subtype", dodge = 0.8)

  stat_test_sig <- stat_test %>% filter(p.adj < 0.05)

  p <- ggplot(enrichment_plot, aes(x = subtype, y = enrichment_score)) +
    geom_boxplot(aes(fill = subtype), outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.2, size = 3, alpha = 1) +
    facet_wrap(~ CM, ncol = 4, scales = "free_y") +
    labs(x = "Subtype", y = "Enrichment Score", title = title_label) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          strip.background = element_rect(fill = "lightgray"),
          axis.text.x = element_text(angle = 0),
          panel.grid = element_blank()) +
    scale_fill_brewer(palette = "Set2")

  if (nrow(stat_test_sig) > 0) {
    p <- p + stat_pvalue_manual(stat_test_sig, label = "p.adj.signif",
                               tip.length = 0.01, step.increase = 0.1)
  }
  p
}

p_meta <- plot_cbea(enrichment_meta, group_j, "Metagenome")
p_sc   <- plot_cbea(enrichment_sc, group_j, "scRNA")

ggsave(j_metagenome_pdf, plot = p_meta, width = 15, height = 5, device = "pdf")
ggsave(j_scRNA_pdf, plot = p_sc, width = 15, height = 5, device = "pdf")



# =========================================================================== #
# Figure K: CM enrichment analysis in public metagenomic cohorts.
# Heatmap of standardized enrichment scores per sample ordered by RECIST;
# boxplot and ROC curve for CM03 distinguishing PR vs PD.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggsignif)
library(effsize)
library(pROC)
library(pheatmap)
library(RColorBrewer)
library(BiocSet)
library(TreeSummarizedExperiment)
library(CBEA)
library(curatedMetagenomicData)

# ------------------------------ File paths ----------------------------------
metadata_path  <- "path/to/metadata.rds"
gini_file      <- "path/to/H_gini.txt"
output_dir     <- "path/to/F4"

# Choose study (FrankelAE_2017, LeeKA_2022, WindTT_2020)
object_name <- "WindTT_2020"

# -------------------- Load public metagenomic data --------------------------
se_relative <- sampleMetadata |>
  filter(study_name == object_name) |>
  returnSamples("relative_abundance", rownames = "short")

abundance_mat <- assay(se_relative, "relative_abundance")
clin_df <- as.data.frame(colData(se_relative))

# Remove specific samples (known outliers)
if (object_name == "WindTT_2020") {
  abundance_mat <- abundance_mat[, colnames(abundance_mat) != "op159_M1", drop = FALSE]
}
if (object_name == "FrankelAE_2017") {
  abundance_mat <- abundance_mat[, colnames(abundance_mat) != "frank_P14", drop = FALSE]
}

# Build TreeSummarizedExperiment
seq <- TreeSummarizedExperiment(
  assays = list(Metagenomic = as.matrix(abundance_mat)),
  colData = DataFrame(sample_id = colnames(abundance_mat)),
  rowData = DataFrame(SPECIES = rownames(abundance_mat))
)

# ------------- Species to CM mapping (from scRNA-seq metadata) --------------
meta_data <- readRDS(metadata_path) %>%
  group_by(Species) %>%
  filter(n() >= 100) %>%
  ungroup()

meta_data$module <- "Unclassified"
cm_clusters <- list(
  CM01 = c(1, 18, 29),
  CM02 = c(10, 2, 3, 39, 46, 49, 6, 7),
  CM03 = c(11, 13, 32, 40, 47, 9),
  CM04 = c(0, 12, 22, 31, 36, 38, 42)
)
for (mod_name in names(cm_clusters)) {
  meta_data$module[meta_data$seurat_clusters %in% cm_clusters[[mod_name]]] <- mod_name
}

# Proportional assignment of species to CM
count_table <- table(meta_data$Species, meta_data$module)
prop_table <- prop.table(count_table, margin = 1)

result_df <- data.frame(
  Species = rownames(prop_table),
  CM = colnames(prop_table)[apply(prop_table, 1, which.max)],
  Dominant_CMProp = apply(prop_table, 1, max),
  stringsAsFactors = FALSE
)
rownames(result_df) <- NULL

# Merge with Gini and keep species with Dominant_CMProp > 0.3
results <- read.table(gini_file, header = TRUE, sep = "\t",
                      quote = "", stringsAsFactors = FALSE)
cm_species_df <- inner_join(result_df, results, by = "Species")

main_group <- cm_species_df %>%
  filter(Dominant_CMProp > 0.3) %>%
  select(Species, MainGroup = CM) %>%
  column_to_rownames("Species")

# Create BiocSet for CBEA
elementset_df <- data.frame(
  element = rownames(main_group),
  set = main_group$MainGroup,
  stringsAsFactors = FALSE
)
set <- BiocSet_from_elementset(elementset_df)

# ------------------------ CBEA enrichment analysis --------------------------
set.seed(123)
mod <- cbea(
  obj = seq, set = set, output = "zscore",
  abund_values = "Metagenomic",
  distr = "norm", parametric = TRUE,
  adj = TRUE, thresh = 0.05, n_perm = 1000
)

enrichment_matrix <- do.call(cbind, mod$R)
colnames(enrichment_matrix) <- names(mod$R)

# Remove Unclassified and scale by column (CM)
enrichment_df <- enrichment_matrix[, setdiff(colnames(enrichment_matrix), "Unclassified")]
enrichment_df <- scale(enrichment_df)

# ---------------------- Heatmap: samples ordered by RECIST ------------------
recist_vec <- clin_df[rownames(enrichment_df), "RECIST"]
recist_factor <- factor(recist_vec, levels = c("CR", "PR", "SD", "PD"))

# Order samples within each RECIST group by hierarchical clustering
all_order <- c()
for (grp in levels(recist_factor)) {
  idx_grp <- which(recist_factor == grp)
  if (length(idx_grp) == 1) {
    all_order <- c(all_order, idx_grp)
  } else {
    sub_df <- enrichment_df[idx_grp, , drop = FALSE]
    hc <- hclust(dist(sub_df), method = "complete")
    all_order <- c(all_order, idx_grp[hc$order])
  }
}

enrichment_df_sorted <- enrichment_df[all_order, , drop = FALSE]
annotation_row <- data.frame(RECIST = recist_factor[all_order])
rownames(annotation_row) <- rownames(enrichment_df_sorted)

# Color palette
spectral_colors <- rev(brewer.pal(11, "Spectral"))
my_colors <- colorRampPalette(spectral_colors)(100)
breaks <- seq(-2, 2, length.out = length(my_colors) + 1)

heatmap_file <- file.path(output_dir, paste0("K_heatmap_", object_name, ".pdf"))
pheatmap(
  enrichment_df_sorted,
  annotation_row = annotation_row,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = my_colors,
  breaks = breaks,
  main = "",
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = NA,
  filename = heatmap_file
)

# --------------- Boxplot: CM03 enrichment across RECIST groups --------------
box_data <- as.data.frame(enrichment_df)
box_data$RECIST <- clin_df[rownames(box_data), "RECIST"]
box_data <- box_data[!is.na(box_data$RECIST), ]
box_data$RECIST <- factor(box_data$RECIST, levels = c("CR", "PR", "SD", "PD"))

# Identify valid groups (>=2 samples with non‑NA CM03)
group_counts <- table(box_data$RECIST[!is.na(box_data$CM03)])
valid_groups <- names(group_counts[group_counts >= 2])

if (length(valid_groups) >= 2) {
  # Pairwise t‑tests (p < 0.1) with Cohen's d
  all_combinations <- combn(valid_groups, 2, simplify = FALSE)
  sig_comparisons <- list()
  sig_labels <- c()
  
  for (comb in all_combinations) {
    g1 <- box_data$CM03[box_data$RECIST == comb[1]]
    g2 <- box_data$CM03[box_data$RECIST == comb[2]]
    g1 <- g1[!is.na(g1)]
    g2 <- g2[!is.na(g2)]
    if (length(g1) >= 2 && length(g2) >= 2) {
      t_res <- t.test(g1, g2)
      if (t_res$p.value < 0.1) {
        d_val <- cohen.d(g1, g2)$estimate
        p_label <- ifelse(t_res$p.value < 0.001, "p < 0.001",
                          paste0("p = ", round(t_res$p.value, 3)))
        d_label <- paste0("d = ", round(d_val, 2))
        label <- paste(p_label, d_label, sep = ", ")
        sig_comparisons <- c(sig_comparisons, list(comb))
        sig_labels <- c(sig_labels, label)
      }
    }
  }
  
  p_box <- ggplot(box_data, aes(x = RECIST, y = CM03, fill = RECIST)) +
    geom_boxplot(alpha = 1, width = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 1) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none", panel.grid = element_blank()) +
    labs(x = "RECIST Classification",
         y = paste(object_name, "CM03 Enrichment Score"))
  
  if (length(sig_comparisons) > 0) {
    p_box <- p_box + geom_signif(
      comparisons = sig_comparisons,
      annotations = sig_labels,
      map_signif_level = FALSE,
      step_increase = 0.1,
      tip_length = 0.02,
      textsize = 3
    )
  }
  
  boxplot_file <- file.path(output_dir, paste0("K_boxplot_", object_name, ".pdf"))
  ggsave(boxplot_file, p_box, width = 4, height = 6, dpi = 300)
}

# ----------- ROC curve: CM03 discriminating PR vs PD ------------------------
roc_data <- box_data %>%
  filter(RECIST %in% c("PR", "PD")) %>%
  mutate(label = ifelse(RECIST == "PR", 1, 0)) %>%
  select(CM03, label) %>%
  na.omit()

if (sum(roc_data$label == 1) >= 1 && sum(roc_data$label == 0) >= 1) {
  roc_obj <- roc(label ~ CM03, data = roc_data, quiet = TRUE)
  auc_val <- auc(roc_obj)
  
  roc_file <- file.path(output_dir, paste0("K_ROC_", object_name, "_CM03_PR_vs_PD.pdf"))
  pdf(roc_file, width = 6, height = 6)
  plot(roc_obj, col = "#f8caaa", lwd = 4,
       xlab = "1 - Specificity (False Positive Rate)",
       ylab = "Sensitivity",
       main = "ROC Curve: CM03 distinguishing PR vs PD")
  legend("bottomright",
         legend = paste("ROC (AUC =", round(auc_val, 2), ")"),
         col = "#f8caaa", lwd = 4, bty = "n")
  dev.off()
}