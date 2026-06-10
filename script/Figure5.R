# =========================================================================== #
# Figure A: Data preparation for cluster composition and functional gene 
# expression. Outputs cluster-subtype proportions, cell counts, and Z-score 
# normalized CPM of selected genes across clusters.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(openxlsx)

# ------------------------------ File paths ----------------------------------
metadata_path   <- "path/to/metadata.rds"
seurat_obj_path <- "path/to/01merge.rds"
marker_file     <- "path/to/04marker_cluster_m01_l00_r05_func.tsv"
output_dir      <- "path/to/F5"
prop_file       <- file.path(output_dir, "A_1cluster_subtype_proportions.xlsx")
count_file      <- file.path(output_dir, "A_2cell_counts.xlsx")
gene_file       <- file.path(output_dir, "A_3gene_expression_by_cluster.xlsx")

# -------------------------- Load data ---------------------------------------
metadata   <- readRDS(metadata_path) %>% filter(treatment_condition == "Pre")
seurat_obj <- readRDS(seurat_obj_path)
seurat_obj@meta.data <- metadata
rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$barcode

# ------------------- 1. Cluster-subtype proportions -------------------------
metadata$cluster_label <- paste0("Cluster ", metadata$seurat_clusters)

cluster_order <- c(
  "Cluster 1", "Cluster 15", "Cluster 18", "Cluster 20", "Cluster 25", 
  "Cluster 29", "Cluster 41", "Cluster 44", "Cluster 2", "Cluster 3", 
  "Cluster 6", "Cluster 7", "Cluster 10", "Cluster 39", "Cluster 46", 
  "Cluster 49", "Cluster 9", "Cluster 11", "Cluster 13", "Cluster 32", "Cluster 33",
  "Cluster 40", "Cluster 47", "Cluster 0", "Cluster 12", "Cluster 22", 
  "Cluster 31", "Cluster 36", "Cluster 38", "Cluster 42"
)
cluster_order_unique <- unique(cluster_order)

result <- metadata %>%
  group_by(cluster_label, subtype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster_label) %>%
  mutate(proportion = count / sum(count)) %>%
  select(-count) %>%
  pivot_wider(names_from = cluster_label,
              values_from = proportion,
              values_fill = 0) %>%
  arrange(subtype)

existing_clusters <- intersect(cluster_order_unique, colnames(result))
result_ordered <- result %>% select(subtype, all_of(existing_clusters))
write.xlsx(result_ordered, prop_file, rowNames = FALSE)

# ----------------------- 2. Cell counts per cluster -------------------------
cluster_counts <- table(metadata$seurat_clusters)
cluster_numbers <- as.numeric(gsub("Cluster ", "", specified_order))

result_df <- data.frame(
  Cluster = specified_order,
  Count = sapply(cluster_numbers, function(x) {
    if (as.character(x) %in% names(cluster_counts)) cluster_counts[as.character(x)] else 0
  }),
  stringsAsFactors = FALSE
)
write.xlsx(result_df, count_file, row.names = FALSE)

# ------------------- 3. Gene expression Z-score matrix ----------------------
genes <- c("clpB", "dnak", "groL", "htpG", "tpx", "katA", "susC-10", "mglB-2",
           "susC-24", "susD", "susF", "susE", "manZ-2", "sugC", "lacE", "cbgA-1",
           "lacF-1", "tpl", "gdh", "alst-5", "gdh-4", "hutu-1", "prdA-4", "grdE",
           "fhs", "nif", "purR-6", "mutB-2", "mutB-1", "pccB-2", "scpA-2",
           "cat1-1", "gcdA", "bcd", "gctA", "mdh", "gap", "pyk", "eno",
           "glpk-1", "sspC2", "spoIVA", "spoVB")

# Average raw counts per cluster (SCT assay, slot = "counts")
raw_avg <- AverageExpression(seurat_obj,
                             assays = "SCT",
                             features = genes,
                             group.by = "seurat_clusters",
                             slot = "counts")$SCT

available_genes <- intersect(genes, rownames(raw_avg))
raw_avg <- raw_avg[available_genes, ]
colnames(raw_avg) <- paste("Cluster", colnames(raw_avg))

# CPM normalization
calculate_cpm <- function(count_matrix) {
  cluster_sums <- colSums(count_matrix)
  t(t(count_matrix) / cluster_sums) * 1e6
}
cpm_data <- calculate_cpm(raw_avg)

# log2(CPM + 1) transformation
logcpm_data <- log2(cpm_data + 1)

# Z-score per gene
zscore_data <- t(scale(t(logcpm_data)))
zscore_data[is.na(zscore_data)] <- 0

# Clip extreme values at 98th percentile of absolute values
data_abs <- abs(zscore_data)
data_range <- quantile(data_abs[data_abs > 0], probs = 0.98, na.rm = TRUE)
if (is.finite(data_range)) {
  zscore_data[zscore_data >  data_range] <-  data_range
  zscore_data[zscore_data < -data_range] <- -data_range
}

# Reorder columns by cluster_order, keep only existing
existing_clusters <- cluster_order[cluster_order %in% colnames(zscore_data)]
zscore_data <- zscore_data[, existing_clusters, drop = FALSE]
zscore_data <- zscore_data[available_genes, , drop = FALSE]

final_data <- as.data.frame(zscore_data)
final_data <- cbind(Gene = rownames(final_data), final_data)

# Merge KEGG annotation
marker <- read.delim(marker_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
marker <- marker[, c("gene", "KEGG")]
marker <- marker %>%
  group_by(gene) %>%
  summarise(KEGG = paste(unique(KEGG), collapse = ","), .groups = "drop")

final_data_with_kegg <- final_data %>%
  left_join(marker, by = c("Gene" = "gene")) %>%
  select(Gene, KEGG, everything())

write.xlsx(final_data_with_kegg, gene_file, row.names = FALSE)



# =========================================================================== #
# Figure B: Heatmap of Z-score normalized mean abundance of selected KEGG
# modules across NMF subtypes in pre-treatment metagenomic samples.
# Modules were selected by Wilcoxon test (p<0.05) and effect size per subtype,
# with additional manually curated modules added.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(tidyverse)
library(readr)
library(pheatmap)
library(effsize)
library(RColorBrewer)

# ------------------------------ File paths ----------------------------------
uniref_ko_path    <- "path/to/map_ko_uniref90.txt"
kegg_modules_path <- "path/to/KEGG_modules.tab"
genefamilies_path <- "path/to/combined_genefamilies.tsv"
clinical_path     <- "path/to/Clinical_information.rds"
output_pdf        <- "path/to/F5/B.pdf"

# ================== Part 1: Build KEGG module abundance matrix ==============

# Read KO-to-UniRef90 mapping
uniref_ko_map <- read_tsv(uniref_ko_path, col_names = FALSE, show_col_types = FALSE) %>%
  rename(ko = X1) %>%
  pivot_longer(-ko, names_to = "col_name", values_to = "uniref90") %>%
  select(-col_name) %>%
  filter(!is.na(uniref90) & uniref90 != "") %>%
  distinct()

# Read KEGG module definitions and split KO list
kegg_modules <- read_tsv(kegg_modules_path,
                         col_names = c("module", "description", "kos"),
                         show_col_types = FALSE) %>%
  separate_rows(kos, sep = ";") %>%
  filter(!is.na(kos) & kos != "") %>%
  distinct()

# Merge to get UniRef90–module mapping
full_map <- uniref_ko_map %>%
  inner_join(kegg_modules, by = c("ko" = "kos")) %>%
  select(uniref90, ko, module)

# Read gene family abundances and reshape to long format
genefamilies <- read_tsv(genefamilies_path, comment = "#", show_col_types = FALSE)
genefamilies_long <- genefamilies %>%
  pivot_longer(-`Gene Family`, names_to = "sample", values_to = "abundance") %>%
  rename(uniref90 = `Gene Family`) %>%
  filter(!grepl("^UNMAPPED$|^UNGROUPED$", uniref90))

# Join with mapping and sum abundance per module per sample
module_abundance <- genefamilies_long %>%
  inner_join(full_map, by = "uniref90") %>%
  group_by(module, sample) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# Pivot to wide matrix and add clean module identifiers
final_table <- module_abundance %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = 0) %>%
  left_join(kegg_modules %>% select(module, description) %>% distinct(), by = "module") %>%
  mutate(
    function_name = str_replace(description, "\\s*\\[.*\\]$", ""),
    function_name = str_replace(function_name, "=>.*$", ""),
    module_id = paste0(module, "_", function_name)
  ) %>%
  select(module_id, everything(), -module, -description, -function_name)

# Convert to data frame with module_id as row names
kegg_module_data <- as.data.frame(final_table)
rownames(kegg_module_data) <- kegg_module_data$module_id
kegg_module_data <- kegg_module_data[, -1]

# Clean sample names (remove suffix)
colnames(kegg_module_data) <- gsub("_2_clean_Abundance-RPKs", "", colnames(kegg_module_data))

# Normalize to relative abundance (sample-wise)
kegg_module_data <- sweep(kegg_module_data, 2, colSums(kegg_module_data), "/")

# Filter low-abundance modules (mean > 0.0001)
kegg_module_data <- kegg_module_data[rowMeans(kegg_module_data) > 0.0001, ]
# Re-normalize after filtering
kegg_module_data <- sweep(kegg_module_data, 2, colSums(kegg_module_data), "/")

# Transpose to samples × modules
kegg_mat <- t(kegg_module_data)

# ================== Part 2: Subtype association analysis ====================

# Load clinical data, keep only Pre-treatment with known subtype
group <- readRDS(clinical_path) %>%
  filter(treatment_condition == "Pre") %>%
  rename(model = subtype) %>%
  filter(model != "Unknow")

# Intersect samples
common_samples <- intersect(rownames(kegg_mat), group$sample)
kegg_mat <- kegg_mat[common_samples, ]
group <- group[group$sample %in% common_samples, ]

# Prepare long-format data
kegg_long <- kegg_mat %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "module", values_to = "abundance") %>%
  left_join(group, by = "sample")

# Wilcoxon test + Cohen's d for each module vs each subtype
clinical_models <- unique(group$model)
modules <- colnames(kegg_mat)

results <- data.frame()
for (m in clinical_models) {
  for (mod in modules) {
    mod_data <- kegg_long %>%
      filter(module == mod) %>%
      mutate(group = ifelse(model == m, "target", "others"))
    wt <- wilcox.test(abundance ~ group, data = mod_data)
    cd <- cohen.d(mod_data$abundance[mod_data$group == "target"],
                  mod_data$abundance[mod_data$group == "others"])$estimate
    results <- rbind(results, data.frame(
      module = mod,
      clinical_model = m,
      p_value = wt$p.value,
      cohen_d = cd,
      stringsAsFactors = FALSE
    ))
  }
}

# FDR correction
results$fdr <- p.adjust(results$p_value, method = "fdr")

# Select top 1 module per subtype by effect size (p < 0.05)
sig_results <- results %>%
  group_by(clinical_model) %>%
  filter(p_value < 0.05) %>%
  slice_max(order_by = cohen_d, n = 1) %>%
  ungroup()

# Manually add extra modules of interest (prefix list)
target_prefixes <- c(
  "M00023","M00089","M00191","M00224","M00256","M00330","M00332","M00444","M00477","M00490",
  "M00500","M00260","M00033","M00088","M00118","M00119","M00170","M00269","M00539","M00543","M00547",
  "M00119","M00024","M00025","M00040","M00113","M00503","M00281","M00002","M00003","M00004",
  "M00005","M00006","M00007","M00008","M00009","M00010","M00011","M00012","M00015","M00019",
  "M00093","M00115","M00120","M00140","M00141","M00172","M00133","M00157","M00177","M00526","M00527",
  "M00019","M00050","M00051","M00083","M00088","M00115","M00120","M00133","M00134","M00140",
  "M00152","M00170","M00176","M00240","M00250","M00260","M00299","M00311","M00377","M00394",
  "M00527","M00570","M00572","M00609","map02060","map04141","M00024","M00025","M00023","M00040",
  "M00136","M00336","M00435","M00545"
)

all_modules <- colnames(kegg_mat)
target_full_names <- all_modules[grepl(paste(target_prefixes, collapse = "|"), all_modules)]
modules_to_add <- setdiff(target_full_names, sig_results$module)

if (length(modules_to_add) > 0) {
  sig_results <- bind_rows(
    sig_results,
    data.frame(module = modules_to_add,
               clinical_model = NA_character_,
               p_value = NA_real_,
               cohen_d = NA_real_,
               fdr = NA_real_,
               stringsAsFactors = FALSE)
  )
}

# ================== Part 3: Heatmap preparation =============================

# Calculate mean abundance per module per subtype
heatmap_mean <- kegg_long %>%
  filter(module %in% sig_results$module) %>%
  group_by(module, model) %>%
  summarise(mean_value = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  right_join(expand.grid(module = unique(sig_results$module),
                         model = unique(group$model),
                         stringsAsFactors = FALSE),
             by = c("module", "model")) %>%
  mutate(mean_value = replace_na(mean_value, 0)) %>%
  pivot_wider(names_from = model, values_from = mean_value, values_fill = 0) %>%
  column_to_rownames("module")

colnames(heatmap_mean) <- paste0("Subtype ", colnames(heatmap_mean))
heatmap_mean <- heatmap_mean[, order(as.numeric(gsub("Subtype ", "", colnames(heatmap_mean))))]

# Z-score per row (module)
mat_scaled <- t(scale(t(as.matrix(heatmap_mean))))
mat_scaled[is.na(mat_scaled)] <- 0

# Color settings
legend_breaks <- seq(-1, 1, length.out = 100)
my_colors <- colorRampPalette(c("#a4bbd7", "#835398"))(length(legend_breaks) - 1)

# Draw and save heatmap
p <- pheatmap(mat_scaled,
              scale = "none",
              color = my_colors,
              breaks = legend_breaks,
              border_color = "grey60",
              show_colnames = TRUE,
              show_rownames = TRUE,
              clustering_method = "ward.D2",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              cluster_cols = FALSE,
              fontsize_row = 8,
              fontsize_col = 10,
              cellwidth = 20,
              angle_col = 45)

ggsave(output_pdf, plot = p, width = 10, height = 10, dpi = 300, bg = "white")



# =========================================================================== #
# Figure C: Heatmap of Z-score normalized mean abundance of selected serum 
# metabolites across NMF subtypes. Metabolites were identified by Wilcoxon 
# test (p < 0.05) and Cohen's d for each subtype, with additional manually 
# curated metabolites included. Row annotations show SuperClass and Class 
# from HMDB.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(tidyverse)
library(readxl)
library(sva)
library(effsize)
library(pheatmap)
library(RColorBrewer)

# ------------------------------ File paths ----------------------------------
mapping1_path    <- "path/to/BQ_JWQ20240911_CJ300_RFB.xlsx"
mapping2_path    <- "path/to/BQ-JWQ20250312-LC-RFB-2-定量结果.xlsx"
metab1_path      <- "path/to/metabolome/01.xlsx"
metab2_path      <- "path/to/metabolome/03.xlsx"
clinical_path    <- "path/to/Clinical_information.rds"
hmdb_class_path  <- "path/to/python_hmdb_metabolites_all.csv"
output_pdf       <- "path/to/F5/C.pdf"

# ================== 1. Load and harmonize metabolite data ===================

# Build HMDB ID mapping from two annotation tables
map1 <- read_excel(mapping1_path) %>% select(1, 9)
colnames(map1) <- c("Metabolome", "HMDBID")
map2 <- read_excel(mapping2_path) %>% select(1, 9)
colnames(map2) <- c("Metabolome", "HMDBID")
mapping_all <- bind_rows(map1, map2) %>% distinct()
# Replace "/" with Unknown IDs
idx <- which(mapping_all$HMDBID == "/")
mapping_all$HMDBID[idx] <- paste0("Unknown_", seq_along(idx))

name_map <- setNames(mapping_all$HMDBID, mapping_all$Metabolome)

# Function to rename columns using the map
rename_to_hmdb <- function(df, nm) {
  old <- names(df)
  new <- nm[old]
  new[is.na(new)] <- old[is.na(new)]
  names(df) <- new
  df
}

# Read raw metabolomics data
met01 <- read_excel(metab1_path)
met02 <- read_excel(metab2_path)
met01 <- rename_to_hmdb(met01, name_map)
met02 <- rename_to_hmdb(met02, name_map)

# Harmonize columns across batches
all_cols <- union(colnames(met01), colnames(met02))
met01 <- met01 %>% mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>% select(all_of(all_cols))
met02 <- met02 %>% mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>% select(all_of(all_cols))

# Merge and batch-correct
dup_samples <- intersect(met01$Samples, met02$Samples)
common_metabs <- intersect(names(met01)[-1], names(met02)[-1])

combined <- bind_rows(
  met01 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 1L),
  met02 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 2L)
) %>%
  mutate(BioID = ifelse(Samples %in% dup_samples, Samples, paste0("uniq_", row_number())))

combat_mat <- combined %>%
  select(all_of(common_metabs)) %>%
  as.matrix() %>%
  {log10(. + 1)} %>%
  t() %>%
  ComBat(batch = combined$Batch, mod = model.matrix(~1, data = combined), par.prior = TRUE) %>%
  t() %>%
  {10^. - 1}
combat_mat[combat_mat < 0] <- 0

# Consolidate duplicates and create final matrix
metab_final <- combat_mat %>%
  as.data.frame() %>%
  mutate(Sample = combined$Samples, BioID = combined$BioID) %>%
  group_by(BioID) %>%
  summarise(Sample = first(Sample), across(all_of(common_metabs), mean), .groups = "drop") %>%
  arrange(factor(Sample, levels = c(met01$Samples, setdiff(met02$Samples, dup_samples)))) %>%
  select(Sample, all_of(common_metabs)) %>%
  column_to_rownames("Sample")

# Filter low-abundance metabolites (present >0.5 in at least 20% of samples)
metab_final <- metab_final[, colSums(metab_final > 0.5) >= max(5, 0.2 * nrow(metab_final))]

# ================== 2. Load clinical data ===================================
group <- readRDS(clinical_path) %>%
  filter(treatment_condition == "Pre") %>%
  rename(model = subtype) %>%
  filter(model != "Unknow")

# Intersect samples
common_samples <- intersect(rownames(metab_final), group$sample)
metab_sub <- metab_final[common_samples, ]
group_sub <- group[group$sample %in% common_samples, ]

# Convert to long format
metab_long <- metab_sub %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "metabolite", values_to = "abundance") %>%
  left_join(group_sub, by = "sample")

# ================== 3. Wilcoxon test and effect size ========================
clinical_models <- unique(group_sub$model)
metabolites <- colnames(metab_sub)

results <- data.frame()
for (m in clinical_models) {
  for (met in metabolites) {
    met_data <- metab_long %>%
      filter(metabolite == met) %>%
      mutate(group = ifelse(model == m, "target", "others"))
    wt <- wilcox.test(abundance ~ group, data = met_data)
    cd <- cohen.d(met_data$abundance[met_data$group == "target"],
                  met_data$abundance[met_data$group == "others"])$estimate
    results <- rbind(results, data.frame(
      metabolite = met,
      clinical_model = m,
      p_value = wt$p.value,
      cohen_d = cd,
      stringsAsFactors = FALSE
    ))
  }
}
results$fdr <- p.adjust(results$p_value, method = "fdr")

# ================== 4. Select significant and forced metabolites ============
sig_results <- results %>%
  group_by(clinical_model) %>%
  filter(p_value < 0.05) %>%
  arrange(desc(cohen_d)) %>%
  ungroup()

# Define metabolites to force include (using original names)
forced_keywords <- unique(c(
  "Butyric acid", "Fructose-6-phosphate", "Dihydroxyacetone phosphate",
  "Succinic acid", "L-Leucine", "L-Phenylalanine", "L-Arginine",
  "Acetylcholine", "Tryptamine", "Tyramine", "Imidazolepropionic acid",
  "7-Ketolithocholic acid", "Ursocholic acid", "β-Muricholic acid", "Xylose",
  "L-(+)-rhamnose", "Methylimidazoleacetic acid", "3β-Cholic Acid",
  "Ursodeoxycholic acid 3-Sulfate", "Homovanillic acid", "Sucrose",
  "7-Ketodeoxycholic acid", "Phenylacetic acid", "2-Methylbenzoic acid",
  "8,11,14_Eicosatrienoic acid", "Docosapentaenoic acid DPA",
  "beta-Hydroxyisovaleric acid", "Nonoic acid", "L-Proline", "L-Valine",
  "L-Methionine", "Phenylpyruvic acid", "L-Serine", "2-Oxohexanoic acid",
  "Indole-3-carboxylic acid", "Dodecanoic acid", "N-Acetyl-L-tyrosine",
  "Apocholic acid", "4-Methyl-2-oxopentanoic acid", "L-Dopa",
  "Glyceric acid", "Ethylmalonic acid"
))
# Map forced keywords to HMDB IDs
forced_hmdb <- mapping_all$HMDBID[match(forced_keywords, mapping_all$Metabolome)]
forced_rows <- results %>% filter(tolower(metabolite) %in% tolower(forced_hmdb))
missing_forced <- forced_rows %>% anti_join(sig_results, by = c("clinical_model", "metabolite"))
sig_results <- bind_rows(sig_results, missing_forced) %>% arrange(clinical_model, desc(cohen_d))

# ================== 5. Prepare heatmap data =================================
sample_models <- group_sub[, c("sample", "model")] %>% distinct()
all_comb <- expand.grid(metabolite = unique(sig_results$metabolite),
                        model = unique(sample_models$model),
                        stringsAsFactors = FALSE)

heatmap_mean <- metab_long %>%
  filter(metabolite %in% sig_results$metabolite) %>%
  group_by(metabolite, model) %>%
  summarise(mean_value = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  right_join(all_comb, by = c("metabolite", "model")) %>%
  mutate(mean_value = replace_na(mean_value, 0)) %>%
  pivot_wider(names_from = model, values_from = mean_value, values_fill = 0) %>%
  column_to_rownames("metabolite")

colnames(heatmap_mean) <- paste0("Subtype ", colnames(heatmap_mean))
heatmap_mean <- heatmap_mean[, order(as.numeric(gsub("Subtype ", "", colnames(heatmap_mean))))]

# Replace HMDB IDs with actual metabolite names and add HMDB classification
hmdb_class <- read.csv(hmdb_class_path, stringsAsFactors = FALSE)
heatmap_mean <- heatmap_mean %>%
  rownames_to_column("HMDBID") %>%
  left_join(mapping_all[, c("HMDBID", "Metabolome")], by = "HMDBID") %>%
  left_join(hmdb_class, by = c("HMDBID" = "HMDB_ID")) %>%
  column_to_rownames("Metabolome")
heatmap_mean[is.na(heatmap_mean)] <- "Undefined"

# Extract annotation columns and heatmap matrix
annotation_row <- heatmap_mean[, c("SuperClass", "Class")]
mat_data <- heatmap_mean[, paste0("Subtype ", 1:3)]

# ================== 6. Build annotation colors ==============================
super_classes <- unique(annotation_row$SuperClass)
class_types   <- unique(annotation_row$Class)
n_super <- length(super_classes)
n_class <- length(class_types)

if (n_super <= 12) {
  super_colors <- brewer.pal(n_super, "Paired")
} else {
  super_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_super)
}
if (n_class <= 12) {
  class_colors <- brewer.pal(n_class, "Paired")
} else {
  class_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_class)
}
ann_colors <- list(
  SuperClass = setNames(super_colors, super_classes),
  Class = setNames(class_colors, class_types)
)

# ================== 7. Draw and save heatmap ================================
mat_scaled <- t(scale(t(as.matrix(mat_data))))
mat_scaled[is.na(mat_scaled)] <- 0

legend_breaks <- seq(-1, 1, length.out = 100)
my_colors <- colorRampPalette(c("#a4bbd7", "#835398"))(length(legend_breaks) - 1)

pheatmap(
  mat = mat_scaled,
  scale = "none",
  color = my_colors,
  breaks = legend_breaks,
  border_color = "grey60",
  show_colnames = TRUE,
  show_rownames = TRUE,
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  cellwidth = 20,
  angle_col = 45,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  annotation_legend = TRUE,
  legend = TRUE,
  filename = output_pdf,
  width = 10,
  height = 10,
  bg = "white"
)



# =========================================================================== #
# Figures D–G: Dot plots showing average expression and percent expressed of
# selected functional gene sets across NMF subtypes in pre-treatment scRNA-seq.
# D: Carbohydrate transport (sus, lac, mgl, sugC, etc.)
# E: Central carbon & fermentation (pckA, frdA, mutB, etc.)
# F: Amino acid metabolism (gdh, aspA, tpl, tdcB, etc.)
# G: Stress response & redox (sodB, katA, clpB, etc.)
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# ------------------------------ File paths ----------------------------------
scRNA_rds_path  <- "path/to/01merge.rds"
metadata_path   <- "path/to/metadata.rds"
output_dir      <- "path/to/F5"

# --------------------------- Load and filter data ---------------------------
scRNA <- readRDS(scRNA_rds_path)
metadata <- readRDS(metadata_path)
scRNA@meta.data <- metadata
rownames(scRNA@meta.data) <- scRNA@meta.data$barcode

# Keep only pre‑treatment cells with known subtype
scRNA <- subset(scRNA, treatment_condition == "Pre")
scRNA <- subset(scRNA, subtype != "Unknow")

# =========================================================================== #
#                        Figure D: Carbohydrate transport
# =========================================================================== #
genes_D <- c("lacE", "lacF", "lacZ", "cbgA", "mglB", "manZ",
             "sugC", "susB", "susC", "susD", "susE", "susF")
valid_D <- intersect(genes_D, rownames(scRNA))

expr_D <- GetAssayData(scRNA, assay = "RNA", slot = "data")[valid_D, , drop = FALSE]
subtype_D <- scRNA@meta.data$subtype

df_D <- data.frame(cell = colnames(expr_D), subtype = subtype_D)
for (g in valid_D) df_D[[g]] <- expr_D[g, ]
df_long_D <- df_D %>%
  pivot_longer(cols = all_of(valid_D), names_to = "gene", values_to = "expression")

plot_data_D <- df_long_D %>%
  group_by(subtype, gene) %>%
  summarise(avg_exp = mean(expression),
            pct_exp = sum(expression > 0) / n() * 100,
            .groups = "drop")

p_D <- ggplot(plot_data_D, aes(x = factor(subtype), y = gene)) +
  geom_point(aes(size = pct_exp, color = avg_exp)) +
  scale_color_gradientn(colours = c("white", "#5D90BA", "#4a1486"),
                        name = "Average expression") +
  scale_size_continuous(range = c(1, 10), name = "Percent expressed") +
  labs(x = "Subtype", y = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "gray60", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm")) +
  coord_flip()

ggsave(file.path(output_dir, "D_DotPlot_genes_by_subtype.pdf"),
       p_D, width = 4 + 0.3 * length(valid_D), height = 3, device = "pdf")

# =========================================================================== #
#           Figure E: Central carbon metabolism and fermentation
# =========================================================================== #
gene_order_E <- c("pckA", "frdA", "mutB", "pccB", "scpA", "pct",
                  "lcdAB", "ppdK", "ldh", "gctA", "por", "bcd", "cat1")
valid_E <- intersect(gene_order_E, rownames(scRNA))

expr_E <- GetAssayData(scRNA, assay = "RNA", slot = "data")[valid_E, , drop = FALSE]
subtype_E <- scRNA@meta.data$subtype

df_E <- data.frame(cell = colnames(expr_E), subtype = subtype_E)
for (g in valid_E) df_E[[g]] <- expr_E[g, ]
df_long_E <- df_E %>%
  pivot_longer(cols = all_of(valid_E), names_to = "gene", values_to = "expression")

plot_data_E <- df_long_E %>%
  group_by(subtype, gene) %>%
  summarise(avg_exp = mean(expression),
            pct_exp = sum(expression > 0) / n() * 100,
            .groups = "drop") %>%
  mutate(gene = factor(gene, levels = valid_E),
         avg_exp_capped = pmax(pmin(avg_exp, 0.1), 0),
         pct_exp_capped = pmin(pct_exp, 7))

p_E <- ggplot(plot_data_E, aes(x = factor(subtype), y = gene)) +
  geom_point(aes(size = pct_exp_capped, color = avg_exp_capped)) +
  scale_color_gradientn(colours = c("white", "#5D90BA", "#4a1486"),
                        name = "Average expression", limits = c(0, 0.1)) +
  scale_size_continuous(range = c(1, 12), name = "Percent expressed",
                        breaks = c(5, 6, 7), labels = c("5","6",">=7"),
                        limits = c(0, 7)) +
  labs(x = "Subtype", y = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "gray60", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm")) +
  coord_flip()

ggsave(file.path(output_dir, "E_DotPlot_genes_by_subtype.pdf"),
       p_E, width = 6 + 0.3 * length(valid_E), height = 3, device = "pdf")

# =========================================================================== #
#                    Figure F: Amino acid metabolism
# =========================================================================== #
genes_F <- c("gdh", "grdE", "fhs", "aspA", "purR", "prdA", "tpl", "tdcB")
valid_F <- intersect(genes_F, rownames(scRNA))

expr_F <- GetAssayData(scRNA, assay = "RNA", slot = "data")[valid_F, , drop = FALSE]
subtype_F <- scRNA@meta.data$subtype

df_F <- data.frame(cell = colnames(expr_F), subtype = subtype_F)
for (g in valid_F) df_F[[g]] <- expr_F[g, ]
df_long_F <- df_F %>%
  pivot_longer(cols = all_of(valid_F), names_to = "gene", values_to = "expression")

plot_data_F <- df_long_F %>%
  group_by(subtype, gene) %>%
  summarise(avg_exp = mean(expression),
            pct_exp = sum(expression > 0) / n() * 100,
            .groups = "drop")

p_F <- ggplot(plot_data_F, aes(x = factor(subtype), y = gene)) +
  geom_point(aes(size = pct_exp, color = avg_exp)) +
  scale_color_gradientn(colours = c("white", "#5D90BA", "#4a1486"),
                        name = "Average expression") +
  scale_size_continuous(range = c(1, 14), name = "Percent expressed") +
  labs(x = "Subtype", y = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "gray60", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm")) +
  coord_flip()

ggsave(file.path(output_dir, "F_DotPlot_genes_by_subtype.pdf"),
       p_F, width = 4 + 0.3 * length(valid_F), height = 3, device = "pdf")

# =========================================================================== #
#                  Figure G: Stress response and redox
# =========================================================================== #
gene_order_G <- c("sodB", "katA", "clpB", "groL", "htpG",
                  "aphC", "tpx", "zwf", "trxB", "dps")
valid_G <- intersect(gene_order_G, rownames(scRNA))

expr_G <- GetAssayData(scRNA, assay = "RNA", slot = "data")[valid_G, , drop = FALSE]
subtype_G <- scRNA@meta.data$subtype

df_G <- data.frame(cell = colnames(expr_G), subtype = subtype_G)
for (g in valid_G) df_G[[g]] <- expr_G[g, ]
df_long_G <- df_G %>%
  pivot_longer(cols = all_of(valid_G), names_to = "gene", values_to = "expression")

plot_data_G <- df_long_G %>%
  group_by(subtype, gene) %>%
  summarise(avg_exp = mean(expression),
            pct_exp = sum(expression > 0) / n() * 100,
            .groups = "drop") %>%
  mutate(gene = factor(gene, levels = valid_G),
         avg_exp_capped = pmax(pmin(avg_exp, 0.2), 0),
         pct_exp_capped = pmin(pct_exp, 20))

p_G <- ggplot(plot_data_G, aes(x = factor(subtype), y = gene)) +
  geom_point(aes(size = pct_exp_capped, color = avg_exp_capped)) +
  scale_color_gradientn(colours = c("white", "#5D90BA", "#4a1486"),
                        name = "Average expression", limits = c(0, 0.2)) +
  scale_size_continuous(range = c(1, 12), name = "Percent expressed",
                        breaks = c(1, 10, 20), labels = c("1","10",">=20"),
                        limits = c(0, 20)) +
  labs(x = "Subtype", y = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "gray60", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm")) +
  coord_flip()

ggsave(file.path(output_dir, "G_DotPlot_genes_by_subtype.pdf"),
       p_G, width = 6 + 0.3 * length(valid_G), height = 3, device = "pdf")