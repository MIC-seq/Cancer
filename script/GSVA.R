# =========================================================================== #
# GSVA analysis pipeline: pathway enrichment scores for microbial clusters.
# Steps: 1) Build gene annotation with KEGG/GO; 2) Merge with cluster markers;
# 3) Retrieve KEGG pathway names; 4) Construct GMT file; 5) Run ssGSEA;
# 6) Identify pathway markers per cluster.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(GSVA)
library(GSEABase)
library(Matrix)
library(KEGGREST)
library(stringr)
library(reshape2)

# ------------------------------ File paths ----------------------------------
input_rds           <- "path/to/01merge.rds"
gene_back_file      <- "path/to/df5_new_adPv.txt"
gene_annot_rds      <- "path/to/df5_new_adPv_adNpKG.data.rds"
output_dir          <- "path/to/GSVA_output"
marker_file         <- file.path(output_dir, "01marker_cluster_m01_l00_r05.txt")
marker_func_file    <- file.path(output_dir, "02marker_cluster_m01_l00_r05_func.tsv")
map_file            <- file.path(output_dir, "03funcmarker_map.tsv")
gmt_file            <- file.path(output_dir, "04gmt.txt")
gsva_rds            <- file.path(output_dir, "gsva_output.rds")
gsva_marker_file    <- file.path(output_dir, "05GSVA_marker.txt")

# ======================= 1. Build gene annotation ===========================
gene_back <- read.table(gene_back_file, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, quote = "")
colnames(gene_back) <- c("gene", "function")

ad_name <- gene_back
ad_name$Name <- str_extract(ad_name$`function`, "(?<=Name=)[^;]+")
ad_name <- ad_name[!is.na(ad_name$Name) & ad_name$Name != "-", ]
ad_name$product <- str_extract(ad_name$`function`, "(?<=product=)[^;]+")
ad_name <- ad_name[!is.na(ad_name$product) & ad_name$product != "-", ]
ad_name$KEGG <- str_extract(ad_name$`function`, "(?<=KEGG=)[^;]+")
ad_name$GO <- str_extract(ad_name$`function`, "(?<=Ontology_term=)[^;]+")
ad_name <- ad_name %>% mutate(across(where(is.character), ~ gsub("_", "-", .)))
saveRDS(ad_name, gene_annot_rds)

# ==================== 2. Merge markers with annotation ======================
ad_name <- readRDS(gene_annot_rds)
marker <- read.table(marker_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
marker_func <- left_join(marker, ad_name, by = c("gene" = "Name"))

# Remove duplicate columns and keep first KEGG entry
marker_func <- marker_func[, !names(marker_func) %in% c("function", "gene.y")]
marker_func <- marker_func %>% distinct()

combine_unique <- function(x) {
  unique_values <- unique(x[!is.na(x)])
  if (length(unique_values) == 0) return(NA)
  paste(unique_values, collapse = ";")
}
marker_func <- marker_func %>%
  group_by(across(1:7)) %>%
  summarise(across(everything(), combine_unique), .groups = "drop")

marker_func$KEGG <- sapply(marker_func$KEGG, function(x) {
  if (is.na(x)) return(NA)
  matches <- regmatches(x, regexpr("ko:K\\d+", x))
  if (length(matches) > 0) matches[1] else NA
})

# ---------- Retrieve KEGG pathways for each gene ----------
kegg_path <- matrix(NA, nrow = nrow(marker_func), ncol = 20)
for (i in seq_len(nrow(marker_func))) {
  kegg_entries <- unlist(strsplit(as.character(marker_func$KEGG[i]), "[,;]"))
  kegg_entries <- gsub("^ko:|-", "", kegg_entries)
  kegg_entries <- kegg_entries[kegg_entries != ""]
  all_pathways <- character(0)
  for (kid in kegg_entries) {
    if (nchar(kid) > 0) {
      query <- tryCatch(keggGet(kid), error = function(e) NULL)
      if (!is.null(query) && !is.null(query[[1]]$PATHWAY)) {
        all_pathways <- c(all_pathways, unname(query[[1]]$PATHWAY))
      }
    }
  }
  if (length(all_pathways) == 0) {
    kegg_path[i, 1] <- "nopath"
  } else {
    unique_pathways <- unique(all_pathways)
    n <- min(length(unique_pathways), 20)
    kegg_path[i, seq_len(n)] <- unique_pathways[seq_len(n)]
  }
}
kegg_path_n <- as.data.frame(kegg_path)
kegg_path_n <- cbind(marker_func, kegg_path_n)
write.table(kegg_path_n, file = marker_func_file, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# ================= 3. Fetch KEGG pathway names (map) ========================
kegg_path_n <- read.delim(marker_func_file)
kegg_path <- matrix(NA, nrow = nrow(kegg_path_n), ncol = 20)

for (i in seq_len(nrow(kegg_path_n))) {
  kegg_id <- kegg_path_n$KEGG[i]
  if (is.na(kegg_id) | kegg_id == "") {
    kegg_path[i, 1] <- "invalid_id"
    next
  }
  tryCatch({
    query <- keggGet(kegg_id)
    if (is.null(query[[1]]$PATHWAY)) {
      kegg_path[i, 1] <- "nopath"
    } else {
      pathways <- names(query[[1]]$PATHWAY)
      kegg_path[i, seq_along(pathways)] <- pathways
    }
  }, error = function(e) {
    kegg_path[i, 1] <- "error"
  })
}

colnames(kegg_path) <- paste0("W", seq_len(ncol(kegg_path)))
kegg_path_selg <- cbind(kegg_path_n, kegg_path)
kegg_path_selg <- kegg_path_selg[, c(1:10, 31:50)]
kegg_path_selg[is.na(kegg_path_selg)] <- ""
colnames(kegg_path_selg) <- gsub("^W(\\d+)$", "V\\1", colnames(kegg_path_selg))
write.table(kegg_path_selg, file = map_file, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# ===================== 4. Build GMT file ====================================
kegg_path_n <- read.delim(map_file)
kegg_path_n <- kegg_path_n[, c(7, 9, 11:30)]
kegg_path_n <- kegg_path_n %>% distinct()
rownames(kegg_path_n) <- kegg_path_n$gene
kegg_path_n[kegg_path_n == ""] <- NA
kegg_path_n <- subset(kegg_path_n, V1 != "nopath")

kegg_path_selg <- kegg_path_n[, 3:min(22, ncol(kegg_path_n))]

cc <- unique(na.omit(unlist(kegg_path_selg)))
list_path <- matrix(NA, nrow = length(cc), ncol = nrow(kegg_path_selg))
rownames(list_path) <- cc

for (i in seq_len(nrow(list_path))) {
  for (si in seq_len(nrow(kegg_path_selg))) {
    non_na_vals <- na.omit(unlist(kegg_path_selg[si, ]))
    if (rownames(list_path)[i] %in% non_na_vals) {
      list_path[i, si] <- rownames(kegg_path_selg)[si]
    }
  }
}
list_path <- as.data.frame(list_path)
colnames(list_path) <- kegg_path_n$KEGG

kegg_pathway_anno <- c()
for (i in seq_len(nrow(list_path))) {
  path1 <- keggGet(rownames(list_path)[i])
  path2 <- path1[[1]]$PATHWAY_MAP
  kegg_pathway_anno <- rbind(kegg_pathway_anno, path2)
}
list_path$pathwayanno <- kegg_pathway_anno[, 1]
list_path_selg <- list_path

data_gmt <- matrix(NA, nrow = nrow(list_path_selg), ncol = 3)
for (i in seq_len(nrow(list_path_selg))) {
  data_gmt[i, 1] <- list_path_selg$pathwayanno[i]
  data_gmt[i, 2] <- rownames(list_path_selg)[i]
  data_f <- na.omit(t(list_path_selg[i, seq_len(ncol(list_path_selg) - 1)]))
  qq <- paste(data_f[, 1], collapse = " ")
  data_gmt[i, 3] <- qq
}
data_gmt <- as.data.frame(data_gmt)
data_gmt$V3 <- gsub(" ", "\\t", data_gmt$V3)
rownames(data_gmt) <- data_gmt$V1
data_gmt <- na.omit(data_gmt)
write.table(data_gmt, gmt_file, sep = "\t", col.names = FALSE,
            quote = FALSE, row.names = FALSE)

# ====================== 5. Run ssGSEA =======================================
crc_matrix <- readRDS(input_rds)
countexp2 <- as.matrix(GetAssayData(crc_matrix, assay = "RNA", slot = "counts"))
geneSets <- getGmt(gmt_file)

ssgseaPar <- ssgseaParam(as.matrix(countexp2), geneSets)
ssgseaScores <- gsva(ssgseaPar)

signature_exp <- data.frame(ssgseaScores)
saveRDS(signature_exp, gsva_rds)

# ================= 6. Cluster-specific pathway markers ======================
gsva_seurat <- CreateSeuratObject(counts = as(as.matrix(signature_exp), "sparseMatrix"),
                                  project = "gsva", min.cells = 1, min.features = 1)
scRNA <- readRDS(input_rds)
gsva_seurat <- AddMetaData(gsva_seurat, scRNA@meta.data$seurat_clusters,
                           col.name = "seurat_clusters")
Idents(gsva_seurat) <- "seurat_clusters"
new_levels <- sort(as.numeric(levels(Idents(gsva_seurat))))
Idents(gsva_seurat) <- factor(Idents(gsva_seurat), levels = new_levels)

marker_pathways <- FindAllMarkers(gsva_seurat, only.pos = TRUE,
                                  min.pct = 0.1, logfc.threshold = 0,
                                  return.thresh = 0.05, verbose = FALSE)
write.table(marker_pathways, file = gsva_marker_file, sep = "\t",
            quote = FALSE, row.names = FALSE)