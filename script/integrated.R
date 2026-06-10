library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(ggplot2)
library(Scillus)
library(paletteer) 
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(clustree)
library(Matrix)
library(ggrepel)
library(dplyr)
library(ggsci)
library(tidydr)
library(gridExtra)
library(viridis)

input_rds <- "/path/to/matrix/"
output_path <- "/path/to/output/"
sample_names <- c("sample1","sample2","sample3")

# 定义一个通用函数，用于读取数据并创建 Seurat 对象
process_sample <- function(sample, base_path) {
  matrix_data <- readRDS(paste0(base_path, sample, "/0func_gene_matrix.rds"))
  seurat_obj <- CreateSeuratObject(counts = matrix_data, project = sample, min.cell = 1, min.features = 10)
  species_file <- paste0(base_path, sample, "/1report_Specie.txt")
  if (file.exists(species_file)) {
    species_info <- read.delim(species_file, header = TRUE, stringsAsFactors = FALSE)
    rownames(species_info) <- species_info$barcode
    meta_data <- seurat_obj@meta.data
    meta_data$species_info <- "Unknown"
    matched_cells <- intersect(rownames(meta_data), rownames(species_info))
    meta_data[matched_cells, "species_info"] <- species_info[matched_cells, "species_info"]
    seurat_obj@meta.data <- meta_data
  } else {
    warning(paste("Warning: File not found:", species_file))
    seurat_obj@meta.data$species_info <- "Unknown"
  }
  return(seurat_obj)
}

# 处理原始文件夹中的样本
scRNAlist <- lapply(sample_names, function(sample) {
  process_sample(sample, input_rds)
})

all_scRNAlist <- c(scRNAlist)
scRNA_merge <- merge(x = all_scRNAlist[[1]], y = all_scRNAlist[-1])
scRNA_merge <- subset(scRNA_merge, subset = species_info != "Homo sapiens")
scRNA_merge@meta.data <- scRNA_merge@meta.data %>% mutate(species_info = str_replace(species_info, fixed("[Ruminococcus] gnavus"), "Mediterraneibacter gnavus"))
scRNA_merge <- SCTransform(scRNA_merge, verbose = FALSE)
scRNA_merge <- RunPCA(scRNA_merge, verbose = FALSE, npcs = 100)
scRNA_merge <- FindNeighbors(scRNA_merge, dims = 1:40, verbose = FALSE)
scRNA_merge <- FindClusters(scRNA_merge, resolution = 1, verbose = FALSE)
scRNA_merge <- RunUMAP(scRNA_merge, dims = 1:40, verbose = FALSE)
scRNA_merge <- RunTSNE(scRNA_merge, dims = 1:40, check_duplicates = FALSE, verbose = FALSE)
#saveRDS(scRNA_merge, file = paste0(output_path, "01merge.rds"))
output_file <- paste0(outpath,"01metadata.txt")
write.table(scRNA_merge@meta.data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)