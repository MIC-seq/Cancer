#!/usr/bin/env Rscript

# =============================================================================
# Functional gene matrix construction and species annotation for 10X data
# Usage: Rscript this_script.R --ref_genome <file> --gene_annot <file> \
#          --input_dir <dir> --output_base <dir> --samples <sample1,sample2,...>
# =============================================================================

library(Seurat)
library(dplyr)
library(stringr)
library(Matrix)
library(matrixStats)

# ---- Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

# Helper function to extract argument value
get_arg <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NA)
  args[idx + 1]
}

ref_genome_file <- get_arg(args, "--ref_genome")
gene_annot_file  <- get_arg(args, "--gene_annot")
input_dir        <- get_arg(args, "--input_dir")
output_base      <- get_arg(args, "--output_base")
samples_str      <- get_arg(args, "--samples")

# Check required arguments
if (any(is.na(c(ref_genome_file, gene_annot_file, input_dir, output_base, samples_str)))) {
  stop("Missing required arguments. Usage:\n",
       "Rscript script.R --ref_genome <file> --gene_annot <file> \\\n",
       "  --input_dir <dir> --output_base <dir> --samples <sample1,sample2,...>")
}

# Parse samples into a vector
samples <- strsplit(samples_str, ",")[[1]]
samples <- trimws(samples)  # remove possible whitespace

# =============================================================================
# Load reference data
# =============================================================================

# Read genome metadata and add Phocaeicola vulgatus
ref_genome <- read.table(ref_genome_file, sep = "\t", header = TRUE)
ref_genome[nrow(ref_genome) + 1, ] <- NA
ref_genome[nrow(ref_genome), c("Genome", "Lineage")] <- list(
  "BVU",
  "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus"
)
rownames(ref_genome) <- ref_genome$Genome

# Read gene annotation file
gene_back <- read.table(gene_annot_file, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, quote = "")
colnames(gene_back) <- c("gene", "function")

# =============================================================================
# Main loop over samples
# =============================================================================
for (sample in samples) {
  outpath <- file.path(output_base, "01each", sample)
  if (!dir.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
  }
  
  # Read 10X data
  input_sample_dir <- file.path(input_dir, sample)
  scRNA.data <- Read10X(data.dir = input_sample_dir, gene.column = 1)
  scrna_combine <- CreateSeuratObject(counts = scRNA.data, project = sample,
                                      min.cell = 1, min.features = 1)
  
  # Extract count matrix
  scrna_combine_count_matrix <- scrna_combine@assays$RNA@counts
  
  # ---- Step 1: Assign each barcode to the most abundant gene's species ----
  mgyg_name <- strsplit(rownames(scrna_combine_count_matrix), '-')
  name_sp <- sapply(mgyg_name, `[`, 1)  # get the first part (genome tag)
  
  barcode_sp_functionumap <- c()
  for (n in seq_len(ncol(scrna_combine_count_matrix))) {
    if (n %% 5000 == 0) cat("Processing barcode", n, "\n")
    not_zero_name <- name_sp[which(scrna_combine_count_matrix[, n] != 0)]
    freq_table <- table(not_zero_name)
    barcode_sp_functionumap[n] <- names(freq_table)[which.max(freq_table)]
  }
  barcode_sp_functionumap <- as.data.frame(barcode_sp_functionumap)
  rownames(barcode_sp_functionumap) <- colnames(scrna_combine_count_matrix)
  
  # Map genome tag to full species information
  barcode_sp_functionumap_addsp <- ref_genome[ref_genome$Genome %in% 
                                               barcode_sp_functionumap$barcode_sp_functionumap, ]
  barcode_sp_functionumap_addsp$Lineage <- strsplit(barcode_sp_functionumap_addsp$Lineage, split = ";")
  barcode_sp_functionumap_addsp$species_info <- barcode_sp_functionumap_addsp$Genome
  for (i in seq_along(barcode_sp_functionumap_addsp$Lineage)) {
    barcode_sp_functionumap_addsp$species_info[i] <- barcode_sp_functionumap_addsp$Lineage[[i]][7]
  }
  barcode_sp_functionumap_addsp$species_info <- gsub("s__", "", barcode_sp_functionumap_addsp$species_info)
  barcode_sp_functionumap_addsp$species_info <- trimws(barcode_sp_functionumap_addsp$species_info)
  barcode_sp_functionumap_addsp <- barcode_sp_functionumap_addsp[barcode_sp_functionumap_addsp$species_info != "", ]
  
  # Merge to get species per barcode
  barcode_gene <- barcode_sp_functionumap
  barcode_gene$barcode <- rownames(barcode_gene)
  gene_barcode_species <- merge(barcode_gene, barcode_sp_functionumap_addsp,
                                by.x = "barcode_sp_functionumap", by.y = "Genome", all.x = TRUE)
  gene_barcode_species <- subset(gene_barcode_species, !is.na(species_info))
  
  # ---- Step 2: Build functional gene annotation and filter matrix ----
  total_back <- rownames(scrna_combine@assays$RNA@counts)
  total_back <- data.frame(gene = total_back, stringsAsFactors = FALSE)
  total_back$gene <- gsub("-", "_", total_back$gene)
  
  # Adjust BVU gene names to match annotation keys
  gene_num <- total_back
  gene_num$gene <- gsub("(BVU)([A-Za-z0-9]+)", "\\1_\\2", gene_num$gene)
  gene_num <- merge(gene_back, gene_num, by = "gene")
  
  # Parse annotation fields
  get_function <- gene_num
  get_function_l <- strsplit(gene_num$`function`, split = ";")
  gene_f <- matrix(NA, nrow = length(get_function_l), ncol = 4)
  
  for (i in seq_along(get_function_l)) {
    for (si in seq_along(get_function_l[[i]])) {
      if (str_detect(get_function_l[[i]][si], "KEGG"))         gene_f[i, 1] <- get_function_l[[i]][si]
      if (str_detect(get_function_l[[i]][si], "Ontology_term")) gene_f[i, 2] <- get_function_l[[i]][si]
      if (str_detect(get_function_l[[i]][si], "Name"))          gene_f[i, 3] <- get_function_l[[i]][si]
      if (str_detect(get_function_l[[i]][si], "product"))       gene_f[i, 4] <- get_function_l[[i]][si]
    }
  }
  
  colnames(gene_f) <- c("KEGGnumber", "Ontology_term", "genename", "product")
  get_function <- cbind(get_function, gene_f)
  get_function <- get_function %>%
    mutate(
      genename = ifelse(genename == "Name=", NA, genename),
      product  = ifelse(product == "product=", NA, product)
    )
  
  # Filter out uninformative genes
  get_function <- get_function[!is.na(get_function$genename), ]
  get_function <- get_function[!is.na(get_function$product), ]
  get_function <- subset(get_function, Ontology_term != "Ontology_term=" | is.na(Ontology_term))
  get_function <- get_function[get_function$KEGGnumber != "KEGG=-", ]
  get_function <- get_function[!duplicated(get_function$gene), ]
  id <- grep("ribosomal", get_function$product)
  if (length(id) > 0) get_function <- get_function[-id, ]
  
  # Clean annotation strings
  get_function$KEGGnumber    <- gsub("KEGG=", "", get_function$KEGGnumber)
  get_function$Ontology_term <- gsub("Ontology_term=", "", get_function$Ontology_term)
  get_function$genename      <- gsub("Name=", "", get_function$genename)
  get_function$product       <- gsub("product=", "", get_function$product)
  get_function$gene          <- gsub("_", "-", get_function$gene)
  
  # Subset count matrix to functional genes and replace IDs with gene names
  func_gene_matrix <- scrna_combine@assays$RNA@counts
  rownames(func_gene_matrix) <- gsub("BVU", "BVU-", rownames(func_gene_matrix))
  func_gene_matrix <- func_gene_matrix[rownames(func_gene_matrix) %in% get_function$gene, ]
  func_gene_matrix <- func_gene_matrix[, colSums(func_gene_matrix) > 0]
  
  gene_map <- setNames(as.character(get_function$genename), as.character(get_function$gene))
  rownames(func_gene_matrix) <- gene_map[rownames(func_gene_matrix)]
  rownames(func_gene_matrix) <- gsub("_", "-", rownames(func_gene_matrix))
  
  # Aggregate by gene name (take maximum per gene)
  dense_mat <- as.matrix(func_gene_matrix)
  row_names <- rownames(func_gene_matrix)
  gene_groups <- split(seq_len(nrow(dense_mat)), row_names)
  unique_genes <- names(gene_groups)
  
  new_mat <- matrix(0, nrow = length(unique_genes), ncol = ncol(dense_mat))
  rownames(new_mat) <- unique_genes
  colnames(new_mat) <- colnames(dense_mat)
  
  for (gene in unique_genes) {
    rows <- gene_groups[[gene]]
    if (length(rows) == 1) {
      new_mat[gene, ] <- dense_mat[rows, ]
    } else {
      new_mat[gene, ] <- colMaxs(dense_mat[rows, , drop = FALSE])
    }
  }
  
  sparse_mat_new <- as(new_mat, "dgCMatrix")
  sparse_mat_new <- sparse_mat_new[, colSums(sparse_mat_new) != 0]
  
  # Save functional gene matrix
  saveRDS(sparse_mat_new, file = file.path(outpath, "0func_gene_matrix.rds"))
  
  # ---- Step 3: Save species annotation and summary ----
  report <- gene_barcode_species[, c("barcode", "species_info")]
  write.table(report, file = file.path(outpath, "1report_specie.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  species_counts <- as.data.frame(table(report$species_info))
  colnames(species_counts) <- c("Species", "Count")
  write.table(species_counts, file = file.path(outpath, "2count_Species.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # ---- Step 4: Create Seurat object with species metadata ----
  scRNA <- CreateSeuratObject(counts = sparse_mat_new, project = sample,
                              min.cell = 1, min.features = 1)
  scRNA@meta.data$barcode <- rownames(scRNA@meta.data)
  scRNA@meta.data <- left_join(scRNA@meta.data, report, by = "barcode",
                               relationship = "many-to-many")
  scRNA@meta.data <- scRNA@meta.data %>% filter(!is.na(species_info))
  
  serial_cellnum_list <- scRNA@meta.data %>%
    group_by(species_info) %>%
    summarise(cellnum = n(), .groups = "drop") %>%
    arrange(species_info) %>%
    mutate(serial = as.integer(factor(species_info)) - 1)
  
  scRNA@meta.data <- left_join(scRNA@meta.data, serial_cellnum_list,
                               by = "species_info", relationship = "many-to-many")
  scRNA@meta.data$serial_cellnum <- paste0(scRNA@meta.data$species_info, " (",
                                           scRNA@meta.data$cellnum, ")")
  rownames(scRNA@meta.data) <- scRNA@meta.data$barcode
  
  saveRDS(scRNA, file = file.path(outpath, "3scobj_withspecies.rds"))
  
  cat("Finished processing sample:", sample, "\n")
}