# Load necessary libraries
library(Seurat)
library(SeuratObject)
library(dplyr)
library(stringr)

# Load reference genome metadata
ref_genome <- read.table("path/to/genomes-all_metadata_2.0.tsv", sep = "\t", header = TRUE)

# Manually add Phocaeicola vulgatus (BVU) entry
ref_genome[nrow(ref_genome) + 1, ] <- NA  
ref_genome[nrow(ref_genome), c("Genome", "Lineage")] <- list(
  "BVU", 
  "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola vulgatus"
)
rownames(ref_genome) <- ref_genome$Genome

# 49 samples
samples <- c("sample1","sample2","sample3")

for (sample in samples) {
  outpath <- paste0("path/to/output/01each/", sample, "/")
  if (!dir.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
  }

  # Read 10X data
  input_file <- paste0("path/to/input/matrix/", sample, "/")
  scRNA.data <- Read10X(data.dir = input_file, gene.column = 1)

  # Create Seurat object (minimal filtering)
  scrna_combine <- CreateSeuratObject(counts = scRNA.data, project = sample, 
                                      min.cell = 1, min.features = 1)

  # Extract filtered count matrix
  scrna_combine_count_matrix <- scrna_combine@assays$RNA@counts

  # Extract genome IDs from rownames (format: "GenomeID-geneName")
  mgyg_name <- strsplit(rownames(scrna_combine_count_matrix), '-')
  name_sp <- sapply(mgyg_name, `[`, 1)

  # For each barcode, find the most frequently expressed genome
  barcode_sp_functionumap <- c()
  for (n in 1:ncol(scrna_combine_count_matrix)) {  
    if (n %% 5000 == 0) { print(n) }
    not_zero_name <- name_sp[which(scrna_combine_count_matrix[, n] != 0)]
    freq_table <- table(not_zero_name)
    barcode_sp_functionumap[n] <- names(freq_table)[which.max(freq_table)]
  }
  barcode_sp_functionumap <- as.data.frame(barcode_sp_functionumap)
  rownames(barcode_sp_functionumap) <- colnames(scrna_combine_count_matrix)

  # Prepare mapping: barcode -> dominant genome
  barcode_gene <- barcode_sp_functionumap
  barcode_gene$barcode <- rownames(barcode_gene)

  # Subset reference to only genomes actually observed
  genomes_used <- intersect(ref_genome$Genome, barcode_sp_functionumap$barcode_sp_functionumap)
  lineage_info <- ref_genome[genomes_used, ]
  lineage_split <- strsplit(lineage_info$Lineage, split = ";")

  # Extract taxonomy levels (Phylum, Family, Genus, Species)
  lineage_info$Phylum  <- sapply(lineage_split, function(x) trimws(gsub("p__", "", x[2])))
  lineage_info$Family  <- sapply(lineage_split, function(x) trimws(gsub("f__", "", x[5])))
  lineage_info$Genus   <- sapply(lineage_split, function(x) trimws(gsub("g__", "", x[6])))
  lineage_info$Species <- sapply(lineage_split, function(x) trimws(gsub("s__", "", x[7])))

  # Keep only records with at least one valid taxonomic annotation
  lineage_info <- lineage_info[lineage_info$Phylum != "" | lineage_info$Family != "" | 
                               lineage_info$Genus != "" | lineage_info$Species != "", ]

  # Merge taxonomy with barcode-genome mapping
  gene_barcode_anno <- merge(barcode_gene, 
                             lineage_info[, c("Genome", "Phylum", "Family", "Genus", "Species")],
                             by.x = "barcode_sp_functionumap", by.y = "Genome",
                             all.x = TRUE)

  # Output annotation tables for each taxonomy level
  report_phylum  <- gene_barcode_anno[!is.na(gene_barcode_anno$Phylum),  c("barcode", "Phylum")]
  report_family  <- gene_barcode_anno[!is.na(gene_barcode_anno$Family),  c("barcode", "Family")]
  report_genus   <- gene_barcode_anno[!is.na(gene_barcode_anno$Genus),   c("barcode", "Genus")]
  report_species <- gene_barcode_anno[!is.na(gene_barcode_anno$Species), c("barcode", "Species")]

  write.table(report_phylum,  file = paste0(outpath, "1report_Phylum.txt"),  sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(report_family,  file = paste0(outpath, "1report_Family.txt"),  sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(report_genus,   file = paste0(outpath, "1report_Genus.txt"),   sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(report_species, file = paste0(outpath, "1report_Specie.txt"),  sep = "\t", quote = FALSE, row.names = FALSE)
}