###############################################################################
# Figure B: Circular phylogenetic tree with family highlights, phylum ring,  #
#           and optionally a species abundance bar chart on the outer ring.  #
###############################################################################

# ----------------------------- Load required packages -----------------------
library(dplyr)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(readr)
library(tidyverse)
library(ggnewscale)
library(tidytree)
library(ggstance)

# ----------------------------- Define file paths (replace with actual paths) --
genome_map_file    <- "path/to/01genome_map_v2.1.tsv"
reference_genomes_dir <- "path/to/reference_genomes"
metadata_file      <- "path/to/metadata.rds"
tree_file          <- "path/to/reference_genomes.tre.treefile"

output_families_pdf  <- "path/to/circular_tree_top10_families.pdf"
output_abundance_pdf <- "path/to/circular_tree_top10_families_abundance.pdf"

# ----------------------------- Prepare genome-to-species mapping -------------
# Read mapping table: column1 = genome ID, column2 = species name
df <- read_tsv(genome_map_file, col_names = FALSE)

# Get genome IDs that actually have a downloaded reference file
file_names <- list.files(reference_genomes_dir, pattern = "\\.fna\\.gz$")
match_strings <- gsub("\\.fna\\.gz$", "", file_names)
df <- df %>% filter(X1 %in% match_strings)

# ----------------------------- Load metadata and compute species abundance ---
metadata <- readRDS(metadata_file)
metadata <- metadata[metadata$Species != "Unknown", ]

# Calculate species abundance as percentage of total occurrences
new_df <- metadata %>%
  group_by(Phylum, Family, Species) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(abundance = (count / sum(count)) * 100) %>%
  select(-count) %>%
  arrange(desc(abundance))

# Merge abundance data with genome IDs
new_df <- new_df %>% left_join(df, by = c("Species" = "X2"))

# ----------------------------- Read and inspect the phylogenetic tree --------
tree <- read.tree(tree_file)
tree$tip.label <- trimws(tree$tip.label)
cat("Number of tips in the tree:", length(tree$tip.label), "\n")

# Clean genome ID column
new_df$X1 <- trimws(new_df$X1)
matched <- intersect(new_df$X1, tree$tip.label)
cat("Number of matched tips:", length(matched), "\n")

# ----------------------------- Identify family-level nodes (MRCA) -----------
phylo_tree <- tree

family_nodes <- new_df %>%
  group_by(Family) %>%
  summarise(tips = list(X1)) %>%
  rowwise() %>%
  mutate(
    node = if (length(tips) > 1) {
      tip_vec <- unlist(tips)
      missing <- setdiff(tip_vec, phylo_tree$tip.label)
      if (length(missing) > 0) {
        warning("Family ", Family, " contains tips not in tree: ",
                paste(missing, collapse = ", "), " - skipping")
        NA_integer_
      } else {
        node_val <- tryCatch(getMRCA(phylo_tree, tip_vec), error = function(e) NULL)
        if (is.null(node_val)) {
          warning("Unable to compute MRCA for family ", Family, " - skipping")
          NA_integer_
        } else {
          node_val
        }
      }
    } else {
      NA_integer_
    }
  ) %>%
  filter(!is.na(node)) %>%
  ungroup()

cat("Number of families with computable MRCA:", nrow(family_nodes), "\n")

# Remove the root node if present
all_parents <- unique(tree$edge[, 1])
all_children <- unique(tree$edge[, 2])
root_node <- setdiff(all_parents, all_children)[1]
family_nodes_filtered <- family_nodes[family_nodes$node != root_node, ]

# Select top 10 families by number of species
tip_counts <- lengths(family_nodes_filtered$tips)
top10_idx <- order(tip_counts, decreasing = TRUE)[1:10]
family_top10 <- family_nodes_filtered[top10_idx, ]
family_top10$Family <- as.factor(family_top10$Family)

# ----------------------------- Prepare tip metadata (phylum) -----------------
tip_phylum <- new_df %>%
  filter(X1 %in% tree$tip.label) %>%
  select(label = X1, Phylum) %>%
  distinct()
tip_phylum$Phylum <- as.factor(tip_phylum$Phylum)

# Generate consistent colour palettes
family_colors <- scales::hue_pal()(nlevels(family_top10$Family))
names(family_colors) <- levels(family_top10$Family)

phylum_colors <- scales::hue_pal()(nlevels(tip_phylum$Phylum))
names(phylum_colors) <- levels(tip_phylum$Phylum)

# ----------------------------- Plot 1: Family highlights + phylum ring -------
p_top10 <- ggtree(tree, layout = "circular", size = 0.3) +
  geom_hilight(data = family_top10, aes(node = node, fill = Family), alpha = 0.2) +
  scale_fill_manual(name = "Family", values = family_colors)

p_top10 <- p_top10 %<+% tip_phylum

p_top10 <- p_top10 +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = Phylum),
    width = 0.05,
    offset = 0.05,
    color = NA
  ) +
  scale_fill_manual(name = "Phylum", values = phylum_colors) +
  theme(legend.position = "right")

ggsave(output_families_pdf, p_top10, width = 16, height = 14, limitsize = FALSE)

# ----------------------------- Prepare abundance data (sqrt transform) ------
abundance_data <- new_df %>%
  select(label = X1, abundance) %>%
  mutate(label = factor(label, levels = tree$tip.label),
         abundance_sqrt = sqrt(abundance)) %>%
  arrange(label) %>%
  filter(!is.na(abundance))

# ----------------------------- Plot 2: Add abundance bar chart on outer ring
p_top10_v2 <- ggtree(tree, layout = "circular", size = 0.3) +
  geom_hilight(data = family_top10, aes(node = node, fill = Family), alpha = 0.2) +
  scale_fill_manual(name = "Family", values = family_colors)

p_top10_v2 <- p_top10_v2 %<+% tip_phylum

p_top10_v2 <- p_top10_v2 +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = Phylum),
    width = 0.05,
    offset = 0.05,
    color = NA
  ) +
  scale_fill_manual(name = "Phylum", values = phylum_colors) +
  theme(legend.position = "right")

p_top10_v2 <- p_top10_v2 +
  new_scale_fill() +
  geom_fruit(
    data = abundance_data,
    geom = geom_bar,
    mapping = aes(x = abundance_sqrt, y = label),
    stat = "identity",
    orientation = "y",
    width = 0.8,
    offset = 0.3,
    fill = "steelblue"
  ) +
  scale_x_continuous(
    name = "Abundance (square root transformed)",
    breaks = pretty(abundance_data$abundance_sqrt, n = 5)
  ) +
  theme(legend.position = "right")

ggsave(output_abundance_pdf, p_top10_v2, width = 16, height = 14, limitsize = FALSE)



###############################################################################
# Figure C: Phylum-level relative abundance boxplots with Kruskal-Wallis test #
###############################################################################

# ----------------------------- Load required packages -----------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ----------------------------- Define file paths (replace with actual paths) --
metadata_path <- "path/to/metadata.rds"
output_path   <- "path/to/C_phylum_boxplots_smRNA.pdf"

# ----------------------------- Load and filter metadata -----------------------
metadata <- readRDS(metadata_path)
metadata <- subset(metadata, PRPDSD != "Unknow")

# Print column names and unique phyla for inspection
cat("Column names:", colnames(metadata), "\n")
cat("Unique phyla:", unique(metadata$Phylum), "\n")

# ----------------------------- Target phyla -----------------------------------
target_phyla <- c("Bacteroidota", "Firmicutes", "Proteobacteria", "Fusobacteriota")

# ----------------------------- Compute relative abundance per sample ----------
# Total cell counts per sample
sample_totals <- metadata %>%
  count(orig.ident) %>%
  rename(total_cells = n)

all_samples <- unique(metadata$orig.ident)

# Ensure all combinations of sample and target phylum exist
all_combinations <- expand.grid(
  orig.ident = all_samples,
  Phylum = target_phyla,
  stringsAsFactors = FALSE
)

# Count cells per sample and phylum, then calculate relative abundance
phylum_counts <- metadata %>%
  filter(Phylum %in% target_phyla) %>%
  count(orig.ident, Phylum) %>%
  right_join(all_combinations, by = c("orig.ident", "Phylum")) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  left_join(sample_totals, by = "orig.ident") %>%
  mutate(relative_abundance = n / total_cells)

# Add group information (PRPDSD)
phylum_abundance <- phylum_counts %>%
  left_join(distinct(metadata, orig.ident, PRPDSD), by = "orig.ident")

# ----------------------------- Subset to target phyla -------------------------
filtered_data <- phylum_abundance %>% filter(Phylum %in% target_phyla)

# ----------------------------- Kruskal-Wallis test per phylum -----------------
p_values <- filtered_data %>%
  group_by(Phylum) %>%
  kruskal_test(relative_abundance ~ PRPDSD) %>%
  mutate(
    p_label = case_when(
      p < 0.001 ~ "p < 0.001",
      p < 0.01  ~ paste0("p = ", round(p, 3)),
      TRUE       ~ paste0("p = ", round(p, 2))
    ),
    Phylum_label = Phylum
  )

# Merge p-value labels back into data
filtered_data <- filtered_data %>%
  left_join(p_values %>% select(Phylum, Phylum_label), by = "Phylum")

# ----------------------------- Set factor levels for plotting -----------------
phylum_label_order <- p_values %>%
  arrange(match(Phylum, target_phyla)) %>%
  pull(Phylum_label)
filtered_data$Phylum_label <- factor(filtered_data$Phylum_label,
                                     levels = phylum_label_order)
filtered_data$PRPDSD <- factor(filtered_data$PRPDSD,
                               levels = c("PR", "SD", "PD"))

# ----------------------------- Build boxplot ----------------------------------
p <- ggplot(filtered_data, aes(x = PRPDSD, y = relative_abundance, fill = PRPDSD)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 1) +
  geom_jitter(width = 0.2, size = 3, alpha = 1, color = "black") +
  facet_wrap(~ Phylum_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("PR" = "#A82E25", "PD" = "#EB7E35", "SD" = "#F5CBBF")) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.01", "0.1", "1")
  ) +
  labs(
    x = "Treatment Group",
    y = "Relative Abundance (log10 scale)",
    title = "smRNA"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ----------------------------- Save and display -------------------------------
ggsave(output_path, plot = p, width = 14, height = 6, dpi = 300)
print(p)




###############################################################################
# Figure D: Combined plots (top bar, heatmap, right bar) for family-level and #
#           genus-level abundance with Ro/e values                            #
###############################################################################

# ----------------------------- Load required packages -----------------------
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)

# ----------------------------- Define file paths (replace with actual paths) --
metadata_path     <- "path/to/metadata.rds"
output_family_pdf <- "path/to/D_combined_family_plot_roe.pdf"
output_genus_pdf  <- "path/to/D_combined_genus_plot_roe.pdf"

# ----------------------------- Load and filter metadata -----------------------
metadata <- readRDS(metadata_path)
metadata <- subset(metadata, PRPDSD != "Unknow")

# Sample identifier column (adjust if necessary)
sample_col <- "orig.ident"

# ============================================================================
# PART 1: Family-level combined plot
# ============================================================================

# ----------------------------- 1.1 Filter families with total cells > 1000 ----
family_totals_global <- metadata %>%
  group_by(Family) %>%
  summarise(total_cells = n(), .groups = 'drop') %>%
  filter(total_cells > 1000)

metadata_family <- metadata %>%
  semi_join(family_totals_global, by = "Family")

# ----------------------------- 1.2 Top bar data: family percentage per group --
prpdsd_totals_family <- metadata_family %>%
  group_by(PRPDSD) %>%
  summarise(total_cells = n(), .groups = 'drop')

family_prpdsd_counts <- metadata_family %>%
  group_by(PRPDSD, Family) %>%
  summarise(count = n(), .groups = 'drop')

family_prpdsd_percentage <- family_prpdsd_counts %>%
  left_join(prpdsd_totals_family, by = "PRPDSD") %>%
  mutate(percentage = count / total_cells * 100)

family_prpdsd_percentage$PRPDSD <- factor(family_prpdsd_percentage$PRPDSD,
                                          levels = c("PR", "SD", "PD"))

# ----------------------------- 1.3 Heatmap data: Ro/e per family -------------
global_total_cells_family <- sum(family_totals_global$total_cells)
global_prop_family <- family_totals_global %>%
  mutate(global_prop = total_cells / global_total_cells_family) %>%
  select(Family, global_prop)

group_obs_prop_family <- family_prpdsd_counts %>%
  left_join(prpdsd_totals_family, by = "PRPDSD") %>%
  mutate(obs_prop = count / total_cells) %>%
  select(PRPDSD, Family, obs_prop)

roe_data_family <- group_obs_prop_family %>%
  left_join(global_prop_family, by = "Family") %>%
  mutate(roe = obs_prop / global_prop)

heatmap_matrix_family <- roe_data_family %>%
  select(PRPDSD, Family, roe) %>%
  pivot_wider(names_from = PRPDSD, values_from = roe, values_fill = 0) %>%
  column_to_rownames("Family") %>%
  as.matrix()

existing_groups <- intersect(c("PR", "SD", "PD"), colnames(heatmap_matrix_family))
heatmap_matrix_family <- heatmap_matrix_family[, existing_groups, drop = FALSE]

# Hierarchical clustering of rows
if (nrow(heatmap_matrix_family) > 1) {
  row_dist <- dist(heatmap_matrix_family)
  row_clust <- hclust(row_dist, method = "complete")
  heatmap_row_order_family <- rownames(heatmap_matrix_family)[row_clust$order]
} else {
  heatmap_row_order_family <- rownames(heatmap_matrix_family)
}

# Convert to long format for ggplot
heatmap_df_family <- as.data.frame(heatmap_matrix_family) %>%
  rownames_to_column(var = "Family") %>%
  pivot_longer(cols = -Family, names_to = "PRPDSD", values_to = "roe")

heatmap_df_family$Family <- factor(heatmap_df_family$Family,
                                   levels = rev(heatmap_row_order_family))
heatmap_df_family$PRPDSD <- factor(heatmap_df_family$PRPDSD,
                                   levels = c("PR", "SD", "PD"))

# ----------------------------- 1.4 Heatmap plot ------------------------------
heatmap_gg_family <- ggplot(heatmap_df_family,
                            aes(x = PRPDSD, y = Family, fill = roe)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(
    low = "navy", mid = "white", high = "#F6944B",
    midpoint = 1,
    limits = c(0.5, 1.5),
    oob = squish,
    name = "Ro/e"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "right"
  )

# ----------------------------- 1.5 Top bar plot ------------------------------
family_prpdsd_percentage$Family <- factor(family_prpdsd_percentage$Family,
                                          levels = heatmap_row_order_family)

top_bar_family <- ggplot(family_prpdsd_percentage,
                         aes(x = PRPDSD, y = percentage, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = c("PR", "SD", "PD")) +
  labs(x = NULL, y = "Family abundance (%)", fill = "") +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(ncol = 2)) +
  theme_minimal() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(1, 1, 0.5, 1, "cm")
  )

# ----------------------------- 1.6 Right bar data: mean abundance per sample --
sample_totals_family <- metadata_family %>%
  group_by(!!sym(sample_col)) %>%
  summarise(sample_total = n(), .groups = 'drop')

sample_family_counts <- metadata_family %>%
  group_by(!!sym(sample_col), Family) %>%
  summarise(count = n(), .groups = 'drop')

sample_family_abundance <- sample_family_counts %>%
  left_join(sample_totals_family, by = sample_col) %>%
  mutate(abundance = count / sample_total * 100) %>%
  left_join(distinct(metadata_family, !!sym(sample_col), PRPDSD),
            by = sample_col)

mean_abundance_family <- sample_family_abundance %>%
  group_by(Family, PRPDSD) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = 'drop')

right_data_family <- mean_abundance_family %>%
  group_by(Family) %>%
  mutate(
    total = sum(mean_abundance),
    percentage = ifelse(total == 0, 0, mean_abundance / total * 100)
  ) %>%
  ungroup() %>%
  select(Family, PRPDSD, percentage)

right_data_family$PRPDSD <- factor(right_data_family$PRPDSD,
                                   levels = c("PD", "SD", "PR"))
right_data_family$Family <- factor(right_data_family$Family,
                                   levels = rev(heatmap_row_order_family))

right_bar_family <- ggplot(right_data_family,
                           aes(x = Family, y = percentage, fill = PRPDSD)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative proportion within family (%)", fill = "") +
  scale_fill_manual(values = c("PR" = "#90162D", "SD" = "#F7A24F", "PD" = "#C6133B")) +
  theme_minimal() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, margin = margin(b = 10)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(1, 1, 0.5, 1, "cm")
  )

# ----------------------------- 1.7 Combine and save family plot --------------
combined_family <- (top_bar_family + plot_spacer()) /
  (heatmap_gg_family + right_bar_family) +
  plot_layout(heights = c(1, 4), widths = c(4, 1.5))

ggsave(output_family_pdf, plot = combined_family,
       width = 12, height = 12, bg = "white", dpi = 300)

# ============================================================================
# PART 2: Genus-level combined plot
# ============================================================================

# ----------------------------- 2.1 Filter genera with total cells > 1000 -----
genus_totals_global <- metadata %>%
  group_by(Genus) %>%
  summarise(total_cells = n(), .groups = 'drop') %>%
  filter(total_cells > 1000)

metadata_genus <- metadata %>%
  semi_join(genus_totals_global, by = "Genus")

# ----------------------------- 2.2 Top bar data: genus percentage per group --
prpdsd_totals_genus <- metadata_genus %>%
  group_by(PRPDSD) %>%
  summarise(total_cells = n(), .groups = 'drop')

genus_prpdsd_counts <- metadata_genus %>%
  group_by(PRPDSD, Genus) %>%
  summarise(count = n(), .groups = 'drop')

genus_prpdsd_percentage <- genus_prpdsd_counts %>%
  left_join(prpdsd_totals_genus, by = "PRPDSD") %>%
  mutate(percentage = count / total_cells * 100)

genus_prpdsd_percentage$PRPDSD <- factor(genus_prpdsd_percentage$PRPDSD,
                                         levels = c("PR", "SD", "PD"))

# ----------------------------- 2.3 Heatmap data: Ro/e per genus --------------
global_total_cells_genus <- sum(genus_totals_global$total_cells)
global_prop_genus <- genus_totals_global %>%
  mutate(global_prop = total_cells / global_total_cells_genus) %>%
  select(Genus, global_prop)

group_obs_prop_genus <- genus_prpdsd_counts %>%
  left_join(prpdsd_totals_genus, by = "PRPDSD") %>%
  mutate(obs_prop = count / total_cells) %>%
  select(PRPDSD, Genus, obs_prop)

roe_data_genus <- group_obs_prop_genus %>%
  left_join(global_prop_genus, by = "Genus") %>%
  mutate(roe = obs_prop / global_prop)

heatmap_matrix_genus <- roe_data_genus %>%
  select(PRPDSD, Genus, roe) %>%
  pivot_wider(names_from = PRPDSD, values_from = roe, values_fill = 0) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

existing_groups_genus <- intersect(c("PR", "SD", "PD"),
                                   colnames(heatmap_matrix_genus))
heatmap_matrix_genus <- heatmap_matrix_genus[, existing_groups_genus, drop = FALSE]

if (nrow(heatmap_matrix_genus) > 1) {
  row_dist_genus <- dist(heatmap_matrix_genus)
  row_clust_genus <- hclust(row_dist_genus, method = "complete")
  heatmap_row_order_genus <- rownames(heatmap_matrix_genus)[row_clust_genus$order]
} else {
  heatmap_row_order_genus <- rownames(heatmap_matrix_genus)
}

heatmap_df_genus <- as.data.frame(heatmap_matrix_genus) %>%
  rownames_to_column(var = "Genus") %>%
  pivot_longer(cols = -Genus, names_to = "PRPDSD", values_to = "roe")

heatmap_df_genus$Genus <- factor(heatmap_df_genus$Genus,
                                 levels = rev(heatmap_row_order_genus))
heatmap_df_genus$PRPDSD <- factor(heatmap_df_genus$PRPDSD,
                                  levels = c("PR", "SD", "PD"))

# ----------------------------- 2.4 Heatmap plot ------------------------------
heatmap_gg_genus <- ggplot(heatmap_df_genus,
                           aes(x = PRPDSD, y = Genus, fill = roe)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(
    low = "navy", mid = "white", high = "#F6944B",
    midpoint = 1,
    limits = c(0.5, 1.5),
    oob = squish,
    name = "Ro/e"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "right"
  )

# ----------------------------- 2.5 Top bar plot ------------------------------
genus_prpdsd_percentage$Genus <- factor(genus_prpdsd_percentage$Genus,
                                        levels = heatmap_row_order_genus)

n_genus_top <- length(unique(genus_prpdsd_percentage$Genus))
top_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_genus_top)

top_bar_genus <- ggplot(genus_prpdsd_percentage,
                        aes(x = PRPDSD, y = percentage, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = c("PR", "SD", "PD")) +
  labs(x = NULL, y = "Genus abundance (%)", fill = "") +
  scale_fill_manual(values = top_colors) +
  guides(fill = guide_legend(ncol = 3)) +
  theme_minimal() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(1, 1, 0.5, 1, "cm")
  )

# ----------------------------- 2.6 Right bar data ----------------------------
sample_totals_genus <- metadata_genus %>%
  group_by(!!sym(sample_col)) %>%
  summarise(sample_total = n(), .groups = 'drop')

sample_genus_counts <- metadata_genus %>%
  group_by(!!sym(sample_col), Genus) %>%
  summarise(count = n(), .groups = 'drop')

sample_genus_abundance <- sample_genus_counts %>%
  left_join(sample_totals_genus, by = sample_col) %>%
  mutate(abundance = count / sample_total * 100) %>%
  left_join(distinct(metadata_genus, !!sym(sample_col), PRPDSD),
            by = sample_col)

mean_abundance_genus <- sample_genus_abundance %>%
  group_by(Genus, PRPDSD) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = 'drop')

right_data_genus <- mean_abundance_genus %>%
  group_by(Genus) %>%
  mutate(
    total = sum(mean_abundance),
    percentage = ifelse(total == 0, 0, mean_abundance / total * 100)
  ) %>%
  ungroup() %>%
  select(Genus, PRPDSD, percentage)

right_data_genus$PRPDSD <- factor(right_data_genus$PRPDSD,
                                  levels = c("PD", "SD", "PR"))
right_data_genus$Genus <- factor(right_data_genus$Genus,
                                 levels = rev(heatmap_row_order_genus))

right_bar_genus <- ggplot(right_data_genus,
                          aes(x = Genus, y = percentage, fill = PRPDSD)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative proportion within genus (%)", fill = "") +
  scale_fill_manual(values = c("PR" = "#90162D", "SD" = "#F7A24F", "PD" = "#C6133B")) +
  theme_minimal() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, margin = margin(b = 10)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(1, 1, 0.5, 1, "cm")
  )

# ----------------------------- 2.7 Combine and save genus plot ---------------
combined_genus <- (top_bar_genus + plot_spacer()) /
  (heatmap_gg_genus + right_bar_genus) +
  plot_layout(heights = c(1, 4), widths = c(4, 1.5))

ggsave(output_genus_pdf, plot = combined_genus,
       width = 12, height = 12, bg = "white", dpi = 300)



###############################################################################
# Figure E: Genus-level relative abundance boxplots with Kruskal-Wallis test  #
#           arranged in a 2x2 facet grid using ggh4x                          #
###############################################################################

# ----------------------------- Load required packages -----------------------
library(tidyverse)
library(ggpubr)
library(rstatix)

# Install and load ggh4x if not already installed
if (!require(ggh4x)) install.packages("ggh4x")
library(ggh4x)

# ----------------------------- Define file paths (replace with actual paths) --
metadata_path <- "path/to/metadata.rds"
output_path   <- "path/to/E_genus_boxplots_smRNA.pdf"

# ----------------------------- Load and filter metadata -----------------------
metadata <- readRDS(metadata_path)
metadata <- subset(metadata, PRPDSD != "Unknow")

# ----------------------------- Define target genera ---------------------------
target_Genera <- c("Bacteroides", "Phocaeicola", "Faecalibacterium", "Fusicatenibacter")

# ----------------------------- Compute relative abundance per sample ----------
sample_totals <- metadata %>%
  count(orig.ident) %>%
  rename(total_cells = n)

all_samples <- unique(metadata$orig.ident)

# Ensure all combinations of sample and target genus exist
all_combinations <- expand.grid(
  orig.ident = all_samples,
  Genus = target_Genera
)

genus_counts <- metadata %>%
  filter(Genus %in% target_Genera) %>%
  count(orig.ident, Genus, .drop = FALSE) %>%
  right_join(all_combinations, by = c("orig.ident", "Genus")) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  left_join(sample_totals, by = "orig.ident") %>%
  mutate(relative_abundance = n / total_cells)

# Add group information (PRPDSD)
genus_abundance <- genus_counts %>%
  left_join(distinct(metadata, orig.ident, PRPDSD), by = "orig.ident")

filtered_data <- genus_abundance %>% filter(Genus %in% target_Genera)

# ----------------------------- Kruskal-Wallis test per genus ------------------
p_values <- filtered_data %>%
  group_by(Genus) %>%
  kruskal_test(relative_abundance ~ PRPDSD) %>%
  mutate(
    p_label = case_when(
      p < 0.001 ~ "p < 0.001",
      p < 0.01  ~ paste0("p = ", round(p, 3)),
      TRUE      ~ paste0("p = ", round(p, 2))
    ),
    Genus_label = paste0(Genus, " (", p_label, ")")
  )

# Merge labels back
filtered_data <- filtered_data %>%
  left_join(p_values %>% select(Genus, Genus_label), by = "Genus")

# ----------------------------- Set factor levels and panel layout -------------
genus_label_order <- p_values %>%
  arrange(match(Genus, target_Genera)) %>%
  pull(Genus_label)

filtered_data$Genus_label <- factor(filtered_data$Genus_label, levels = genus_label_order)
filtered_data$PRPDSD <- factor(filtered_data$PRPDSD, levels = c("PR", "SD", "PD"))

# Create row/column mapping for 2x2 layout
genus_levels <- levels(filtered_data$Genus_label)
row_assign <- c(1, 1, 2, 2)   # rows: 1,1,2,2
col_assign <- c(1, 2, 1, 2)   # cols: 1,2,1,2

row_col_map <- data.frame(
  Genus_label = genus_levels,
  row = factor(paste0("Row", row_assign), levels = c("Row1", "Row2")),
  col = factor(paste0("Col", col_assign), levels = c("Col1", "Col2"))
)

filtered_data <- filtered_data %>%
  left_join(row_col_map, by = "Genus_label")

# Prepare label positions for each panel (upper-left corner)
label_data <- filtered_data %>%
  group_by(Genus_label, row, col) %>%
  summarise(
    x_pos = 1,
    y_pos = ifelse(max(relative_abundance) == 0, 1e-4, max(relative_abundance) * 1.2),
    .groups = "drop"
  )

# ----------------------------- Build boxplot ----------------------------------
p <- ggplot(filtered_data, aes(x = PRPDSD, y = relative_abundance, fill = PRPDSD)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 1) +
  geom_jitter(width = 0.2, size = 2, alpha = 1, color = "black", shape = 16) +
  # facet_grid2 from ggh4x with axes = "all" to show axis on every panel
  facet_grid2(row ~ col, scales = "free_y", axes = "all") +
  geom_text(data = label_data,
            aes(x = x_pos, y = y_pos, label = Genus_label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("PR" = "#A82E25", "PD" = "#EB7E35", "SD" = "#F5CBBF")) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.01", "0.1", "1")
  ) +
  labs(
    x = "Treatment Group",
    y = "Relative Abundance (log10 scale)",
    title = "smRNA (Genus level)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ----------------------------- Save and display -------------------------------
ggsave(output_path, plot = p, width = 6, height = 7, dpi = 300)

# =========================================================================== #
# Figure: Species-level odds ratios for PR vs PD                              #
#         Lollipop plots from univariable logistic regression                 #
# =========================================================================== #

# ====================== Paths ======================
prop_path        <- "path/to/prop_sample1_species.csv"
metagenome_path  <- "path/to/DNA_merged_abundance_table.txt"
group_path       <- "path/to/Clinical_information.rds"
output_dir       <- "path/to/F1"
output_file1     <- file.path(output_dir, "F_PRPD_lollipop.pdf")
output_file2     <- file.path(output_dir, "F_PRSDPD_Metagenome_lollipop.pdf")

# ====================== Packages ======================
library(tidyverse)
library(broom)
library(scales)

# ======================================================
# Part 1: Analysis from species proportion matrix (CSV)
# ======================================================

# ---------- Load data ----------
group <- readRDS(group_path)
prop_matrix <- read.csv(prop_path, row.names = 1, check.names = FALSE)

# ---------- Filter clinical data ----------
group <- group %>%
  filter(treatment_condition == "Pre", PRPDSD != "Unknow") %>%
  filter(PRPDSD %in% c("PR", "PD")) %>%
  select(sample, PRPDSD)

# ---------- Process species proportion matrix ----------
original_names <- colnames(prop_matrix)
colnames(prop_matrix) <- make.names(original_names, unique = TRUE)

prop_matrix <- prop_matrix[, colMeans(prop_matrix) > 0.01]
prop_matrix <- prop_matrix / rowSums(prop_matrix)

merged_data <- prop_matrix %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  inner_join(group, by = "sample")

# ---------- Prepare analysis dataset ----------
analysis_data <- merged_data %>%
  mutate(response = if_else(PRPDSD == "PR", 1, 0)) %>%
  select(-PRPDSD, -sample)

species_cols <- setdiff(names(analysis_data), "response")
analysis_data[species_cols] <- log10(analysis_data[species_cols] + 1e-5)

# ---------- Univariable logistic regression ----------
results_list <- map(species_cols, ~ {
  formula <- reformulate(.x, response = "response")
  glm(formula, data = analysis_data, family = binomial(link = "logit")) %>%
    tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    slice(1) %>%
    mutate(species = .x)
})

results_df <- bind_rows(results_list) %>%
  select(species, estimate, conf.low, conf.high, p.value) %>%
  rename(OR = estimate, CI_low = conf.low, CI_high = conf.high)

# ---------- Prepare forest-plot data ----------
forest_data <- results_df %>%
  mutate(CI = sprintf("%.2f (%.2f-%.2f)", OR, CI_low, CI_high),
         p.value = signif(p.value, 3))

forest_data <- forest_data %>%
  mutate(
    special_case = case_when(
      OR < 1e-5 ~ "extreme_low",
      OR > 1e5 ~ "extreme_high",
      is.infinite(CI_high) ~ "inf_high",
      is.na(CI_high) ~ "na_high",
      TRUE ~ "normal"
    ),
    OR_adj = case_when(
      special_case == "extreme_low" ~ 1e-5,
      special_case == "extreme_high" ~ 1e5,
      TRUE ~ OR
    ),
    CI_low_adj = case_when(
      special_case == "extreme_low" ~ CI_low,
      special_case == "extreme_high" ~ CI_low,
      TRUE ~ CI_low
    ),
    CI_high_adj = case_when(
      special_case == "inf_high" ~ 100,
      special_case == "na_high" ~ 100,
      special_case == "extreme_high" ~ 1e5,
      TRUE ~ CI_high
    ),
    display_label = case_when(
      special_case == "extreme_low" ~ sprintf("%.1e (%.1e-%.1e)", OR, CI_low, CI_high),
      special_case == "extreme_high" ~ sprintf("%.1e (%.1e-%.1e)", OR, CI_low, CI_high),
      special_case == "inf_high" ~ sprintf("%.2f (%.2f-Inf)", OR, CI_low),
      special_case == "na_high" ~ sprintf("%.1e (%.1e-NA)", OR, CI_low),
      TRUE ~ CI
    ),
    color_group = ifelse(OR_adj < 1, "below_one", "above_one")
  ) %>%
  arrange(OR_adj) %>%
  mutate(species = factor(species, levels = unique(species)))

# ---------- Lollipop plot ----------
forest_pic1 <- ggplot(forest_data, aes(x = OR_adj, y = species)) +
  geom_vline(xintercept = 1, linetype = "solid", color = "gray50") +
  geom_segment(aes(x = 1, xend = OR_adj, y = species, yend = species, color = color_group),
               linewidth = 1.2) +
  geom_point(aes(fill = color_group), shape = 21, size = 5, stroke = 0.8) +
  scale_x_log10(
    limits = c(1e-6, 1e6),
    breaks = c(0.000001, 0.0001, 0.01, 1, 100, 10000, 1000000),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_color_manual(
    values = c("below_one" = "darkgoldenrod", "above_one" = "skyblue"),
    labels = c("below_one" = "PD-associated (OR < 1)", "above_one" = "PR-associated (OR > 1)"),
    name = "Association"
  ) +
  scale_fill_manual(
    values = c("below_one" = "darkgoldenrod", "above_one" = "skyblue"),
    labels = c("below_one" = "PD-associated (OR < 1)", "above_one" = "PR-associated (OR > 1)"),
    name = "Association"
  ) +
  labs(
    x = "Odds Ratio (PR vs PD)",
    y = "",
    subtitle = "Response: PR = 1 (responders), PD = 0 (non-responders)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray85", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray85", linetype = "dashed"),
    axis.text.x = element_text(angle = 0, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.subtitle = element_text(hjust = 0.5, size = 9)
  )

# ---------- Save output ----------
ggsave(filename = output_file1, plot = forest_pic1, bg = "white",
       width = 9, height = 7, dpi = 300)

# ======================================================
# Part 2: Analysis from metagenomic abundance table
# ======================================================

# ---------- Load metagenome data ----------
meta_df <- read.delim(file = metagenome_path, check.names = FALSE, row.names = 1)

# ---------- Filter species-level annotations ----------
matrix <- meta_df %>%
  rownames_to_column(var = "Taxonomy") %>%
  filter(str_count(Taxonomy, "\\|") == 6, str_detect(Taxonomy, "\\|s__")) %>%
  mutate(Species = str_replace(Taxonomy, ".*\\|s__", "")) %>%
  select(-Taxonomy) %>%
  column_to_rownames("Species") %>%
  t() %>%
  as.data.frame()

rownames(matrix) <- gsub("_metaphlan", "", rownames(matrix))

# ---------- Retain prevalent species ----------
prop_matrix2 <- matrix %>%
  select(names(.)[colMeans(.) > 0.2])

# ---------- Reload and filter clinical data (SD-PD included this time) ----------
group2 <- readRDS(group_path) %>%
  filter(treatment_condition == "Pre", PRPDSD != "Unknow") %>%
  filter(PRPDSD %in% c("PR", "SD-PD", "PD")) %>%
  select(sample, PRPDSD)

# ---------- Merge ----------
merged_data2 <- prop_matrix2 %>%
  rownames_to_column("sample") %>%
  inner_join(group2, by = "sample")

# ---------- Prepare analysis dataset (PR vs PD only) ----------
analysis_data2 <- merged_data2 %>%
  filter(PRPDSD %in% c("PR", "PD")) %>%
  mutate(response = if_else(PRPDSD == "PR", 1, 0)) %>%
  select(-PRPDSD, -sample)

species_cols2 <- setdiff(names(analysis_data2), "response")
analysis_data2 <- analysis_data2 %>%
  rename_with(~ make.names(.x), all_of(species_cols2))
species_cols_renamed <- make.names(species_cols2)

analysis_data2[species_cols_renamed] <- log10(analysis_data2[species_cols_renamed] + 1e-5)

# ---------- Univariable logistic regression ----------
results_list2 <- map(species_cols_renamed, ~ {
  species_var <- .x
  pr_data <- analysis_data2[analysis_data2$response == 1, species_var]
  pd_data <- analysis_data2[analysis_data2$response == 0, species_var]
  
  if (sd(pr_data) == 0 || sd(pd_data) == 0) {
    return(tibble(
      term = species_var,
      estimate = ifelse(mean(pr_data) > mean(pd_data), 1e10, 1e-10),
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      species = species_var
    ))
  }
  
  formula <- reformulate(species_var, response = "response")
  tryCatch({
    model <- glm(formula, data = analysis_data2, family = binomial(link = "logit"))
    if (!model$converged) stop("Non-converged model")
    result <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
    result <- result %>% filter(term == species_var)
    result$species <- species_var
    return(result)
  }, error = function(e) {
    return(tibble(
      term = species_var,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      species = species_var
    ))
  })
})

results_df2 <- bind_rows(results_list2) %>%
  select(species, estimate, conf.low, conf.high, p.value) %>%
  rename(OR = estimate, CI_low = conf.low, CI_high = conf.high) %>%
  mutate(species_original = species)

results_df_filtered2 <- results_df2 %>%
  filter(!is.na(OR), OR <= 1e5, OR >= 1e-5)

# ---------- Prepare forest-plot data ----------
forest_data2 <- results_df_filtered2 %>%
  mutate(CI = sprintf("%.2f (%.2f-%.2f)", OR, CI_low, CI_high),
         p.value = signif(p.value, 3)) %>%
  select(species_original, OR, CI_low, CI_high, CI, p.value)

if (nrow(forest_data2) == 0) stop("No valid data for forest plot.")

# ---------- Lollipop plot ----------
forest_data2 <- forest_data2 %>%
  mutate(
    OR_adj = OR,
    color_group = case_when(
      OR_adj > 1 ~ "above_one",
      OR_adj < 1 ~ "below_one",
      TRUE ~ "equal_one"
    )
  ) %>%
  filter(OR_adj >= 0.1 & OR_adj <= 5) %>%
  arrange(OR_adj) %>%
  mutate(species_original = factor(species_original, levels = species_original))

forest_pic2 <- ggplot(forest_data2, aes(x = OR_adj, y = species_original)) +
  geom_vline(xintercept = 1, linetype = "solid", color = "gray50") +
  geom_segment(aes(x = 1, xend = OR_adj, y = species_original, yend = species_original,
                   color = color_group),
               linewidth = 1.2) +
  geom_point(aes(fill = color_group), shape = 21, size = 5, stroke = 0.8) +
  scale_x_log10(
    limits = c(0.1, 5),
    breaks = c(0.1, 0.5, 1, 5, 10),
    labels = function(x) format(x, scientific = FALSE, trim = TRUE)
  ) +
  scale_color_manual(
    values = c("below_one" = "darkgoldenrod", "above_one" = "skyblue"),
    labels = c("below_one" = "PD-associated (OR < 1)",
               "above_one" = "PR-associated (OR > 1)"),
    name = "Association"
  ) +
  scale_fill_manual(
    values = c("below_one" = "darkgoldenrod", "above_one" = "skyblue"),
    labels = c("below_one" = "PD-associated (OR < 1)",
               "above_one" = "PR-associated (OR > 1)"),
    name = "Association"
  ) +
  labs(
    x = "Odds Ratio (PR vs PD)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray85", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray85", linetype = "dashed"),
    axis.text.x = element_text(angle = 0, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

# ---------- Save output ----------
ggsave(filename = output_file2, plot = forest_pic2, bg = "white",
       width = 9, height = 7, dpi = 300)



# =========================================================================== #
# Figure G: Single-phylum relative abundance boxplots across Pre/Post/Post2     #
#         with Kruskal–Wallis test (smRNA or Metagenome)                      #
# =========================================================================== #

# ====================== Paths ======================
rds_path          <- "path/to/Clinical_information"
metagenome_path   <- "path/to/DNA_merged_abundance_table.txt"
output_dir        <- "path/to/F1"

# ====================== Parameters (modify as needed) ======================
PRPDSD        <- "PR"              # "PR", "SD", "PD"
target_phylum <- "Firmicutes"      # e.g. "Bacteroidota", "Proteobacteria"
data_to_plot  <- "smRNA"           # "smRNA" or "Metagenome"

# ====================== Packages ======================
library(tidyverse)
library(ggpubr)
library(rstatix)

# ====================== Load and filter metadata (smRNA source) ======================
metadata <- readRDS(file.path(rds_path, "metadata.rds")) %>%
  filter(treatment_condition != "Unknow",
         name != "PanZhengfu",
         name != "WangWei",
         PRPDSD == PRPDSD)

# Prepare sample-level info for merging
metadata_samples <- metadata %>%
  distinct(orig.ident, treatment_condition, name) %>%
  filter(treatment_condition %in% c("Pre", "Post", "Post2"))

# Compute relative abundance of the target phylum from cell counts (smRNA)
phylum_counts_meta <- metadata %>%
  count(orig.ident, Phylum) %>%
  group_by(orig.ident) %>%
  mutate(total_cells = sum(n),
         relative_abundance = n / total_cells) %>%
  ungroup() %>%
  left_join(metadata_samples, by = "orig.ident") %>%
  filter(Phylum == target_phylum,
         !is.na(treatment_condition))

# ====================== Load metagenomic phylum-level abundance ======================
meta_df <- read.delim(file = metagenome_path, check.names = FALSE, row.names = 1)

# Keep phylum-level annotations (rank depth 2, pattern "|p__")
df2 <- meta_df %>%
  rownames_to_column(var = "Taxonomy") %>%
  filter(str_count(Taxonomy, "\\|") == 1, str_detect(Taxonomy, "\\|p__")) %>%
  mutate(Phylum = str_replace(Taxonomy, ".*\\|p__", "")) %>%
  select(-Taxonomy) %>%
  column_to_rownames("Phylum") %>%
  t() %>%
  as.data.frame() %>%
  select(where(~ any(. != 0)))

# Clean sample names
rownames(df2) <- gsub("_metaphlan", "", rownames(df2))

# Harmonize phylum names
df2 <- df2 %>%
  rename(
    Actinobacteriota   = Actinobacteria,
    Bacteroidota       = Bacteroidetes,
    Methanobacteriota  = Euryarchaeota,
    Lentisphaerota     = Lentisphaerae,
    Verrucomicrobiota  = Verrucomicrobia,
    Synergistota       = Synergistetes,
    Fusobacteriota     = Fusobacteria
  ) %>%
  select(order(colnames(.)))

# Convert to long format and filter target phylum
df2_long <- df2 %>%
  rownames_to_column(var = "orig.ident") %>%
  filter(orig.ident %in% metadata_samples$orig.ident) %>%
  pivot_longer(cols = -orig.ident, names_to = "Phylum", values_to = "relative_abundance") %>%
  left_join(metadata_samples, by = "orig.ident") %>%
  filter(Phylum == target_phylum,
         !is.na(treatment_condition))

# ====================== Kruskal–Wallis tests ======================
# smRNA
p_meta <- phylum_counts_meta %>%
  kruskal_test(relative_abundance ~ treatment_condition) %>%
  mutate(
    p_label = case_when(
      p < 0.001 ~ "p < 0.001",
      p < 0.01  ~ paste0("p = ", round(p, 3)),
      TRUE      ~ paste0("p = ", round(p, 2))
    ),
    Phylum_label = paste0("smRNA", " (", p_label, ")")
  )

# Metagenome
p_df2 <- df2_long %>%
  kruskal_test(relative_abundance ~ treatment_condition) %>%
  mutate(
    p_label = case_when(
      p < 0.001 ~ "p < 0.001",
      p < 0.01  ~ paste0("p = ", round(p, 3)),
      TRUE      ~ paste0("p = ", round(p, 2))
    ),
    Phylum_label = paste0("Metagenome", " (", p_label, ")")
  )

# Attach labels (not shown in the plot but available for future faceting)
phylum_counts_meta <- phylum_counts_meta %>%
  mutate(Phylum_label = p_meta$Phylum_label[1],
         data_source  = "smRNA")

df2_long <- df2_long %>%
  mutate(Phylum_label = p_df2$Phylum_label[1],
         data_source  = "Metagenome")

# ====================== Select data source to plot ======================
if (data_to_plot == "smRNA") {
  plot_data <- phylum_counts_meta
} else if (data_to_plot == "Metagenome") {
  plot_data <- df2_long
} else {
  stop("data_to_plot must be either 'smRNA' or 'Metagenome'")
}

# Ensure correct timepoint order
plot_data$treatment_condition <- factor(plot_data$treatment_condition,
                                        levels = c("Pre", "Post", "Post2"))

# Add per-patient jitter for consistent point positions across time
plot_data <- plot_data %>%
  group_by(name) %>%
  mutate(jitter_offset = runif(1, -0.1, 0.1)) %>%
  ungroup() %>%
  mutate(x_numeric = as.numeric(treatment_condition) + jitter_offset)

# ====================== Boxplot ======================
p <- ggplot(plot_data, aes(x = treatment_condition, y = relative_abundance,
                           fill = treatment_condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 1) +
  geom_line(aes(x = x_numeric, group = name),
            alpha = 0.3, linewidth = 0.5, color = "black") +
  geom_point(aes(x = x_numeric),
             size = 3, alpha = 1, color = "black", shape = 16) +
  scale_fill_manual(values = c("Pre"  = "#A82E25",
                               "Post" = "#EB7E35",
                               "Post2"= "#F5CBBF")) +
  scale_y_continuous(
    breaks = c(0.0, 0.3, 0.6, 0.9),
    labels = c("0.0", "0.3", "0.6", "0.9")
  ) +
  labs(x = "",
       y = "Relative Abundance",
       title = paste(target_phylum, " (", PRPDSD, ")")) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text          = element_text(face = "bold", size = 15),
    axis.title          = element_text(face = "bold", size = 15),
    axis.text.x         = element_text(angle = 0, hjust = 0.5),
    plot.title          = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position     = "none",
    panel.border        = element_blank(),
    axis.line           = element_line(color = "black", size = 0.5),
    axis.ticks          = element_line(color = "black", size = 0.5),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank()
  )

# ====================== Save output ======================
ggsave(file.path(output_dir,
                 paste0("G_", data_to_plot, "_", target_phylum, "_", PRPDSD, ".pdf")),
       plot = p, width = 3, height = 4, dpi = 300)


# =========================================================================== #
# Figure H: Family-level and Genus-level Ro/e heatmaps with stacked bar plots   #
#         for treatment conditions Pre, Post, Post2                           #
# =========================================================================== #

# ====================== Paths ======================
metadata_path   <- "path/to/metadata.rds"
output_dir      <- "path/to/F1"
output_family   <- file.path(output_dir, "H_combined_family_treatment_roe.pdf")
output_genus    <- file.path(output_dir, "H_combined_genus_treatment_roe.pdf")

# ====================== Packages ======================
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(scales)

# ====================== Load metadata ======================
metadata <- readRDS(metadata_path) %>%
  filter(treatment_condition != "Unknow")

# ======================================================
# Part 1: Family-level analysis
# ======================================================

# ---------- Filter families with total cells > 1000 ----------
family_totals_global <- metadata %>%
  group_by(Family) %>%
  summarise(total_cells = n(), .groups = 'drop') %>%
  filter(total_cells > 1000)

metadata_fam <- metadata %>%
  semi_join(family_totals_global, by = "Family")

# ---------- Top bar: family composition per condition ----------
condition_totals_fam <- metadata_fam %>%
  group_by(treatment_condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

family_condition_counts <- metadata_fam %>%
  group_by(treatment_condition, Family) %>%
  summarise(count = n(), .groups = 'drop')

family_condition_percentage <- family_condition_counts %>%
  left_join(condition_totals_fam, by = "treatment_condition") %>%
  mutate(percentage = count / total_cells * 100) %>%
  mutate(treatment_condition = factor(treatment_condition, 
                                      levels = c("Pre", "Post", "Post2")))

# ---------- Right bar: intra-family composition by condition ----------
sample_col <- "orig.ident"

sample_totals_fam <- metadata_fam %>%
  group_by(!!sym(sample_col)) %>%
  summarise(sample_total = n(), .groups = 'drop')

sample_family_counts <- metadata_fam %>%
  group_by(!!sym(sample_col), Family) %>%
  summarise(count = n(), .groups = 'drop')

sample_family_abundance <- sample_family_counts %>%
  left_join(sample_totals_fam, by = sample_col) %>%
  mutate(abundance = count / sample_total * 100) %>%
  left_join(metadata_fam %>% 
              select(!!sym(sample_col), treatment_condition) %>% 
              distinct(),
            by = sample_col)

mean_abundance_fam <- sample_family_abundance %>%
  group_by(Family, treatment_condition) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = 'drop')

right_data_fam <- mean_abundance_fam %>%
  group_by(Family) %>%
  mutate(total = sum(mean_abundance),
         percentage = ifelse(total == 0, 0, mean_abundance / total * 100)) %>%
  ungroup() %>%
  select(Family, treatment_condition, percentage) %>%
  mutate(treatment_condition = factor(treatment_condition,
                                      levels = c("Post2", "Post", "Pre")))

# ---------- Heatmap: Ro/e values ----------
global_total_cells_fam <- sum(family_totals_global$total_cells)
global_prop_fam <- family_totals_global %>%
  mutate(global_prop = total_cells / global_total_cells_fam) %>%
  select(Family, global_prop)

group_obs_prop_fam <- family_condition_counts %>%
  left_join(condition_totals_fam, by = "treatment_condition") %>%
  mutate(obs_prop = count / total_cells) %>%
  select(treatment_condition, Family, obs_prop)

roe_data_fam <- group_obs_prop_fam %>%
  left_join(global_prop_fam, by = "Family") %>%
  mutate(roe = obs_prop / global_prop)

heatmap_matrix_fam <- roe_data_fam %>%
  select(treatment_condition, Family, roe) %>%
  pivot_wider(names_from = treatment_condition, values_from = roe, values_fill = 0) %>%
  column_to_rownames("Family") %>%
  as.matrix()

existing_groups <- intersect(c("Pre", "Post", "Post2"), colnames(heatmap_matrix_fam))
heatmap_matrix_fam <- heatmap_matrix_fam[, existing_groups, drop = FALSE]

if (nrow(heatmap_matrix_fam) > 1) {
  row_dist <- dist(heatmap_matrix_fam)
  row_clust <- hclust(row_dist, method = "complete")
  heatmap_row_order_fam <- rownames(heatmap_matrix_fam)[row_clust$order]
} else {
  heatmap_row_order_fam <- rownames(heatmap_matrix_fam)
}

heatmap_df_fam <- as.data.frame(heatmap_matrix_fam) %>%
  rownames_to_column(var = "Family") %>%
  pivot_longer(cols = -Family, names_to = "treatment_condition", values_to = "roe") %>%
  mutate(Family = factor(Family, levels = rev(heatmap_row_order_fam)),
         treatment_condition = factor(treatment_condition, 
                                      levels = c("Pre", "Post", "Post2")))

# ---------- Family heatmap ----------
color_low <- 0.5
color_high <- 1.5

heatmap_fam <- ggplot(heatmap_df_fam, aes(x = treatment_condition, y = Family, fill = roe)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(
    low = "navy", mid = "white", high = "#F6944B",
    midpoint = 1,
    limits = c(color_low, color_high),
    oob = squish,
    name = "Ro/e"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(0, 0, 10, 0),
    legend.position = "right"
  )

# ---------- Top bar (family composition) ----------
family_condition_percentage$Family <- factor(family_condition_percentage$Family,
                                             levels = heatmap_row_order_fam)
n_family <- length(unique(family_condition_percentage$Family))
top_colors_fam <- colorRampPalette(brewer.pal(12, "Paired"))(n_family)

top_bar_fam <- ggplot(family_condition_percentage,
                       aes(x = treatment_condition, y = percentage, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = c("Pre", "Post", "Post2")) +
  labs(x = NULL, y = "Family abundance (%)", fill = "") +
  scale_fill_manual(values = top_colors_fam) +
  guides(fill = guide_legend(ncol = 2)) +
  theme_minimal() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(1, 1, 0.5, 1, "cm")
  )

# ---------- Right bar (intra-family proportion) ----------
right_data_fam$Family <- factor(right_data_fam$Family,
                                levels = rev(heatmap_row_order_fam))

right_bar_fam <- ggplot(right_data_fam,
                         aes(x = Family, y = percentage, fill = treatment_condition)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative proportion within family (%)", fill = "") +
  scale_fill_manual(values = c("Pre" = "#90162D", "Post" = "#F7A24F", "Post2" = "#C6133B")) +
  theme_minimal() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, margin = margin(b = 10)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(1, 1, 0.5, 1, "cm")
  )

# ---------- Assemble and save family plot ----------
combined_family <- (top_bar_fam + plot_spacer()) / (heatmap_fam + right_bar_fam) +
  plot_layout(heights = c(1, 4), widths = c(4, 1.5))

ggsave(output_family, plot = combined_family,
       width = 12, height = 12, bg = "white", dpi = 300)

# ======================================================
# Part 2: Genus-level analysis (same workflow)
# ======================================================

# ---------- Filter genera with total cells > 1000 ----------
genus_totals_global <- metadata %>%
  group_by(Genus) %>%
  summarise(total_cells = n(), .groups = 'drop') %>%
  filter(total_cells > 1000)

metadata_gen <- metadata %>%
  semi_join(genus_totals_global, by = "Genus")

# ---------- Top bar: genus composition ----------
condition_totals_gen <- metadata_gen %>%
  group_by(treatment_condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

genus_condition_counts <- metadata_gen %>%
  group_by(treatment_condition, Genus) %>%
  summarise(count = n(), .groups = 'drop')

genus_condition_percentage <- genus_condition_counts %>%
  left_join(condition_totals_gen, by = "treatment_condition") %>%
  mutate(percentage = count / total_cells * 100) %>%
  mutate(treatment_condition = factor(treatment_condition,
                                      levels = c("Pre", "Post", "Post2")))

# ---------- Right bar: intra-genus composition ----------
sample_totals_gen <- metadata_gen %>%
  group_by(!!sym(sample_col)) %>%
  summarise(sample_total = n(), .groups = 'drop')

sample_genus_counts <- metadata_gen %>%
  group_by(!!sym(sample_col), Genus) %>%
  summarise(count = n(), .groups = 'drop')

sample_genus_abundance <- sample_genus_counts %>%
  left_join(sample_totals_gen, by = sample_col) %>%
  mutate(abundance = count / sample_total * 100) %>%
  left_join(metadata_gen %>%
              select(!!sym(sample_col), treatment_condition) %>%
              distinct(),
            by = sample_col)

mean_abundance_gen <- sample_genus_abundance %>%
  group_by(Genus, treatment_condition) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = 'drop')

right_data_gen <- mean_abundance_gen %>%
  group_by(Genus) %>%
  mutate(total = sum(mean_abundance),
         percentage = ifelse(total == 0, 0, mean_abundance / total * 100)) %>%
  ungroup() %>%
  select(Genus, treatment_condition, percentage) %>%
  mutate(treatment_condition = factor(treatment_condition,
                                      levels = c("Post2", "Post", "Pre")))

# ---------- Heatmap: Ro/e ----------
global_total_cells_gen <- sum(genus_totals_global$total_cells)
global_prop_gen <- genus_totals_global %>%
  mutate(global_prop = total_cells / global_total_cells_gen) %>%
  select(Genus, global_prop)

group_obs_prop_gen <- genus_condition_counts %>%
  left_join(condition_totals_gen, by = "treatment_condition") %>%
  mutate(obs_prop = count / total_cells) %>%
  select(treatment_condition, Genus, obs_prop)

roe_data_gen <- group_obs_prop_gen %>%
  left_join(global_prop_gen, by = "Genus") %>%
  mutate(roe = obs_prop / global_prop)

heatmap_matrix_gen <- roe_data_gen %>%
  select(treatment_condition, Genus, roe) %>%
  pivot_wider(names_from = treatment_condition, values_from = roe, values_fill = 0) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

existing_groups_gen <- intersect(c("Pre", "Post", "Post2"), colnames(heatmap_matrix_gen))
heatmap_matrix_gen <- heatmap_matrix_gen[, existing_groups_gen, drop = FALSE]

if (nrow(heatmap_matrix_gen) > 1) {
  row_dist <- dist(heatmap_matrix_gen)
  row_clust <- hclust(row_dist, method = "complete")
  heatmap_row_order_gen <- rownames(heatmap_matrix_gen)[row_clust$order]
} else {
  heatmap_row_order_gen <- rownames(heatmap_matrix_gen)
}

heatmap_df_gen <- as.data.frame(heatmap_matrix_gen) %>%
  rownames_to_column(var = "Genus") %>%
  pivot_longer(cols = -Genus, names_to = "treatment_condition", values_to = "roe") %>%
  mutate(Genus = factor(Genus, levels = rev(heatmap_row_order_gen)),
         treatment_condition = factor(treatment_condition,
                                      levels = c("Pre", "Post", "Post2")))

# ---------- Genus heatmap ----------
heatmap_gen <- ggplot(heatmap_df_gen, aes(x = treatment_condition, y = Genus, fill = roe)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(
    low = "navy", mid = "white", high = "#F6944B",
    midpoint = 1,
    limits = c(color_low, color_high),
    oob = squish,
    name = "Ro/e"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(0, 0, 10, 0),
    legend.position = "right"
  )

# ---------- Top bar (genus composition) ----------
genus_condition_percentage$Genus <- factor(genus_condition_percentage$Genus,
                                           levels = heatmap_row_order_gen)
n_genus <- length(unique(genus_condition_percentage$Genus))
top_colors_gen <- colorRampPalette(brewer.pal(12, "Paired"))(n_genus)

top_bar_gen <- ggplot(genus_condition_percentage,
                       aes(x = treatment_condition, y = percentage, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = c("Pre", "Post", "Post2")) +
  labs(x = NULL, y = "Genus abundance (%)", fill = "") +
  scale_fill_manual(values = top_colors_gen) +
  guides(fill = guide_legend(ncol = 3)) +
  theme_minimal() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(1, 1, 0.5, 1, "cm")
  )

# ---------- Right bar (intra-genus proportion) ----------
right_data_gen$Genus <- factor(right_data_gen$Genus,
                               levels = rev(heatmap_row_order_gen))

right_bar_gen <- ggplot(right_data_gen,
                         aes(x = Genus, y = percentage, fill = treatment_condition)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative proportion within genus (%)", fill = "") +
  scale_fill_manual(values = c("Pre" = "#90162D", "Post" = "#F7A24F", "Post2" = "#C6133B")) +
  theme_minimal() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, margin = margin(b = 10)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(1, 1, 0.5, 1, "cm")
  )

# ---------- Assemble and save genus plot ----------
combined_genus <- (top_bar_gen + plot_spacer()) / (heatmap_gen + right_bar_gen) +
  plot_layout(heights = c(1, 4), widths = c(4, 1.5))

ggsave(output_genus, plot = combined_genus,
       width = 12, height = 12, bg = "white", dpi = 300)

# =========================================================================== #
# Figure I: Genus-level marker gene dot plot with expression heatmap,           #
#         log2FC, and Cox hazard ratios (PR vs PD)                           #
# =========================================================================== #

# ====================== Paths ======================
seurat_rds      <- "path/to/01merge.rds"
metadata_rds    <- "path/to/metadata.rds"
marker_file     <- "path/to/I_marker_inGenus_gt1000.txt"
output_dir      <- "path/to/F1"
output_pdf      <- file.path(output_dir, "I_DotPlot_publication_genus_heatmap.pdf")

# ====================== Packages ======================
library(Seurat)
library(tidyverse)
library(survival)
library(scales)
library(RColorBrewer)

# ====================== Load data ======================
scRNA <- readRDS(seurat_rds)
metadata <- readRDS(metadata_rds)
scRNA@meta.data <- metadata
rownames(scRNA@meta.data) <- scRNA@meta.data$barcode

DefaultAssay(scRNA) <- "RNA"

# ====================== Filter genera with >1000 cells ======================
genus_counts <- table(scRNA@meta.data$Genus)
top_genus <- names(genus_counts[genus_counts > 1000])
scRNA_subset <- subset(scRNA, subset = Genus %in% top_genus)

# ====================== Find marker genes per genus ======================
orig_idents <- unique(scRNA_subset@meta.data$Genus)
all_markers <- data.frame()

for (ident in orig_idents) {
  ident_markers <- FindMarkers(
    scRNA_subset,
    ident.1 = ident,
    group.by = 'Genus',
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0,
    return.thresh = 0.05
  )
  ident_markers$orig.ident <- ident
  all_markers <- rbind(all_markers, ident_markers)
}

all_markers <- all_markers %>%
  rownames_to_column('gene') %>%
  select(gene, everything())

cluster_gene <- all_markers %>%
  group_by(orig.ident) %>%
  rename(cluster = orig.ident)

write.table(cluster_gene, file = marker_file, sep = "\t", quote = FALSE, row.names = FALSE)

# ====================== Subset to functional gene list ======================
func_genes_vec <- c(
  "sodB","eno","gnd","tpx","zwf","tufA","por","gapA",
  "pckA","gapA","fda","frdA","infB","uvrA","deaD",
  "cpxP","rpoH","pckA","dfx","nagB-1","nif","fusA",
  "gdhA","dnaK","groL","clpB1","htpG","gdh","gcdA",
  "bcd","gctA","gap","mglB-2","ppdk","susC-10","susD",
  "susF","susE","manZ-2","lagD","nagE","agaC","bhsa-2",
  "aceE","bla-2","bluF-2","tufB","por-1","por-2","glgC-1",
  "mdh","mutB-2","pyk","eno","tpx","zwf","gnd","trxA","hag",
  "bcd","sugC","flgK","pckA","dxs-1","aspA","susC-24",
  "por-2","grdB-4","grdE-5","por-1","sigH-3","cshA",
  "cspL-2","cutC-3","tpa","ald2","dsvA","katA","por",
  "tufA","mdh-2","pnp-2","gapA","por-2","oppA-1","tuf",
  "tpl","alst-5","por-1","clpB","dnak-2","hutu-1","gdh-4",
  "fhs","purR-6","tuf1","por-2","glgP","thlA","dnak2","rpoC",
  "rpoB","ftsZ-2","aprA-1","aprB-1","dhaT","sat","katA",
  "sodB","grdE","grdB","gdh","trxB-1","glpK-3","glpF",
  "lhgo-2","hag-1","gdh","gcdA","bcd","lacE","cbgA-1",
  "lacF-1","mleN-2","thlA-1","por-4","mdh","mutB-2","pyk",
  "gap","tuf","sspC2","spoIVA","spoVB","mutB-1","scpA-2",
  "cat1-1","gsiB","braC","glgC-1","fldC","clpB","dnak",
  "groL","htpG","hbd","thlA","crt","cutC-2","por","hadB",
  "carE","carD-1","gdh","cutC-3","dsva","gadB","pccB-2",
  "sodB","acnA","prdA-4","adhE-2","por","aprA","mop",
  "dhaT-2","dsVB","glpk-1","f1aB","cbiQ","dapL","speF-1","potE-3","lysS","por"
)

df_func <- data.frame(gene = func_genes_vec, stringsAsFactors = FALSE)

marker <- read.delim(marker_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
marker_common <- marker[marker$gene %in% df_func$gene, ]

# Keep top 300 markers by avg_log2FC per genus
marker_top2 <- marker_common %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 300) %>%
  ungroup()

# Remove genes that are difficult to interpret
marker_top2 <- marker_top2[!(marker_top2$gene %in% c("cspL-2", "sigH-3")), ]

# ====================== Prepare Seurat object for plots ======================
seurat_obj <- scRNA

# Correct PFS_situation for specific samples
seurat_obj@meta.data$PFS_situation[seurat_obj@meta.data$orig.ident %in% c("A1536", "A1593", "A1787")] <- 0

# Subset to Pre treatment only and remove unknown PRPDSD
seurat_obj <- subset(seurat_obj, subset = treatment_condition %in% c("Post", "Post2"), invert = TRUE)
seurat_obj <- subset(seurat_obj, subset = PRPDSD %in% c("Unknow"), invert = TRUE)

# ====================== Gene list ======================
marker_genes <- unique(marker_top2$gene)
genes_to_plot <- intersect(marker_genes, rownames(seurat_obj))

# ====================== Plot 1: Dot plot (genus x gene) ======================
genus_info <- seurat_obj$Genus
expr_matrix <- GetAssayData(seurat_obj, slot = "data")

genus_to_include <- names(genus_counts[genus_counts >= 1000])

dotplot_data <- data.frame()
for (genus in genus_to_include) {
  cell_idx <- which(genus_info == genus)
  for (gene in genes_to_plot) {
    expr_values <- expr_matrix[gene, ]
    avg_expr <- mean(expr_values[cell_idx])
    pct_expressed <- sum(expr_values[cell_idx] > 0) / length(cell_idx) * 100
    dotplot_data <- rbind(dotplot_data, data.frame(
      Genus = genus,
      Gene = gene,
      AvgExpr = avg_expr,
      PctExpr = pct_expressed
    ))
  }
}

# Order genera by phylum
genus_phylum <- seurat_obj@meta.data %>%
  as.data.frame() %>%
  filter(Genus %in% genus_to_include) %>%
  select(Genus, Phylum) %>%
  group_by(Genus) %>%
  summarise(Phylum = names(sort(table(Phylum), decreasing = TRUE))[1], .groups = 'drop')

genus_order <- genus_phylum %>%
  arrange(Phylum, Genus) %>%
  pull(Genus)

missing_genus <- setdiff(genus_to_include, genus_order)
genus_order <- c(genus_order, missing_genus)

dotplot_data$Genus <- factor(dotplot_data$Genus, levels = genus_order)

p_dot <- ggplot(dotplot_data, aes(x = Gene, y = Genus, color = AvgExpr, size = PctExpr)) +
  geom_point() +
  scale_color_gradientn(
    colors = c("white", "#5D90BA", "#4a1486"),
    limits = c(0, 0.1), oob = squish, name = "Average Expression"
  ) +
  scale_size_continuous(limits = c(0, 100), range = c(5, 13), name = "% Expressed") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90", size = 0.25),
    panel.grid.minor = element_line(color = "gray90", size = 0.25)
  ) +
  labs(x = "", y = "")

# ====================== Prepare survival metadata ======================
meta <- seurat_obj@meta.data %>%
  filter(!orig.ident %in% c("A1741", "A1869")) %>%
  mutate(
    time_numeric = as.numeric(time),
    event = case_when(
      PFS_situation == "1"      ~ 1,
      PFS_situation == "0"      ~ 0,
      PFS_situation == "Unknow" ~ 0
    )
  ) %>%
  filter(!is.na(time_numeric) & time_numeric > 0)

sample_surv <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    time  = unique(time_numeric)[1],
    event = unique(event)[1],
    PRPDSD = unique(PRPDSD)[1],
    .groups = 'drop'
  ) %>%
  filter(!is.na(time) & !is.na(event) & !is.na(PRPDSD),
         PRPDSD %in% c("PR", "SD", "PD"))

# Average expression per sample
expr_per_sample <- AverageExpression(seurat_obj,
                                     features = genes_to_plot,
                                     group.by = "orig.ident",
                                     slot = "data")[[1]]
expr_per_sample <- as.data.frame(t(expr_per_sample)) %>%
  rownames_to_column(var = "orig.ident")

cox_data <- expr_per_sample %>%
  left_join(sample_surv, by = "orig.ident") %>%
  filter(!is.na(time) & !is.na(event) & !is.na(PRPDSD))

# ====================== Plot 2: Expression heatmap by response ======================
sample_ids <- cox_data$orig.ident
cells_in_samples <- rownames(meta)[meta$orig.ident %in% sample_ids]

expr_sub <- expr_matrix[genes_to_plot, cells_in_samples, drop = FALSE]
group <- as.character(meta[cells_in_samples, "PRPDSD"])

response_levels <- unique(cox_data$PRPDSD)

avg_expr_list <- list()
for (resp in response_levels) {
  cells_in_resp <- which(group == resp)
  if (length(cells_in_resp) == 0) {
    avg_expr_list[[resp]] <- rep(NA, length(genes_to_plot))
  } else {
    avg_expr_list[[resp]] <- rowMeans(expr_sub[, cells_in_resp, drop = FALSE])
  }
}
avg_expr_matrix <- do.call(rbind, avg_expr_list)
colnames(avg_expr_matrix) <- genes_to_plot

# Z-score per gene
avg_expr_matrix_scaled <- apply(avg_expr_matrix, 2, function(x) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  x_no_na <- x[!is.na(x)]
  if (length(x_no_na) > 1) {
    (x - mean(x_no_na)) / sd(x_no_na)
  } else {
    rep(0, length(x))
  }
})
rownames(avg_expr_matrix_scaled) <- rownames(avg_expr_matrix)

avg_expr_long <- as.data.frame(avg_expr_matrix_scaled) %>%
  rownames_to_column(var = "Response") %>%
  pivot_longer(cols = -Response, names_to = "Gene", values_to = "ScaledExpr") %>%
  mutate(Response = factor(Response, levels = response_levels))

expr_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

p_expr_heat <- ggplot(avg_expr_long, aes(x = Gene, y = Response, fill = ScaledExpr)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradientn(
    colors = expr_colors,
    limits = c(-1, 1),
    oob = squish,
    name = "Scaled Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ====================== Plot 3: log2FC (PR/PD) ======================
if ("PR" %in% response_levels && "PD" %in% response_levels) {
  pr_mean <- avg_expr_matrix["PR", genes_to_plot, drop = TRUE]
  pd_mean <- avg_expr_matrix["PD", genes_to_plot, drop = TRUE]
  log2fc <- log2((pr_mean + 1e-6) / (pd_mean + 1e-6))
  genes_sorted <- genes_to_plot[order(log2fc, decreasing = TRUE)]
  log2fc_df <- data.frame(
    Group = "PR vs PD",
    Gene = genes_to_plot,
    log2FC = log2fc,
    stringsAsFactors = FALSE
  ) %>%
    mutate(Group = factor(Group, levels = "PR vs PD"),
           label = sprintf("%.2f", log2FC))
} else {
  genes_sorted <- genes_to_plot
  log2fc_df <- data.frame(
    Group = "PR vs PD", Gene = genes_to_plot,
    log2FC = NA, label = "NA", stringsAsFactors = FALSE
  )
}

p_log2fc <- ggplot(log2fc_df, aes(x = Gene, y = Group, fill = log2FC)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_distiller(
    palette = "BrBG", direction = -1,
    limits = c(-2, 2),
    oob = squish,
    name = "log2FC (PR/PD)",
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ====================== Plot 4: Cox HR (PR vs PD only) ======================
cox_data_hr <- cox_data %>% filter(PRPDSD %in% c("PR", "PD"))
genes <- genes_to_plot
hr_vec <- setNames(rep(NA, length(genes)), genes)
p_vec <- setNames(rep(NA, length(genes)), genes)

for (gene in genes) {
  form <- as.formula(paste("Surv(time, event) ~ `", gene, "`", sep = ""))
  cox_fit <- tryCatch(coxph(form, data = cox_data_hr), error = function(e) NULL)
  if (!is.null(cox_fit)) {
    s <- summary(cox_fit)
    hr_vec[gene] <- exp(s$coefficients[1, "coef"])
    p_vec[gene]  <- s$coefficients[1, "Pr(>|z|)"]
  }
}

star_vec <- rep("", length(genes))
names(star_vec) <- genes
star_vec[p_vec < 0.001 & !is.na(p_vec)] <- "***"
star_vec[p_vec < 0.01  & !is.na(p_vec) & p_vec >= 0.001] <- "**"
star_vec[p_vec < 0.05  & !is.na(p_vec) & p_vec >= 0.01] <- "*"

hr_long <- data.frame(
  Group = "Overall",
  Gene = genes,
  HR = hr_vec,
  Star = star_vec,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Group = factor(Group, levels = "Overall"),
    HR_fill = ifelse(is.infinite(HR), 2, HR),
    label = ifelse(is.infinite(HR), "Inf",
                   ifelse(is.na(HR), "NA", Star))
  )

p_heat <- ggplot(hr_long, aes(x = Gene, y = Group, fill = HR_fill)) +
  geom_tile(color = "white", size = 0.1) +
  geom_text(aes(label = label), size = 4, vjust = 0.8) +
  scale_fill_distiller(
    palette = "PiYG",
    direction = 1,
    limits = c(0, 2),
    oob = squish,
    name = "Hazard Ratio",
    na.value = "grey100"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ====================== Align factor levels across plots ======================
dotplot_data$Gene <- factor(dotplot_data$Gene, levels = genes_sorted)
avg_expr_long$Gene <- factor(avg_expr_long$Gene, levels = genes_sorted)
log2fc_df$Gene <- factor(log2fc_df$Gene, levels = genes_sorted)
hr_long$Gene <- factor(hr_long$Gene, levels = genes_sorted)

# ====================== Combine and save ======================
combined_plot <- (p_dot / p_expr_heat / p_log2fc / p_heat) +
  plot_layout(heights = c(4, 0.5, 0.2, 0.2), guides = "collect") &
  theme(legend.justification = "bottom")

ggsave(filename = output_pdf,
       plot = combined_plot,
       width = 35,
       height = 15,
       units = "in",
       dpi = 300,
       device = "pdf")