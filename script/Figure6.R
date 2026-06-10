# =========================================================================== #
# Figure A: Boxplots of functional module (CM) activity comparing Pre vs Post
# treatment across NMF subtypes. Four panels show Wilcoxon or paired t-test
# with/without extreme value trimming. Paired tests use patient name.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(CoVarNet)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ------------------------------ File paths ----------------------------------
clinical_path       <- "path/to/Clinical_information.rds"
prop_csv_path       <- "path/to/prop_sample1_clusters.csv"
output_dir          <- "path/to/F6"
wilcox_orig_pdf     <- file.path(output_dir, "A_CM_wilcox_original.pdf")
wilcox_trim_pdf     <- file.path(output_dir, "A_CM_wilcox_trimmed.pdf")
pairedt_orig_pdf    <- file.path(output_dir, "A_CM_pairedt_original.pdf")
pairedt_trim_pdf    <- file.path(output_dir, "A_CM_pairedt_trimmed.pdf")

# --------------------- Load clinical data -----------------------------------
group <- readRDS(clinical_path) %>%
  filter(name != "Unknow")

# ------------------ Load and prepare cluster proportion matrix --------------
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

# Normalize (min-max) and run nsNMF (K = 4)
mat_fq_norm <- freq_normalize(mat_fq_raw, normalize = "minmax")
K <- 4
NMF_K <- nmf(mat_fq_norm, K, method = "nsNMF", seed = rep(77, 6), nrun = 30)
colnames(basis(NMF_K)) <- c("CM01","CM04","CM03","CM02")
rownames(coef(NMF_K)) <- c("CM01","CM04","CM03","CM02")

# Compute pairwise correlation and network (needed for scoef extraction)
cor_pair <- pair_correlation(mat_fq_raw, method = "pearson")
topnum <- floor(nrow(basis(NMF_K)) / ncol(basis(NMF_K)))

# Network construction helper function (silent, simplified)
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

  listA <- Corres[, c("subCluster1","majorCluster1")]; colnames(listA) <- c("subCluster","majorCluster")
  listB <- Corres[, c("subCluster2","majorCluster2")]; colnames(listB) <- c("subCluster","majorCluster")
  ann <- rbind(listA, listB); ann <- ann[!duplicated(ann$subCluster), ]; rownames(ann) <- ann$subCluster

  w <- basis(NMFres)
  if (is.null(colnames(w))) {
    colnames(w) <- sprintf("CM%02d", 1:ncol(w))
    rownames(coef(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
    colnames(basis(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
  }
  w <- w[, order(colnames(w))]

  weight <- reshape2::melt(w)
  colnames(weight) <- c("subCluster","cm","weight")
  weight <- weight %>% group_by(cm) %>% arrange(cm, desc(weight)) %>% ungroup()
  if (apply_top_n) {
    weight_filtered <- weight %>% group_by(cm) %>% arrange(desc(weight)) %>% slice(1:top_n) %>% ungroup()
  } else {
    weight_filtered <- weight
  }

  meta_sc <- data.frame(subCluster = c(Corres$subCluster1, Corres$subCluster2),
                        majorCluster = c(Corres$majorCluster1, Corres$majorCluster2))
  meta_sc <- meta_sc[!duplicated(meta_sc$subCluster), ]; rownames(meta_sc) <- meta_sc$subCluster

  pl_df <- Corres
  if (apply_corr) pl_df <- pl_df[pl_df$correlation > corr, ]
  if (apply_fdr)  pl_df <- pl_df[pl_df$pval_fdr <= fdr, ]
  if (apply_spe) {
    n_all <- nrow(meta_sc)
    spe_cutoff <- 1 - ((top_n - 1) * 2 - 1)/((n_all - 1) * 2 - 1)
    pl_df <- pl_df[pl_df$spe >= spe_cutoff, ]
  }

  graph <- igraph::graph_from_data_frame(pl_df, directed = FALSE)
  node_global <- data.frame(subCluster = igraph::V(graph)$name,
                            majorCluster = meta_sc$majorCluster[match(igraph::V(graph)$name, meta_sc$subCluster)])
  node_each <- data.frame(); edge_each <- data.frame()

  for (cm in unique(weight_filtered$cm)) {
    cm_nodes <- weight_filtered$subCluster[weight_filtered$cm == cm]
    sub_graph <- igraph::subgraph(graph, vids = intersect(igraph::V(graph)$name, cm_nodes))
    if (remove_isolated && length(igraph::V(sub_graph)) > 0)
      sub_graph <- igraph::delete.vertices(sub_graph, igraph::V(sub_graph)[igraph::degree(sub_graph) == 0])
    if (length(igraph::V(sub_graph)) != 0) {
      tmp1 <- data.frame(cm = cm, subCluster = igraph::V(sub_graph)$name)
      node_each <- rbind(node_each, tmp1)
      tmp2 <- pl_df[(pl_df$subCluster1 %in% tmp1$subCluster) &
                    (pl_df$subCluster2 %in% tmp1$subCluster), ]
      edge_each <- rbind(edge_each, data.frame(cm = cm, tmp2))
    }
  }
  if (nrow(node_each) > 0) {
    node_each$majorCluster <- meta_sc$majorCluster[match(node_each$subCluster, meta_sc$subCluster)]
    weight_final <- merge(node_each, weight_filtered)
    weight_final <- weight_final[order(weight_final$cm, -weight_final$weight), ]
    rownames(weight_final) <- 1:nrow(weight_final)
  } else {
    weight_final <- data.frame()
  }
  list(global = list(node = node_global, edge = pl_df),
       each = list(node = node_each, edge = edge_each),
       raw = NMFres, filter = weight_final, ann = ann)
}

network <- sm_cm_network(
  NMF_K, cor_pair,
  corr = 0.2, top_n = topnum,
  apply_fdr = FALSE, apply_corr = TRUE,
  apply_top_n = TRUE, apply_spe = TRUE, remove_isolated = TRUE
)

# Extract CM activity scores (samples × CMs)
cm_activity <- t(scoef(network$raw))

# --------------------------- Plotting parameters ----------------------------
cm_list <- c("CM01","CM02","CM03","CM04")
subtype_list <- c(1,2,3)
treatment_colors <- c("Pre" = "#83a3d5", "Post" = "#dc7081")

# Helper: trim top and bottom values within each treatment group
trim_extremes <- function(df) {
  df %>%
    group_by(treatment_condition) %>%
    arrange(activity, .by_group = TRUE) %>%
    slice(if (n() >= 3) 2:(n() - 1) else 1:n()) %>%
    ungroup()
}

# Helper: compute y-axis range per CM
compute_ylim <- function(cm_activity, group, cm_list, subtype_list, use_trimmed = FALSE) {
  cm_ylim <- list()
  for(cm in cm_list) {
    all_values <- numeric(0)
    for(st in subtype_list) {
      selected_samples <- group %>% filter(subtype == st) %>% pull(sample)
      current_cm <- cm_activity[rownames(cm_activity) %in% selected_samples, cm, drop = FALSE]
      temp_data <- data.frame(sample = rownames(current_cm),
                              activity = current_cm[,1],
                              stringsAsFactors = FALSE) %>%
        left_join(group %>% filter(subtype == st) %>% select(sample, treatment_condition), by = "sample") %>%
        filter(treatment_condition %in% c("Pre","Post")) %>%
        filter(!is.na(activity))
      if(nrow(temp_data) > 0) {
        if(use_trimmed) temp_data <- trim_extremes(temp_data)
        if(nrow(temp_data) > 0) all_values <- c(all_values, temp_data$activity)
      }
    }
    if(length(all_values) > 0) {
      rng <- range(all_values, na.rm = TRUE)
      diff_range <- diff(rng)
      cm_ylim[[cm]] <- c(rng[1] - 0.05 * diff_range, rng[2] + 0.15 * diff_range)
    } else {
      cm_ylim[[cm]] <- c(0,1)
    }
  }
  cm_ylim
}

# Main panel drawing function
draw_CM_panel <- function(cm_activity, group, cm_list, subtype_list,
                          test_method = c("wilcox","paired_t"),
                          trim_extremes_flag = FALSE) {
  test_method <- match.arg(test_method)
  cm_list <- intersect(cm_list, colnames(cm_activity))
  subtype_list <- intersect(subtype_list, unique(group$subtype))
  cm_ylim <- compute_ylim(cm_activity, group, cm_list, subtype_list, use_trimmed = trim_extremes_flag)
  plot_list <- list(); plot_index <- 1

  for(cm in cm_list) {
    for(st in subtype_list) {
      selected_samples <- group %>% filter(subtype == st) %>% pull(sample)
      current_cm <- cm_activity[rownames(cm_activity) %in% selected_samples, cm, drop = FALSE]

      plot_data <- data.frame(sample = rownames(current_cm),
                              activity = current_cm[,1],
                              stringsAsFactors = FALSE) %>%
        left_join(group %>% filter(subtype == st) %>% select(sample, treatment_condition, name), by = "sample") %>%
        filter(treatment_condition %in% c("Pre","Post")) %>%
        filter(!is.na(activity))

      if(trim_extremes_flag) plot_data_trimmed <- trim_extremes(plot_data) else plot_data_trimmed <- plot_data

      if(nrow(plot_data_trimmed) == 0) {
        p <- ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=4) + theme_void() +
          labs(title = ifelse(cm == cm_list[1], paste("Subtype", st), ""),
               y = ifelse(st == 1, paste(cm, "activity"), "")) +
          scale_y_continuous(breaks = c(0,0.5,1), limits = cm_ylim[[cm]])
        if(cm != cm_list[1]) p <- p + theme(plot.title = element_blank())
        if(cm != "CM04") p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        if(st != 1) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
        plot_list[[plot_index]] <- p; plot_index <- plot_index + 1; next
      }

      plot_data_trimmed$treatment_condition <- factor(plot_data_trimmed$treatment_condition, levels = c("Pre","Post"))
      p <- ggplot(plot_data_trimmed, aes(x = treatment_condition, y = activity, fill = treatment_condition)) +
        geom_boxplot(color = "black", outlier.shape = NA, width = 0.4) +
        geom_jitter(width = 0.1, alpha = 0.6, size = 1.2, color = "black") +
        scale_fill_manual(values = treatment_colors) +
        scale_y_continuous(breaks = c(0,0.5,1), limits = cm_ylim[[cm]]) +
        labs(title = ifelse(cm == cm_list[1], paste("Subtype", st), ""),
             x = "", y = ifelse(st == 1, paste(cm, "activity"), "")) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
              axis.title.y = element_text(size = 9, face = "bold"),
              axis.text = element_text(size = 8),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              axis.line = element_line(color = "black", linewidth = 0.5),
              legend.position = "none")
      if(cm != cm_list[1]) p <- p + theme(plot.title = element_blank())
      if(cm != "CM04") p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      if(st != 1) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())

      # Statistical annotation
      pval <- NA
      pre_vals <- plot_data_trimmed$activity[plot_data_trimmed$treatment_condition == "Pre"]
      post_vals <- plot_data_trimmed$activity[plot_data_trimmed$treatment_condition == "Post"]
      if(test_method == "wilcox") {
        if(length(pre_vals) > 0 && length(post_vals) > 0) {
          test_res <- wilcox.test(pre_vals, post_vals)
          pval <- test_res$p.value
        }
        if(!is.na(pval)) {
          label_text <- if(pval < 0.001) "italic(P) < 0.001" else paste0("italic(P) == ", format(pval, digits=3))
          y_pos <- max(plot_data_trimmed$activity, na.rm=TRUE) + 0.05 * diff(range(plot_data_trimmed$activity, na.rm=TRUE))
          p <- p + annotate("text", x = 1.5, y = y_pos, label = label_text, size = 3, hjust = 0.5, parse = TRUE)
        }
      } else {
        paired_wide <- plot_data_trimmed %>%
          filter(!is.na(name)) %>%
          select(name, treatment_condition, activity) %>%
          pivot_wider(names_from = treatment_condition, values_from = activity) %>%
          filter(!is.na(Pre) & !is.na(Post))
        if(nrow(paired_wide) >= 2) {
          test_res <- t.test(paired_wide$Pre, paired_wide$Post, paired = TRUE)
          pval <- test_res$p.value
          label_text <- if(pval < 0.001) "italic(P) < 0.001" else paste0("italic(P) == ", format(pval, digits=3))
        } else {
          label_text <- "N<2 pairs"
        }
        y_pos <- max(plot_data_trimmed$activity, na.rm=TRUE) + 0.05 * diff(range(plot_data_trimmed$activity, na.rm=TRUE))
        p <- p + annotate("text", x = 1.5, y = y_pos, label = label_text, size = 3, hjust = 0.5,
                          parse = (label_text != "N<2 pairs"))
      }
      plot_list[[plot_index]] <- p; plot_index <- plot_index + 1
    }
  }
  wrap_plots(plot_list, ncol = 3)
}

# ---------------------- Generate four panels ---------------------------------
p1 <- draw_CM_panel(cm_activity, group, cm_list, subtype_list, test_method = "wilcox", trim_extremes_flag = FALSE)
p1 <- p1 + plot_annotation(title = "Wilcoxon, original data")
p2 <- draw_CM_panel(cm_activity, group, cm_list, subtype_list, test_method = "wilcox", trim_extremes_flag = TRUE)
p2 <- p2 + plot_annotation(title = "Wilcoxon, trimmed extremes")
p3 <- draw_CM_panel(cm_activity, group, cm_list, subtype_list, test_method = "paired_t", trim_extremes_flag = FALSE)
p3 <- p3 + plot_annotation(title = "Paired t-test, original data")
p4 <- draw_CM_panel(cm_activity, group, cm_list, subtype_list, test_method = "paired_t", trim_extremes_flag = TRUE)
p4 <- p4 + plot_annotation(title = "Paired t-test, trimmed extremes")

ggsave(wilcox_orig_pdf,  p1, width = 5.5, height = 8, dpi = 500)
ggsave(wilcox_trim_pdf,  p2, width = 5.5, height = 8, dpi = 500)
ggsave(pairedt_orig_pdf, p3, width = 5.5, height = 8, dpi = 500)
ggsave(pairedt_trim_pdf, p4, width = 5.5, height = 8, dpi = 500)



# =========================================================================== #
# Figure B: Circular heatmap of integrated cluster activity across subtypes
# and treatment conditions (Pre/Post). Outer rings show CM group membership,
# treatment condition, and per-subtype mean activity for each cluster.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(CoVarNet)
library(dplyr)
library(stringr)
library(tidyr)
library(circlize)
library(plotrix)

# ------------------------------ File paths ----------------------------------
clinical_path <- "path/to/Clinical_information.rds"
prop_csv_path <- "path/to/prop_sample1_clusters.csv"
output_pdf    <- "path/to/F6/B_circlize_heatmap.pdf"

# --------------------- Load and prepare clinical data -----------------------
group <- readRDS(clinical_path) %>% filter(name != "Unknow")

# ----------------- Load cluster proportion matrix ---------------------------
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

# Normalize and run nsNMF (K = 4)
mat_fq_norm <- freq_normalize(mat_fq_raw, normalize = "minmax")
K <- 4
NMF_K <- nmf(mat_fq_norm, K, method = "nsNMF", seed = rep(77, 6), nrun = 30)
colnames(basis(NMF_K)) <- c("CM01","CM04","CM03","CM02")
rownames(coef(NMF_K)) <- c("CM01","CM04","CM03","CM02")

# Extract CM activity and integrated activity matrix
cm_activity <- t(scoef(NMF_K))
all_cluster_weights <- basis(NMF_K)
integrated_activity_matrix <- cm_activity %*% t(all_cluster_weights)

# ------------------ Prepare plotting data -----------------------------------
activity_df <- as.data.frame(integrated_activity_matrix)
activity_df$sample <- rownames(activity_df)

activity_long <- activity_df %>%
  pivot_longer(cols = starts_with("cluster"),
               names_to = "cluster",
               values_to = "activity")

group <- group %>% filter(treatment_condition %in% c("Pre", "Post"))
activity_merged <- merge(activity_long, group, by = "sample")

# Compute mean activity per cluster, treatment and subtype
pre_activity <- activity_merged %>%
  filter(treatment_condition == "Pre") %>%
  group_by(cluster, subtype) %>%
  summarise(mean_activity = mean(activity, na.rm = TRUE), .groups = "drop")

post_activity <- activity_merged %>%
  filter(treatment_condition == "Post") %>%
  group_by(cluster, subtype) %>%
  summarise(mean_activity = mean(activity, na.rm = TRUE), .groups = "drop")

pre_wide <- pre_activity %>%
  pivot_wider(names_from = subtype, values_from = mean_activity,
              names_prefix = "subtype_")
post_wide <- post_activity %>%
  pivot_wider(names_from = subtype, values_from = mean_activity,
              names_prefix = "subtype_")

all_clusters <- paste0("cluster", 0:49)

pre_complete <- data.frame(cluster = all_clusters) %>%
  left_join(pre_wide, by = "cluster")
post_complete <- data.frame(cluster = all_clusters) %>%
  left_join(post_wide, by = "cluster")
pre_complete[is.na(pre_complete)] <- 0
post_complete[is.na(post_complete)] <- 0

# Build per-subtype matrices (50 clusters x 2 treatments)
subtype_matrices <- list()
for (i in 1:3) {
  subtype_name <- paste0("subtype_", i)
  mat <- matrix(0, nrow = 50, ncol = 2)
  rownames(mat) <- all_clusters
  colnames(mat) <- c("Pre", "Post")
  mat[, 1] <- as.numeric(pre_complete[[subtype_name]])
  mat[, 2] <- as.numeric(post_complete[[subtype_name]])
  subtype_matrices[[i]] <- mat
}

# CM colours
cm_colors <- c("#9DD0C7", "#9180AC", "#D9BDD8", "#E58579")  # CM01..CM04

# Define cluster groups (with duplicates allowed for shared clusters)
CM01 <- c("cluster1","cluster15","cluster18","cluster20","cluster25",
          "cluster29","cluster30","cluster41","cluster44")
CM02 <- c("cluster2","cluster3","cluster6","cluster7","cluster8",
          "cluster10","cluster39","cluster46","cluster49")
CM03 <- c("cluster8","cluster9","cluster11","cluster13","cluster30",
          "cluster32","cluster40","cluster47")
CM04 <- c("cluster0","cluster12","cluster22","cluster30","cluster31",
          "cluster36","cluster38","cluster42")

selected_clusters <- c(CM01, CM02, CM03, CM04)

# Create unique sector names (append suffix for duplicates)
selected_clusters_unique <- character(length(selected_clusters))
cluster_count <- list()
for (i in seq_along(selected_clusters)) {
  cl <- selected_clusters[i]
  if (is.null(cluster_count[[cl]])) cluster_count[[cl]] <- 1
  else cluster_count[[cl]] <- cluster_count[[cl]] + 1
  selected_clusters_unique[i] <- paste0(cl, "_", cluster_count[[cl]])
}

# Group boundaries for gaps
group_sizes <- c(length(CM01), length(CM02), length(CM03), length(CM04))
group_ends <- cumsum(group_sizes)
n_clusters <- length(selected_clusters_unique)
gap_degrees <- rep(1, n_clusters)
gap_degrees[n_clusters] <- 15
cm_group_vector <- rep(c("CM01","CM02","CM03","CM04"), times = group_sizes)

# Activity colour map
col_fun <- colorRamp2(seq(0, 0.05, length.out = 100),
                      colorRampPalette(c("#3583B4","white","#E84D94"))(100))

# ------------------------ Draw circular heatmap -----------------------------
pdf(output_pdf, width = 5000/300, height = 4000/300)
par(mar = c(1, 1, 3, 1))

circos.clear()
circos.par(gap.degree = gap_degrees, cell.padding = c(0, 0, 0, 0),
           start.degree = 90)
circos.initialize(factors = factor(selected_clusters_unique,
                                   levels = selected_clusters_unique),
                  xlim = c(0, 2))

# Track 1: CM group background
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {},
             track.height = 0.05, bg.border = NA)

for (cm_grp in c("CM01","CM02","CM03","CM04")) {
  idx <- which(cm_group_vector == cm_grp)
  sectors <- selected_clusters_unique[idx]
  col <- cm_colors[which(c("CM01","CM02","CM03","CM04") == cm_grp)]
  highlight.sector(sectors, track.index = 1, col = col,
                   border = col, lwd = 0.1)
}

# Track 2: Treatment condition (Pre/Post)
treatment_colors <- c("#FFA500", "#8B4513")
circos.track(ylim = c(0, 1),
  panel.fun = function(x, y) {
    circos.rect(0, 0, 1, 1, col = treatment_colors[1], border = "white", lwd = 0.5)
    circos.rect(1, 0, 2, 1, col = treatment_colors[2], border = "white", lwd = 0.5)
  }, track.height = 0.05, bg.border = NA)

# Track 3: Subtype activity (three rows combined)
circos.track(ylim = c(0, 3),
  panel.fun = function(x, y) {
    sector.index <- CELL_META$sector.index
    original_cluster <- sub("_[0-9]+$", "", sector.index)
    cluster_index <- which(all_clusters == original_cluster)
    for (i in 1:3) {
      y_bottom <- 3 - i
      y_top <- y_bottom + 1
      mat <- subtype_matrices[[i]]
      if (is.matrix(mat) && nrow(mat) >= cluster_index) {
        pre_val  <- mat[cluster_index, 1]
        post_val <- mat[cluster_index, 2]
        circos.rect(0, y_bottom, 1, y_top, col = col_fun(pre_val),
                    border = "white", lwd = 0.5)
        circos.rect(1, y_bottom, 2, y_top, col = col_fun(post_val),
                    border = "white", lwd = 0.5)
      }
    }
  }, track.height = 0.35, bg.border = NA)

# Track 4: Cluster labels
circos.track(ylim = c(0, 1),
  panel.fun = function(x, y) {
    sector.index <- CELL_META$sector.index
    original_cluster <- sub("_[0-9]+$", "", sector.index)
    cl_num <- gsub("cluster", "", original_cluster)
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] - mm_y(20),
                paste0("Cluster ", cl_num), facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), cex = 1.5)
  }, track.height = 0.05, bg.border = NA)

# External legends
legend_x <- 1.05
legend_y <- 0.9
spacing  <- 0.25

legend(legend_x, legend_y,
       legend = c("Pre", "Post"), fill = treatment_colors,
       title = "Treatment", cex = 1.5, xpd = TRUE, bty = "n")
legend(legend_x, legend_y - spacing,
       legend = c("CM01","CM02","CM03","CM04"), fill = cm_colors,
       title = "CM Groups", cex = 1.5, xpd = TRUE, bty = "n")
legend(legend_x, legend_y - 2*spacing,
       legend = c("Subtype 1","Subtype 2","Subtype 3"),
       title = "", cex = 1.5, xpd = TRUE, bty = "n")

color.legend(legend_x, legend_y - 3*spacing - 0.03,
             legend_x + 0.15, legend_y - 3*spacing,
             legend = round(c(0, 0.05), 2),
             rect.col = colorRampPalette(c("#3583B4","white","#E84D94"))(100),
             gradient = "x", align = "rb")
text(legend_x + 0.075, legend_y - 3*spacing + 0.01,
     "Activity", cex = 1.5, xpd = TRUE, adj = c(0.5, 0))

circos.clear()
dev.off()



# =========================================================================== #
# Figure F: Heatmap of CM04-associated cluster frequencies along pseudotime.
# Figure G: PHATE plots coloured by treatment, pseudotime, and CM04 activity.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(CoVarNet)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(reticulate)

# Python environment (adjust if needed)
use_python("~/anaconda3/envs/screnv/bin/python", required = TRUE)
sc <- reticulate::import("scanpy")
pa <- reticulate::import("palantir", convert = FALSE)

# ------------------------------ File paths ----------------------------------
clinical_path  <- "path/to/Clinical_information.rds"
prop_csv_path  <- "path/to/prop_sample1_clusters.csv"
output_dir     <- "path/to/F6"
heatmap_pdf    <- file.path(output_dir, "F_heatmap_CM04.pdf")
phate_pdf      <- file.path(output_dir, "G_combined_Clusters_treatment_subtype_Pseudotime.pdf")

# --------------------- Load and prepare clinical data -----------------------
group <- readRDS(clinical_path) %>% filter(name != "Unknow")

# ----------------- Load cluster proportion matrix ---------------------------
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

# ---------------- NMF decomposition and network construction ----------------
mat_fq_norm <- freq_normalize(mat_fq_raw, normalize = "minmax")
NMF_K <- nmf(mat_fq_norm, 4, method = "nsNMF", seed = rep(77, 6), nrun = 30)
colnames(basis(NMF_K)) <- c("CM01","CM04","CM03","CM02")
rownames(coef(NMF_K)) <- c("CM01","CM04","CM03","CM02")

cor_pair <- pair_correlation(mat_fq_raw, method = "pearson")
topnum <- floor(nrow(basis(NMF_K)) / ncol(basis(NMF_K)))

# Network construction function (silent, used to extract CM activity)
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

  listA <- Corres[, c("subCluster1","majorCluster1")]; colnames(listA) <- c("subCluster","majorCluster")
  listB <- Corres[, c("subCluster2","majorCluster2")]; colnames(listB) <- c("subCluster","majorCluster")
  ann <- rbind(listA, listB); ann <- ann[!duplicated(ann$subCluster), ]; rownames(ann) <- ann$subCluster

  w <- basis(NMFres)
  if (is.null(colnames(w))) {
    colnames(w) <- sprintf("CM%02d", 1:ncol(w))
    rownames(coef(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
    colnames(basis(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
  }
  w <- w[, order(colnames(w))]

  weight <- reshape2::melt(w)
  colnames(weight) <- c("subCluster","cm","weight")
  weight <- weight %>% group_by(cm) %>% arrange(cm, desc(weight)) %>% ungroup()
  weight_filtered <- if (apply_top_n) {
    weight %>% group_by(cm) %>% slice(1:top_n) %>% ungroup()
  } else weight

  meta_sc <- data.frame(subCluster = c(Corres$subCluster1, Corres$subCluster2),
                        majorCluster = c(Corres$majorCluster1, Corres$majorCluster2))
  meta_sc <- meta_sc[!duplicated(meta_sc$subCluster), ]; rownames(meta_sc) <- meta_sc$subCluster

  pl_df <- Corres
  if (apply_corr) pl_df <- pl_df[pl_df$correlation > corr, ]
  if (apply_fdr)  pl_df <- pl_df[pl_df$pval_fdr <= fdr, ]
  if (apply_spe) {
    n_all <- nrow(meta_sc)
    spe_cutoff <- 1 - ((top_n - 1) * 2 - 1) / ((n_all - 1) * 2 - 1)
    pl_df <- pl_df[pl_df$spe >= spe_cutoff, ]
  }

  graph <- igraph::graph_from_data_frame(pl_df, directed = FALSE)
  node_global <- data.frame(subCluster = igraph::V(graph)$name,
                            majorCluster = meta_sc$majorCluster[match(igraph::V(graph)$name, meta_sc$subCluster)])
  node_each <- data.frame(); edge_each <- data.frame()

  for (cm in unique(weight_filtered$cm)) {
    cm_nodes <- weight_filtered$subCluster[weight_filtered$cm == cm]
    sub_graph <- igraph::subgraph(graph, vids = intersect(igraph::V(graph)$name, cm_nodes))
    if (remove_isolated && length(igraph::V(sub_graph)) > 0)
      sub_graph <- igraph::delete.vertices(sub_graph, igraph::V(sub_graph)[igraph::degree(sub_graph) == 0])
    if (length(igraph::V(sub_graph)) != 0) {
      tmp1 <- data.frame(cm = cm, subCluster = igraph::V(sub_graph)$name)
      node_each <- rbind(node_each, tmp1)
      tmp2 <- pl_df[(pl_df$subCluster1 %in% tmp1$subCluster) &
                    (pl_df$subCluster2 %in% tmp1$subCluster), ]
      edge_each <- rbind(edge_each, data.frame(cm = cm, tmp2))
    }
  }
  if (nrow(node_each) > 0) {
    node_each$majorCluster <- meta_sc$majorCluster[match(node_each$subCluster, meta_sc$subCluster)]
    weight_final <- merge(node_each, weight_filtered)
    weight_final <- weight_final[order(weight_final$cm, -weight_final$weight), ]
    rownames(weight_final) <- 1:nrow(weight_final)
  } else weight_final <- data.frame()

  list(global = list(node = node_global, edge = pl_df),
       each = list(node = node_each, edge = edge_each),
       raw = NMFres, filter = weight_final, ann = ann)
}

network <- sm_cm_network(
  NMF_K, cor_pair, corr = 0.2, top_n = topnum,
  apply_fdr = FALSE, apply_corr = TRUE,
  apply_top_n = TRUE, apply_spe = TRUE, remove_isolated = TRUE
)

# Extract CM activity and CM04 subsets
cm_activity <- t(scoef(network$raw))
subsets <- network$filter[network$filter$cm == "CM04", "subCluster"]

# ----------------------- Prepare metadata for PHATE -------------------------
ann <- group
rownames(ann) <- ann$sample
sample <- unique(ann[ann$treatment_condition != "Unknow", "sample"])
sample <- sort(sample)

meta <- cbind(t(mat_fq_raw[, sample]), ann[sample, ], cm_activity[sample, ])

# Reproducible PHATE reduction
freq_reduction_reproducible <- function(meta, subsets, seed = 123) {
  set.seed(seed)
  sc$settings$seed <- as.integer(seed)
  mat_gb_norm <- meta[, subsets]
  mat_gb_norm <- t(mat_gb_norm)
  mat_gb_norm <- mat_gb_norm - apply(mat_gb_norm, 1, mean)
  mat_ggb_norm <- mat_gb_norm / apply(mat_gb_norm, 1, sd)
  mat_gb_norm <- t(mat_gb_norm)
  data <- anndata::AnnData(X = mat_gb_norm, obsm = list(meta.data = meta))
  sc$pp$neighbors(data, random_state = as.integer(seed))
  sc$tl$leiden(data, random_state = as.integer(seed))
  data$obs["Clusters"] <- plyr::mapvalues(
    data$obs$leiden, levels(data$obs$leiden), seq_along(levels(data$obs$leiden))
  )
  sc$external$tl$phate(data, n_components = 2L, a = 40L, random_state = as.integer(seed))
  data$obsm$meta.data["phate1"] <- data$obsm$X_phate[, 1]
  data$obsm$meta.data["phate2"] <- data$obsm$X_phate[, 2]
  data$obsm$meta.data["Clusters"] <- data$obs$Clusters
  data
}

set.seed(1234)
data <- freq_reduction_reproducible(meta, subsets)

# --------------- Figure G: PHATE plots --------------------------------------
p1 <- gr.phate(data, "treatment_condition", color = c("#EA738D", "#89ABE3", "#FFC000"))

root <- "1"
data <- pseudotime(data, root)
p2 <- gr.phate(data, "Pseudotime", color = RColorBrewer::brewer.pal(9, "Purples")[3:7])

CMnum <- "CM04"
spectral_tail <- tail(RColorBrewer::brewer.pal(11, "Spectral"), 6)
p3_colors <- c(spectral_tail, "#3F007D")
p3 <- gr.phate(data, CMnum, color = rev(p3_colors))

combined_g <- p1 + p2 + p3 + plot_layout(ncol = 3)
ggsave(phate_pdf, combined_g, width = 6, height = 1.5, units = "in", limitsize = FALSE)

# --------------- Figure F: Trajectory heatmap -------------------------------
gr.trajectorynew <- function(data, subset, rowann, ann_colors = list()) {
  res <- data$obsm$meta.data
  res <- res[order(res$Pseudotime), ]
  ann_col <- ComplexHeatmap::HeatmapAnnotation(
    df = res[, rowann],
    col = ann_colors,
    border = FALSE,
    annotation_name_gp = grid::gpar(fontsize = 6),
    annotation_name_side = "left",
    simple_anno_size = unit(2, "mm"),
    annotation_legend_param = list(
      labels_gp = grid::gpar(fontsize = 6),
      title_gp = grid::gpar(fontsize = 6),
      grid_width = unit(3, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  X <- as.data.frame(scale(res[, subset]))
  ComplexHeatmap::Heatmap(
    matrix = t(X),
    col = circlize::colorRamp2(seq(-1, 2, length.out = 4),
                               RColorBrewer::brewer.pal(7, "RdBu")[4:1]),
    border = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    width = unit(6, "cm"),
    height = unit(12, "cm"),
    column_title_gp = grid::gpar(fontsize = 7),
    top_annotation = ann_col,
    show_column_names = FALSE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 6),
    column_dend_height = unit(5, "mm"),
    row_dend_width = unit(5, "mm"),
    heatmap_legend_param = list(
      title = "Subset freq.",
      at = c(-1, 0, 1, 2),
      legend_direction = "horizontal",
      title_position = "leftcenter",
      grid_height = unit(3, "mm"),
      legend_width = unit(1, "cm"),
      labels_gp = grid::gpar(fontsize = 6),
      title_gp = grid::gpar(fontsize = 6)
    )
  )
}

# Colour mapping for CM04
cm_vals <- data$obsm$meta.data[[CMnum]]
breaks <- seq(min(cm_vals, na.rm = TRUE), max(cm_vals, na.rm = TRUE),
              length.out = length(p3_colors))
ann_colors <- setNames(list(circlize::colorRamp2(breaks, rev(p3_colors))), CMnum)

sorted_clusters <- subsets[order(as.numeric(gsub("cluster", "", subsets)))]

pdf(heatmap_pdf, width = 5, height = 6)
set.seed(1234)
gr.trajectorynew(
  data       = data,
  subset     = sorted_clusters,
  rowann     = c("Pseudotime", CMnum, "treatment_condition"),
  ann_colors = ann_colors
)
dev.off()



# =========================================================================== #
# Figure H: Stacked bar charts showing the number of Pre/Post differentially
# expressed genes (DEGs) per cluster within each NMF subtype, and a heatmap
# of overlap ratios with cluster marker genes. Analyses run per subtype.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(RColorBrewer)

# ------------------------------ File paths ----------------------------------
clinical_path      <- "path/to/Clinical_information.rds"
scRNA_path         <- "path/to/01merge.rds"
marker_all_path    <- "path/to/04marker_cluster_m01_l00_r05.txt"
output_pdf         <- "path/to/F6/H_marker.pdf"

# --------------------------- Load data --------------------------------------
group <- readRDS(clinical_path) %>%
  filter(name != "Unknow",
         subtype != "Unknow",
         treatment_condition %in% c("Pre", "Post"))

scRNA <- readRDS(scRNA_path)
scRNA <- subset(scRNA, subset = orig.ident %in% group$sample)

marker_all <- read.table(marker_all_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

# ====================== Function: DEGs per cluster ==========================
run_subtype_degs <- function(scRNA, group, subtype_val) {
  # Subset cells of the given subtype
  cells <- subset(scRNA,
                  subset = orig.ident %in% group$sample[group$subtype == subtype_val])
  cells$treatment_condition <- group$treatment_condition[
    match(cells$orig.ident, group$sample)
  ]
  Idents(cells) <- "seurat_clusters"
  clusters <- as.character(unique(cells$seurat_clusters))
  
  marker_results <- list()
  for (cl in clusters) {
    cluster_cells <- subset(cells, idents = cl)
    Idents(cluster_cells) <- "treatment_condition"
    cell_counts <- table(cluster_cells$treatment_condition)
    
    if (all(c("Pre", "Post") %in% names(cell_counts)) &&
        all(cell_counts[c("Pre", "Post")] >= 3)) {
      markers <- tryCatch(
        FindMarkers(cluster_cells, ident.1 = "Post", ident.2 = "Pre",
                    min.pct = 0.1, logfc.threshold = 0, return.thresh = 0.05),
        error = function(e) NULL
      )
      if (!is.null(markers) && nrow(markers) > 0) {
        markers$cluster <- cl
        markers$gene <- rownames(markers)
        rownames(markers) <- NULL
        marker_results[[cl]] <- markers
      }
    }
  }
  
  res <- if (length(marker_results) > 0) {
    do.call(rbind, marker_results) %>%
      arrange(as.numeric(cluster)) %>%
      filter(p_val_adj < 0.05)
  } else {
    data.frame()
  }
  res
}

# ====================== Function: overlap ratio =============================
calc_overlap_ratio <- function(subtype_markers, marker_all) {
  CM01 <- c("cluster1","cluster15","cluster18","cluster20","cluster25",
            "cluster29","cluster30","cluster41","cluster44")
  CM02 <- c("cluster2","cluster3","cluster6","cluster7","cluster8",
            "cluster10","cluster39","cluster46","cluster49")
  CM03 <- c("cluster8","cluster9","cluster11","cluster13","cluster30",
            "cluster32","cluster40","cluster47")
  CM04 <- c("cluster0","cluster12","cluster22","cluster30","cluster31",
            "cluster36","cluster38","cluster42")
  
  extract_number <- function(x) as.numeric(gsub("cluster", "", x))
  full_order <- c(extract_number(CM01), extract_number(CM02),
                  extract_number(CM03), extract_number(CM04))
  
  sapply(full_order, function(cl) {
    cl_chr <- as.character(cl)
    deg_genes <- subtype_markers %>%
      filter(cluster == cl_chr) %>%
      pull(gene)
    marker_genes <- marker_all %>%
      filter(cluster == cl_chr) %>%
      pull(gene)
    if (length(marker_genes) > 0) {
      length(intersect(deg_genes, marker_genes)) / length(marker_genes)
    } else 0
  })
}

# ====================== Function: combined plot =============================
plot_subtype_markers <- function(markers_data, plot_title, overlap_vec,
                                 show_x_axis = FALSE, show_y_axis = FALSE) {
  if (nrow(markers_data) == 0) return(NULL)
  
  CM01 <- c("cluster1","cluster15","cluster18","cluster20","cluster25",
            "cluster29","cluster30","cluster41","cluster44")
  CM02 <- c("cluster2","cluster3","cluster6","cluster7","cluster8",
            "cluster10","cluster39","cluster46","cluster49")
  CM03 <- c("cluster8","cluster9","cluster11","cluster13","cluster30",
            "cluster32","cluster40","cluster47")
  CM04 <- c("cluster0","cluster12","cluster22","cluster30","cluster31",
            "cluster36","cluster38","cluster42")
  
  extract_cluster_number <- function(x) as.numeric(gsub("cluster", "", x))
  full_order <- c(extract_cluster_number(CM01), extract_cluster_number(CM02),
                  extract_cluster_number(CM03), extract_cluster_number(CM04))
  
  # Build counts
  cluster_counts <- lapply(seq_along(full_order), function(i) {
    cl <- as.character(full_order[i])
    if (cl %in% markers_data$cluster) {
      sub <- markers_data[markers_data$cluster == cl, ]
      up <- sum(sub$avg_log2FC > 0, na.rm = TRUE)
      down <- sum(sub$avg_log2FC < 0, na.rm = TRUE)
    } else { up <- 0; down <- 0 }
    data.frame(cluster_label = cl, position = i, up = up, down = down)
  }) %>% bind_rows()
  
  cluster_long <- cluster_counts %>%
    pivot_longer(cols = c(up, down), names_to = "direction", values_to = "count") %>%
    mutate(position = factor(position, levels = seq_along(full_order)),
           direction = factor(direction, levels = c("down", "up"),
                              labels = c("Pre", "Post")))
  
  x_labels <- cluster_counts$cluster_label[!duplicated(cluster_counts$position)]
  
  # Bar plot
  p_bar <- ggplot(cluster_long, aes(x = position, y = count, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7, position = "stack") +
    scale_fill_manual(values = c("Pre" = "#377EB8", "Post" = "#E41A1C"),
                      name = "Expression in Post") +
    scale_x_discrete(labels = x_labels, expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    labs(title = plot_title,
         x = if (show_x_axis) "Cluster" else "",
         y = if (show_y_axis) "Number of DEGs" else "") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, face = "bold", color = "black"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      plot.margin = unit(c(5.5, 5.5, 0, 5.5), "points")
    )
  if (!show_y_axis) p_bar <- p_bar + theme(axis.title.y = element_blank())
  
  # Heatmap of overlap
  heat_data <- data.frame(
    position = factor(seq_along(full_order), levels = seq_along(full_order)),
    overlap = overlap_vec
  )
  rdylgn_colors <- rev(brewer.pal(11, "PRGn")[2:8])
  
  p_heat <- ggplot(heat_data, aes(x = position, y = 1, fill = overlap)) +
    geom_tile(color = "white", linewidth = 0.5, width = 0.7, height = 0.8) +
    scale_fill_gradientn(colors = rdylgn_colors, limits = c(0, 0.5),
                         oob = scales::squish, name = "Overlap\nratio") +
    scale_x_discrete(labels = x_labels, expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme_minimal() +
    theme(
      axis.text.x = if (show_x_axis) element_text(angle = 0, hjust = 0.5,
                                                  vjust = 0.5, size = 9,
                                                  color = "black") else element_blank(),
      axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
      axis.title.x = if (show_x_axis) element_text(size = 12, face = "bold") else element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = unit(c(0, 5.5, 5.5, 5.5), "points")
    )
  
  p_bar / p_heat + plot_layout(heights = c(6, 1))
}

# ==================== Run DEG analysis per subtype ==========================
subtype1markers <- run_subtype_degs(scRNA, group, 1)
subtype2markers <- run_subtype_degs(scRNA, group, 2)
subtype3markers <- run_subtype_degs(scRNA, group, 3)

# Compute overlap ratios
ov1 <- calc_overlap_ratio(subtype1markers, marker_all)
ov2 <- calc_overlap_ratio(subtype2markers, marker_all)
ov3 <- calc_overlap_ratio(subtype3markers, marker_all)

# Generate plots
p1 <- plot_subtype_markers(subtype1markers, "Subtype 1", overlap_vec = ov1,
                           show_x_axis = FALSE, show_y_axis = FALSE)
p2 <- plot_subtype_markers(subtype2markers, "Subtype 2", overlap_vec = ov2,
                           show_x_axis = FALSE, show_y_axis = TRUE)
p3 <- plot_subtype_markers(subtype3markers, "Subtype 3", overlap_vec = ov3,
                           show_x_axis = TRUE, show_y_axis = FALSE)

# Combine and save
plots_list <- list(p1, p2, p3)
valid_plots <- !sapply(plots_list, is.null)

if (sum(valid_plots) > 0) {
  combined_plot <- wrap_plots(plots_list[valid_plots], ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
  ggsave(output_pdf, plot = combined_plot,
         width = 12, height = 2 * sum(valid_plots), dpi = 300)
}



# =========================================================================== #
# Figure I: GSVA score trends along pseudotime for CM04-associated clusters
# (clusters 1, 6, 9). Loess smoothing with 95% confidence bands.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(CoVarNet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GSVA)
library(reticulate)
use_python("~/anaconda3/envs/screnv/bin/python", required = TRUE)
sc <- reticulate::import("scanpy")

# ------------------------------ File paths ----------------------------------
clinical_path   <- "path/to/Clinical_information.rds"
prop_csv_path   <- "path/to/prop_sample1_clusters.csv"
gsva_rds_path   <- "path/to/gsva_output.rds"
metadata_path   <- "path/to/metadata.rds"
output_pdf      <- "path/to/F6/I_GSVA_trend_clusters_1_6_9.pdf"

# ------------------ 1. Load and prepare clinical data -----------------------
group <- readRDS(clinical_path) %>% filter(name != "Unknow")

# ------------------ 2. Load cluster proportion matrix -----------------------
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

# ------------------ 3. NMF and network construction -------------------------
mat_fq_norm <- freq_normalize(mat_fq_raw, normalize = "minmax")
NMF_K <- nmf(mat_fq_norm, 4, method = "nsNMF", seed = rep(77, 6), nrun = 30)
colnames(basis(NMF_K)) <- c("CM01","CM04","CM03","CM02")
rownames(coef(NMF_K)) <- c("CM01","CM04","CM03","CM02")

cor_pair <- pair_correlation(mat_fq_raw, method = "pearson")
topnum <- floor(nrow(basis(NMF_K)) / ncol(basis(NMF_K)))

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

  listA <- Corres[, c("subCluster1","majorCluster1")]; colnames(listA) <- c("subCluster","majorCluster")
  listB <- Corres[, c("subCluster2","majorCluster2")]; colnames(listB) <- c("subCluster","majorCluster")
  ann <- rbind(listA, listB); ann <- ann[!duplicated(ann$subCluster), ]; rownames(ann) <- ann$subCluster

  w <- basis(NMFres)
  if (is.null(colnames(w))) {
    colnames(w) <- sprintf("CM%02d", 1:ncol(w))
    rownames(coef(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
    colnames(basis(NMFres)) <- sprintf("CM%02d", 1:ncol(w))
  }
  w <- w[, order(colnames(w))]

  weight <- reshape2::melt(w)
  colnames(weight) <- c("subCluster","cm","weight")
  weight <- weight %>% group_by(cm) %>% arrange(cm, desc(weight)) %>% ungroup()
  if (apply_top_n) weight_filtered <- weight %>% group_by(cm) %>% slice(1:max(1, top_n)) %>% ungroup() else weight_filtered <- weight

  meta_sc <- data.frame(subCluster = c(Corres$subCluster1, Corres$subCluster2),
                        majorCluster = c(Corres$majorCluster1, Corres$majorCluster2))
  meta_sc <- meta_sc[!duplicated(meta_sc$subCluster), ]; rownames(meta_sc) <- meta_sc$subCluster

  pl_df <- Corres
  if (apply_corr) pl_df <- pl_df[pl_df$correlation > corr, ]
  if (apply_fdr)  pl_df <- pl_df[pl_df$pval_fdr <= fdr, ]
  if (apply_spe) {
    n_all <- nrow(meta_sc)
    spe_cutoff <- 1 - ((top_n - 1) * 2 - 1) / ((n_all - 1) * 2 - 1)
    pl_df <- pl_df[pl_df$spe >= spe_cutoff, ]
  }

  graph <- igraph::graph_from_data_frame(pl_df, directed = FALSE)
  node_global <- data.frame(subCluster = igraph::V(graph)$name,
                            majorCluster = meta_sc$majorCluster[match(igraph::V(graph)$name, meta_sc$subCluster)])
  node_each <- data.frame(); edge_each <- data.frame()

  for (cm in unique(weight_filtered$cm)) {
    cm_nodes <- weight_filtered$subCluster[weight_filtered$cm == cm]
    sub_graph <- igraph::subgraph(graph, vids = intersect(igraph::V(graph)$name, cm_nodes))
    if (remove_isolated && length(igraph::V(sub_graph)) > 0)
      sub_graph <- igraph::delete.vertices(sub_graph, igraph::V(sub_graph)[igraph::degree(sub_graph) == 0])
    if (length(igraph::V(sub_graph)) != 0) {
      tmp1 <- data.frame(cm = cm, subCluster = igraph::V(sub_graph)$name)
      node_each <- rbind(node_each, tmp1)
      tmp2 <- pl_df[(pl_df$subCluster1 %in% tmp1$subCluster) &
                    (pl_df$subCluster2 %in% tmp1$subCluster), ]
      edge_each <- rbind(edge_each, data.frame(cm = cm, tmp2))
    }
  }
  if (nrow(node_each) > 0) {
    node_each$majorCluster <- meta_sc$majorCluster[match(node_each$subCluster, meta_sc$subCluster)]
    weight_final <- merge(node_each, weight_filtered)
    weight_final <- weight_final[order(weight_final$cm, -weight_final$weight), ]
    rownames(weight_final) <- 1:nrow(weight_final)
  } else weight_final <- data.frame()

  list(global = list(node = node_global, edge = pl_df),
       each = list(node = node_each, edge = edge_each),
       raw = NMFres, filter = weight_final, ann = ann)
}

network <- sm_cm_network(
  NMF_K, cor_pair,
  corr = 0.2, top_n = topnum,
  apply_fdr = FALSE, apply_corr = TRUE,
  apply_top_n = TRUE, apply_spe = TRUE, remove_isolated = TRUE
)

# CM04 subsets
subsets <- network$filter[network$filter$cm == "CM04", "subCluster"]

# ------------------ 4. Metadata and PHATE/pseudotime ------------------------
ann <- group; rownames(ann) <- ann$sample
sample <- sort(unique(ann[ann$treatment_condition != "Unknow", "sample"]))
cm_activity <- t(scoef(network$raw))
meta <- cbind(t(mat_fq_raw[, sample]), ann[sample, ], cm_activity[sample, ])

freq_reduction_reproducible <- function(meta, subsets, seed = 123) {
  set.seed(seed)
  sc$settings$seed <- as.integer(seed)
  mat_gb_norm <- meta[, subsets]
  mat_gb_norm <- t(mat_gb_norm)
  mat_gb_norm <- mat_gb_norm - apply(mat_gb_norm, 1, mean)
  mat_gb_norm <- mat_gb_norm / apply(mat_gb_norm, 1, sd)
  mat_gb_norm <- t(mat_gb_norm)
  data <- anndata::AnnData(X = mat_gb_norm, obsm = list(meta.data = meta))
  sc$pp$neighbors(data, random_state = as.integer(seed))
  sc$tl$leiden(data, random_state = as.integer(seed))
  data$obs["Clusters"] <- plyr::mapvalues(
    data$obs$leiden, levels(data$obs$leiden), seq_along(levels(data$obs$leiden)))
  sc$external$tl$phate(data, n_components = 2L, a = 40L, random_state = as.integer(seed))
  data$obsm$meta.data["phate1"] <- data$obsm$X_phate[, 1]
  data$obsm$meta.data["phate2"] <- data$obsm$X_phate[, 2]
  data$obsm$meta.data["Clusters"] <- data$obs$Clusters
  data
}

data <- freq_reduction_reproducible(meta, subsets)
root <- "1"
data <- pseudotime(data, root)
sample_Pseudotime <- data.frame(
  SampleID = rownames(data$obsm$meta.data),
  Pseudotime = data$obsm$meta.data$Pseudotime,
  stringsAsFactors = FALSE
)

# ------------------ 5. GSVA trend per cluster -------------------------------
gsva_matrix <- readRDS(gsva_rds_path)
metadata <- readRDS(metadata_path)

plot_cluster_gsva <- function(cluster_id) {
  gsva_long <- as.data.frame(gsva_matrix) %>%
    rownames_to_column("pathway") %>%
    pivot_longer(-pathway, names_to = "cell_id", values_to = "gsva_score")
  merged_data <- gsva_long %>%
    left_join(metadata %>% rownames_to_column("cell_id"), by = "cell_id")
  result_df <- merged_data %>%
    filter(seurat_clusters == cluster_id) %>%
    group_by(orig.ident, pathway) %>%
    summarise(mean_gsva = mean(gsva_score), .groups = "drop")
  result_matrix <- result_df %>%
    pivot_wider(names_from = orig.ident, values_from = mean_gsva) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  sample_avg <- colMeans(result_matrix, na.rm = TRUE)
  sample_gsvascore <- data.frame(SampleID = names(sample_avg),
                                 gsva_score = sample_avg)
  sample_merged <- inner_join(sample_Pseudotime, sample_gsvascore, by = "SampleID") %>%
    mutate(gsva_score = (gsva_score - min(gsva_score, na.rm = TRUE)) /
             (max(gsva_score, na.rm = TRUE) - min(gsva_score, na.rm = TRUE)))
  ggplot(sample_merged, aes(x = Pseudotime, y = gsva_score)) +
    geom_smooth(method = "loess", se = TRUE, level = 0.95,
                color = "blue", fill = "lightblue", linewidth = 5) +
    labs(title = paste0("Cluster ", cluster_id, " (CM04)"),
         x = "Pseudotime", y = "GSVA score") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.text = element_text(color = "black", size = 25),
          axis.title = element_text(color = "black", size = 25),
          plot.title = element_text(color = "black", size = 25, hjust = 0.5))
}

plot_list <- lapply(c(1, 6, 9), plot_cluster_gsva)
combined_plot <- wrap_plots(plot_list, ncol = 3)
ggsave(output_pdf, combined_plot, width = 24, height = 6, dpi = 300)



# =========================================================================== #
# Figure J–K: Mirrored density plots of significant serum metabolites 
# comparing Pre vs Post treatment, and heatmap of log2(Pre/Post) ratios 
# across NMF subtypes with significance stars.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(tidyverse)
library(readxl)
library(sva)
library(broom)
library(pheatmap)

# ------------------------------ File paths ----------------------------------
metab1_path        <- "path/to/metabolome/01.xlsx"
metab2_path        <- "path/to/metabolome/03.xlsx"
clinical_path      <- "path/to/Clinical_information.rds"
output_violin_pdf  <- "path/to/F6/J_violin.pdf"
output_heat_pdf    <- "path/to/F6/J_heat.pdf"

# ====================== 1. Metabolome preprocessing =========================
met01 <- read_excel(metab1_path)
met02 <- read_excel(metab2_path)

# Harmonize columns
all_cols <- union(colnames(met01), colnames(met02))
met01 <- met01 %>%
  mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>%
  select(all_of(all_cols))
met02 <- met02 %>%
  mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>%
  select(all_of(all_cols))

# Batch correction
dup_samples   <- intersect(met01$Samples, met02$Samples)
common_metabs <- intersect(names(met01)[-1], names(met02)[-1])

combined <- bind_rows(
  met01 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 1L),
  met02 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 2L)
) %>%
  mutate(BioID = ifelse(Samples %in% dup_samples, Samples,
                        paste0("uniq_", row_number())))

combat_mat <- combined %>%
  select(all_of(common_metabs)) %>%
  as.matrix() %>%
  {log10(. + 1)} %>%
  t() %>%
  ComBat(batch = combined$Batch, mod = model.matrix(~1, data = combined),
         par.prior = TRUE) %>%
  t() %>%
  {10^. - 1}
combat_mat[combat_mat < 0] <- 0

# Consolidate duplicate samples
metab_final <- combat_mat %>%
  as.data.frame() %>%
  mutate(Sample = combined$Samples, BioID = combined$BioID) %>%
  group_by(BioID) %>%
  summarise(Sample = first(Sample), across(all_of(common_metabs), mean),
            .groups = "drop") %>%
  arrange(factor(Sample, levels = c(met01$Samples,
                                    setdiff(met02$Samples, dup_samples)))) %>%
  select(Sample, all_of(common_metabs)) %>%
  column_to_rownames("Sample")

# Filter low-abundance metabolites
metab_data <- metab_final[, colSums(metab_final > 0.5) >=
                            max(5, 0.2 * nrow(metab_final))]
metab_data <- metab_data %>% mutate(across(where(is.numeric), ~ .x + 1))

# Build name mapping (original <-> coded)
mapping_list <- data.frame(
  original_name = colnames(metab_data),
  new_name = paste0("metabolite_", seq_len(ncol(metab_data)))
)
colnames(metab_data) <- mapping_list$new_name

# ====================== 2. Clinical data ====================================
group <- readRDS(clinical_path) %>%
  filter(treatment_condition %in% c("Pre", "Post"), subtype != "Unknow")

# ====================== 3. Prepare analysis data ============================
metab_data <- metab_data %>% rownames_to_column("Sample_ID")
group <- group %>% rownames_to_column("Sample_ID")

common_samples <- intersect(metab_data$Sample_ID, group$Sample_ID)

metab_data_common <- metab_data %>%
  filter(Sample_ID %in% common_samples) %>%
  arrange(Sample_ID)

group_common <- group %>%
  filter(Sample_ID %in% common_samples) %>%
  arrange(Sample_ID)

# Long format with original metabolite names
metab_long <- metab_data_common %>%
  pivot_longer(-Sample_ID, names_to = "Metabolite", values_to = "Abundance") %>%
  left_join(mapping_list %>% select(new_name, original_name),
            by = c("Metabolite" = "new_name")) %>%
  mutate(Metabolite = original_name) %>%
  select(-original_name) %>%
  left_join(group_common %>% select(Sample_ID, treatment_condition, subtype),
            by = "Sample_ID")

# Overall Pre vs Post t-test
diff_results_all <- metab_long %>%
  group_by(Metabolite) %>%
  do(tidy(t.test(Abundance ~ treatment_condition, data = .))) %>%
  ungroup() %>%
  rename(p_value = p.value) %>%
  left_join(
    metab_long %>% filter(treatment_condition == "Pre") %>%
      group_by(Metabolite) %>%
      summarise(Pre_mean = mean(Abundance, na.rm = TRUE), .groups = "drop"),
    by = "Metabolite"
  ) %>%
  left_join(
    metab_long %>% filter(treatment_condition == "Post") %>%
      group_by(Metabolite) %>%
      summarise(Post_mean = mean(Abundance, na.rm = TRUE), .groups = "drop"),
    by = "Metabolite"
  ) %>%
  mutate(
    log2FC = log2(Post_mean / Pre_mean),
    abs_log2FC = abs(log2FC),
    p_adjust = p.adjust(p_value, method = "fdr")
  )

# ====================== 4. Select metabolites ===============================
# Subtype-specific t-test to find significant metabolites
diff_by_subtype <- metab_long %>%
  group_by(Metabolite, subtype) %>%
  filter(n_distinct(treatment_condition) == 2,
         sum(treatment_condition == "Pre") >= 2,
         sum(treatment_condition == "Post") >= 2) %>%
  do(tidy(t.test(Abundance ~ treatment_condition, data = .))) %>%
  ungroup() %>%
  rename(p_value = p.value)

sig_mets <- diff_by_subtype %>%
  filter(p_value < 0.05) %>%
  distinct(Metabolite) %>%
  pull(Metabolite)

# Force include three SCFAs
target_metabolites <- c("Acetic acid", "Propanoic acid", "Butyric acid")
final_selected_mets <- unique(c(sig_mets, target_metabolites))

diff_results_filtered <- diff_results_all %>%
  filter(Metabolite %in% final_selected_mets)

target_met_names <- unique(diff_results_filtered$Metabolite)
idx <- match(target_met_names, mapping_list$original_name)
target_met_new <- data.frame(
  original_name = target_met_names,
  new_name = mapping_list$new_name[idx]
)

# ====================== 5. Mirrored density plot (Figure J) =================
metab_long_filtered <- metab_data_common %>%
  pivot_longer(-Sample_ID, names_to = "Metabolite", values_to = "Concentration") %>%
  filter(Metabolite %in% target_met_new$new_name) %>%
  left_join(group_common %>% select(Sample_ID, treatment_condition), by = "Sample_ID") %>%
  left_join(target_met_new, by = c("Metabolite" = "new_name")) %>%
  mutate(Metabolite = original_name) %>%
  select(-original_name)

# Keep metabolites with >=2 samples per group
valid_mets <- metab_long_filtered %>%
  group_by(Metabolite, treatment_condition) %>%
  summarise(n = sum(!is.na(Concentration)), .groups = "drop") %>%
  group_by(Metabolite) %>%
  filter(all(n >= 2)) %>%
  pull(Metabolite) %>% unique()

metab_long_filtered <- metab_long_filtered %>%
  filter(Metabolite %in% valid_mets)

pval_labels <- diff_results_filtered %>%
  filter(Metabolite %in% valid_mets) %>%
  mutate(p_label = ifelse(is.na(p_value), "p = NA",
                          paste0("p = ", format(p_value, digits = 3))))

peak_data <- metab_long_filtered %>%
  group_by(Metabolite, treatment_condition) %>%
  summarise(peak_x = {
    vals <- Concentration[!is.na(Concentration)]
    if (length(vals) > 1) { d <- density(vals); d$x[which.max(d$y)] }
    else NA_real_
  }, .groups = "drop") %>%
  pivot_wider(names_from = treatment_condition, values_from = peak_x) %>%
  rename(Pre_peak = Pre, Post_peak = Post)

plot_data <- metab_long_filtered %>%
  left_join(pval_labels, by = "Metabolite") %>%
  left_join(peak_data, by = "Metabolite")

p_violin <- ggplot(plot_data, aes(fill = treatment_condition)) +
  geom_density(data = subset(plot_data, treatment_condition == "Pre"),
               aes(x = Concentration, y = after_stat(density)),
               color = "black", alpha = 0.6) +
  geom_density(data = subset(plot_data, treatment_condition == "Post"),
               aes(x = Concentration, y = -after_stat(density)),
               color = "black", alpha = 0.6) +
  geom_vline(data = peak_data, aes(xintercept = Pre_peak),
             linetype = "dashed", color = "black") +
  geom_vline(data = peak_data, aes(xintercept = Post_peak),
             linetype = "dashed", color = "black") +
  geom_text(data = pval_labels,
            aes(x = Inf, y = Inf, label = p_label),
            hjust = 1.1, vjust = 2, size = 5, color = "black",
            inherit.aes = FALSE) +
  scale_fill_manual(
    name = "Condition",
    values = c("Pre" = "#dc6f84", "Post" = "#83a3d5"),
    breaks = c("Pre", "Post"),
    labels = c("Pre", "Post")
  ) +
  facet_wrap(~ Metabolite, ncol = 3, scales = "free") +
  theme_minimal() +
  theme(
    axis.line.y.left   = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y.left  = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length.y.left = unit(0.2, "cm"),
    axis.text.y.left   = element_text(color = "black", size = 15, margin = margin(r = 2)),
    axis.title.y.left  = element_text(size = 15, margin = margin(r = 5)),
    axis.line.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right  = element_blank(),
    axis.line.x        = element_line(color = "black"),
    axis.ticks.x       = element_line(color = "black"),
    axis.text.x        = element_text(color = "black", size = 15),
    axis.title.x       = element_text(size = 15),
    strip.text         = element_text(face = "bold", size = 13),
    legend.position    = "bottom",
    legend.title       = element_text(size = 13, face = "bold"),
    legend.text        = element_text(size = 12),
    legend.key.size    = unit(1.2, "lines"),
    panel.grid         = element_blank(),
    panel.spacing      = unit(1, "lines")
  ) +
  labs(x = "Concentration", y = "Density (Pre: positive, Post: negative)") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(labels = function(x) abs(x), expand = expansion(mult = 0.05))

# Save mirrored density plot
if (length(valid_mets) > 0) {
  n_rows <- ceiling(length(valid_mets) / 3)
  height_in <- max(5, 2.5 * n_rows)
  ggsave(output_violin_pdf, p_violin,
         width = 10, height = height_in, units = "in", limitsize = FALSE)
}

# ====================== 6. Pre/Post ratio heatmap (Figure K) ================
valid_new_names <- mapping_list %>%
  filter(original_name %in% valid_mets) %>% pull(new_name)

ratio_data <- metab_data_common %>%
  select(Sample_ID, all_of(valid_new_names)) %>%
  pivot_longer(-Sample_ID, names_to = "Metabolite_new", values_to = "Concentration") %>%
  left_join(group_common %>% select(Sample_ID, treatment_condition, subtype), by = "Sample_ID") %>%
  left_join(mapping_list %>% select(new_name, original_name),
            by = c("Metabolite_new" = "new_name")) %>%
  mutate(Metabolite = original_name) %>%
  select(-original_name, -Metabolite_new)

ratio_data_filtered <- ratio_data %>%
  group_by(Metabolite, subtype, treatment_condition) %>%
  filter(sum(!is.na(Concentration)) >= 2) %>%
  ungroup()

mean_by_subtype <- ratio_data_filtered %>%
  group_by(Metabolite, subtype, treatment_condition) %>%
  summarise(mean_conc = mean(Concentration, na.rm = TRUE), .groups = "drop")

mean_wide <- mean_by_subtype %>%
  pivot_wider(names_from = treatment_condition, values_from = mean_conc,
              names_prefix = "mean_") %>%
  mutate(ratio = mean_Pre / mean_Post)

ratio_matrix <- mean_wide %>%
  select(Metabolite, subtype, ratio) %>%
  pivot_wider(names_from = subtype, values_from = ratio) %>%
  column_to_rownames("Metabolite") %>%
  as.matrix()

ratio_matrix_clean <- ratio_matrix[is.finite(rowSums(ratio_matrix)), , drop = FALSE]
ratio_matrix_log2 <- log2(ratio_matrix_clean)

# Per-metabolite per-subtype t-test
pvalue_matrix <- ratio_data_filtered %>%
  group_by(Metabolite, subtype) %>%
  do({
    data <- .
    pre_vals <- data$Concentration[data$treatment_condition == "Pre"]
    post_vals <- data$Concentration[data$treatment_condition == "Post"]
    p_val <- tryCatch(t.test(pre_vals, post_vals)$p.value,
                      error = function(e) NA_real_)
    data.frame(p_value = p_val)
  }) %>%
  ungroup() %>%
  pivot_wider(names_from = subtype, values_from = p_value) %>%
  column_to_rownames("Metabolite") %>%
  as.matrix()

common_mets <- intersect(rownames(ratio_matrix_log2), rownames(pvalue_matrix))
common_subtypes <- intersect(colnames(ratio_matrix_log2), colnames(pvalue_matrix))
ratio_matrix_log2 <- ratio_matrix_log2[common_mets, common_subtypes, drop = FALSE]
pvalue_matrix <- pvalue_matrix[common_mets, common_subtypes, drop = FALSE]

star_matrix <- matrix("", nrow = nrow(pvalue_matrix), ncol = ncol(pvalue_matrix),
                      dimnames = dimnames(pvalue_matrix))
star_matrix[pvalue_matrix < 0.1]  <- "*"
star_matrix[pvalue_matrix < 0.05] <- "**"
star_matrix[pvalue_matrix < 0.01] <- "***"

col_palette <- colorRampPalette(c("#d73027", "#f7f7f7", "#4575b4"))(100)
breaks <- seq(-1.5, 1.5, length.out = 101)

pheatmap(ratio_matrix_log2,
         color = col_palette,
         breaks = breaks,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         display_numbers = star_matrix,
         fontsize_number = 10,
         number_color = "black",
         fontsize_row = 7,
         fontsize_col = 10,
         angle_col = 0,
         main = "Pre/Post Ratio (log2) Across Subtypes\n(* p<0.1, ** p<0.05, *** p<0.01)",
         filename = output_heat_pdf,
         width = 8,
         height = max(6, nrow(ratio_matrix_log2) * 0.25))



# =========================================================================== #
# Figure L: Correlation between treatment-induced changes (Post – Pre) in 
# serum metabolite levels and functional module (CM) activity, stratified 
# by NMF subtype. One multi-panel PDF is generated per metabolite.
# =========================================================================== #

# ----------------------------- Load packages --------------------------------
library(tidyverse)
library(readxl)
library(sva)
library(CoVarNet)
library(patchwork)
library(reticulate)

use_python("~/anaconda3/envs/screnv/bin/python", required = TRUE)
sc <- reticulate::import("scanpy")

# ------------------------------ File paths ----------------------------------
metab1_path        <- "path/to/metabolome/01.xlsx"
metab2_path        <- "path/to/metabolome/03.xlsx"
clinical_path      <- "path/to/Clinical_information.rds"
prop_csv_path      <- "path/to/prop_sample1_clusters.csv"
output_dir         <- "path/to/F6/L_allCM_allsubtypes"

# ====================== 1. Metabolome preprocessing =========================
met01 <- read_excel(metab1_path)
met02 <- read_excel(metab2_path)

all_cols <- union(colnames(met01), colnames(met02))
met01 <- met01 %>%
  mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>%
  select(all_of(all_cols))
met02 <- met02 %>%
  mutate(across(setdiff(all_cols, colnames(.)), ~ 0)) %>%
  select(all_of(all_cols))

dup_samples   <- intersect(met01$Samples, met02$Samples)
common_metabs <- intersect(names(met01)[-1], names(met02)[-1])

combined <- bind_rows(
  met01 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 1L),
  met02 %>% select(Samples, all_of(common_metabs)) %>% mutate(Batch = 2L)
) %>%
  mutate(BioID = ifelse(Samples %in% dup_samples, Samples,
                        paste0("uniq_", row_number())))

combat_mat <- combined %>%
  select(all_of(common_metabs)) %>%
  as.matrix() %>%
  {log10(. + 1)} %>%
  t() %>%
  ComBat(batch = combined$Batch, mod = model.matrix(~1, data = combined),
         par.prior = TRUE) %>%
  t() %>%
  {10^. - 1}
combat_mat[combat_mat < 0] <- 0

metab_final <- combat_mat %>%
  as.data.frame() %>%
  mutate(Sample = combined$Samples, BioID = combined$BioID) %>%
  group_by(BioID) %>%
  summarise(Sample = first(Sample), across(all_of(common_metabs), mean),
            .groups = "drop") %>%
  arrange(factor(Sample, levels = c(met01$Samples,
                                    setdiff(met02$Samples, dup_samples)))) %>%
  select(Sample, all_of(common_metabs)) %>%
  column_to_rownames("Sample")

# Filter low-abundance metabolites
metab_data <- metab_final[, colSums(metab_final > 0.5) >=
                            max(5, 0.2 * nrow(metab_final))]
metab_data <- metab_data %>% mutate(across(where(is.numeric), ~ .x + 1))

# Build name mapping (original <-> coded)
mapping_list <- data.frame(
  original_name = colnames(metab_data),
  new_name = paste0("metabolite_", seq_len(ncol(metab_data)))
)
colnames(metab_data) <- mapping_list$new_name

# ====================== 2. Clinical data ====================================
group <- readRDS(clinical_path) %>%
  filter(treatment_condition %in% c("Pre", "Post"), subtype != "Unknow")

# ====================== 3. NMF and CM activity ==============================
cluster_sample <- read.csv(prop_csv_path, row.names = 1)
colnames(cluster_sample) <- gsub("^X", "cluster", colnames(cluster_sample))
cluster_sample <- t(cluster_sample)

common_samples <- intersect(colnames(cluster_sample), group$sample)
cluster_sample <- cluster_sample[, common_samples, drop = FALSE]
mat_fq_raw <- as.data.frame(cluster_sample)

# Annotate functional categories
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

mat_fq_norm <- freq_normalize(mat_fq_raw, normalize = "minmax")
NMF_K <- nmf(mat_fq_norm, 4, method = "nsNMF", seed = rep(77, 6), nrun = 30)
colnames(basis(NMF_K)) <- c("CM01","CM04","CM03","CM02")
rownames(coef(NMF_K)) <- c("CM01","CM04","CM03","CM02")

cor_pair <- pair_correlation(mat_fq_raw, method = "pearson")
topnum <- floor(nrow(basis(NMF_K)) / ncol(basis(NMF_K)))

# Network construction helper (returns CM activity via scoef)
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

  listA <- Corres[, c("subCluster1","majorCluster1")]; colnames(listA) <- c("subCluster","majorCluster")
  listB <- Corres[, c("subCluster2","majorCluster2")]; colnames(listB) <- c("subCluster","majorCluster")
  ann <- rbind(listA, listB); ann <- ann[!duplicated(ann$subCluster), ]; rownames(ann) <- ann$subCluster

  w <- basis(NMFres)
  if (is.null(colnames(w))) {
    colnames(w) <- sprintf("CM%02d", seq_len(ncol(w)))
    rownames(coef(NMFres)) <- sprintf("CM%02d", seq_len(ncol(w)))
    colnames(basis(NMFres)) <- sprintf("CM%02d", seq_len(ncol(w)))
  }
  w <- w[, order(colnames(w))]

  weight <- reshape2::melt(w)
  colnames(weight) <- c("subCluster","cm","weight")
  weight <- weight %>% group_by(cm) %>% arrange(cm, desc(weight)) %>% ungroup()
  weight_filtered <- if (apply_top_n) {
    weight %>% group_by(cm) %>% slice(seq_len(max(1, top_n))) %>% ungroup()
  } else weight

  meta_sc <- data.frame(subCluster = c(Corres$subCluster1, Corres$subCluster2),
                        majorCluster = c(Corres$majorCluster1, Corres$majorCluster2))
  meta_sc <- meta_sc[!duplicated(meta_sc$subCluster), ]; rownames(meta_sc) <- meta_sc$subCluster

  pl_df <- Corres
  if (apply_corr) pl_df <- pl_df[pl_df$correlation > corr, ]
  if (apply_fdr)  pl_df <- pl_df[pl_df$pval_fdr <= fdr, ]
  if (apply_spe) {
    n_all <- nrow(meta_sc)
    spe_cutoff <- 1 - ((top_n - 1) * 2 - 1) / ((n_all - 1) * 2 - 1)
    pl_df <- pl_df[pl_df$spe >= spe_cutoff, ]
  }

  graph <- igraph::graph_from_data_frame(pl_df, directed = FALSE)
  node_global <- data.frame(subCluster = igraph::V(graph)$name,
                            majorCluster = meta_sc$majorCluster[match(igraph::V(graph)$name, meta_sc$subCluster)])
  node_each <- data.frame(); edge_each <- data.frame()

  for (cm in unique(weight_filtered$cm)) {
    cm_nodes <- weight_filtered$subCluster[weight_filtered$cm == cm]
    sub_graph <- igraph::subgraph(graph, vids = intersect(igraph::V(graph)$name, cm_nodes))
    if (remove_isolated && length(igraph::V(sub_graph)) > 0)
      sub_graph <- igraph::delete.vertices(sub_graph, igraph::V(sub_graph)[igraph::degree(sub_graph) == 0])
    if (length(igraph::V(sub_graph)) != 0) {
      tmp1 <- data.frame(cm = cm, subCluster = igraph::V(sub_graph)$name)
      node_each <- rbind(node_each, tmp1)
      tmp2 <- pl_df[(pl_df$subCluster1 %in% tmp1$subCluster) &
                    (pl_df$subCluster2 %in% tmp1$subCluster), ]
      edge_each <- rbind(edge_each, data.frame(cm = cm, tmp2))
    }
  }
  if (nrow(node_each) > 0) {
    node_each$majorCluster <- meta_sc$majorCluster[match(node_each$subCluster, meta_sc$subCluster)]
    weight_final <- merge(node_each, weight_filtered)
    weight_final <- weight_final[order(weight_final$cm, -weight_final$weight), ]
    rownames(weight_final) <- seq_len(nrow(weight_final))
  } else weight_final <- data.frame()

  list(global = list(node = node_global, edge = pl_df),
       each = list(node = node_each, edge = edge_each),
       raw = NMFres, filter = weight_final, ann = ann)
}

network <- sm_cm_network(
  NMF_K, cor_pair,
  corr = 0.2, top_n = topnum,
  apply_fdr = FALSE, apply_corr = TRUE,
  apply_top_n = TRUE, apply_spe = TRUE, remove_isolated = TRUE
)

cm_activity <- t(scoef(network$raw))

# ====================== 4. Metabolite list and loop =========================
met_list <- c("metabolite_214", "metabolite_247", "metabolite_84",
              "metabolite_278", "metabolite_259", "metabolite_288",
              "metabolite_320", "metabolite_257", "metabolite_4",
              "metabolite_81", "metabolite_181", "metabolite_2",
              "metabolite_190", "metabolite_80", "metabolite_17",
              "metabolite_221", "metabolite_129", "metabolite_254",
              "metabolite_188", "metabolite_189", "metabolite_287",
              "metabolite_40", "metabolite_134", "metabolite_329",
              "metabolite_3", "metabolite_227", "metabolite_118",
              "metabolite_207", "metabolite_61", "metabolite_192",
              "metabolite_194")

name_map <- setNames(mapping_list$original_name, mapping_list$new_name)
cm_list <- c("CM01", "CM02", "CM03", "CM04")
subtype_list <- c(1, 2, 3)

for (met in met_list) {
  real_name <- name_map[met]
  if (is.na(real_name)) next

  metabolite_data <- metab_data[, met, drop = FALSE]
  colnames(metabolite_data) <- "Metabolite"

  subtype_row_plots <- list()
  for (st in subtype_list) {
    group_filtered <- group[group$subtype == st, , drop = FALSE]
    common_ids <- Reduce(intersect, list(rownames(metabolite_data),
                                         rownames(cm_activity),
                                         rownames(group_filtered)))

    cm_plots <- list()
    for (cm in cm_list) {
      df <- data.frame(
        SampleID   = common_ids,
        Metabolite = metabolite_data[common_ids, 1],
        CM         = cm_activity[common_ids, cm],
        Time       = group_filtered[common_ids, "treatment_condition"],
        Subject    = group_filtered[common_ids, "name"]
      )

      df_filter <- df %>% filter(Time %in% c("Pre", "Post"))
      paired_subjects <- df_filter %>%
        group_by(Subject) %>%
        filter(all(c("Pre", "Post") %in% Time)) %>%
        pull(Subject) %>% unique()

      df_paired <- df_filter %>% filter(Subject %in% paired_subjects)

      p <- ggplot() + annotate("text", x = 0, y = 0,
                               label = "Not enough finite observations",
                               size = 5, color = "gray50") +
        theme_void() +
        labs(title = paste(cm, "(subtype =", st, ")")) +
        theme(plot.title = element_text(hjust = 0.5, size = 12))

      if (nrow(df_paired) > 0) {
        delta_df <- tryCatch({
          df_paired %>%
            select(Subject, Time, Metabolite, CM) %>%
            pivot_wider(names_from = Time, values_from = c(Metabolite, CM)) %>%
            mutate(Delta_Metabolite = Metabolite_Post - Metabolite_Pre,
                   Delta_CM         = CM_Post - CM_Pre) %>%
            drop_na()
        }, error = function(e) NULL)

        if (!is.null(delta_df) && nrow(delta_df) >= 3) {
          cor_delta <- cor.test(delta_df$Delta_Metabolite,
                                delta_df$Delta_CM, method = "pearson")
          label_text <- sprintf("r = %.3f, p = %.3f",
                                cor_delta$estimate, cor_delta$p.value)
          p <- ggplot(delta_df, aes(x = Delta_Metabolite, y = Delta_CM)) +
            geom_point(size = 2.5, alpha = 0.7, color = "blue") +
            geom_smooth(method = "lm", se = TRUE, level = 0.95,
                        formula = y ~ x, fill = "gray", color = "black") +
            labs(x = paste0("\u0394 ", real_name, " (Post \u2013 Pre)"),
                 y = paste0("\u0394 ", cm, " Activity (Post \u2013 Pre)"),
                 title = paste(cm, "(subtype =", st, ")")) +
            theme_classic(base_size = 12) +
            theme(panel.grid = element_blank()) +
            annotate("text", x = Inf, y = Inf,
                     label = label_text,
                     hjust = 1.1, vjust = 1.5, size = 4.5,
                     color = "black", fontface = "plain")
        }
      }
      cm_plots[[cm]] <- p
    }

    row_plot <- wrap_plots(cm_plots, ncol = 4) +
      plot_annotation(title = paste("Subtype", st),
                      theme = theme(plot.title = element_text(hjust = 0.5,
                                                              size = 14, face = "bold")))
    subtype_row_plots[[paste0("Subtype_", st)]] <- row_plot
  }

  final_plot <- wrap_plots(subtype_row_plots, ncol = 1)
  ggsave(file.path(output_dir, paste0("L_Delta_Correlation_", real_name, ".pdf")),
         plot = final_plot, device = "pdf", width = 14, height = 12)
}