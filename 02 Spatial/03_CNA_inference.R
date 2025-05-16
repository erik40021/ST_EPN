library(infercnv)
library(Seurat)
library(ggplot2)
library(parallel)
library(ComplexHeatmap)
library(igraph)
library(patchwork)
library(circlize)
library(RColorBrewer)

source("utils/r_utils.R")
source("utils/seurat_utils.R")


# ------ PART 1: infer CNAs ------

gene_positions <- read.table("gene_ordering.txt", header = T) # load human GRCh38 reference genome
cna_cols = c("#F2B701","#3969AC","#EF4868","#11A579","#6d2c6e","#ffe39e","#bd6908","#66C5CC","#ff918a","cyan2","#ffea03","#80BA5A","#D4D915","#d5bdaf",
             "#CF1C90","#4b4b8f","#B95FBB","#748cab","yellow","grey","#ebc5ae","brown","#bab475","#4e81a3","#967bad", "#542f3d","#f5bfd3","#adf590","#e3d536","#cb1adb")
names(cna_cols) = 1:length(cna_cols)

# prepare spatial cortex reference data
ref_obj1 = readRDS(file = "sstobj_#UKF256_C_ST.rds"); ref_obj2 = readRDS(file = "sstobj_#UKF265_C_ST.rds")
ref_mtx1 = GetAssayData(ref_obj1, layer = "counts"); colnames(ref_mtx1) = paste0("ref1_", colnames(ref_mtx1))
ref_mtx2 = GetAssayData(ref_obj2, layer = "counts"); colnames(ref_mtx2) = paste0("ref2_", colnames(ref_mtx2))
all_genes <- union(rownames(ref_mtx1), rownames(ref_mtx2))
ref_mtx <- matrix(0, nrow = length(all_genes), ncol = ncol(ref_mtx1) + ncol(ref_mtx2))
rownames(ref_mtx) <- all_genes; colnames(ref_mtx) <- c(colnames(ref_mtx1), colnames(ref_mtx2))
ref_mtx[rownames(ref_mtx1), colnames(ref_mtx1)] = as.matrix(ref_mtx1); ref_mtx[rownames(ref_mtx2), colnames(ref_mtx2)] = as.matrix(ref_mtx2)

options(scipen = 100)
options(bitmapType='cairo')
ht_opt$message = FALSE

input_dir = ""
out_dir = ""
sample_ids = gsub("^sstobj_|\\.rds$", "", list.files(input_dir))

# recommended to be run in parallel on server
n_samples = 12; inner_cores = 16

message("starting ", n_samples, " in parallel, using ", inner_cores, " cores each | start time: ", format(Sys.time(), "%d.%m. %X"))
res = mclapply(sample_ids, mc.cores = n_samples, function(s) { 
  run_infercnv_spatial(s, input_dir, out_dir, inner_cores, gene_positions, prepare_plot_data = T, do_plots = F, cols = cna_cols, louvain_res = 1)
  gc() 
})

# execute function below before running
# function is the same as for single-nucleus data, but additionally plots CNA scores in spatial objects
run_infercnv_spatial = function(s, input_dir, out_dir, inner_cores, gene_positions, prepare_plot_data = T, do_plots = F, cols = NULL, louvain_res = 1) {
  message("starting inferCNV run for ", s, " (", format(Sys.time(), "%d.%m. %X"), ")")
  out_dir <- file.path(out_dir, s); if (!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
  sstobj = readRDS(file = file.path(input_dir, "sstobj_", s, ".rds")) # load Seurat object of sample
  obs_mtx = GetAssayData(sstobj, layer = "counts")
  all_genes <- union(rownames(ref_mtx), rownames(obs_mtx))
  combined_mtx <- matrix(0, nrow = length(all_genes), ncol = ncol(ref_mtx) + ncol(obs_mtx))
  rownames(combined_mtx) <- all_genes; colnames(combined_mtx) <- c(colnames(ref_mtx), colnames(obs_mtx))
  combined_mtx[rownames(ref_mtx), colnames(ref_mtx)] = as.matrix(ref_mtx); combined_mtx[rownames(obs_mtx), colnames(obs_mtx)] = as.matrix(obs_mtx)
  message("Combined reference matrix (", dim(ref_mtx)[1], " x ", dim(ref_mtx)[2], ") and sample matrix (", dim(obs_mtx)[1], " x ", dim(obs_mtx)[2], ") to matrix of dimension: ", dim(combined_mtx)[1], " x ", dim(combined_mtx)[2])
  cell_annots <- data.frame(cell_id = colnames(combined_mtx), group = c(rep("ref", ncol(ref_mtx)), rep("sample", ncol(obs_mtx))))
  subset_gene_positions <- gene_positions[which(gene_positions$gene %in% rownames(combined_mtx)), ]
  write.table(cell_annots, file = file.path(out_dir, "cell_annots.txt"), sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(subset_gene_positions, file = file.path(out_dir, "gene_positions_filtered.txt"), sep = "\t", col.names = FALSE, row.names = FALSE)
  
  message("creating inferCNV object")
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = combined_mtx, ref_group_names = "ref", annotations_file = file.path(out_dir, "cell_annots.txt"), 
                                      delim="\t", gene_order_file = file.path(out_dir, "gene_positions_filtered.txt"))
  message('running infercnv')
  infercnv_obj = infercnv::run(infercnv_obj, out_dir = out_dir, cluster_by_groups=F, denoise=F, HMM=F, no_prelim_plot=T, write_phylo=F, cluster_references=F, save_rds=F,
                              min_cells_per_gene = 10, cutoff = 0.1, analysis_mode = 'subclusters', tumor_subcluster_pval = 0.05, leiden_resolution = 0.001,
                              num_threads = inner_cores, output_format = "png", window_length = 101)
  message('finished infercnv run')

  if (prepare_plot_data) { 
    cna_mtx = infercnv_obj@expr.data # cna prediction matrix
    cna_mtx = cna_mtx[, -unlist(infercnv_obj@reference_grouped_cell_indices)] # remove all reference cells

    malignancy_score = calculate_malignancy_score_by_CNAs(cna_mtx, print_stats = TRUE, normalise = F)
    saveRDS(malignancy_score, file = file.path(out_dir, "malignancy_score.rds"))
    malignancy_score_norm = calculate_malignancy_score_by_CNAs(cna_mtx, print_stats = TRUE, normalise = T)
    saveRDS(malignancy_score_norm, file = file.path(out_dir, "malignancy_score_norm.rds"))

    leiden_clusters = lapply(infercnv_obj@tumor_subclusters$subclusters$all_observations, function(clust) colnames(infercnv_obj@expr.data)[clust])
    names(leiden_clusters) = sub("all_observations_s", "", names(leiden_clusters))
    ord_leiden_clusters = order_clusters_by_CNA_malignancy(malignancy_score, leiden_clusters)
    saveRDS(ord_leiden_clusters, file = file.path(out_dir, "ordered_cna_leiden_clusters.rds"))

    message("starting Louvain clustering for matrix [", dim(cna_mtx)[1], " x ", dim(cna_mtx)[2], "]")
    col_dist = dist(t(cna_mtx))
    cell_similarity <- 1 / as.matrix(col_dist)  # Inverting distances to similarities
    cell_similarity[is.infinite(cell_similarity)] <- 0 # Replace Inf (which may result from division by zero) with 0
    g <- graph_from_adjacency_matrix(cell_similarity, mode = "undirected", weighted = TRUE, diag = FALSE) # Create graph from similarity matrix
    set.seed(42); louvain_result <- cluster_louvain(g, resolution = louvain_res)
    l_clusters = louvain_result$membership; names(l_clusters) = louvain_result$names
    message(length(unique(l_clusters)), " Louvain clusters identified, with distribution: ", paste0(table(l_clusters), sep = " "))
    ord_louvain_clusters = order_clusters_by_CNA_malignancy(malignancy_score, l_clusters)
    saveRDS(ord_louvain_clusters, file = paste0(out_dir, "/ordered_cna_louvain-res", louvain_res, "_clusters.rds")) # save clusters ordered by malignancy score for later 

    # Transfer identified CNA-derived clusters to tissue space (same as in 'st_annotate_malignancy.R')
    sstobj$cna_malignancy = malignancy_score
    sstobj$leiden_clusters = ord_leiden_clusters
    sstobj$louvain_clusters = ord_louvain_clusters
    ratio = get_spatial_aspect_ratio(sstobj); pt_size = get_spatial_point_size(sstobj, scale_factor = 3.5) # spatial plot params
    p1 = ggplot(data.frame(Value = malignancy_score), aes(x = Value)) + geom_density(alpha = 0.5) + labs(x = "Malignancy Score", y = "Spot Density") + 
        annotate("text", x = mean(malignancy_score), y = 0, label = sprintf("SD: %.2f", sqrt(var(malignancy_score))), hjust = 0.5, vjust = -2, color = "black", size = 4) + theme_minimal()
    p2 = SpatialFeaturePlot(sstobj, features = "cna_malignancy", pt.size.factor = pt_size, image.alpha = 0) + theme(legend.position = "right", aspect.ratio = ratio)
    ggsave("04a_cna_malignancy_in_tissue_+distribution.png", plot = p2 + p1 + plot_layout(widths = c(2, 1)), path = out_dir, width = 13, height = 7)
    p1 = SpatialDimPlot(sstobj, group.by = "leiden_clusters", pt.size.factor = pt_size, image.alpha = 0, cols = cols) + theme(aspect.ratio = ratio)
    p2 = SpatialDimPlot(sstobj, group.by = "louvain_clusters", pt.size.factor = pt_size, image.alpha = 0, cols = cols) + theme(aspect.ratio = ratio)
    ggsave("04b_cna_clusters_in_tissue.png", plot = p1 + p2, path = out_dir, width = 13, height = 7)
    if (do_plots) {
      # 1. Default order
      png(filename = file.path(out_dir, "01_heatmap_CNAs_default-order.png"), width = 2600, height = 1600, res = 300)
      print(heatmap_cna(t(cna_mtx), gene_positions, row_title = "Spots in default order"))
      dev.off()

      # 2. Order by CNA-derived malignancy score
      reordered_expr_data = cna_mtx[, order(malignancy_score, decreasing = TRUE)]
      png(filename = file.path(out_dir, "02_heatmap_CNAs_ordered-by-mscore.png"), width = 2600, height = 1600, res = 300)
      print(heatmap_cna(t(reordered_expr_data), gene_positions, row_title = "Spots ordered by CNA expression"))
      dev.off()

      # 3. Order by and within clusters
      # a) using existing leiden clusters
      by_clusters_reordered_expr_data = cna_mtx[, names(ord_leiden_clusters)]
      ord_leiden_clusters = factor(ord_leiden_clusters, levels = unique(ord_leiden_clusters))
      cluster_anno <- HeatmapAnnotation(seurat_cluster = ord_leiden_clusters, col = list(seurat_cluster = cna_cols), which = "row", show_annotation_name = F)
      png(filename = file.path(out_dir, "03a_heatmap_CNAs_leiden_clusters_ordered-by-mscore.png"), width = 2600, height = 1600, res = 300)
      print(heatmap_cna(t(by_clusters_reordered_expr_data), gene_positions, cluster_splits = ord_leiden_clusters, row_anno = cluster_anno, row_title = "Spots ordered by CNA expression"))
      dev.off()
      # b) using new unbiased Louvain clustering (better distribution of clusters, no singletons!)
      by_clusters_reordered_expr_data = cna_mtx[, names(ord_louvain_clusters)]
      ord_louvain_clusters = factor(ord_louvain_clusters, levels = unique(ord_louvain_clusters))
      cluster_anno <- HeatmapAnnotation(cna_cluster = ord_louvain_clusters, col = list(cna_cluster = cols), which = "row", show_annotation_name = F)
      png(filename = file.path(out_dir, "03b_heatmap_CNAs_louvain-res", louvain_res, "_ordered-by-mscore.png"), width = 2600, height = 1600, res = 300)
      print(heatmap_cna(t(by_clusters_reordered_expr_data), gene_positions, cluster_splits = ord_louvain_clusters, row_anno = cluster_anno, row_title = "Spots ordered by CNA expression"))
      dev.off()
    }
  }
}


# ------ PART 2: correlate CNAs to NMF programs ------

metaprograms = readxl::read_excel("spatial_NMF_metaprograms_main.xlsx", sheet = "MP genes (final)")
# extract and prepare data from every inferCNV object (only first time)
raw_mscore_cors_per_mp = list()
norm_mscore_cors_per_mp = list()
avg_profile_cors_per_mp = list()
for (s in sample_ids) {
  message(">>> loading sample '", s, "' at [", format(Sys.time(), "%d.%m. %X"), "] <<<")
  sstobj = readRDS(paste0("sstobj_", s, ".rds"))
  infercnv_obj = readRDS(file.path(s, "run.final.infercnv_obj"))
  norm_malignancy_score = readRDS(file.path(s, "malignancy_score_norm.rds"))
  raw_malignancy_score = readRDS(file.path(s, "malignancy_score.rds"))
  cna_mtx = infercnv_obj@expr.data # cna prediction matrix
  cna_mtx = cna_mtx[, colnames(cna_mtx) %in% colnames(sstobj)]
  # calculate correlation of MP scores and CNA malignancy score
  norm_mscore_cors_per_mp[[s]] = sapply(1:ncol(metaprograms), function(i) cor(sstobj[[paste0("MP", i, "_mp-norm")]][names(norm_malignancy_score), ], norm_malignancy_score))
}
saveRDS(norm_mscore_cors_per_mp, file = "Data/norm_malignancy_vs_mp_scores_correlation_per_sample.rds")

# relate normalised malignancy score to MPs
norm_mscore_cors_per_mp = readRDS("norm_malignancy_vs_mp_scores_correlation_per_sample.rds")
data <- do.call(rbind, norm_mscore_cors_per_mp); colnames(data) = colnames(metaprograms)
data <- reshape2::melt(data); colnames(data) <- c("sample", "metaprogram", "correlation")
mean_scores = data %>% group_by(metaprogram) %>% summarise(mean_cor = mean(correlation, na.rm = TRUE)) %>% arrange(desc(mean_cor))
data <- data %>% mutate(metaprogram = factor(metaprogram, levels = mean_scores$metaprogram))
data$subtype = metadata$Subtype[match(data$sample, metadata$`Study ID`)]

# a. ggpubr version:
ggviolin(data, x = "metaprogram", y = "correlation", add = "boxplot", fill = "metaprogram", palette = st_mp_cols, add.params = list(fill = "white")) + labs(x = "metaprogram", y = "correlation to malignancy score") # + stat_compare_means(method = "t.test", label = "p.format", comparisons = NULL) # Global significance across groups
ggsave(filename = "01a_violin_norm-malignancy-score_vs_MP_score_correlation.png", width = 9, height = 5, path = out_dir)
# # split by subtype:
# ggviolin(data, x = "metaprogram", y = "correlation", add = "boxplot", fill = "metaprogram", palette = cna_cols, add.params = list(fill = "white")) + 
#   labs(x = "metaprogram", y = "correlation to malignancy score") + # stat_compare_means(method = "t.test", label = "p.format", comparisons = NULL) # Global significance across groups
#   facet_wrap(~subtype)
# ggsave(filename = "01b_violin_malignancy-score_vs_MP_score_correlation_split-by-subtype.png", width = 12, height = 8, path = out_dir)
# b. ggbetweenstats version:
ggbetweenstats(data, metaprogram, correlation, type = "parametric", pairwise.display = "none", p.adjust.method = "holm", mean.plotting = TRUE,
               mean.ci = TRUE, messages = FALSE, violin.args = list(width = 1, alpha = 0), boxplot.args = list(alpha = 0, width = 0.2), 
               point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 1), alpha = 1, size = 1.5)) + theme_classic() +
  scale_color_manual(values = st_mp_cols) + labs(title = "Violin Plot of Mscore by Metaprogram", x = "metaprogram", y = "correlation to malignancy score")
ggsave(filename = "01b_violin_norm-malignancy-score_vs_MP_score_correlation.png", width = 12, height = 8, path = out_dir)





                             
