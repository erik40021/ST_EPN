library(infercnv)
library(Seurat)
library(parallel)
library(ComplexHeatmap)
library(igraph)

source("utils/r_utils.R")

sobj_path = "ST_EPN.rds"
out_dir = "inferCNV"
sobj = readRDS(file = sobj_path)
sample_ids = unique(sobj@meta.data[["orig.ident"]])[order(unique(sobj@meta.data[["orig.ident"]]))] # ordered alphabetically
gene_positions = read.table("input_data/misc/gene_ordering.txt", header = T)
cna_cols = c("#F2B701","#3969AC","#EF4868","#11A579","#6d2c6e","#ffe39e","#bd6908","#66C5CC","#ff918a","cyan2","#ffea03","#80BA5A","#D4D915","#d5bdaf",
             "#CF1C90","#4b4b8f","#B95FBB","#748cab","yellow","grey","#ebc5ae","brown","#bab475","#4e81a3","#967bad", "#542f3d","#f5bfd3","#adf590","#e3d536","#cb1adb")
names(cna_cols) = 0:(length(cna_cols)-1)

options(scipen = 100)
options(bitmapType='cairo')
ht_opt$message = FALSE


ref_clusters <- c(1, 20, 25, 28) # high confidence TME clusters

# recommended to be run in parallel on server
n_samples = 12; inner_cores = 16

message(length(sample_ids), " samples found in object (", toString(sample_ids), ")")
message("calculating ", n_samples, " samples in parallel, using ", inner_cores, " cores each | start time: ", format(Sys.time(), "%d.%m. %X"))
res = mclapply(sample_ids, mc.cores = n_samples, function(s) {
  subset_obj = subset(sobj, subset = orig.ident == s)
  run_infercnv_sc(s, subset_obj, out_dir, ref_clusters, inner_cores, gene_positions, prepare_plot_data = T, do_plots = F, cols = cna_cols, louvain_res = 1)
  gc()
})
message("finished all samples | end time: ", format(Sys.time(), "%d.%m. %X"))

# execute function below before running
# runs inferCNV on a given single-cell subset obj representing one sample. Optionally plots additional CNA-derived clusterings. 
run_infercnv_sc = function(s, subset_sobj, out_dir, ref_clusters, inner_cores, gene_positions, prepare_plot_data = T, do_plots = F, cols = NULL, louvain_res = 1) {
  message("starting inferCNV run for ", s, " (", format(Sys.time(), "%d.%m. %X"), ")")
  out_dir <- file.path(out_dir, s); if(!dir.exists(out_dir)) { dir.create(out_dir, r = T) }

  valid_clusters <- names(table(Idents(subset_sobj)))[which(table(Idents(subset_sobj)) > 10)] # at least 10 cells per cluster
  subset_sobj <- subset(subset_sobj, idents = valid_clusters)
  
  cell_annots <- data.frame(cell_id = rownames(subset_sobj[[]]), annotation = Idents(subset_sobj))
  rownames(cell_annots) <- 1:nrow(cell_annots)
  subset_gene_positions <- gene_positions[which(gene_positions$gene %in% Features(subset_sobj)),]
  write.table(cell_annots, file = file.path(out_dir,"cell_annots.txt"), sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(subset_gene_positions, file = file.path(out_dir,"gene_positions_filtered.txt"), sep = "\t", col.names = FALSE, row.names = FALSE)
  ref_clusters <- names(table(Idents(subset_sobj)))[which(table(Idents(subset_sobj)) > 10 & names(table(Idents(subset_sobj))) %in% ref_clusters)]
  
  message("creating inferCNV object")
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = subset_sobj[["RNA"]]$counts, annotations_file = file.path(out_dir, "cell_annots.txt"), delim="\t",
                                      gene_order_file = file.path(out_dir, "gene_positions_filtered.txt"), ref_group_names = ref_clusters)
  message('running infercnv')
  infercnv_obj = infercnv::run(infercnv_obj, out_dir = out_dir, denoise=F, HMM=F, no_prelim_plot=T, write_phylo=F, cluster_references=F, save_rds=F,
                              min_cells_per_gene = 10, cutoff = 0.1, analysis_mode = 'subclusters', tumor_subcluster_pval = 0.05, leiden_resolution = 0.001,
                              num_threads = inner_cores, output_format = "png", window_length = 101, cluster_by_groups=F)

  message('finished infercnv run')
  if (prepare_plot_data) {
    cna_mtx = infercnv_obj@expr.data # cna prediction matrix
    cna_mtx = cna_mtx[, -unlist(infercnv_obj@reference_grouped_cell_indices)] # remove all reference cells

    malignancy_score = calculate_malignancy_score_by_CNAs(cna_mtx, print_stats = TRUE, normalise = F)
    saveRDS(malignancy_score, file = file.path(out_dir, "malignancy_score.rds"))
    malignancy_score_norm = calculate_malignancy_score_by_CNAs(cna_mtx, print_stats = TRUE, normalise = T)
    saveRDS(malignancy_score_norm, file = file.path(out_dir, "malignancy_score_norm.rds"))

    leiden_clusters = lapply(infercnv_obj@tumor_subclusters$subclusters$all_observations, function(clust) colnames(infercnv_obj@expr.data)[clust]) # makes sense only when using cluster_by_groups = F
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
      png(filename = paste0(out_dir, "/03b_heatmap_CNAs_louvain-res", louvain_res, "_ordered-by-mscore.png"), width = 2600, height = 1600, res = 300)
      print(heatmap_cna(t(by_clusters_reordered_expr_data), gene_positions, cluster_splits = ord_louvain_clusters, row_anno = cluster_anno, row_title = "Spots ordered by CNA expression"))
      dev.off()
    }
  }
}

