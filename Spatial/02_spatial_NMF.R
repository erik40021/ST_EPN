library(Seurat)
library(writexl)
library(NMF)
library(patchwork)

source("utils/nmf_additional_functions.R")
source("utils/spatial_utils.R")

s = commandArgs(trailingOnly = TRUE)[1] # ID of sample to run NMF for (first command line argument)

in_path = "input_data/sst_objects_stepn"
out_path = "output_data/NMF/spatial"

options(bitmapType = 'cairo')
range = 2:20
n_features = 5000

message("Started spatial NMF pipeline for sample ", s)
out_dir = file.path(out_path, s); if (!dir.exists(out_dir)) { dir.create(out_dir, r = T) }

# 1. ---- load and prepare data -----
sstobj = readRDS(file = paste0(in_path, "/sstobj_", s, ".rds"))
data = as.matrix(GetAssayData(sstobj, layer = 'scale.data'))
data = data[VariableFeatures(sstobj)[1:n_features], ]
data[data < 0] = 0  # set negative values to zero
data = data[apply(data, 1, var) > 0,] # subset on genes with a variance above zero
message("[", s, "] Extracted matrix of dimensions ", dim(data)[1], " (genes) x ", dim(data)[2], " (spots)")

# 2. ---- run NMF -----
message("running NMF for ranks ", range[1], "-", tail(range, 1), " | start time: ", format(Sys.time(), "%d.%m. %X"))
res.list = lapply(range, function(r) {
    message(r, " (", format(Sys.time(), "%X"), ")")
    nmf(data, rank = r, nrun = 1, seed = "ica", method = "nsNMF")
})
names(res.list) = range
saveRDS(res.list, file = paste0(out_dir, "/raw_res-list_", s, ".rds"))

# 3. ---- select optimal rank -----
modules.list = lapply(res.list, NMFToModules, gmin = 5)
comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
print(comp)
mi = min(comp)
r = names(which(comp == mi))
r = r[length(r)] # chooses the rank for which the next bigger rank does not add another "useful" (gmin >= 5) module
message("selecting rank ", r)
res = res.list[[r]]
saveRDS(res, file = paste0(out_dir, "/res_final-rank_", s, ".rds")) # saves selected result

# 4. ---- calculate signature genes of each module, and add modules to object -----
tryCatch(expr = {
  modules = NMFToModules(res, gmin = 5) # finds the top associated genes of each module
  gene_scores = basis(res); colnames(gene_scores) = names(modules)
  cell_scores = coefficients(res); rownames(cell_scores) = names(modules)
  # order modules by euclidean distance
  h = Heatmap(cell_scores, clustering_distance_columns = "euclidean") # using heatmap only to cluster modules
  o = row_order(h)
  gene_scores = gene_scores[, o]; cell_scores = cell_scores[o, ]; modules = modules[o] # reorder all data
  sstobj = AddMetaData(sstobj, t(cell_scores), col.name = rownames(cell_scores))
  # assign each cell to a module using clustering
  h = Heatmap(cell_scores, clustering_distance_columns = "euclidean") 
  hcl = as.hclust(column_dend(h))
  sig = cutree(hcl, k = length(modules))
  nmf = c(by(t(cell_scores), INDICES = sig, FUN = function(x){names(modules)[which.max(colMeans(x))]}))[sig] # assigns NMF programs based on maximal expression per cell/spot
  sstobj@meta.data$nmf = factor(nmf, levels = names(modules))
  col_nmf = c(brewer.pal(12, "Set3"), brewer.pal(8, "Set1"))[1:nlevels(sstobj$nmf)]
  names(col_nmf) = levels(sstobj$nmf)
}, error = function(e){c()})
# saveRDS(sstobj, file = paste0(out_dir, "/post_NMF_sobj_r", r, "_", s, ".rds")) # save seurat object with nmf modules

# saving the modules to an excel file
modules_df <- as.data.frame(t(do.call(rbind, lapply(modules, `length<-`, max(sapply(modules, length)))))) # convert list to dataframe
write_xlsx(modules_df, path=paste0(out_dir, "/modules_of_rank", r, "_", s, ".xlsx"))

# 5. ---- plot modules vs. clusters of sobj -----
col_cluster = hue_pal()(nlevels(sstobj@meta.data$seurat_clusters)) # saving the cluster colors right away
names(col_cluster) = levels(sstobj@meta.data$seurat_clusters)
for (reduction in c("pca", "umap")) {
  p1 = DimPlot(sstobj, pt.size = 2, reduction = reduction, group.by = "seurat_clusters", cols = col_cluster, label = TRUE, label.size=6)
  p2 = DimPlot(sstobj, pt.size = 2, reduction = reduction, group.by = "nmf", cols = col_nmf)
  p3 = FeaturePlot(sstobj, reduction = reduction, features = names(modules))
  ggsave(paste0(out_dir, "/", reduction, "_", s, "_for_rank", r, ".png"), plot = wrap_plots(p1, p2, p3, nrow = 1), height = 15, width = 30)
}

# 6. ---- plot heatmaps with top DEGs and save GO-annotations ----
pdf(paste0(out_dir, "/heatmaps_", s, "_for_rank", r, ".pdf"), height = 20, width = 10)
# Clusters
markers_clusters = FindAllMarkers2(sstobj, do.plot = TRUE, enrichment.type = "GO", group.by = "seurat_clusters", cols = col_cluster, print.bar = FALSE, do.enrichment=F)
# NMF
markers_nmf = FindAllMarkers2(sstobj, do.plot = TRUE, enrichment.type = "GO", group.by = "nmf", cols = col_nmf, print.bar = FALSE, do.enrichment=F)
# Modules
top_ann = HeatmapAnnotation(df = data.frame("seurat_clusters" = sstobj@meta.data$seurat_clusters, "nmf" = sstobj@meta.data$nmf), col = list("seurat_clusters" = col_cluster, "nmf" = col_nmf), which = "column")
h = Heatmap(cell_scores, name = "cell_scores", top_ann = top_ann, show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, 
            column_order = order(sstobj@meta.data$nmf, sstobj@meta.data$seurat_clusters), breaks = c(0,max(cell_scores)/2), colors = c("white","red"))
print(h)
h = Heatmap(gene_scores[unlist(modules),], name = "gene_scores", show_row_names = TRUE, row_names_gp = gpar(cex = 0.5), cluster_rows = FALSE, cluster_columns = FALSE,
            split = factor(unlist(mapply(rep, names(modules), sapply(modules, length))), levels = names(modules)))
print(h)
dev.off()

q("no") # terminate without saving
