library(Seurat)
library(ggplot2)
library(patchwork)
library(pals)
library(RColorBrewer)
library(readxl)

utils_dir = ""
cloud_path = ""
source(file.path(utils_dir, "r_utils.R"))
source(file.path(utils_dir, "Single cell/seurat_utils.R"))
source(file.path(utils_dir, "Spatial/spatial_utils.R"))


metadata = readxl::read_excel("masterlist_EPN-ST.xlsx", skip = 1); metadata = metadata[grepl("2", metadata$`Analyses*`), ]
sample_ids = metadata$`Study ID`

# --------------------------------------- PART 1: create seurat objects and filter spots --------------------------------------------

# ------ key processing & filtering parameters ------
min_counts_per_spot = 300                 # minimal number of counts a spot must have (filtering)
min_counts_per_spot_necrotic = 100        # same but for necrotic spots
max_mt_pct = 20                           # maximal mitochondrial percentage a spot can have (filtering)
max_mt_pct_necrotic = 50                  # same but for necrotic spots
dims = 30                                 # dimensions to use for PCA and UMAP
res = 0.5                                 # resolution to use for clustering
n_var_features = 7000                     # number of most variable features to compute
necrosis_genes = c("VEGFA", "CXCL8", "HILPDA", "ADM", "CA12", "CA9", "NDRG1", "HK2", "TREM1", "ANKRD37", "ANGPTL4", "C15orf48", "IGFBP5", "EGLN3", "SERPINE1", "PLOD2") # genes taken from http://glioblastoma.alleninstitute.org ("Pseudopalisading cells around necrosis")
necrosis_threshold = 0.01
smoothing_win_size = 3
# --------------------------------------------------

options(future.globals.maxSize = 1000 * 1024^2) # increase RAM limit for SCTransform function for big samples
for (s in sample_ids) {
  message(">>> starting Seurat pipeline for sample '", s, "' at [", format(Sys.time(), "%d.%m. %X"), "] <<<")
  raw_data_dir = file.path(cloud_path, metadata[metadata$`Study ID` == s, ]$`KK code`, "outs")
  out_dir = file.path("Samples", s)
  degs_dir = file.path(out_dir, "DEGs"); if(!file.exists(degs_dir)) dir.create(degs_dir, recursive = T)
  
  sstobj = Load10X_Spatial(raw_data_dir, filename = "filtered_feature_bc_matrix.h5", slice = "tissue_lowres_image")
  sstobj@project.name = s; sstobj$orig.ident = s; Idents(sstobj) = sstobj$orig.ident
  
  ratio = get_spatial_aspect_ratio(sstobj); pt_size = get_spatial_point_size(sstobj, scale_factor = 3.5) # spatial plot params
  p = SpatialFeaturePlot(sstobj, features = NULL, pt.size.factor = pt_size, alpha = 0.2) + theme(aspect.ratio = ratio) + NoLegend() # plot H&E image without spots
  ggsave("00_spots_in_tissue.png", plot = p, path = out_dir, width = 10, height = 10)
  
  # 1. QC + filtering spots by counts and mt_pct
  p1 = VlnPlot(sstobj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  p2 = SpatialFeaturePlot(sstobj, features = "nCount_Spatial", pt.size.factor = pt_size) + theme(aspect.ratio = ratio)
  ggsave("01a_QC_counts.png", plot = wrap_plots(p1, p2), path = out_dir, width = 15, height = 7)
  p1 = VlnPlot(sstobj, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
  p2 = SpatialFeaturePlot(sstobj, features = "nFeature_Spatial", pt.size.factor = pt_size) + theme(aspect.ratio = ratio)
  ggsave("01b_QC_features.png", plot = wrap_plots(p1, p2), path = out_dir, width = 15, height = 7)
  sstobj$percent_mt = PercentageFeatureSet(sstobj, pattern = "^MT-") # filter spots with high percent (max_mt_pct) mitochondrial genes
  sstobj$percent_mt[is.na(sstobj$percent_mt)] = 0
  p1 = VlnPlot(sstobj, features = "percent_mt", pt.size = 0.1) + NoLegend()
  p2 = SpatialFeaturePlot(sstobj, features = "percent_mt", pt.size.factor = pt_size) + theme(aspect.ratio = ratio)
  ggsave("01c_QC_mt-pct.png", plot = wrap_plots(p1, p2), path = out_dir, width = 15, height = 7)
  
  # score for necrotic marker genes before filtering. Apply weaker thresholds for high-necrosis-scoring spots
  sstobj = subset(sstobj, subset = nCount_Spatial > 0)
  sstobj = SCTransform(sstobj, assay = "Spatial") # only run for scoring necrosis: normalise genes, find most variable genes, and scale those genes
  sstobj = try_add_module_score(sstobj, features = list(necrosis_genes), name = "necrosis", ctrl = 100, min_ctrl = 50, nbin = 24, min_bin = 14, verbose = T)
  necrotic <- ifelse(sstobj$necrosis1 >= necrosis_threshold, "pot.necrotic", "non-necrotic")
  
  # smoothen necrotic regions:
  message("smoothing of zones in ", s, " using 'smoothing_win_size' = ", smoothing_win_size)
  spots_positions = read.csv(list.files(file.path(raw_data_dir, "spatial"), pattern = "tissue_positions.*\\.csv$", full.names = T), header = F, skip = 1)
  neighbors_table = neighbors_table_funcV2(spots_positions, data.frame(necrosis=necrotic,barcodes=names(necrotic))); neighbors_table[neighbors_table== "NaN"] = NA
  necrotic_smo = necrotic
  for (spot in colnames(sstobj)) {
    win_spots = c(spot); for (j in 1:smoothing_win_size) win_spots = unique(c(win_spots, unique(na.omit(as.character(neighbors_table[win_spots, ])))))
    win_necrosis = necrotic[names(necrotic) %in% win_spots]
    necrotic_smo[spot] = ifelse(length(win_necrosis[win_necrosis == "pot.necrotic"])/length(win_necrosis) >= 0.9, "necrotic", 
                                   ifelse(length(win_necrosis[win_necrosis == "pot.necrotic"])/length(win_necrosis) >= 0.5, "pot.necrotic", "non-necrotic"))
  }
  sstobj$necrotic = necrotic_smo
  
  p1a = SpatialFeaturePlot(sstobj, features = "necrosis1", pt.size.factor = pt_size) + theme(aspect.ratio = ratio, legend.position = "right") + labs(fill = "necrosis")
  p1b = SpatialDimPlot(sstobj, group.by = "necrotic", pt.size.factor = pt_size) + theme(aspect.ratio = ratio) + scale_fill_manual(values = c(necrotic="#a82203", pot.necrotic="#f1af3a", "non-necrotic"="#208cc0"))
  p2 = SpatialFeaturePlot(sstobj, features = necrosis_genes, pt.size.factor = pt_size) & theme(aspect.ratio = ratio)
  ggsave("01d_necrosis.png", plot = ((p1a+p1b+plot_layout(ncol=1))|(p2+plot_layout(nrow = 3))) + plot_layout(widths = c(1,2)), path = out_dir, width = 15, height = 7)
  
  p1 = SpatialFeaturePlot(sstobj, features = "nCount_Spatial", pt.size.factor = pt_size, max.cutoff = 2000) + theme(aspect.ratio = ratio) + ggtitle("Before (max counts capped at 2000)")
  p2 = SpatialFeaturePlot(sstobj, features = "percent_mt", pt.size.factor = pt_size, max.cutoff = max_mt_pct_necrotic) + theme(aspect.ratio = ratio) + ggtitle(paste0("Before (max mt.pct capped at ", max_mt_pct_necrotic, "%)"))
  sstobj = subset(sstobj, subset = (necrotic == "non-necrotic" & percent_mt < max_mt_pct & nCount_Spatial >= min_counts_per_spot) | 
                                   (necrotic %in% c("necrotic", "pot.necrotic") & percent_mt < max_mt_pct_necrotic & nCount_Spatial >= min_counts_per_spot_necrotic))
  ratio = get_spatial_aspect_ratio(sstobj); pt_size = get_spatial_point_size(sstobj, scale_factor = 3.5) # update after filtering
  p3 = SpatialFeaturePlot(sstobj, features = NULL, pt.size.factor = pt_size) + NoLegend() + theme(aspect.ratio = ratio) + 
  ggtitle(paste0("After filtering\n(features >= ", min_counts_per_spot, "/", min_counts_per_spot_necrotic, " & mt_pct < ", max_mt_pct, "/", max_mt_pct_necrotic, "%)"))
  ggsave("01e_QC_post_filtering.png", plot = wrap_plots(p1, p2, p3), path = out_dir, width = 15, height = 7)
  
  # 2. Standard processing (transformation + clustering)
  sstobj = SCTransform(sstobj, assay = "Spatial", variable.features.n = n_var_features) # final run: normalise genes, find most variable genes, and scale those genes
  sstobj = RunPCA(sstobj, verbose = FALSE, features = VariableFeatures(sstobj))
  
  sstobj = FindNeighbors(sstobj, reduction = "pca", dims = 1:dims)
  sstobj = FindClusters(sstobj, verbose = FALSE, resolution = res)
  sstobj = RunUMAP(sstobj, reduction = "pca", dims = 1:dims)
  p1 = DimPlot(sstobj, reduction = "umap", label = TRUE, cols = "Set1")
  p2 = SpatialDimPlot(sstobj, label = TRUE, label.size = 5, pt.size.factor = pt_size, cols = "Set1") + NoLegend() + theme(aspect.ratio = ratio)
  ggsave(paste0("02a_clusters-res", res, "_in_umap_and_tissue.png"), plot = p1 + p2, path = out_dir, width = 20, height = 10)
  p = SpatialDimPlot(sstobj, cells.highlight = CellsByIdentities(object = sstobj), facet.highlight = TRUE, pt.size.factor = pt_size, 
                     ncol = c(1,2,3,3,3,3,4,4,5,5,rep(6, 20))[length(levels(sstobj$seurat_clusters))]) & theme(aspect.ratio = ratio)
  ggsave(paste0("02b_clusters-res", res, "_in_tissue_split.png"), plot = p, path = out_dir, width = 20, height = 10)
  p1 = VlnPlot(sstobj, features = "nFeature_SCT", pt.size = 0.1) + NoLegend()
  p2 = SpatialDimPlot(sstobj, label = TRUE, label.size = 5, pt.size.factor = pt_size, cols = "Set1") + NoLegend() + theme(aspect.ratio = ratio)
  ggsave(paste0("02c_clusters-res", res, "_vs_features.png"), plot = wrap_plots(p1, p2), path = out_dir, width = 20, height = 10)

  # 4. Save fully processed object
  saveRDS(sstobj, file = paste0("Objects/sstobj_", s, ".rds"))
}

