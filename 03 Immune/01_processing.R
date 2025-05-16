library(Seurat)
library(ggplot2)

source("utils/r_utils.R")
source("utils/seurat_utils.R")

sobj = readRDS(file = "ST_EPN.rds")
immu_sobj = subset(sobj, idents = c(1,20,25,28)) # immune cluster idents

base_dir = ""


# --- 1. re-do default pre-processing on new subset of cells ---

immu_sobj <- NormalizeData(immu_sobj)
immu_sobj <- FindVariableFeatures(immu_sobj)
immu_sobj <- ScaleData(immu_sobj)
immu_sobj <- RunPCA(immu_sobj, verbose = FALSE)

pcs_chosen = 30; cluster_resolution = 0.6
immu_sobj <- FindNeighbors(immu_sobj, reduction = "pca", dims = 1:pcs_chosen)
immu_sobj <- FindClusters(immu_sobj, resolution = cluster_resolution, cluster.name = "unintegrated_clusters")
immu_sobj = RunUMAP(immu_sobj, reduction = "pca", dims = 1:pcs_chosen, reduction.name = "umap.unintegrated")


# --- 2. plot new NON-integration based clustering and UMAP ---

# plot clusters vs samples
p = DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "orig.ident", label = T, raster = F)
ggsave(filename = paste0("01a_umap_unintegrated_clusters-vs-samples_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(non_int_dir, "Clustering"))
# plot clusters vs patients
p = DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "patients", label = T, raster = F)
ggsave(filename = paste0("01b_umap_unintegrated_clusters-vs-patients_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(non_int_dir, "Clustering"))
# plot clusters vs batches
p = DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "batch", raster = F)
ggsave(filename = paste0("02_umap_unintegrated_clusters-vs-batches_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(non_int_dir, "Clustering"))
# plot clusters vs subtypes
p = DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "subtype", raster = F)
ggsave(filename = paste0("03a_umap_unintegrated_clusters-vs-subtypes_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(non_int_dir, "Clustering"))
# plot clusters vs subclasses
p = DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "subclass", raster = F)
ggsave(filename = paste0("03b_umap_unintegrated_clusters-vs-subclass_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(non_int_dir, "Clustering"))
# plot clusters vs subclasses (bc)
p = DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "subclass_bc", raster = F)
ggsave(filename = paste0("03c_umap_unintegrated_clusters-vs-subclass-bc_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(non_int_dir, "Clustering"))
# plot clusters vs ag
p = DimPlot(immu_sobj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  FeaturePlot(immu_sobj, reduction = "umap.unintegrated", features = "age_at_resection", pt.size = 0.5, cols = c("#00CCFF", "#CC0000"), raster = F)
ggsave(filename = paste0("04_umap_unintegrated_clusters-vs-age_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(non_int_dir, "Clustering"))
p = VlnPlot(immu_sobj, features = c("nFeature_RNA", "percent.mt", "nCount_RNA"), group.by = "unintegrated_clusters", sort = T, stack = T, flip = T, raster = F) + theme(legend.position = "none")
ggsave(filename = paste0("01_violin-stacked_unintegrated_qc-metrics_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 5, path = file.path(non_int_dir, "QC"))
p = dimplot_standard_qc_metrics(immu_sobj, "umap.unintegrated", raster = F)
ggsave(filename = paste0("02_umap_unintegrated_qc-metrics_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 15, path = file.path(non_int_dir, "QC"))


# --- 3. re-do batch integration and integration-based clustering and UMAP ---
table(immu_sobj$batch)
immu_sobj$RNA <- split(immu_sobj$RNA, f = immu_sobj$batch) # split into technical processing batches
anchors = 5
reduction_name = "bat.integrated.rpca" # name of new batch integration
options(future.globals.maxSize = 28000 * 1024^2)
immu_sobj = IntegrateLayers(immu_sobj, method = RPCAIntegration, k.weight = 90, k.anchor = anchors, orig.reduction = "pca", new.reduction = reduction_name, verbose = T)

pcs_chosen = 30; cluster_resolution = 0.4
immu_sobj <- FindNeighbors(immu_sobj, reduction = "bat.integrated.rpca", dims = 1:pcs_chosen)
immu_sobj <- FindClusters(immu_sobj, resolution = cluster_resolution, cluster.name = "bat.rpca_clusters")
immu_sobj = RunUMAP(immu_sobj, reduction = "bat.integrated.rpca", dims = 1:pcs_chosen, reduction.name = "umap.bat.rpca")
immu_sobj$seurat_clusters = immu_sobj$bat.rpca_clusters


# --- 4. plot new batch integration-based clustering and UMAP ---

# plot clusters vs samples
p = DimPlot(immu_sobj, reduction = "umap.bat.rpca", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.bat.rpca", group.by = "orig.ident", raster = F)
ggsave(filename = paste0("01a_umap_integrated_clusters-vs-samples_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(bat_int_dir, "Clustering"))
# plot clusters vs patients
p = DimPlot(immu_sobj, reduction = "umap.bat.rpca", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.bat.rpca", group.by = "patients", raster = F)
ggsave(filename = paste0("01b_umap_integrated_clusters-vs-patients_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(bat_int_dir, "Clustering"))
# plot clusters vs batches
p = DimPlot(immu_sobj, reduction = "umap.bat.rpca", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.bat.rpca", group.by = "batch", raster = F)
ggsave(filename = paste0("02_umap_integrated_clusters-vs-batches_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(bat_int_dir, "Clustering"))
# plot clusters vs subtypes
p = DimPlot(immu_sobj, reduction = "umap.bat.rpca", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.bat.rpca", group.by = "subtype", raster = F)
ggsave(filename = paste0("03a_umap_integrated_clusters-vs-subtypes_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(bat_int_dir, "Clustering"))
# plot clusters vs subclasses
p = DimPlot(immu_sobj, reduction = "umap.bat.rpca", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.bat.rpca", group.by = "subclass", raster = F)
ggsave(filename = paste0("03b_umap_integrated_clusters-vs-subclass_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(bat_int_dir, "Clustering"))
# plot clusters vs subclasses (bc)
p = DimPlot(immu_sobj, reduction = "umap.bat.rpca", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(immu_sobj, reduction = "umap.bat.rpca", group.by = "subclass_bc", raster = F)
ggsave(filename = paste0("03c_umap_integrated_clusters-vs-subclass-bc_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = file.path(bat_int_dir, "Clustering"))
p = VlnPlot(immu_sobj, features = c("nFeature_RNA", "percent.mt", "nCount_RNA"), group.by = "bat.rpca_clusters", sort = T, stack = T, flip = T, raster = F) + theme(legend.position = "none")
ggsave(filename = paste0("01_violin-stacked_integrated_qc-metrics_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 5, path = file.path(bat_int_dir, "QC"))
p = dimplot_standard_qc_metrics(immu_sobj, "umap.bat.rpca", raster = F)
ggsave(filename = paste0("02_umap_integrated_qc-metrics_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 15, path = file.path(bat_int_dir, "QC"))

# save immune subset object
immu_sobj = JoinLayers(immu_sobj)
saveRDS(immu_sobj, paste0("ST_EPN_immune.rds"))





