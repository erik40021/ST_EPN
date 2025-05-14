library(Seurat)
library(patchwork)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

source("utils/seurat_utils.R")

metadata = readxl::read_excel("masterlist_EPN-ST.xlsx", skip = 1); metadata = metadata[grepl("1", metadata$`Analyses*`), ]
sample_ids = metadata$`Study ID`

preprocessing_dir = "Preprocessing/"
raw_obj_dir = "Objects/Raw objects by sample/"
unintegrated_clust_dir = paste0(preprocessing_dir, "/Clustering/01 Before integration/"); if(!dir.exists(unintegrated_clust_dir)) dir.create(unintegrated_clust_dir, r=T)
unintegrated_qc_dir = paste0(preprocessing_dir, "QC/02 Before integration/"); if(!dir.exists(unintegrated_qc_dir)) dir.create(unintegrated_qc_dir, r=T)

### ----- check cohort quality (sequencing QC metrics) ----- ###
# plot number of cells specifically
data = metadata; data$Samples = data$`Study ID`
p = ggplot(data, aes(x=Samples, y=Cells)) + geom_bar(stat="identity") + theme_classic() + theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle("Number of cells per sample")
ggsave(filename = paste0("03_barplot_n-cells.png"), plot = p, width = 15, height = 5, path = "Preprocessing/QC/01 Before merging/")

### ---------- exclusion & reduction of samples -------- ###
# -> decide on which samples to reduce if necessary <-
sample_ids = metadata$`Study ID` # update to new/filtered sample_ids
cell_limit = 6000 # limit will be applied to all samples

### ------------------------- MERGING of raw objects --------------------------- ###
sobjects = list()
for (i in 1:length(sample_ids)) { # load raw, but pre-filtered objects
  sobjects[[sample_ids[i]]] = readRDS(file = paste0(raw_obj_dir, "filtered_", sample_ids[i], ".rds")) 
}
sobjects = reduce_cell_numbers(sobjects, cell_limit)
sobj <- merge(x = sobjects[[1]], y = sobjects[2:length(sample_ids)], project = "ST-EPN_merged") # merge sample objects into a single Seurat object containing individual layers for each sample
### ---------------------------------------------------------------------------- ###

# run streamline processing functions (Normalization, FindVariableFeatures, ScaleData and RunPCA) for each layer individually (layers are not yet integrated)
sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj, verbose = FALSE)
sobj <- RunPCA(sobj, verbose = FALSE)

# calculate clustering and UMAP reduction to visualise the merged, but still unintegrated data
pcs_chosen = 50
cluster_resolution = 0.6
sobj <- FindNeighbors(sobj, dims = 1:pcs_chosen, reduction = "pca")
sobj <- FindClusters(sobj, resolution = cluster_resolution, cluster.name = "unintegrated_clusters")
sobj <- RunUMAP(sobj, dims = 1:pcs_chosen, reduction = "pca", reduction.name = "umap.unintegrated")


# plot clusters vs samples
p = DimPlot(sobj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(sobj, reduction = "umap.unintegrated", group.by = "orig.ident", raster = F, label = T)
ggsave(filename = paste0("umap_unintegrated_clusters-vs-samples_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = unintegrated_clust_dir)

# plot clusters vs batches
sobj = add_annotation_from_excelsheet(sobj, metadata, "Batch", "batch")
p = DimPlot(sobj, reduction = "umap.unintegrated", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(sobj, reduction = "umap.unintegrated", group.by = "batch", raster = F)
ggsave(filename = paste0("umap_unintegrated_clusters-vs-batches_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = unintegrated_clust_dir)

# plot clusters vs subtypes
sobj = add_annotation_from_excelsheet(sobj, metadata, "Subtype", "subtype")
p = DimPlot(sobj, group.by = "unintegrated_clusters", reduction = "umap.unintegrated", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(sobj, reduction = "umap.unintegrated", group.by = "subtype", raster = F)
ggsave(filename = paste0("umap_unintegrated_clusters-vs-subtypes_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = unintegrated_clust_dir)

# plot clusters vs subclass
sobj = add_annotation_from_excelsheet(sobj, metadata, "Subclass", "subclass")
p = DimPlot(sobj, group.by = "unintegrated_clusters", reduction = "umap.unintegrated", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(sobj, reduction = "umap.unintegrated", group.by = "subclass", raster = F)
ggsave(filename = paste0("umap_unintegrated_clusters-vs-subclass_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = unintegrated_clust_dir)

# plot clusters vs subclass (classifier)
sobj = add_annotation_from_excelsheet(sobj, metadata, "Subclass (classifier)", "subclass_bc")
p = DimPlot(sobj, group.by = "unintegrated_clusters", reduction = "umap.unintegrated", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(sobj, reduction = "umap.unintegrated", group.by = "subclass_bc", raster = F)
ggsave(filename = paste0("umap_unintegrated_clusters-vs-subclass-bc_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = unintegrated_clust_dir)

# plot clusters vs age
sobj = add_annotation_from_excelsheet(sobj, metadata, "Age at resection", "age_at_resection")
p = DimPlot(sobj, group.by = "unintegrated_clusters", reduction = "umap.unintegrated", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  FeaturePlot(sobj, reduction = "umap.unintegrated", features = "age_at_resection", pt.size = 0.5, cols = c("#00CCFF", "#CC0000"), raster = F)
ggsave(filename = paste0("umap_unintegrated_clusters-vs-age_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = unintegrated_clust_dir)

# plot clusters vs recurrence status
sobj = add_annotation_from_excelsheet(sobj, metadata, "Rec", "recurrence")
p = DimPlot(sobj, group.by = "unintegrated_clusters", reduction = "umap.unintegrated", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(sobj, reduction = "umap.unintegrated", group.by = "recurrence", raster = F)
ggsave(filename = paste0("umap_unintegrated_clusters-vs-recurrence_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = unintegrated_clust_dir)

# plot clusters vs patients
metadata$patient_id = sub("(\\d+)[A-Za-z].*$", "\\1", metadata$`Study ID`)
sobj = add_annotation_from_excelsheet(sobj, metadata, "patient_id", "patients")
p = DimPlot(sobj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + 
  DimPlot(sobj, reduction = "umap.unintegrated", group.by = "patients", raster = F, label = T)
ggsave(filename = paste0("umap_unintegrated_clusters-vs-patients_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 20, height = 9, path = unintegrated_clust_dir)


# -- Quality control BEFORE integration --
# 1. investigate distribution of high expressing cells. Expected for proliferating and undiff. tumor cells. Otherwise an indication for technical sequencing doublets!
p = FeaturePlot(sobj, reduction = "umap.unintegrated", features = "nFeature_RNA", raster = F)
ggsave(filename = paste0("umap_unintegrated_feature-expression_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 7, path = unintegrated_qc_dir)

p = VlnPlot(sobj, features = "nFeature_RNA", raster = F, alpha = 0.2)
ggsave(filename = paste0("violin_unintegrated_feature-expression_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 7, path = unintegrated_qc_dir)

# 2. investigate mitochondrial gene expression: If strong mitochondrial expression-dependent clustering is visible: Check for the three signs of low quality cells and optionally perform more stringent filtering
p = FeaturePlot(sobj, reduction = "umap.unintegrated", features = "percent.mt", raster = F)
ggsave(filename = paste0("umap_unintegrated_mt-percent-expression_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 7, path = unintegrated_qc_dir)

p = VlnPlot(sobj, features = "percent.mt", raster = F, alpha = 0.2)
ggsave(filename = paste0("violin_unintegrated_mt-percent-expression_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 7, path = unintegrated_qc_dir)

# 2.1 investigate also the total expression of mitochondrial genes
sobj = add_absolute_feature_expression(sobj, pattern = "^MT-", col.name = "total.mt")
p = FeaturePlot(sobj, reduction = "umap.unintegrated", features = "total.mt", raster = F)
ggsave(filename = paste0("umap_unintegrated_mt-total-expression_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 7, path = unintegrated_qc_dir)

p = VlnPlot(sobj, features = "total.mt", raster = F, alpha = 0.2)
ggsave(filename = paste0("violin_unintegrated_mt-totalexpression_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 7, path = unintegrated_qc_dir)

# 3. investigate counts detected
p = FeaturePlot(sobj, reduction = "umap.unintegrated", features = "nCount_RNA", raster = F)
ggsave(filename = paste0("umap_unintegrated_counts-detected_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 7, path = unintegrated_qc_dir)

p = VlnPlot(sobj, features = "nCount_RNA", raster = F, alpha = 0.2)
ggsave(filename = paste0("violin_unintegrated_counts-detected_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 10, height = 7, path = unintegrated_qc_dir)

# 4. plot all single-cell qc metrics together as stacked violin plot
p = VlnPlot(sobj, features = c("nFeature_RNA", "percent.mt", "nCount_RNA"), sort = T, stack = T, flip = T, group.by = "orig.ident", raster = F, alpha = 0.2) + theme(legend.position = "none")
ggsave(filename = paste0("violin-stacked_qc-metrics_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 15, height = 5, path = unintegrated_qc_dir)

# ... and all together as umaps
p = dimplot_standard_qc_metrics(sobj, "umap.unintegrated", raster = F)
ggsave(filename = paste0("umaps_qc-metrics_res", cluster_resolution, "_pcs", pcs_chosen, ".png"), plot = p, width = 25, height = 13, path = unintegrated_qc_dir)

### ----------------------------- INTEGRATION ---------------------------------- ###
sobj$RNA <- split(sobj$RNA, f = sobj$batch) # split into those groups of samples that you want to integrate (e.g. technical batches, or samples)
anchors = 5
reduction_name = "integrated.rpca"
options(future.globals.maxSize = 28000 * 1024^2)
sobj = IntegrateLayers(sobj, method = RPCAIntegration, k.anchor = anchors, orig.reduction = "pca", new.reduction = reduction_name, verbose = T)
### ---------------------------------------------------------------------------- ###

umap_name = "umap.rpca"
pcs_chosen = 50
cluster_resolution = 0.5
sobj@misc = list(cluster_res_rpca = cluster_resolution, pcs_chosen_rpca = pcs_chosen) # store the chosen parameters for traceability
sobj <- FindNeighbors(sobj, reduction = reduction_name, dims = 1:pcs_chosen)
sobj <- FindClusters(sobj, resolution = cluster_resolution, cluster.name = "rpca_clusters")
sobj = RunUMAP(sobj, reduction = reduction_name, dims = 1:pcs_chosen, reduction.name = umap_name)


# ---- Quality control AFTER integration ----
int_qc_dir = paste0(preprocessing_dir, "QC/03 Integrated/"); if (!dir.exists(int_qc_dir)) dir.create(int_qc_dir)
# 1. investigate distribution of high expressing cells. Expected for proliferating and undiff. tumor cells. Otherwise an indication for technical sequencing doublets!
p = FeaturePlot(sobj, reduction = umap_name, features = "nFeature_RNA", raster = F)
ggsave(filename = paste0("umap_rpca-int_feature-expression_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 10, height = 5, path = int_qc_dir)
p = VlnPlot(sobj, features = "nFeature_RNA", raster = F, alpha = 0.2)
ggsave(filename = paste0("violin_rpca-int_feature-expression.png"), plot = p, width = 10, height = 5, path = int_qc_dir)

# 2. investigate counts
p = FeaturePlot(sobj, reduction = umap_name, features = "nCount_RNA", raster = F)
ggsave(filename = paste0("umap_rpca-int_count-expression_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 10, height = 5, path = int_qc_dir)

# 3. investigate mitochondrial gene expression: Usually, all clusters should have similar percent.mt, if not: perform more stringent filtering?
p = FeaturePlot(sobj, reduction = umap_name, features = "percent.mt", raster = F)
ggsave(filename = paste0("umap_rpca-int_mt-expression_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 10, height = 5, path = int_qc_dir)
p = VlnPlot(sobj, features = "percent.mt", raster = F, alpha = 0.2)
ggsave(filename = paste0("violin_rpca-int_mt-expression.png"), plot = p, width = 10, height = 5, path = int_qc_dir)


# ---- Investigate integrated clustering in more detail with biological and technical variables ----
int_dir = paste0(preprocessing_dir, "Clustering/02 rpca-integrated by batch/"); if(!dir.exists(int_dir)) dir.create(int_dir, r=T)

# assess clusters vs samples
p = DimPlot(sobj, reduction = umap_name, label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + DimPlot(sobj, reduction = umap_name, group.by = "orig.ident", raster = F)
ggsave(filename = paste0("01a_umap_rpca-int_clusters-vs-samples_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 20, height = 9, path = int_dir)

# assess clusters vs patients
p = DimPlot(sobj, reduction = umap_name, label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + DimPlot(sobj, reduction = umap_name, group.by = "patients", raster = F)
ggsave(filename = paste0("01b_umap_rpca-int_clusters-vs-patients_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 20, height = 9, path = int_dir)

# assess clusters vs batches
p = DimPlot(sobj, reduction = umap_name, label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + DimPlot(sobj, reduction = umap_name, group.by = "batch", raster = F)
ggsave(filename = paste0("02_umap_rpca-int_clusters-vs-batches_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 20, height = 9, path = int_dir)

# assess clusters vs subtypes
p = DimPlot(sobj, reduction = umap_name, label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + DimPlot(sobj, reduction = umap_name, group.by = "subtype", raster = F)
ggsave(filename = paste0("03a_umap_rpca-int_clusters-vs-subtypes_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 20, height = 9, path = int_dir)

# assess clusters vs subclasses
p = DimPlot(sobj, reduction = umap_name, label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + DimPlot(sobj, reduction = umap_name, group.by = "subclass", raster = F)
ggsave(filename = paste0("03b_umap_rpca-int_clusters-vs-subclass_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 20, height = 9, path = int_dir)

# assess clusters vs subclasses (classifier)
p = DimPlot(sobj, reduction = umap_name, label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + DimPlot(sobj, reduction = umap_name, group.by = "subclass_bc", raster = F)
ggsave(filename = paste0("03c_umap_rpca-int_clusters-vs-subclass-bc_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 20, height = 9, path = int_dir)

# assess clusters vs age
p = DimPlot(sobj, reduction = umap_name, label = T, raster = F) + NoLegend() + ggtitle("clusters.ident") + FeaturePlot(sobj, reduction = umap_name, features = "age_at_resection", raster = F, cols = c("lightblue", "#CC0000"))
ggsave(filename = paste0("04_umap_rpca-int_clusters-vs-age_res", cluster_resolution, "_k-anchors", anchors, ".png"), plot = p, width = 20, height = 9, path = int_dir)


# visualise more precisely each sample's distribution in UMAP
out_dir = file.path(int_dir, "split by sample"); if (!dir.exists(out_dir)) dir.create(out_dir)
lapply(1:(length(unique(sobj$orig.ident))/4), function(i) {
  samples = unique(sobj$orig.ident)[(4*(i-1)+1):(4*i)]
  subset_sobj = subset(sobj, subset = orig.ident %in% samples)
  p = DimPlot(subset_sobj, reduction = "umap.rpca", split.by = "orig.ident")
  ggsave(filename = paste0("umap_split-by-samples_", paste(samples, collapse = ","), ".png"), plot = p, width = 20, height = 9, path = out_dir)
})

sobj = JoinLayers(sobj) # <- run before saving if splitted before
saveRDS(sobj, "Objects/ST_EPN.rds")

