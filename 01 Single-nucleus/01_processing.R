library(Seurat)
library(patchwork)
library(RColorBrewer)
library(ggplot2)

cloud_path = ""
source("utils/seurat_utils.R")

preprocessing_dir = "Preprocessing/"; if(!dir.exists(preprocessing_dir)) dir.create(preprocessing_dir)
plot_dir = paste0(preprocessing_dir, "QC/01 Before merging/"); if(!dir.exists(plot_dir)) dir.create(plot_dir, r=T)
raw_obj_dir = "Objects/Raw objects by sample/"; if(!dir.exists(raw_obj_dir)) dir.create(raw_obj_dir, r=T)

metadata = readxl::read_excel("masterlist_EPN-ST.xlsx", skip = 1); metadata = metadata[grepl("1", metadata$`Analyses*`), ]
sample_ids = metadata$`Study ID`
sample_kks = metadata$`KK code` # internal ID


# ---------------------- PART 0: OBJECT CREATION ---------------------
# 1. collect all raw .h5 files to prepare for input for CellBender:
for (i in 1:length(sample_ids)) {
  message("collecting raw .h5 file for sample '", sample_ids[i], "' [", format(Sys.time(), "%d.%m. %X"), "]")
  h5_file = paste0(cloud_path, metadata[metadata$`KK code` == sample_kks[i], ]$`SC data path`, sample_kks[i], "/count/sample_raw_feature_bc_matrix.h5")
  target_path = paste0(preprocessing_dir, "Raw h5 files/", sample_kks[i], "_sample_raw_feature_bc_matrix.h5")
  file.copy(h5_file, target_path)
}

# 2. >>> run CellBender with default parameters in Python on CUDA-compatible GPU linux server <<<

# 3. create objects based on consensus CellBender and CellRanger-filtered .h5 count matrices
objects = lapply(1:length(sample_ids), function(i) {
  message("creating seurat object for sample '", sample_ids[i], "' [", format(Sys.time(), "%d.%m. %X"), "]")
  filtered_h5_file = paste0(cloud_path, "CB_filtered_h5/", sample_kks[i], "/processed_filtered.h5")
  mtx = Read_CellBender_h5_Mat(filtered_h5_file) # default implementation (Read10X_h5) is bugged
  obj_cb = CreateSeuratObject(counts = mtx, min.cells = 3, min.features = 200, project = sample_ids[i])
  mtx = Read10X_h5(paste0(cloud_path, metadata[metadata$`KK code` == sample_kks[i], ]$`SC data path`, sample_kks[i], "/count/sample_filtered_feature_bc_matrix.h5"))
  obj_cr = CreateSeuratObject(counts = mtx, min.cells = 3, min.features = 200, project = sample_ids[i])
  obj = obj_cb[rownames(obj_cr), colnames(obj_cr)] # constructs the consensus (intersection) matrix of CellBender and CellRanger. Combines the CB-fitered matrix with CR-filtered cells and features
  return(obj)
})


# ---------------------- PART 1: PLOT QCs PRE-FILTERING -------------------- 

# 1. check standard QC metrics: features, counts, mt.pct
for (feature in c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
  plots = lapply(1:length(sample_ids), function(i) {
    if (feature == "percent.mt") objects[[i]][["percent.mt"]] <<- PercentageFeatureSet(objects[[i]], pattern = "^MT-")
    return(VlnPlot(objects[[i]], features = feature) + NoLegend() + ggtitle(""))
  })
  ggsave(filename = paste0("01_unfiltered_qc_", feature, ".png"), plot = wrap_plots(plots, ncol = 20), width = 40, height = 20, path = plot_dir)
}
plots = lapply(1:length(sample_ids), function(i) VlnPlot(objects[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
ggsave(filename = paste0("01_unfiltered_qc_combined.png"), plot = wrap_plots(plots), width = 40, height = 20, path = plot_dir)

# 2. check correlation of counts vs. mt.pct and counts vs. features
plots = lapply(1:length(sample_ids), function(i){ FeatureScatter(objects[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend() + 
    FeatureScatter(objects[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") })
ggsave(filename = paste0("01_unfiltered_qc_correlations_combined.png"), plot = wrap_plots(plots), width = 40, height = 18, path = plot_dir)


# ---------------------- PART 2: FILTERING -------------------- 

# based on manual inspection of all relevant quality metrics, filter out low quality cells (doublets, empty droplets, ...)
# --------------
min_features = 200
max_features = 8000
max_mt_percentage = 15
# --------------
objects = lapply(objects, subset, subset = nFeature_RNA >= min_features & nFeature_RNA <= max_features & percent.mt <= max_mt_percentage)


# ---------------------- PART 3: PLOT QCs POST-FILTERING -------------------- 
write.table(data.frame(min_features, max_features, max_mt_percentage), file = paste0(plot_dir, "filters_used.txt"), sep = "\t", row.names = FALSE)

# 1. check standard QC metrics: features, counts, mt.pct
for (feature in c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
  plots = lapply(1:length(sample_ids), function(i) VlnPlot(objects[[i]], features = feature) + NoLegend() + ggtitle(""))
  ggsave(filename = paste0("02_filtered_qc_", feature, ".png"), plot = wrap_plots(plots, ncol = 20), width = 40, height = 20, path = plot_dir)
}
plots = lapply(1:length(sample_ids), function(i) VlnPlot(objects[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
ggsave(filename = paste0("02_filtered_qc_combined.png"), plot = wrap_plots(plots), width = 40, height = 20, path = plot_dir)

# 2. check correlation of counts vs. mt.pct and counts vs. features
plots = lapply(1:length(sample_ids), function(i){ FeatureScatter(objects[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend() + 
    FeatureScatter(objects[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") })
ggsave(filename = paste0("02_filtered_qc_correlations_combined.png"), plot = wrap_plots(plots), width = 40, height = 18, path = plot_dir)


# ---------------------- PART 4: AUTOMATED DOUBLET/MULTIPLET REMOVAL -------------------- 
library(DoubletFinder)
plot_dir = paste0(plot_dir, "Doublet removal/"); if (!dir.exists(plot_dir)) dir.create(plot_dir)
# Set the estimated the number of doublets (e.g., 4% of cells)
exp_doublet_pct = 0.04 # expected percentage of doublets
objects_doub_filt = c()

for (i in 1:length(sample_ids)) { # takes < 1 hour
  message("removing doublets in sample '", sample_ids[i], "' [", format(Sys.time(), "%d.%m. %X"), "]")
  obj = objects[[i]]
  obj = obj %>% NormalizeData(verbose = F) %>% FindVariableFeatures(verbose = F) %>% ScaleData(verbose = F) %>% 
    RunPCA(features = VariableFeatures(obj), verbose = F) %>% RunUMAP(dims = 1:10, verbose = F)
  
  # 1. run a parameter sweep to find the optimal pK value for DoubletFinder
  {sink("null"); sweep_res <- paramSweep(obj, PCs = 1:10); sink();} # parameter sweep (sink just for message suppresion)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)  # summarize results
  pK <- find.pK(sweep_stats)  # find the optimal pK
  optimal_pK <- as.character(pK$pK[which.max(pK$BCmetric)])  # select pK with highest BC metric
  message("Optimal pK value found: ", optimal_pK)
  
  # 2. run DoubletFinder to classify doublets
  n_doublets <- round(ncol(obj) * exp_doublet_pct)
  message("Estimated number of doublets in sample: ", n_doublets)
  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = as.numeric(optimal_pK), nExp = n_doublets, sct = F)
  
  # 3. visualize UMAP with doublet status
  doublet_status <- obj@meta.data[[paste0("DF.classifications_0.25_", optimal_pK, "_", n_doublets)]]  # add doublet status
  obj$doublet_status = doublet_status
  p = DimPlot(obj, reduction = "umap", group.by = "doublet_status") + ggtitle("DoubletFinder Results")
  ggsave(filename = paste0("umap_doublet_finder_result_", sample_ids[i], ".png"), plot = p, width = 8, height = 4, path = plot_dir)
  
  # 4. filter out doublets (in unprocessed initial object)
  obj = objects[[i]]; obj$doublet_status = doublet_status
  obj_filtered <- subset(obj, subset = doublet_status == "Singlet")  # keep only singlets
  message("Removed ", ncol(obj) - ncol(obj_filtered), " cells (before: ", ncol(obj), ", after: ", ncol(obj_filtered), ")\n")
  objects_doub_filt[[i]] = obj_filtered
  gc()
}
objects = objects_doub_filt # overwrites previous objects

# save filtered objects
for (i in 1:length(sample_ids)) saveRDS(objects[[i]], file = paste0(raw_obj_dir, "filtered_", sample_ids[i], ".rds"))


# ---------------------- PART 5: MERGE FILTERED OBJECTS --------------------            
raw_obj_dir = "Objects/Raw objects by sample/"
unintegrated_clust_dir = paste0(preprocessing_dir, "/Clustering/01 Before integration/"); if(!dir.exists(unintegrated_clust_dir)) dir.create(unintegrated_clust_dir, r=T)
unintegrated_qc_dir = paste0(preprocessing_dir, "QC/02 Before integration/"); if(!dir.exists(unintegrated_qc_dir)) dir.create(unintegrated_qc_dir, r=T)

# 1. check cohort quality (sequencing QC metrics)
# plot number of cells specifically
data = metadata; data$Samples = data$`Study ID`
p = ggplot(data, aes(x=Samples, y=Cells)) + geom_bar(stat="identity") + theme_classic() + theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle("Number of cells per sample")
ggsave(filename = paste0("03_barplot_n-cells.png"), plot = p, width = 15, height = 5, path = "Preprocessing/QC/01 Before merging/")

# 2. exclusion & reduction of samples
# -> decide on which samples to reduce if necessary <-
sample_ids = metadata$`Study ID` # update to new/filtered sample_ids
cell_limit = 6000 # limit will be applied to all samples

# 3. merging
sobjects = list()
for (i in 1:length(sample_ids)) { # load raw, but pre-filtered objects
  sobjects[[sample_ids[i]]] = readRDS(file = paste0(raw_obj_dir, "filtered_", sample_ids[i], ".rds")) 
}
sobjects = reduce_cell_numbers(sobjects, cell_limit)
sobj <- merge(x = sobjects[[1]], y = sobjects[2:length(sample_ids)], project = "ST-EPN_merged") # merge sample objects into a single Seurat object containing individual layers for each sample

# 4. processing
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

# 5. plot biological and technical variables
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


# ---------------------- PART 5: INTEGRATE OBJECTS -------------------- 
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

