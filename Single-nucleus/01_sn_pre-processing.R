library(Seurat)
library(patchwork)
library(ggpubr)

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


# ---------------------- PART 5: SAVE FILTERED OBJECTS -------------------- 
for (i in 1:length(sample_ids)) saveRDS(objects[[i]], file = paste0(raw_obj_dir, "filtered_", sample_ids[i], ".rds"))


