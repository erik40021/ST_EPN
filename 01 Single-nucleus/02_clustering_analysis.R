library(Seurat)
library(ggpubr)

source("utils/seurat_utils.R")
source("utils/r_utils.R")

sobj = readRDS(file = "Objects/ST_EPN.rds")
out_dir = "Main analysis/Clustering"; if (!dir.exists(out_dir)) dir.create(out_dir, r = T)



# ---------------- PART 1: cluster annotation, malignancy classification, and labeling of integrated object -------------------

p = DimPlot(sobj, reduction = "umap.rpca", label = TRUE, label.size = 7, raster = F) + NoLegend()
ggsave(filename = paste0("00_umap_rpca_numbered.png"), plot = p, width = 14, height = 10, path = out_dir)

# i) calculate DEGs and annotate with GSEA and differential GO terms
sobj@active.ident = sobj$rpca_clusters
sobj$rpca_clusters = factor(sobj$rpca_clusters, levels = 0:(length(unique(sobj$rpca_clusters))-1))
degs = FindAllMarkers(sobj, logfc.threshold = 0.25)
save_degs_sheetwise(degs, paste0("Main analysis/Annotation/DEGs_ST-EPN_wilcox_", format(Sys.Date(), "%d-%m-%Y"), ".xlsx"))

annotate_clusters_GSEA(list.files("Main analysis/Annotation", pattern = "DEGs_*", full.names = T), "Main analysis/Annotation/GSEA/")
degs = read_excel_allsheets(list.files("Main analysis/Annotation", pattern = "DEGs_*", full.names = T)) # read in to get into proper format
degs = degs[order(as.numeric(names(degs)))]
top_degs_list = lapply(degs, function(x) head(x, 100)$gene)
enrich_signatures_GO_differential(top_degs_list, "Main analysis/Annotation/GO", save_combined = T, save_individually = T)

# ii) plot expression of top DEGs as violins
p = plot_violin_signature(sobj, top_degs_list)
ggsave(filename = paste0("01a_vln_all_top100-degs.png"), plot = p, width = 12, height = 12, path = "Main analysis/Annotation")
p = plot_violin_signature(sobj, sapply(degs, function(x) head(x, 2)$gene), sort = F)
ggsave(filename = paste0("01b_vln_all_top2-degs.png"), plot = p, width = 12, height = 16, path = "Main analysis/Annotation")

# iii) label non-malignant/TME clusters based on the previous annotations
cluster_labels = c("(0)", "Myeloid (1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)", "(10)", "(11)", "(12)", 
                   "Endothelial (13)", "(14)", "Stromal (15)", "(16)", "(17)", "(18)","(19)", "Lymphoid.2 (20)", "Oligodendrocyte (21)", 
                   " OPC (22)", "(23)", "Neuron (24)", "Lymphoid (25)", "(26)", "Astrocyte (27)", "Myeloid.SE (28)")
sobj$rpca_labels = factor(sobj$rpca_clusters, labels = cluster_labels)
p = DimPlot(sobj, reduction = "umap.rpca", group.by = "rpca_labels", label = TRUE, label.size = 5, raster = F) + NoLegend()
ggsave(filename = paste0("05_umap_rpca_TME-annotated.png"), plot = p, width = 14, height = 10, path = out_dir)

# iv) select reference cells for CNA inference
ref_candidates = c(1, 20, 25, 28)
p = plot_cluster_composition(sobj, "subtype", ref_candidates) # check how many cells each sample contributes to these clusters (important for inferCNV reference cells)
ggsave(filename = paste0("06a_barplot_TME-cluster_compositions_by_subtype.png"), plot = p, width = 5, height = 6, path = out_dir)
p = plot_cluster_composition(sobj, "orig.ident", ref_candidates)
ggsave(filename = paste0("06b_barplot_TME-cluster_compositions_by_sample.png"), plot = p, width = 10, height = 8, path = out_dir)

# v) infer CNAs per cell, derive malignancy score, and cluster CNAs
# >>> follow '03_CNA_inference.R' <<<

# vi) transfer CNA inference results to object
non_malignant_clusters = c(1, 13, 15, 17, 20, 21, 22, 24, 25, 27, 28) # result of CNA inference and cluster annotation based malignancy classification 
malignant_clusters = as.numeric(setdiff(unique(sobj@meta.data$rpca_clusters), non_malignant_clusters))
sobj = store_malignancy_annotation(sobj, malignant_clusters, "rpca_clusters", "malignant") # update malignancy label to final decision

# vii) add new labels and plot 'malignant' and TME labels
sobj$rpca_clusters = factor(sobj$rpca_clusters, levels = 0:(length(unique(sobj$rpca_clusters))-1))
malig_tme_labels = rep("malignant", length(levels(sobj$rpca_clusters)))
malig_tme_labels[non_malignant_clusters+1] = c("Myeloid", "Endothelial", "Stromal", "Unassigned", "Lymphoid", "Oligodendrocyte", "OPC", "Neuron", "Lymphoid", "Astrocyte", "Myeloid")
sobj$rpca_labels_malignant_TME = factor(sobj$rpca_clusters, labels = malig_tme_labels)
sobj$rpca_labels_malignant_TME[sobj$malignant == "malignant"] = "malignant" # make sure also TME-cluster malig cells are labeled correctly
p = DimPlot(sobj, reduction = "umap.rpca", group.by = "rpca_labels_malignant_TME", cols = c("grey", gg_color_hue(length(non_malignant_clusters))), label.size = 7, label = F, raster = F)
ggsave(filename = "08_umap_malignancy_annotation_with-TME.png", plot = p, width = 15, height = 10, path = out_dir)

custom_cols = c(malignant = "grey80", Myeloid="#ff918a", Lymphoid = "#F2B701", Endothelial = "#EF4868", Stromal = "#bab475", Unassigned = "#d5bdaf", 
             Oligodendrocyte = "#4e81a3", OPC = "#4b4b8f", Neuron = "#e3d536", Astrocyte = "#11A579")
p = DimPlot(sobj, reduction = "umap.rpca", group.by = "rpca_labels_malignant_TME", cols = custom_cols, label.size = 7, label = F, raster = F)
ggsave(filename = "08_umap_malignancy_annotation_with-TME_custom.png", plot = p, width = 10, height = 8, path = out_dir)

sobj$rpca_labels_malignant_TME = factor(sobj$rpca_labels_malignant_TME, levels = c("malignant", "Unassigned", "Oligodendrocyte", "OPC", "Astrocyte", "Neuron", "Endothelial", "Stromal", "Myeloid", "Myeloid.2", "Lymphoid", "Lymphoid.2")) # change order
top_TME_genes = sapply(degs[c(22, 23, 28, 25, 14, 16, 2, 29, 21, 26)], function(x) head(x, 2)$gene)
p = plot_violin_signature(sobj, top_TME_genes, group.by = "rpca_labels_malignant_TME", sort = F, idents = c(21,22,27,24,13,15,1,28,20,25))
ggsave(filename = paste0("02_vln_TME_top2-degs.png"), plot = p, width = 6, height = 5, path = "Main analysis/Annotation")

custom_cols = c("#4e81a3", "#4e81a3", "#4b4b8f","#4b4b8f", "#11A579","#11A579", "#e3d536","#e3d536" ,"#EF4868","#EF4868", "#bab475","#bab475", "#ff918a","#ff918a", "#F2B701","#F2B701")
top_TME_genes = sapply(degs[c(22, 23, 28, 25, 14, 16, 2, 21)], function(x) head(x, 2)$gene)
p = plot_violin_signature(sobj, top_TME_genes, cols = custom_cols, group.by = "rpca_labels_malignant_TME", sort = F, idents = c(21,22,27,24,13,15,1,20))
ggsave(filename = paste0("02_vln_TME_top2-degs_custom.png"), plot = p, width = 6, height = 5, path = "Main analysis/Annotation")

top_TME_genes = lapply(degs[c(22, 23, 28, 25, 14, 16, 2, 21)], function(x) head(x, 50)$gene)
p = plot_violin_signature(sobj, top_TME_genes, as_heatmap = T, group.by = "rpca_labels_malignant_TME", cluster_rows = F, cluster_cols = F, sort = F)
ggsave(filename = paste0("02_heatmap_TME_top100-degs_custom.png"), plot = p, width = 4, height = 4, path = "Main analysis/Annotation")

# top_TME_genes = lapply(degs[c(22, 23, 28, 25, 14, 16, 2, 21)], function(x) head(x, 50)$gene)
# plot_violin_signature(sobj, top_TME_genes, as_heatmap = T, group.by = "rpca_labels_malignant_TME", cluster_rows = F, cluster_cols = F, sort = F)
# ggsave(filename = paste0("02_heatmap_TME_top100-degs_custom.png"), plot = p, width = 4, height = 4, path = "Main analysis/Annotation")

# viii) save integrated object after adding all previous labels
saveRDS(sobj, "Objects/ST_EPN.rds")


# ---------------- PART 2: compare subtype differences -------------------

# ----------- compare the variance of gene expression among subtypes -------------
data = calculate_variance_across_groups(sobj, "subtype")
data_all = data[!data$group %in% c("NOS", "PLAGL1"), ]; data_all$group = factor(data_all$group, levels = c("ZFTA", "SE", "YAP1"))
data_all$log_vars = log(data_all$gene_vars) + 1
ggviolin(data_all, x = "group", y = "log_vars", add = "boxplot", fill = "group", palette = subtype_cols, add.params = list(fill = "white")) + 
  labs(x = "subtype", y = "log-scaled variance (all genes)") + stat_compare_means(method = "t.test", label = "p.format", comparisons = list(c("YAP1", "ZFTA"), c("ZFTA", "SE"))) # Global significance across groups
ggsave("01a_barplot_variance_among_subtypes_all-genes.png", width = 6, height = 4, path = "Main analysis/Subtype comparison")

data_hvg = data_all[data_all$gene %in% VariableFeatures(sobj), ]
ggviolin(data_hvg, x = "group", y = "log_vars", add = "boxplot", fill = "group", palette = subtype_cols, add.params = list(fill = "white")) + 
  labs(x = "subtype", y = "log-scaled variance (only hvgs)") + stat_compare_means(method = "t.test", label = "p.format", comparisons = list(c("YAP1", "ZFTA"), c("ZFTA", "SE"))) # Global significance across groups
  # ylim(0, quantile(data_hvg$gene_vars, 0.99))
ggsave("01b_barplot_variance_among_subtypes_highly-variable-genes.png", width = 4, height = 4, path = "Main analysis/Subtype comparison")

# compare averages:
sapply(c("ZFTA", "SE", "YAP1"), function(st) mean(data_hvg[data_hvg$group == st, ]$gene_vars))
sapply(c("ZFTA", "SE", "YAP1"), function(st) mean(data_hvg[data_hvg$group == st, ]$log_vars))



       
# ---------------- PART 3: subset malignant cells, and process and cluster -------------------
       
# subset malignant object according to prior classification
malig_sobj = subset(sobj, cells = rownames(sobj@meta.data[sobj$malignant == "malignant", ]))

# --- 1. re-do default pre-processing on new subset of cells ---
malig_sobj <- NormalizeData(malig_sobj)
malig_sobj <- FindVariableFeatures(malig_sobj)
malig_sobj <- ScaleData(malig_sobj)
malig_sobj <- RunPCA(malig_sobj, verbose = FALSE)
pcs_chosen = 30; cluster_resolution = 0.6
malig_sobj <- FindNeighbors(malig_sobj, reduction = "pca", dims = 1:pcs_chosen)
malig_sobj <- FindClusters(malig_sobj, resolution = cluster_resolution, cluster.name = "unintegrated_clusters")
malig_sobj = RunUMAP(malig_sobj, reduction = "pca", dims = 1:pcs_chosen, reduction.name = "umap.unintegrated")

# --- 2. re-do batch integration and integration-based clustering and UMAP ---
malig_sobj$RNA <- split(malig_sobj$RNA, f = malig_sobj$batch) # split into technical processing batches
anchors = 5
reduction_name = "bat.integrated.rpca" # name of new batch integration
options(future.globals.maxSize = 28000 * 1024^2)
malig_sobj = IntegrateLayers(malig_sobj, method = RPCAIntegration, k.anchor = anchors, orig.reduction = "pca", new.reduction = reduction_name, verbose = T)
pcs_chosen = 30; cluster_resolution = 0.6
malig_sobj <- FindNeighbors(malig_sobj, reduction = "bat.integrated.rpca", dims = 1:pcs_chosen)
malig_sobj <- FindClusters(malig_sobj, resolution = cluster_resolution, cluster.name = "bat.rpca_clusters")
malig_sobj = RunUMAP(malig_sobj, reduction = "bat.integrated.rpca", dims = 1:pcs_chosen, reduction.name = "umap.bat.rpca")
malig_sobj$seurat_clusters = malig_sobj$bat.rpca_clusters

malig_sobj = JoinLayers(malig_sobj)
saveRDS(malig_sobj, "ST_EPN_malig.rds")



# ---------------- PART 4: Transfer metaprograms to Seurat object and score cells ---------------
       
# requires finished NMF program identification
# # >>> for that follow '04_NMF_programs.R' <<<
  
# assign each cell to one MP
mps = openxlsx::read.xlsx("NMF_metaprograms_main.xlsx"), sheet = "MP genes (final)")
mp_list = lapply(mps, function(col) col)
malig_sobj = Seurat::AddModuleScore(malig_sobj, features = mp_list, name = "MP") # add expression scores of each MP to each cell
malig_sobj = assign_cells_by_max_score(malig_sobj, score_names = paste0("MP", 1:length(mp_list)), meta_data_entry_name = "metaprogram")
levels(malig_sobj$metaprogram) = names(mp_list) # rename 'metaprogram' metadata entry

p = DimPlot(malig_sobj, reduction = "umap.bat.rpca", label = T, raster = F) + NoLegend() + DimPlot(malig_sobj, reduction = "umap.bat.rpca", group.by = "metaprogram", cols = mp_cols, raster = F)
ggsave(filename = "01_umap_metaprograms.png", plot = p, width = 20, height = 10, path = seurat_out_dir)



