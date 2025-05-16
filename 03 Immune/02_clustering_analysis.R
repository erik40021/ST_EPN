library(Seurat)
library(ggplot2)

source("utils/r_utils.R")
source("utils/seurat_utils.R")

immu_sobj = readRDS(file = "ST_EPN_immune.rds")
base_dir = ""


# ------------------------ PART 1: analyse integrated object cluster-wise -------------------------

pcs_chosen = 30; cluster_resolution = 0.4
immu_sobj <- FindClusters(immu_sobj, resolution = cluster_resolution, cluster.name = "bat.rpca_clusters")
immu_sobj$seurat_clusters = immu_sobj$bat.rpca_clusters

# --- annotate clusters ---
out_dir = ""
p = DimPlot(immu_sobj, reduction = "umap.bat.rpca", label = TRUE, label.size = 7, raster = F) + NoLegend()
ggsave(filename = paste0("00_umap_rpca_numbered.png"), plot = p, width = 6, height = 6, path = out_dir)

# i) calculate DEGs and annotate with GSEA and differential GO terms
immu_sobj@active.ident = immu_sobj$bat.rpca_clusters
immu_sobj$bat.rpca_clusters = factor(immu_sobj$bat.rpca_clusters, levels = 0:(length(unique(immu_sobj$bat.rpca_clusters))-1))
degs = FindAllMarkers(immu_sobj, logfc.threshold = 0.25)
save_degs_sheetwise(degs, paste0("DEGs_ST-EPN_wilcox_", format(Sys.Date(), "%d-%m-%Y"), ".xlsx"))

annotate_clusters_GSEA(list.files("", pattern = "DEGs_*", full.names = T), "")
degs = read_excel_allsheets(list.files("", pattern = "DEGs_*", full.names = T)) # read in to get into proper format
degs = degs[order(as.numeric(names(degs)))]
top_degs_list = lapply(degs, function(x) head(x, 100)$gene)
enrich_signatures_GO_differential(top_degs_list, "", save_combined = T, save_individually = T)

p = plot_violin_signature(immu_sobj, top_degs_list)
ggsave(filename = paste0("01a_vln_all_top100-degs.png"), plot = p, width = 12, height = 12, path = out_dir)
p = plot_violin_signature(immu_sobj, sapply(degs, function(x) head(x, 2)$gene), sort = F)
ggsave(filename = paste0("01b_vln_all_top2-degs.png"), plot = p, width = 12, height = 16, path = out_dir)

# >> exclude clusters 5 and 12 because of high tumour and/or non-immune TME gene expression (5 = tumour cells, 12 = non-immune TME (or multiplets)) <<
immu_sobj = subset(immu_sobj, idents = setdiff(unique(immu_sobj$bat.rpca_clusters), c(12,5)))

# iii) label clusters using the previous annotations
cluster_labels = c('0'="M-nos", "1"="M-nos", "2"="M-nos", "3"="M-nos", "4"="M-nos", "5"="M-Interferon", "6"="M-Neutrophil", "7"="L-NKcell", "8"="M-hypoxic", 
                   "9"="L-Tcell", "10"="M-protfold", "11"="Cell-cycle", "12"="L-plasma", "13"="L-Bcell", "14"="Ion-trans", "15"="M-mast")
immu_sobj$rpca_labels = factor(immu_sobj$bat.rpca_clusters, labels = cluster_labels)
cluster_order = c("M-nos", "M-Interferon", "M-Neutrophil", "M-hypoxic","M-protfold", "M-mast", "L-NKcell", "L-Tcell", "L-plasma", "L-Bcell", "Cell-cycle", "Ion-trans")
immu_sobj$rpca_labels = factor(immu_sobj$rpca_labels, levels = cluster_order)
p = DimPlot(immu_sobj, reduction = "umap.bat.rpca", group.by = "rpca_labels", label = TRUE, label.size = 5, raster = F)
ggsave(filename = paste0("02_umap_rpca_TME-annotated.png"), plot = p, width = 9, height = 7, path = "Immune analysis/Annotation")



# ------------------------ PART 2: transfer NMF metaprograms to object -------------------------
# requires finished identification of immune NMF programs (same as in '02 Spatial/02_NMF_programs.R')
mps = openxlsx::read.xlsx(paste0("immune_NMF_metaprograms_main.xlsx"), sheet = "MP genes (final)")
seurat_out_dir = ""

# 1. assign each cell to one MP
immu_sobj = Seurat::AddModuleScore(immu_sobj, features = as.list(mps), name = "MP") # add expression scores of each MP to each cell
immu_sobj = assign_cells_by_max_score(immu_sobj, score_names = paste0("MP", 1:length(as.list(mps))), meta_data_entry_name = "metaprogram")
levels(immu_sobj$metaprogram) = names(as.list(mps)) # rename 'metaprogram' metadata entry
saveRDS(immu_sobj, "ST_EPN_immune.rds")

p = DimPlot(immu_sobj, reduction = "umap.bat.rpca", label = T, raster = F) + NoLegend() + DimPlot(immu_sobj, reduction = "umap.bat.rpca", group.by = "metaprogram", cols = my_cols, raster = F)
ggsave(filename = "01a_umap_metaprograms.png", plot = p, width = 15, height = 7, path = seurat_out_dir)



