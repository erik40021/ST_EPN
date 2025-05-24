library(Seurat)
library(ggpubr)

source("utils/r_utils.R")
source("utils/seurat_utils.R")

out_dir = ""

# ---------------- PART 1: prepare single-nucleus data ----------------
# 1. create object with cell assignment that differentiates individual immune programs and hypoxic and non-hypoxic tumour cells
sobj = readRDS("ST_EPN.rds")
malig_sobj = readRDS("ST_EPN_malig.rds")
sobj$immu_tumour = "NA"
sobj@meta.data[colnames(immu_sobj), ]$immu_tumour = as.character(immu_sobj$metaprogram)
sobj@meta.data[colnames(malig_sobj), ]$immu_tumour = ifelse(malig_sobj$metaprogram == "Hypoxia", "Tumour-hyp", "Tumour")
sobj$immu_tumour = factor(sobj$immu_tumour, levels = unique(sobj$immu_tumour))
sobj@active.ident = sobj$immu_tumour
immu_tumour_sobj = subset(sobj, idents = "NA", invert = T) # remove all non-tumour, non-immune cells
table(immu_tumour_sobj$immu_tumour)
immu_tumour_sobj$immu_tumour = droplevels(immu_tumour_sobj$immu_tumour) # drop empty 'non-malignant' level
saveRDS(immu_tumour_sobj, "ST_EPN_immu_tumour.rds")

# 2. downsample non-hypoxic tumour cells to balance difference in size compared to hypoxic tumour cells
tumor_cell_ids = WhichCells(immu_tumour_sobj, idents = "Tumour")
downsampled_cells <- sample(tumor_cell_ids, 4932) # Randomly sample cells from all non-hypoxic tumour cells, matching the size of hyp-tumour
selected_cells <- c(downsampled_cells, setdiff(Cells(immu_tumour_sobj), tumor_cell_ids))
immu_tumour_sobj_ds = subset(immu_tumour_sobj, cells = selected_cells)
table(Idents(immu_tumour_sobj_ds))
writeMM(immu_tumour_sobj_ds@assays$RNA$data, file = "matrix.mtx") # save normalized counts
# save gene and cell names
write(x = rownames(immu_tumour_sobj_ds@assays$RNA$data), file = "features.tsv")
write(x = colnames(immu_tumour_sobj_ds@assays$RNA$data), file = "barcodes.tsv")

# 3. generate metadata
immu_tumour_sobj_ds@meta.data$cell = rownames(immu_tumour_sobj_ds@meta.data)
metadata = immu_tumour_sobj_ds@meta.data[, c("cell", "immu_tumour")]
write.table(metadata, file = file.path(save_dir, "meta.tsv"), sep="\t", quote=F, row.names=F)


# ---------------- PART 2: CCI analysis using CellPhoneDB v5 ----------------
# 1. run CellPhoneDB v5.0.0 with default parameters (function 'cpdb_statistical_analysis_method') in Python

# 2. (optional) visualise results with InterCellar
InterCellar::run_app(reproducible = TRUE) # starts shiny app in browser

# 3. customise interaction plots
interactions = read.csv(file.path(out_dir, "preprocessed_table.csv")) # using InterCellar output of CellPhoneDB results, but CellPhoneDB results can be used directly too

dotplot_interactions = function(interactions, senders = NULL, receivers = NULL, score_threshold = 0, differential = F, min_diffs = c(s=0.5,p=0.1), max_point_size = 5) {
  if (is.null(senders)) senders = unique(interactions$clustA)
  if (is.null(receivers)) receivers = unique(interactions$clustB)
  sel_ints <- interactions %>% filter(clustA %in% senders & clustB %in% receivers & score >= score_threshold)
  sel_ints$clust_pair = paste0(sel_ints$clustA, "::", sel_ints$clustB)
  sel_ints <- sel_ints %>% arrange(desc(score)); sel_ints$int_pair = factor(sel_ints$int_pair, levels = unique(sel_ints$int_pair))
  sel_ints <- sel_ints %>% arrange(clustA)
  sel_ints$clust_pair = factor(sel_ints$clust_pair, levels = paste0(rep(unique(sel_ints$clustA), each = length(unique(sel_ints$clustB))), "::", rev(receivers)))
  if (differential) {
    # currently only works for two groups probably
    diff_ints = apply(sel_ints, 1, function(int) { # check which interactions are differential
      int_to_compare = sel_ints[sel_ints$int_pair == int["int_pair"] & sel_ints$clustB == int["clustB"], ]
      if (nrow(int_to_compare) > 1) {
        diff(int_to_compare$score) > min_diffs[1] | diff(int_to_compare$p_value) > min_diffs[2]
      } else TRUE
    })
    sel_ints = sel_ints[diff_ints, ]
    diffs <- sel_ints %>% group_by(int_pair, clustA) %>% dplyr::summarize(avg_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = clustA, values_from = avg_score)
    diffs$is_class1 = as.logical(!is.na(diffs[, 2]) & (is.na(diffs[, 3]) | diffs[, 2] > diffs[, 3]))
    diffs$max_score = apply(diffs, 1, function(row) max(row[2:3], na.rm = T))
    diffs = diffs %>% arrange(is_class1, desc(max_score))
    sel_ints$int_pair = factor(sel_ints$int_pair, levels = diffs$int_pair)
  }
  p = ggplot(sel_ints, aes(x = int_pair, y = clust_pair)) + geom_point(aes(size = score, color = p_value)) +  # Dot size by score, color by -log10(p_value)
    scale_color_gradient2(low = "black", mid = "#003967", high = "grey95", midpoint = 0.01, name = "p value") +
    scale_size_continuous(range = c(1, max_point_size), name = "Score") +
    labs(title = paste0("Interactions (score > ", score_threshold, ")"), x = "Target", y = "Interaction Pair", size = "Score") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, size = 16))
  return(p)
}

# plot all interactions from tumour -> immune and tumour_hyp -> immune, with score > min_score
dotplot_interactions(interactions, sender = "Tumour", receivers = setdiff(unique(interactions$clustB), c("Tumour", "Tumour-hyp")), min_score = 1)

senders = c("Tumour", "Tumour-hyp")
receivers = c("Hypoxia-re", "Sys.inflam", "HS-UPR", "T.cell-naive", "Macrophage", "T.cell-cyt", "Monocyte", "Microglia", "B.cell", "Cell-cycle")
dotplot_interactions(interactions, senders = senders, receivers = receivers, score_threshold = 0)

dotplot_interactions(interactions, senders = senders, receivers = receivers, score_threshold = 0, differential = T, 
                     min_diffs = c(0.2, 0.1), max_point_size = 4.5) # + coord_flip()
ggsave("02_cci_dotplot_tumour-hyp_vs_tumour.png", width = 9, height = 5, path = out_dir)


# plot differential interactions from both tumour and tumour-hyp as circle plot
cols = c(Tumour="black", 'Hypoxia-re'="grey30", Sys.inflam="brown", 'HS-UPR'="#d5bdaf", "T.cell-naive"="#4e81a3", 
         Macrophage="#EF4868", "T.cell-cyt"="#66C5CC", Monocyte="#ff918a", Microglia="#967bad", B.cell="#e3d536", "Cell-cycle"="#F2B701")

png(file.path(out_dir, "03a_cci_circleplot_tumour-immune.png"), width = 2200, height = 2200, res = 300)
circlePlot(sel_ints[sel_ints$clustA == "Tumour", ], cols, score_cols = c("grey90", "#3969AC", "#00203b"))
dev.off()

png(file.path(out_dir, "03b_cci_circleplot_tumour-immune.png"), width = 2200, height = 2200, res = 300)
circlePlot(sel_ints[sel_ints$clustA == "Tumour-hyp", ], cols, score_cols = c("grey90", "#3969AC", "#00203b"))
dev.off()

