
# does GO enrichment with all three GO domains (BP, MF, CC) for gene signatures (as columns of a dataframe)
enrich_signatures_GO = function(sig_df, out_dir, do_plots=T) {
  library(openxlsx); library(org.Hs.eg.db); library(clusterProfiler); library(patchwork)
  plot_dir = file.path(out_dir, "GO dotplots"); if (!dir.exists(plot_dir)) dir.create(plot_dir, r = T)
  workbook = createWorkbook()
  for (i in 1:ncol(sig_df)) {
    sig = names(sig_df)[i]
    message("running GO enrichment for signature '", sig, "' (", ncol(sig_df)-(i-1), " left)")
    go_results <- enrichGO(gene= sig_df[[i]], OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL")
    if (do_plots) {
      go_bp = go_results; go_bp@result = go_bp@result[go_bp@result[[1]] == "BP",]
      if (length(go_bp@result[[1]]) > 0) p1 = clusterProfiler::dotplot(go_bp, x = "GeneRatio", color = "p.adjust", showCategory = 15, font.size = 10, title = "GO:BP")
      go_mf = go_results; go_mf@result = go_mf@result[go_mf@result[[1]] == "MF",]
      if (length(go_mf@result[[1]]) > 0) p2 = clusterProfiler::dotplot(go_mf, x = "GeneRatio", color = "p.adjust", showCategory = 15, font.size = 10, title = "GO:MF")
      go_cc = go_results; go_cc@result = go_cc@result[go_cc@result[[1]] == "CC",]
      if (length(go_cc@result[[1]]) > 0) p3 = clusterProfiler::dotplot(go_cc, x = "GeneRatio", color = "p.adjust", showCategory = 15, font.size = 10, title = "GO:CC")
      if (!any(!sapply(list(p1,p2,p3), is.null))) next
      ggsave(paste0(i, "_dotplots_go-enrichment_", sig, ".png"), plot = wrap_plots(Filter(Negate(is.null), list(p1,p2,p3)), nrow = 1), path = plot_dir, width = 20, height = 8)
    }
    addWorksheet(workbook, sig); writeData(workbook, sig, go_results@result, startCol = 1, colNames = TRUE)
    gc()
    p1 = NULL; p2 = NULL; p3 = NULL
  }
  saveWorkbook(workbook, file = paste0(out_dir, "/GO-enrichment_", format(Sys.Date(), "%d-%m-%Y"), ".xlsx"), overwrite = T)
}

# same as above but using differential GO enrichment (shows GO terms which are different across the signatures)
enrich_signatures_GO_differential = function(sig_df, out_dir = NULL, save_combined = T, save_individually = F, GO_label_width = 50, width = 15, height = 13, x_lab_size = 10) {
  library(openxlsx); library(org.Hs.eg.db); library(clusterProfiler); library(patchwork); library(stringr)
  message("comparing signatures regarding differential GO-BP terms ..."); comp_bp = compareCluster(sig_df, fun = "enrichGO", keyType = "SYMBOL", ont = "BP", OrgDb = "org.Hs.eg.db")
  message("comparing signatures regarding differential GO-MF terms ..."); comp_mf = compareCluster(sig_df, fun = "enrichGO", keyType = "SYMBOL", ont = "MF", OrgDb = "org.Hs.eg.db")
  message("comparing signatures regarding differential GO-CC terms ..."); comp_cc = compareCluster(sig_df, fun = "enrichGO", keyType = "SYMBOL", ont = "CC", OrgDb = "org.Hs.eg.db")
  p1 = clusterProfiler::dotplot(comp_bp, title = "Differential enrichment of GO-BP terms", showCategory = 4) + theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1, size = x_lab_size))
  p2 = clusterProfiler::dotplot(comp_mf, title = "Differential enrichment of GO-MF terms", showCategory = 4) + theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1, size = x_lab_size))
  p3 = clusterProfiler::dotplot(comp_cc, title = "Differential enrichment of GO-CC terms", showCategory = 4) + theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1, size = x_lab_size))
  if (is.null(out_dir)) { return(Filter(Negate(is.null), list(p1, p2, p3))) } else {
    if (!dir.exists(out_dir)) dir.create(out_dir)
    if (save_combined) ggsave(paste0("dotplot_differential_combined_GO-enrichment.png"), plot = wrap_plots(Filter(Negate(is.null), list(p1,p2,p3)), nrow = 1), path = out_dir, width = 32, height = 16)
    if (save_individually) { # save individual ontologies, allowing one more term per signature, and more space for GO term labels
      ggsave(paste0("dotplot_differential_GO-BP-enrichment.png"), plot = p1 + scale_y_discrete(labels = function(x) str_wrap(x, width = GO_label_width)), path = out_dir, width = width, height = height)
      ggsave(paste0("dotplot_differential_GO-MF-enrichment.png"), plot = p2 + scale_y_discrete(labels = function(x) str_wrap(x, width = GO_label_width)), path = out_dir, width = width, height = height)
      ggsave(paste0("dotplot_differential_GO-CC-enrichment.png"), plot = p3 + scale_y_discrete(labels = function(x) str_wrap(x, width = GO_label_width)), path = out_dir, width = width, height = height)
    }
  }
}


annotate_clusters_GSEA = function(deg_input_file, out_dir, make_plots = TRUE, type_readable = F) {
  library(clusterProfiler); library(msigdbr); library(org.Hs.eg.db); library(openxlsx); library(readxl); library(patchwork)
  if (!dir.exists(out_dir)) dir.create(out_dir, r = T)
  celltype_collection <- msigdbr(species = "Homo sapiens", category = "C8") # celltype signature
  hallmark_collection <- msigdbr(species = "Homo sapiens", category = "H") # hallmark signatures
  go_bp_collection <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") # GO biological process
  collections = list("celltype" = celltype_collection, "hallmark" = hallmark_collection, "GO-BP" = go_bp_collection)
  
  workbooks = list(); for (i in 1:length(collections)) workbooks = c(workbooks, createWorkbook())
  sheets = getSheetNames(deg_input_file)
  for (cluster in sheets) {
    df_deg <- data.frame(read_excel(path = deg_input_file, sheet = cluster))
    if (type_readable) df_deg = dplyr::bind_rows(df_deg[, 1:7], df_deg[, 9:15]) # bind both DEG parts (previously separated for readability) together
    df_deg = df_deg[!is.na(df_deg$gene) & df_deg$p_val_adj < 0.05, ]
    lfc <- df_deg$avg_log2FC; names(lfc) <- df_deg$gene; lfc_sorted_genes <- sort(lfc, decreasing = TRUE)
    message("--- collecting data of ", cluster, " (", length(sheets), " clusters total) for GSEA ---")
    res_plots = list()
    for (i in 1:length(collections)) {
      collection = data.frame(collections[[i]])
      message("starting GSEA for ", cluster, " using collection ", names(collections)[i])
      gsea_res <- clusterProfiler::GSEA(geneList = lfc_sorted_genes, minGSSize = 10, maxGSSize = 800, pvalueCutoff = 0.05, seed = F, 
        pAdjustMethod = "BH", TERM2GENE = dplyr::select(collection, gs_name, gene_symbol))
      gsea_df = data.frame(gsea_res, row.names = NULL); gsea_df$Description = NULL # extract data frame from result object and get in readable shape
      if (make_plots & length(gsea_df$ID) > 0) {
        p = clusterProfiler::dotplot(gsea_res, x = "NES", color = "p.adjust", showCategory = 15, font.size = 10, size = NULL)
        res_plots[[i]] = p
      }
      addWorksheet(workbooks[[i]], cluster)
      writeData(workbooks[[i]], cluster, gsea_df, startCol = 1, colNames = TRUE)
    }
    if (length(res_plots) > 0){ # save plots of the cluster for each collection together
      wrap_plots(Filter(Negate(is.null), res_plots), nrow = 1)
      ggsave(filename = paste0("gsea-dotplot_", cluster, "_combined.png"), width = 20, height = 10, path = out_dir)
    }
  }
  for (i in 1:length(collections)) saveWorkbook(workbooks[[i]], file = paste0(out_dir, "/gsea_", names(collections)[i], "_rpca-int_ST-EPN_integrated_", format(Sys.Date(), "%d-%m-%Y"), ".xlsx"), overwrite = TRUE) # save all workbooks
}


jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union*100)
}

compare_jaccard_ref_vs_sig = function(sig, reference_sig) {
  m = matrix(ncol = length(sig), nrow = length(reference_sig))
  colnames(m) = if (is.data.frame(sig)) { colnames(sig) } else { names(sig) }
  rownames(m) = if (is.data.frame(reference_sig)) { colnames(reference_sig) } else { names(reference_sig) }
  i = 1; j = 1;
  for (s in sig) {
    for (ref_s in reference_sig) {
      m[j, i] = jaccard(s, ref_s)
      j = j + 1
    }
    j = 1; i = i + 1
  }
  return(m)
}

reorder_excel_sheet_columns = function(file, sheets, order) {
  wb <- loadWorkbook(file)
  for (sheet in sheets) {
    sheet_data <- readWorkbook(wb, sheet = sheet)
    sheet_data <- sheet_data[, order]
    writeData(wb, sheet = sheet, x = sheet_data, startCol = 1, startRow = 1, colNames = TRUE)
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}

save_degs_sheetwise = function(degs, filename, order_signs_separately = F, order_by_log2fc = T) {
  library(openxlsx)
  wb = createWorkbook()
  for (i in unique(degs$cluster)){
    addWorksheet(wb, i)
    cluster_markers = degs[degs$cluster == i,]
    if (order_signs_separately) { cluster_markers$`x` = NA
      cluster_markers = cbind_fill_na(cluster_markers[cluster_markers$avg_log2FC > 0, ], cluster_markers[cluster_markers$avg_log2FC <= 0, ])
    } else if (order_by_log2fc) cluster_markers = cluster_markers[order(cluster_markers$avg_log2FC, decreasing = T), ]
    writeData(wb, i, cluster_markers, startCol = 1, rowNames = F, colNames = TRUE)
  }
  saveWorkbook(wb, filename, overwrite = TRUE)
}

read_excel_allsheets <- function(filename, tibble = FALSE, specific_col = NULL, skip = 0) {
  library(readxl)
  sheets <- readxl::excel_sheets(filename)
  if (is.null(specific_col)) {
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, skip = skip))
  } else {
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, skip = skip)[[specific_col]])
  }
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}

calculate_malignancy_score_by_CNAs = function(cna_mtx, print_stats = TRUE, normalise = T) {
  # 1. calculate "CNA signal" 
  cna_mtx = cna_mtx - 1 # inferCNV default output is centered around 1, but needs to be centered around 0 for applying abs()
  threshold <- apply(abs(cna_mtx), 1, quantile, probs = 0.75) # Defining threshold for high CNAs (e.g., above 75th percentile)
  high_cna_bool <- abs(cna_mtx) > threshold # boolean matrix marking high CNA regions per cell
  high_cna = abs(cna_mtx) * high_cna_bool # matrix with CNA values left only at the high CNA regions
  CNAsignal <- colMeans(high_cna)
  # 2. calculate "CNA correlation"
  average_profile <- rowMeans(cna_mtx, na.rm = TRUE)
  CNAcorrelation <- numeric(ncol(cna_mtx))
  for (i in 1:length(CNAcorrelation)) {
    cna_cor = cor(cna_mtx[, i], average_profile, use = "complete.obs", method = "pearson")
    CNAcorrelation[i] = ifelse(is.na(cna_cor), 0, cna_cor)
  }
  # 3. combine both scores (see Greenwald 2024, page 24 for reference)
  cna_score = CNAcorrelation + (CNAsignal * (max(CNAcorrelation) / max(CNAsignal)))
  cna_score_final = if (normalise) normalize01(cna_score) else cna_score
  if (print_stats) message("variance before normalisation: ", var(cna_score), " (sd = ", sd(cna_score), ")")
  return(cna_score_final)
}


heatmap_cna = function(cna_mtx, gene_positions, cluster_splits = NULL, row_anno = NULL, row_title = NULL, cluster_row_slices = F, quantiles = c(0.01, 0.99)) {
  library(ComplexHeatmap); library(circlize); library(RColorBrewer)
  quantiles <- quantile(cna_mtx, probs = quantiles)
  distances <- abs(cna_mtx - 1); distances[cna_mtx < quantiles[1] | cna_mtx > quantiles[2]] <- NA
  max_distance <- distances[which.max(distances)]
  color_palette = rev(brewer.pal(11, "RdBu")); color_mapping = colorRamp2(seq(1-max_distance, 1+max_distance, length.out = length(color_palette)), color_palette)
  gene_positions_subset <- gene_positions[match(colnames(cna_mtx), gene_positions$gene), ]
  chromosomes = sub("^chr", "", gene_positions_subset$chr)
  chr_splits = factor(chromosomes, levels = as.character(1:22))
  row_title = ifelse(is.null(row_title), ifelse(is.null(cluster_splits), "", "%s"), row_title)
  cluster_rows = cluster_row_slices
  chr_arm_anno = HeatmapAnnotation(arms = gene_positions_subset$arm, col = list(arms = c(p = "#ededed", q = "#474747")), which = "col", 
                                   show_annotation_name = F, show_legend = F, annotation_height = unit(0.5, "mm"))
  h = Heatmap(cna_mtx, name = "CNA expr.", show_column_names = FALSE, show_row_names = FALSE, cluster_rows = cluster_rows, cluster_columns = F, 
              left_annotation = row_anno, top_annotation = chr_arm_anno, cluster_row_slices = cluster_row_slices, row_title = row_title, column_split = chr_splits, 
              row_split = cluster_splits, col = color_mapping, border = T, column_gap = unit(0, "mm"), row_gap = unit(0, "mm"))
  return(h)
}

order_clusters_by_CNA_malignancy = function(malignancy_score, clusters) {
  cluster_list <- if (is.list(clusters)) {clusters} else {split(names(clusters), clusters)}
  mean_scores <- sapply(cluster_list, function(cells) mean(malignancy_score[cells]))
  ordered_cluster_names <- names(mean_scores)[order(mean_scores, decreasing = TRUE)]
  sorted_clusters <- cluster_list[ordered_cluster_names]
  sorted_clusters <- lapply(sorted_clusters, function(cells) {cells[order(malignancy_score[cells], decreasing = TRUE)]})
  reordered_clusters = unlist(lapply(names(sorted_clusters), FUN = function(cluster) rep(cluster, length(sorted_clusters[[cluster]]))))
  names(reordered_clusters) = unlist(sorted_clusters)
  return(reordered_clusters)
}
                                                                                               
# colors used in the study:
my_cols = c("#F2B701","#3969AC","#EF4868","#11A579","#6d2c6e","#ffe39e","#66C5CC","#ff918a","#748cab","#bd6908","#d5bdaf","#80BA5A","#CF1C90","yellow2",
             "#B95FBB","cyan2","#4b4b8f","#D4D915","grey","#ebc5ae","brown","#bab475","#4e81a3","#967bad","#542f3d","#f5bfd3","#adf590","#e3d536","#cb1adb")
subtype_cols = c(ZFTA = "#3969AC", SE = "#c73210", YAP1 = "#e3d536", PLAGL1 = "#b0977b", NOS = "#ebc5ae")
subclass_cols = c(ZFTA_RELA="#003967", ZFTA_NC="#66C5CC", ST_SE="#F97B72", SP_SE="#bd6908", YAP1="#e3d536", PLAGL1="#b0977b")
mp_cols = c(Stress="#F2B701", Hypoxia="#3969AC", Cilia="#ff918a", 'Cilia-reg'="#11A579", Astroglial="#6d2c6e", 'Cell-cycle'="#ffe39e", 
            'ZFTA-fus'="#66C5CC", NPC="#EF4868", Respiration="#748cab", 'ECM-remod'="#bd6908", Wounded="#d5bdaf", 'Interferon-re'="#80BA5A")
st_mp_cols = c(Stress="#F2B701", Hypoxia="#3969AC", Cilia="#ff918a", Glial="#11A579", Astroglial="#6d2c6e", Fibroblast="#bd6908", 
               'ZFTA-fus'="#66C5CC", 'Cell-cycle'="#ffe39e", Myeloid="#f5bfd3", Vasc="#d5bdaf", Inflammatory="#EF4868", 'Interferon-re'="#80BA5A")

