
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
