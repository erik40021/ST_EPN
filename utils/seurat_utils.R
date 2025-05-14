
Read_CellBender_h5_Mat <- function(file_name, use.names = TRUE, unique.features = TRUE) {
  # Check hdf5r installed
  if (!requireNamespace('hdf5r', quietly = TRUE)) cli_abort(message = c("Please install hdf5r to read HDF5 files", "i" = "`install.packages('hdf5r')`"))
  # Check file
  if (!file.exists(file_name)) stop("File not found")  
  if (use.names) {
    feature_slot <- 'features/name'
  } else {
    feature_slot <- 'features/id'
  }
  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")
  counts <- infile[["matrix/data"]]
  indices <- infile[["matrix/indices"]]
  indptr <- infile[["matrix/indptr"]]
  shp <- infile[["matrix/shape"]]
  features <- infile[[paste0("matrix/", feature_slot)]][]
  barcodes <- infile[["matrix/barcodes"]]  
  sparse.mat <- Matrix::sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(x = counts[]), dims = shp[], repr = "T")
  if (unique.features) features <- make.unique(names = features)  
  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as(object = sparse.mat, Class = "CsparseMatrix")
  infile$close_all()
  return(sparse.mat)
}


plot_features_signature = function(sobj, sig_name, features, reduction = "umap.unintegrated", do_plot = T, raster = F){
  gene.set <- rownames(GetAssay(sobj)$data[rownames(GetAssay(sobj)$data) %in% features,])
  diff = setdiff(features, gene.set)
  if (length(diff) > 0) {
    cat(" Warning! Non-matching gene symbols: \n", diff, "\n")
    features = intersect(features, gene.set)
  }
  # get mean expression of features of interest per cell
  mean.exp <- colMeans(GetAssay(sobj)$data[gene.set, ], na.rm = TRUE)
  # add mean expression values in 'sobj@meta.data$gene.set.score'
  if (all(names(x = mean.exp) == rownames(x = sobj@meta.data))) {
    sobj@meta.data$gene.set.score <- mean.exp
  }
  plot = FeaturePlot(sobj, features = "gene.set.score", cols = c("lightgrey", "red"), reduction = reduction, split.by = NULL, raster = raster) + 
    theme(plot.subtitle = element_text(hjust = 0.5))
  title = if (length(features) > 9) {
    if (length(features) <= 20) { ggtitle(sig_name, subtitle = paste0("[", toString(features[1:(length(features)/2)]), ",\n", toString(features[floor(length(features)/2+1):length(features)]), "]")) 
      } else { ggtitle(sig_name, subtitle = paste0("[", toString(features[1:10]), ",\n", toString(features[11:19]), ", ...]")) }
    } else { ggtitle(sig_name, subtitle = paste0("[", toString(features), "]")) }
  plot = plot + title
  if (do_plot) print(plot)
  return(plot)
}

plot_violin_signature = function(sobj, features, sig_name = NA, stacked = T, as_heatmap = F, x_lab = NULL, y_lab = NULL, cols = NULL,
                                 cluster_rows = T, cluster_cols = T, group.by = NULL, sort = T, scale_heatmap = "none", idents = NULL) {
  if (stacked) {
    if (is.list(features)) {
      if (is.null(y_lab)) { warning("using names of features list as y-axis labels"); y_lab = names(features) }
      for (i in 1:length(features)){
        gene.set <- rownames(GetAssay(sobj)$data[rownames(GetAssay(sobj)$data) %in% features[[i]],])
        diff = setdiff(features[[i]], gene.set)
        if (length(diff) > 0) {
          cat(" Warning! Non-matching gene symbols: \n", diff, "\n")
          features[[i]] = intersect(features[[i]], gene.set)
        }
        # get mean expression of features of interest per cell
        mean.exp <- colMeans(GetAssay(sobj)$data[gene.set, ], na.rm = TRUE)
        if (all(names(x = mean.exp) == rownames(x = sobj@meta.data))) {
          sobj@meta.data[[y_lab[i]]] <- mean.exp
        }
      }
      # if (is.null(x_lab)) x_lab = 0:length(features)-1
      plot = VlnPlot(sobj, stack = TRUE, flip = TRUE, features = y_lab, group.by = group.by, idents = idents, cols = cols) + NoLegend() + ggtitle(sig_name) + labs(x = group.by) # + scale_x_discrete(labels=x_lab)
      if (as_heatmap) {
        library(pheatmap); library(dplyr)
        data_raw = plot$data %>% group_by(feature, ident) %>% dplyr::summarize(mean_expression = mean(expression))
        data_raw$ident = factor(data_raw$ident, levels = unique(data_raw$ident))
        data_formatted = data.frame(lapply(split(data_raw, f = data_raw$ident), function(x) x$mean_expression))
        colnames(data_formatted) = if (is.null(x_lab)) { unique(data_raw$ident) } else x_lab
        rownames(data_formatted) = y_lab
        return(pheatmap::pheatmap(as.matrix(data_formatted), angle_col = "45", fontsize = 16, main = sig_name, cluster_rows = cluster_rows, cluster_cols = cluster_cols, scale = scale_heatmap, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)))
      }
    } else {
      plot = VlnPlot(sobj, stack = TRUE, flip = TRUE, features = features, group.by = group.by, idents = idents, sort = sort, cols = cols) + NoLegend() + ggtitle(sig_name) + labs(x = group.by) # + scale_x_discrete(labels=x_lab)
      if (as_heatmap) {
        library(pheatmap); library(dplyr)
        data_raw = plot$data %>% group_by(feature, ident) %>% dplyr::summarize(mean_expression = mean(expression))
        data_raw$ident = factor(data_raw$ident, levels = unique(data_raw$ident))
        data_formatted = data.frame(lapply(split(data_raw, f = data_raw$ident), function(x) x$mean_expression))
        if (!is.null(x_lab)) colnames(data_formatted) = x_lab
        rownames(data_formatted) = unique(data_raw$feature)
        plot = pheatmap::pheatmap(as.matrix(data_formatted), angle_col = "45", fontsize = 11, main = sig_name, cluster_rows = cluster_rows, cluster_cols = cluster_cols, scale = scale_heatmap, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100))
      }
    }
  } else {
    gene.set <- rownames(GetAssay(sobj)$data[rownames(GetAssay(sobj)$data) %in% features,])
    diff = setdiff(features, gene.set)
    if (length(diff) > 0) {
      cat(" Warning! Non-matching gene symbols: \n", diff, "\n")
      features = intersect(features, gene.set)
    }
    # get mean expression of features of interest per cell
    mean.exp <- colMeans(GetAssay(sobj)$data[gene.set, ], na.rm = TRUE)
    # add mean expression values in 'sobj@meta.data$gene.set.score'
    if (all(names(x = mean.exp) == rownames(x = sobj@meta.data))) {
      sobj@meta.data$gene.set.score <- mean.exp
    }
    if (is.null(x_lab)) x_lab = if (is.null(group.by)) { levels(Idents(sobj)) } else x_lab = unique(sobj@meta.data[[group.by]])
    plot = VlnPlot(sobj, stack = F, flip = TRUE, features = "gene.set.score", pt.size = 0, group.by = group.by, idents = idents, cols = cols) + NoLegend() + ggtitle(sig_name) + labs(x = group.by) # + scale_x_discrete(labels=x_lab)
  }
  return(plot)
}


calculate_variance_across_groups = function(sobj, group.by = "subtype") {
  data = lapply(unique(sobj[[group.by]])[[group.by]], function(group) {
    if (is.na(group)) next
    expr <- FetchData(sobj, vars = group.by)
    subset_obj = sobj[, which(expr == group)]
    data = AggregateExpression(subset_obj, group.by = "orig.ident")$RNA # pseudobulk per sample
    gene_vars = apply(data, 1, var) # variance among all samples of the group, per gene
    # mean(gene_vars) # average variance considering all genes
    data = data.frame(gene = names(gene_vars), group = group, gene_vars = gene_vars)
    data$group = group; data$gene <- rownames(data)
    return(data)
  })
  data = do.call(rbind, data)
  return(data)
}

reduce_cell_numbers = function(sobjects, cell_limit) {
  for (i in 1:length(sobjects)){
    if (length(colnames(sobjects[[i]])) > cell_limit) {
      message("Reducing cell number in ", sobjects[[i]]@project.name, " from ", 
              length(colnames(sobjects[[i]])), " to ", cell_limit, " cells by random sampling")
      sobjects[[i]] = sobjects[[i]][, sample(colnames(sobjects[[i]]), size = cell_limit, replace = F)]
    }
  }
  return(sobjects)
}

dimplot_standard_qc_metrics = function(sobj, reduction, cells = NULL, pt.size = NULL, raster = NULL, feature_names = NULL) {
  library(patchwork)
  if (is.null(feature_names)) feature_names = c("nFeature_RNA", "nCount_RNA", "percent.mt", "total.mt")
  p1 = FeaturePlot(sobj, reduction = reduction, cells = cells, features = feature_names[1], pt.size = pt.size, raster = raster)
  p2 = FeaturePlot(sobj, reduction = reduction, cells = cells, features = feature_names[2], pt.size = pt.size, raster = raster)
  p3 = FeaturePlot(sobj, reduction = reduction, cells = cells, features = feature_names[3], pt.size = pt.size, raster = raster, cols = c("lightgrey", "red"))
  p4 = FeaturePlot(sobj, reduction = reduction, cells = cells, features = feature_names[4], pt.size = pt.size, raster = raster, cols = c("lightgrey", "red"))
  return(wrap_plots(p1, p2, p3, p4))
}

vlnplot_standard_qc_metrics = function(sobj, reduction) {
  library(patchwork)
  p1 = VlnPlot(sobj, features = "nFeature_RNA", pt.size = 0.5) + NoLegend()
  p2 = VlnPlot(sobj, features = "nCount_RNA", pt.size = 0.5) + NoLegend()
  p3 = VlnPlot(sobj, features = "percent.mt", pt.size = 0.5) + NoLegend()
  p4 = VlnPlot(sobj, features = "total.mt", pt.size = 0.5) + NoLegend()
  return(wrap_plots(p1, p2, p3, p4))
}

# creates new slot in sobj meta.data to annotate malignant and non-malignant clusters for the specified clustering
# EXAMPLE use ||| sobj = store_malignancy_annotation(sobj, c(0,1,2,3,4,7,9,10,12,13,14,17,26), "rpca_clusters", "rpca_malignancy")
store_malignancy_annotation = function(sobj, malig_cluster_ids, cluster_ident, malig_ident_name){
  cluster_labels = sobj@meta.data[[cluster_ident]]
  all_ids = levels(cluster_labels)
  malignancy_labels = rep("non-malignant", length(all_ids))
  malignancy_labels[malig_cluster_ids + 1] = "malignant"
  levels(cluster_labels) = malignancy_labels
  sobj@meta.data[[malig_ident_name]] = cluster_labels
  return(sobj)
}

add_annotation_from_excelsheet = function(sobj, metadata, col_name, annotation_id, sample_id_col_name = "Study ID") {
  subtype_annotation = data.frame(metadata[[col_name]], "s_ids" = metadata[[sample_id_col_name]])
  ljoined_anno = merge(subtype_annotation, data.frame("s_ids" = sobj@meta.data$orig.ident), by = "s_ids", all.x = T, sort = F)
  sobj@meta.data[[annotation_id]] = ljoined_anno[[2]]
  return(sobj)
}

# EXAMPLE use:   sobj = store_cluster_labels(sobj, c("label1", "label2", ...), "rpca_labels", "rpca_clusters")
store_cluster_labels = function(sobj, labels, label_name, clusters_to_annotate = "rpca_clusters") {
  # cluster_ids = unique(sobj@meta.data[[clusters_to_annotate]]) # 1. copy meta.data slot of clusters to be annotated
  # if (length(labels) != length(levels(cluster_ids))) stop("Error: Number of given labels (n=", length(labels), ") does not match number of clusters to be annotated (n=", length(levels(cluster_ids)), ")!")
  # levels(cluster_ids) = labels # 2. rename levels of this factor type
  # sobj@meta.data[[label_name]] = cluster_ids # 3. store in newly create meta.data slot
  # return(sobj)
  message("use: 'sobj$rpca_labels = factor(sobj$rpca_clusters, labels = cluster_labels)' instead")
}

add_absolute_feature_expression = function (object, pattern = NULL, features = NULL, col.name = NULL, assay = NULL) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = features) && !is.null(x = pattern)) {
    warn(message = "Both pattern and features provided. Pattern is being ignored.")
  }
  absolute.featureset <- list()
  layers <- Layers(object = object, search = "counts")
  for (i in seq_along(along.with = layers)) {
    layer <- layers[i]
    features.layer <- features %||% grep(pattern = pattern, x = rownames(x = object[[assay]][layer]), value = TRUE)
    layer.data <- LayerData(object = object, assay = assay, layer = layer)
    layer.sums <- colSums(x = layer.data[features.layer, , drop = FALSE])
    absolute.featureset[[i]] <- layer.sums
  }
  absolute.featureset <- unlist(absolute.featureset)
  if (!is.null(x = col.name)) {
    object <- AddMetaData(object = object, metadata = absolute.featureset, col.name = col.name)
    return(object)
  }
  return(absolute.featureset)
}

compare_jaccard_degs_vs_ref = function(all_degs_observed, reference_sig, min_log2fc = 0, max_p = 0.05, same_length = F){
  j_results = c()
  cat("Calculating Jaccard index for reference set of length:", length(reference_sig), "\n")
  for (cluster_degs in all_degs_observed){
    if (same_length){
      x = cluster_degs[order(cluster_degs$avg_log2FC, decreasing = T),]
      x = x[x$avg_log2FC > 0, ]$gene[1:length(reference_sig)]
    } else {
      x = cluster_degs[cluster_degs$avg_log2FC > min_log2fc & cluster_degs$p_val_adj < max_p, ]$gene
    }
    j = jaccard(reference_sig, x)
    cat("Jaccard index for cluster ", cluster_degs$cluster[1], " (n=", length(x), "): ", j, "\n", sep = "")
    j_results = c(j_results, j)
  }
  cat("Highest match found in cluster", which(j_results == max(j_results))-1, "with Jaccard index",max(j_results), "\n")
  names(j_results) = 0:(length(j_results)-1)
  return(j_results)
}


plot_cluster_composition = function(sobj, group.by, clusters, clustering = "rpca_clusters") {
  data = list()
  meta_data = sobj@meta.data[sobj@meta.data[[clustering]] %in% clusters, ] # meta.data subset of only the requested clusters
  selected_group = unique(sobj@meta.data[[group.by]])
  for (i in selected_group) {
    d = c()
    for (j in clusters) {
      y = length(meta_data[meta_data[[group.by]] == i & meta_data[[clustering]] == j, 1])
      d = c(d, y)
    }
    data[[i]] = d
  }
  df = data.frame(rep(names(data), each = length(clusters))); colnames(df)[1] = "group"; df$n_cells = unlist(data); df$clusters = rep(as.character(clusters), length(unique(df$group)))
  df$clusters <- factor(df$clusters, levels = unique(df$clusters)) # to fixate the given order of clusters
  return(ggplot(df, aes(x=group,y=n_cells, fill = clusters)) + geom_bar(stat="identity", position="stack") + labs(x=group.by) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
}

# plots percentage of cells of specified group contributing to the specified clusters, relative to the total number of cells of the group
plot_relative_subtype_contributions = function(sobj, group.by, clusters = NULL, clustering = "rpca_clusters", exclude_group = NULL) {
  if (is.null(clusters)) clusters = levels(sobj@meta.data[[clustering]])
  data = list()
  meta_data = sobj@meta.data[sobj@meta.data[[clustering]] %in% clusters, ] # meta.data subset of only the requested clusters
  selected_group = unique(sobj@meta.data[[group.by]])
  if (!is.null(exclude_group)) selected_group = setdiff(selected_group, exclude_group)
  for (i in selected_group) {
    d = c()
    for (j in clusters) {
      y = length(meta_data[meta_data[[group.by]] == i & meta_data[[clustering]] == j, 1])
      d = c(d, y)
    }
    data[[i]] = d/length(sobj@meta.data[sobj@meta.data[[group.by]] == i, 1])
  }
  df = data.frame(rep(names(data), each = length(clusters))); colnames(df)[1] = "group"; df$percent_cells = unlist(data); df$clusters = rep(as.character(clusters), length(unique(df$group)))
  df$clusters <- factor(df$clusters, levels = unique(df$clusters)) # to fixate the given order of clusters
  return(ggplot(df, aes(x=group,y=percent_cells, fill = clusters)) + geom_bar(stat="identity", position="stack") + labs(x=group.by) + 
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle("Relative contribution to clusters"))
}

# to compare similarity of reference signatures (e.g. Gojo) and metaprograms
plot_heatmap_mp_vs_ref = function(refs, mps, min_log2fc = 1, max_p = 0.05, cluster_rows = T, scale = "none", cluster_cols = F) {
  library(pheatmap)
  all_intersects = data.frame()
  for (ref_sig in refs) {
    ref_sig = ref_sig[ref_sig$p_val < max_p & ref_sig$avg_logFC > min_log2fc, ]$gene
    j_intersects = c()
    for (mp_sig in mps) {
      j_intersects = c(j_intersects, jaccard(ref_sig, mp_sig))
    }
    all_intersects = rbind(all_intersects, j_intersects)
  }
  data = data.frame(all_intersects)
  colnames(data) = colnames(mps); rownames(data) = names(refs)
  return(pheatmap(data, angle_col = "45", fontsize = 16, main = "Metaprograms vs. Gojo signatures - Jaccard Index", 
                  cluster_rows = cluster_rows, cluster_cols = cluster_cols, scale = scale))
}

plot_heatmap_mps_vs_cells = function(sobj, mps, reduced_cell_n = NA) {
  library(pheatmap)
  # sobj = ScaleData(sobj, features = mps)
  if (!is.na(reduced_cell_n)) { # randomly takes out the specified number of cells
    data = (FetchData(sobj, vars = mps, layer = "data", cells = sample.int(length(colnames(sobj)), reduced_cell_n)))
  } else { data = (FetchData(sobj, vars = mps, layer = "data")) }
  data = t(scale_per_gene(data))
  for (gene in unique(mps[duplicated(mps)])) { # manually adds duplicated genes back to data that were left out by FetchData
    duplicated_genes = which(mps == gene)
    add_index = duplicated_genes[-1]
    for (ind in add_index) {
      data = rbind(data[1:(ind-1),], data[duplicated_genes[1],] ,data[-(1:(ind-1)),])
      rownames(data)[ind] = gene
    }
  }
  if (length(mps) <= 150) {
    show_rownames = T
  } else { show_rownames = F }
  return(pheatmap(data, angle_col = 45, fontsize = 11, main = "Metaprograms", border_color = NA, treeheight_col = 0, 
                  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)[20:100], 
                  show_colnames = F, show_rownames = show_rownames, cutree_rows = 14, cluster_rows = F, silent = T))
}


try_add_module_score <- function(sobj, features, name, ctrl = 100, ctrl_decrement = 5, min_ctrl = 50, nbin = 24, min_bin = 18, verbose = F) {
  current_ctrl <- ctrl  # Startwert für ctrl
  repeat {
    if (current_ctrl < min_ctrl) break
    tryCatch({
      sobj <- Seurat::AddModuleScore(sobj, slot = "data", features = features, name = name, ctrl = current_ctrl)
      break # wenn erfolgreich, beende die Schleife
    }, error = function(e) { current_ctrl <<- current_ctrl - ctrl_decrement; if (verbose) message(e) }) # try with reduced number of ctrl features
  }
  if (current_ctrl >= min_ctrl) {
    if (current_ctrl != ctrl) message("[WARNING] Highest working number of ctrl features lower than desired: ", current_ctrl, " instead of ", ctrl)
    return(sobj)
  } else {
    message("[WARNING] No working number of ctrl features found for min_ctrl = ", min_ctrl, "!")
    current_bin = nbin
    repeat {
      if (current_bin < min_bin) break
      tryCatch({
        sobj <- Seurat::AddModuleScore(sobj, slot = "data", features = features, name = name, ctrl = ctrl, nbin = current_bin)
        break # wenn erfolgreich, beende die Schleife
      }, error = function(e) { current_bin <<- current_bin - 1; if(verbose) message(e) }) # try with reduced number of bins
    }
    if (current_bin >= min_bin) { message("[WARNING] Highest working number of bins lower than desired: ", current_bin, " instead of ", nbin, " (low variance in data?)")
    } else { message("[ERROR] No working number of bins found for min_bin = ", min_bin, "!") }
    return(sobj)
  }
}

# score_names = vector of metadata entries after running AddModuleScores
assign_cells_by_max_score = function(sobj, score_names, meta_data_entry_name) {
  best_fitting_sig = FetchData(sobj, vars = score_names)
  best_fitting_sig = apply(best_fitting_sig, 1, FUN = function(x){colnames(best_fitting_sig)[which.max(x)]})
  sobj@meta.data[[meta_data_entry_name]] = factor(best_fitting_sig, levels = score_names)
  return(sobj)
}

transform_modscores_by_cell = function(sobj, score_names, transformation = c("min_max", "min>0_max", "softmax", "normalisation"), min_score_threshold = 0, transformed_score_names = score_names) {
  scores = sobj@meta.data[, score_names]
  transformed_scores = t(apply(scores, 1, FUN = function(cell) { 
    score = rep(0, length(cell))
    present_MPs = which(cell > min_score_threshold) # consider only MPs as present which exceed the minimum score threshold
    if (transformation == "softmax") {
      warning("not recommended!! softmax is sensitive to the input scores' magnitudes. For the usual range of module scores, all exponentials will be almost the same")
      score[present_MPs] = exp(cell[present_MPs]) / sum(exp(cell[present_MPs])) # softmax on filtered scores
      return(unlist(score))
    }
    if (transformation == "min_max") {
      score[present_MPs] = normalize01(unlist(cell[present_MPs])) # scales all scores of present MPs into [0,1] across cells
      return(score)
    }
    if (transformation == "min>0_max") {
      score[present_MPs] = normalize01(c(0, unlist(cell[present_MPs])))[-1] # scales all scores of present MPs into (0,1] across cells
      return(score)
    }
    if (transformation == "normalisation") {
      score[present_MPs] = cell[present_MPs] / sum(cell[present_MPs])
      return(score)
    }
  }))
  sobj@meta.data[, transformed_score_names] = transformed_scores
  return(sobj)
}

transform_modscores_by_program = function(sobj, score_names, transformation = c("min_max"), min_score_threshold = NULL, transformed_score_names = score_names) {
  scores = sobj@meta.data[, score_names]
  transformed_scores = apply(scores, 2, FUN = function(program) {
    if (!is.null(min_score_threshold)) {
      score = rep(0, length(program))
      pos_cells = which(program > min_score_threshold) # consider only cells as present which exceed the minimum score threshold
      if (transformation == "min_max") {
        score[pos_cells] = normalize01(unlist(program[pos_cells])) # scales all scores of present MPs into [0,1] across programs
        return(score)
      }
    }
    else { # don't filter by threshold here
      if (transformation == "min_max") return(normalize01(program)) # scales all scores of all MPs into [0,1] across programs
    }
  })
  sobj@meta.data[, transformed_score_names] = transformed_scores
  return(sobj)
}


















