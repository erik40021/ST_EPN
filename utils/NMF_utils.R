run_NMF = function(sobj_path, out_path, range, n_features) {
  message("Running NMF for sample ", s)
  out_dir = file.path(out_path, s); if (!dir.exists(out_dir)) { dir.create(out_dir, r = T) }
  
  # 1. ---- load and prepare data -----
  sobj = readRDS(file = paste0(sobj_path, "/sobj", s, ".rds"))
  data = as.matrix(Seurat::GetAssayData(sobj, layer = 'scale.data'))
  data = data[Seurat::VariableFeatures(sobj)[1:n_features], ]
  data[data < 0] = 0  # set negative values to zero
  data = data[apply(data, 1, var) > 0,] # subset on genes with a variance above zero
  message("[", s, "] Extracted matrix of dimensions ", dim(data)[1], " (genes) x ", dim(data)[2], " (spots)")
  
  # 2. ---- run NMF -----
  message("running NMF for ranks ", range[1], "-", tail(range, 1), " | start time: ", format(Sys.time(), "%d.%m. %X"))
  res.list = lapply(range, function(r) {
      message(r, " (", format(Sys.time(), "%X"), ")")
      nmf(data, rank = r, nrun = 1, seed = "ica", method = "nsNMF")
  })
  names(res.list) = range
  saveRDS(res.list, file = paste0(out_dir, "/raw_res-list_", s, ".rds"))
}





plot_NMF_heatmap = function(data, xlab = NULL, ylab = NULL, limits = c(2, 25), cols = custom_magma) {
  midpoint = limits[1] + (limits[2] - limits[1]) / 2
  p = ggplot(data, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + geom_tile() + 
    scale_color_gradient2(limits=limits, low=cols[1:111], mid=cols[112:222], high=cols[223:333], midpoint = midpoint, oob=squish, name="Similarity\n(Jaccard index)") +                                
    scale_fill_gradient2(limits=limits, low=cols[1:111], mid=cols[112:222], high=cols[223:333], midpoint = midpoint, oob=squish, name="Similarity\n(Jaccard index)")  +
    theme(axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), 
          axis.title = element_text(size = 16), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom",
          axis.title.x = element_text(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) + xlab(xlab) + ylab(ylab)
  return(p)
}


calculate_dominant_genes = function(programs) {
  if (length(programs) == 1) return(robust_programs_original[, programs])
  genes_mp_freq = sort(table(robust_programs_original[, programs]), decreasing = TRUE)
  Genes_at_border = genes_mp_freq[which(genes_mp_freq == genes_mp_freq[50])] # genes with overlap equal to the 50th gene
  if (length(Genes_at_border) > 1) {
    Genes_curr_robust_program_score = c()
    for (i in programs) {
      current_sample = paste0(sub("_rank.*", "", i), "_ranks")
      list_index = grep(current_sample, names(aggregated_w_matrices))
      Q = aggregated_w_matrices[[list_index]][ match(toupper(names(Genes_at_border)),toupper(rownames(aggregated_w_matrices[[list_index]])))[!is.na(match(toupper(names(Genes_at_border)),toupper(rownames(aggregated_w_matrices[[list_index]]))))] ,i]
      names(Q) = names(Genes_at_border[!is.na(match(toupper(names(Genes_at_border)),toupper(rownames(aggregated_w_matrices[[list_index]]))))])  ### sometimes when adding genes the names do not appear
      Genes_curr_robust_program_score = c(Genes_curr_robust_program_score, Q)
    }
    Genes_curr_robust_program_score_sort = sort(Genes_curr_robust_program_score, decreasing = TRUE)
    Genes_curr_robust_program_score_sort = Genes_curr_robust_program_score_sort[unique(names(Genes_curr_robust_program_score_sort))]
    final_genes = c(names(genes_mp_freq[which(genes_mp_freq > genes_mp_freq[50])]), names(Genes_curr_robust_program_score_sort))[1:50]
  } else { final_genes = names(genes_mp_freq)[1:50] }
  if (length(final_genes[is.na(final_genes)]) > 0) stop("NA included in genes of MP")
  return(final_genes)
}

# original from Gavish et al. 2023
NMFToModules = function(res, gmin = 5){
  scores = basis(res)
  coefs = coefficients(res)
  # Remove if fewer than gmin genes
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)) ranks_y[ranks_x[,i] > 1,i] = Inf
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  l = sapply(modules, length)
  keep = (l >= gmin)
  if (sum(keep) < 2) { message("Less than 2 modules found in NMF res of rank ", dim(res)[3], " with sufficient number of genes (> ", gmin, ")"); return(NULL)}
  scores = scores[, keep]
  coefs = coefs[keep, ]
  
  # Find modules
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  names(modules) = sapply(modules, '[', 1)
  names(modules) = paste('m', names(modules), sep = '_')
  names(modules) = gsub('-','_',names(modules))
  return(modules)
}

# original from Gavish et al. 2023
# Function for selecting robust non-negative matrix factorization (NMF) programs
# - nmf_programs = a list; each element contains a matrix with NMF programs (top 50 genes) generated for a specific sample using different NMF factorization ranks. 
# - intra_min = minimum overlap with a program from the same sample (for selecting robust programs)
# - intra_max = maximum overlap with a program from the same sample (for removing redundant programs)
# - inter_filter = logical; indicates whether programs should be filtered based on their similarity to programs of other cell lines
# - inter_min = minimum overlap with a program from another sample 
# -> returns a character vector with the NMF programs selected
robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10) {
  message("started filtering of robust programs with initially ", sum(sapply(nmf_programs, ncol)), " programs")
  # Select NMF programs based on the minimum overlap with other NMF programs from the same sample
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
  names(nmf_sel) <- names(nmf_programs)
  message("programs remaining after filtering for intra_min >= ", intra_min, ": ", sum(sapply(nmf_sel, ncol)))
  
  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same sample and
  # ii) the minimum overlap with programs from another sample
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs
  
  final_filter <- NULL
  for(i in names(nmf_sel)) {
    sample_id = sub("_rank.*", "", i)
    a <- inter_intersect[grep(sample_id, colnames(inter_intersect), invert = T), grep(sample_id, colnames(inter_intersect))]
    b <- sort(apply(a, 2, max), decreasing = T) # for each sample, ranks programs based on their maximum overlap with programs of other samples
    if(inter_filter==T) b <- b[b>=inter_min] # selects programs with an intersection of at least inter_min (e.g. 10) to other samples' programs, all others are thrown out
    if(length(b) > 1) {
      c <- names(b[1]) 
      for(y in 2:length(b)) {
        if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have an intersection smaller than or equal to intra_max (e.g. 10) with a previously selected program of the same sample
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  message("programs remaining after filtering for intra_max <= ", intra_max, " and inter_min >= ", inter_min ,": ", length(unique(final_filter)))
  return(unique(final_filter))                                                      
}
