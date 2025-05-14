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
