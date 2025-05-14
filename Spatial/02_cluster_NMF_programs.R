library(reshape2)
library(ggplot2)
library(scales)
library(writexl)
library(readxl)
library(pals)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(pheatmap)
library(patchwork)
library(dplyr)
library(Polychrome)

utils_dir = "utils"
source(file.path(utils_dir, "r_utils.R"))
source(file.path(utils_dir, "NMF_utils.R"))

base_out_dir = "NMF"; if(!dir.exists(base_out_dir)) dir.create(base_out_dir, r=T)
input_dir = "NMF/programs by sample"
metadata = read_excel("masterlist_EPN-ST.xlsx", skip = 1); metadata = metadata[grepl("2", metadata$`Analyses*`), ] # subset to only spatial samples
sample_ids = list.files(input_dir); sample_ids = sample_ids[sample_ids %in% metadata$`Study ID`]
sample_colors = as.character(createPalette(N = length(sample_ids), seedcolors = alphabet.colors(26)))
cluster_colors = my_cols

# adapted from Gavish et al. 2023 and Muenter et al. 2025

### ------------------ part 1: Aggregate W matrices (genes x modules) of all samples and all its ranks --------------------------
max_rank = 20 # limit highest rank considered to not inflate number of programs
aggregated_w_matrices = list()
for (s in sample_ids) {
  message("Collecting W matrices of sample ", s)
  in_dir = file.path(input_dir, s)
  
  res.list = readRDS(paste0(in_dir, "/raw_res-list_", s, ".rds")) # loads NMF res for all calculated ranks
  modules.list = Filter(Negate(is.null), lapply(res.list, NMFToModules, gmin = 5))[1:(max_rank-1)]
  comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
  r = names(which(comp == min(comp)))
  r = r[length(r)] # chooses the rank for which the next bigger rank does not add another "useful" (gmin >= 5) module
  message("selecting rank ", r)
  res.list = res.list[1:(as.numeric(r)-1)] # subset the list to keep only the highest rank and all lower ranks
  
  w_matrices_of_sample = data.frame(matrix(nrow = nrow(res.list[[1]])))
  for(j in 1:length(res.list)) {
    w = basis(res.list[[j]]) # extract basis matrix W for current rank
    cols_pre = ncol(w_matrices_of_sample) # get number of columns of dataframe in order to later be able to change the correct column names
    w_matrices_of_sample = cbind(w_matrices_of_sample, w) # add W to the dataframe for the respective sample
    cols = ncol(w_matrices_of_sample)
    n_modules = cols - cols_pre
    colnames(w_matrices_of_sample)[(cols_pre+1):cols] = paste0(s, "_rank", j+1, "_mod", 1:n_modules)
  }
  w_matrices_of_sample = w_matrices_of_sample[,-1]
  aggregated_w_matrices = c(aggregated_w_matrices, list(w_matrices_of_sample)) # add dataframe to the main list aggregating all samples' W matrices
  names(aggregated_w_matrices)[length(aggregated_w_matrices)] = paste0(s, "_ranks_2-", highest_rank)
}
message("Aggregating done for in total ", sum(unlist(lapply(aggregated_w_matrices, length))), " modules of ", length(aggregated_w_matrices), " different samples")
saveRDS(aggregated_w_matrices, file = file.path(base_out_dir, "aggregated_w_matrices.rds"))



### ------------------- part 2: Generate robust programs and metaprograms ----------------------------

# ----------------------------------
# >>>>>>> a) key ROBUST PROGRAM parameters <<<<<<<
intra_min_parameter = 25 # minimal gene overlap between programs of the same sample needed to call them 'robust programs'
intra_max_parameter = 10 # programs of the same sample with HIGHER gene overlap are considered redundant (and therefore removed)
inter_min_parameter = 8 # minimal gene overlap between programs of two different samples, only used if inter_filter = TRUE
# >>>>>>> b) key METAPROGRAM parameters <<<<<<<
min_intersect = 10 # minimal number of genes intersecting between robust programs to be aggregated to metaprogram
min_group = 2 # minimal number of (significant) intersections a robust program must have to qualify for metaprogram initiation
# >>>>>>> c) metaprogram evaluation and re-assignment parameters (optional) <<<<<<<
min_intersect_final = 8 # minimal mean intersecting genes between a robust program and all genes of the initial metaprogram programs to be finally assigned to it
min_mp_size_final = 5 # minimal number of robust programs a final metaprogram must contain
# ----------------------------------

custom_magma = c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10), rev(viridis::magma(323, begin = 0.18)))
out_dir = paste0(base_out_dir, "/MPs ", intra_min_parameter, ".", intra_max_parameter, ".", inter_min_parameter, ".", min_intersect, ".", 
                 min_group, ".", min_intersect_final, ".", min_mp_size_final); if (!dir.exists(out_dir)) dir.create(out_dir)
mp_gen_dir = file.path(out_dir, "MP generation"); if (!dir.exists(mp_gen_dir)) dir.create(mp_gen_dir)
anno_out_dir = file.path(out_dir, "Annotation"); if (!dir.exists(anno_out_dir)) dir.create(anno_out_dir)


# 2.1 --- Select robust NMF programs ---
aggregated_w_matrices = readRDS(file.path(base_out_dir, "aggregated_w_matrices.rds"))
message("Loaded ", sum(unlist(lapply(aggregated_w_matrices, length))), " programs of ", length(aggregated_w_matrices), " samples")

# - nmf_modules = a list; each element is a matrix of all modules of a sample and their top 50 genes 
# - inter_filter = logical; indicates whether programs should be filtered based on their similarity to programs of other samples
nmf_modules = lapply(aggregated_w_matrices, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50])) # get top 50 genes for each NMF program 

# for each sample, select robust NMF programs (i.e. observed using different ranks in the same sample), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
robust_program_names = robust_nmf_programs(nmf_modules, intra_min = intra_min_parameter, intra_max = intra_max_parameter, inter_filter = T, inter_min = inter_min_parameter)  
robust_programs = lapply(nmf_modules, function(x) x[, colnames(x) %in% robust_program_names]) # attach genes back to the identified robust programs
robust_programs = do.call(cbind, robust_programs); colnames(robust_programs) = robust_program_names
robust_programs_original <- robust_programs
if (length(robust_programs) == 0) stop("No robust NMF programs found for the selected filters")
message("Identified ", ncol(robust_programs), " different robust programs (of initially ", sum(unlist(lapply(aggregated_w_matrices, length))), " programs)")


# 2.2 --- Cluster robust NMF programs to generate MPs ---

### Parameters for clustering
Min_intersect_initial = min_intersect         # the minimal intersection required to define the first robust program in a metaprogram (set above)
Min_intersect_metaprogram = min_intersect     # the minimal intersection required to add a new robust program to the forming metaprogram (set above)
Min_group_size = min_group                    # the minimal group size to consider for defining the first robust program in a metaprogram (set above)

robust_program_intersects = apply(robust_programs, 2, function(x) apply(robust_programs, 2, function(y) length(intersect(x,y)))) # calculate gene overlap between programs
# hierarchical clustering of the similarity matrix 
robust_program_intersects_hc = hclust(as.dist(50-robust_program_intersects), method="average") 
robust_program_intersects_hc = reorder(as.dendrogram(robust_program_intersects_hc), colMeans(robust_program_intersects))
robust_program_intersects = robust_program_intersects[order.dendrogram(robust_program_intersects_hc), order.dendrogram(robust_program_intersects_hc)]
Sorted_intersection = sort(apply(robust_program_intersects, 2, function(x) (length(which(x>=Min_intersect_initial))-1)), decreasing = TRUE) # number of significant intersections for every robust programs, sorted
MP_programs = list()   ### Every entry contains the robust programs of a chosen metaprogram
MP_genes = list()
k = 1
robust_program_intersects_original = robust_program_intersects # store for later 

while (Sorted_intersection[1] >= Min_group_size) {
  curr_metaprogram = c(names(Sorted_intersection[1]))
  ### intersection between all remaining robust programs and Genes in MP 
  Genes_MP                    = robust_programs[,names(Sorted_intersection[1])] # Genes in the forming MP are first chosen to be those in the first robust program. Genes_MP always has only 50 genes and evolves during the formation of the metaprogram
  robust_programs             = robust_programs[,-match(names(Sorted_intersection[1]) , colnames(robust_programs))]  # remove selected robust program
  Intersection_with_Genes_MP  = sort(apply(robust_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other robust programs and Genes_MP  
  robust_program_history      = Genes_MP  # has genes in all robust programs in the current metaprogram, for redefining Genes_MP after adding a new robust program 
  ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span beyond the 50th genes, they are sorted according to their robust program score.    
  while (Intersection_with_Genes_MP[1] >= Min_intersect_metaprogram) {
    curr_metaprogram  = c(curr_metaprogram , names(Intersection_with_Genes_MP)[1])
    Genes_MP_temp   = sort(table(c(robust_program_history , robust_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all robust programs in the current metaprogram 
    Genes_at_border = Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene
    if (length(Genes_at_border)>1) {
      ### Sort last genes in Genes_at_border according to maximal robust program gene scores
      ### Run across all robust program in curr_metaprogram and extract robust program scores for each gene
      Genes_curr_robust_program_score = c()
      for (i in curr_metaprogram) {
        current_sample = paste0(sub("_rank.*", "", i), "_ranks")
        list_index = grep(current_sample, names(aggregated_w_matrices))
        Q = aggregated_w_matrices[[list_index]][ match(names(Genes_at_border),toupper(rownames(aggregated_w_matrices[[list_index]])))[!is.na(match(names(Genes_at_border),toupper(rownames(aggregated_w_matrices[[list_index]]))))] ,i] 
        names(Q) = names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(aggregated_w_matrices[[list_index]]))))])  ### sometimes when adding genes the names do not appear 
        Genes_curr_robust_program_score = c(Genes_curr_robust_program_score,  Q )
      }
      Genes_curr_robust_program_score_sort = sort(Genes_curr_robust_program_score , decreasing = TRUE)
      Genes_curr_robust_program_score_sort = Genes_curr_robust_program_score_sort[unique(names(Genes_curr_robust_program_score_sort))]   
      Genes_MP_temp             = c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_robust_program_score_sort))
    } else {
      Genes_MP_temp = names(Genes_MP_temp)[1:50] 
    }
    robust_program_history     = c(robust_program_history , robust_programs[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP        = Genes_MP_temp[1:50]
    robust_programs    = robust_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(robust_programs))]  # remove selected NMF
    Intersection_with_Genes_MP = sort(apply(robust_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  }
  MP_programs[[paste0("MP_",k)]] = curr_metaprogram
  MP_genes[[paste0("MP_",k)]] = Genes_MP
  k = k+1
  robust_program_intersects = robust_program_intersects[-match(curr_metaprogram, rownames(robust_program_intersects)), -match(curr_metaprogram, colnames(robust_program_intersects))]  # Remove current chosen metaprogram
  Sorted_intersection = sort(apply(robust_program_intersects , 2, function(x) (length(which(x>=Min_intersect_initial))-1)), decreasing = TRUE)   # Sort intersection of remaining NMFs not included in any of the previous metaprograms
  message("new metaprogram formed containing ", length(curr_metaprogram), " robust programs")
  # print(paste0("sorted intersection[1] now: ", Sorted_intersection[1]))
}

# 2.3 --- Re-sort robust programs to metaprograms (OPTIONAL) ---
# a) re-assign all robust programs to the best matching MP (if different to initial assignment), then re-calculate genes
final_MP_programs = replicate(length(MP_programs), vector(), simplify = FALSE); names(final_MP_programs) = names(MP_programs)
for (prog in colnames(robust_programs_original)) {
  # mp_intersects = sapply(MP_genes, function(mp) length(intersect(mp, robust_programs_original[, prog]))) # intersection to MP regarding MP genes
  mp_intersects <- sapply(MP_programs, function(mp_progs) mean(robust_program_intersects_original[prog, mp_progs[prog != mp_progs]])) # intersection to MP regarding MP programs
  candidate_mps = mp_intersects[mp_intersects >= min_intersect_final]
  if (length(candidate_mps) > 0) {
    final_mp = names(candidate_mps[which.max(candidate_mps)])
    final_MP_programs[[final_mp]] = c(final_MP_programs[[final_mp]], prog)
  }
}
final_MP_programs = Filter(Negate(is.null), lapply(final_MP_programs, function(mp) if (length(mp) >= min_mp_size_final) mp)) # filter out metaprograms that are too small
# b) re-calculate final MP genes
final_MP_genes = lapply(final_MP_programs, calculate_dominant_genes)
# c) re-order robust programs within metaprograms by intersection with each other, starting with the program with the highest gene intersection to the MP genes
for (mp in names(final_MP_programs)) {
  mp_mat = robust_program_intersects_original[final_MP_programs[[mp]], final_MP_programs[[mp]]]
  prog_similarity = sapply(rownames(mp_mat), function(p) length(intersect(final_MP_genes[[mp]], robust_programs_original[, p])))
  final_MP_programs[[mp]] = names(prog_similarity)[order(prog_similarity, decreasing = T)] # order by gene similarity to MP genes
}
MP_programs = final_MP_programs; MP_genes = final_MP_genes

# 2.4 --- save all robust programs, intersections, metaprograms and genes (as .rds and in excel) ---
saveRDS(robust_programs_original, file.path(mp_gen_dir, "robust_programs.RDS"))
saveRDS(robust_program_intersects_original, file.path(mp_gen_dir, "robust_programs_intersects.RDS"))
saveRDS(MP_genes, file.path(mp_gen_dir, "metaprogram_genes.RDS"))
saveRDS(MP_programs, file.path(mp_gen_dir, "metaprogram_programs.RDS"))

robust_program_intersects_excel = data.frame(robust_programs = rownames(robust_program_intersects_original), robust_program_intersects_original)
MP_excel = as.data.frame(MP_genes)
n.clus = sapply(MP_programs, length)
n.max = max(n.clus)
MP_programs_temp = MP_programs
for (i in 1:length(MP_programs_temp)){
  MP_programs_temp[[i]] = c(MP_programs_temp[[i]], rep(NA, (n.max-n.clus)[i]))
}
remaining_NMF_excel = data.frame(robust_programs[, !colnames(robust_programs) %in% unlist(MP_programs)])
write_xlsx(list("MP genes" = MP_excel, "All robust NMF Programs" = as.data.frame(robust_programs_original), "Robust programs intersects" = robust_program_intersects_excel, 
                "MP programs" = as.data.frame(MP_programs_temp), "Robust programs not in MPs" = remaining_NMF_excel), 
           path = file.path(out_dir, "NMF_metaprograms_main.xlsx"))


### ------------------- part 3: Plot heatmaps of metaprograms ----------------------------
# a) main MP heatmap but with reordered metaprograms (optional)
mps_reordered = MP_programs[order(sapply(MP_programs, length), decreasing = T)]
reorder = 1:length(mps_reordered)
for (i in 1:(length(mps_reordered)-1)) {
  max_intersect = 0; max_index = i+1
  for (j in (i+1):length(mps_reordered)) {
    intersect = sum(unlist(robust_program_intersects_original[mps_reordered[[reorder[i]]], mps_reordered[[reorder[j]]]]))
    intersect = intersect / length(mps_reordered[[reorder[j]]])
    if (intersect > max_intersect) { max_intersect = intersect; max_index = j }
  }
  reorder[c(i+1, max_index)] <- reorder[c(max_index, i+1)] # swap the index of the top intersecting MP to the (relative) front
}
mps_reordered = mps_reordered[reorder]
inds_sorted = c()
for (j in 1:length(mps_reordered)) inds_sorted = c(inds_sorted , match(mps_reordered[[j]] , colnames(robust_program_intersects_original)))
robust_program_intersects_meltI = reshape2::melt(robust_program_intersects_original[inds_sorted,inds_sorted])
g = plot_NMF_heatmap(robust_program_intersects_meltI, xlab = paste0("robust NMF program (n = ",length(inds_sorted), ")"), limits = c(5,25))
colors = cluster_colors[1:length(mps_reordered)]; names(colors) = names(mps_reordered)
color_data = data.frame(y = unique(robust_program_intersects_meltI$Var2), Color = unlist(sapply(names(mps_reordered), function(x) rep(colors[[x]], length(mps_reordered[[x]])))), name = unlist(sapply(names(mps_reordered), function(x) rep(x, length(mps_reordered[[x]])))))
color_data <- color_data %>% group_by(Color) %>% dplyr::mutate(name = if_else(row_number() == ceiling(dplyr::n()/2), name, "")) %>% ungroup()
anno = ggplot(data = color_data, aes(x = 0, y = y)) + geom_tile(aes(fill = Color)) + scale_fill_identity() + theme_void() +
  theme(axis.title = element_blank(), axis.text.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(1,0,1,1), "cm")) + labs(x = NULL, y = NULL) + scale_y_discrete(breaks = color_data$y, labels = color_data$name)
g <- anno + g + plot_layout(widths = c(1, 14))
ggsave(g, filename = "03d_heatmap_MPs_reordered-by-similarity.png", width = 10.5, height = 8, path = mp_gen_dir)


# b) main heatmap and heatmap with all robust programs
inds_sorted = c()
for (j in 1:length(MP_programs)) inds_sorted = c(inds_sorted , match(MP_programs[[j]] , colnames(robust_program_intersects_original)))
inds_new = c(inds_sorted, which(is.na( match(1:dim(robust_program_intersects_original)[2],inds_sorted)))) ### clustered robust programs will appear first, and the latter are the robust programs that were not clustered
robust_program_intersects_meltI_NEW = reshape2::melt(robust_program_intersects_original[inds_new,inds_new])
g1 = plot_NMF_heatmap(robust_program_intersects_meltI_NEW, xlab = paste0("All robust NMF programs (n = ",ncol(robust_program_intersects_original), ")"), ylab = "robust NMF program", limits = c(5,25))
ggsave(g1, filename = "02_heatmap_all_robust_programs.png", width = 9, height = 8, path = mp_gen_dir)

# --> main heatmap <--
robust_program_intersects_meltI = reshape2::melt(robust_program_intersects_original[inds_sorted,inds_sorted])
g2 = plot_NMF_heatmap(robust_program_intersects_meltI, xlab = paste0("robust NMF program (n = ",length(inds_sorted), ")"), ylab = "robust NMF program", limits = c(5,25))
colors = cluster_colors[1:length(MP_programs)]; names(colors) = names(MP_programs)
color_data = data.frame(y = unique(robust_program_intersects_meltI$Var2), Color = unlist(sapply(names(MP_programs), function(x) rep(colors[[x]], length(MP_programs[[x]])))), name = unlist(sapply(names(MP_programs), function(x) rep(x, length(MP_programs[[x]])))))
color_data <- color_data %>% group_by(Color) %>% dplyr::mutate(name = if_else(row_number() == ceiling(dplyr::n()/2), name, "")) %>% ungroup()
anno = ggplot(data = color_data, aes(x = 0, y = y)) + geom_tile(aes(fill = Color)) + scale_fill_identity() + theme_void() +
  theme(axis.title = element_blank(), axis.text.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(1,0,1,1), "cm")) + labs(x = NULL, y = NULL) + scale_y_discrete(breaks = color_data$y, labels = color_data$name)
ggsave(anno + g2 + labs(y = NULL) + plot_layout(widths = c(1, 14)), filename = "01_heatmap_MPs_main.png", width = 10, height = 8, path = mp_gen_dir)

cols = c(colorRampPalette(c("white", rev(viridis::mako(323, begin = 0.13))[1]))(20), rev(viridis::mako(323, begin = 0.1)))
plot_NMF_heatmap(robust_program_intersects_meltI, xlab = paste0("robust NMF program (n = ",length(inds_sorted), ")"), ylab = "robust NMF program", limits = c(5,25), cols = cols)
ggsave(anno + g2 + labs(y = NULL) + plot_layout(widths = c(1, 14)), filename = "01_heatmap_MPs_main_alternative-color.png", width = 10, height = 8, path = mp_gen_dir)


# c) labeled main heatmaps
ncols = length(MP_programs) # create color vector for the axis labels to be able to differentiate between the MPs
my_cols = if (ncols <= 18) cluster_colors else cols25(n = ncols)
color_vector = c()
for (i in 1:ncols) color_vector = c(color_vector, rep(my_cols[i], length(MP_programs[[i]])))  # repeat the respective color from the color list as many times as there are samples in the respective MP
color_vector_long = c(color_vector, rep("black", (length(inds_new) - length(color_vector)))) # create a longer vector for the heatmap including all NMF program
g = g1 + theme(axis.text.y = element_text(color = color_vector_long, size = 8), axis.title.x = element_text()) + xlab(paste0("Number of MPs: ", length(MP_programs), "   -   Robust NMF Programs: ", length(inds_new)))
ggsave(g, filename = "03a_heatmap_MPs_labeled+unselected-MPs.png", width = 15, height = 13, path = mp_gen_dir)

# by sample
sample_names = sub("_rank.*$", "", unlist(MP_programs)); names(sample_names) = NULL
color_palette = sample_colors[1:length(unique(sample_names))]; names(color_palette) = unique(sample_names); sample_colors = color_palette[sample_names]
color_data = data.frame(x = unique(robust_program_intersects_meltI$Var2), Color = sample_colors)
anno = ggplot(data = color_data, aes(x = x, y = 0)) + geom_tile(aes(fill = Color)) + scale_fill_identity() + theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(1,1,0,1), "cm")) + labs(x = NULL, y = NULL)
g <- anno + g2 + plot_layout(heights = c(1, 14)) # combine annotation bar and heatmap
ggsave(g, filename = "03b_heatmap_MPs_labeled-by-sample.png", width = 14.6, height = 14, path = mp_gen_dir)

# by patient
patients_names = sub("_rank.*$", "", unlist(MP_programs)); names(patients_names) = NULL
patients_names = ifelse(grepl("[0-9][a-zA-Z]$", patients_names), sub("([0-9])[a-zA-Z]$", "\\1", patients_names), patients_names)
color_palette = sample_colors[1:length(unique(patients_names))]; names(color_palette) = unique(patients_names); patient_colors = color_palette[patients_names]
color_data = data.frame(x = unique(robust_program_intersects_meltI$Var2), Color = patient_colors)
anno = ggplot(data = color_data, aes(x = x, y = 0)) + geom_tile(aes(fill = Color)) + scale_fill_identity() + theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(1,1,0,1), "cm")) + labs(x=NULL, y=NULL)
g <- anno + g2 + plot_layout(heights = c(1, 14)) # combine annotation bar and heatmap
ggsave(g, filename = "03c_heatmap_MPs_labeled-by-patient.png", width = 14.6, height = 14, path = mp_gen_dir)



### ------------------- part 4: Plot composition of metaprograms ----------------------------
# i. composition of MPs regarding samples
samples_per_metaprogram = lapply(MP_programs, function(x){sub("_rank.*$", "",x)})
data = data.frame(samples = rep(unique(sample_ids), length(samples_per_metaprogram)))
data$metaprograms = as.factor(rep(1:length(samples_per_metaprogram), each = length(unique(sample_ids))))
data$n_percent = unlist(lapply(samples_per_metaprogram, function(y) { unlist(lapply(unique(sample_ids), function(x) { sum(x == y)/length(y) })) }))
p = ggplot(data, aes(x=metaprograms, y=n_percent, fill = samples)) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.key.size = unit(4, "mm")) + 
  scale_fill_manual(values = sample_colors) + geom_col() + scale_y_continuous(breaks = c(0, 0.5, 1)) + coord_cartesian(ylim = c(0, 1))
ggsave(p, filename = "04a_barplot_MPs_by-samples.png", width = 6, height = 4, path = mp_gen_dir)

# ii. composition of MPs regarding subtypes and subclass
data$subtype = metadata[match(data$samples, metadata$`Study ID`), ]$Subtype
data$subclass = metadata[match(data$samples, metadata$`Study ID`), ]$`Subclass (classifier)`
p1 = ggplot(data, aes(x=metaprograms, y=n_percent, fill = subtype)) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.key.size = unit(4, "mm")) +
  scale_fill_manual(values = cluster_colors) + geom_col() + scale_y_continuous(breaks = c(0, 0.5, 1)) + coord_cartesian(ylim = c(0, 1))
p2 = ggplot(data, aes(x=metaprograms, y=n_percent, fill = subclass)) + geom_bar(stat="identity", position="stack") + theme_classic() + theme(legend.key.size = unit(4, "mm")) +
  scale_fill_manual(values = cluster_colors) + geom_col() + scale_y_continuous(breaks = c(0, 0.5, 1)) + coord_cartesian(ylim = c(0, 1))
ggsave(p1 + p2, filename = "04b_barplot_MPs_by-subtype.png", width = 7, height = 2.5, path = mp_gen_dir)

# iii. contributions of samples and patients to MPs
mps_per_sample = lapply(unique(sample_ids), FUN = function(x) sapply(samples_per_metaprogram, FUN = function(mp) sum(x == mp)))
data$programs_contributed = as.vector(t(do.call(cbind, mps_per_sample))) # transform to matrix and transpose to match data order
p = ggplot(data, aes(x=samples, y=programs_contributed, fill = metaprograms)) + geom_bar(stat="identity", position="stack") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.key.size = unit(4, "mm")) + scale_fill_manual(values = cluster_colors)
ggsave(filename = "04c_barplot_samples_by_MP-contributions.png", plot = p, width = 9, height = 2.5, path = mp_gen_dir)

# iv. contributions of subtypes and subtypes_precise to MPs (biological variables)
subtype_sizes = sapply(unique(data$subtype), FUN = function(x) sum(x == metadata[match(data$samples, metadata$`Study ID`), ]$Subtype, na.rm = T))
subtype_precise_sizes = sapply(unique(data$subclass), FUN = function(x) sum(x == metadata[match(data$samples, metadata$`Study ID`), ]$`Subclass (classifier)`, na.rm = T))
data$rel_programs_contributed = data$programs_contributed / subtype_sizes[match(data$subtype, unique(data$subtype))]
p1 = ggplot(na.omit(data), aes(x=subtype, y=rel_programs_contributed, fill = metaprograms)) + geom_bar(stat="identity", position="stack") + theme_classic() +
  scale_fill_manual(values = cluster_colors) + ylab("programs contributed / subtype size") + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.key.size = unit(4, "mm"))
data$rel_programs_contributed = data$programs_contributed / subtype_precise_sizes[match(data$subclass, unique(data$subclass))]
p2 = ggplot(na.omit(data), aes(x=subclass, y=rel_programs_contributed, fill = metaprograms)) + geom_bar(stat="identity", position="stack") + theme_classic() +
  scale_fill_manual(values = cluster_colors) + ylab("programs contributed / subclass size") + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.key.size = unit(4, "mm"))
ggsave(p1 + p2, filename = "04d_barplot_subtypes_by_MP-contributions.png", width = 7, height = 4, path = mp_gen_dir)

# v. composition of MPs regarding sex/gender and age (biological variables)
data$sex = metadata[match(data$samples, metadata$`Study ID`), ]$Sex; data$age = metadata[match(data$samples, metadata$`Study ID`), ]$`Age at resection`
p1 = ggplot(data, aes(x=metaprograms, y=n_percent, fill = sex)) + geom_bar(stat="identity", position="stack") + theme_classic() +
  scale_fill_manual(values = cluster_colors) + geom_col() + scale_y_continuous(breaks = c(0, 0.5, 1)) + coord_cartesian(ylim = c(0, 1))
p2 = ggplot(data, aes(x=metaprograms, y=n_percent, fill = age)) + geom_bar(stat="identity", position="stack") + theme_classic() +
  scale_fill_manual(values = cluster_colors) + geom_col() + scale_y_continuous(breaks = c(0, 0.5, 1)) + coord_cartesian(ylim = c(0, 1))
ggsave(p1, filename = "04e_barplot_MPs_by-sex-and-age.png", width = 6, height = 2.5, path = mp_gen_dir)

                               
# ------------------------ part 5: Annotate MPs ----------------------------
metaprograms = read_excel(file.path(out_dir, "NMF_metaprograms_main.xlsx"), sheet = "MP genes")

# a) annotate with EPN-ST signatures
# heatmap_col1 = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
heatmap_col2 = colorRampPalette(c("white", rev(magma(20))[1:17]))(100)
curated_ST_EPN_signatures = as.list(read.xlsx("curated_ST-EPN_signatures.xlsx"))
jac_mtx = compare_jaccard_ref_vs_sig(metaprograms, curated_ST_EPN_signatures); jac_mtx = reorder_similarity_matrix_diagonally(jac_mtx, reorder_rows = T)
pheatmap::pheatmap(jac_mtx, color = heatmap_col1, angle_col = "45", fontsize = 22, scale = "none", main = paste0("MPs vs. EPN-ST signatures by jaccard ind."),
                   cluster_cols = F, cluster_rows = F, filename = file.path(anno_out_dir, "01_heatmap_MPs_vs_EPN-ST-sigs_jaccard_unscaled.png"), width = 10, height = 8)
# customised:
heatmap_col1 = colorRampPalette(c("white", "#e0e1dd", "#c3cbd6", "#9eb8d9", "#778da9", "#415a77", "#1b263b", "black"))(100)
jac_mtx = jac_mtx[rowSums(jac_mtx) > 5, ]
pheatmap::pheatmap(t(jac_mtx), color = heatmap_col1, angle_col = "90", fontsize = 22, scale = "none", main = paste0("MPs vs. EPN-ST signatures by jaccard ind."),

# b) annotate with GO terms
enrich_signatures_GO(metaprograms, anno_out_dir)
enrich_signatures_GO_differential(metaprograms, anno_out_dir, save_combined = T, save_individually = T, GO_label_width = 50, width = 9, height = 10, x_lab_size = 10)

# c) (optional) annotate with reference signatures: Tirosh 3CA metaprograms
anno_out_dir = file.path(anno_out_dir, "Tirosh ref"); if (!dir.exists(anno_out_dir)) dir.create(anno_out_dir, r = T)
heatmap_col1 = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100); heatmap_col2 = colorRampPalette(c("white", rev(magma(20))[1:17]))(100)
ref_all_tirosh = read_tirosh()
for (i in 1:length(ref_all_tirosh)) {
  jac_mtx = compare_jaccard_ref_vs_sig(metaprograms, ref_all_tirosh[[i]]); jac_mtx = reorder_similarity_matrix_diagonally(jac_mtx, reorder_rows = T)
  mps_name = names(ref_all_tirosh)[i]
  pheatmap::pheatmap(jac_mtx, color = heatmap_col2, angle_col = "45", fontsize = 22, scale = "none", main = paste0("MPs vs. Tirosh ", mps_name, " MPs by jaccard ind."),
                     cluster_cols = F, cluster_rows = F, filename = paste0(anno_out_dir, "/01_heatmap_tirosh_", mps_name, "-mps_jaccard_unscaled.png"), width = 14, height = 15)
  pheatmap::pheatmap(jac_mtx, color = heatmap_col1, angle_col = "45", fontsize = 22, scale = "column", main = paste0("MPs vs. Tirosh ", mps_name, " MPs by jaccard ind. (cs)"),
                     cluster_cols = F, cluster_rows = F, filename = paste0(anno_out_dir, "/02_heatmap_tirosh_", mps_name, "-mps_jaccard_col-scaled.png"), width = 14, height = 15)
}








# ------------------------ part 6: Plot customly ordered final heatmap ----------------------------

new_order = c(1,13,4,11, 8,5,3,6,7,9,10,12, 2) # define custom order of MPs here, based on biological similarity and to give better overview

MP_programs_reordered = MP_programs[new_order]; MP_genes_reordered = MP_genes[new_order]

inds_sorted = c()
for (j in 1:length(MP_programs)) inds_sorted = c(inds_sorted , match(MP_programs_reordered[[j]] , colnames(robust_program_intersects_original)))
robust_program_intersects_meltI = reshape2::melt(robust_program_intersects_original[inds_sorted,inds_sorted])
g = plot_NMF_heatmap(robust_program_intersects_meltI, xlab = paste0("All robust NMF programs (n = ", length(inds_sorted), ")"), ylab = NULL, limits = c(5,25))
colors = cluster_colors[1:length(MP_programs_reordered)]; names(colors) = names(MP_programs_reordered)
color_data = data.frame(y = unique(robust_program_intersects_meltI$Var2), Color = unlist(sapply(names(MP_programs_reordered), function(x) rep(colors[[x]], length(MP_programs_reordered[[x]])))), name = unlist(sapply(names(MP_programs_reordered), function(x) rep(x, length(MP_programs_reordered[[x]])))))
color_data <- color_data %>% group_by(Color) %>% dplyr::mutate(name = if_else(row_number() == ceiling(dplyr::n()/2), name, "")) %>% ungroup()
anno = ggplot(data = color_data, aes(x = 0, y = y)) + geom_tile(aes(fill = Color)) + scale_fill_identity() + theme_void() +
  theme(axis.title = element_blank(), axis.text.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(1,0,1,1), "cm")) + labs(x = NULL, y = NULL) + scale_y_discrete(breaks = color_data$y, labels = color_data$name)
g <- anno + g + plot_layout(widths = c(1, 14))
ggsave(g, filename = "heatmap_MPs_custom_order.png", width = 10.5, height = 8, path = mp_gen_dir)


# --- if certain to use the new order, replot all previous plots and update order of MPs in excel ---
# 1. change order in labeled excel table
reorder_excel_sheet_columns(file.path(out_dir, "NMF_metaprograms_main.xlsx"), sheets = c("MP genes", "MP programs"), new_order)

# 2. apply labels in new order from excel to MP_genes and MP_programs
metaprograms = read_excel(file.path(out_dir, "NMF_metaprograms_main.xlsx"), sheet = "MP genes")
MP_programs = MP_programs_reordered; names(MP_programs) = colnames(metaprograms)
MP_genes = MP_genes_reordered; names(MP_genes) = colnames(metaprograms)
saveRDS(MP_genes, file.path(mp_gen_dir, "metaprogram_genes.RDS"))
saveRDS(MP_programs, file.path(mp_gen_dir, "metaprogram_programs.RDS"))

# 3. adapt dirs to save new order-based plots in new folder (change only the two subfolders)
mp_gen_dir = file.path(out_dir, "Custom order/MP generation"); if (!dir.exists(mp_gen_dir)) dir.create(mp_gen_dir, r = T)
anno_out_dir = file.path(out_dir, "Custom order/Annotation"); if (!dir.exists(anno_out_dir)) dir.create(anno_out_dir)
# ------

