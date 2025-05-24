library(Seurat)
library(patchwork)
library(readxl)
library(ggplot2)
library(dplyr)

source("utils/spatial_utils.R")
source("utils/spatial_score_utils.R")

cloud_path = ""
metadata = read_excel("masterlist_EPN-ST.xlsx", skip = 1); metadata = metadata[grepl("2", metadata$`Analyses*`), ]
sample_ids = metadata$`Study ID`

base_dir = ""
metaprograms = readxl::read_excel("spatial_NMF_metaprograms_main.xlsx", sheet = "MP genes (final)")
mp_names = colnames(metaprograms)
win_sizes = c(5, 8, 11)


# ---------------------------- PART 0: extract data (only one time) -------------------------

all_spots_positions = list(); all_spots_programs_comp_norm = list(); all_spots_programs_comp = list(); all_spots_programs_comp_spot_norm = list()
all_norm_malignancy_scores = list(); all_raw_malignancy_scores = list()
for (s in sample_ids) {
  message("extracting data from ", s)
  sstobj = readRDS(file = paste0("sstobj_", s, ".rds"))
  raw_data_dir = file.path(cloud_path, metadata[metadata$`Study ID` == s, ]$`KK code`, "outs")
  spots_positions = read.csv(list.files(file.path(raw_data_dir, "spatial"), pattern = "tissue_positions.*\\.csv$", full.names = T), header = F, skip = 1)
  row.names(spots_positions) <- spots_positions$V1
  spots_programs_comp = FetchData(sstobj, vars = paste0("MP", 1:ncol(metaprograms)))
  spots_programs_comp$barcodes = rownames(spots_programs_comp); colnames(spots_programs_comp) = c(mp_names, "barcodes")
  spots_programs_comp_spot_norm = FetchData(sstobj, vars = paste0("MP", 1:ncol(metaprograms), "_spot-norm")); colnames(spots_programs_comp_spot_norm) = mp_names
  spots_programs_comp_norm = FetchData(sstobj, vars = paste0("MP", 1:ncol(metaprograms), "_mp-norm")) # normalised by program
  spots_programs_comp_norm$barcodes = rownames(spots_programs_comp_norm); colnames(spots_programs_comp_norm) = c(mp_names, "barcodes")
  all_spots_positions[[s]] = spots_positions; all_spots_programs_comp_norm[[s]] = spots_programs_comp_norm; 
  all_spots_programs_comp[[s]] = spots_programs_comp; all_spots_programs_comp_spot_norm[[s]] = spots_programs_comp_spot_norm
  all_norm_malignancy_scores[[s]] = norm_malignancy_score; all_raw_malignancy_scores[[s]] = raw_malignancy_score
}
saveRDS(all_spots_positions, file = file.path(base_dir, "all_spots_positions.rds")) # <- used for all spatial scores
saveRDS(all_spots_programs_comp_norm, file = file.path(base_dir, "all_spots_programs_comp_norm.rds")) # <- used for complexity score and downstream analyses
saveRDS(all_spots_programs_comp, file = file.path(base_dir, "all_spots_programs_comp.rds")) # <- used for association and zone abundances (PART 3)
saveRDS(all_spots_programs_comp_spot_norm, file = file.path(base_dir, "all_spots_programs_comp_spot_norm.rds")) # <- used to define hyp-classes
saveRDS(all_norm_malignancy_scores, file =  file.path(base_dir, "../all_norm_malignancy_scores.rds")) # <- used to differentiate malignant from non-malignant spots ('CNA' score)
saveRDS(all_raw_malignancy_scores, file =  file.path(base_dir, "../all_raw_malignancy_scores.rds")) # <- same but not normalised, to differentiate magnitudes of the score across samples better


# ------------------------------- PART 1: calculate spatial coherence and complexity -----------------------

# 1.1 calculate spatial coherence, based on normalised or raw MP scores (optional)
coherence = lapply(sample_ids, function(s) calculate_spatial_coherence(s, mp_names, use_norm_scores = T))
names(coherence) = sample_ids
saveRDS(coherence, file.path(base_dir, "coherence.rds"))
          
# 1.2 calculate spatial complexity 
# recommended to be run in parallel on server
# to run locally:
win_sizes = c(5, 8, 11)
rand_num = 100
for (w in win_sizes) {
  all_scores = lapply(sample_ids, function(s) calculate_spatial_complexity(s, w, rand_num, mp_names)) # runs sequentially
  saveRDS(all_scores, paste0(base_dir, "/01_scores_w", w, ".rds"))
}

# load complexity scores of all samples
complexity = lapply(sample_ids, function(s) readRDS(paste0("complexity_scores_", s, ".rds")))
names(complexity) = sample_ids

# average across win_sizes
complexity_win_combined = lapply(sample_ids, FUN = function(s) {
  df1 = complexity[[s]][[1]]; df2 = complexity[[s]][[2]]; df3 = complexity[[s]][[3]]
  array_of_dfs = array(c(as.matrix(df1), as.matrix(df2), as.matrix(df3)), dim = c(nrow(df1), ncol(df1), 3))
  mean_scores = as.data.frame(apply(array_of_dfs, c(1, 2), function(x) mean(x, na.rm = TRUE))) # averaging the three win_sizes
  rownames(mean_scores) = rownames(df1); colnames(mean_scores) = colnames(df1)
  return(mean_scores)
})
names(complexity_win_combined) = sample_ids
saveRDS(complexity_win_combined, file.path(base_dir, "complexity_win_combined.rds"))
                              
# 1.3 plot coherence and complexity scores in tissue. first by program, then merged
out_dir = ""
for (s in sample_ids) {
  message("plotting scores in tissue of ", s)
  sstobj = readRDS(file = paste0("Objects/sstobj_", s, ".rds")) # takes a while bc. every object is loaded here...
  ratio = get_spatial_aspect_ratio(sstobj); pt_size = get_spatial_point_size(sstobj, scale_factor = 1.15) # spatial plot params
  co_plots = list(); gr_plots = list() # individual program scores
  for (mp in mp_names) {
    sstobj$coherence = coherence[[s]][, mp]
    co_plots[[mp]] = SpatialPlot(sstobj, pt.size.factor = pt_size, features = "coherence", image.scale = NULL) + ggtitle(mp) + theme(plot.title = element_text(hjust = 0.5), aspect.ratio = ratio)
    sstobj$complexity = complexity_win_combined[[s]][[mp]]
    gr_plots[[mp]] = SpatialPlot(sstobj, pt.size.factor = pt_size, features = "complexity", image.scale = NULL) + ggtitle(mp) + theme(plot.title = element_text(hjust = 0.5), aspect.ratio = ratio)
  }
  global_co_max <- max(sapply(co_plots, function(p) max(p$data$coherence, na.rm = TRUE)))
  co_plots = lapply(co_plots, function(p) p + scale_fill_gradientn(colors = Seurat:::SpatialColors(n = 100), limits = c(0, global_co_max)))
  ggsave(filename = paste0(s, "_a_spatial_coherence.png"), plot = wrap_plots(co_plots, ncol = 5), width = 20, height = 10, path = out_dir)
  global_gr_max <- max(sapply(gr_plots, function(p) max(p$data$complexity, na.rm = TRUE)))
  gr_plots = lapply(gr_plots, function(p) p + scale_fill_gradientn(colors = Seurat:::SpatialColors(n = 100), limits = c(0, global_gr_max)))
  ggsave(filename = paste0(s, "_b_spatial_complexity.png"), plot = wrap_plots(gr_plots, ncol = 5), width = 20, height = 10, path = out_dir)

  # 'merge' scores across programs by taking the maximal complexity score
  sstobj$coherence = apply(coherence[[s]], 1, function(x) max(x, na.rm = T))
  co = SpatialPlot(sstobj, pt.size.factor = pt_size, features = "coherence", image.scale = NULL) + ggtitle("max coherence") + theme(legend.position = "right", aspect.ratio = ratio)
  sstobj$complexity = apply(complexity_win_combined[[s]], 1, function(x) max(x, na.rm = T))
  cplx = SpatialPlot(sstobj, pt.size.factor = pt_size, features = "complexity", image.scale = NULL) + ggtitle("max complexity") + theme(aspect.ratio = ratio)
  ggsave(filename = paste0(s, "_c_spatial_scores_summary.png"), plot = (co|cplx), width = 20, height = 10, path = out_dir)
}
coherences_prog_max_combined = lapply(coherence, function(s_df) { apply(s_df, 1, max) })
saveRDS(coherences_prog_max_combined, file.path(base_dir, "coherence_max_combined.rds"))

complexity_prog_max_combined_wins = lapply(complexity, function(s) { lapply(c(1,2,3), function(i) { apply(s[[i]], 1, max, na.rm = T) }) }); names(complexity_prog_max_combined_wins) = sample_ids
complexity_prog_max_combined_wins = lapply(complexity_prog_max_combined_wins, function(s) { names(s) = win_sizes; s })
saveRDS(complexity_prog_max_combined_wins, file.path(base_dir, "complexity_prog_max_combined_wins.rds"))
complexity_prog_max_combined = lapply(complexity_win_combined, function(s_df) { apply(s_df, 1, max, na.rm = T) })
saveRDS(complexity_prog_max_combined, file.path(base_dir, "complexity_prog_max_combined.rds"))
complexity_prog_mean_combined = lapply(complexity_win_combined, function(s_df) { apply(s_df, 1, mean, na.rm = T) })
saveRDS(complexity_prog_mean_combined, file.path(base_dir, "complexity_prog_mean_combined.rds"))
                    
# 1.4 plot averages
# i) plot average coherence and gradients of programs across samples, corrected for the mean program abundance in a sample
p1 = plot_mean_score(coherence, "coherence", adjust_by_abund = F) + ggtitle("Coherence")
p2 = plot_mean_score(complexity_win_combined, "complexity") + ggtitle("Complexity")
ggsave("01a_mean_spatial_scores_by_MPs_coh-unadjusted.png", plot = p1 + p2, width = 7, height = 3, path = base_dir)
p = plot_mean_score(complexity, "complexity", separate_wins = T) + ggtitle("Complexity (by windows)")
ggsave("01b_mean_spatial_complexity_by-wins.png", plot = p, width = 7, height = 3, path = base_dir)

# split by subtype
subtype = metadata$Subtype[match(sample_ids, metadata$`Study ID`)]; subtype = factor(subtype, levels = c("ZFTA", "SE", "YAP1"))
p = plot_mean_score(coherence, score_name = "coherence", split.by = subtype, adjust_by_abund = T)
ggsave("01c_mean_coherence_by_MPs_split-by-subtype.png", plot = p, width = 10, height = 3, path = base_dir)
p = plot_mean_score(complexity_win_combined, score_name = "complexity", split.by = subtype)
ggsave("01d_mean_complexity_by_MPs_split-by-subtype.png", plot = p, width = 10, height = 3, path = base_dir)

# split by hyp_class
hyp_norm_score = lapply(all_spots_programs_comp_spot_norm, function(s) s[, "Hypoxia"])
hyp_threshold = 0.3 # defining hypoxic presence as > 30% hypoxia expression in a spot (based on spot-norm scores)
hyp_pct_norm_scores = sapply(hyp_norm_score, function(s) sum(s > hyp_threshold) / length(s))
hyp_pct_norm_scores[order(hyp_pct_norm_scores)]
quantile(hyp_pct_norm_scores)
# divide cohort into 3 classes based on quantiles (< 40%, 40-60%, > 60%)
hyp_class = ifelse(hyp_pct_norm_scores < quantile(hyp_pct_norm_scores, seq(0, 1, 0.2))[3], "hyp_low", 
                    ifelse(hyp_pct_norm_scores >= quantile(hyp_pct_norm_scores, seq(0, 1, 0.2))[4], "hyp_high", "intermediate"))
hyp_class = factor(hyp_class, levels = c("hyp_high", "hyp_low", "intermediate"))
saveRDS(hyp_class, "Data/hyp_class.rds")
saveRDS(hyp_pct_norm_scores, "Data/hyp_score.rds")
p = plot_mean_score(coherence, score_name = "coherence", split.by = hyp_class, adjust_by_abund = T)
ggsave("01e_mean_coherence_by_MPs_split-by-hyp-status.png", plot = p, width = 10, height = 3, path = base_dir)
p = plot_mean_score(complexity_win_combined, score_name = "complexity", split.by = hyp_class)
ggsave("01f_mean_complexity_by_MPs_split-by-hyp-status.png", plot = p, width = 10, height = 3, path = base_dir)

# ii) plot coherence and complexity of each program per sample, all together and split by subtype
data_coh = program_coherence_per_sample_data()
data_cplx = program_complexity_per_sample_data()
p = plot_score_per_sample(data_coh, ylab = "mean coherence / abundance") + labs(fill = "program", title = "Spatial coherence per sample (meaned across spots and win_sizes)")
ggsave("02a_mean_MP_coherence_by_sample.png", plot = p, width = 10, height = 5, path = base_dir)
p = plot_score_per_sample(data_cplx, ylab = "spatial complexity") + labs(title = "Spatial complexity per sample (meaned across spots and win_sizes)")
ggsave("02b_mean_MP_complexity_by_sample.png", plot = p, width = 10, height = 5, path = base_dir)

p = plot_score_per_sample(data_coh, split.by = subtype, ylab = "mean coherence / abundance")
ggsave("02c_mean_MP_coherence_by_sample_split-by-subtype.png", plot = p, width = 10, height = 5, path = base_dir)
p = plot_score_per_sample(data_cplx, split.by = subtype, ylab = "spatial complexity", anno = hyp_score, anno_name = "hyp score")
ggsave("02d_mean_MP_complexity_by_sample_split-by-subtype.png", plot = p, width = 10, height = 5, path = base_dir)

p = plot_score_per_sample(data_coh, split.by = hyp_class, ylab = "mean coherence / abundance")
ggsave("02e_mean_MP_coherence_by_sample_split-by-hyp-status.png", plot = p, width = 10, height = 5, path = base_dir)
p = plot_score_per_sample(data_cplx, split.by = hyp_class, ylab = "spatial complexity", anno = hyp_score, anno_name = "hyp score")
ggsave("02f_mean_MP_complexity_by_sample_split-by-hyp-status.png", plot = p, width = 10, height = 5, path = base_dir)

# ordered by hyp score
hyp_order = names(hyp_score)[order(hyp_score, decreasing = T)]
p = plot_score_per_sample(data_cplx, ylab = "spatial complexity", anno = hyp_score, anno_name = "hyp score", order_manual = hyp_order)
ggsave("03a_mean_MP_complexity_by_sample_ordered-by-hyp-score.png", plot = p, width = 10, height = 5, path = base_dir)
p = plot_score_per_sample(data_cplx, split.by = hyp_class, ylab = "spatial complexity", anno = hyp_score, anno_name = "hyp score", order_manual = hyp_order)
ggsave("03b_mean_MP_complexity_by_sample_ordered-by-hyp-score-split-by-hyp-class.png", plot = p, width = 10, height = 5, path = base_dir)

# iii) plot mean complexity per program as a direct comparison of two groups
p = plot_mean_score_comparative(complexity_win_combined, split.by = hyp_class, pair = c("hyp_high", "hyp_low"), order_by_difference = T, cols = c(hyp_high = "black", hyp_low = "#e3d536"))
ggsave("04a_mean_MP_complexity_by_sample_comparing-hyp-class_ordered-by-diff.png", plot = p, width = 4.5, height = 3, path = base_dir)
p = plot_mean_score_comparative(complexity_win_combined, split.by = hyp_class, pair = c("hyp_high", "hyp_low"), order_by_difference = F, cols = c(hyp_high = "black", hyp_low = "#e3d536"))
ggsave("04b_mean_MP_complexity_by_sample_comparing-hyp-class_ordered-by-hyp-high-cplx.png", plot = p, width = 4.5, height = 3, path = base_dir)

# iv) statistics on observed differences between hyp classes
data_cplx$hyp_class = hyp_class[match(data_cplx$sample, names(hyp_class))]; data_cplx$hyp_score = hyp_score[match(data_cplx$sample, names(hyp_score))]
# condensed data (one row per sample)
data = data_cplx %>% distinct(sample, .keep_all = T); data = data[, c(3,4,6,7)]
# Mann-Whitney U test: hyp-class (present or absent) vs. average sample complexity score
wilcox.test(average_score ~ hyp_class, data = data[data$hyp_class %in% c("hyp_low", "hyp_high"), ])
cols = c(hyp_high = "black", intermediate = "#4e81a3", hyp_low = "#e3d536")
data$hyp_class = factor(data$hyp_class, levels = unique(data$hyp_class))
ggplot(data, aes(x = average_score, y = hyp_score*100, color = hyp_class)) + geom_point(size = 2) + # geom_smooth(method = "lm", se = F, color = "grey50") +
  scale_color_manual(values = cols) + labs(x = "Average complexity", y = "Percent hypoxic spots") + theme_classic()
ggsave("05_hypoxic-expr_vs_avg-cplx_statistic.png", width = 4, height = 2, path = base_dir)

# find every MPs complexity vs. hyp-score correlation
mp_hyp_pvals = c()
mp_hyp_cors = sapply(mp_names, function(mp) {
  data = data_cplx[data_cplx$metaprogram == mp, ]
  mp_hyp_pvals[mp] <<- summary(lm(data$score ~ data$hyp_score))$coefficients[2,4]
  cor(data$score, data$hyp_score)
}); names(mp_hyp_cors) = mp_names

# find every MPs complexity vs. hyp-class significance
mp_hyp_pvals = sapply(mp_names, function(mp) {
  data = data_cplx[data_cplx$metaprogram == mp & data_cplx$hyp_class %in% c("hyp_low", "hyp_high"), ]
  wilcox.test(score ~ hyp_class, data = data)$p.value
}); names(mp_hyp_pvals) = mp_names
mp_hyp_pvals[order(mp_hyp_pvals)]



# -------------------------------- PART 2: define structural zones ----------------------------                    
out_dir = ""
                             
# 2.1 define structural zones based on max-combined spatial complexity scores per sample
quantiles = as.numeric(quantile(na.omit(unlist(complexity_prog_max_combined)), probs = (seq(0, 1, 0.05))))
# assign labels based on quantiles 
complexity_zones = sapply(sample_ids, function(s) {ifelse(complexity_prog_max_combined[[s]] <= quantiles[6], "disorganised_high",         # 0-25%
                                                          ifelse(complexity_prog_max_combined[[s]] <= quantiles[11], "disorganised_low",         # 26-50%
                                                                 ifelse(complexity_prog_max_combined[[s]] <= quantiles[16], "structured_low",    # 51-75%
                                                                        "structured_high")))})                                               # 76-100%
smoothing_win_size = 4
# smoothes by majority of neighbouring spots: 'clear' structured or disorganised zones: > 60% of spots classified accordingly, mixed zones: > 40% of both structured and disorganised zones
# smoothing approach adapted from Greenwald et al. 24
all_zones_smoothened = lapply(sample_ids, function(s) {
  message("smoothing of zones in ", s, " using 'smoothing_win_size' = ", smoothing_win_size)
  neighbors_table = neighbors_table_funcV2(all_spots_positions[[s]], all_spots_programs_comp_norm[[s]]); neighbors_table[neighbors_table== "NaN"] = NA
  spots_smo = complexity_zones[[s]]
  for (spot in names(spots_smo)) {
    win_spots = c(spot); for (j in 1:smoothing_win_size) win_spots = unique(c(win_spots, unique(na.omit(as.character(neighbors_table[win_spots, ])))))
    win_zones = spots_smo[names(spots_smo) %in% win_spots]
    spots_smo[spot] = ifelse(length(win_zones[win_zones %in% c("disorganised_high", "disorganised_low")])/length(win_zones) >= 0.6, 
                             ifelse(length(win_zones[win_zones == "disorganised_high"]) > length(win_zones[win_zones == "disorganised_low"]), "disorganised_high", "disorganised_low"), 
                             ifelse(length(win_zones[win_zones %in% c("structured_low", "structured_high")])/length(win_zones) >= 0.6, 
                                    ifelse(length(win_zones[win_zones == "structured_high"]) > length(win_zones[win_zones == "structured_low"]), "structured_high", "structured_low"),
                                    "intermediate")) # case when none of the two major zones (disorg. or struct) reached at least 60 % agreement
  }
  return(spots_smo)
})
complexity_zones = all_zones_smoothened; names(complexity_zones) = sample_ids
saveRDS(complexity_zones, file.path(base_dir, "complexity_zones.rds"))

# 2.2 plot structural zones based on complexity or coherence (together with discrete malignancy) in every sample
zone_type = "complexity"       # choose 'coherence' or 'complexity'-based zones to run this with
malignancy_type = "raw"        # choose 'norm' or 'raw' for normalised or non-normalised maligancy score
all_zones = if (zone_type == "coherence") all_coherence_zones else complexity_zones
all_mscores = if (malignancy_type == "norm") all_norm_malignancy_scores else all_raw_malignancy_scores
mscore_threshold = if (malignancy_type == "norm") 0.5 else 1 # threshold to decide if spot is considered malignant or not
zone.colors = c(structured_high = "#a82203", structured_low = "#cf5e4e", disorganised_high = "#003967", disorganised_low = "#208cc0", non_malignant = "#11A579", intermediate ="#f1af3a")

spot_malignancy = lapply(sample_ids, function(s) ifelse(all_mscores[[s]] >= mscore_threshold, "mal", "non_mal")); names(spot_malignancy) = sample_ids
for (s in sample_ids) {
  message("plotting structural zone in ", s)
  spots_filt = all_spots_positions[[s]][all_spots_positions[[s]]$V1 %in% rownames(all_spots_programs_comp_norm[[s]]), ]
  spots_filt$zone = all_zones[[s]][spots_filt$V1]
  spots_filt$malignancy = spot_malignancy[[s]][spots_filt$V1]
  if (s %in% c("399", "415A")) { spots_filt$tmp = spots_filt$V3; spots_filt$V3 = -spots_filt$V4; spots_filt$V4 = spots_filt$tmp } # adjust for different orientation in Visium v1
  p = ggplot(spots_filt, aes(x=V4, y=-V3, color = zone, alpha = malignancy)) + geom_point(shape = 16, size=2.3) + labs(title=paste(s, "spatial", zone_type, "zones")) + 
    scale_alpha_manual(values = c("non_mal" = 0.3, "mal" = 1)) + scale_color_manual(values = zone.colors, name = "zone") + 
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(s, "_structural_zones_", zone_type, ".png"), plot = p, width = 7, height = 6, path = out_dir)
}


# ---------------------- PART 3: investigate MP abundance in structural zones --------------------
out_dir = ""

# 3.1 use raw scores and a hard threshold to consider only spots with significant expression of a signature
mp_zone_abund = lapply(sample_ids, function(s) { # simpler variant without significance testing
  mp_scores = all_spots_programs_comp[[s]]
  zone_abund = sapply(mp_names, function(mp) {
    tapply(mp_scores[, mp], as.character(all_zones[[s]][rownames(mp_scores)]), function(zone_spot_scores) sum(zone_spot_scores > 0.6))
  })
  if (is.null(dim(zone_abund))) { zone_abund = data.frame(t(zone_abund)); rownames(zone_abund) = unique(all_zones[[s]][rownames(mp_scores)]); colnames(zone_abund) = mp_names }
  for (row_name in c("structured_high", "structured_low", "intermediate", "disorganised_low", "disorganised_high")) { # add missing zones to allow later comparisons!
    if (!row_name %in% rownames(zone_abund)) { zone_abund = rbind(zone_abund, rep(0, ncol(zone_abund))); rownames(zone_abund)[nrow(zone_abund)] = row_name }
  }
  zone_abund = zone_abund[c("structured_high", "structured_low", "intermediate", "disorganised_low", "disorganised_high"), ] # bring in consistent order and remove 'intermediate'
  return(zone_abund)
}); names(mp_zone_abund) = sample_ids

# 3.2 plot relative abundances of MPs in zones as barplots
plots = plot_mp_zone_abundance(mp_zone_abund, zone.colors, ylim = c(0.8, 1.2))
ggsave("01_mean_abundance_in_zones.png", plot = plots[[1]], width = 5, height = 3, path = out_dir)

