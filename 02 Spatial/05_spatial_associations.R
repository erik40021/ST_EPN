library(stringr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(pheatmap)
library(Seurat)
library(circlize)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(patchwork)

source("utils/spatial_score_utils.R")

base_dir = ""

metadata = read_excel("masterlist_EPN-ST.xlsx", skip = 1); metadata = metadata[grepl("2", metadata$`Analyses*`), ]
sample_ids = metadata$`Study ID`
metaprograms = readxl::read_excel("spatial_NMF_metaprograms_main.xlsx", sheet = "MP genes (final)")
mp_names = colnames(metaprograms)


# ------------------------------ PART 1: Pairwise program associations across distances ------------------------------
out_dir = ""

# 1.1 run association calculation
for (s in sample_ids) { 
  calculate_association(s, mp_names, pairs, pairs_names, all_spots_positions[[s]], all_spots_programs_comp[[s]], out_dir, max_win_size) 
}

# 1.2 extract scores from all samples
association_results <- lapply(sample_ids, function(s) readRDS(paste0(out_dir, "/", s, "_obs_association_w", max_win_size, ".rds"))); names(association_results) <- sample_ids
association = lapply(association_results, function(s) s[[1]]); names(association) <- sample_ids
saveRDS(association, file.path(base_dir, "association.rds"))


# ------------------------------ PART 2: Plot association scores ------------------------------
out_dir = ""

# 2.1 plot all associations as heatmaps
# a) individually by sample
for (s in sample_ids) {
  mtx = as.matrix(association[[s]])
  png(paste0(out_dir, "/heatmap_", s, "_association.png"), width = 2200, height = 2400, res = 300)
  print(plot_association_heatmap(mtx, order_mode = "association_type", title = paste0("Association ", s)))
  dev.off()
}

# b) averaged across all samples
avg_association <- Reduce(function(x, y) Map(`+`, x, y), association); avg_association <- as.data.frame(lapply(avg_association, function(x) x / length(sample_ids)))
colnames(avg_association) = 0:max_win_size; rownames(avg_association) = pairs_names
saveRDS(avg_association, file.path(base_dir, "avg_association.rds"))#; saveRDS(avg_pvals, file.path(base_dir, "avg_pvals.rds"))
png(file.path(base_dir, "01a_heatmap_avg_association.png"), width = 2200, height = 2400, res = 300)
plot_association_heatmap(avg_association, order_mode = "association_type", title = paste0("Avg. association across all samples"))
dev.off()
png(file.path(base_dir, "01a_heatmap_avg_association_fixed-scale.png"), width = 2000, height = 2400, res = 300)
plot_association_heatmap(avg_association, order_mode = "association_type", title = paste0("Avg. association across all samples"), fixed_scale = T)
dev.off()

# c) averaged across subtypes
heatmaps = lapply(unique(metadata$Subtype), function(subtype) {
  samples = metadata$`Study ID`[metadata$Subtype == subtype]
  subt_association <- Reduce(function(x, y) Map(`+`, x, y), association[samples]); subt_association <- as.data.frame(lapply(subt_association, function(x) x / length(samples)))
  rownames(subt_association) = pairs_names
  heatmap = plot_association_heatmap(subt_association, order_mode = "association_type", title = subtype, silent = T, fixed_scale = T)
  return(heatmap$gtable)
})
ggsave("01b_heatmap_avg_association_by_subtypes_fixed-scale.png", plot = grid.arrange(grobs = heatmaps, nrow = 1), width = 15, height = 8, path = base_dir)

# d) averaged across hyp_classes
heatmaps = lapply(unique(hyp_class), function(hyp) {
  samples = names(hyp_class)[hyp_class == hyp]
  subt_association <- Reduce(function(x, y) Map(`+`, x, y), association[samples]); subt_association <- as.data.frame(lapply(subt_association, function(x) x / length(samples)))
  rownames(subt_association) = pairs_names
  heatmap = plot_association_heatmap(subt_association, order_mode = "association_type", title = hyp, silent = T, fixed_scale = T)
  return(heatmap$gtable)
})
ggsave("01c_heatmap_avg_association_by_hyp-status_fixed-scal.png", plot = grid.arrange(grobs = heatmaps, nrow = 1), width = 15, height = 8, path = base_dir)


# 2.2 plot proximities as 'target' plots
# visualisation was inspired by Greenwald et al. 24
# i. all pairs individually
out_dir = ""
for (pair in rownames(avg_association)) { 
  data <- data.frame(score = unlist(avg_association[pair, 16:1]), radius = rev(seq(from = 0.5, to = 7.5, length.out = 16)))
  png(paste0(out_dir, "/target_all_proximities_", pair, ".png"), width = 2000, height = 1500, res = 300)
  target_plot(data, radii = "all")
  dev.off()
}
# ii. stacked pairs (39 per plot) for overview
avg_association_ordered = avg_association[order(rowMeans(avg_association, na.rm = T), decreasing = T), ]
for (i in 1:ceiling(nrow(avg_association_ordered)/39)) {
  png(paste0(base_dir, "/03_stacked_overview_targets_", (i-1)*39+1, "-", min(nrow(avg_association_ordered), (i*39)), ".png"), width = 5000, height = 2500, res = 400)
  par(mfrow=c(5, 8), mar=c(0, 0, 1.5, 0.2))  # Setting up a grid of 5x8 plots with minimal margins
  for (pair in rownames(avg_association_ordered)[((i-1)*39+1):min(nrow(avg_association_ordered), (i*39))]) {
    data <- data.frame(score = unlist(avg_association_ordered[pair, 16:1]), radius = rev(seq(from = 0.5, to = 7.5, length.out = 16)))
    target_plot(data, radii = "all", add_legend = F)
    split_title = strsplit(pair, " x ")[[1]]
    title(main=paste(split_title[1], "\n", split_title[2], sep=""), cex.main = 1, line=-0.5) # Re-add ")" to split title
  }
  plot(0, 0, type="n", xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="", xaxt='n', yaxt='n', bty='n') # plot also the legend
  color_bar <- matrix(1:100, ncol=1)
  rasterImage(as.raster(matrix(colorRampPalette(brewer.pal(11, "RdBu"))(100)[color_bar], ncol=1)), xleft=0.4, ybottom=0, xright=0.6, ytop=1, interpolate=TRUE)
  tick_marks <- c(1, 0.5, 0, -0.5, -1); tick_positions <- (tick_marks + 1) / 2  # Normalize between [0, 1]
  axis(2, at=tick_positions, labels = F, cex.axis=0.7, pos = c(0.6, 0)); axis(4, at=tick_positions, labels = F, cex.axis=0.7, pos = c(0.4, 0))
  text(0.6, tick_positions, labels=tick_marks, cex=0.9, pos=4); text(0.5, 1.1, labels="association", xpd=TRUE, cex = 1.2)
  dev.off()
}


# ------------------------------ PART 3: Define functional zones ------------------------------
out_dir = ""
mscore_threshold = 1
spot_malignancy = lapply(sample_ids, function(s) ifelse(all_raw_malignancy_scores[[s]] >= mscore_threshold, "mal", "non_mal")); names(spot_malignancy) = sample_ids

# define functional zones based on higher-order model
zones = list(hyp = c("Hypoxia", "Stress"), 'hyp-adj' = c("Inflammatory"), 'angio-adj' = c("Myeloid", "Vasc", "Fibroblast"),
             Tdiff = c("Astroglial", "Cilia", "Interferon-re"), Tundiff = c("ZFTA-fus", "NSC"))
library(ggforce)
library(scatterpie)
all_spots_fzone_score_norm = lapply(all_spots_programs_comp, function(s) {
  as.data.frame(t(apply(s, 1, function(spot) { 
    scores = sapply(zones, function(zone) mean(as.numeric(unlist(spot[zone]))))
    scores[scores < 0] = 0
    scores / sum(scores) # normalise by spot
  })))
})
saveRDS(all_spots_fzone_score_norm, file.path(base_dir, "functional_zones_norm.rds"))

# plot functional zones in tissues
for (s in sample_ids) {
  message("plotting functional zone in ", s)
  spots_filt = all_spots_positions[[s]][all_spots_positions[[s]]$V1 %in% rownames(all_spots_programs_comp[[s]]), ]
  spots_filt = cbind(spots_filt, all_spots_fzone_score_norm[[s]][spots_filt$V1, ])
  spots_filt$malignancy = spot_malignancy[[s]][spots_filt$V1]
  if (s %in% c("399", "415A")) { spots_filt$tmp = spots_filt$V3; spots_filt$V3 = -spots_filt$V4; spots_filt$V4 = spots_filt$tmp } # adjust for different orientation in Visium v1
  p = spatial_feature_pies(spots_filt, pie_size = 0.85, fzone_cols = fzone_cols)
  ggsave(paste0(s, "_functional_zones_scatterpie.png"), plot = p, width = 7, height = 6, path = file.path(out_dir, "02 per sample (model-based)"))
}


# ------------------------------ PART 4: Composition of spatial zones ------------------------------
out_dir = ""
# use higher-resolution single-cell derived programs and determine their abundance in structural and functional zones to precisely locate
# them spatially (despite their undetectability in spot-resolution data)
sig_cols = as.character(my_cols)
szone_order = c("structured_high", "structured_low", "intermediate", "disorganised_low", "disorganised_high")
fzone_order = c("hyp", "hyp-adj", "angio-adj", "Tdiff", "Tundiff")
complexity_zones = readRDS("complexity_zones.rds")
functional_zones = readRDS(file.path(base_dir, "functional_zones_norm.rds"))

# 4.1 as a base reference, plot the average spot composition regarding s and fzones
# a. average for structural zones
all_szone_comps = readRDS("complexity_rel_sample_zone_abund.rds")
comps_by_subtype = sapply(unique(subtype), function(subt) rowMeans(all_szone_comps[, subtype == subt, drop = FALSE])); colnames(comps_by_subtype) = unique(subtype)
comps_by_hypclass = sapply(unique(hyp_class), function(hyp) rowMeans(all_szone_comps[, hyp_class == hyp, drop = FALSE])); colnames(comps_by_hypclass) = unique(hyp_class)
data = rbind(all = rowMeans(all_szone_comps), t(comps_by_subtype), t(comps_by_hypclass))
data = reshape2::melt(data); colnames(data) = c("group", "zone", "value")
data$zone = factor(data$zone, levels = c("structured_high", "structured_low", "disorganised_high", "disorganised_low", "intermediate"))
data[data$zone %in% c("intermediate", "disorganised_low", "disorganised_high"), ]$value = -data[data$zone %in% c("intermediate", "disorganised_low", "disorganised_high"), ]$value
p = ggplot(data, aes(x = group, y = value, fill = zone)) + geom_bar(stat = "identity", position = "stack") + geom_hline(yintercept=0, color="black") +
  labs(x = NULL, y = "relative abundance") + scale_fill_manual(values = szone_cols) + 
  scale_y_continuous(limits = c(-1,1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c(1, 0.5, 0, 0.5, 1)) + theme_minimal() + #ylim(-ylim[1], ylim[2]) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid.major.x = element_blank())
ggsave("avg_composition_of_szones_across_groups.png", plot = p, width = 7, height = 6, path = out_dir)
# b. average for functional zones
all_fzone_comps = sapply(functional_zones, function(s) colMeans(s, na.rm = T))
comps_by_subtype = sapply(unique(subtype), function(subt) rowMeans(all_fzone_comps[, subtype == subt, drop = FALSE])); colnames(comps_by_subtype) = unique(subtype)
comps_by_hypclass = sapply(unique(hyp_class), function(hyp) rowMeans(all_fzone_comps[, hyp_class == hyp, drop = FALSE])); colnames(comps_by_hypclass) = unique(hyp_class)
data = rbind(all = rowMeans(all_fzone_comps), t(comps_by_subtype), t(comps_by_hypclass))
data = reshape2::melt(data); colnames(data) = c("group", "zone", "value")
data$zone = factor(data$zone, levels = c("hyp", "hyp-adj", "Tundiff", "Tdiff", "angio-adj"))
data[data$zone %in% c("angio-adj", "Tdiff", "Tundiff"), ]$value = -data[data$zone %in% c("angio-adj", "Tdiff", "Tundiff"), ]$value
p = ggplot(data, aes(x = group, y = value, fill = zone)) + geom_bar(stat = "identity", position = "stack") + geom_hline(yintercept=0, color="black") +
  labs(x = NULL, y = "relative abundance") + scale_fill_manual(values = fzone_cols) + 
  scale_y_continuous(limits = c(-1,1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c(1, 0.5, 0, 0.5, 1)) + theme_minimal() + #ylim(-ylim[1], ylim[2]) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid.major.x = element_blank())
ggsave("avg_composition_of_fzones_across_groups.png", plot = p, width = 7, height = 6, path = out_dir)


# 4.2 investigate specific reference myeloid programs (from Miller et al. 2025)
out_dir = ""
myeloid_sigs = as.list(openxlsx::read.xlsx("Miller 2025 cNMF genes.xlsx", sheet = "Myeloid Top 100 genes per prog."))
min_score_threshold = 0.3 # threshold to decide if a signatures is 'significantly expressed' in a spot
myeloid_scores = get_sig_scores(myeloid_sigs, plot_sigs = T, out_dir = out_dir)[[1]] # only need raw scores here
saveRDS(myeloid_scores, file.path(out_dir, "myeloid_sig_scores.rds"))
myeloid_scores = readRDS(file.path(out_dir, "myeloid_sig_scores.rds"))
szone_abund = get_sig_szone_abund(myeloid_scores, complexity_zones, min_score_threshold)
fzone_abund = get_sig_fzone_abund_by_threshold(myeloid_scores, functional_zones, min_score_threshold)
plot_all_sig_zone_abundances(szone_abund, fzone_abund, out_dir)


# 4.3 investigate the specific single-nucleus immune programs from this study
out_dir = ""
immune_sigs = as.list(openxlsx::read.xlsx("immune_NMF_metaprograms_main.xlsx", sheet = "MP genes (final)"))
min_score_threshold = 0.3 # threshold to decide if a signatures is 'significantly expressed' in a spot
immune_scores = get_sig_scores(immune_sigs, plot_sigs = T, out_dir = out_dir)[[1]] # only need raw scores here
saveRDS(immune_scores, file.path(out_dir, "immune_sig_scores.rds"))
immune_scores = readRDS(file.path(out_dir, "immune_sig_scores.rds"))
szone_abund = get_sig_szone_abund(immune_scores, complexity_zones, min_score_threshold)
fzone_abund = get_sig_fzone_abund_by_threshold(immune_scores, functional_zones, min_score_threshold)
plot_all_sig_zone_abundances(szone_abund, fzone_abund, out_dir)



