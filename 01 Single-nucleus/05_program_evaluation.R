library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

source("utils/seurat_utils.R")
source("utils/r_utils.R")

malig_sobj = readRDS(file = "Objects/ST_EPN_malig.rds")
mps = openxlsx::read.xlsx("NMF_metaprograms_main.xlsx", sheet = "MP genes (final)")
out_dir = "MP evaluation"; if (!dir.exists(out_dir)) dir.create(out_dir)


# ------ PART 1: evaluate metaprograms regarding favourable and unfavourable genes ------

# favourability of single genes was downloaded from R2-Browser with Korshunov et al. bulk ZFTA transcriptomes (via PAGE enrichment in fav. vs. unfav. samples)
gene_favs = read.table(file.path(out_dir, "/all_genes_enrichment.txt"))[, 2:4]; colnames(gene_favs) = c("gene", "pval_adj", "log2fc")
gene_favs$gene = gsub("ORF", "orf", gene_favs$gene)
data <- mps %>% pivot_longer(cols = everything(), names_to = "mp", values_to = "gene") %>% filter(!is.na(gene))
data = cbind(gene_favs, metaprogram = data$mp[match(gene_favs$gene, data$gene)])
ggplot(data, aes(x = log2fc, y = -log10(pval_adj), color = metaprogram)) + geom_point(size = 3, alpha = 0.9, shape = 16) + 
  scale_color_manual(values = mp_cols) + # geom_smooth(method = "lm", se = FALSE, aes(group = metaprogram), linetype = "dashed", size = 0.5, alpha = 0.5) +  # Linear trendlines
  theme_minimal() + labs(x = "log2 fold change", y = "-log10(p-value)") + # coord_cartesian(xlim = c(-2, 2), ylim = c(0, 6)) + 
  theme(legend.position = "right", panel.grid.major.x = element_line(size = 0.2), panel.grid.major.y = element_line(size = 0.2), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) # + geom_vline(xintercept = 0, linetype = "dashed", color = "black")
ggsave("volcano_favourability_by_mp-genes.png", width = 8, height = 4, path = file.path(out_dir, "Favourability"))
# re-structure by MP, so that every MP's gene is counted for the average, even if it is duplicated!
data_by_mp <- mps %>% pivot_longer(cols = everything(), names_to = "metaprogram", values_to = "gene") %>% filter(!is.na(gene))
data_by_mp = cbind(data_by_mp, gene_favs[match(data_by_mp$gene, gene_favs$gene), 2:3])
avg_fc <- data_by_mp %>% group_by(metaprogram) %>% dplyr::summarise(mean_lfc = mean(log2fc, na.rm = TRUE))
ggplot(avg_fc, aes(x = mean_lfc, y = reorder(metaprogram, mean_lfc), fill = metaprogram)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = mp_cols) + # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
  geom_vline(xintercept = 0, color = "white") + labs(x = NULL, y = NULL, fill = "average\nlog2fc") + # coord_cartesian(xlim = c(-2, 2)) +
  theme_minimal() + theme(legend.position = "right", panel.grid.major.x = element_line(size = 0.2), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
ggsave("avg_favourability_of_MPs.png", width = 6, height = 1.5, path = file.path(out_dir, "Favourability"))

# plot all MP genes separatly, with gene labels and trendline
library(ggrepel)
plots = lapply(unique(data_by_mp$metaprogram), function(mp) {
  data = data_by_mp[data_by_mp$metaprogram == mp, ]
  data$label <- ifelse(data$log2fc > 0.5 & data$pval_adj < 0.05, data$gene, NA)
  max_x = max(abs(data$log2fc), na.rm = T)
  ggplot(data, aes(x = log2fc, y = -log10(pval_adj), color = metaprogram)) + geom_point(size = 3, alpha = 0.9, shape = 16) + 
    scale_color_manual(values = mp_cols) + geom_text_repel(aes(label = label), vjust = -0.5, size = 3, color = "black") + theme_minimal() + 
    labs(x = "log2 fold change", y = "-log10(p-value)", title = mp) + coord_cartesian(xlim = c(-max_x, max_x)) +
    theme(legend.position = "none", panel.grid.major.x = element_line(size = 0.2), panel.grid.major.y = element_line(size = 0.2), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) # + geom_vline(xintercept = 0, linetype = "dashed", color = "black")
})
ggsave("volcanos_favourability_by_mp_genes_separated-by-mp.png", plot = wrap_plots(plots, ncol = 4), width = 20, height = 10, path = file.path(out_dir, "Favourability"))
                   
mps_fav = c("Cilia", "Astroglial", "Inteferon-re", "Cilia-reg"); mps_unfav = c("ECM-remod", "ZFTA-fus", "Stress", "Hypoxia", "Cell-cycle", "Respiration")

# do statistics on the average enrichment per MP
shapiro.test(data_by_mp[data_by_mp$metaprogram == "ECM-remod", "log2fc"]) # p < 0.05, so the enrichment is not normally distributed for ECM-remod
qqnorm(data_by_mp[data_by_mp$metaprogram == "ECM-remod", "log2fc"])
shapiro.test(data_by_mp[, "log2fc"]) # p < 0.05, so enrichment in general is not normally distributed
qqnorm(data_by_mp[, "log2fc"])

# choose non-parametric test (data is not normally distributed)
wilcox.test(data_by_mp[data_by_mp$metaprogram == "ECM-remod", "log2fc"], mu = 0)  # Test against no enrichment (0), tells how different the data is from 0
all_pvals = sapply(unique(data_by_mp$metaprogram), function(mp) wilcox.test(data_by_mp[data_by_mp$metaprogram == mp, "log2fc"], mu = 0)$p.value)
x = all_pvals[order(all_pvals, decreasing = T)]; x[x < 0.05] # show only significant


                   
# ------ PART 2: evaluate here the MP-composition across samples, subtypes and subclasses ------
                   
# 1.1 composition by SAMPLE
# a. ordered by subclass
mp_cells_df = malig_sobj@meta.data[, c("orig.ident", "metaprogram", "subtype")]; colnames(mp_cells_df) = c("sample", "MP", "subtype")
mp_cells_df = mp_cells_df %>% group_by(sample, MP, subtype) %>% dplyr::summarise(cell_count = n(), .groups = "drop") %>% group_by(sample, subtype) %>% dplyr::mutate(pct = cell_count / sum(cell_count)) %>% ungroup()
low_cell_count_samples = names(table(malig_sobj$orig.ident))[table(malig_sobj$orig.ident) < 100] # exclude samples with too few malignant cells (unreliable)
mp_cells_df = mp_cells_df[!mp_cells_df$sample %in% low_cell_count_samples, ]
mp_cells_df$subclass = malig_sobj$subclass[match(mp_cells_df$sample, malig_sobj$orig.ident)]
mp_cells_df$subtype = factor(mp_cells_df$subtype, levels = c("ZFTA", "PLAGL1", "YAP1", "SE", "NOS")) # custom order of subtypes
mp_cells_df <- mp_cells_df %>% # group_by(subtype) %>% mutate(subclass = factor(subclass, levels = mp_cells_df %>% group_by(subtype, subclass) %>% summarise(size = n(), .groups = "drop") %>% arrange(subtype, desc(size)) %>% pull(subclass))) %>% ungroup() %>%
  group_by(subtype, subclass) %>% dplyr::mutate(sample = { 
    wide_data <- tidyr::pivot_wider(cur_data(), names_from = sample, values_from = pct); wide_data[is.na(wide_data)] <- 0
    if (nrow(wide_data) > 0 && ncol(wide_data) > 2) { mat <- as.matrix(wide_data[, -1, drop = FALSE]); if (ncol(mat) > 1) { hc <- hclust(dist(t(mat))); factor(sample, levels = hc$labels[hc$order]) 
      } else factor(sample, levels = unique(sample)) } else factor(sample, levels = unique(sample)) }) %>% ungroup()
ggplot(data = mp_cells_df, aes(x = sample, y = pct, fill = MP)) + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = mp_cols, name = "MP") +
  ggnewscale::new_scale_fill() + geom_tile(aes(x = sample, y = 1.04, fill = subtype), width = 1, height = 0.03, inherit.aes = FALSE) + scale_fill_manual(values = subtype_cols, name = "Subtype") +
  ggnewscale::new_scale_fill() + geom_tile(aes(x = sample, y = 1.08, fill = subclass), width = 1, height = 0.03, inherit.aes = FALSE) + scale_fill_manual(values = subclass_cols, name = "Subclass") + 
  labs(x = "", y = "proportion of cells") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), legend.key.size = unit(4, "mm"))
ggsave("01a_barplot_cell-proportions_by_sample_subclass-hclust-ordered.png", width = 10, height = 6, path = out_dir)

# b. only show significantly favourable vs. non-favourable programs
fav_mps_df = mp_cells_df[!mp_cells_df$MP %in% c("NPC", "Wounded", "Respiration"), ]
mps_fav = c("Cilia", "Astroglial", "Interferon-re", "Cilia-reg"); mps_unfav = c("ECM-remod", "ZFTA-fus", "Stress", "Hypoxia", "Cell-cycle")
fav_mps_df = fav_mps_df %>% group_by(sample) %>% dplyr::mutate(pct = cell_count / sum(cell_count)) %>% ungroup() # <- calculate percentage again based on only the selected MPs (makes no difference actually to ->) # mutate(pct = pct / sum(pct)) %>% ungroup()
fav_mps_df[fav_mps_df$MP %in% mps_fav, ]$pct = -fav_mps_df[fav_mps_df$MP %in% mps_fav, ]$pct
undiff_sums = sapply(unique(fav_mps_df$sample), function(s) sum(fav_mps_df[!(fav_mps_df$MP %in% mps_fav) & fav_mps_df$sample == s, ]$pct))
fav_mps_df$sample = factor(fav_mps_df$sample, levels = unique(fav_mps_df$sample)[order(undiff_sums, decreasing = T)])
fav_mps_df$MP = factor(fav_mps_df$MP, levels = c(mps_unfav, mps_fav))
min_y = min(sapply(unique(fav_mps_df$sample), function(s) sum(fav_mps_df[fav_mps_df$MP %in% mps_fav & fav_mps_df$sample == s, ]$pct)))
p = ggplot(data = fav_mps_df, aes(x = sample, y = pct, fill = MP)) + geom_bar(stat = "identity") + scale_fill_manual(values = mp_cols, name = "MP") + geom_hline(yintercept=0, color="black") +
  ggnewscale::new_scale_fill() + geom_tile(aes(x = sample, y = max(undiff_sums) + 0.06, fill = subtype), width = 1, height = 0.03, inherit.aes = FALSE) + scale_fill_manual(values = subtype_cols, name = "Subtype") + 
  ggnewscale::new_scale_fill() + geom_tile(aes(x = sample, y = max(undiff_sums) + 0.1, fill = subclass), width = 1, height = 0.03, inherit.aes = FALSE) + scale_fill_manual(values = subclass_cols, name = "Subclass") + 
  scale_y_continuous(limits = c(min_y, max(undiff_sums) + 0.12), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c(1, 0.5, 0, 0.5, 1)) +
  labs(x = "", y = "proportion of cells") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), legend.key.size = unit(4, "mm"))
ggsave("01b_barplot_cell-proportions_by_sample_fav-vs-nonfav_scaled_horizontal.png", plot = p, width = 8, height = 5, path = out_dir)

# statistics
# number of ZFTA samples:
mp_cells_df %>% group_by(sample) %>% dplyr::filter(MP == "ZFTA-fus") %>% dplyr::summarise(count = n_distinct(sample)) %>% dplyr::summarise(sum(count))
# number of ZFTA samples where ZFTA-fus is the most dominant program:
mp_cells_df[mp_cells_df$subtype == "ZFTA", ] %>% group_by(sample) %>% slice_max(pct, n = 1, with_ties = FALSE) %>%  # Select the program with the max percentage
  dplyr::filter(MP == "ZFTA-fus") %>% dplyr::summarise(count = n()) # %>% dplyr::summarise(sum(count))
# number of ZFTA samples where ZFTA-fus is expressed by more than 10/30%:
mp_cells_df[mp_cells_df$subtype == "ZFTA", ] %>% group_by(sample) %>% dplyr::filter(MP == "ZFTA-fus" & pct >= 0.1) # %>% ungroup() %>% dplyr::summarise(length(sample))
mp_cells_df[mp_cells_df$subtype == "ZFTA", ] %>% group_by(sample) %>% dplyr::filter(MP == "ZFTA-fus" & pct >= 0.3) # %>% ungroup() %>% dplyr::summarise(length(sample))

# mean & median expression of ZFTA-fus:
x = mp_cells_df[mp_cells_df$subtype == "ZFTA" & mp_cells_df$MP == "ZFTA-fus", ]$pct
print("mean/median expression of ZFTA-fus among ZFTA samples (with sd): "); mean(x); median(x); sd(x)
# same for ECM-remod:
x = mp_cells_df[mp_cells_df$subtype == "ZFTA" & mp_cells_df$MP == "ECM-remod", ]$pct
print("mean/median expression of ECM-remod among ZFTA samples (with sd): "); mean(x); median(x); sd(x)
# number of ZFTA samples where ECM-remod is expressed by more than 8/10%:
mp_cells_df[mp_cells_df$subtype == "ZFTA", ] %>% group_by(sample) %>% dplyr::filter(MP == "ECM-remod" & pct >= 0.08) # %>% ungroup() %>% dplyr::summarise(length(sample))
mp_cells_df[mp_cells_df$subtype == "ZFTA", ] %>% group_by(sample) %>% dplyr::filter(MP == "ECM-remod" & pct >= 0.1) # %>% ungroup() %>% dplyr::summarise(length(sample))

# 1.2 composition by SUBTYPES
# a. ordered by subtype
data = mp_cells_df %>% group_by(subtype, MP) %>% dplyr::mutate(pct = mean(pct, na.rm = T)) %>% distinct(subtype, MP, pct) %>% ungroup()
ggplot(data, aes(x = subtype, y = pct, fill = MP)) + geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = mp_cols, name = "MP") +
  labs(x = "", y = "proportion of cells") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), legend.key.size = unit(4, "mm"))
ggsave("02a_barplot_cell-proportions_by_subtype.png", width = 3, height = 2.5, path = out_dir)

# b. fav. vs. unfav. MPs
fav_mps_df_avg = fav_mps_df %>% group_by(subtype, MP) %>% dplyr::mutate(pct = mean(pct, na.rm = T)) %>% distinct(subtype, MP, pct) %>% ungroup()
undiff_sums = sapply(unique(fav_mps_df_avg$subtype), function(s) sum(fav_mps_df_avg[fav_mps_df_avg$MP %in% mps_unfav & fav_mps_df_avg$subtype == s, ]$pct))
fav_mps_df_avg$subtype = factor(fav_mps_df_avg$subtype, levels = unique(fav_mps_df_avg$subtype)[order(undiff_sums, decreasing = T)])
fav_mps_df_avg$MP = factor(fav_mps_df_avg$MP, levels = c(mps_unfav, mps_fav))
min_y = min(sapply(unique(fav_mps_df_avg$subtype), function(s) sum(fav_mps_df_avg[!fav_mps_df_avg$MP %in% mps_unfav & fav_mps_df_avg$subtype == s, ]$pct)))
ggplot(data = fav_mps_df_avg, aes(x = subtype, y = pct, fill = MP)) + geom_bar(stat = "identity") + scale_fill_manual(values = mp_cols, name = "MP") + geom_hline(yintercept=0, color="black") +
  scale_y_continuous(limits = c(min_y, max(undiff_sums) + 0.12), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c(1, 0.5, 0, 0.5, 1)) + labs(x = "", y = "proportion of cells") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
ggsave("02b_barplot_cell-proportions_by_subtype_fav-vs-unfav_scaled.png", width = 3, height = 2.5, path = out_dir)






                   
