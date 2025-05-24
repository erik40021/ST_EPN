library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # reference annotation for EPIC methylation data
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # reference annotation for 450k methylation data
library(minfi)
library(Seurat)
library(dplyr)

setwd("")
source("utils/methylation_utils.R")


metadata = readxl::read_excel("masterlist_EPN-ST.xlsx", skip = 1); metadata = metadata[!is.na(metadata$IDAT) & metadata$IDAT != "see primary", ]
metadata$patiend_id = sub("(\\d+)[A-Za-z].*$", "\\1", metadata$`Study ID`)


out_dir = ""
pohl_metadata = as.data.frame(readxl::read_excel(file.path(out_dir, "Pohl_metadata.xls"), skip = 22)) # cohort metadata from Pohl et al. 2024
pohl_metadata = pohl_metadata[grepl("ZFTA|YAP|ST_SE", pohl_metadata$`characteristics: Classification result v12`), ]
rownames(pohl_metadata) = sub("_Grn\\.idat$", "", pohl_metadata$`idat file 1`)
ref_dir = ""

# --- PART 1: prepare meth. data (only needed once) ---
# 1. prepare ref meth. data
ref_files = list.files(ref_dir, pattern = ".idat$", recursive = T, full.names = F)
idat_red = sub(".*GSM[0-9]+[_/]", "", ref_files); idat_red = idat_red[grepl("Red", idat_red)]
message("samples not found in metadata: ", idat_red[!idat_red %in% pohl_metadata$`idat file 2`]) # check that all ref files have metadata (should be empty here!)
platform = pohl_metadata[match(idat_red, pohl_metadata$`idat file 2`), ]$platform
basenames <- unique(sub("_(Grn|Red)\\.idat$", "", list.files(ref_dir, pattern = ".idat$", recursive = T, full.names = T)))
pohl_metadata = pohl_metadata[pohl_metadata$`idat file 2` %in% idat_red, ] # subset metadata to only those samples where IDATs exist (41 couldn't be downloaded!)
ref_450k = basenames[platform == "GPL13534"]; ref_EPIC = basenames[platform == "GPL21145"]
# 2. prepare our meth. data
ms_450k_full = list.files("450k", pattern = ".idat$", recursive = T, full.names = T)
ms_EPIC_full = list.files("EPIC (850K)", pattern = ".idat$", recursive = T, full.names = T)
ms_450k = unique(sub("_(Grn|Red)\\.idat$", "", ms_450k_full)); ms_EPIC = unique(sub("_(Grn|Red)\\.idat$", "", ms_EPIC_full))
# 3. read in all samples together (MS + ref), split by platforms
all_450k = c(ref_450k, ms_450k); all_EPIC = c(ref_EPIC, ms_EPIC)
sum(duplicated(all_450k)); sum(duplicated(all_EPIC)) # check if there are duplicates (samples in both cohorts) and remove from REF if true
RGset_450k = read.metharray(basenames = all_450k, force = T)
RGset_EPIC = read.metharray(basenames = all_EPIC, force = T)
RGset_450k <- preprocessNoob(RGset_450k); RGset_EPIC <- preprocessNoob(RGset_EPIC) # normalise independently
common_probes <- intersect(rownames(getBeta(RGset_450k)), rownames(getBeta(RGset_EPIC))) # find common probes
# subset methylation sets to common probes
RGset_450k <- RGset_450k[common_probes, ]; RGset_EPIC <- RGset_EPIC[common_probes, ]
# 4. merge back into one common-probes (450k) dataset
RGset_combined <- combineArrays(RGset_450k, RGset_EPIC) # casts EPIC arrays into 450k type too
colnames(RGset_combined) = sub("^GSM[0-9]+_", "", colnames(RGset_combined))
saveRDS(RGset_combined, file.path(out_dir, "RGset_combined.rds"))
# 5. combine metadata dfs of MS and REF cohorts
colnames(metadata)[colnames(metadata) %in% c("Subclass (classifier)", "Age at resection", "Sex", "Meth. classif. score")] = c("sex", "age", "meth_class", "class_score")
metadata = metadata[metadata$IDAT %in% colnames(RGset_combined), ]
combined_metadata = data.frame(ID = c(metadata$patiend_id, paste0("ref_", 1:length(basenames)))); rownames(combined_metadata) = combined_metadata$ID
combined_metadata = cbind(combined_metadata, IDAT = c(metadata$IDAT, sub("_(Grn|Red)\\.idat$", "", pohl_metadata$`idat file 1`)),
                          age = c(metadata$age, pohl_metadata$`characteristics: Age`),
                          sex = c(metadata$sex, pohl_metadata$`characteristics: Sex`),
                          meth_class = c(metadata$meth_class, pohl_metadata$`characteristics: Classification result v12`),
                          class_score = c(metadata$class_score, as.numeric(sub("'", "", pohl_metadata$`characteristics: Classification score`))),
                          os = c(rep(NA, nrow(metadata)), pohl_metadata$`characteristics: OS [months]`),
                          pfs = c(rep(NA, nrow(metadata)), pohl_metadata$`characteristics: PFS [months]`),
                          os_status = c(rep(NA, nrow(metadata)), as.numeric(pohl_metadata$`characteristics: OS status`)),
                          pfs_status = c(rep(NA, nrow(metadata)), as.numeric(pohl_metadata$`characteristics: PFS status`)))
combined_metadata$meth_class = stringr::str_replace_all(combined_metadata$meth_class, c(EPN_ST_ZFTA_FUS_C = "ZFTA_FUS_C", EPN_ST_ZFTA_FUS_D = "ZFTA_FUS_D",
    EPN_ST_ZFTA_FUS_E = "ZFTA_FUS_E", EPN_ST_ZFTA_RELA_A = "ZFTA_RELA_A", EPN_ST_ZFTA_RELA_B = "ZFTA_RELA_B", EPN_YAP = "YAP1", EPN_ST_SE = "ST_SE"))
combined_metadata$sex = toupper(combined_metadata$sex)
combined_metadata$cohort = "REF"; combined_metadata[!grepl("ref", rownames(combined_metadata)), ]$cohort = "KK"
saveRDS(combined_metadata, file.path(out_dir, "combined_metadata.rds"))
# 6. reorder RGset and extract beta values
RGset_combined = RGset_combined[, combined_metadata$IDAT] # reorder samples (cols) to match metadata order
saveRDS(RGset_combined, file.path(out_dir, "RGset_combined.rds")) # (optional) save to load faster next time
message(dim(RGset_combined)[2], " samples collected in RGset_combined containing ", dim(RGset_combined)[1], " probes")
beta_values = getBeta(RGset_combined); beta_values = beta_values[, !duplicated(colnames(beta_values))] # remove ref samples that are in our cohort too
saveRDS(beta_values, file.path(out_dir, "beta_values.rds"))


# --- PART 2: analyse combined data set ---
beta_values = readRDS(file.path(out_dir, "beta_values.rds")); beta_values = beta_values[, !duplicated(colnames(beta_values))]
combined_metadata = readRDS(file.path(out_dir, "combined_metadata.rds")); combined_metadata = combined_metadata[!duplicated(combined_metadata$IDAT), ]

# 1. calculate the top most variable features for downstream analysis
colnames(beta_values) = combined_metadata$ID
feature_variances <- apply(beta_values, 1, var)
top_features <- names(sort(feature_variances, decreasing = TRUE))[1:100000]
saveRDS(top_features, file = file.path(out_dir, "top_var_features.rds"))
top_betas = beta_values[top_features, ]

# 2. PCA
pca_res <- prcomp(t(top_betas), scale. = TRUE) # perform PCA on the top 100,000 features
explained_variance <- pca_res$sdev^2 / sum(pca_res$sdev^2); cumulative_variance <- cumsum(explained_variance)
png(file.path(out_dir, "elbow_plot_pca.png"), width = 1000, height = 1000)
plot(cumulative_variance, type = "b", pch = 19, frame = FALSE, xlab = "Number of Principal Components", ylab = "Cumulative Proportion of Variance Explained", main = "Elbow Plot of Principal Components")
abline(h = 0.9, col = "red", lty = 2) # Add a horizontal line at 90% cumulative variance
dev.off()
pca_data <- pca_res$x[, 1:40] # Use top principal components for clustering

tsne = Rtsne::Rtsne(pca_data[, 1:40], perplexity = 50, verbose = TRUE, theta = 0, pca = F, iterations = 2000) # option 2: t-SNE
meth_data = cbind(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2], combined_metadata)

# 3. multi-iteration clustering with confidence calculation
meth_class = meth_data$meth_class; meth_class[meth_data$class_score <= 0.8] = NA
meth_class[grepl("ZFTA", meth_class)] = "ZFTA"; meth_class[grepl("ST_SE", meth_class)] = "ST_SE"; meth_class[grepl("YAP1", meth_class)] = "YAP1"

# 3.1 louvain clustering based on top 40 PCs distance matrix (using the graph built above)
clusters_list <- run_louvain_iterations(graph, res = 1.1, n_iter = 1000)
assigned_classes <- lapply(clusters_list, assign_classes, meth_class) # Map clusters to classes for each iteration
class_matrix_louvain_pc <- do.call(cbind, assigned_classes)
# 3.2 kmeans clustering
clusters_list <- run_kmeans_iterations(pca_data, k = 10, n_iter = 1000)
assigned_classes <- lapply(clusters_list, assign_classes, meth_class)
class_matrix_kmeans <- do.call(cbind, assigned_classes)
class_matrix_combined = cbind(class_matrix_louvain_pc, class_matrix_kmeans)
confidence_scores <- calculate_class_consistency(class_matrix_combined) # calculate consistency scores based on assigned classes

