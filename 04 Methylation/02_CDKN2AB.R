library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # reference annotation for EPIC methylation data
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # reference annotation for 450k methylation data
library(minfi)
library(conumee2)
library(ggplot2)
library(dplyr)

metadata = readxl::read_excel("masterlist_EPN-ST.xlsx", skip = 1); metadata = metadata[!is.na(metadata$IDAT) & metadata$IDAT != "see primary", ]
metadata$patiend_id = sub("(\\d+)[A-Za-z].*$", "\\1", metadata$`Study ID`)


# --- load and pre-process raw data ---
# 1. study samples (n=49)
idat_dir_450k = "450K"
idat_dir_epic = "EPIC (850K)"
RGset_450k = read.metharray.exp(base = idat_dir_450k, recursive = T, force = T)
RGset_EPIC = read.metharray.exp(base = idat_dir_epic, recursive = T, force = T)
RGset_450k <- preprocessNoob(RGset_450k); RGset_EPIC <- preprocessNoob(RGset_EPIC) # 1. normalise independently
common_probes <- intersect(rownames(getBeta(RGset_450k)), rownames(getBeta(RGset_EPIC))) # 2. find common probes
RGset_450k <- RGset_450k[common_probes, ]; RGset_EPIC <- RGset_EPIC[common_probes, ] # 3. subset methylation sets to common probes
RGset <- combineArrays(RGset_450k, RGset_EPIC) # 4. merge back into one common-probes (~450k) dataset
RGset = RGset[, metadata$IDAT] # reorder to match metadata+
message(dim(RGset)[2], " samples collected in combined RGset, each containing ", dim(RGset)[1], " probes")

# 2. normal reference brain samples (n=4+12)
# source 1: TCGA adult normal brain methylomes (n = 4, see https://portal.gdc.cancer.gov/, case IDs: C3N-03450, C3N-03026, TCGA-74-6573, C3N-03446)
# source 2: children brain methylomes from Franklin et al. 25 (n = 12, see https://www.biorxiv.org/content/10.1101/2025.02.21.639467v1, all cases with age > 0)
idat_dir_ref_450k = "450k"
idat_dir_ref_epic = "EPIC"
RGset_ref_450k = read.metharray.exp(base = idat_dir_ref_450k, recursive = TRUE, targets = NULL, force = T) # use "recursive" = TRUE to read in all files in the dir
RGset_ref_epic = read.metharray.exp(base = idat_dir_ref_epic, recursive = TRUE, targets = NULL, force = T) # use "recursive" = TRUE to read in all files in the dir
RGset_ref_450k <- preprocessNoob(RGset_ref_450k); RGset_ref_epic <- preprocessNoob(RGset_ref_epic) # 1. normalise independently
common_probes <- intersect(rownames(getBeta(RGset_ref_450k)), rownames(getBeta(RGset_ref_epic))) # 2. find common probes
RGset_ref_450k <- RGset_ref_450k[common_probes, ]; RGset_ref_epic <- RGset_ref_epic[common_probes, ] # 3. subset methylation sets to common probes
RGset_ref <- combineArrays(RGset_ref_450k, RGset_ref_epic) # 4. merge back into one common-probes (~450k) dataset
message(dim(RGset_ref)[2], " samples collected in combined reference RGset, each containing ", dim(RGset_ref)[1], " probes")

common_probes <- intersect(rownames(getBeta(RGset)), rownames(getBeta(RGset_ref)))
RGset <- RGset[common_probes, ]; RGset_ref <- RGset_ref[common_probes, ] # subset methylation sets to common probes
# RGset_all = combineArrays(RGset2, RGset_ref2)
cnv_obj <- CNV.load(RGset)
cnv_ref_obj = CNV.load(RGset_ref)
data("detail_regions")
anno <- CNV.create_anno(detail_regions = detail_regions) # must be based on hg19 because of 450k data
gr <- anno@probes[anno@probes@ranges@NAMES %in% rownames(cnv_obj@intensity)]
anno@probes <- subsetByOverlaps(anno@probes, gr) # subset anno to match probes available

cdkn2_ranges = detail_regions[detail_regions$name == "CDKN2A/B"]
cdkn2_probes <- subsetByOverlaps(anno@probes, GRanges(seqnames = "chr9", cdkn2_ranges@ranges))
cdkn2_thick_probes = subsetByOverlaps(anno@probes, GRanges(seqnames = "chr9", cdkn2_ranges$thick))
cdkn2_neighbours = setdiff(names(cdkn2_thick_probes), names(cdkn2_probes))

# specify exact probes of CDKN2A and CDKN2B:
cdkn2a_probes = names(cdkn2_probes[1:8])
cdkn2b_probes = setdiff(names(cdkn2_probes), cdkn2a_probes)

cdkn2ab_genmomic_loci = c(21968233, 21995735, 22009595)
cnv_cdkn2 = lapply(metadata$`Study ID`, function(s) {
  message("running ", s)
  x <- CNV.fit(cnv_obj[as.character(metadata[metadata$`Study ID` == s, "IDAT"])], cnv_ref_obj, anno)
  x <- CNV.bin(x); x <- CNV.detail(x); x <- CNV.segment(x)
  x@fit$ratio[names(cdkn2_probes), ]
})

cnv_cdkn2 = as.data.frame(cnv_cdkn2); colnames(cnv_cdkn2) = metadata$`Study ID`
avg_cdkn2_ratio = colMeans(cnv_cdkn2) # mean copy-number signal (ratio) of all CDKN2A/B probes 

