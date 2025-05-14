library(Seurat)
library(writexl)
library(NMF)
library(patchwork)

# recommended to be run in parallel on server, with one process for each sample s
s = commandArgs(trailingOnly = TRUE)[1] # ID of sample to run NMF for (first command line argument when run on server)

in_path = ""
out_path = ""

options(bitmapType = 'cairo')
range = 2:20
n_features = 5000

message("Running NMF for sample ", s)
out_dir = file.path(out_path, s); if (!dir.exists(out_dir)) { dir.create(out_dir, r = T) }

# 1. ---- load and prepare data -----
sstobj = readRDS(file = paste0(in_path, "/sstobj_", s, ".rds"))
data = as.matrix(GetAssayData(sstobj, layer = 'scale.data'))
data = data[VariableFeatures(sstobj)[1:n_features], ]
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
