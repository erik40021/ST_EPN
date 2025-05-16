library(RCy3)

source("utils/r_utils.R")
source("utils/spatial_utils.R")
source("utils/spatial_score_utils.R")


base_dir = ""

metadata = readxl::read_excel("masterlist_EPN-ST.xlsx", skip = 1); metadata = metadata[grepl("2", metadata$`Analyses*`), ]
sample_ids = metadata$`Study ID`
metaprograms = readxl::read_excel("spatial_NMF_metaprograms_main.xlsx", sheet = "MP genes (final)")
mp_names = colnames(metaprograms)
type_cols = c(coloc = "#C43D3C", contact = "#208cc0", dist_contact = "#ffd60a")
layer_cols = c("#99582a", "#fa7764", "#f1af3a", "#82a5b5", "#b5c3c9", "#617680", "grey60", "grey80", "grey90", "grey95")
association = readRDS(file.path(base_dir, "association.rds"))



# 1. unbiased model construction
# using all-sample-averaged associations
avg_asso = readRDS(file.path(base_dir, "avg_association.rds"))
avg_pvals = combine_pvals_dfs(adj_pvals)
all_data = prepare_graph_data(avg_asso, avg_pvals)
graph = construct_graph(all_data, root_node = "Hypoxia", verticality_threshold = 0.5)
graph = graph[match(mp_names, graph$id), ]

# 2. visualisation in Cytoscape
# reduce (over-)crowdedness of associations by showing edges only between nodes of the same layer for coloc, of neighb. layers for contact, or of non-neighbouring layers for proximal and distant
data = all_data[all_data$sign == "pos.", ]
data$is_same_layer = apply(data, 1, function(edge) graph[graph$id == edge[1], 2] == graph[graph$id == edge[2], 2])
data$is_neighb_layer = apply(data, 1, function(edge) abs(graph[graph$id == edge[1], 2] - graph[graph$id == edge[2], 2]) == 1)
data$is_dist_layer = apply(data, 1, function(edge) abs(graph[graph$id == edge[1], 2] - graph[graph$id == edge[2], 2]) > 1)
data = data[(data$type == "coloc" & data$is_same_layer) | (data$type == "contact" & data$is_neighb_layer) | (data$type == "dist_contact" & data$is_dist_layer), ]

createNetworkFromDataFrames(edges = data, nodes = graph, title = "diff_interactions", collection = "all_ints_networks")
setVisualStyle("default"); setCurrentNetwork("diff_interactions")
setEdgeLineWidthMapping("weight", c(0, max(data$weight)), c(0, 10)) # Anpassung der Kantendicke
setEdgeColorMapping("type", mapping.type = "d", table.column.values = names(type_cols), colors = type_cols)
setNodeBorderColorDefault("#000000"); setNodeBorderWidthDefault(1)
setNodeColorBypass(graph$id, layer_cols[graph$layer])
# manual coloring:
cols = c(Hypoxia="#99402a", Stress="#99402a",'Inflammatory'="#fa7764", Myeloid="#f1c33a", Vasc="#f1c33a", Fibroblast="#f1c33a",
         Glial="grey90", 'Interferon-re'="#b5c3c9", Astroglial="#82a5b5", 'ZFTA-fus'="grey75", Cilia="#788d96", 'Cell-cycle'="#617680")
setNodeColorBypass(names(cols), cols)



