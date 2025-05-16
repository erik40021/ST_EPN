get_spatial_aspect_ratio = function(sstobj) {
  coord <- GetTissueCoordinates(sstobj@images[[1]]) # grabs either 'tissue_lowres_image' or 'slice1'
  ratio = (max(coord[, 1]) - min(coord[, 1])) / (max(coord[, 2]) - min(coord[, 2]))
  return(ratio)
}

# define scale_factor as the size a spot should have in a full-size sample (spanning the full 5.5 x 5.5 mm area)
get_spatial_point_size = function(sstobj, scale_factor = 2) {
  coord <- GetTissueCoordinates(sstobj@images[[1]])
  area = (max(max(coord[, 1]) - min(coord[, 1]), max(coord[, 2]) - min(coord[, 2])))^2
  pt_size = max(scale_factor, scale_factor * 1/sqrt(area/2000000))
  return(pt_size)
}


# --- graph model functions ---
prepare_graph_data = function(asso, pvals, max_type_only = F, p_threshold = 0.05) {
  type_inds = list(coloc = 1, contact = 2:6, dist_contact = 7:15)
  if (max_type_only) {
    type_names = c("coloc", rep("contact", 5), rep("dist_contact", 10))
    abs_max = apply(asso, 1, function(pair) pair[which.max(abs(pair))])
    sign = ifelse(abs_max > 0, "pos.", ifelse(abs_max < 0, "neg.", "ns"))
    type = apply(asso, 1, function(pair) type_names[which.max(abs(pair))])
    avg_asso_in_type = sapply(rownames(asso), function(pair_name) mean(unlist(asso[pair_name, type_inds[[type[pair_name]]]])))
    avg_p_in_type = sapply(rownames(pvals), function(pair_name) mean(unlist(pvals[pair_name, type_inds[[type[pair_name]]]])))
    sign = sapply(names(sign), function(pair_name) if (avg_p_in_type[pair_name] < p_threshold) sign[pair_name] else "ns")
    all_data = data.frame(source = sub(" x .*", "", rownames(asso)), target = sub("^.*? x ", "", rownames(asso)), 
                          weight = avg_asso_in_type, type = type, sign = sign, pval = avg_p_in_type)
  } else {
    types = rep(names(type_inds), each = nrow(asso))
    avg_asso_per_type = unlist(lapply(names(type_inds), function(type) apply(asso, 1, function(pair) mean(pair[type_inds[[type]]]))))
    avg_p_per_type = unlist(lapply(names(type_inds), function(type) apply(pvals, 1, function(pair) mean(pair[type_inds[[type]]]))))
    sign = ifelse(avg_p_per_type < p_threshold, ifelse(avg_asso_per_type > 0, "pos.", "neg."), "ns")
    all_data = data.frame(source = rep(sub(" x .*", "", rownames(asso)), 3), target = rep(sub("^.*? x ", "", rownames(asso)), 3),
                          weight = avg_asso_per_type, type = types, sign = sign, pval = avg_p_per_type)
  }
  return(all_data)
}

# stacks contact associations horizontally, co-localisations vertically
construct_graph <- function(data, root_node, verticality_threshold = 0.5) {
  layers <- list(root_node)
  visited <- c(root_node) # Keep track of visited nodes
  layer_index <- 1 # Current layer index
  while (TRUE) {
    current_layer_nodes = layers[[layer_index]]
    # 1. vertical: add co-localising node(s) to current layer
    while (TRUE) {
      coloc_nodes = data %>% filter((source %in% current_layer_nodes | target %in% current_layer_nodes) & !((source %in% visited & target %in% visited)))
      coloc_nodes <- coloc_nodes %>% mutate(adj_node = if_else(source %in% current_layer_nodes, target, source)) %>% filter(!adj_node %in% visited)
      vertical_weights <- coloc_nodes %>% filter(type == "coloc") %>% group_by(adj_node) %>% dplyr::summarize(vertical_strength = mean(weight), .groups="drop")
      vertical_nodes = vertical_weights[vertical_weights$adj_node %in% coloc_nodes[coloc_nodes$type == "coloc", ]$adj_node & 
                                          vertical_weights$vertical_strength > verticality_threshold, ]$adj_node
      if (length(vertical_nodes) == 0) break
      current_layer_nodes = c(current_layer_nodes, vertical_nodes)
      layers[[layer_index]] = current_layer_nodes; visited <- c(visited, vertical_nodes)
      coloc_nodes = data %>% filter((source %in% current_layer_nodes | target %in% current_layer_nodes) & !((source %in% visited & target %in% visited)))
    }
    # 2. horizontal: add the highest associated contact node to form the next layer
    grouped_edges <- data %>% filter((source %in% current_layer_nodes | target %in% current_layer_nodes) & !((source %in% visited & target %in% visited)))
    if (nrow(grouped_edges) == 0) { message("no more associations found to layer ", layer_index); break }
    grouped_edges <- grouped_edges %>% mutate(adj_node = if_else(source %in% current_layer_nodes, target, source)) %>% filter(!adj_node %in% visited & type == "contact")
    horizontal_weight <- grouped_edges %>% filter(type == "contact") %>% group_by(adj_node) %>% dplyr::summarize(total_strength = mean(weight), .groups = "drop")
    if (nrow(horizontal_weight) == 0) { message("no more associations found to layer ", layer_index); break }
    next_layer_node = horizontal_weight[which.max(horizontal_weight$total_strength), ]$adj_node
    visited <- c(visited, next_layer_node)
    layers[[layer_index + 1]] <- next_layer_node
    layer_index <- layer_index + 1
  }
  node_positions <- integer(0); names(node_positions) <- character(0)
  node_positions = sapply(1:length(layers), function(i) layers[[i]])
  nodes_layers = unlist(sapply(1:length(node_positions), function(i) rep(i, length.out = length(node_positions[[i]]))))
  names(nodes_layers) = unlist(node_positions)
  nodes_df <- data.frame(id = mp_names); nodes_df$layer = nodes_layers[match(nodes_df$id, names(nodes_layers))]
  nodes_df = nodes_df[order(nodes_df$layer),]
  nodes_df$layer[is.na(nodes_df$layer)] = layer_index + 1
  return(nodes_df)
}
