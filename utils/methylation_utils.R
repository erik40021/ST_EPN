
# assign majority known tumor class to each cluster
assign_classes <- function(clusters, meth_class) {
  cluster_classes <- sapply(unique(clusters)[order(unique(clusters))], function(cluster) {
    samples_in_cluster <- which(clusters == cluster)
    known_labels <- meth_class[samples_in_cluster]
    if (all(is.na(known_labels))) {
      return(NA)  # No known labels in this cluster
    }
    return(names(sort(table(known_labels), decreasing = TRUE))[1])  # Majority vote
  })
  # Map each sample's cluster number to the corresponding class
  assigned_classes <- cluster_classes[clusters]
  return(assigned_classes)
}

run_louvain_iterations <- function(graph, res = 1, n_iter = 100) {
  library(igraph)
  clusters_list <- list()
  for (i in seq_len(n_iter)) clusters_list[[i]] = membership(cluster_louvain(graph, weights = E(graph)$weight, resolution = res))
  return(clusters_list)
}
run_kmeans_iterations <- function(pca_data, k = 3, n_iter = 100) {
  clusters_list <- list()
  for (i in seq_len(n_iter)) clusters_list[[i]] = kmeans <- kmeans(pca_data, centers = k)$cluster
  return(clusters_list)
}

calculate_class_consistency <- function(class_matrix) {
  confidence_scores <- apply(class_matrix, 1, function(row) {
    valid_classes <- row[!is.na(row)]  # Exclude NA values
    if (length(valid_classes) == 0) return(NA)  # No consistent class
    max(table(valid_classes)) / length(valid_classes)  # Fraction of consistent assignments
  })
  return(confidence_scores)
}
