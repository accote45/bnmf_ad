# figure1_panel_d_cluster_freq.R
# Panel D: Horizontal bar plot showing how many ancestries each cluster
# appears in (based on cosine-similarity matching across ancestries).

# Assumes figure_utils.R has been sourced.


#' Create Panel D: Cluster frequency across ancestries
#'
#' Assigns canonical cluster IDs based on cosine similarity matching,
#' then counts how many ancestries each canonical cluster appears in.
#'
#' @param h_list Named list of H-matrix tibbles (from load_h_matrices)
#' @param cluster_matches tibble from match_all_clusters()
#' @return ggplot object
plot_panel_d <- function(h_list, cluster_matches) {

  # Build a mapping: each (ancestry, cluster) -> canonical cluster ID
  # Strategy: use union-find style grouping from pairwise matches

  # Collect all (ancestry, cluster) nodes
  all_nodes <- map_dfr(names(h_list), function(anc) {
    tibble(ancestry = anc, cluster = h_list[[anc]]$Cluster)
  }) %>%
    mutate(node_id = paste(ancestry, cluster, sep = "_"))

  # Build adjacency from matched pairs (only real matches, not NA)
  valid_matches <- cluster_matches %>% filter_out(when_any(is.na(cluster_1), is.na(cluster_2)))

  # Simple union-find
  parent <- set_names(all_nodes$node_id, all_nodes$node_id)

  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]  # path compression
      x <- parent[x]
    }
    x
  }

  union_nodes <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) parent[ra] <<- rb
  }

  for (i in seq_len(nrow(valid_matches))) {
    node_a <- paste(valid_matches$anc1[i], valid_matches$cluster_1[i], sep = "_")
    node_b <- paste(valid_matches$anc2[i], valid_matches$cluster_2[i], sep = "_")
    if (node_a %in% names(parent) && node_b %in% names(parent)) {
      union_nodes(node_a, node_b)
    }
  }

  # Assign canonical IDs
  all_nodes <- all_nodes %>%
    mutate(canonical = map_chr(node_id, find_root))

  # Relabel canonical clusters as "Cluster 1", "Cluster 2", etc.
  canonical_labels <- unique(all_nodes$canonical)
  label_map <- set_names(paste("Cluster", seq_along(canonical_labels)), canonical_labels)
  all_nodes <- all_nodes %>%
    mutate(canonical_label = label_map[canonical])

  # Count ancestries per canonical cluster
  freq_df <- all_nodes %>%
    group_by(canonical_label) %>%
    summarise(count = n_distinct(ancestry), .groups = "drop")

  max_ancestries <- length(h_list)

  ggplot(freq_df, aes(x = count, y = reorder(canonical_label, count))) +
    geom_col(fill = "#5B8FA8", width = 0.6) +
    scale_x_continuous(
      breaks = seq(0, max_ancestries),
      limits = c(0, max_ancestries + 0.5)
    ) +
    coord_flip() +
    theme_big_text() +
    labs(x = "Number of Ancestries", y = "Cluster")
}
