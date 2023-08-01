# Adapted from https://github.com/samuel-marsh/scCustomize/tree/master

Cluster_Stats <- function(sobj, group.by.var = "orig.ident") {
  
  # Check on meta data column
  possible_meta_col <- colnames(sobj@meta.data)
  if (!group.by.var %in% possible_meta_col) {
    stop(message = "{.val {group.by.var}} was not found in meta.data slot of Seurat Object.")
  }
  
  # Extract total percents
  total_percent <- prop.table(x = table(sobj@active.ident)) * 100
  total_percent <- data.frame(total_percent) %>%
    dplyr::rename(Cluster = .data[["Var1"]])
  
  # Extract total cell number per cluster across all samples
  total_cells <- table(sobj@active.ident) %>%
    data.frame() %>%
    dplyr::rename(Cluster = .data[["Var1"]], Number = .data[["Freq"]])
  
  # Cluster overall stats across all animals
  cluster_stats <- suppressMessages(dplyr::left_join(total_cells, total_percent))
  
  # Extract cells per metadata column per cluster
  cells_per_cluster_2 <- table(sobj@active.ident, sobj@meta.data[, group.by.var])
  cells_per_cluster_2 <- data.frame(cells_per_cluster_2) %>%
    dplyr::rename(Cluster = .data[["Var1"]], group.by.var = .data[["Var2"]], cell_number = .data[["Freq"]])
  
  cells_per_cluster_2 <- cells_per_cluster_2 %>%
    tidyr::pivot_wider(names_from = group.by.var, values_from = .data[["cell_number"]])
  
  # Merge cells per metadata column per cluster with cluster stats
  cluster_stats_2 <- suppressMessages(dplyr::left_join(cluster_stats, cells_per_cluster_2))
  
  # Calculate and extract percents of cells per cluster per
  percent_per_cluster_2 <- prop.table(x = table(sobj@active.ident, sobj@meta.data[, group.by.var]), margin = 2) * 100
  percent_per_cluster_2 <- data.frame(percent_per_cluster_2) %>%
    dplyr::rename(Cluster = .data[["Var1"]], group.by.var = .data[["Var2"]], percent = .data[["Freq"]])
  percent_per_cluster_2 <- percent_per_cluster_2 %>%
    tidyr::pivot_wider(names_from = group.by.var, values_from = .data[["percent"]]) %>%
    tibble::column_to_rownames("Cluster")
  colnames(percent_per_cluster_2) <- paste(colnames(percent_per_cluster_2), "%", sep = "_")
  
  percent_per_cluster_2 <- percent_per_cluster_2 %>%
    tibble::rownames_to_column(var = "Cluster")
  
  # Merge percent cells per metadata column per cluster with cluster stats and add Totals column
  cluster_stats <- suppressMessages(dplyr::left_join(cluster_stats_2, percent_per_cluster_2)) %>%
    janitor::adorn_totals("row")
  
  return(cluster_stats)
}
