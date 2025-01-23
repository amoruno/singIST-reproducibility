# Description: function to select top N genes within each cell type and ensure common genes across all cell types
select_top_genes <- function(data, pathway_full_name = "JAK-STAT signaling pathway [KEGG]", top_n = 10) {
  # Initialize an empty data frame to store results
  top_genes <- data.frame()
  # Filter out other pathways
  data <- data[data$Pathway_full_name == pathway_full_name, ]
  # List to store top genes from each cell type
  top_gene_list <- list()
  
  # Loop over each cell type
  for (cell_type in unique(data$Cell_type)) {
    # Filter data for the current cell type
    cell_data <- data %>% filter(Cell_type == cell_type)
    
    # Select top N genes based on the absolute value of Gene_contribution
    top_cell_genes <- cell_data %>%
      arrange(desc(abs(Gene_contribution))) %>%
      head(top_n)
    
    # Append the selected genes to the list
    top_gene_list[[cell_type]] <- top_cell_genes$Gene_name
    
    # Append the selected genes to the result
    top_genes <- bind_rows(top_genes, top_cell_genes)
  }
  
  # Get a unique list of top genes across all cell types
  unique_top_genes <- unique(unlist(top_gene_list))
  
  # Filter the original data to include only the top genes
  final_top_genes <- data %>% filter(Gene_name %in% unique_top_genes)
  
  # Ensure that all top genes are present for each cell type by joining the missing genes
  final_result <- data.frame()
  for (cell_type in unique(data$Cell_type)) {
    cell_data <- final_top_genes %>% filter(Cell_type == cell_type)
    missing_genes <- setdiff(unique_top_genes, cell_data$Gene_name)
    if (length(missing_genes) > 0) {
      missing_data <- data %>%
        filter(Gene_name %in% missing_genes, Cell_type == cell_type) %>%
        arrange(desc(abs(Gene_contribution)))
      cell_data <- bind_rows(cell_data, missing_data)
    }
    final_result <- bind_rows(final_result, cell_data)
  }
  
  return(final_result)
}
