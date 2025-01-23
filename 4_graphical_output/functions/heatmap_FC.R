# Description: plot FC for each gene, cell type, disease model and pathway
heatmap_FC <- function(df, pathway_full_name = "Dendritic cells in Th1/Th2 Development [BIOCARTA]",
                         value_size = 5, gene_size = 8, height = NULL, width = NULL){
  df$Sample <- paste0(df$Disease_model, "_", df$cell_type)
  df_ordered <- df[order(df$Pathway_full_name, df$Disease_model, df$cell_type, df$gene_name), ]
  df_2 <- df_ordered[df_ordered$Pathway_full_name == pathway_full_name, 
                     colnames(df_ordered) %in% c("Sample", "logFC", "gene_name")]
  df_2$logFC <- df_2$logFC
  # Pivot to wider format
  df_wide <- df_2 %>%
    pivot_wider(names_from = Sample, values_from = logFC)
  # Transpose the dataframe
  df_transposed <- as.data.frame(df_wide)
  genes <- df_transposed$gene_name
  df_transposed$gene_name <- NULL
  df_transposed <- as.matrix(df_transposed)
  rownames(df_transposed) <- genes
  # Metadata
  clusters <- rep(c("Dendritic Cells", "Keratinocytes", "Langerhans Cells", "Melanocytes", "T-cell"), times = 3)
  stimulation <- rep(c("IMQ", "OVA", "OXA"), each = 5)
  
  # Function to create cell borders and remove "NaN" text
  # Matrix with in silico sum values
  # Pivot to wider format
  df_2 <- df_ordered[df_ordered$Pathway_full_name == pathway_full_name, 
                     colnames(df_ordered) %in% c("Sample", "InSilicoSum", "gene_name")]
  border_df <- df_2 %>%
    pivot_wider(names_from = Sample, values_from = InSilicoSum)
  # Transpose the dataframe
  border_matrix <- as.data.frame(border_df)
  genes <- border_matrix$gene_name
  border_matrix$gene_name <- NULL
  border_matrix <- as.matrix(border_matrix)
  rownames(border_matrix) <- genes
  
  # Create matrix to check if the row should be black for each disease model
  is_row_black_per_disease <- apply(border_matrix, 2, function(column) border_matrix[, colnames(border_matrix) %in% c("IMQ", "OVA", "OXA")] == -1)
                                      
  # Updated cell_fun function
  library(grid)
  cell_fun <- function(j, i, x, y, width, height, fill) {
    disease_model <- sub("_.*", "", colnames(df_transposed)[j])
    if (any(border_matrix[i, grep(disease_model, colnames(border_matrix))] == -1, na.rm = TRUE)) {
      grid.rect(x, y, width, height, gp = gpar(fill = "grey", col = "grey"))
    } else if (is.na(df_transposed[i, j])) {
      grid.rect(x, y, width, height, gp = gpar(fill = "grey", col = "grey"))
    } else if (border_matrix[i, j] == 2) {
      grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 1))
    }
    # Only show text value for non-zero FC
    if(abs(df_transposed[i,j]) > .Machine$double.eps && !is.na(df_transposed[i,j])){
        grid.text(sprintf("%.1f", df_transposed[i,j]), x, y, 
                  gp = gpar(fontsize = value_size))
        
      }
    }
  
  # Annotation object
  ha_col <- HeatmapAnnotation(
    "Cell type" = clusters,
    "Disease model" = stimulation,
    col = list("Cell type" = c("Dendritic Cells" = "#1f77b4", "Keratinocytes" = "#ff7f0e", "Langerhans Cells" = "#2ca02c", "Melanocytes" =  "#17becf", "T-cell" = "#d62728"),
               "Disease model" = c("IMQ" = "#9467bd", "OVA" = "#8c564b", "OXA" = "#e377c2")),
    annotation_name_gp = gpar(fontsize = 8),  # Change the font size of the annotation labels
    annotation_legend_param = list(
      "Cell type" = list(title_gp = gpar(fontsize = 8, fontface = 2), labels_gp = gpar(fontsize = 8)),
      "Disease model" = list(title_gp = gpar(fontsize = 8, fontface = 2), labels_gp = gpar(fontsize = 8))
    )
  )
  
  # Custom color function
  color_function <- function(x){
    colors <- colorRamp2(c(-100, 0, 100), c("#1f77b4", "#ffffff", "#d62728"))
    col <- colors(x)
    col[is.na(x)] <- "grey"
    return(col)
  }
  
  # Create heatmap
  ht <- Heatmap(df_transposed, top_annotation = ha_col, cluster_rows = FALSE, cluster_columns = FALSE,
                height = height, width = width, column_title = "B",            
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),  # Customize title appearance
                show_column_names = FALSE, row_names_gp = gpar(fontsize = gene_size, fontface = 2),
                col = color_function,   
                cell_fun = cell_fun,
                heatmap_legend_param = list(title = expression(r[tilde(g)]^{b}), title_gp = gpar(fontface = 2))
                #layer_fun = function(j, i, x, y, width, height, fill) {
                # since grid.text can also be vectorized
                #  grid.text(sprintf("%.1f", pindex(df_transposed, i, j)), x, y, 
                #            gp = gpar(fontsize = value_size))
                #}
  )
  return(ht)
}
