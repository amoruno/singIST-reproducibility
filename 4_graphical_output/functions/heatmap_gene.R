# Description: create complexheatmap to plot gene contribution to cell type recapitulation for each disease model and pathway
heatmap_gene <- function(df, pathway_full_name = "Dendritic cells in Th1/Th2 Development [BIOCARTA]",
                         value_size = 5, gene_size = 8, height = NULL, width = NULL){
  df$Sample <- paste0(df$disease_model, "_", df$Cell_type)
  df_ordered <- df[order(df$Pathway_full_name, df$disease_model, df$Cell_type, df$Gene_name), ]
  df_2 <- df_ordered[df_ordered$Pathway_full_name == pathway_full_name, 
                      colnames(df_ordered) %in% c("Sample", "Gene_contribution", "Gene_name")]
  df_2$Gene_contribution <- 100*df_2$Gene_contribution
  # Pivot to wider format
  df_wide <- df_2 %>%
    pivot_wider(names_from = Sample, values_from = Gene_contribution)
  # Transpose the dataframe
  df_transposed <- as.data.frame(df_wide)
  genes <- df_transposed$Gene_name
  df_transposed$Gene_name <- NULL
  df_transposed <- as.matrix(df_transposed)
  rownames(df_transposed) <- genes
  # Metadata
  clusters <- rep(c("Dendritic Cells", "Keratinocytes", "Langerhans Cells", "Melanocytes", "T-cell"), times = 3)
  stimulation <- rep(c("IMQ", "OVA", "OXA"), each = 5)

  # Annotation object
  ha_col <- HeatmapAnnotation(
    "Cell type" = clusters,
    "Disease model" = stimulation,
    col = list("Cell type" = c("Dendritic Cells" = "#1f77b4", "Keratinocytes" = "#ff7f0e", "Langerhans Cells" = "#2ca02c", "Melanocytes" =  "#17becf", "T-cell" = "#d62728"),
               "Disease model" = c("IMQ" = "#9467bd", "OVA" = "#8c564b", "OXA" = "#e377c2")),
    annotation_name_gp = gpar(fontsize = 0),  # Change the font size of the annotation labels
    annotation_legend_param = list(
      "Cell type" = list(title_gp = gpar(fontsize = 8, fontface = 2), labels_gp = gpar(fontsize = 8)),
      "Disease model" = list(title_gp = gpar(fontsize = 8, fontface = 2), labels_gp = gpar(fontsize = 8))
    )
  )

  # Updated cell_fun function
  library(grid)
  cell_fun <- function(j, i, x, y, width, height, fill) {
    if(abs(df_transposed[i,j]) > .Machine$double.eps){
      grid.text(sprintf("%.1f", df_transposed[i,j]), x, y, 
                            gp = gpar(fontsize = value_size))
      
    }
  }
  
  # Create heatmap
  ht <- Heatmap(df_transposed, top_annotation = ha_col, cluster_rows = FALSE, cluster_columns = FALSE,
                height = height, width = width, column_title = "A",          
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),  # Customize title appearance
                show_column_names = FALSE, row_names_gp = gpar(fontsize = gene_size, fontface = 2),
                col = colorRamp2(c(-100, 0, 100), c("#1f77b4", "#ffffff", "#d62728")),   
                heatmap_legend_param = list(title = expression(hat(Delta)*f["g,AD"]^b), title_gp = gpar(fontface = 2)),
                cell_fun = cell_fun
                #layer_fun = function(j, i, x, y, width, height, fill) {
                  # since grid.text can also be vectorized
                #  grid.text(sprintf("%.1f", pindex(df_transposed, i, j)), x, y, 
                #            gp = gpar(fontsize = value_size))
                #  }
  )
  return(ht)
}
