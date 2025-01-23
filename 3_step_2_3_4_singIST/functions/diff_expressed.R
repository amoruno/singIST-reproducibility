diff_expressed <- function(SeuratObject, condition_1 = c(), condition_2 = c(), adjpval.treshold = 0.05,
                           logfc.treshold = 0.05){
  # FindMarkers function by row of dataset
  apply_function <- function(row, data = SeuratObject) {
    logFC <- FindMarkers(object = data, ident.1 = row[1], ident.2 = row[2], 
                         assay = "RNA", slot = "data", logfc.threshold = logfc.treshold)
    # Filter only genes with adjpval.treshold
    logFC <- logFC[logFC$p_val_adj <= adjpval.treshold, ]
    return(logFC)
  }
  # Combinations to test
  combinations <- outer(condition_1, condition_2 , paste, sep = "_")
  
  # Apply the function to each row of the data frame
  output <- apply(combinations, 1, apply_function)
  names(output) <- condition_1
  return(output)
}
