# Description: compute gene contribution to cell type recapitulation for a disease model and cell type
gene_diff_contribution <- function(disease_model_contribution, human_contribution, cell_ideal_recapitulation){
  cell_types <- rownames(cell_ideal_recapitulation)
  output <- data.frame()
  for(b in cell_types){
    # Select gene contributions of cell type b
    cell_to_pick <- paste0("^", b)
    delta_prime <- disease_model_contribution[,colnames(disease_model_contribution)[grep(cell_to_pick, colnames(disease_model_contribution))]]
    delta <- human_contribution[,colnames(human_contribution)[grep(cell_to_pick, colnames(human_contribution))]]
    delta_tilde <- delta_prime-delta # As all delta_tilde are equal pick the first sample
    Gamma <- cell_ideal_recapitulation[rownames(cell_ideal_recapitulation) == b]
    contribution <- data.frame("Cell_type" = rep(b, nrow(delta_tilde)),
                               "Gene_name" = rownames(delta_tilde),
                               "Gene_contribution" = delta_tilde[,1]/Gamma)
    output <- rbind(output, contribution)
  }
  return(output)
}
