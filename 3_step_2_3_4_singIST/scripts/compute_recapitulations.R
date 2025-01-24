# Load logFC in disease model
~/singIST-reproducibility/3_step_2_3_4_singIST/exported_results/
all_input_disease_model <- paste0("C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/input_disease_model/", list.files("C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/input_disease_model/", pattern = "\\.RData$"))
all_input_human <- paste0("C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/output/", list.files("C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/output/", pattern = "\\.RData$"))

Pathway_recap <- lapply(seq_along(orthologs), function(i, input_disease_model = all_input_disease_model, 
                                               input_human = all_input_human){
  # Load human objects
  load(all_input_human[i])
  human <- object
  
  # Load disease model objects
  load(all_input_disease_model[i])
  disease_model <- output
  
  # Set output object to return
  
  # Select samples from base class and target class
  base_class <- which(object$Y.matrix == 0)
  target_class <- which(object$Y.matrix == 1)
  
  # Format human data block
  X.matrix.human <- human_data(human$X.matrix)

  # Output logFC with cell type/orthology status
  FC_oxa <- FC(X.matrix.human, disease_model$Oxazolone, model_name = "OXA", Pathway_name = disease_model$Pathway_name)  # Oxa
  FC_imq <- FC(X.matrix.human, disease_model$Imiquimod, model_name = "IMQ", Pathway_name = disease_model$Pathway_name)  # imq
  FC_ova <- FC(X.matrix.human, disease_model$Ovalbumine, model_name = "OVA", Pathway_name = disease_model$Pathway_name) # Ova

  # Generate In Silico treated human samples with disease model logFC
  InSilico_oxa <- InSilico(human$X.matrix[base_class, ], FC_oxa, human$asmbPLSDA.fit)
  InSilico_imq <- InSilico(human$X.matrix[base_class, ], FC_imq, human$asmbPLSDA.fit)
  InSilico_ova <- InSilico(human$X.matrix[base_class, ], FC_ova, human$asmbPLSDA.fit)
  
  ## RECAPITULATION CALCULATION
  # Pathway ideal recapitulation
  cell_type_new <- c("Dendritic Cells","Keratinocytes", "Langerhans Cells", "T-cell", "Melanocytes")
  pred_human <- response_recapitulation(human$X.matrix, 
                                          human$asmbPLSDA.fit, 
                                          type = "pathway",
                                          cell_type_new)
  
  Omega <- median(pred_human[target_class])-median(pred_human[base_class])
  
  # Disease models
  pred_oxa <- response_recapitulation(InSilico_oxa, 
                                      human$asmbPLSDA.fit, 
                                      type = "pathway",
                                      cell_type_new)
  pred_imq <- response_recapitulation(InSilico_imq, 
                                      human$asmbPLSDA.fit, 
                                      type = "pathway",
                                      cell_type_new)
  pred_ova <- response_recapitulation(InSilico_ova, 
                                      human$asmbPLSDA.fit, 
                                      type = "pathway",
                                      cell_type_new)
  pathway_recapitulation_oxa <- pathway_recapitulation(pred_oxa, pred_human, Omega)
  pathway_recapitulation_imq <- pathway_recapitulation(pred_imq, pred_human, Omega)
  pathway_recapitulation_ova <- pathway_recapitulation(pred_ova, pred_human, Omega)

  recapitulation_pathway_output <- data.frame("Pathway_name" = human$Pathway_name, "OXA" = pathway_recapitulation_oxa, 
                                      "IMQ" = pathway_recapitulation_imq, "OVA" = pathway_recapitulation_ova,
                                      "OXA_orth" = length(intersect(disease_model$Oxazolone$Orthologs_observed, human$Gene_set))/length(human$Gene_set), 
                                      "IMQ_orth" = length(intersect(disease_model$Imiquimod$Orthologs_observed, human$Gene_set))/length(human$Gene_set), 
                                      "OVA_orth" = length(intersect(disease_model$Ovalbumine$Orthologs_observed, human$Gene_set))/length(human$Gene_set)
  )
  
  ## CELL TYPE RECAPITULATION
  cell_pathway <- response_recapitulation(human$X.matrix, human$asmbPLSDA.fit, type = "cell", cell_type_new)
  gamma_pathway <- cell_ideal_recapitulation(cell_pathway)
  
  pred_oxa_cell <- response_recapitulation(InSilico_oxa, human$asmbPLSDA.fit, type = "cell", cell_type_new)
  pred_imq_cell <- response_recapitulation(InSilico_imq, human$asmbPLSDA.fit, type = "cell", cell_type_new)
  pred_ova_cell <- response_recapitulation(InSilico_ova, human$asmbPLSDA.fit, type = "cell", cell_type_new)
  
  cell_recapitulation_oxa <- cell_observed_recapitulation(cell_pathway, pred_oxa_cell, gamma_pathway)
  cell_recapitulation_imq <- cell_observed_recapitulation(cell_pathway, pred_imq_cell, gamma_pathway)
  cell_recapitulation_ova <- cell_observed_recapitulation(cell_pathway, pred_ova_cell, gamma_pathway)
  
  recapitulation_cell_output <- data.frame("Pathway_name" = rep(human$Pathway_name, length(cell_type_new)),
                                           "cell_type" = rownames(cell_recapitulation_imq),
                                           "OXA" = cell_recapitulation_oxa$gamma,
                                           "IMQ" = cell_recapitulation_imq$gamma,
                                           "OVA" = cell_recapitulation_ova$gamma) 
  
  ## GENE CONTRIBUTION TO CELL TYPE RECAPITULATION
  pred_human_gene <- response_recapitulation(human$X.matrix[base_class,], human$asmbPLSDA.fit, type = "gene", cell_type_new)
  
  # Compute gene contributions
  pred_oxa_gene <- response_recapitulation(InSilico_oxa, human$asmbPLSDA.fit, type = "gene", cell_type_new)
  pred_imq_gene <- response_recapitulation(InSilico_imq, human$asmbPLSDA.fit, type = "gene", cell_type_new)
  pred_ova_gene <- response_recapitulation(InSilico_ova, human$asmbPLSDA.fit, type = "gene", cell_type_new)
  
  # Compute difference in gene contribution
  contribution_gene_oxa <- gene_diff_contribution(pred_oxa_gene, pred_human_gene, gamma_pathway)
  contribution_gene_imq <- gene_diff_contribution(pred_imq_gene, pred_human_gene, gamma_pathway)
  contribution_gene_ova <- gene_diff_contribution(pred_ova_gene, pred_human_gene, gamma_pathway)
  
  contribution_gene_oxa$disease_model <- "OXA"
  contribution_gene_imq$disease_model <- "IMQ"
  contribution_gene_ova$disease_model <- "OVA"
  contribution_gene_output <- rbind(contribution_gene_oxa, contribution_gene_ova, contribution_gene_imq)
  contribution_gene_output$Pathway_name <- human$Pathway_name
  
  return(list("Pathway" = recapitulation_pathway_output, "Cell" = recapitulation_cell_output,
              "Gene" = contribution_gene_output, "FC" = list("OXA" = FC_oxa, "IMQ" = FC_imq, "OVA" = FC_ova)))
})
