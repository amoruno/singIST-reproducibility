# Description: compute Fold Change for target class vs. base class of a disease model given the logFC from FindMarkers
 # Calculates log Fold Change between target class and base of disease model
  FC <- function(human, disease_model, model_name = NULL, Pathway_name = NULL){
    log_FC <- numeric(nrow(human))
    case <- numeric(nrow(human)) # case = 0 not ortholog, case = -1 no cell type, case = 1 not significant, 
                                 # case = 2 significant
    rownames_output <- c()
    # epsilon <- .Machine$double.eps
    cell_disease_model <- unique(sub("\\*.*", "", colnames(disease_model$X.matrix)))
    gene_disease_model <- unique(sub(".*\\*", "", colnames(disease_model$X.matrix)))
    logFC_FindMarkers <- disease_model$logFC
    for(i in 1:nrow(human)){
      cell <- sub("\\*.*", "", rownames(human)[i])
      gene <- sub(".*\\*", "", rownames(human)[i])
      rownames_output[i] <- paste0(cell, "*", gene)
      # Cell-type exists
      if(cell %in% cell_disease_model){
        # There exists a one-to-one ortholog gene
        if(gene %in% gene_disease_model){
          element <- paste0(cell, "*", gene)
          index <- which(rownames(logFC_FindMarkers) == element)
          if(length(index) == 1){ # Look if gene is in logFC from FindMarkers 
            if(logFC_FindMarkers[index, "p_val_adj"] <= 0.05){ # significant
              log_FC[i] <- logFC_FindMarkers[index, "avg_log2FC"]
              case[i] <- 2
            }else{ # not significant
              log_FC[i] <- 0
              case[i] <- 1
            }
          }else{ # If not in logFC_FindMarkers but ortholog then it is non significant
            log_FC[i] <- 0
            case[i] <- 1
          }
        }else{
          log_FC[i] <- NaN # Not ortholog set to NaN
          case[i] <- -1
        }
      }else{
        # Cell-type non-existent
        log_FC[i] <- 0
        case[i] <- 0
      }
    }
    output <- data.frame(logFC = log_FC, InSilicoSum = case)
    rownames(output) <- rownames_output
    # Save genes and cell type onto variable
    output$cell_type <- sub("\\*.*", "", rownames_output)
    output$gene_name <- sub(".*\\*", "", rownames_output)
    # Add pathway name and disease model
    output$Pathway_name <- Pathway_name
    output$Disease_model <- model_name
    return(output)
  }
