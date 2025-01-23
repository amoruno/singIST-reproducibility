 # Description: compute recapitulations in the response space from the superpathway, cell type and gene contributions 
  response_recapitulation <- function(data, fit_asmb, type = c("pathway", "cell", "gene"), 
                                      cell_type_new){
    response_recapitulation <- c()
    for(i in 1:nrow(data)){
      output <- cell_gene_latent_recapitulations(data, i, fit_asmb, cell_type_new, 
                                                 center = TRUE, scale = TRUE) # change to center = fit_asmb$center and scale = fit_asmb$scale
      # Response pathway recapitulation
      if(type == "pathway"){
        response_recapitulation <- cbind(response_recapitulation, output$Y_pred_num)
      }
      else if(type == "cell"){
        # Response cell recapitulation
        response_recapitulation <- cbind(response_recapitulation, output$gamma %*% as.matrix(t(output$Y_weight)))
      }
      else if(type == "gene"){
        # Response gene recapitulation
        j_num <- 1
        for(j in rownames(output$gamma)){
          index <- which(sub("_.*", "", colnames(output$delta)) == j)
          w_super <- as.matrix(fit_asmb$X_super_weight[j_num,])
          Y_weight <- t(as.matrix(output$Y_weight))
          aux <- as.matrix(w_super*Y_weight)
          delta <- as.matrix(output$delta[, index]) %*% aux
          colnames(delta) <- paste0(j, "_", i)
          response_recapitulation <- cbind(response_recapitulation, delta)
          j_num <- j_num + 1
        }
      }
    }
    return(response_recapitulation)
  }
