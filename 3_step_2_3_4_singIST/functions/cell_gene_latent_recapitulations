# Description: compute gene contributions, cell type contributions and superpathway contributions from asmbPLS-DA scores and weights
  cell_gene_latent_recapitulations <- function(data, sample, fit_asmb, cell_type_new, center = TRUE, scale = TRUE){
    delta_cbind <- c()
    Delta <- c()
    Y_fit <- 0
    gamma <- c()
    Gamma <- c()
    n.PLS.binary <- ncol(fit_asmb$Y_weight)
    FC_applied_test <- as.matrix(t(data[sample, colnames(data)]))
    
    X.dim <- fit_asmb$X_dim
    Y_col_mean <- fit_asmb$Y_col_mean
    X_col_mean <- fit_asmb$X_col_mean
    X_col_sd <- fit_asmb$X_col_sd
    
    # Newdata scale 
    if(center){
      if(scale){
        # center = TRUE and scale = TRUE
        for(i in 1:sum(X.dim)){
          FC_applied_test[, i] <- (FC_applied_test[, i] - X_col_mean[i])/X_col_sd[i]
        }
      }else{
        # center = TRUE and scale = FALSE
        for(i in 1:sum(X.dim)){
          FC_applied_test[, i] <- FC_applied_test[, i] - X_col_mean[i]
        }
      }
    }else{
      if(scale){
        # center = FALSE and scale = CENTER
        for(i in 1:sum(X.dim)){
          FC_applied_test[, i] <- FC_applied_test[i, ]/X_col_sd[i]
        }
      }
    } # center = FALSE and scale = TRUE, no action 
    
    for(i in 1:n.PLS.binary){
      for(j in 1:length(X.dim)){
        C_b <- as.matrix(FC_applied_test[, (1+X.dim[j]*(j-1)):(X.dim[j]*j)])
        value_aux <- unique(sub("\\*.*", "", rownames(C_b)))
        index_aux <- which(cell_type_new == value_aux) # Pick the index of cell value in C_b
        omega <- as.matrix(fit_asmb$X_weight[[j]][,i])
        delta <- C_b*omega/(sqrt(X.dim[j]))
        rownames(delta) <- sub(".*\\*", "", rownames(delta))
        colnames(delta) <- paste0(cell_type_new[index_aux], "_",i)
        delta_cbind <- cbind(delta_cbind, delta)
      }
      Delta <- cbind(Delta, colSums(delta_cbind[, (1+length(X.dim)*(i-1)):(length(X.dim)*i)]))
      Delta <- as.matrix(Delta)
      rownames(Delta) <- sub("\\_.*", "", rownames(Delta)) 
      omega_super <- as.matrix(fit_asmb$X_super_weight[,i])
      aux <- cbind(Delta[,i], omega_super)
      gamma <- cbind(gamma, aux[,1]*aux[,2])
      Gamma <- cbind(Gamma, sum(gamma[,i]))
      Gamma <- as.matrix(Gamma)
      #Deflation of block data
      for(j in 1:length(X.dim)){
        C_b <- as.matrix(FC_applied_test[, (1+X.dim[j]*(j-1)):(X.dim[j]*j)])
        # Deflation by block loadings 
        p_temp <- as.matrix(fit_asmb$X_loading[[j]][,i])
        p_aux <- as.matrix(delta_cbind[, (1+length(X.dim)*(i-1)):(length(X.dim)*i)])
        p <- as.matrix(p_aux[, j])
        FC_applied_test[, (1+X.dim[j]*(j-1)):(X.dim[j]*j)] <- C_b - p*p_temp
      }
      # Response space recapitulation
      q <- fit_asmb$Y_weight[,i]
      Y_fit <- Y_fit + Gamma[,i]*q
    }
    
    # center = TRUE, add the mean of Y back
    if(center){
      for(i in 1:length(Y_col_mean)){
        Y_fit <-  Y_fit + Y_col_mean[i]
      }
    }
    
    output <- list(Delta = Delta, delta = delta_cbind, gamma = gamma, Gamma = Gamma,
                   Y_pred_num = Y_fit, Y_weight = fit_asmb$Y_weight)
    return(output)
  } 
