jackk_asmbplsda <- function(object, quantile.comb, 
                            n.PLS, X.dim, Y.matrix, X.matrix, center = TRUE,
                            scale = TRUE, outcome.type = c("binary", "multiclass"),
                            expected.measure.increase = 0.005, maxiter = 100, npermut = 100){
  # Number of samples
  K = nrow(Y.matrix)
  # Initialize jackknife parameter object
  jackk_param <- list()
  VAR_INF <- list()
  
  # Estimate optimal model by jackknife
  for(j in 1:K) {
    aux = 1:K
    training_index =  aux[aux != j]
    
    # obtain validation and training sets
    E_matrix_training = as.matrix(X.matrix[training_index,])
    F_matrix_training = as.matrix(Y.matrix[training_index,])
    
    # fit model using training set
    asmbPLSDA_fit_results = asmbPLSDA.fit(E_matrix_training, F_matrix_training, n.PLS, 
                                          X.dim, quantile.comb, outcome.type, center, 
                                          scale, maxiter)
    
    # Compute CIP/GIP 
    aux <- CIP_GIP(asmbPLSDA_fit_results)
    if(j == 1){
      jackk_param <- aux
    }else{
      for(b in 1:length(X.dim)){
        jackk_param$GIP[[b]] <- cbind(jackk_param$GIP[[b]], aux$GIP[[b]])
      }
      jackk_param$CIP <- cbind(jackk_param$CIP, aux$CIP)
    }
  }
  
  # Generate null CIP and GIP distributions
  for (j in 1:(npermut + 1)) {
    set.seed(seed = j)
    X_perm_aux <- matrix(nrow = K, ncol =  ncol(X.matrix))
    colnames(X_perm_aux) <- colnames(X.matrix)
    rownames(X_perm_aux) <- rownames(X.matrix)
    for(b in 1:length(X.dim)){
      X_perm_aux[1:K, (1+X.dim[b]*(b-1)):(X.dim[b]*b)] <- X.matrix[sample(K), sample((1+X.dim[b]*(b-1)):(X.dim[b]*b))] 
                                                                          
    }
    colnames(X_perm_aux) <- colnames(X.matrix)
    
    # Fit optimal model with permutated blocks 
    Modelpermut <- asmbPLSDA.fit(X.matrix = X_perm_aux, 
                                 Y.matrix = Y.matrix, 
                                 PLS.comp = n.PLS, 
                                 X.dim = X.dim, 
                                 quantile.comb = quantile.comb,
                                 center = center,
                                 scale = scale,
                                 outcome.type = outcome.type)
    if(j == 1){
      NULL_VAR_INF <- CIP_GIP(Modelpermut)
    }else{
      aux <- CIP_GIP(Modelpermut)
      for(b in 1:length(X.dim)){
        NULL_VAR_INF$GIP[[b]] <- cbind(NULL_VAR_INF$GIP[[b]], aux$GIP[[b]])
      }
      NULL_VAR_INF$CIP <- cbind(NULL_VAR_INF$CIP, aux$CIP)
    }
  }
  
  # Compute Mann-Whitney U test between H0 and Jackknife distribution
  # Compute test for each gene within cell
  GIP_pvalue  <- lapply(seq_along(jackk_param$GIP), function(i){
    jack_cell <- jackk_param$GIP[[i]]
    null_cell <- NULL_VAR_INF$GIP[[i]]
    tests <- mapply(function(jack_row, null_row) {
                      wilcox.CIP.GIP(jack_row, null_row)
                    }, split(jack_cell, row(jack_cell)), split(null_cell, row(null_cell)), SIMPLIFY = TRUE)
    output <- data.frame("pvalue" = tests)
    rownames(output) <- rownames(jack_cell)
    return(output)
  })
  
  CIP_pvalue <- lapply(seq_along(jackk_param$GIP), function(i){ # Number of blocks = seq_along(jackk_param$GIP)
    jack_cell_CIP <- jackk_param$CIP[i,]
    null_cell_CIP <- NULL_VAR_INF$CIP[i,]
    tests <- wilcox.CIP.GIP(jack_cell_CIP, null_cell_CIP)
    output <- tests
    return(output)
  })

  return(list("jack_param" = jackk_param, "NULL_CIP_GIP" = NULL_VAR_INF, "CIP_pvalue" = CIP_pvalue,"GIP_pvalue" = GIP_pvalue))
}
