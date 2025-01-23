  # Description: calculates in silico treated human samples with FC disease model. This is the biological link function.
  InSilico <- function(human, logFC, asmbpls){
    rownames(logFC) <- colnames(human)
    aux_logFC <- t(logFC[1])
    case <- logFC[2]
    X_prime <- matrix(ncol = ncol(human), nrow = nrow(human))
    # Estimated average gene expression by asmbPLS
    mean <- asmbpls$X_col_mean 
    for(i in 1:ncol(human)){
      # Case with no ortholog and/or no cell type present
      if(case[i,] <= 0){
        X_prime[1:nrow(human), i] <- human[1:nrow(human), i]
        # Case with ortholog, cell type present and not negative logFC
      }else if(case[i,] >= 1 && min(human[, i]+mean[i]*aux_logFC[i]) >= 0){ 
        X_prime[1:nrow(human), i] <- human[1:nrow(human), i] + mean[i]*aux_logFC[i] 
      }else if(case[i,] >= 1 && min(human[, i]+mean[i]*aux_logFC[i]) < 0){ 
        index <- which.min(human[, i])
        X_prime[1:nrow(human), i] <- human[1:nrow(human), i] - human[index, i]
      }
    }
    colnames(X_prime) <- colnames(human)
    return(X_prime)
  }
