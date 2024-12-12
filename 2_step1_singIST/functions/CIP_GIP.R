 # Function to compute CIP and GIP metrics
  CIP_GIP <- function(object){
    # Assuming binary class
    # Loadings
    w_super <- object$X_super_weight^2
    q <- object$Y_weight
    nblocks <- nrow(w_super)
    
    # Compute CIP
    CIP <- (w_super %*% t(q))/sum(q)
    
    # Compute GIP
    GIP <- vector("list", nblocks)
    for(b in 1:nblocks){
      w <- object$X_weight[[b]]^2
      GIP[[b]] <- (w %*% t(q))/sum(q)
    }
    return(list("GIP" = GIP, "CIP" = CIP))
  }
