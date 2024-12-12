# Description: Function to perform permutation testing with 2-fold
# cross-validation for adaptive sparse multi-block partial least squares
# discriminant analysis, in order to evaluate model validity

# object: asmbPLSDA object to be assessed
#
permut_asmbplsda <- function(object, quantile.comb, n.PLS, Y.matrix, X.matrix, 
                             npermut = 100, nbObsPermut = NULL, Nc = 1, 
                             CV_error, measure = "B_accuracy"){
  library(FactoMineR)
  # Target class
  q <- ncol(object$Y_group)
  # Number samples
  nr <- nrow(object$Y_group)
  # Model options
  center <- object$center
  scale <- object$scale
  outcome.type <- object$Outcome_type
  X.dim <- object$X_dim
  # Selected validation measure
  if(measure == "accuracy") {
    measure_selected = 1
  }
  if(measure == "B_accuracy") {
    measure_selected = 2
  }
  if(measure == "precision") {
    measure_selected = 3
  }
  if(measure == "recall") {
    measure_selected = 4
  }
  if(measure == "F1") {
    measure_selected = 5
  }
  
  # GIP of optimal model
  VAR_INF <- CIP_GIP(object)
  # Lists that will store the permutations
  dimlabP <- c("NoPermut", paste("permut", (1:npermut), sep = ""))
  res <- list()
  res$RV.YYpermut.values <- data.frame(dimlabP, rep(NA, (npermut + 
                                                           1)), stringsAsFactors = TRUE)
  colnames(res$RV.YYpermut.values) <- c("dimlabP", "RV.YYpermut.values")
  res$cor.YYpermut.values <- data.frame(dimlabP, matrix(NA, 
                                                        ncol = q, nrow = (npermut + 1)), stringsAsFactors = TRUE)
  colnames(res$cor.YYpermut.values) <- c("dimlabP", "target class")
  res$prctGlob.Ychange.values <- data.frame(dimlabP, rep(NA, 
                                                         (npermut + 1)), stringsAsFactors = TRUE)
  colnames(res$prctGlob.Ychange.values) <- c("dimlabP", "prctGlob.Ychange.values")
  res$prct.Ychange.values <- data.frame(dimlabP, matrix(NA, 
                                                        ncol = q+1, nrow = (npermut + 1)), stringsAsFactors = TRUE)
  colnames(res$prct.Ychange.values) <- c("dimlabP", "target", measure)
  
  # Compute LOOCV F1 distribution
  
  ######################################
  # EVALUATION METRICS TO BE SELECTED ##
  ######################################
  for (j in 1:(npermut + 1)) {
    set.seed(seed = j)
    nObsP <- NA
    if (j == 1) {
      nObsP <- 0
    }
    if (j > 1) {
      if (is.null(nbObsPermut) == TRUE) {
        nObsP <- sample(x = nr, size = 1, replace = FALSE)
      }
      else {
        nObsP <- nbObsPermut
      }
    }
    Ypermut <- Y.matrix
    if ((nObsP > 0) == TRUE) {
      for (o in 1:nObsP) {
        # Permute response
        indObsPermut_Y <- sample(x = 1:nr, size = 2, replace = FALSE)
        Y1 <- Ypermut[indObsPermut_Y[1], ]
        Y2 <- Ypermut[indObsPermut_Y[2], ]
        Ypermut[indObsPermut_Y[1], ] <- Y2
        Ypermut[indObsPermut_Y[2], ] <- Y1
      }
    }
    rownames(Ypermut) <- rownames(Y.matrix)
    
    # Save permutation statistics for Y block
    res$cor.YYpermut.values[j, 2:(q + 1)] <- sapply(1:q, 
                                                    function(Q) cor(Y.matrix[, Q], Ypermut[, Q]))
    res$prct.Ychange.values[j, 2:(q + 1)] <- sapply(1:q, 
                                                    function(Q) (sum(Y.matrix[, Q] != Ypermut[, Q])/nr))
    res$RV.YYpermut.values[j, 2] <- coeffRV(Y.matrix, Ypermut)$rv
    res$prctGlob.Ychange.values[j, 2] <- sum((Y.matrix - Ypermut)^2)/(nr * q)
    
    # Save permutation statistics for X block
    
    # Number of CV
    s <- sample(x = nr, size = Nc)
    
    # Train and validation sets of permutated blocks
    # Tranpose X block if Nc = 1
    #Xpermut_val <- t(as.matrix(Xpermut[s, ]))
    #Xpermut_train <- as.matrix(Xpermut[-s, ])
    X_val <- t(as.matrix(X.matrix[s, ]))
    X_train <- as.matrix(X.matrix[-s, ])
    Ypermut_val <- t(as.matrix(Ypermut[s, ]))
    Ypermut_train <- as.matrix(Ypermut[-s, ])
    
    # Fit optimal model with permutated blocks 
    Modelpermut <- asmbPLSDA.fit(X.matrix = X_train, 
                                     Y.matrix = Ypermut_train, 
                                     PLS.comp = n.PLS, 
                                     X.dim = X.dim, 
                                     quantile.comb = quantile.comb,
                                     center = center,
                                     scale = scale,
                                     outcome.type = outcome.type)
    # Compute CIP and GIP
    aux <- CIP_GIP(Modelpermut)
    
    for(b in 1:length(X.dim)){
      VAR_INF$GIP[[b]] <- cbind(VAR_INF$GIP[[b]], aux$GIP[[b]])
    }
    VAR_INF$CIP <- cbind(VAR_INF$CIP, aux$CIP)
    # Predict response with permuted model and validation permutation
    ModelPredict <- asmbPLSDA.predict(Modelpermut, 
                                      X_train, 
                                      PLS.comp = n.PLS)$Y_pred
    ModelPredict_Val <- asmbPLSDA.predict(Modelpermut, 
                                          X_val, 
                                          PLS.comp = n.PLS)$Y_pred
    Yperm_pred <- append(ModelPredict, ModelPredict_Val, after=s-1)
    measures <- Results_comparison_measure(Yperm_pred, Y.matrix, 
                                           outcome.type = outcome.type)
    selected_measure <- measures[measure_selected]
    res$prct.Ychange.values[j, 2:(q + 2)] <- selected_measure
  }
  # IC95 of null distribution
  IC95 <- function(m) {
    rt <- c(round(quantile(m, 0.05), 5), round(quantile(m, 0.95), 5))
    return(rt)
  }
  # Wilcox-test to differentiate between null and CV error distribution 
  null_errors <- as.vector(res$prct.Ychange.values[, measure])
  if(measure == "F1"){
    # If F1 is NaN => Recall + Precision = 0, the model performance is poor impute to 0
    null_errors[is.nan(null_errors)] <- 0 
  }
  ecdf_errors <- ecdf(null_errors)
  res$IC <- IC95(null_errors)
  res$pvalue <- 1-ecdf_errors(CV_error)
  # Add CIP/GIP null distribution to output
  res$VAR_INF <- VAR_INF
  return(res)
}
