asmbPLSDA.cv.loo <- function(X.matrix,
                             Y.matrix,
                             PLS_term = 1,
                             X.dim, 
                             quantile.comb.table, 
                             outcome.type = c("binary", "multiclass"),
                             Method = NULL,
                             measure = "B_accuracy",
                             parallel = FALSE,
                             cores = NULL,
                             expected.measure.increase = 0.005,
                             center = TRUE,
                             scale = TRUE,
                             maxiter = 100){
  
  # Libraries 
  library(asmbPLS)
  if(parallel){
    library(parallel)
    library(furrr)
  }
  
  n_group = ncol(Y.matrix)
  quantile_table_CV <- matrix(data = rep(0, PLS_term), nrow = PLS_term, 
                              ncol = (length(X.dim) + 5)) # Table to save the best quantile combination and the corresponding measures
  
  
  #if (outcome.type == "multiclass") {
  #  for (int n = 0; n < ncv; ++n) {
  #   CV_index_results[n] = CV_index_multiclass(F_matrix, K);
  #  }
  #}
  
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
  
  # Number of LOO samples
  K = nrow(Y.matrix)
  
  # Number of quantile combinations
  n_quantile_comb = nrow(quantile.comb.table)
  
  for (i in 1:PLS_term) {
    results_CV_summary_n = matrix(data = rep(0, n_quantile_comb), nrow = n_quantile_comb,
                                   ncol = K, byrow = TRUE)
    F_matrix_validation_bind = matrix(data = rep(0, n_quantile_comb), nrow = n_quantile_comb,
                                      ncol = K, byrow = TRUE)
    if(parallel){
      if(is.null(cores)){
        workers <- detectCores(logical=FALSE)-1
        # Start cluster
        plan(multisession, workers = workers)
        # Object to execute parallel on
        j <- data.frame("j" = 1:K)
        output <- furrr::future_pmap(j, quantile_computation, results_CV_summary_n = results_CV_summary_n, 
                                     F_matrix_validation_bind = F_matrix_validation_bind, X.matrix = X.matrix,
                                     Y.matrix = Y.matrix, PLS_term = i, X.dim = X.dim, 
                                     quantile.comb.table = quantile.comb.table, outcome.type = outcome.type,
                                     quantile_table_CV = quantile_table_CV, K = K, n_quantile_comb = n_quantile_comb,
                                     Method = Method, measure = measure, expected.measure.increase = expected.measure.increase,
                                     center = center, scale = center, maxiter = maxiter, .progress = TRUE, 
                                     .options = furrr_options(globals = FALSE, seed = TRUE))
        
        # Joint results from the lists
        results_CV_summary_n <- output[[1]]$results_CV_summary_n
        F_matrix_validation_bind <- output[[1]]$F_matrix_validation_bind
        for(ncol in 2:K){
          results_CV_summary_n[,ncol] <- output[[ncol]]$results_CV_summary_n[,ncol]
          F_matrix_validation_bind[,ncol] <- output[[ncol]]$F_matrix_validation_bind[,ncol]
        }
        
        # Stop the cluster
        plan(sequential)
      }else{
        workers <- cores
        # Start cluster
        plan(multisession, workers = workers)
        # Object to execute parallel on
        j <- data.frame("j" = 1:K)
        output <- furrr::future_pmap(j, quantile_computation, results_CV_summary_n = results_CV_summary_n, 
                                     F_matrix_validation_bind = F_matrix_validation_bind, X.matrix = X.matrix,
                                     Y.matrix = Y.matrix, PLS_term = i, X.dim = X.dim, 
                                     quantile.comb.table = quantile.comb.table, outcome.type = outcome.type,
                                     quantile_table_CV = quantile_table_CV, K = K, n_quantile_comb = n_quantile_comb,
                                     Method = Method, measure = measure, expected.measure.increase = expected.measure.increase,
                                     center = center, scale = center, maxiter = maxiter, .progress = TRUE, 
                                     .options = furrr_options(globals = FALSE, seed = TRUE))
        
        # Joint results from the lists
        results_CV_summary_n <- output[[1]]$results_CV_summary_n
        F_matrix_validation_bind <- output[[1]]$F_matrix_validation_bind
        for(ncol in 2:K){
          results_CV_summary_n[,ncol] <- output[[ncol]]$results_CV_summary_n[,ncol]
          F_matrix_validation_bind[,ncol] <- output[[ncol]]$F_matrix_validation_bind[,ncol]
        }
        # Stop the cluster
        plan(sequential)
      }
    }else{
      for(j in 1:K) {
        validation_index = j
        aux = 1:K
        training_index =  aux[aux != j]
        temp = 0
        
        # obtain validation and training sets
        E_matrix_validation = as.matrix(t(X.matrix[validation_index,]))
        F_matrix_validation = as.matrix(t(Y.matrix[validation_index,]))
        E_matrix_training = as.matrix(X.matrix[training_index,])
        F_matrix_training = as.matrix(Y.matrix[training_index,])
        
        # calculate overall/balanced accuracy using different quantile combinations
        for (l in 1:n_quantile_comb) {
          quantile_table_CV[i, 1:length(X.dim)] = quantile.comb.table[l, 1:length(X.dim)]
          if(i == 1){
            quantile_temp = t(as.matrix(quantile_table_CV[1:i, 1:length(X.dim)]))
          }else{
            quantile_temp = as.matrix(quantile_table_CV[1:i, 1:length(X.dim)])
          }
          # fit model using training set
          asmbPLSDA_fit_results = asmbPLSDA.fit(E_matrix_training, F_matrix_training, i, X.dim, quantile_temp, outcome.type, center, scale, maxiter)
          asmbPLSDA_predict_results = asmbPLSDA.predict(asmbPLSDA_fit_results, E_matrix_validation, i, Method)
          Y_pred = as.numeric(asmbPLSDA_predict_results["Y_pred"])
          results_CV_summary_n[l, j] = Y_pred
          F_matrix_validation_bind[l, j] = F_matrix_validation
        }
      }
    }
    
    # calculate the mean accuracy for each quantile combination
    measure_acc = c()
    for(l in 1:n_quantile_comb){
      Y_pred_ = as.vector(results_CV_summary_n[l,])
      F_matrix_validation_ = as.vector(F_matrix_validation_bind[l,])
      measure_new = Results_comparison_measure(Y_pred_, F_matrix_validation_, outcome.type)
      measure_acc = c(measure_acc, measure_new[measure_selected])
    }
    
    # find the quantile combination with the highest accuracy
    index_max_measure = which.max(measure_acc)
    quantile_table_CV[i, 1:length(X.dim)] = quantile.comb.table[index_max_measure,]
    
    Y_pred_bind = matrix()
    F_matrix_validation_bind = matrix()
    
    # obtain corresponding measures for the selected quantile combination
      for (j in 1:K) {
        validation_index = j
        aux = 1:K
        training_index =  aux[aux != j]
        temp = 0
        
        # obtain validation and training sets
        E_matrix_validation = as.matrix(t(X.matrix[validation_index,]))
        F_matrix_validation = as.matrix(t(Y.matrix[validation_index,]))
        E_matrix_training = as.matrix(X.matrix[training_index,])
        F_matrix_training = as.matrix(Y.matrix[training_index,])
        
        if(i == 1){
          quantile_temp = t(as.matrix(quantile_table_CV[1:i, 1:length(X.dim)]))
        }else{
          quantile_temp = as.matrix(quantile_table_CV[1:i, 1:length(X.dim)])
        }
        # fit model using training set
        asmbPLSDA_fit_results = asmbPLSDA.fit(E_matrix_training, F_matrix_training, i, X.dim, quantile_temp, outcome.type, center, scale, maxiter)
        asmbPLSDA_predict_results = asmbPLSDA.predict(asmbPLSDA_fit_results, E_matrix_validation, i, Method)
        Y_pred = as.matrix(asmbPLSDA_predict_results["Y_pred"])
        rownames(Y_pred) = NULL
        Y_pred_bind = rbind(Y_pred_bind, Y_pred)
        F_matrix_validation_bind = rbind(F_matrix_validation_bind, F_matrix_validation)
      }
    measure = Results_comparison_measure(Y_pred_bind, F_matrix_validation_bind, outcome.type)
    quantile_table_CV[i, (length(X.dim)+1):ncol(quantile_table_CV)] = measure
    #Results_comparison_measure(Y_pred, F_matrix_validation, outcome.type)
  }
  
  optimal_nPLS = 1
  if(PLS_term > 1) {
    for (i in 1:(PLS_term-1)) {
      current_measure =  quantile_table_CV[i, (length(X.dim) + measure_selected)]
      next_measure =  quantile_table_CV[i+1, (length(X.dim) + measure_selected)]
      if(next_measure > current_measure + expected.measure.increase) {
        optimal_nPLS = optimal_nPLS + 1
      } else {break}
    }
  }
  
  output = list("quantile_table_CV"= quantile_table_CV,
                "optimal_nPLS" = optimal_nPLS)
  
  return(output)
}

