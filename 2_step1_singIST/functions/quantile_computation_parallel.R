 # Function to parallelize model fit for each quantile
    quantile_computation <- 
      function(j, ..., results_CV_summary_n, F_matrix_validation_bind, 
               X.matrix, Y.matrix, PLS_term = 1, X.dim, 
               quantile.comb.table, outcome.type = c("binary", "multiclass"), quantile_table_CV, K, n_quantile_comb,
               Method = NULL, measure = "B_accuracy", expected.measure.increase = 0.005,
               center = TRUE, scale = TRUE, maxiter = 100){
        suppressWarnings(library(asmbPLS))
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
          quantile_table_CV[PLS_term, 1:length(X.dim)] <- quantile.comb.table[l, 1:length(X.dim)]
          if(PLS_term == 1){
            quantile_temp = t(as.matrix(quantile_table_CV[1:PLS_term, 1:length(X.dim)]))
          }else{
            quantile_temp = as.matrix(quantile_table_CV[1:PLS_term, 1:length(X.dim)])
          }
          # fit model using training set
          asmbPLSDA_fit_results = asmbPLSDA.fit(E_matrix_training, F_matrix_training, PLS_term, X.dim, quantile_temp, outcome.type, center, scale, maxiter)
          asmbPLSDA_predict_results = asmbPLSDA.predict(asmbPLSDA_fit_results, E_matrix_validation, PLS_term, Method)
          Y_pred = as.numeric(asmbPLSDA_predict_results["Y_pred"])
          results_CV_summary_n[l, j] <- Y_pred
          F_matrix_validation_bind[l, j] <- F_matrix_validation                              
        }
        return(list("results_CV_summary_n" = results_CV_summary_n,
                    "F_matrix_validation_bind" = F_matrix_validation_bind,
                    "obs" = j)
        )
      }
