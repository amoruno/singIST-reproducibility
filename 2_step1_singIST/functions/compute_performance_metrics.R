# Function to compute metrics
  Results_comparison_measure <- function(Y_predict, 
                                         Y_true, 
                                         outcome.type = c("binary", "multiclass")){
    Y_col = length(as.vector(Y_true)) # Number of groups of samples
    n_match = 0
    n_TP = 0 # Number of true positive
    n_TN = 0 # Number of true negative
    n_FP = 0 # Number of false positive
    n_FN = 0 #Number of false negative
    balanced_accuracy_multicalss = 0
    n_recall_multiclass = 0
    
    for (i in 1:Y_col) {
      temp_accu = which(as.vector(Y_predict)[i] == as.vector(Y_true)[i])
      temp_TP = which(as.vector(Y_predict)[i] == 1 && as.vector(Y_true)[i] == 1)
      temp_TN = which(as.vector(Y_predict)[i] == 0 && as.vector(Y_true)[i] == 0)
      temp_FP = which(as.vector(Y_predict)[i] == 1 && as.vector(Y_true)[i] == 0)
      temp_FN = which(as.vector(Y_predict)[i] == 0 && as.vector(Y_true)[i] == 1)
      
      n_temp_TP = length(temp_TP)
      n_temp_TN = length(temp_TN)
      n_temp_FP = length(temp_FP)
      n_temp_FN = length(temp_FN)
      
      n_match = n_match + temp_accu
      n_TP = n_TP + n_temp_TP
      n_TN = n_TN + n_temp_TN
      n_FP = n_FP + n_temp_FP
      n_FN = n_FN + n_temp_FN
      
      n_recall_multiclass = n_temp_TP/(n_temp_TP + n_temp_FN)
      balanced_accuracy_multicalss = balanced_accuracy_multicalss + n_recall_multiclass
    }
    accuracy = (n_TP + n_TN)/(n_TP + n_TN + n_FP + n_FN)
    precision = n_TP/(n_TP + n_FP)
    recall = n_TP/(n_TP + n_FN) # Sensitivity also
    specificity = n_TN/(n_TN + n_FP)
    F1 = 2 * precision * recall/(precision + recall)
    
    balanced_accuracy = 0
    balanced_accuracy_binary = (recall + specificity)/2
    balanced_accuracy_multicalss = balanced_accuracy_multicalss/Y_col
    if (outcome.type == "binary") {
      balanced_accuracy = balanced_accuracy_binary
    }
    if (outcome.type == "multiclass") {
      balanced_accuracy = balanced_accuracy_multicalss
    }
    
    output <- c("accuracy" = accuracy, 
                "balanced_accuracy" = balanced_accuracy,
                "precision" = precision,
                "recall" = recall,
                "F1" = F1)
    return(output)
  }
