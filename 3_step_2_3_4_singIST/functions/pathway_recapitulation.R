  # Description: superpathway recapitulation function
  pathway_recapitulation <- function(pred_scist, pred_human, omega_ideal){
    # pred_scist contains the In Silico Treated healthy samples 1:4
    # while pred_human contains the true bsae class
    aux <- median(pred_scist[1:4])-median(pred_human[6:9])
    return(100*(aux/omega_ideal))
  }
