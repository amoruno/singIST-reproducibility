pval <- function(val){
  # If p-val < 0.01
  if(val <= 0.001){
    return(paste0("p <= 0.001"))
  }else if(0.001 < val & val <= 0.01){
    return(paste0("p <= 0.01"))
  }else if(0.01 < val & val <= 0.05){
    return(paste0("p <= 0.05"))
  }else if(val > 0.05 & val < 1){
    return(paste0("p = ", round(val, 3)))
  }else if(val >= 1){
    return(paste0("p = 1"))
  }
}
