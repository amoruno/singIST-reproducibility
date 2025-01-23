# Description: format human data blocks 
  human_data <- function(data){
    target_class <- "AD" # AD
    X.matrix.new <- as.data.frame(data)
    class <- c(rep("AD", 5), rep("HC", 4)) # Number of observations
    X.matrix.new <- cbind(class, X.matrix.new)
    target_sum <- aggregate(. ~  class, X.matrix.new[X.matrix.new$class == target_class,], sum)
    rownames(target_sum) <- target_sum$class
    target_sum <- as.data.frame(t(target_sum[, 2:ncol(target_sum)])) 
    return(target_sum)
  }
