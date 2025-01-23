# Description: cell ideal and observed recapitulation
  cell_ideal_recapitulation <- function(cell_response_recapitulation){
    output <- data.frame()
    for(i in 1:nrow(cell_response_recapitulation)){
      aux <- median(cell_response_recapitulation[i, 1:4])-median(cell_response_recapitulation[i, 6:9])
      output <- rbind(output, aux)
    }
    rownames(output) <- rownames(cell_response_recapitulation)
    colnames(output) <- "gamma"
    return(output)
  }
