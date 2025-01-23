# Description: compute cell type observed recapitulation
cell_observed_recapitulation <- function(cell_response_recapitulation, observed, ideal_recapitulation){
    output <- data.frame()
    for(i in 1:nrow(cell_response_recapitulation)){
      aux <- 100*(median(observed[i, 1:4])-median(cell_response_recapitulation[i, 6:9]))/ideal_recapitulation[i,]
      output <- rbind(output, aux)
    }
    rownames(output) <- rownames(cell_response_recapitulation)
    colnames(output) <- "gamma"
    return(output)
  }
