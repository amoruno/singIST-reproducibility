  # Wilcox test function to return pvalue
  wilcox.CIP.GIP <- function(jack, null){
    return(wilcox.test(jack, null, alternative = "greater")$p.value)
  }
