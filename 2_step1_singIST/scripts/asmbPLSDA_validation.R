# libraries
library(BiocGenerics)
library(BiocFileCache)
library(AnnotationHub)
library(ExperimentHub)
library(msigdb)

# Store global significance test pvalues
pvalues <- c()
for(i in all_output){ # change to all_output
  load(paste0("C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/output", "/", i)) # change to i
  # Pathway name
  summary_table[[i]]$Pathway_name <- object$Pathway_name
  # Optimal number of PLS
  summary_table[[i]]$Optimal_nPLS <- object$asmbPLSDA.cv$optimal_nPLS
  # Gene set size
  summary_table[[i]]$Gene_set_size <- sum(object$Gene_set %in%
                                            sub(".*\\*(.*)", "\\1", rownames(object$asmbPLSDA.fit$X_weight[[1]])))
  # Cell type lambda
  summary_table[[i]]$lambdas <- matrix(object$asmbPLSDA.cv$quantile_table_CV[1:object$asmbPLSDA.cv$optimal_nPLS, 1:length(object$X.dim)], nrow = object$asmbPLSDA.cv$optimal_nPLS, ncol = length(object$X.dim))
  colnames(summary_table[[i]]$lambdas) <- c("DC", "KC", "LC", "MC", "TC")
  
  # Quantitative validation of optimal model
  # Permutation test of optimal value
  model_validation <- 
    lapply(1, function(i, npermut = 10000){
      CV_error <- object$asmbPLSDA.cv$quantile_table_CV[object$asmbPLSDA.cv$optimal_nPLS, 10]
      model_statistics <-  permut_asmbplsda(object$asmbPLSDA.fit, object$asmbPLSDA.cv$quantile_table_CV, 
                                            object$asmbPLSDA.cv$optimal_nPLS, object$Y.matrix, object$X.matrix, 
                                            npermut = npermut, nbObsPermut = 9, Nc = 1, 
                                            CV_error = CV_error, measure = "F1")
      return(model_statistics)}
    )
  
  # Save pvalue for correcting later on
  pvalues <- c(pvalues, model_validation[[1]]$pvalue)

  # Mann-Whitney U test: null distribution GIP vs Jackknife
  npermut <- 10000
  CIP_GIP_mannwhitney <- jackk_asmbplsda(object$asmbPLSDA.fit.binary, 
                                     object$asmbPLSDA.cv$quantile_table_CV, 
                                     object$asmbPLSDA.cv$optimal_nPLS, object$X.dim, 
                                     object$Y.matrix, object$X.matrix, center = object$asmbPLSDA.fit$center,
                                     scale = object$asmbPLSDA.fit$scale, 
                                     outcome.type = object$asmbPLSDA.fit$Outcome_type,
                                     expected.measure.increase = 0.005, maxiter = 100, npermut = npermut)
  
  # Bonferroni correction with m0 = floor(lambda*number genes in cell)
  CIP_GIP_adj_pval <- lapply(seq_along(CIP_GIP_mannwhitney$GIP_pvalue), function(j, pathway = summary_table[[i]]){
    # Correction factor
    lambda <- apply(pathway$lambdas, 2, function(x) prod(x, na.rm = TRUE))[j]
    m_cb <- pathway$Gene_set_size 
    m_0 <- ifelse(floor(lambda*m_cb) == 0, 1, floor(lambda*m_cb)) # Avoid case where floor = 0, at least
                                                                  # we consider 1 true null hypothesis
    # Adjust pvalue
    pval <- CIP_GIP_mannwhitney$GIP_pvalue[[j]][,1]
    adj_pval <- sapply(pval*m_0, pval)
    output <- data.frame("adj_pval" = adj_pval)
    rownames(output) <- rownames(CIP_GIP_mannwhitney$pvalue[[j]])
    return(output)
  })
  summary_table[[i]]$GIP_adj_pval <- CIP_GIP_adj_pval
  # CIP pvalue
  summary_table[[i]]$CIP_pvalue <- CIP_GIP_mannwhitney$CIP_pvalue
  # Count number of significant GIP within cell type
  DC <- sum(summary_table[[i]]$GIP_adj_pval[[1]][,1] %in% 
                                                          c("p <= 0.001", "p <= 0.01", "p <= 0.05"))
  KC <- sum(summary_table[[i]]$GIP_adj_pval[[2]][,1] %in% 
                                                          c("p <= 0.001", "p <= 0.01", "p <= 0.05"))
  LC <- sum(summary_table[[i]]$GIP_adj_pval[[3]][,1] %in% 
                                                          c("p <= 0.001", "p <= 0.01", "p <= 0.05"))
  MC <- sum(summary_table[[i]]$GIP_adj_pval[[4]][,1] %in% 
                                                          c("p <= 0.001", "p <= 0.01", "p <= 0.05"))
  TC <- sum(summary_table[[i]]$GIP_adj_pval[[5]][,1] %in% 
                                                          c("p <= 0.001", "p <= 0.01", "p <= 0.05"))
  
  summary_table[[i]]$GIP_significant <- matrix(c(DC, KC, LC, MC, TC), ncol = length(object$X.dim), nrow = 1)
  colnames(summary_table[[i]]$GIP_significant) <- c("DC", "KC", "LC", "MC", "TC")
  # remove object
  #rm(object)
  #rm(model_validation)
}
