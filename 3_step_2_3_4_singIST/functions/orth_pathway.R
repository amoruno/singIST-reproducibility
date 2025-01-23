# Description: orthology mapping between mus musculus organism and humans
orth_pathway <- function(gene_set){
  # Get mouse orthologs
  mart <- biomaRt::useMart(
    biomart = "ensembl", 
    dataset = paste0("hsapiens", "_gene_ensembl"))
  
  orthologs <- getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", 
                                    "mmusculus_homolog_orthology_confidence", 
                                    "mmusculus_homolog_orthology_type"),
                     filters    = list(hgnc_symbol = gene_set, with_mmusculus_homolog = TRUE),
                     mart       = mart)
  
  orthologs <- data.table::data.table(orthologs)
  names(orthologs) <- c("gene", "ortholog", "confidence", "type")
  
  orthologs <- orthologs[type == "ortholog_one2one", ]
  
  # filter out m to m relationship
  orthologs <- orthologs[
    !gene %in% orthologs[, .N, by = gene][N > 1, gene] & 
      !ortholog %in% orthologs[, .N, by = ortholog][N > 1, ortholog],
  ]
  
  orthologs_ensembl <- as.vector(orthologs$ortholog)
  
  # Extract mouse gene symbol
  mart <- biomaRt::useMart(
    biomart = "ensembl", 
    dataset = paste0("mmusculus", "_gene_ensembl"))
  
  rat_gene_symbol <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                           filters    = list(ensembl_gene_id = orthologs_ensembl),
                           mart       = mart)
  
  colnames(rat_gene_symbol)[1] <- "ortholog"
  
  # Extract human gene symbol
  mart <- biomaRt::useMart(
    biomart = "ensembl", 
    dataset = paste0("hsapiens", "_gene_ensembl"))
  
  gene_ensembl <- orthologs$gene
  
  human_gene_symbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                             filters    = list(ensembl_gene_id = gene_ensembl),
                             mart       = mart)
  
  colnames(human_gene_symbol)[1] <- "gene"
  
  # Final ortholog dataset, human gene symbol, rat gene symbol
  aux <- merge(orthologs, rat_gene_symbol, by = "ortholog")
  final <- merge(aux, human_gene_symbol, by = "gene")
  return(final)
}
