# libraries
library(BiocGenerics)
library(BiocFileCache)
library(AnnotationHub)
library(ExperimentHub)
library(msigdb)

# Load list with lognorm counts, gene set and pathway name
all_input <- list.files("C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/input")

for(i in all_input){
  # Load input
  load(paste0("C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/input", "/", i))

  training_pathway <- 
    lapply(1, function(void){
           # Filter only HC and AD lesional samples
           pathway <- output_dataset
           pathway$Lognorm_counts <- pathway$Lognorm_counts[pathway$Lognorm_counts$treatment_new %in% c("HC", "Baseline"),]
           # Homogenize cell type annotation
           homogenize_cell_type <- function(celltype) {
             celltype <- ifelse(grepl("T", celltype, fixed = TRUE), "T-cell", celltype)
             celltype <- ifelse(grepl("KC", celltype, fixed = TRUE), "Keratinocytes", celltype)
             celltype <- ifelse(grepl("Melanocytes", celltype, fixed = TRUE), "Melanocytes", celltype)
             celltype <- ifelse(grepl("DC", celltype, fixed = TRUE), "Dendritic Cells", celltype)
             celltype <- ifelse(grepl("LC", celltype, fixed = TRUE), "Langerhan Cells", celltype)
             return(celltype)
           }
           pathway$Lognorm_counts$CELLTYPE_new <- apply(pathway$Lognorm_counts["CELLTYPE_new"], 1, homogenize_cell_type)
           # Filter only cell types to analyze
           pathway$Lognorm_counts <- pathway$Lognorm_counts[pathway$Lognorm_counts$CELLTYPE_new %in% c("T-cell", "Keratinocytes",
                                                                                                       "Dendritic Cells", "Langerhan Cells",
                                                                                                        "Melanocytes"),]
           # Create block of matrices X with the pseudobulk of log-normalized counts
           aux_X.matrix <- aggregate(. ~ CELLTYPE_new + Sample_id,
                                     pathway$Lognorm_counts[, colnames(pathway$Lognorm_counts) 
                                                            %in% c(pathway$Gene_set, "Sample_id", "CELLTYPE_new")],
                                     sum
                                     )
           df_X.matrix <- data.frame("Sample_id" = unique(pathway$Lognorm_count$Sample_id))
           cell_type_new <- unique(pathway$Lognorm_counts$CELLTYPE_new)
           for(i in cell_type_new){
             blockcells <- aux_X.matrix[aux_X.matrix$CELLTYPE_new == i, ]
             colnames(blockcells)[(which(colnames(aux_X.matrix) == "Sample_id")+1):ncol(blockcells)] <- paste0(i, "*", colnames(blockcells)[(which(colnames(aux_X.matrix) == "Sample_id")+1):ncol(blockcells)])
             df_X.matrix <- merge(df_X.matrix, blockcells[, colnames(blockcells) != "CELLTYPE_new"], 
                                  by = "Sample_id")
           }
           # Train asmbPLS-DA
           library(asmbPLS)
           library(RcppAlgos) # permuteGeneral function
           
           ## Load parameters to train asmbPLS-DA
           pathway$Y.matrix <- as.matrix(ifelse(grepl("HC", df_X.matrix$Sample_id, fixed = TRUE), 0, 1))
           pathway$X.matrix <- as.matrix(df_X.matrix[, colnames(df_X.matrix) != "Sample_id"])
           pathway$X.dim <- colSums(sapply(cell_type_new, function(cell_type){ grepl(cell_type, colnames(df_X.matrix))}))
           pathway$quantile.comb.table.cv <- as.matrix(permuteGeneral(seq(0.05, 0.95, by = 0.10), 
                                                              m = length(cell_type_new), 
                                                              TRUE))
           colnames(pathway$quantile.comb.table.cv) <- paste0("block_", rep(1:length(cell_type_new),1))
           
           # Set binary outcome
           pathway$outcome.type <- "binary"
           # LOOCV
           pathway$asmbPLSDA.cv <- asmbPLSDA.cv.loo(X.matrix = pathway$X.matrix, 
                                                    Y.matrix = pathway$Y.matrix, 
                                                    PLS_term = 3,
                                                    X.dim = pathway$X.dim, 
                                                    quantile.comb.table = pathway$quantile.comb.table.cv, 
                                                    outcome.type = pathway$outcome.type,
                                                    center = TRUE,
                                                    scale = TRUE,
                                                    measure = "F1",
                                                    parallel = TRUE,
                                                    cores = NULL)
           # Fit
           pathway$asmbPLSDA.fit <- asmbPLSDA.fit(X.matrix = pathway$X.matrix, 
                                                  Y.matrix = pathway$Y.matrix, 
                                                  PLS.comp = pathway$asmbPLSDA.cv$optimal_nPLS, 
                                                  X.dim = pathway$X.dim, 
                                                  quantile.comb = pathway$asmbPLSDA.cv$quantile_table_CV,
                                                  center = TRUE,
                                                  scale = TRUE,
                                                  outcome.type = "binary")
           # Save objects
           filepath <- paste0("C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/output/",
                              pathway$Pathway_name, ".RData")
           object <- pathway
           save(object , file = filepath)
           
           return(pathway)}
           )
  rm(output_dataset) # remove this 
}
