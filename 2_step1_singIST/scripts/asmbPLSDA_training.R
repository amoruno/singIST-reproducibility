# Load list with lognorm counts, gene set and pathway name
all_input <- list.files("~/singIST-reproducibility/1_rawdata/exported_results/", pattern = "\\.RData$")

for(i in all_input){
  # Load input
  load(paste0("~/singIST-reproducibility/1_rawdata/exported_results", "/", i))

  training_pathway <- 
    lapply(1, function(void){
           output_dataset$Lognorm_counts <- as.data.frame(output_dataset$Lognorm_counts)
           pathway <- output_dataset

           # Split rownames into sample id and cell type 
           split <- do.call(rbind, strsplit(rownames(pathway$Lognorm_counts), "_"))
           pathway$Lognorm_counts$CELLTYPE_new <- split[,1]
           pathway$Lognorm_counts$Sample_id <- split[,2]
           
           # Filter only cell types to analyze
           pathway$Lognorm_counts <- pathway$Lognorm_counts[pathway$Lognorm_counts$CELLTYPE_new %in% c("T-cell", "Keratinocytes",
                                                                                                       "Dendritic Cells", "Langerhans Cells",
                                                                                                        "Melanocytes"),]
           # Create X.matrix block
           df_X.matrix <- data.frame("Sample_id" = unique(pathway$Lognorm_count$Sample_id))
           aux_X.matrix <- pathway$Lognorm_counts
           cell_type_new <- unique(pathway$Lognorm_counts$CELLTYPE_new)
           for(i in cell_type_new){
             blockcells <- aux_X.matrix[aux_X.matrix$CELLTYPE_new == i, ]
             colnames(blockcells)[1:(which(colnames(aux_X.matrix) == "CELLTYPE_new")-1)] <- paste0(i, "*", colnames(blockcells)[1:(which(colnames(aux_X.matrix) == "CELLTYPE_new")-1)])
             df_X.matrix <- merge(df_X.matrix, blockcells[, colnames(blockcells) != "CELLTYPE_new"], 
                                  by = "Sample_id")
           }
             
           # Train asmbPLS-DA
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
           filepath <- paste0("~/singIST-reproducibility/2_step1_singIST/exported_results/asmbPLSDA_training/",
                              pathway$Pathway_name, ".RData")
           object <- pathway
           save(object , file = filepath)
           return(pathway)}
           )
}
