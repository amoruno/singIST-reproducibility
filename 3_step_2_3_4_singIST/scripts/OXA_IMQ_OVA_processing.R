# Load asmbPLS-DA fitted pathway object
all_input <- list.files("~/singIST-reproducibility/2_step1_singIST/exported_results/asmbPLSDA_training/", pattern = "\\.RData$")

# IMQ/OXA
imq_oxa <- readRDS(file = "~/OXA_IMQ.rds")
imq_oxa_ <- UpdateSeuratObject(imq_oxa)

# OVA
ova <- readRDS(file = "~/OVA.rds")
ova_ <- UpdateSeuratObject(ova)

# Identify clusters according to original annotations
imq_oxa_clust <- imq_oxa_[,imq_oxa_$integrated_snn_res.0.35 %in% c(1, 3, 4, 6, 12, 14, 15, 16)]
imq_oxa_clust$stim <- factor(imq_oxa_clust$stim, levels = c("OXA", "ETOH", "IMQ", "VEH"), ordered = TRUE)

ova_clust <- ova_[, ova_$celltype %in% c("T cells", "Keratinocytes",
                                                    "Dendritic cells")] # Langerhans Cells are not considered since the disease model has < 100 cells
                                                                         
ova_clust$group <- factor(ova_clust$group, levels = c("OVA", "SAL"), ordered = TRUE)

# Homogenize cell type annotation
# IMQ/OXA
imq_oxa_metadata <- imq_oxa_clust@meta.data
imq_oxa_metadata[imq_oxa_metadata$integrated_snn_res.0.35 %in% c(1,3,6), "cell_type_new"] <- "T-cell"
imq_oxa_metadata[imq_oxa_metadata$integrated_snn_res.0.35 %in% c(4,14,16), "cell_type_new"] <- "Dendritic Cells"
imq_oxa_metadata[imq_oxa_metadata$integrated_snn_res.0.35 %in% c(12), "cell_type_new"] <- "Keratinocytes"
imq_oxa_metadata[imq_oxa_metadata$integrated_snn_res.0.35 %in% c(15), "cell_type_new"] <- "Langerhans Cells"

imq_oxa_clust@meta.data <- imq_oxa_metadata
imq_oxa_pseudobulk <- AggregateExpression(imq_oxa_clust, assays = "RNA", group.by = c("cell_type_new", "orig.ident", "stim"),
                                          normalization.method = "LogNormalize", return.seurat = TRUE)
imq_oxa_pseudobulk_matrix <- t(as.data.frame(LayerData(imq_oxa_pseudobulk, "data")))

# OVA
# Harmonize cell type annotation
ova_metadata <- ova_clust@meta.data
ova_metadata[ova_metadata$celltype == "T cells", "cell_type_new"] <- "T-cell"
ova_metadata[ova_metadata$celltype == "Dendritic cells", "cell_type_new"] <- "Dendritic Cells"
ova_metadata[ova_metadata$celltype == "Keratinocytes", "cell_type_new"] <- "Keratinocytes"
#ova_metadata[ova_metadata$celltype == "Langerhans cells", "cell_type_new"] <- "Langerhans Cells"

ova_clust@meta.data <- ova_metadata
ova_pseudobulk <- AggregateExpression(ova_clust, assays = "RNA", group.by = c("cell_type_new", "orig.ident", "group"), 
                                          normalization.method = "LogNormalize", return.seurat = TRUE)
ova_pseudobulk_matrix <- t(as.data.frame(LayerData(ova_pseudobulk, "data")))

# Compute DEGs
imq_oxa_clust$celltype.stim <- paste0(imq_oxa_clust$cell_type_new, "_", imq_oxa_clust$stim)
Idents(imq_oxa_clust) <- "celltype.stim" 

logFC_OXA <- diff_expressed(imq_oxa_clust, condition_1 = unique(imq_oxa_clust$cell_type_new), 
                            condition_2 = c("OXA", "ETOH"), logfc.treshold = 0.25)
logFC_IMQ <- diff_expressed(imq_oxa_clust, condition_1 = unique(imq_oxa_clust$cell_type_new), 
                            condition_2 = c("IMQ", "VEH"), logfc.treshold = 0.25)

ova_clust$celltype.stim <- paste0(ova_clust$cell_type_new, "_", ova_clust$group)
Idents(ova_clust) <- "celltype.stim" 
logFC_OVA <- diff_expressed(ova_clust, condition_1 = unique(ova_clust$cell_type_new), 
                            condition_2 = c("OVA", "SAL"), logfc.treshold = 0.1)

# Group results
logFC_FindMarkers <- list("Oxazolone" = logFC_OXA, "Imiquimod" = logFC_IMQ, "Ovalbumine" = logFC_OVA)

# Select orthologs for each pathway for mus musculus organism
orthologs <- lapply(seq_along(all_input), function(i, path_input = "/singIST-reproducibility/2_step1_singIST/exported_results/asmbPLSDA_training"){
                                              # Initiliaze output object to return
                                              output <- list("Pathway_name" = NULL, "Gene_set_observed" = c(),
                                                             "Gene_set_ortholog" = c(list("external_gene_name" = c(),
                                                                                          "Hgnc_symbol" = c())))
                                              # Human trained asmbPLS-DA
                                              load(paste0(path_input, "/", all_input[i]))
                                              # Assign pathway name
                                              output$Pathway_name <- object$Pathway_name
                                              # Obtain gene set only of observed genes
                                              aux <- object$Gene_set %in%
                                                        sub(".*\\*(.*)", "\\1", rownames(object$asmbPLSDA.fit$X_weight[[1]]))
                                              output$Gene_set_observed <- object$Gene_set[aux]
                                              print(output$Pathway_name)
                                              # Map ortholog genes
                                              gene_set <- output$Gene_set_observed
                                              aux <- orth_pathway(gene_set)
                                              output$Gene_set_ortholog$external_gene_name <- aux$external_gene_name
                                              output$Gene_set_ortholog$Hgnc_symbol <- aux$hgnc_symbol

                                              return(output)
                                            })

save(orthologs, file = "C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/input_disease_model/Ortholog object/orthologs.RData")

# Obtain superpathway dataset of each disease model
pathway_disease_model <- lapply(seq_along(orthologs), function(i, dataset_oxa_imq = imq_oxa_pseudobulk_matrix, dataset_ova = ova_pseudobulk_matrix, 
                                                              dataset_orthologs = orthologs, logFC = logFC_FindMarkers){
  # Initialize output object
  output <- list(Pathway_name = NULL, Imiquimod = list("Lognorm_counts" = data.frame(), "X.matrix" = matrix(),
                                                       "Orthologs_observed" = c(), "logFC" = data.frame()), 
                 Ovalbumine = list("Lognorm_counts" = data.frame(), "X.matrix" = matrix(),
                                     "Orthologs_observed" = c(), "logFC" = data.frame()),
                 Oxazolone = list("Lognorm_counts" = data.frame(), "X.matrix" = matrix(),
                                    "Orthologs_observed" = c(), "logFC" = data.frame()))
  # Set pathway
  pathway_orthologs <- dataset_orthologs[[i]]
  
  ## IMQ-OXA
  genes_position <- colnames(dataset_oxa_imq) %in% pathway_orthologs$Gene_set_ortholog$external_gene_name
  genes_to_pick <- colnames(dataset_oxa_imq)[genes_position]
  genes_to_pick_hgnc <- pathway_orthologs$Gene_set_ortholog$Hgnc_symbol[match(genes_to_pick, pathway_orthologs$Gene_set_ortholog$external_gene_name)]
  
  data <- as.data.frame(dataset_oxa_imq[, intersect(colnames(dataset_oxa_imq), c(genes_to_pick))]) # change to only "stim"
  colnames(data)[1:length(genes_to_pick)] <- genes_to_pick_hgnc
  
  # Split rownames into sample id and cell type 
  split <- do.call(rbind, strsplit(rownames(data), "_"))
  data$cell_type_new <- split[,1]
  data$orig.ident <- split[,2]
  data$stim <- split[,3]
  
  # Average the pseudobulk samples to calculate logFC later on
  data <- aggregate(. ~ cell_type_new + stim, data[, colnames(data) != "orig.ident"], mean)
  
  # Select experimental units
  data_imq <- data[data$stim == "IMQ" | data$stim == "VEH", ]
  data_oxa <- data[data$stim == "OXA" | data$stim == "ETOH", ]
  
  # Generate X.matrix block 
  cell_type_new <- unique(data_imq$cell_type_new) # cell types are imq = oxa
  X.matrix_imq <- data.frame(stim = factor(c("IMQ", "VEH"), levels = c("IMQ", "VEH"), ordered = TRUE))
  X.matrix_oxa <- data.frame(stim = factor(c("OXA", "ETOH"), levels = c("OXA", "ETOH"), ordered = TRUE))

  for(i in cell_type_new){
    # X.matrix block
    blockcells_imq <- data_imq[data_imq$cell_type_new == i, ]
    colnames(blockcells_imq)[3:ncol(blockcells_imq)] <- paste0(i, "*", colnames(blockcells_imq)[3:ncol(blockcells_imq)])
    X.matrix_imq <- merge(X.matrix_imq, blockcells_imq[, colnames(blockcells_imq) != "cell_type_new"], 
                          by = "stim")
    
    blockcells_oxa <- data_oxa[data_oxa$cell_type_new == i, ]
    colnames(blockcells_oxa)[3:ncol(blockcells_oxa)] <- paste0(i, "*", colnames(blockcells_oxa)[3:ncol(blockcells_oxa)])
    X.matrix_oxa <- merge(X.matrix_oxa, blockcells_oxa[, colnames(blockcells_oxa) != "cell_type_new"], 
                          by = "stim")
    
    ## logFC from FindMarkers
    if(i %in% names(logFC$Oxazolone)){
      intersect_genes <- intersect(genes_to_pick, rownames(logFC$Oxazolone[[i]]))
      if(length(intersect_genes) > 0){
        logFC$Oxazolone[[i]] <- logFC$Oxazolone[[i]][intersect_genes, ]
        # Compute sign(log2FC)*2^log2FC for statistically significant genes
        significant_genes <- logFC$Oxazolone[[i]]$p_val_adj <= 0.05
        logFC$Oxazolone[[i]][significant_genes, "avg_log2FC"] <- sign(logFC$Oxazolone[[i]][significant_genes, "avg_log2FC"])*2^logFC$Oxazolone[[i]][significant_genes, "avg_log2FC"]
        # Otherwise assign to 0
        logFC$Oxazolone[[i]][!significant_genes, "avg_log2FC"] <- 0
        # Set rownames to HGNC symbol
        rownames(logFC$Oxazolone[[i]]) <- genes_to_pick_hgnc[match(rownames(logFC$Oxazolone[[i]]), genes_to_pick)]
      }
    }
    if(i %in% names(logFC$Imiquimod)){
      intersect_genes <- intersect(genes_to_pick, rownames(logFC$Imiquimod[[i]]))
      if(length(intersect_genes) > 0){
        logFC$Imiquimod[[i]] <- logFC$Imiquimod[[i]][intersect_genes, ]
        # Compute sign(log2FC)*2^log2FC for statistically significant genes
        significant_genes <- logFC$Imiquimod[[i]]$p_val_adj <= 0.05
        logFC$Imiquimod[[i]][significant_genes, "avg_log2FC"] <- sign(logFC$Imiquimod[[i]][significant_genes, "avg_log2FC"])*2^logFC$Imiquimod[[i]][significant_genes, "avg_log2FC"]
        # Otherwise assign to 0
        logFC$Imiquimod[[i]][!significant_genes, "avg_log2FC"] <- 0
        # Set rownames to HGNC symbol
        rownames(logFC$Imiquimod[[i]]) <- genes_to_pick_hgnc[match(rownames(logFC$Imiquimod[[i]]), genes_to_pick)]
      }
    }
    output$Oxazolone$logFC <- do.call(rbind, logFC$Oxazolone)
    rownames(output$Oxazolone$logFC) <- gsub("\\.", "*", rownames(output$Oxazolone$logFC))
    output$Imiquimod$logFC <- do.call(rbind, logFC$Imiquimod)
    rownames(output$Imiquimod$logFC) <- gsub("\\.", "*", rownames(output$Imiquimod$logFC))
  }
  
  # Sort the merged data frame by the ID column
  X.matrix_oxa <- X.matrix_oxa[order(X.matrix_oxa$stim), ]
  X.matrix_imq <- X.matrix_imq[order(X.matrix_imq$stim), ]

  # Save to output object
  output$Pathway_name <- pathway_orthologs$Pathway_name
  output$Imiquimod$Lognorm_counts <- data_imq
  output$Oxazolone$Lognorm_counts <- data_oxa
  output$Imiquimod$X.matrix <- X.matrix_imq
  output$Oxazolone$X.matrix <- X.matrix_oxa
  output$Imiquimod$Orthologs_observed <- output$Oxazolone$Orthologs_observed <- genes_to_pick_hgnc
  
  ## OVA
  genes_position <- colnames(dataset_ova) %in% pathway_orthologs$Gene_set_ortholog$external_gene_name
  genes_to_pick <- colnames(dataset_ova)[genes_position]
  genes_to_pick_hgnc <- pathway_orthologs$Gene_set_ortholog$Hgnc_symbol[match(genes_to_pick, pathway_orthologs$Gene_set_ortholog$external_gene_name)]
  
  data <- as.data.frame(dataset_ova[, intersect(colnames(dataset_ova), c(genes_to_pick))]) # change to grooup and cell type new
  colnames(data)[1:length(genes_to_pick)] <- genes_to_pick_hgnc
  
  # Split rownames into sample id and cell type 
  split <- do.call(rbind, strsplit(rownames(data), "_"))
  data$cell_type_new <- split[,1]
  data$orig.ident <- split[,2]
  data$group <- split[,3]
  
  # Average the pseudobulk samples to calculate logFC later on
  data <- aggregate(. ~ cell_type_new + group, data[, colnames(data) != "orig.ident"], mean)
  
  # Select experimental units
  data_ova <- data[data$group == "OVA" | data$group == "SAL", ]
  
  # Generate X.matrix block 
  cell_type_new <- unique(data_ova$cell_type_new)
  X.matrix_ova <- data.frame(group = factor(c("OVA", "SAL"), levels = c("OVA", "SAL"), ordered = TRUE))
  
  for(i in cell_type_new){
    blockcells_ova <- data_ova[data_ova$cell_type_new == i, ]
    colnames(blockcells_ova)[3:ncol(blockcells_ova)] <- paste0(i, "*", colnames(blockcells_ova)[3:ncol(blockcells_ova)])
    X.matrix_ova <- merge(X.matrix_ova, blockcells_ova[, colnames(blockcells_ova) != "cell_type_new"], 
                          by = "group")
    ## logFC from FindMarkers
    if(i %in% names(logFC$Ovalbumine)){
      intersect_genes <- intersect(genes_to_pick, rownames(logFC$Ovalbumine[[i]]))
      if(length(intersect_genes) > 0){
        logFC$Ovalbumine[[i]] <- logFC$Ovalbumine[[i]][intersect_genes, ]
        # Compute sign(log2FC)*2^log2FC for statistically significant genes
        significant_genes <- logFC$Ovalbumine[[i]]$p_val_adj <= 0.05
        logFC$Ovalbumine[[i]][significant_genes, "avg_log2FC"] <- sign(logFC$Ovalbumine[[i]][significant_genes, "avg_log2FC"])*2^logFC$Ovalbumine[[i]][significant_genes, "avg_log2FC"]
        # Otherwise assign to 0
        logFC$Ovalbumine[[i]][!significant_genes, "avg_log2FC"] <- 0
        # Set rownames to HGNC symbol
        rownames(logFC$Ovalbumine[[i]]) <- genes_to_pick_hgnc[match(rownames(logFC$Ovalbumine[[i]]), genes_to_pick)]
      }
    }
    output$Ovalbumine$logFC <- do.call(rbind, logFC$Ovalbumine)
    rownames(output$Ovalbumine$logFC) <- gsub("\\.", "*", rownames(output$Ovalbumine$logFC))
  }
  
  # Sort the merged data frame by the ID column
  X.matrix_ova <- X.matrix_ova[order(X.matrix_ova$group), ]
  
  # Save to output object
  output$Ovalbumine$Lognorm_counts <- data_ova
  output$Ovalbumine$X.matrix <- X.matrix_ova
  output$Ovalbumine$Orthologs_observed <- genes_to_pick_hgnc
  
  save_path <- paste0("C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/input_disease_model/Disease_Models_", output$Pathway_name, ".RData")
  save(output, file = save_path)
  return(output)
})
