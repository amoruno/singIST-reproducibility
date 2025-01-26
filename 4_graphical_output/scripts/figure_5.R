# Generate data frame with gene contributions to cell type recapitulation, each gene, cell type, disease model and pathway
df_gene <- do.call(rbind, lapply(Pathway_recap, function(x){return(x$Gene)}))
df_gene <- merge(df_pvalues, df_gene, by = "Pathway_name")

# Generate data frame with FC for each gene, cell type, disease model and pathway
df_FC <- do.call(rbind, lapply(Pathway_recap, function(x){
  return(rbind(x$FC$OXA, x$FC$OVA, x$FC$IMQ))
  })
  )
df_FC <- merge(df_pvalues, df_FC, by = "Pathway_name")

save(df_gene, file = "C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/output_recapitulation/df_gene.RData")
save(df_FC, file = "C:/Users/amoruno/OneDrive - Almirall S.A/Doctorat/Publication 1/All pathways analysis/output_recapitulation/df_FC.RData")


# Add a custom legend for NaN values
lgd <- Legend(
  labels = "Not ortholog",
  legend_gp = gpar(fill = "grey"),
  title = "Orthology status"
)

# Define the legend for black borders
border_legend <- Legend(
  labels = c(expression("FDR <= 0.05")),
  title = "FC significance",
  type = "lines",
  legend_gp = gpar(col = "black", fill = NA, lwd = 2),
  nrow = 1, by_row = TRUE, labels_gp = gpar(fontsize = 8)
)


# Dendritic Cells
df_dc <- df_gene[df_gene$Pathway_full_name == "Dendritic cells in Th1/Th2 Development [BIOCARTA]", ]
pathway_dc <- heatmap_gene(df_dc, pathway_full_name = "Dendritic cells in Th1/Th2 Development [BIOCARTA]",
                            value_size = 5, gene_size = 6, height = 1000, width = 300)

# Plot applied fold changes
df_FC_dc <- df_FC[df_FC$Pathway_full_name ==  "Dendritic cells in Th1/Th2 Development [BIOCARTA]", ]
df_FC_dc <- df_FC_dc[df_FC_dc$gene_name %in% unique(df_dc$Gene_name),]
pathway_dc_FC <- heatmap_FC(df_FC_dc, pathway_full_name = "Dendritic cells in Th1/Th2 Development [BIOCARTA]",
                             value_size = 6, gene_size = 6, height = 1000, width = 300)

draw(pathway_dc+pathway_dc_FC, 
     show_heatmap_legend = FALSE, 
     show_annotation_legend = FALSE,
     #annotation_legend_list = list(lgd, border_legend), 
     column_title = "Dendritic cells in Th1/Th2 Development [BIOCARTA]")

# JAK-STAT
# Apply the function to select top X genes
df_jak <- select_top_genes(df_gene, pathway_full_name = "JAK-STAT signaling pathway [KEGG]", top_n = 5)
pathway_jak <- heatmap_gene(df_jak, pathway_full_name = "JAK-STAT signaling pathway [KEGG]",
                            value_size = 5, gene_size = 5, height = 1000, width = 300)

# Plot applied fold changes
df_FC_jak <- df_FC[df_FC$Pathway_full_name ==  "JAK-STAT signaling pathway [KEGG]", ]
df_FC_jak <- df_FC_jak[df_FC_jak$gene_name %in% unique(df_jak$Gene_name),]
pathway_jak_FC <- heatmap_FC(df_FC_jak, pathway_full_name = "JAK-STAT signaling pathway [KEGG]",
                             value_size = 5, gene_size = 5, height = 1000, width = 300)
draw(pathway_jak+pathway_jak_FC, 
     show_heatmap_legend = FALSE, 
     show_annotation_legend = FALSE,
     #annotation_legend_list = list(lgd, border_legend), 
     column_title = "JAK-STAT signaling pathway [KEGG]")

# Cytokine-Cytokine
# Apply the function to select top X genes
df_cyt <- select_top_genes(df_gene, pathway_full_name = "Cytokine-Cytokine receptor interaction [KEGG]", top_n = 5)

pathway_cyt <- heatmap_gene(df_cyt, pathway_full_name = "Cytokine-Cytokine receptor interaction [KEGG]",
                          value_size = 5, gene_size = 5, height = 1000, width = 300)

# Plot applied fold changes
df_FC_cyt <- df_FC[df_FC$Pathway_full_name ==  "Cytokine-Cytokine receptor interaction [KEGG]", ]
df_FC_cyt <- df_FC_cyt[df_FC_cyt$gene_name %in% unique(df_cyt$Gene_name),]
pathway_cyt_FC <- heatmap_FC(df_FC_cyt, pathway_full_name = "Cytokine-Cytokine receptor interaction [KEGG]",
                          value_size = 5, gene_size = 5, height = 1000, width = 300)
draw(pathway_cyt+pathway_cyt_FC, 
     show_heatmap_legend = FALSE, 
     show_annotation_legend = FALSE,
     # annotation_legend_list = list(lgd, border_legend),
     column_title = "Cytokine-Cytokine receptor interaction [KEGG]")


# Chemokine 
# Apply the function to select top X genes
df_che <- select_top_genes(df_gene, pathway_full_name = "Chemokine receptors bind chemokines [REACTOME]", top_n = 5)
pathway_che <- heatmap_gene(df_che, pathway_full_name = "Chemokine receptors bind chemokines [REACTOME]",
                            value_size = 5, gene_size = 5, height = 1000, width = 300)

# Plot applied fold changes
df_FC_che <- df_FC[df_FC$Pathway_full_name ==  "Chemokine receptors bind chemokines [REACTOME]", ]
df_FC_che <- df_FC_che[df_FC_che$gene_name %in% unique(df_che$Gene_name),]
pathway_che_FC <- heatmap_FC(df_FC_che, pathway_full_name = "Chemokine receptors bind chemokines [REACTOME]",
                             value_size = 5, gene_size = 5, height = 1000, width = 300)
draw(pathway_che+pathway_che_FC, 
     show_heatmap_legend = FALSE, 
     show_annotation_legend = FALSE,
     #annotation_legend_list = list(lgd, border_legend), 
     column_title = "Chemokine receptors bind chemokines [REACTOME]")
