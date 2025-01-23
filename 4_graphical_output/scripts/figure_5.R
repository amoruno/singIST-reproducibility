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

