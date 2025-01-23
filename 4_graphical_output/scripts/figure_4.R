# Plot cell type recapitulations by pathway and disease model
df_cell <- do.call(rbind, lapply(Pathway_recap, function(x){return(x$Cell)}))

# Add adjusted global significance pvalues of superpathway asmbsPLS-DA classifier
df_cell <- merge(df_pvalues, df_cell, by = "Pathway_name")
df_cell <- df_cell[order(df_cell$pvalues, decreasing = FALSE), ]# Order by pvalues

# Data frame containing cell type recapitulations 
df_cell_long <- df_cell[, c("OXA", "OVA", "IMQ", "Pathway_full_name", "cell_type", "adj.pvalues")] %>% 
  pivot_longer(
    cols = c("OXA", "OVA", "IMQ"), 
    names_to = "DiseaseModel", 
    values_to = "Recapitulation"
  ) %>%  mutate(DiseaseModel = str_extract(DiseaseModel, "^[^_]*")) %>% as.data.frame

# Make custom theme
theme_heat <- theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size=5.7),
        axis.text.x = element_text(angle = -65, hjust = 0, size = 6),
        axis.text.y = element_text(size = 6))

# Default 40 character target width.
swr = function(string, nwrap=40) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

# Create line breaks in Pathway name
df_cell_long$Pathway_full_name = swr(df_cell_long$Pathway_full_name)

# Plot for each pathway, cell type and disease model
cell_type_recapitulation_plot <- ggplot(df_cell_long, aes(x = DiseaseModel, y = cell_type)) +
  geom_tile(aes(fill = Recapitulation), color = "white") +
  facet_wrap(~Pathway_full_name, scales = "free") + theme_heat +
  scale_fill_gradient2(low = "blue", mid = "white", high = "green", 
                       midpoint = 0, limits=c(-200, 200))

# plot with text overlay and viridis color palette
cell_type_recapitulation_plot.A <- cell_type_recapitulation_plot + geom_text(aes(label = paste0(round(Recapitulation, 1), "%")), 
                       color = "black", size = 2) +
  labs(x = "Disease model", y= "Cell type",
    fill = TeX("$\\widehat{\\Gamma} f_{AD}^b$")) +
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.subtitle = element_text(face = "bold", color = "grey35")) +
  theme(plot.caption = element_text(color = "grey68"))
