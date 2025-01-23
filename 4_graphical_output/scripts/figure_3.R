# Generate data frame with the superpathway recapitulations
df_pathway <- do.call(rbind, lapply(Pathway_recap, function(x){return(x$Pathway)}))

# Add superpathway pvalues
df_pathway <- merge(df_pvalues, df_pathway, by = "Pathway_name")
df_pathway <- df_pathway[order(df_pathway$pvalues, decreasing = TRUE), ]# Order by pvalues

# Plots and descriptive statistics

theme_heat <- theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5)) # Common theme

# Orthology statistics of disease model organisms
df_long <- df_pathway[, c("OXA_orth", "OVA_orth", "IMQ_orth", "Pathway_full_name")] %>% 
  pivot_longer(
    cols = ends_with("_orth"), 
    names_to = "DiseaseModel", 
    values_to = "Orthology"
  ) %>%  mutate(DiseaseModel = str_extract(DiseaseModel, "^[^_]*")) %>% as.data.frame

df_long$Pathway_full_name <- factor(df_long$Pathway_full_name, levels = unique(df_long$Pathway_full_name),
                                    ordered = TRUE) # Factor to preserve order

orthology.B <- ggplot(data = df_long, aes(fill=DiseaseModel, x = Pathway_full_name, y = Orthology)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.05), labels = percent) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "", 
       y = "Observed one-to-one orthology coverage", #title = "Pathway orthology coverage by disease model",
       fill = "Disease model", tag = "B") +
  theme_heat +
  coord_flip()


# Superpathway recapitulation plot
df_long.A <- df_pathway[, c("OXA", "OVA", "IMQ", "Pathway_full_name", "adj.pvalues")] %>% 
  pivot_longer(
    cols = c("OXA", "OVA", "IMQ"), 
    names_to = "DiseaseModel", 
    values_to = "Recapitulation"
  ) %>%  mutate(DiseaseModel = str_extract(DiseaseModel, "^[^_]*")) %>% #%>% 
  #mutate(Pathway_full_name = case_when(
  #adj.pvalues == "p <= 0.001" ~ paste0(Pathway_full_name, "***"),
  #adj.pvalues == "p <= 0.01" ~ paste0(Pathway_full_name, "**"),
  #adj.pvalues == "p <= 0.05" ~ paste0(Pathway_full_name, "*"),
  #TRUE ~ Pathway_full_name
  #)) 
  as.data.frame

df_long.A$Pathway_full_name <- factor(df_long.A$Pathway_full_name, levels = unique(df_long.A$Pathway_full_name),
                                      ordered = TRUE) # Factor to preserve order

superpathway_recapitulation.A <- ggplot(df_long.A, aes(x = DiseaseModel, y = Pathway_full_name)) +
  geom_tile(aes(fill = Recapitulation), color = "white") +
  theme_heat + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "green", midpoint = 0, limits=c(-200,200)) +
  geom_text(aes(label = paste0(round(Recapitulation, 1), "%")), 
            color = "black") +
  labs(x = "Disease Model", y = "Pathway name", #title = "Superpathway observed recapitulation",
       #subtitle = "Atopic Dermatitis disease models",
       fill = TeX("$\\widehat{\\Omega} f_{AD}$"), tag = "A") +
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.subtitle = element_text(face = "bold", color = "grey35")) +
  theme(plot.caption = element_text(color = "grey68"))
