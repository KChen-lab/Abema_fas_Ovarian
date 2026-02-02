## Source data of this plot is in supplementary files of the paper.
## GSEA of on-treatment vs. pre-treatment samples in patients with long PFS (tumor compartment) using published metaprograms (Barkley, et. al Nat Genet 2020)

library(ggplot2)

Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=coverage, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title="pathways NES from GSEA between on and pre treatment in Tumor with Long PFS") + 
  theme_minimal()+ scale_color_gradient(low="red", high="blue")+
  scale_size(range = c(4, 6), name="Coverage, [%]")+
  theme(axis.text = element_text(size=18),strip.text = element_text(size=18))+
  theme(legend.position  = "right", legend.text = element_text(size=18),title = element_text(size=18))
Bubble


