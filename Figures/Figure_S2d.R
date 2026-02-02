## Source data of this plot is in supplementary files of the paper.
## GSEA of on-treatment vs. pre-treatment samples in all patients (stromal compartment)

set.seed(2023)
library(ggplot2)


Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=covered, color=padj))+ #
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between on and pre treatment in TME") + 
  theme_minimal()+ scale_color_gradient(low="red", high="blue")+
  scale_size(range = c(4, 6), name="Coverage, [%]")+
  theme(axis.text = element_text(size=18),strip.text = element_text(size=18))+
  theme(legend.position  = "right", legend.text = element_text(size=18),title = element_text(size=18))
Bubble
