## Source data of this plot is in supplementary files of the paper.

library(ggplot2)
library(ggrepel)
Tumor_long_post_vs_pre_DEG <- read_csv("path/to/Figure_5d.csv")
colnames(Tumor_long_post_vs_pre_DEG)[1] <- "Gene"

Gene_highlight <- c("UPK1B","WNK4","MIF","NT5E","H1-5","MYCN","THY1","IFITM1","IL32","KLK5","H3C2","TAGLN","CLDN1",
                    "CLDN16","SPP1","UQCR10","TP53","KLK10","WNT7A","CCND1","CXCL12","LIF","CDK18","NDUFB6","HMGA1",
                    "HMGB2","CDK1","IL4R","PCNA","PKM","GREB1","HLA-C","IFITM3","ISG20","MKI67","IGHG4","APOD","IGKC",
                    "IGHG3","CXCL14","IGHG2","IL6","CCL4L2","MGP","NR4A1","FOS","DUSP1","JUNB","HSP90B1","EGR1","ATF3")

Tumor_long_post_vs_pre_DEG$invert_P <- (-log10(Tumor_long_post_vs_pre_DEG$adj.P.Val)) 
colnames(Tumor_long_post_vs_pre_DEG)[1] <- "Gene"
Tumor_long_post_vs_pre_DEG$Color <- "NS"
Tumor_long_post_vs_pre_DEG$Color[Tumor_long_post_vs_pre_DEG$adj.P.Val < 0.1 & Tumor_long_post_vs_pre_DEG$logFC >0] <- "FDR < 0.1, pos"
Tumor_long_post_vs_pre_DEG$Color[Tumor_long_post_vs_pre_DEG$adj.P.Val < 0.05 & Tumor_long_post_vs_pre_DEG$logFC >0] <- "FDR < 0.05, pos"
Tumor_long_post_vs_pre_DEG$Color[Tumor_long_post_vs_pre_DEG$adj.P.Val < 0.01 & Tumor_long_post_vs_pre_DEG$logFC >0] <- "FDR < 0.01, pos"
Tumor_long_post_vs_pre_DEG$Color[Tumor_long_post_vs_pre_DEG$adj.P.Val < 0.001 & Tumor_long_post_vs_pre_DEG$logFC >0] <- "FDR < 0.001, pos"
Tumor_long_post_vs_pre_DEG$Color[Tumor_long_post_vs_pre_DEG$adj.P.Val < 0.1 & Tumor_long_post_vs_pre_DEG$logFC < 0] <- "FDR < 0.1, neg"
Tumor_long_post_vs_pre_DEG$Color[Tumor_long_post_vs_pre_DEG$adj.P.Val < 0.05 & Tumor_long_post_vs_pre_DEG$logFC < 0] <- "FDR < 0.05, neg"
Tumor_long_post_vs_pre_DEG$Color[Tumor_long_post_vs_pre_DEG$adj.P.Val < 0.01 & Tumor_long_post_vs_pre_DEG$logFC < 0] <- "FDR < 0.01, neg"
Tumor_long_post_vs_pre_DEG$Color[Tumor_long_post_vs_pre_DEG$adj.P.Val < 0.001 & Tumor_long_post_vs_pre_DEG$logFC < 0] <- "FDR < 0.001, neg"
Tumor_long_post_vs_pre_DEG$Color <- factor(Tumor_long_post_vs_pre_DEG$Color, levels = c("NS","FDR < 0.1, pos","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos",
                                                                                        "FDR < 0.1, neg","FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg"))




ggplot(Tumor_long_post_vs_pre_DEG, aes(x = logFC, y = invert_P, color = Color,label = Gene)) +
  #geom_vline(xintercept = c(1, -1), lty = "dashed") +
  #geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point(size=0.3) +
  labs(x = "Pre-treatment <- log2(FC) -> On-treatment",
       y = "-log10(p.adj)",
       color = "Significance") +
  scale_color_manual(values = c(
    `FDR < 0.001, pos` = "#bd0026",
    `FDR < 0.01, pos` = "#fc4e2a",
    `FDR < 0.05, pos` = "#feb24c",
    `FDR < 0.1, pos` = "#ffeda0",
    `FDR < 0.001, neg` = "#08306b",
    `FDR < 0.01, neg` = "#2171b5",
    `FDR < 0.05, neg` = "#6baed6",
    `FDR < 0.1, neg` = "#c6dbef",
    `NS` = "#f0f0f0"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(Tumor_long_post_vs_pre_DEG, Gene %in% Gene_highlight),
                  size = 3, point.padding = 0, color = "black",
                  min.segment.length = 0, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right") 
