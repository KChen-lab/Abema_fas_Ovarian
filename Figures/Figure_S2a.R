


library(ggplot2)
library(ggrepel)
Tumor_long_post_vs_pre_DEG <- read_csv("path-to-folder/Figure S2a.csv")
colnames(Tumor_long_post_vs_pre_DEG)[1] <- "Gene"

Gene_highlight <- c("CDK18", 'CLIC3', 'HLA-C','IL32','IL4R', 'UQCR10','HLA-B','COX5B','HLA-A', 'HLA-F','FHL2',"CLDN1",
                    "IDO1","TGM2","CLDN16",'SEMA3B',"TUBB4B","ATP1B1","MDK","WNT7A",
                    "MUC5B",'FGFR3',"CCND1","IGFBP4","RHOB",'JUN',"BHLHE40","CEBPD","INHBA","IL6","JUNB","SOCS3",
                    "EGR2","SGK1","NR4A1","ZFP36","GADD45B","EGR1","DUSP1","FKBP5","CCN1","FOS","ATF3","SERPINE1",
                    'APOD',"CXCL14","IGHG4","CCN2","IGHG2","IGKG","IGHG3")

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

