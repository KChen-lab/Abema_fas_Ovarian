
## Source data of this plot is in supplementary files of the paper.

set.seed(2023)
library("EnhancedVolcano")
library(ggrepel)
library(ggplot2)
symbol <- rownames(top.table_Tumor_PFS) 
top.table_Tumor_PFS$symbol <- symbol

                                                     
keyvals <- ifelse((top.table_Tumor_PFS$logFC < 0 )& (top.table_Tumor_PFS$adj.P.Val < 0.001), '#203569',
                  ifelse((top.table_Tumor_PFS$logFC < 0 )& (top.table_Tumor_PFS$adj.P.Val < 0.01), '#2572b5',
                         ifelse((top.table_Tumor_PFS$logFC < 0 )& (top.table_Tumor_PFS$adj.P.Val < 0.05), '#6baed6',
                                ifelse((top.table_Tumor_PFS$logFC < 0 )& (top.table_Tumor_PFS$adj.P.Val < 0.1), '#c7dbee',
                                       ifelse((top.table_Tumor_PFS$logFC< 0) & (top.table_Tumor_PFS$adj.P.Val< 0.05), '#f1f0f0',
                                              
                                                      ifelse((top.table_Tumor_PFS$logFC > 0 )& (top.table_Tumor_PFS$adj.P.Val < 0.001), '#be2029',
                                                            ifelse((top.table_Tumor_PFS$logFC > 0 )& (top.table_Tumor_PFS$adj.P.Val < 0.01), '#f05232',
                                                                  ifelse((top.table_Tumor_PFS$logFC > 0 )& (top.table_Tumor_PFS$adj.P.Val < 0.05), '#fbb14e',
                                                                        ifelse((top.table_Tumor_PFS$logFC > 0 )& (top.table_Tumor_PFS$adj.P.Val < 0.1), '#feec9f',
                                                     '#f1f0f0')))))))))

names(keyvals)[keyvals == '#203569'] <- 'p.adj < 0.001, Neg'
names(keyvals)[keyvals == '#2572b5'] <- 'p.adj < 0.01, Neg'
names(keyvals)[keyvals == '#6baed6'] <- 'p.adj < 0.05, Neg'
names(keyvals)[keyvals == '#c7dbee'] <- 'p.adj < 0.1, Neg'
names(keyvals)[keyvals == '#f1f0f0'] <- 'NS'
names(keyvals)[keyvals == '#be2029'] <- 'p.adj < 0.001, Pos'
names(keyvals)[keyvals == '#f05232'] <- 'p.adj < 0.01, Pos'
names(keyvals)[keyvals == '#fbb14e'] <- 'p.adj < 0.05, Pos'
names(keyvals)[keyvals == '#feec9f'] <- 'p.adj < 0.1, Pos'

Volcano <- EnhancedVolcano(top.table_Tumor_PFS,
                           lab =top.table_Tumor_PFS$symbol,
                    
                           selectLab=c("IDO1","CFAP157","FOXJ1","PIFO","HLA−F","IFI27","TAP1",'HLA−A','HLA−B','STAT1',"HLA−DRB1","B2M",'ISG15',
                           'CDK4','FH','NDUFA6','ATP2B4','GREB1','UQCR11','SDHD','TOP2A','CDK1','SULT2B1','TUBB'),
                           x = 'logFC',
                           y = 'adj.P.Val',
                           xlim = c(-4.8,6.6),
                           ylim=c(0,30),
                           ylab=bquote(~-Log[10]~Adjusted~italic(P)),
                           xlab = bquote(~Log~ 'FC'),
                           labSize = 5,
                           title = 'Differential expressed genes\nBetween short and long PFS at pre treatment in Tumor',
                           pCutoff = 0.05,
                           FCcutoff = 1.0,
                           colCustom = keyvals,
                           colAlpha = 1,
                           drawConnectors = T,
                           legendPosition='bottom')

Volcano
