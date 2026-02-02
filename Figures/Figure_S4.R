## Spatial overlay plots

set.seed(2025)
options(java.parameters = "-Xmx8000m")
##https://github.com/Nanostring-Biostats/SpatialOmicsOverlay

library(SpatialOmicsOverlay)
library(GeomxTools)

# annotation file
## This data is labworksheet file. It is under data folder of github website.
anno <-list.files(path="/path-to-file", pattern="*.txt", full.names=T)

## Normalized data 
## This file is in the GEO submission
Ovary <- read.csv("Ovary_cancer_DSP_WTA_normalized_data_corrected.csv",row.names=1) 

## Create objects to visualize gene expressions on top of ome-tiff files
## Patient 13 Pre treatment
tiff_1 <-list.files(path="path-to-ometiff-files", pattern="Lee4T1S1.ome.tiff", full.names=T)

Lee4T1S1 <- readSpatialOverlay(ometiff = tiff_1, annots = anno, 
                               slideName = "Lee4T1S1", image = T,
                               saveFile = FALSE, outline = FALSE)
Anno_1 <- readLabWorksheet(lw = anno, slideName = "Lee4T1S1")




## Selected genes
Lee4T1S1 <- addPlottingFactor(overlay = Lee4T1S1, annots = Ovary, 
                              plottingFactor = "CDK4")

Lee4T1S1 <- addPlottingFactor(overlay = Lee4T1S1, annots = Ovary, 
                                   plottingFactor = "GREB1")

Lee4T1S1 <- addPlottingFactor(overlay = Lee4T1S1, annots = Ovary, 
                              plottingFactor = "NDUFA6")
Lee4T1S1 <- addPlottingFactor(overlay = Lee4T1S1, annots = Ovary, 
                              plottingFactor = "TOP2A")

Lee4T1S1<- addPlottingFactor(overlay = Lee4T1S1, annots = Ovary, 
                              plottingFactor = "UQCR11")

Lee4T1S1<- addPlottingFactor(overlay = Lee4T1S1, annots = Ovary, 
                                  plottingFactor = "FH")
Lee4T1S1<- addPlottingFactor(overlay = Lee4T1S1, annots = Ovary, 
                                  plottingFactor = "ATP2B4")

samps_1 <- Anno_1$Sample_ID[Anno_1$segment == "Tumor" & 
                              Anno_1$slide.name == slideName(Lee4T1S1)]
# Crop object only for Tumor cells
Lee4T1S1_crop <- cropSamples(overlay = Lee4T1S1, sampleIDs = samps_1, sampsOnly = TRUE)


### Pt8 Pre treatment


tiff_2 <-list.files(path="path-to-ometiff-files", pattern="Lee1T1S4.ome.tiff", full.names=T)

Lee1T1S4 <- readSpatialOverlay(ometiff = tiff_2, annots = anno, 
                               slideName = "Lee1T1S4", image = T,
                               saveFile = FALSE, outline = FALSE)
Anno_2 <- readLabWorksheet(lw = anno, slideName = "Lee1T1S4")




## Selected genes
Lee1T1S4<- addPlottingFactor(overlay = Lee1T1S4, annots = Ovary, 
                              plottingFactor = "CDK4")

Lee1T1S4 <- addPlottingFactor(overlay = Lee1T1S4, annots = Ovary, 
                              plottingFactor = "GREB1")

Lee1T1S4 <- addPlottingFactor(overlay = Lee1T1S4, annots = Ovary, 
                              plottingFactor = "NDUFA6")
Lee1T1S4 <- addPlottingFactor(overlay = Lee1T1S4, annots = Ovary, 
                              plottingFactor = "TOP2A")

Lee1T1S4<- addPlottingFactor(overlay = Lee1T1S4, annots = Ovary, 
                             plottingFactor = "UQCR11")

Lee1T1S4<- addPlottingFactor(overlay = Lee1T1S4, annots = Ovary, 
                             plottingFactor = "FH")
Lee1T1S4<- addPlottingFactor(overlay = Lee1T1S4, annots = Ovary, 
                             plottingFactor = "ATP2B4")

samps_2 <- Anno_2$Sample_ID[Anno_2$segment == "Tumor" & 
                              Anno_2$slide.name == slideName(Lee1T1S4)]
# Crop object only for Tumor cells
Lee1T1S4_crop <- cropSamples(overlay = Lee1T1S4, sampleIDs = samps_2, sampsOnly = TRUE)

######## Pt11 Pre treatment


tiff_3 <-list.files(path="path-to-ometiff-files", pattern="Lee3T1S1.ome.tiff", full.names=T)

Lee3T1S1<- readSpatialOverlay(ometiff = tiff_3, annots = anno, 
                               slideName = "Lee3T1S1", image = T,
                               saveFile = FALSE, outline = FALSE)
Anno_3 <- readLabWorksheet(lw = anno, slideName = "Lee3T1S1")




## Selected genes
Lee3T1S1<- addPlottingFactor(overlay = Lee3T1S1, annots = Ovary, 
                             plottingFactor = "CDK4")

Lee3T1S1 <- addPlottingFactor(overlay = Lee3T1S1, annots = Ovary, 
                              plottingFactor = "GREB1")

Lee3T1S1 <- addPlottingFactor(overlay = Lee3T1S1, annots = Ovary, 
                              plottingFactor = "NDUFA6")
Lee3T1S1 <- addPlottingFactor(overlay = Lee3T1S1, annots = Ovary, 
                              plottingFactor = "TOP2A")

Lee3T1S1<- addPlottingFactor(overlay = Lee3T1S1, annots = Ovary, 
                             plottingFactor = "UQCR11")

Lee3T1S1<- addPlottingFactor(overlay = Lee3T1S1, annots = Ovary, 
                             plottingFactor = "FH")
Lee3T1S1<- addPlottingFactor(overlay = Lee3T1S1, annots = Ovary, 
                             plottingFactor = "ATP2B4")

samps_3 <- Anno_3$Sample_ID[Anno_3$segment == "Tumor" & 
                              Anno_3$slide.name == slideName(Lee3T1S1)]

# Crop object only for Tumor cells
Lee3T1S1_crop <- cropSamples(overlay = Lee3T1S1, sampleIDs = samps_3, sampsOnly = TRUE)

######################################

### generate images 
library(ggpubr)
library("scales")

#CDK4

#Pt 13 pre treatment
A <- plotSpatialOverlay(overlay = Lee4T1S1_crop, hiRes = T, colorBy = "CDK4", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  
  ggplot2::labs(title = "CDK4 Expression Pre treatment Pt 13 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(5.4, 7.6)),
                       guide = "colorbar", limits=c(5.4, 7.6))
A

pdf(file='Spatial_Overlay_PT_13_pretreatment_CDK4.pdf')
A
dev.off()

#Pt 8 pre treatment
B <- plotSpatialOverlay(overlay = Lee1T1S4_crop, hiRes = T, colorBy = "CDK4", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "CDK4 Expression Pre treatment Pt 8 (Short PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(5.4, 7.6)),
                       guide = "colorbar", limits=c(5.4, 7.6))
B

pdf(file='Spatial_Overlay_PT_8_pretreatment_CDK4.pdf')
B
dev.off()

#Pt 11 pre treatment
C <- plotSpatialOverlay(overlay = Lee3T1S1_crop, hiRes = T, colorBy = "CDK4", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "CDK4 Expression Pre treatment Pt 11 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(5.4, 7.6)),
                       guide = "colorbar", limits=c(5.4, 7.6))
C

pdf(file='Spatial_Overlay_PT_11_pretreatment_CDK4.pdf')
C
dev.off()

#GREB1

#Pt 13 pre treatment
A <- plotSpatialOverlay(overlay = Lee4T1S1_crop, hiRes = T, colorBy = "GREB1", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  #viridis::scale_color_viridis()+
  ggplot2::labs(title = "GREB1 Expression Pre treatment Pt 13 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(3.8, 6.7)),
                       guide = "colorbar", limits=c(3.8, 6.7))
A

pdf(file='Spatial_Overlay_PT_13_pretreatment_GREB1.pdf')
A
dev.off()

#Pt 8 pre treatment
B <- plotSpatialOverlay(overlay = Lee1T1S4_crop, hiRes = T, colorBy = "GREB1", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "GREB1 Expression Pre treatment Pt 8 (Short PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(3.8, 6.7)),
                       guide = "colorbar", limits=c(3.8, 6.7))
B

pdf(file='Spatial_Overlay_PT_8_pretreatment_GREB1.pdf')
B
dev.off()

#Pt 11 pre treatment
C <- plotSpatialOverlay(overlay = Lee3T1S1_crop, hiRes = T, colorBy = "GREB1", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "GREB1 Expression Pre treatment Pt 11 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(3.8, 6.7)),
                       guide = "colorbar", limits=c(3.8, 6.7))
C

pdf(file='Spatial_Overlay_PT_11_pretreatment_GREB1.pdf')
C
dev.off()

##NDUFA6

#Pt 13 pre treatment
A <- plotSpatialOverlay(overlay = Lee4T1S1_crop, hiRes = T, colorBy = "NDUFA6", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "NDUFA6 Expression Pre treatment Pt 13 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(4.3, 6.5)),
                       guide = "colorbar", limits=c(4.3, 6.5))
A

pdf(file='Spatial_Overlay_PT_13_pretreatment_NDUFA6.pdf')
A
dev.off()

#Pt 8 pre treatment
B <- plotSpatialOverlay(overlay = Lee1T1S4_crop, hiRes = T, colorBy = "NDUFA6", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "NDUFA6 Expression Pre treatment Pt 8 (Short PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(4.3, 6.5)),
                       guide = "colorbar", limits=c(4.3, 6.5))
B

pdf(file='Spatial_Overlay_PT_8_pretreatment_NDUFA6.pdf')
B
dev.off()

#Pt 11 pre treatment
C <- plotSpatialOverlay(overlay = Lee3T1S1_crop, hiRes = T, colorBy = "NDUFA6", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "NDUFA6 Expression Pre treatment Pt 11 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(4.3, 6.5)),
                       guide = "colorbar", limits=c(4.3, 6.5))
C

pdf(file='Spatial_Overlay_PT_11_pretreatment_NDUFA6.pdf')
C
dev.off()

## TOP2A

#Pt 13 pre treatment
A <- plotSpatialOverlay(overlay = Lee4T1S1_crop, hiRes = T, colorBy = "TOP2A", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "TOP2A Expression Pre treatment Pt 13 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(3.15, 6.5)),
                       guide = "colorbar", limits=c(3.15, 6.5))
A

pdf(file='Spatial_Overlay_PT_13_pretreatment_TOP2A.pdf')
A
dev.off()

#Pt 8 pre treatment
B <- plotSpatialOverlay(overlay = Lee1T1S4_crop, hiRes = T, colorBy = "TOP2A", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "TOP2A Expression Pre treatment Pt 8 (Short PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(3.15, 6.5)),
                       guide = "colorbar", limits=c(3.15, 6.5))
B

pdf(file='Spatial_Overlay_PT_8_pretreatment_TOP2A.pdf')
B
dev.off()

#Pt 11 pre treatment
C <- plotSpatialOverlay(overlay = Lee3T1S1_crop, hiRes = T, colorBy = "TOP2A", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "TOP2A Expression Pre treatment Pt 11 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(3.15, 6.5)),
                       guide = "colorbar", limits=c(3.15, 6.5))
C

pdf(file='Spatial_Overlay_PT_11_pretreatment_TOP2A.pdf')
C
dev.off()

##UQCR11

#Pt 13 pre treatment
A <- plotSpatialOverlay(overlay = Lee4T1S1_crop, hiRes = T, colorBy = "UQCR11", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "UQCR11 Expression Pre treatment Pt 13 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(5.6, 8.5)),
                       guide = "colorbar", limits=c(5.6, 8.5))
A

pdf(file='Spatial_Overlay_PT_13_pretreatment_UQCR11.pdf')
A
dev.off()

#Pt 8 pre treatment
B <- plotSpatialOverlay(overlay = Lee1T1S4_crop, hiRes = T, colorBy = "UQCR11", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "UQCR11 Expression Pre treatment Pt 8 (Short PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(5.6, 8.5)),
                       guide = "colorbar", limits=c(5.6, 8.5))
B

pdf(file='Spatial_Overlay_PT_8_pretreatment_UQCR11.pdf')
B
dev.off()

#Pt 11 pre treatment
C <- plotSpatialOverlay(overlay = Lee3T1S1_crop, hiRes = T, colorBy = "UQCR11", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "UQCR11 Expression Pre treatment Pt 11 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(5.6, 8.5)),
                       guide = "colorbar", limits=c(5.6, 8.5))
C

pdf(file='Spatial_Overlay_PT_11_pretreatment_UQCR11.pdf')
C
dev.off()

##FH

#Pt 13 pre treatment
A <- plotSpatialOverlay(overlay = Lee4T1S1_crop, hiRes = T, colorBy = "FH", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "FH Expression Pre treatment Pt 13 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(3.9, 7.1)),
                       guide = "colorbar", limits=c(3.9, 7.1))
A

pdf(file='Spatial_Overlay_PT_13_pretreatment_FH.pdf')
A
dev.off()

#Pt 8 pre treatment
B <- plotSpatialOverlay(overlay = Lee1T1S4_crop, hiRes = T, colorBy = "FH", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "FH Expression Pre treatment Pt 8 (Short PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(3.9, 7.1)),
                       guide = "colorbar", limits=c(3.9, 7.1))
B

pdf(file='Spatial_Overlay_PT_8_pretreatment_FH.pdf')
B
dev.off()

#Pt 11 pre treatment
C <- plotSpatialOverlay(overlay = Lee3T1S1_crop, hiRes = T, colorBy = "FH", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "FH Expression Pre treatment Pt 11 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(3.9, 7.1)),
                       guide = "colorbar", limits=c(3.9, 7.1))
C

pdf(file='Spatial_Overlay_PT_11_pretreatment_FH.pdf')
C
dev.off()

### ATP2B4
#Pt 13 pre treatment
A <- plotSpatialOverlay(overlay = Lee4T1S1_crop, hiRes = T, colorBy = "ATP2B4", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "ATP2B4 Expression Pre treatment Pt 13 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(5.9, 8.7)),
                       guide = "colorbar", limits=c(5.9, 8.7))
A

pdf(file='Spatial_Overlay_PT_13_pretreatment_ATP2B4.pdf')
A
dev.off()

#Pt 8 pre treatment
B <- plotSpatialOverlay(overlay = Lee1T1S4_crop, hiRes = T, colorBy = "ATP2B4", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "ATP2B4 Expression Pre treatment Pt 8 (Short PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(5.9, 8.7)),
                       guide = "colorbar", limits=c(5.9, 8.7))
B

pdf(file='Spatial_Overlay_PT_8_pretreatment_ATP2B4.pdf')
B
dev.off()

#Pt 11 pre treatment
C <- plotSpatialOverlay(overlay = Lee3T1S1_crop, hiRes = T, colorBy = "ATP2B4", 
                        scaleBarWidth = 0.3, scaleBarColor = "green",image = T,
                        fluorLegend = T) +
  ggplot2::labs(title = "ATP2B4 Expression Pre treatment Pt 11 (Long PFS)")+
  scale_fill_gradientn(colours = c("blue","white","red"), 
                       values = rescale(c(5.9, 8.7)),
                       guide = "colorbar", limits=c(5.9, 8.7))
C

pdf(file='Spatial_Overlay_PT_11_pretreatment_ATP2B4.pdf')
C
dev.off()

