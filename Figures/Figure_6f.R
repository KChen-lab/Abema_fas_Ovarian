## Source data of this plot is in supplementary files of the paper.

set.seed(2023)

library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

mat.z = read.csv("path-to-folder/Figure_6_f.csv")
meta_data = read.csv("path-to-folder/Figure_6_f_metadata.csv")

col_fun = colorRamp2(c(-1,0,1), c("purple", "black", "yellow3"))

# Create an initial data-frame of the annotation that we want to use
ann <- data.frame(
  Patients = meta_data$Patient_ID,
  stringsAsFactors = FALSE)

# create the colour mapping
colours <- list(
  Patients = c("8"="#fbb463","13"="#8ed1c6"))


# now create the ComplexHeatmap annotation object
###############

colAnn <- HeatmapAnnotation(
  df = ann,
  which = 'row', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'black', # default colour for any NA values in the annotation data-frame, 'ann'
  col = colours,
  
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Patients = list(
      nrow = 4, # number of rows across which the legend will be arranged
      title = 'Sample ID',
      title_position = 'topleft',
      
      legend_direction = 'horizontal',
      title_gp = gpar(fontsize = 18, fontface = 'bold'),
      labels_gp = gpar(fontsize = 18))))
    
    


split <- factor(ann$Patients, levels=c("8","13"))
ht_list = Heatmap(t(mat.z),name="z score", col=col_fun,cluster_columns = T,cluster_rows = T,
                  heatmap_legend_param = list(
                    color_bar = 'continuous',
                    title_position = 'topleft',
                    title_gp=gpar(fontsize = 18, fontface = 'bold'),
                    labels_gp=gpar(fontsize = 18, fontface = 'bold')),
                  clustering_distance_columns = "euclidean",
                  clustering_distance_rows = "euclidean",
                  row_names_side = "right",left_annotation = colAnn,
                  show_row_names = F, show_column_names = T,
                  row_names_gp = gpar(fontsize=18),
                  column_names_rot = 45,
                  column_title_gp =gpar(fontsize = 18),
                  row_split = split,
                  column_title = "Significant differential gene expressions in Tumor Pre treatment between short and long PFS")
ht_list
heatmap <- draw(ht_list,
                heatmap_legend_side = 'bottom',
                annotation_legend_side = 'bottom',
                row_sub_title_side = 'left')

heatmap

