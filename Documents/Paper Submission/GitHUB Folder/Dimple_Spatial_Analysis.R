###############################################################################
### Estimating pairwise distance metrics using DIMPLE
### Author : Jesuchristopher Joseph (jesu.joseph@duke.edu)
### Date : 
### Environment : R 4.4.2
### Packages  : DIMPLE 1.0.0, tidyverse 2.0.0, viridis 0.6.5, ggplot2 4.0.1
###############################################################################


# Loading Library
library(DIMPLE)
library(tidyverse)
library(viridis)
library("ggplot2")


# Import Metadata of sample 046 containing XY coordinates, Leiden cell type labels 

Cellannotataion <- (read.csv("C:/tools/Spatial_Analysis/tutorial_data/046metadataUniform32125.csv",sep=','))



cell_x_coords = Cellannotataion$x_centroid %>% as.numeric()
cell_y_coords = Cellannotataion$y_centroid %>% as.numeric()
cell_marks = Cellannotataion$leiden

slide_ids = 1 # default one slide at a time

for (i in 1:length(cell_x_coords)) 
{
  
  
  slide_ids[i] = ": Spatial plot"
 
}


lung_experiment = DIMPLE::new_MltplxExperiment(
  x = cell_x_coords,
  y = cell_y_coords,
  marks = cell_marks,
  slide_id = slide_ids,
  
)

lung_experiment


plot(lung_experiment[[1]]) + scale_y_reverse() + theme_classic() +
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        legend.key.size = unit(1.0, "cm"),legend.key.width = unit(0.5,"cm") , height = unit(1.5,"cm"),
        legend.text = element_text(size = 20),legend.title=element_text(size=20))

#ggsave("Spatialplot_046.png", dpi=300)  # Save Figure

lung_experiment = update_intensity(lung_experiment,
                                   ps = 10,
                                   bw = 30) 

plot(lung_experiment[[1]]$mltplx_intensity) + scale_y_reverse() + theme_classic() +
       theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        legend.text = element_text(size = 15),legend.title=element_text(size=15),
        plot.title = element_text(size = 15))

#ggsave("Intensityplot_046.png", dpi=300)


lung_experiment = update_dist(lung_experiment,
                              dist_metric = jsd)

plot_dist_matrix(lung_experiment[[1]]) + theme_classic() + 
                   theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=15), 
                         legend.text = element_text(size = 15),legend.title=element_text(size=15), 
                         plot.title = element_text(size = 15))  + geom_text(aes(label = round(dist, 2)),size = 5) +
                         scale_fill_viridis(direction = -1) + scale_color_viridis(direction = -1)
                 
#ggsave("Distanceplot_046.png", dpi=300)



