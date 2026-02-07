###############################################################################
### Distance based Cell type Enrichment analysis
### Author : Jesuchristopher Joseph (jesu.joseph@duke.edu)
### Date : 
### Environment : R 4.4.1
### Packages  : dplyr 1.1.4, tidyr 1.3.1, spatstat 3.3-3, svMisc 1.4.3
###############################################################################


library(dplyr)
library(tidyr)
library(spatstat)
library(svMisc)


XeniummetdataLabelled <- (read.csv("C:/tools/Spatial_Analysis/tutorial_data/NicheAnalysis/Integrated/046MetadataUniform_Integrated_73125.csv", header = TRUE, sep = ","))

head(XeniummetdataLabelled)

print(unique(sort(XeniummetdataLabelled$leiden)))


# [1] "Adipocytes"                "B Cells"                   "Basal/Myoepithelial Cells" "CD8 T Cells"              
# [5] "Dendritic Cells"           "Endothelial Cells"         "Fibroblasts"               "Macrophages"              
# [9] "Mast Cells"                "Plasma Cells"              "Plasmacytoid DCs"          "Stroma"                   
# [13] "T Cells"                   "Tregs"                     "Tumor Cells"              
 

Tumor_Cells <- filter(XeniummetdataLabelled, leiden == "Tumor Cells")

Adipocytes <- filter(XeniummetdataLabelled, leiden == "Adipocytes")

Bcells <- filter(XeniummetdataLabelled, leiden == "B Cells")

Basal_MyoepitelialCells <- filter(XeniummetdataLabelled, leiden == "Basal/Myoepithelial Cells")

CD8Tcells <- filter(XeniummetdataLabelled, leiden == "CD8 T Cells")

Dendritic_cells <- filter(XeniummetdataLabelled, leiden == "Dendritic Cells")

Endothelial_cells <- filter(XeniummetdataLabelled, leiden == "Endothelial Cells")

Fibroblasts <- filter(XeniummetdataLabelled, leiden == "Fibroblasts")

Macrophages <- filter(XeniummetdataLabelled, leiden == "Macrophages")

Mast_cells <- filter(XeniummetdataLabelled, leiden == "Mast Cells")

Plasma_cells <- filter(XeniummetdataLabelled, leiden == "Plasma Cells")

Plasmacytoid_DCs <- filter(XeniummetdataLabelled, leiden == "Plasmacytoid DCs")

TCells <- filter(XeniummetdataLabelled, leiden == "T Cells")

Tregs <- filter(XeniummetdataLabelled, leiden == "Tregs")

Stroma <- filter(XeniummetdataLabelled, leiden == "Stroma")


threshold <- 26
##################################  Near Cells  #####################################

#Adipocytes

# Calculate pairwise distances
pairwise_distancesTumvsadipo <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Adipocytes$X, Adipocytes$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsadipo < threshold, arr.ind = TRUE)

indicesadipo <- as.data.frame(indices)

# Print the indices
print(indicesadipo)


indicesadipocol <- unique(indicesadipo$col)


AdipoNEAR <- data.frame(matrix(NA, nrow = length(indicesadipocol), ncol = dim(Adipocytes)[2]))


j=1
for (i in 1:length(indicesadipocol)) 
{
  AdipoNEAR[i,] <- Adipocytes[indicesadipocol[j],]
  j=j+1
  progress(i)
}


colnames(AdipoNEAR) <- colnames(XeniummetdataLabelled)


###############################################################################

#Bcells

# Calculate pairwise distances
pairwise_distancesTumvsbcells <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Bcells$X, Bcells$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsbcells < threshold, arr.ind = TRUE)

indicesbcells <- as.data.frame(indices)

# Print the indices
print(indicesbcells)


indicesbcellscol <- unique(indicesbcells$col)


BcellsNEAR <- data.frame(matrix(NA, nrow = length(indicesbcellscol), ncol = dim(Bcells)[2]))


j=1
for (i in 1:length(indicesbcellscol)) 
{
  BcellsNEAR[i,] <- Bcells[indicesbcellscol[j],]
  j=j+1
  progress(i)
}


colnames(BcellsNEAR) <- colnames(XeniummetdataLabelled)



###############################################################################

#Basal_MyoepitelialCells

# Calculate pairwise distances
pairwise_distancesTumvsbasal <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Basal_MyoepitelialCells$X, Basal_MyoepitelialCells$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsbasal < threshold, arr.ind = TRUE)

indicesbasal <- as.data.frame(indices)

# Print the indices
print(indicesbasal)


indicesbasalcol <- unique(indicesbasal$col)


BasalNEAR <- data.frame(matrix(NA, nrow = length(indicesbasalcol), ncol = dim(Basal_MyoepitelialCells)[2]))


j=1
for (i in 1:length(indicesbasalcol)) 
{
  BasalNEAR[i,] <- Basal_MyoepitelialCells[indicesbasalcol[j],]
  j=j+1
  progress(i)
}


colnames(BasalNEAR) <- colnames(XeniummetdataLabelled)



##############################################################################


#CD8

# Calculate pairwise distances
pairwise_distancesTumvsCD8 <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, CD8Tcells$X, CD8Tcells$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsCD8 < threshold, arr.ind = TRUE)

indicescd8 <- as.data.frame(indices)

# Print the indices
print(indicescd8)


indicescd8col <- unique(indicescd8$col)


CD8NEAR <- data.frame(matrix(NA, nrow = length(indicescd8col), ncol = dim(CD8Tcells)[2]))


j=1
for (i in 1:length(indicescd8col)) 
{
  CD8NEAR[i,] <- CD8Tcells[indicescd8col[j],]
  j=j+1
  progress(i)
}


colnames(CD8NEAR) <- colnames(XeniummetdataLabelled)


###############################################################################

#Dendritic_cells

# Calculate pairwise distances
pairwise_distancesTumvsdendritic <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Dendritic_cells$X, Dendritic_cells$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsdendritic < threshold, arr.ind = TRUE)

indicesdendritic <- as.data.frame(indices)

# Print the indices
print(indicesdendritic)


indicesdendriticcol <- unique(indicesdendritic$col)


DendriticNEAR <- data.frame(matrix(NA, nrow = length(indicesdendriticcol), ncol = dim(Dendritic_cells)[2]))


j=1
for (i in 1:length(indicesdendriticcol)) 
{
  DendriticNEAR[i,] <- Dendritic_cells[indicesdendriticcol[j],]
  j=j+1
  progress(i)
}


colnames(DendriticNEAR) <- colnames(XeniummetdataLabelled)


###############################################################################

#Endothelial_cells

# Calculate pairwise distances
pairwise_distancesTumvsEndo <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Endothelial_cells$X, Endothelial_cells$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsEndo < threshold, arr.ind = TRUE)

indicesendo <- as.data.frame(indices)

# Print the indices
print(indicesendo)


indicesendocol <- unique(indicesendo$col)


EndoNEAR <- data.frame(matrix(NA, nrow = length(indicesendocol), ncol = dim(Endothelial_cells)[2]))


j=1
for (i in 1:length(indicesendocol)) 
{
  EndoNEAR[i,] <- Endothelial_cells[indicesendocol[j],]
  j=j+1
  progress(i)
}


colnames(EndoNEAR) <- colnames(XeniummetdataLabelled)


###############################################################################

#Fibroblasts

# Calculate pairwise distances
pairwise_distancesTumvsFibro <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Fibroblasts$X, Fibroblasts$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsFibro < threshold, arr.ind = TRUE)

indicesfibro<- as.data.frame(indices)

# Print the indices
print(indicesfibro)


indicesfibrocol <- unique(indicesfibro$col)


FibroNEAR <- data.frame(matrix(NA, nrow = length(indicesfibrocol), ncol = dim(Fibroblasts)[2]))


j=1
for (i in 1:length(indicesfibrocol)) 
{
  FibroNEAR[i,] <- Fibroblasts[indicesfibrocol[j],]
  j=j+1
  progress(i)
}


colnames(FibroNEAR) <- colnames(XeniummetdataLabelled)


################################################################################

#Macrophages

# Calculate pairwise distances
pairwise_distancesTumvsMacrophages <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Macrophages$X, Macrophages$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsMacrophages < threshold, arr.ind = TRUE)

indicesMacrophages <- as.data.frame(indices)

# Print the indices
print(indicesMacrophages)


indicesMacrophagescol <- unique(indicesMacrophages$col)


MacrophagesNEAR <- data.frame(matrix(NA, nrow = length(indicesMacrophagescol), ncol = dim(Macrophages)[2]))


j=1
for (i in 1:length(indicesMacrophagescol)) 
{
  MacrophagesNEAR[i,] <- Macrophages[indicesMacrophagescol[j],]
  j=j+1
  progress(i)
}


colnames(MacrophagesNEAR) <- colnames(XeniummetdataLabelled)


###############################################################################

#Mast_cells

# Calculate pairwise distances
pairwise_distancesTumvsMast <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Mast_cells$X, Mast_cells$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsMast < threshold, arr.ind = TRUE)

indicesmast <- as.data.frame(indices)

# Print the indices
print(indicesmast)


indicesmastcol <- unique(indicesmast$col)


MastNEAR <- data.frame(matrix(NA, nrow = length(indicesmastcol), ncol = dim(Mast_cells)[2]))


j=1
for (i in 1:length(indicesmastcol)) 
{
  MastNEAR[i,] <- Mast_cells[indicesmastcol[j],]
  j=j+1
  progress(i)
}


colnames(MastNEAR) <- colnames(XeniummetdataLabelled)



###############################################################################

#Plasma_cells

# Calculate pairwise distances
pairwise_distancesTumvsPlasma <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Plasma_cells$X, Plasma_cells$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsPlasma < threshold, arr.ind = TRUE)

indicesplasma <- as.data.frame(indices)

# Print the indices
print(indicesplasma)


indicesplasmacol <- unique(indicesplasma$col)


PlasmaNEAR <- data.frame(matrix(NA, nrow = length(indicesplasmacol), ncol = dim(Plasma_cells)[2]))


j=1
for (i in 1:length(indicesplasmacol)) 
{
  PlasmaNEAR[i,] <- Plasma_cells[indicesplasmacol[j],]
  j=j+1
  progress(i)
}


colnames(PlasmaNEAR) <- colnames(XeniummetdataLabelled)


###############################################################################

#Plasmacytoid_DCs

# Calculate pairwise distances
pairwise_distancesTumvsPlasmacytoid <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Plasmacytoid_DCs$X, Plasmacytoid_DCs$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsPlasmacytoid < threshold, arr.ind = TRUE)

indicesplasmacytoid <- as.data.frame(indices)

# Print the indices
print(indicesplasmacytoid)


indicesplasmacytoidcol <- unique(indicesplasmacytoid$col)


PlasmacytoidNEAR <- data.frame(matrix(NA, nrow = length(indicesplasmacytoidcol), ncol = dim(Plasmacytoid_DCs)[2]))


j=1
for (i in 1:length(indicesplasmacytoidcol)) 
{
  PlasmacytoidNEAR[i,] <- Plasmacytoid_DCs[indicesplasmacytoidcol[j],]
  j=j+1
  progress(i)
}


colnames(PlasmacytoidNEAR) <- colnames(XeniummetdataLabelled)


################################################################################

#Tcells

# Calculate pairwise distances
pairwise_distancesTumvsTcells <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, TCells$X, TCells$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsTcells < threshold, arr.ind = TRUE)

indicesTcells <- as.data.frame(indices)

# Print the indices
print(indicesTcells)


indicesTcellscol <- unique(indicesTcells$col)


TcellsNEAR <- data.frame(matrix(NA, nrow = length(indicesTcellscol), ncol = dim(TCells)[2]))


j=1
for (i in 1:length(indicesTcellscol)) 
{
  TcellsNEAR[i,] <- TCells[indicesTcellscol[j],]
  j=j+1
  progress(i)
}


colnames(TcellsNEAR) <- colnames(XeniummetdataLabelled)


################################################################################

#Tregs

# Calculate pairwise distances
pairwise_distancesTumvsTregs <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Tregs$X, Tregs$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsTregs < threshold, arr.ind = TRUE)

indicesTregs <- as.data.frame(indices)

# Print the indices
print(indicesTregs)


indicesTregscol <- unique(indicesTregs$col)


TregsNEAR <- data.frame(matrix(NA, nrow = length(indicesTregscol), ncol = dim(Tregs)[2]))


j=1
for (i in 1:length(indicesTregscol)) 
{
  TregsNEAR[i,] <- Tregs[indicesTregscol[j],]
  j=j+1
  progress(i)
}


colnames(TregsNEAR) <- colnames(XeniummetdataLabelled)



###############################################################################
#Stroma

# Calculate pairwise distances
pairwise_distancesTumvsStroma <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Stroma$X, Stroma$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsStroma < threshold, arr.ind = TRUE)

indicesStroma <- as.data.frame(indices)

# Print the indices
print(indicesStroma)


indicesStromacol <- unique(indicesStroma$col)


StromaNEAR <- data.frame(matrix(NA, nrow = length(indicesStromacol), ncol = dim(Stroma)[2]))


j=1
for (i in 1:length(indicesStromacol)) 
{
  StromaNEAR[i,] <- Stroma[indicesStromacol[j],]
  j=j+1
  progress(i)
}


colnames(StromaNEAR) <- colnames(XeniummetdataLabelled)



###############################################################################

combined1 <- rbind(Tumor_Cells,AdipoNEAR,BcellsNEAR,BasalNEAR,CD8NEAR,DendriticNEAR,EndoNEAR,FibroNEAR,MacrophagesNEAR,MastNEAR,PlasmaNEAR,PlasmacytoidNEAR,TcellsNEAR,TregsNEAR,StromaNEAR )


write.csv(combined1, "C:\\Users\\jj374\\Documents\\RCode\\Brestcancer\\Integrated\\Xeniummetadata_Distance_25_046_81425.csv",row.names=FALSE)


##############################################################################




threshold <- 75
##################################  Far Cells  #################################

#Adipocytes

# Calculate pairwise distances
pairwise_distancesTumvsadipo <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Adipocytes$X, Adipocytes$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsadipo < threshold, arr.ind = TRUE)

indicesadipo <- as.data.frame(indices)

# Print the indices
print(indicesadipo)


indicesadipocol <- unique(indicesadipo$col)


AdipoFar <- data.frame(matrix(NA, nrow = length(indicesadipocol), ncol = dim(Adipocytes)[2]))


j=1
for (i in 1:length(indicesadipocol)) 
{
  AdipoFar[i,] <- Adipocytes[indicesadipocol[j],]
  j=j+1
  progress(i)
}


colnames(AdipoFar) <- colnames(XeniummetdataLabelled)


###############################################################################

#Bcells

# Calculate pairwise distances
pairwise_distancesTumvsbcells <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Bcells$X, Bcells$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsbcells < threshold, arr.ind = TRUE)

indicesbcells <- as.data.frame(indices)

# Print the indices
print(indicesbcells)


indicesbcellscol <- unique(indicesbcells$col)


BcellsFar <- data.frame(matrix(NA, nrow = length(indicesbcellscol), ncol = dim(Bcells)[2]))


j=1
for (i in 1:length(indicesbcellscol)) 
{
  BcellsFar[i,] <- Bcells[indicesbcellscol[j],]
  j=j+1
  progress(i)
}


colnames(BcellsFar) <- colnames(XeniummetdataLabelled)


###############################################################################

#Basal_MyoepitelialCells

# Calculate pairwise distances
pairwise_distancesTumvsbasal <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Basal_MyoepitelialCells$X, Basal_MyoepitelialCells$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsbasal < threshold, arr.ind = TRUE)

indicesbasal <- as.data.frame(indices)

# Print the indices
print(indicesbasal)


indicesbasalcol <- unique(indicesbasal$col)


BasalFar <- data.frame(matrix(NA, nrow = length(indicesbasalcol), ncol = dim(Basal_MyoepitelialCells)[2]))


j=1
for (i in 1:length(indicesbasalcol)) 
{
  BasalFar[i,] <- Basal_MyoepitelialCells[indicesbasalcol[j],]
  j=j+1
  progress(i)
}


colnames(BasalFar) <- colnames(XeniummetdataLabelled)



##############################################################################


#CD8

# Calculate pairwise distances
pairwise_distancesTumvsCD8 <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, CD8Tcells$X, CD8Tcells$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsCD8 < threshold, arr.ind = TRUE)

indicescd8 <- as.data.frame(indices)

# Print the indices
print(indicescd8)


indicescd8col <- unique(indicescd8$col)


CD8Far <- data.frame(matrix(NA, nrow = length(indicescd8col), ncol = dim(CD8Tcells)[2]))


j=1
for (i in 1:length(indicescd8col)) 
{
  CD8Far[i,] <- CD8Tcells[indicescd8col[j],]
  j=j+1
  progress(i)
}


colnames(CD8Far) <- colnames(XeniummetdataLabelled)


###############################################################################

#Dendritic_cells

# Calculate pairwise distances
pairwise_distancesTumvsdendritic <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Dendritic_cells$X, Dendritic_cells$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsdendritic < threshold, arr.ind = TRUE)

indicesdendritic <- as.data.frame(indices)

# Print the indices
print(indicesdendritic)


indicesdendriticcol <- unique(indicesdendritic$col)


DendriticFar <- data.frame(matrix(NA, nrow = length(indicesdendriticcol), ncol = dim(Dendritic_cells)[2]))


j=1
for (i in 1:length(indicesdendriticcol)) 
{
  DendriticFar[i,] <- Dendritic_cells[indicesdendriticcol[j],]
  j=j+1
  progress(i)
}


colnames(DendriticFar) <- colnames(XeniummetdataLabelled)


###############################################################################

#Endothelial_cells

# Calculate pairwise distances
pairwise_distancesTumvsEndo <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Endothelial_cells$X, Endothelial_cells$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsEndo < threshold, arr.ind = TRUE)

indicesendo <- as.data.frame(indices)

# Print the indices
print(indicesendo)


indicesendocol <- unique(indicesendo$col)


EndoFar <- data.frame(matrix(NA, nrow = length(indicesendocol), ncol = dim(Endothelial_cells)[2]))

j=1
for (i in 1:length(indicesendocol)) 
{
  EndoFar[i,] <- Endothelial_cells[indicesendocol[j],]
  j=j+1
  progress(i)
}


colnames(EndoFar) <- colnames(XeniummetdataLabelled)


###############################################################################

#Fibroblasts

# Calculate pairwise distances
pairwise_distancesTumvsFibro <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Fibroblasts$X, Fibroblasts$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsFibro < threshold, arr.ind = TRUE)

indicesfibro<- as.data.frame(indices)

# Print the indices
print(indicesfibro)


indicesfibrocol <- unique(indicesfibro$col)


FibroFar <- data.frame(matrix(NA, nrow = length(indicesfibrocol), ncol = dim(Fibroblasts)[2]))


j=1
for (i in 1:length(indicesfibrocol)) 
{
  FibroFar[i,] <- Fibroblasts[indicesfibrocol[j],]
  j=j+1
  progress(i)
}


colnames(FibroFar) <- colnames(XeniummetdataLabelled)


################################################################################

#Macrophages

# Calculate pairwise distances
pairwise_distancesTumvsMacrophages <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Macrophages$X, Macrophages$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsMacrophages < threshold, arr.ind = TRUE)

indicesMacrophages <- as.data.frame(indices)

# Print the indices
print(indicesMacrophages)


indicesMacrophagescol <- unique(indicesMacrophages$col)


MacrophagesFar <- data.frame(matrix(NA, nrow = length(indicesMacrophagescol), ncol = dim(Macrophages)[2]))


j=1
for (i in 1:length(indicesMacrophagescol)) 
{
  MacrophagesFar[i,] <- Macrophages[indicesMacrophagescol[j],]
  j=j+1
  progress(i)
}


colnames(MacrophagesFar) <- colnames(XeniummetdataLabelled)

################################################################################


#Mast_cells

# Calculate pairwise distances
pairwise_distancesTumvsMast <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Mast_cells$X, Mast_cells$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsMast < threshold, arr.ind = TRUE)

indicesmast <- as.data.frame(indices)

# Print the indices
print(indicesmast)


indicesmastcol <- unique(indicesmast$col)


MastFar <- data.frame(matrix(NA, nrow = length(indicesmastcol), ncol = dim(Mast_cells)[2]))


j=1
for (i in 1:length(indicesmastcol)) 
{
  MastFar[i,] <- Mast_cells[indicesmastcol[j],]
  j=j+1
  progress(i)
}


colnames(MastFar) <- colnames(XeniummetdataLabelled)


###############################################################################

#Plasma_cells

# Calculate pairwise distances
pairwise_distancesTumvsPlasma <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Plasma_cells$X, Plasma_cells$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsPlasma < threshold, arr.ind = TRUE)

indicesplasma <- as.data.frame(indices)

# Print the indices
print(indicesplasma)


indicesplasmacol <- unique(indicesplasma$col)


PlasmaFar <- data.frame(matrix(NA, nrow = length(indicesplasmacol), ncol = dim(Plasma_cells)[2]))


j=1
for (i in 1:length(indicesplasmacol)) 
{
  PlasmaFar[i,] <- Plasma_cells[indicesplasmacol[j],]
  j=j+1
  progress(i)
}


colnames(PlasmaFar) <- colnames(XeniummetdataLabelled)



###############################################################################

#Plasmacytoid_DCs

# Calculate pairwise distances
pairwise_distancesTumvsPlasmacytoid <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Plasmacytoid_DCs$X, Plasmacytoid_DCs$Y)



# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsPlasmacytoid < threshold, arr.ind = TRUE)

indicesplasmacytoid <- as.data.frame(indices)

# Print the indices
print(indicesplasmacytoid)


indicesplasmacytoidcol <- unique(indicesplasmacytoid$col)


PlasmacytoidFar <- data.frame(matrix(NA, nrow = length(indicesplasmacytoidcol), ncol = dim(Plasmacytoid_DCs)[2]))


j=1
for (i in 1:length(indicesplasmacytoidcol)) 
{
  PlasmacytoidFar[i,] <- Plasmacytoid_DCs[indicesplasmacytoidcol[j],]
  j=j+1
  progress(i)
}


colnames(PlasmacytoidFar) <- colnames(XeniummetdataLabelled)



################################################################################
#Tcells

# Calculate pairwise distances
pairwise_distancesTumvsTcells <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, TCells$X, TCells$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsTcells < threshold, arr.ind = TRUE)

indicesTcells <- as.data.frame(indices)

# Print the indices
print(indicesTcells)


indicesTcellscol <- unique(indicesTcells$col)


TcellsFar <- data.frame(matrix(NA, nrow = length(indicesTcellscol), ncol = dim(TCells)[2]))


j=1
for (i in 1:length(indicesTcellscol)) 
{
  TcellsFar[i,] <- TCells[indicesTcellscol[j],]
  j=j+1
  progress(i)
}


colnames(TcellsFar) <- colnames(XeniummetdataLabelled)



################################################################################

#Tregs

# Calculate pairwise distances
pairwise_distancesTumvsTregs <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Tregs$X, Tregs$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsTregs < threshold, arr.ind = TRUE)

indicesTregs <- as.data.frame(indices)

# Print the indices
print(indicesTregs)


indicesTregscol <- unique(indicesTregs$col)


TregsFar <- data.frame(matrix(NA, nrow = length(indicesTregscol), ncol = dim(Tregs)[2]))


j=1
for (i in 1:length(indicesTregscol)) 
{
  TregsFar[i,] <- Tregs[indicesTregscol[j],]
  j=j+1
  progress(i)
}


colnames(TregsFar) <- colnames(XeniummetdataLabelled)



################################################################################


#Stroma

# Calculate pairwise distances
pairwise_distancesTumvsStroma <- crossdist(Tumor_Cells$X, Tumor_Cells$Y, Stroma$X, Stroma$Y)


# Find indices where distance is less than the threshold
indices <- which(pairwise_distancesTumvsStroma < threshold, arr.ind = TRUE)

indicesStroma <- as.data.frame(indices)

# Print the indices
print(indicesStroma)


indicesStromacol <- unique(indicesStroma$col)


StromaFar <- data.frame(matrix(NA, nrow = length(indicesStromacol), ncol = dim(Stroma)[2]))


j=1
for (i in 1:length(indicesStromacol)) 
{
  StromaFar[i,] <- Stroma[indicesStromacol[j],]
  j=j+1
  progress(i)
}


colnames(StromaFar) <- colnames(XeniummetdataLabelled)

################################################################################


combined2 <- rbind(Tumor_Cells,AdipoFar,BcellsFar,BasalFar,CD8Far,DendriticFar,EndoFar,FibroFar,MacrophagesFar,MastFar,PlasmaFar,PlasmacytoidFar,TcellsFar,TregsFar,StromaFar )


write.csv(combined2, "C:\\Users\\jj374\\Documents\\RCode\\Brestcancer\\Integrated\\Xeniummetadata_Distance_75_046_81425.csv",row.names=FALSE)

#####################################Finding Far Cell Continue (cell > 75 um) ###################

Filteredcell75 <- (read.csv("C:\\Users\\jj374\\Documents\\RCode\\Brestcancer\\Integrated\\Xeniummetadata_Distance_75_046_81425.csv", header = TRUE, sep = ","))


Filteredcell75_Colnames = data.frame(Column1 = Filteredcell75$Cell_ID)

Allcells_filter = XeniummetdataLabelled

Cell_ANN_Colnames= data.frame(Column2 = Allcells_filter$Cell_ID)

matches_indices <- match(Filteredcell75_Colnames$Column1, Cell_ANN_Colnames$Column2 )


array1 <- 1:dim(Allcells_filter)[1]
array2 <- matches_indices


diff <- setdiff(array1, array2)


Farcells <- data.frame(matrix(NA, nrow = length(diff), ncol = dim(Allcells_filter)[2]))


for (i in 1: length(diff))
  {
  Farcells[i,] <- Allcells_filter[diff[i],]
    
    }

colnames(Farcells) <- colnames(XeniummetdataLabelled)

Farcells ['distance'] <- "Far"


Farcells <- Farcells[order(Farcells$Cell_ID),]


write.csv(Farcells, "C:\\Users\\jj374\\Documents\\RCode\\Brestcancer\\Integrated\\Xeniummetdata_Fardistance_046_81425.csv",row.names=FALSE)


#####################################Selecting Cells Between 25 um and 75um , the Invasive margin###################

Filteredcell25 <- (read.csv("C:\\Users\\jj374\\Documents\\RCode\\Brestcancer\\Integrated\\Xeniummetadata_Distance_25_046_81425.csv", header = TRUE, sep = ","))


Filteredcell25_Colnames = data.frame(Column1 = Filteredcell25$Cell_ID)
Filteredcell75_Colnames = data.frame(Column2 = Filteredcell75$Cell_ID)


matches_indices <- match(Filteredcell25_Colnames$Column1, Filteredcell75_Colnames$Column2 )


array1 <- 1:dim(Filteredcell75_Colnames)[1]
array2 <- matches_indices


diff <- setdiff(array1, array2)


Margincells <- data.frame(matrix(NA, nrow = length(diff), ncol = dim(Filteredcell75)[2]))


for (i in 1: length(diff))
{
  Margincells[i,] <- Filteredcell75[diff[i],]
  
}

colnames(Margincells) <- colnames(XeniummetdataLabelled)

Margincells ['distance'] <- "Margin"


Margincells <- Margincells[order(Margincells$Cell_ID),]


write.csv(Margincells, "C:\\Users\\jj374\\Documents\\RCode\\Brestcancer\\Integrated\\Xeniummetdata_Margindistance_046_81425.csv",row.names=FALSE)



Nearcells <- Filteredcell25
Nearcells ['distance'] <- "Near"

write.csv(Nearcells, "C:\\Users\\jj374\\Documents\\RCode\\Brestcancer\\Integrated\\Xeniummetdata_Neardistance_046_81425.csv",row.names=FALSE)


combined3 <- rbind(Nearcells,Margincells,Farcells )

combined3 <- combined3[order(combined3$Cell_ID),]


write.csv(combined3, "C:\\Users\\jj374\\Documents\\RCode\\Brestcancer\\Integrated\\XeniummetdataCombined_Integrated_046_12926.csv",row.names=FALSE)





