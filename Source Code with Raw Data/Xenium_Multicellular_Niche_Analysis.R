###############################################################################
### Estimating Multicellular Niche using BuildNicheAssay
### Author : Jesuchristopher Joseph (jesu.joseph@duke.edu)
### Date : 
### Environment : R 4.4.1
### Packages  : Seurat 5.0.3, future 1.68.0,  ggplot2 4.0.1
###############################################################################


library(Seurat)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = +Inf)

library(future)
plan("multisession", workers = 5)

library(ggplot2)

### Input folder containing Xenium Output

xenium_input<- ReadXenium(
  data.dir = "C:/Users/jj374/Documents/RCode/Brestcancer/output-XETG00054__0006530__046__20230623__154919/",
  outs = "matrix",
  type = "centroids",
  mols.qv.threshold = 20
)


segmentations.data <- list(
  "centroids" = CreateCentroids(xenium_input$centroids)
  #"segmentation" = CreateSegmentation(data$segmentations)
)

coords <- CreateFOV(
  coords = segmentations.data,
  type = "centroids",
  molecules = NULL,
  assay = "Xenium"
)

xenium.obj <- CreateSeuratObject(counts = xenium_input$matrix[["Gene Expression"]], assay = "Xenium")

xenium.obj[["NegativeControlCodeword"]] <- CreateAssayObject(counts = xenium_input$matrix[["Negative Control Codeword"]])

xenium.obj[["NegativeControlProbe"]] <- CreateAssayObject(counts = xenium_input$matrix[["Negative Control Probe"]])

xenium.obj[["fov"]] <- coords



xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)


xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")

xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))

ElbowPlot(xenium.obj)


xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)

xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)

xenium.obj <- FindClusters(xenium.obj, resolution = 0.5)


DimPlot(xenium.obj)


plot <- DimPlot(xenium.obj, reduction = "umap") + NoLegend()
LabelClusters(plot = plot, id = "ident")



Cell_ID <- xenium.obj@images$fov@boundaries$centroids@cells
Cell_Coor <- xenium.obj@images$fov@boundaries$centroids@coords

#######################Read CellID  from the metadata ################


data_table  = as.data.frame (xenium.obj@meta.data)

data_table_Colnames= data.frame(Column2 = rownames(data_table))


xenium.obj@meta.data

###########################Reading Cell annotations###########################


ReadCellannotataion <- (read.csv("C:/tools/Spatial_Analysis/tutorial_data/046metadataUniform32125.csv",sep=','))

Cellannotataions =data.frame((ReadCellannotataion[,1]))

Cell_ANN_Colnames= data.frame(Column1 = Cellannotataions[,1])

Cell_ledien = unique(ReadCellannotataion$leiden)


###################################################################################

###############################Indexing and extracting common Columns between Cell data and Cell annotation##################
xenium.obj@meta.data ['Cellordinate'] <- Cell_Coor
xenium.obj@meta.data ['leiden'] <- "notassigned"
xenium.obj@meta.data ['Cell_ID'] <- data_table_Colnames


dim(xenium.obj@meta.data)
dim(ReadCellannotataion)

matches_indices <- match( Cell_ANN_Colnames$Column1,data_table_Colnames$Column2)


# Print the indices of matched strings in the longer column
Index4= print(matches_indices)

xenium.obj@meta.data


# Adding Cell type annotation from Scanpy Output to Metadatfile for Niche analysis

library(svMisc)


j=1
for (i in 1:length(Index4)) 
{
  
  xenium.obj@meta.data[Index4[i],13] = ReadCellannotataion[j,21]  #Assign the corresponding interger values for leiden in Metadata
  j=j+1                                                           #and Choose the correct column values from the cell annotation
  progress(i)
}



xenium.obj@meta.data


Metadata_Subset <-  subset(xenium.obj, subset = leiden %in% c("notassigned"), invert = TRUE)

Metadata_Subset@meta.data

dim(xenium.obj@meta.data)
dim(Metadata_Subset@meta.data)
##################################################################################


Metadata_Subset <- BuildNicheAssay(object = Metadata_Subset, fov = "fov", group.by = "leiden",
                                   niches.k = 6, neighbors.k = 25)


celltype.plot <- ImageDimPlot(Metadata_Subset, fov = "fov", group.by = "leiden", size = 2, cols = "polychrome",
                              dark.background = F) + ggtitle("Cell type")


niche.plot <- ImageDimPlot(Metadata_Subset, fov = "fov",  group.by = "niches", size = 2, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B","#50ab4b"))


celltype.plot | niche.plot


table(Metadata_Subset$leiden, Metadata_Subset$niches)

#saving Niche analysis Output

write.csv(Metadata_Subset@meta.data, "C:\\Users\\jj374\\Documents\\RCode\\Brestcancer\\XeniumdataNiche_Test_046_06_25_12625.csv",row.names=FALSE)



