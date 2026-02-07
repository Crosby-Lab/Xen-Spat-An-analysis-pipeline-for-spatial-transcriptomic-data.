###############################################################################
### Estimating Hotspot detection using Geospatial Analysis Pipeline
### Author : Jesuchristopher Joseph (jesu.joseph@duke.edu)
### Date : 
### Environment : R 4.3.3
### Packages  : dplyr: 1.14, plotly: 4.11.0,  spatstat: 3.3-3, maptools: 1.1-8
###             rgeos: 0.6-4, ggplot2 : 3.5.2 , sf: 1.0-21 , sp: 2.2-0
###############################################################################


library(dplyr) 
library(plotly)
library(spatstat)
library(maptools) 
library(rgeos)
library(ggplot2)
library(sf)
library(sp)  
load("data/stations.RData")
load("data/lnd.RData")
class(stations)



personaldata1 <- data.frame(read.csv("046metadataUniform32125.csv"))

max(personaldata1$Leiden)  # Assign unique integer values to the Cell annotations and name the column as Leiden

unique(personaldata1$Cell_Label)

#"Macrophages"    "Tumor cells"   "Tregs"   "CD8 TCells"    "Endothelial Cells"     "Fibroblasts"    


mydata_1 <- personaldata1 %>% filter(Leiden  == 1)
mydata_2 <- personaldata1 %>% filter(Leiden  == 3)
mydata_3 <- personaldata1 %>% filter(Leiden  == 4)
mydata_4 <- personaldata1 %>% filter(Leiden  == 6)
mydata_5 <- personaldata1 %>% filter(Leiden  == 8)
mydata_6 <- personaldata1 %>% filter(Leiden  == 10)

###############################################################################

my_colors <- c(             
  'sandybrown',  # Tumor cells
  'red',         # Endothelial Cells
  "#556B2F",     # Fibroblasts
  'lightgreen',  # Macrophages
  '#800000',     # CD8 TCells
  "tomato"       # Tregs
 
) 


p <- plot_ly(data = personaldata1, x = ~Cell_X_Position, y = ~Cell_Y_Position, type = "scatter", mode = "markers",  marker = list(size = 3), color = ~Leiden, colors = my_colors,marker = list(size = 6)) %>% 
  layout( xaxis = list(title = "Cell_X_Position"), yaxis = list(title = "Cell_Y_Position"))

p <- p %>% layout(yaxis = list(autorange = "reversed"))

p

################################################################################
# Method to find the Max X-Y coordinate value

x1_range <- range(mydata_1$Cell_X_Position) ; y1_range <- range(mydata_1$Cell_Y_Position)
X2_range <- range(mydata_2$Cell_X_Position) ; y2_range <- range(mydata_2$Cell_Y_Position)
X3_range <- range(mydata_3$Cell_X_Position) ; y3_range <- range(mydata_3$Cell_Y_Position)
x4_range <- range(mydata_4$Cell_X_Position) ; y4_range <- range(mydata_4$Cell_Y_Position)
x5_range <- range(mydata_5$Cell_X_Position) ; y5_range <- range(mydata_5$Cell_Y_Position)
x6_range <- range(mydata_6$Cell_X_Position) ; y6_range <- range(mydata_6$Cell_Y_Position)

X_stack <- matrix(c(x1_range, X2_range, X3_range, x4_range, x5_range, x6_range),nrow = 6, ncol =2, byrow = TRUE)
Y_stack <- matrix(c(y1_range, y2_range, y3_range, y4_range, y5_range, y6_range),nrow = 6, ncol =2, byrow = TRUE)
XY_Stack <- cbind(X_stack,Y_stack)

colnames(XY_Stack) <- c('A', 'B', 'C', 'D')

X_low <- round(min(XY_Stack[,1]))
X_high <- round(max(XY_Stack[,2]))
Y_low <- round(min(XY_Stack[,3]))
Y_high <- round(max(XY_Stack[,4]))

win1 <- owin(xrange = c(X_low-100,X_high+100), yrange = c(Y_low-100,Y_high+100))
XY_low_High <- as.data.frame(win1)



x1 <- mydata_1$Cell_X_Position
y1 <- mydata_1$Cell_Y_Position

x2 <- mydata_2$Cell_X_Position
y2 <- mydata_2$Cell_Y_Position

x3 <- mydata_3$Cell_X_Position
y3 <- mydata_3$Cell_Y_Position


x4 <- mydata_4$Cell_X_Position
y4 <- mydata_4$Cell_Y_Position


x5 <- mydata_5$Cell_X_Position
y5 <- mydata_5$Cell_Y_Position


x6 <- mydata_6$Cell_X_Position
y6 <- mydata_6$Cell_Y_Position

# Generating PPP object for individual cell types

X1 <- ppp(x = x1, y = y1, window = win1)
X1

X2 <- ppp(x = x2, y = y2, window = win1)
X2

X3 <- ppp(x = x3, y = y3, window = win1)
X3

X4 <- ppp(x = x4, y = y4, window = win1)
X4

X5 <- ppp(x = x5, y = y5, window = win1)
X5

X6 <- ppp(x = x6, y = y6, window = win1)
X6


y1_range <- range(X1$y)
y2_range <- range(X2$y)
y3_range <- range(X3$y)
y4_range <- range(X4$y)
y5_range <- range(X5$y)
y6_range <- range(X6$y)


png(filename = "046_Tumor_cells_Plot.png", width = 1000, height = 746)


plot(X1, add = FALSE, cex = 1,main = "Tumor Cells")

plot(X1, add = TRUE, cex = 1,main = "Tumor Cells", col = "sandybrown", axes = TRUE)


dev.off()
dev.new() 

plot(X2, ylim = rev(y2_range),add = FALSE, cex = 1,main = "Endothelial Cells")

plot(X3,ylim = rev(y3_range), add = FALSE, cex = 1,main = "Fibroblasts")

plot(X4, ylim = rev(y4_range),add = FALSE, cex = 1,main = "Macrophages")


png(filename = "046_CD8 T cells_Plot.png", width = 1000, height = 746)


plot(X5, add = FALSE, cex = 1,main = "CD8 TCells")


plot(X5, add = TRUE, cex = 1,main = "CD8 TCells", col = "#800000", axes = TRUE)

dev.off()
dev.new() 

plot(X6, ylim = rev(y6_range), add = FALSE, cex = 1,main = "Tregs")


#################################################################################
# Generating Density Plot, Contour plot and Contour subsetting fro Tumor Samples


png(filename = "046_Tumor_cells_Density_Plot.png", width = 1000, height = 746)

Dens1 <- density(X1, adjust = 0.1)  
class(Dens1)
plot(Dens1, main="Tumor Cells")

dev.off()
dev.new() 


png(filename = "046_Tumor_Cells_Contour.png", width = 1000, height = 746)


Dsg1 <- as(Dens1, "SpatialGridDataFrame")  
Dim1 <- as.image.SpatialGridDataFrame(Dsg1)  
Dcl1 <- contourLines(Dim1, nlevels = 9)  
SLDF1 <- ContourLines2SLDF(Dcl1, CRS(proj4string(lnd)))  
plot(SLDF1, col = terrain.colors(8), main="Tumor_Cells")

dev.off()
dev.new() 


png(filename = "046_Tumor_Cells_Contour_Selected.png", width = 1000, height = 746)

Polyclust1 <- gPolygonize(SLDF1[3, ])
gas1 <- gArea(Polyclust1, byid = T)/10000
median(gas1)
Polyclust1 <- SpatialPolygonsDataFrame(Polyclust1, data = data.frame(gas1), match.ID = F)
plot(Polyclust1,main="Tumor_Cells")

dev.off()
dev.new() 

png(filename = "046_Tumor_Cells_Contour_Selected_Subset.png", width = 1000, height = 746)

Polyclustsub1 <- subset(Polyclust1, gas1 > 1)  #Select Contours

plot(Polyclustsub1,main="Tumor_Cells")

dev.off()
dev.new()

area(Polyclust1)
length(Polyclust1)

###############################################################################

# CD8 TCells

dev.new() 

png(filename = "046_CD8_Tcells_Density_Plot.png", width = 1000, height = 746)

Dens3 <- density(X5, adjust = 0.1)  # create density object
class(Dens3)
plot(Dens3,main="CD8 TCells")

dev.off()
dev.new() 

png(filename = "046_CD8_TCells_Contour.png", width = 1000, height = 746)

Dsg3 <- as(Dens3, "SpatialGridDataFrame")  
Dim3 <- as.image.SpatialGridDataFrame(Dsg3) 
Dcl3 <- contourLines(Dim3, nlevels = 9)  
SLDF3 <- ContourLines2SLDF(Dcl3, CRS(proj4string(lnd)))  
plot(SLDF3, col = terrain.colors(8),main="CD8 TCells")

dev.off()
dev.new() 


png(filename = "046_CD8_TCells_Contour_Selected.png", width = 1000, height = 746)

Polyclust3 <- gPolygonize(SLDF3[3, ])
gas3 <- gArea(Polyclust3, byid = T)/10000
median(gas3)
Polyclust3 <- SpatialPolygonsDataFrame(Polyclust3, data = data.frame(gas3), match.ID = F)
plot(Polyclust3,main="CD8 T Cells")

dev.off()
dev.new() 

png(filename = "046_CD8_TCells_Contour_Selected_Subset.png", width = 1000, height = 746)

Polyclustsub3 <- subset(Polyclust3, gas3 > 1)

plot(Polyclustsub3,main="CD8 T Cells")

dev.off()
dev.new()


area(Polyclustsub3)
length(Polyclustsub3)

###############################################################################

# Generating Overlapping contours between Tumor Cells and CD8 T cells


inter<- gIntersection(Polyclustsub1,Polyclustsub3)


png(filename = "046_Tumor_Cells_Contour.png", width = 1000, height = 746)

plot(Polyclustsub1, pch = 19, col = rgb(1, 0, 0, 0.5), cex = 3, xlim = c(X_low-100, X_high+200) , ylim = c(Y_low-100, Y_high+100), lwd =1)

dev.off()
dev.new()

png(filename = "046_Tumor_Cells_CD8 T cellContour.png", width = 1000, height = 746)


plot(Polyclustsub1, pch = 19, col = rgb(1, 0, 0, 0.5), cex = 3, xlim = c(X_low-100, X_high+200) , ylim = c(Y_low-100, Y_high+100), lwd =1)
plot(Polyclustsub3, add = TRUE, pch = 19, col = rgb(0, 1, 0, 0.5), cex =3, lwd =1,main="Luminal Bc Epithelium vs CD8 T Cells")

dev.off()
dev.new()

png(filename = "046_Tumor_Cells_CD8 T cellContour_overlap.png", width = 1000, height = 746)

plot(Polyclustsub1, pch = 19, col = rgb(1, 0, 0, 0.5), cex = 3, xlim = c(X_low-100, X_high+200) , ylim = c(Y_low-100, Y_high+100), lwd =1)

plot(Polyclustsub3, add = TRUE, pch = 19, col = rgb(0, 1, 0, 0.5), cex =3, lwd =1,main="Luminal Bc Epithelium vs CD8 T Cells")

plot(inter, add=TRUE, border="yellow", lwd =2)

dev.off()

area(inter)  # Measure Area
length(inter@polygons[[1]]@Polygons)   # Measure Contour count

