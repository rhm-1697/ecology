setwd("P:/GEOG70922/Assessment2")


#Install the required libraries

install.packages(c("raster","rgdal","rgeos","gdistance","igraph","fitdistrplus", "RColorBrewer"))
install.packages(c("rasterVis"))

#Load the installed libraries

library(raster) #for raster covariate data

library(rgdal) #for reading different types of GIS files

library(rgeos) #for centroids of polygons

library(gdistance) #for least-cost paths/circuit theory

library(igraph) #for patch-based graphs

library(fitdistrplus) #for fitting kernels

library(RColorBrewer) #for discrete color palette for the plot

library(maptools)

library(rasterVis)


#Read the raster file for Corine Land Cover 2018 data set

corine <- raster("CLC_UK.tif")


#Read in the national park shape files

Cairngorms <- readOGR("SG_CairngormsNationalPark_2010.shp")

LochLomond <- readOGR("SG_LochLomondTrossachsNationalPark_2002.shp")


#Transform the park shape files to the same CRS as the land cover raster layer

Cairngorms <- spTransform(Cairngorms,crs(corine))

LochLomond <- spTransform(LochLomond,crs(corine))


#Get the coordinates of the Cairngorms park

CCoords <- coordinates(Cairngorms)



#Extend the study area by using a buffer

x.min <- min(CCoords[,1]) - 120000

x.max <- max(CCoords[,1]) + 80000

y.min <- min(CCoords[,2]) - 120000

y.max <- max(CCoords[,2]) + 50000



#Get the extent of the buffered area

extent.new <- extent(x.min, x.max, y.min, y.max)



#Crop land cover to the new extent

corine <- crop(corine, extent.new)


#Get the centroids of both parks

Cairngorms_centroid <- gCentroid(Cairngorms, byid = T)

LochLomond_centroid <- gCentroid(LochLomond, byid = T)


#Plot the study area

plot(corine)

plot(Cairngorms, add = T)

plot(LochLomond, add = T)

plot(Cairngorms_centroid, add = T)

plot(LochLomond_centroid, add = T)


legendclc<-read.csv("clc_legend.csv")#read in spreadsheet with legend key



legendclc<-legendclc[,c(1,4)] #select the relevant columns (ID and class names) 

head(legendclc)

landCoverclc<-levels(corine)[[1]]#create object of landcover classes

legendScot <- legendclc[legendclc$CLC_CODE %in% landCoverclc$ID,]#use %in% command to match classes in the raster to those in the legend


par(mar=c(10, 10, 10, 10)) #set margins

##Use colourRampPalette to create custom palette
coul <- colorRampPalette(brewer.pal(8, "Accent"))(33) #select 8 colours from the "Accent" palette and create 33 combinations

plot(corine,legend = FALSE,col = coul)#plot again with the new palette

legend(406000,895000,paste(legendScot[,2]), fill = coul,
       cex = 0.4,bty="n")

#Finally add an entry to symbolize the betweenness layer.
legend(x=65000,y=780000,legend=c("Patch betweenness"),pch=21,col="black",bg="white",cex=0.6,bty="n")


################################################################################################

#Parameterizing a cost (resistance) layer between the Cairngorms and Loch Lomond National Parks

################################################################################################



#Raster categories as factors

corine <- as.factor(corine)



#Inspect levels

levels(corine)



#Resistance values for each land cover type (Stringer et al., 2018)

cost <- c(rep(1000, 3), 200, rep(1000, 7), rep(60, 2), 30, rep(60, 2), rep(1,3), 30, rep(5,2), rep(40,2), rep(5,2), rep(100,4), rep(1000,5))



#Combine the land cover types and the resistance values using the cbind function

RCmat <- cbind(levels(corine)[[1]], cost)


cost_layer <- reclassify(corine, RCmat)



#set breaks based on values of the raster layer

cuts=c(1,5,30,40,60,100,200,1000) 

legend_names <- c("1 - Forest", "5 - Bushes", "30 - Pasture and Grassland", "40 - Rocks", "60 - Crops", "100 - Wetlands", "200 - National Roads", "1000 - Urban Areas and Waterbodies")

#set the color palette for the defined breaks

par(mar=c(5, 5, 5, 5)) #set margins


pal <- colorRampPalette(c("#254117","#347C17","#12AD2B","#77DD77", "#FBB917", "#FFA600", "#EB5406","#FF0000" ))



plot(cost_layer, breaks=cuts, col = pal(7), legend = F, box = T) #plot with defined breaks


legend(400000, 890000, 
       title = "Resistance values",  #title for relevant graph
       inset=c(-0.3,0),
       legend = legend_names, #variable that holds legend names 
       fill = pal(8),
       bty = "n",
       cex = 0.5)


############################################################################

#Plotting a potential wildlife corridor based on known habitat preference

############################################################################


#Create a conductance-based transition layer, with mean of neighboring 8 cells as the transition function

land_cond <- transition(1/cost_layer, transitionFunction=mean, 8)


#Create cumulative cost layers originating from each centroid

C1.cost <- accCost(land_cond, Cairngorms_centroid)

C2.cost <- accCost(land_cond, LochLomond_centroid)



#Get least-cost corridor - use the overlay() function of the raster package to add both rasters together

leastcost_corridor <- overlay(C1.cost, C2.cost, fun=function(x, y){return(x + y)})



#Get the lower quantile by clipping the corridor based on the lowest 7% cost of the overlaid centroids

quantile7 <- quantile(leastcost_corridor, probs=0.07, na.rm=TRUE)



#Make a new truncated layer so only values less than the quantile value are given a value

leastcost_corridor7 <- leastcost_corridor

values(leastcost_corridor7) <- NA

leastcost_corridor7[leastcost_corridor < quantile7] <- 1 # Truncate to identify corridor



# Plot the least cost corridor

plot(corine, legend = F)

plot(leastcost_corridor7, legend=F,axes=F, col = "red", add = T)

points(Cairngorms_centroid, col="grey30")

points(LochLomond_centroid, col="grey30")


#########################################################################################################

#Highlight the relative importance of habitat patches as stepping stones between the two protected areas.

#########################################################################################################


#Read the conifer patch data

conifer <- readOGR("Conifer_Patch.shp")



#Extract the centroids of the conifer patches

conifer_centroids <- gCentroid(conifer, byid = T)



#Create a variable which will store the coordinates of each patch centroid

coords_conifer_centroids <- cbind(conifer_centroids$x, conifer_centroids$y)



#Create a matrix of distance values between the centroids

distmat <- pointDistance(coords_conifer_centroids, lonlat=F)



# Change the distance unit to km for easier interpretability

distmat <- distmat/1000



#Set the mean distance (8 km) as based on Stringer et al. (2018)

mean.dist <- 8



#Create an adjacency matrix which uses a negative exponential function

A.prob <- matrix(0, nrow=nrow(distmat), ncol=ncol(distmat))

alpha <- 1/mean.dist

A.prob <- exp(-alpha*distmat) #Implement negative exponential function

diag(A.prob) <- 0 #Ensures circular connectivity is not considered



#Create an igraph object

graph.Aprob <- graph.adjacency(A.prob, mode="undirected", weighted=T)


#Create a variable that will store betweeness values as an igraph object

Aprob.between <- betweenness(graph.Aprob, weights=1/E(graph.Aprob)$weight)



#Plot betweeness centrality

plot.new()

plot(corine, axes = T, box = F, legend = F)

plot(Cairngorms, add = T)

plot(LochLomond, add = T)

#plot(conifer, add = T, col = "grey")

plot(graph.Aprob, layout=coords_conifer_centroids,vertex.size=Aprob.between*200, vertex.color = "blue", edge.color = "transparent", main= "BetweennessProb", vertex.label=NA, frame=T, rescale = F, add = T)


##########################################################


study_area <- as(raster::extent(corine),"SpatialPolygons")

study_area <- gArea(study_area)/1000000 # Total area of study in km^2



# Area of each conifer patch in km^2

area <- conifer$Shape_Area/1000000



# Calculate landscape-scale connectivity

pstar.mat <- shortest.paths(graph.Aprob, weights= -log(E(graph.Aprob)$weight)) # Calculate all shortest paths between nodes

pstar.mat <- exp(-pstar.mat) # Back-transform to probabilities of connectedness

PCnum <- outer(area, area)*pstar.mat # Get product of all patch areas ij and multiply by probabilities above

PC <- sum(PCnum)/study_area^2 # Divide by total area of the study squared to get the PC metric



# Function for calculating dPC - contribution of each patch to PC.

# Remove each patch using a for loop and repeat the above function

# Adding the difference between PC above and each new calculation (minus patch "i") to an empty vector (dPC)

prob.connectivity <- function(prob.matrix, area,landarea){
  
  
  # dPC
  
  N <- nrow(prob.matrix) # Each row is a patch
  
  dPC <- rep(NA, N) # Empty vector to store results for each patch
  
  
  
  
  for (i in 1:N) {
    
    prob.matrix.i <- prob.matrix[-i,-i]
    
    area.i <-area[-i]
    
    pc.graph.i <- graph.adjacency(prob.matrix.i, mode="undirected", weighted=TRUE)
    
    pstar.mat.i <- shortest.paths(pc.graph.i, weights= -log(E(pc.graph.i)$weight))
    
    pstar.mat.i <- exp(-pstar.mat.i)
    
    PCmat.i <- outer(area.i, area.i)*pstar.mat.i
    
    PC.i <- sum(PCmat.i)/landarea^2
    
    dPC[i] <- (PC-PC.i)/PC*100
    
  }
  
  
  return(dPC)
  
}



# Store dPC value in a variable as an igraph object for plotting

Aprob.PC <- prob.connectivity(prob.matrix=A.prob, area=area, landarea=study_area)

#Get the coordinates of the Cairngorms park

CCoords <- coordinates(Cairngorms)



#Extend the study area by using a buffer

x.min <- min(CCoords[,1]) - 150000

x.max <- max(CCoords[,1]) + 100000

y.min <- min(CCoords[,2]) - 150000

y.max <- max(CCoords[,2]) + 100000



#Get the extent of the buffered area

extent.new <- extent(x.min, x.max, y.min, y.max)



#Crop land cover to the new extent

corine1 <- crop(corine, extent.new)


# Plot dPC probability
plot.new()
plot(corine1, axes = F, box = F)

plot(Cairngorms, add = T)

plot(LochLomond, add = T)

#plot(conifer, add = T, col = "red")

plot(graph.Aprob,layout=coords_conifer_centroids,vertex.size=Aprob.PC*1.5, main= "dPC", vertex.color = "blue", vertex.label.cex = Aprob.PC*0.1 ,vertex.label.color="white", edge.color = "red", rescale = T, add = F)


plot.new()

