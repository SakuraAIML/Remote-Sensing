# Group 1: all 9 bands


#1. Load packages ----
library(raster)     # raster data manipulation tools
library(ggplot2)    # data visualization tools
library(RStoolbox) # raster visualisation tools
library(rasterVis)  # raster data visualization 
library(corrplot)   # correlation plot 
library(ggpubr)     # arranging plots
library(cluster)    # cluster analysis


#2. Load data ----
wd <- "S2A_20220319_subset.tif"

# raster brick objects - combined layers

b <- brick(wd)
b


#3. Create individual band layers ----
# 9 bands (quality scene classification,B2, B3,B4,B5,B6,B7,B8,B11,B12)


b2 <- raster(wd, band = 1) # blue
b3 <- raster(wd, band = 2) # green
b4 <- raster(wd, band = 3) # red
b5 <- raster(wd, band = 4) # veg red edge1
b6 <- raster(wd, band = 5) # veg red edge 2
b7 <- raster(wd, band = 6) # veg red edge 3
b8 <- raster(wd, band = 7) # nir
b11 <- raster(wd, band = 8) # swir1
b12 <- raster(wd, band = 9) # swir2

plot(b12)
#4. Create RasterStack ----

rst <- stack(b2, b3, b4, b5, b6, b7, b8, b11, b12)


# rename each layer within the RasterStack object
names(rst)<- c('Blue','Green','Red','RE-1','RE-2','RE-3','NIR','SWIR-1','SWIR-2')
rst


#5. RGB True & FCC colour composite -----


sentinelRGB <- stack(list(b4, b3, b2))
sentinelFCC <- stack(list(b8, b4, b3))

# create a folder to store plots
if(!dir.exists("./Plots")) dir.create("./Plots")

png('./Plots/true_false_composite.png', width = 2000, height = 500, units = "cm", res = 300)
nf <- layout(matrix(c(1,0,2), 1, 3, byrow = TRUE), width = c(1,0.01,1), respect = TRUE)
plotRGB(sentinelRGB, stretch = "lin", main = "Sentinel True Color Composite", axes = TRUE, margins = FALSE)
plotRGB(sentinelFCC, stretch = "lin", main = "Sentinel False Color Composite", axes = TRUE, margins = FALSE)
dev.off()


#6. Individual bands ----
options(repr.plot.width = 10, repr.plot.height =12)
gplot(rst) +
  geom_raster(aes(x = x, y = y, fill = value))+
  scale_fill_viridis_c() +
  facet_wrap(~variable) +
  coord_quickmap()+
  ggtitle("Sentinel-2 bands") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_classic() +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("./Plots/allbands_tutorial5.png", scale = 1.5, dpi = 300) # to save plot

#6. Clustering ----
# (i) cluster image data with K-MEANS and CLARA for a number of clusters between 2:10
# (ii) assess each clustering solution performance through the average silhouette index

# Extract all values from the raster into a data frame
rstDF <- values(rst)
head(rstDF)
sum(is.na(rstDF))

# check percentage of missing values in each band
p <- function(x) {sum(is.na(x))/length(x)*100}
apply(rstDF, 2, p)

# Check NA's in the data
idx <- complete.cases(rstDF)

# Initiate the raster datasets that will hold all clustering solutions 
# from 2 groups/clusters up to 10
rstKM <- raster(rst[[1]])
rstCLARA <- raster(rst[[1]])


for(nClust in 2:10){
  
  cat("-> Clustering data for nClust =",nClust,"......")
  
  set.seed(1)
  # Perform K-means clustering
  km <- kmeans(rstDF[idx,], centers = nClust, iter.max = 500, nstart = 3)
  
  # Perform CLARA's clustering (using manhattan distance)
  cla <- clara(rstDF[idx, ], k = nClust, metric = "manhattan")
  
  # Create a temporary integer vector for holding cluster numbers
  kmClust <- vector(mode = "integer", length = ncell(rst))
  claClust <- vector(mode = "integer", length = ncell(rst))
  
  # Generate the temporary clustering vector for K-means (keeps track of NA's)
  kmClust[!idx] <- NA
  kmClust[idx] <- km$cluster
  
  # Generate the temporary clustering vector for CLARA (keeps track of NA's too ;-)
  claClust[!idx] <- NA
  claClust[idx] <- cla$clustering
  
  # Create a temporary raster for holding the new clustering solution
  # K-means
  tmpRstKM <- raster(rst[[1]])
  # CLARA
  tmpRstCLARA <- raster(rst[[1]])
  
  # Set raster values with the cluster vector
  # K-means
  values(tmpRstKM) <- kmClust
  # CLARA
  values(tmpRstCLARA) <- claClust
  
  # Stack the temporary rasters onto the final ones
  if(nClust==2){
    rstKM    <- tmpRstKM
    rstCLARA <- tmpRstCLARA
  }else{
    rstKM    <- stack(rstKM, tmpRstKM)
    rstCLARA <- stack(rstCLARA, tmpRstCLARA)
  }
  
  cat(" done!\n\n")
}

if(!dir.exists("./Clustering-Group1")) dir.create("./Clustering-Group1")
# Write the clustering solutions for each algorithm
writeRaster(rstKM,"./Clustering-Group1/KMeans.tif", overwrite=TRUE)
writeRaster(rstCLARA,"./Clustering-Group1/CLARA.tif", overwrite=TRUE)

#4. Silhoutte analysis -----

rstKM <- brick("./Clustering-Group1/KMeans.tif")
rstCLARA <- brick("./Clustering-Group1/CLARA.tif")

# Start a data frame that will store all silhouette values
# for k-means and CLARA   
clustPerfSI <- data.frame(nClust = 2:10, SI_KM = NA, SI_CLARA = NA)


for(i in 1:nlayers(rstKM)){ # Iterate through each layer
  
  cat("-> Evaluating clustering performance for nClust =",(2:10)[i],"......")
  set.seed(1)
  # Extract random cell samples stratified by cluster
  cellIdx_RstKM <- sampleStratified(rstKM[[i]], size = 10000)
  cellIdx_rstCLARA <- sampleStratified(rstCLARA[[i]], size = 10000)
  
  # Get cell values from the Stratified Random Sample from the raster 
  # data frame object (rstDF)
  rstDFStRS_KM <- rstDF[cellIdx_RstKM[,1], ]
  rstDFStRS_CLARA <- rstDF[cellIdx_rstCLARA[,1], ]
  
  # Make sure all columns are numeric (intCriteria function is picky on this)
  rstDFStRS_KM[] <- sapply(rstDFStRS_KM, as.numeric)
  rstDFStRS_CLARA[] <- sapply(rstDFStRS_CLARA, as.numeric)
  
  # Compute the sample-based Silhouette index for: 
  #    
  # K-means
  clCritKM <- intCriteria(traj = rstDFStRS_KM, 
                          part = as.integer(cellIdx_RstKM[,2]), 
                          crit = "Silhouette")
  # and CLARA
  clCritCLARA <- intCriteria(traj = rstDFStRS_CLARA, 
                             part = as.integer(cellIdx_rstCLARA[,2]), 
                             crit = "Silhouette")
  
  # Write the silhouette index value to clustPerfSI data frame holding 
  # all results
  clustPerfSI[i, "SI_KM"]    <- clCritKM[[1]][1]
  clustPerfSI[i, "SI_CLARA"] <- clCritCLARA[[1]][1]
  
  cat(" done!\n\n")
  
}

write.csv(clustPerfSI, file = "./Clustering-Group1/clustPerfSI.csv", row.names = FALSE)

knitr::kable(clustPerfSI, digits = 3, align = "c", 
             col.names = c("#clusters","Avg. Silhouette (k-means)","Avg. Silhouette (CLARA)"))

png('./Clustering-Group1/SI.png', width = 5, height = 5, units = "in", res = 300)
par(mfrow = c(1, 1))
plot(clustPerfSI[,1], clustPerfSI[,2], 
     xlim = c(2,10), ylim = range(clustPerfSI[,2:3]), type = "n", 
     ylab="Average Silhouette Index", xlab="Number of clusters",
     main="Silhouette Index by Number of Clusters")

# Plot Avg Silhouette values across # of clusters for K-means
lines(clustPerfSI[,1], clustPerfSI[,2], col="red")
# Plot Avg Silhouette values across # of clusters for CLARA
lines(clustPerfSI[,1], clustPerfSI[,3], col="blue")

# Grid lines
abline(v = 1:9, lty=2, col="light grey")
abline(h = seq(0.30,0.60,0.05), lty=2, col="light grey")


legend("topright", legend=c("K-means","CLARA"), col=c("red","blue"), lty=1, lwd=1)
dev.off()

# plot the solutions

names(rstKM) <- c("k_2", "k_3", "k_4", "k_5", "k_6",
                  "k_7", "k_8", "k_9", "k_10")
rstKM

names(rstCLARA) <- c("k_2", "k_3", "k_4", "k_5", "k_6",
                     "k_7", "k_8", "k_9", "k_10")
rstCLARA


# clusters 2 to 7
ggR(rstKM, 1:6, geom_raster = TRUE, stretch = "lin") +
  ggtitle("K-means clustering ")+
  scale_fill_gradientn("Cluster", colours = c("blue", "yellow", "green", "red"), guide = FALSE)+
  xlab("Longitude") +
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size=5), 
        axis.text.y = element_text(angle=90),
        axis.title=element_blank())


ggR(rstCLARA, 1:6, geom_raster = TRUE, stretch = "lin") +
  ggtitle("CLARA clustering ")+
  scale_fill_gradientn("Cluster", colours = topo.colors(4), guide = FALSE)+
  xlab("Longitude") +
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size=5), 
        axis.text.y = element_text(angle=90),
        axis.title=element_blank())


