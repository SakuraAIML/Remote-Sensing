# Group 2: Run PCA on 9 bands


#1. Load packages ----
library(raster)     # raster data manipulation tools
library(viridis)    # plotting color palettes
library(rasterVis)  # raster data visualization 
library(corrplot)   # correlation plot 
library(RColorBrewer) # color plots
library(ggpubr)     # arranging plots
library(cluster)    # cluster analysis
library(factoextra) # PCA visualisation
library(FactoMineR) # PCA visualisation
library(clValid) # cluster validation technique

#2. Load data ----
wd <- "S2A_20220319_subset.tif"
# raster brick objects - combined layers

b <- brick(wd)
b
dim(b)

names(b)<- c('Blue','Green','Red','RE-1','RE-2','RE-3','NIR','SWIR-1','SWIR-2')

#3. Create individual band layers ----
# 10 bands (quality scene classification,B2, B3,B4,B5,B6,B7,B8,B11,B12)


b2 <- raster(wd, band = 1) # blue
b3 <- raster(wd, band = 2) # green
b4 <- raster(wd, band = 3) # red
b5 <- raster(wd, band = 4) # veg red edge1
b6 <- raster(wd, band = 5) # veg red edge 2
b7 <- raster(wd, band = 6) # veg red edge 3
b8 <- raster(wd, band = 7) # nir
b11 <- raster(wd, band = 8) # swir1
b12 <- raster(wd, band = 9) # swir2



#4. PCA ----

# Using factoextra ----
values <- getValues(b) # return a matrix with rows representing cells,
# and columns representing layers 
head(values)

set.seed(1)
pr.out <- prcomp(values,
                 scale = TRUE)
names(pr.out)
summary(pr.out)

# Get eigenvalues ----
# Eigenvalues measure the amount of variation retained by each PC
# Eigenvalues are larger for first PCs and small for the subsequent PCs
# That is, the first PCs corresponds to the directions with the max amount of 
# variation in the data set

eig.val <- get_eigenvalue(pr.out)
eig.val # the first 2 PCs have the highest eig values, which cumulatively explain 93.5% of variation

#PC loadings ----
pr.out$rotation

#PC score ----
head(pr.out$x)
dim(pr.out$x) # 2193357       9

# sd of each PC ----
pr.out$sdev

# scree plot ----

fviz_screeplot(pr.out, addlabels = TRUE, barfill = c("#225ea8", "#1d91c0","#41b6c4",
                                                     "#7fcdbb", "#c7e9b4", "#EDF8B1",
                                                     '#556B2F', '#FF7256', '#CAFF70'
                                                     ), 
               barcolor = c("#0c2c84", "#1d91c0","#41b6c4",
                            "#7fcdbb", "#c7e9b4", "#edf8b1",
                            '#556B2F', '#FF7256', "#CAFF70"
                            )) + 
  theme(plot.title = element_text(face = "bold", h = 0.5)) + 
  labs(x = "Principal Components", y = "Percentage of Explained Variances")
ggsave('./Clustering-Group5/PCA/scree.png', scale = 1.5, dpi = 300)

# Graph of variables ----

var <- get_pca_var(pr.out)

head(var$coord, 2) # coordinates
head(var$cos2) # quality of representation
head(var$contrib) # contribution of each var to PCs

fviz_pca_var(pr.out, col.var = "black") # variables correlation
fviz_cos2(pr.out, choice = "var", axes = 1:2) # quality of representation

# color by cos2 values: quality on the factor map
fviz_pca_var(pr.out, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE) # Avoid text overlapping)
ggsave('./Clustering-Group5/PCA/vars_plt.png', scale = 1.5, dpi = 300)

# Contributions of variables to PC1
g1 <- fviz_contrib(pr.out, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
g2 <- fviz_contrib(pr.out, choice = "var", axes = 2, top = 10)
#Total contribution to PC1 and PC2
g3 <- fviz_contrib(pr.out, choice = "var", axes = 1:2, top = 10)

library(cowplot)

ggdraw() +
  draw_plot(g1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(g2, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(g3, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))

ggsave('./Clustering-Group5/PCA/vars_contrib.png', scale = 1.5, dpi = 300)

# Plot most contributing variables
fviz_pca_var(pr.out, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE
) +  scale_fill_discrete(name = "Contribution")

# Color variables by groups
fviz_pca_var(pr.out, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster",
             repel = TRUE) 


#7. Clustering analysis -----
# (i) We will take PC1 and PC2 as input data for clustering
# (ii) Determine the optimal number of clusters using Elbow method, Silhouette
# (iii)Perform clustering on PC1 and PC2
memory.limit(20102*1024)

head(pr.out$x[,1:2])

b_transform <- as.data.frame(pr.out$x[, 1:2]) # create a new data frame which stores the values of PC1 and PC2 to be used as inputs for clustering
head(b_transform)
# create a randomised selection of rows to be used for checking optimal no of clusters
# set the seed, and choose 10,000 samples randomly from the transformed PCA
# use this new data as an input for the checking the optimal number of clusters
# repeat the step 5 times to compare, each time with different seed as we want to create
# a different set of data 

# (i)
set.seed(1)
dat1 <- b_transform[sample(nrow(b_transform), 10000), ]


# draw plots
set.seed(1)
g1 <- fviz_nbclust(dat1, kmeans, method = "wss") + ggtitle("Elbow Method") +
  theme(plot.title = element_text(h = 0.5)) # 3 or 4

set.seed(1)
g2 <- fviz_nbclust(dat1, kmeans, method = 'silhouette')+  ggtitle("Silhouette Method") +
  theme(plot.title = element_text(h = 0.5))# 4

plt <- ggarrange(g1, g2,
                 labels = c("A", "B"),
                 ncol = 2, nrow = 1)
ggsave("./Clustering-Group5/PCA/elbow_sil.png", scale = 1.5, dpi = 300)

# (ii)
set.seed(2)
dat2 <- b_transform[sample(nrow(b_transform), 10000), ]

dim(dat2)

set.seed(2)
fviz_nbclust(dat2, kmeans, method = "wss") # 3 or 4

set.seed(2)
fviz_nbclust(dat2, kmeans, method = 'silhouette') # 4

# (iii)
set.seed(3)
dat3 <- b_transform[sample(nrow(b_transform), 10000), ]

dim(dat3)

set.seed(3)
fviz_nbclust(dat3, kmeans, method = "wss") # 3 or 4

set.seed(3)
fviz_nbclust(dat3, kmeans, method = 'silhouette') # 4

# (iv)
set.seed(4)
dat4 <- b_transform[sample(nrow(b_transform), 10000), ]

dim(dat4)

set.seed(4)
fviz_nbclust(dat4, kmeans, method = "wss") # 3 or 4

set.seed(4)
fviz_nbclust(dat4, kmeans, method = 'silhouette') # 4

# (iiv)
set.seed(5)
dat5 <- b_transform[sample(nrow(b_transform), 10000), ]

dim(dat5)

set.seed(5)
fviz_nbclust(dat5, kmeans, method = "wss") # 3 or 4

set.seed(5)
fviz_nbclust(dat5, kmeans, method = 'silhouette') # 4

# So, the most common no of clusters for kmeans is 3 or 4, will choose 4

# kmeans ----

# Run Kmeans k = 4
set.seed(123)
kmeans_pca = kmeans(dat1, centers = 4, nstart = 50)

# kmeans on transformed PCs
set.seed(1234)
kmeans_pca2 = kmeans(b_transform, centers = 4, nstart = 50)

#clusters silhoutte plot for k=4
sil1 <- silhouette(kmeans_pca$cluster, dist(dat1))

sil1_plt <- fviz_silhouette(sil1, palette = "jco", ggtheme = theme_classic()) + theme(plot.title = element_text(h = 0.5))

# run kmeans k = 3
set.seed(123)
kmeans_pca2 = kmeans(dat1, centers = 3, nstart = 50)
#clusters silhoutte plot for k=3
sil2 <- silhouette(kmeans_pca2$cluster, dist(dat1))

head(sil[, 1:3], 5)

g4 <- fviz_silhouette(sil2, palette = "jco", ggtheme = theme_classic()) + theme(plot.title = element_text(h = 0.5))

# combining graphs
plt <- ggarrange(g1, g2, 
                 g4, g3
)

ggsave('./Clustering-Group5/PCA/optimal_clusters.png', scale = 1.5, dpi = 300)

# clustering plots on PC1 and PC2 for k=4
fviz_cluster(kmeans_pca, data = dat1,
             geom = "point",
             palette = "jco",
             ellipse.type = "norm",
             addEllipses = TRUE, legend.title = "Cluster"

             ) + ggtitle("Visualisation of clustered data using K-means clustering algorithm")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave('./Clustering-Group5/PCA/kmeans_pca.png', scale = 1.5, dpi = 300)

# Plot the image
# function to restrict prediction to the
# first two principal components
pca_predict2 <- function(model, data, ...) {
  predict(model, data, ...)[,1:2]
}

pci <- predict(b, pr.out, fun=pca_predict2)
plot(pci)

head(pr.out)

# clara -----
library(cluster)
?cluster::clara()

# (i)
set.seed(1)
dat1 <- b_transform[sample(nrow(b_transform), 10000), ]

dim(dat1)

set.seed(1)
fviz_nbclust(dat1, clara, method = "wss") # 3 or 4

set.seed(1)
fviz_nbclust(dat1, clara, method = 'silhouette') # 3

# (ii)
set.seed(2)
dat2 <- b_transform[sample(nrow(b_transform), 10000), ]

dim(dat2)

set.seed(2)
fviz_nbclust(dat2, clara, method = "wss") # 3 or 4

set.seed(2)
fviz_nbclust(dat2, clara, method = 'silhouette') # 3

# (iii)
set.seed(3)
dat3 <- b_transform[sample(nrow(b_transform), 10000), ]

dim(dat3)

set.seed(3)
fviz_nbclust(dat3, clara, method = "wss") # 3 

set.seed(3)
fviz_nbclust(dat3, clara, method = 'silhouette') # 4

# (iv)
set.seed(4)
dat4 <- b_transform[sample(nrow(b_transform), 10000), ]

dim(dat4)

set.seed(4)
fviz_nbclust(dat4, clara, method = "wss") # 3 

set.seed(4)
fviz_nbclust(dat4, clara, method = 'silhouette') # 3

# (iiv)
set.seed(5)
dat5 <- b_transform[sample(nrow(b_transform), 10000), ]

dim(dat5)

set.seed(5)
fviz_nbclust(dat5, clara, method = "wss") #  4

set.seed(5)
fviz_nbclust(dat5, clara, method = 'silhouette') # 3

# clara algorithm also chooses 3 or 4 clusters as optimum
# plot both clara and kmeans algorithms on PCs to compare the nb clusters analysis

g1 <- fviz_nbclust(dat1, kmeans, method = "wss") + ggtitle("Elbow method") +
  theme(plot.title = element_text(hjust = 0.5))


g2 <- fviz_nbclust(dat1, kmeans, method = 'silhouette') + ggtitle("Silhouette method") +
  theme(plot.title = element_text(hjust = 0.5))

kmeans_plt <- fviz_cluster(kmeans_pca, data = dat1,
             geom = "point",
             palette = "jco",
             ellipse.type = "norm",
             addEllipses = TRUE, legend.title = "Cluster")+
  theme(plot.title = element_text(hjust = 0.5))+ggtitle("K-means clustering plot") + theme(legend.position = "bottom")

ggdraw() +
  draw_plot(g1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(g2, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(g3, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))


ggsave('./Clustering-Group5/PCA/kmeans_cluster analysis.png', scale = 1.5, dpi = 300)

# clara -----

g1 <- fviz_nbclust(dat5, clara, method = "wss") + ggtitle("Elbow method") +
  theme(plot.title = element_text(hjust = 0.5))


g2 <- fviz_nbclust(dat5, clara, method = 'silhouette') + ggtitle("Silhouette method") +
  theme(plot.title = element_text(hjust = 0.5))

#clusters silhoutte plot for k=4
set.seed(1234)
clara_pca <- clara(dat5, k = 4, metric = "manhattan", samples = 500)

# using full transformed data
set.seed(12345)
clara_pca2 <- clara(b_transform, k = 4, metric = "manhattan", samples = 500)
sil1_cl <- silhouette(clara_pca$cluster, dist(dat5))

sil1_cl_plt <- fviz_silhouette(sil1_cl, palette = "jco", ggtheme = theme_classic()) + theme(plot.title = element_text(h = 0.5))

#clusters silhoutte plot for k=3
set.seed(1234)
clara_pca2 <- clara(dat5, k = 3, metric = "manhattan", samples = 500)
sil_cl2 <- silhouette(clara_pca2$cluster, dist(dat5))


sil1_cl_plt2 <- fviz_silhouette(sil_cl2, palette = "jco", ggtheme = theme_classic()) + theme(plot.title = element_text(h = 0.5))

plt <- ggarrange(g1, g2,
                 sil1_cl_plt2, sil1_cl_plt,
                 nrow = 2, ncol = 2)
ggsave('./Clustering-Group5/PCA/clara_cluster analysis.png', scale = 1.5, dpi = 300)



dd <- cbind(dat5, cluster = clara_pca$clustering)
head(dd, n =4)

clara_plt <- fviz_cluster(clara_pca, data = dat5,
                   geom = "point",
                   palette = "jco",
                   ellipse.type = "norm",
                   addEllipses = TRUE, legend.title = "Cluster")+ ggtitle("CLARA clustering plot")+
  theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "bottom")

plt <- ggarrange(kmeans_plt, clara_plt,
                 nrow=1, ncol = 2)
#ggdraw() +
  #draw_plot(g1, x = 0, y = .5, width = .5, height = .5) +
  #draw_plot(g2, x = .5, y = .5, width = .5, height = .5) +
  #draw_plot(g3, x = 0, y = 0, width = 1, height = 0.5) +
  #draw_plot_label(label = c("A", "B", "C"), size = 15,
                  #x = c(0, 0.5, 0), y = c(1, 1, 0.5))


ggsave('./Clustering-Group5/PCA/clara_cluster analysis.png', scale = 1.5, dpi = 300)


# clustering using pc -----
raster <- b # copy raster brick

raster_df <- as.data.frame(values(b)) # save values in dataframe
head(raster_df)


pca <- prcomp(raster_df, scale. = TRUE)

# save the principal components
raster.pc <- pca$x
raster$PC1 <- raster.pc[, 1]
raster$PC2 <- raster.pc[, 2]
raster$PC3 <- raster.pc[, 3]

# plot PC1:3
ggR(raster,10:12, geom_raster = TRUE) +
  scale_fill_gradientn( colours = rainbow(250), guide = FALSE)+
  ggtitle("PCA transformation")+
  xlab("Longitude") +
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))


# Train k-means on the first and second principle component
# 4 classes
set.seed(123456)
raster$class.pc2 <- kmeans(raster.pc[,1:2], 4)$cluster
# Plot the result
plot(
  raster$class.pc2, 
  col = cm.colors(4), 
  axes = FALSE,
  legend = FALSE
)

plt1 <- ggR(raster$class.pc2, geom_raster = TRUE) +
  scale_fill_gradientn( colours = c("yellow", "blue", "red", "#7CFC00"), guide = FALSE)+
  ggtitle("K-means clustering (k = 4)")+
  xlab("Longitude") +
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

# k means = 3
set.seed(123456)
raster$class.k3 <- kmeans(raster.pc[,1:2], 3)$cluster
plt2 <- ggR(raster$class.k3, geom_raster = TRUE) +
  scale_fill_gradientn( colours = c("yellow", "blue", "green"), guide = FALSE)+
  ggtitle("K-means clustering (k = 3)")+
  xlab("Longitude") +
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

## clara
# k = 4
set.seed(123456)
raster$class.clara <- clara(raster.pc[,1:2], 4, metric = "manhattan")$cluster

plt3 <- ggR(raster$class.clara, geom_raster = TRUE) +
  scale_fill_gradientn( colours = c("#7CFC00", "yellow", "blue", "red"), guide = FALSE)+
  ggtitle("CLARA clustering (k = 4)")+
  xlab("Longitude") +
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

# k = 3

set.seed(123456)
raster$class.clara.k3 <- clara(raster.pc[,1:2], 3, metric = "manhattan")$cluster

plt4 <- ggR(raster$class.clara.k3, geom_raster = TRUE) +
  scale_fill_gradientn( colours = c( "green", "blue", "yellow"), guide = FALSE)+
  ggtitle("CLARA clustering (k = 3)")+ 
  xlab("Longitude") +
  ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))



gridExtra::grid.arrange(plt2, plt1, plt4, plt3, ncol = 2, nrow = 2)
ggsave('./Clustering-Group5/PCA/final_solution.png', scale = 1.5, dpi = 300)

