library (clValid)
b_transform <- as.data.frame(pr.out$x[, 1:2]) # create a new data frame which stores the values of PC1 and PC2 to be used as inputs for clustering
head(b_transform)


set.seed(1)
dat1 <- b_transform[sample(nrow(b_transform), 10000), ]

# Compute clValid
clmethods <- c("kmeans","clara")
intern <- clValid(dat1, nClust = 2:7,
                  clMethods = clmethods, validation = "internal",
                  maxitems = 600, metric = "euclidean")
# Summary
summary(intern)

layout(matrix(c(1,1,3,3), 2, 2, byrow = TRUE))
plot(intern, cex = 0.6, legend = TRUE, main = NULL,  bty = "n")


# Stability measures
clmethods <- c("kmeans","clara")
stab <- clValid(getValues(b), nClust = 2:6, clMethods = clmethods,
                validation = "stability", maxitems = 600)

# Display only optimal Scores
optimalScores(stab)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow = c(2, 2))
plot(stab)

par(mfrow = c(2, 2))
plot(intern)
