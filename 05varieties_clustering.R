library(mgcv)
library(fossil)
library(cluster)
library(dplyr)
library(factoextra)

source("utils.R")
source("plots.R")

## Prepare data ##
##################

Douglas <- get_Douglas_data()

# only presences - unscaled
DougPres <- Douglas[Douglas$PRES==1,]

# data set complete - scaled
DougScaled <- Douglas
DougScaled[,c(10:19)] <- scale(Douglas[,c(10:19)])
summary(DougScaled)

DougScaled$x <- (DougScaled$x - min(DougScaled$x))/1000
DougScaled$y <- (DougScaled$y - min(DougScaled$y))/1000

# data set presences - scaled
DougScaledPres <- DougScaled[DougScaled$PRES==1,]

# separate by roughly the two varieties based on six large-scale DNA classes.
Ecotypes(DougScaledPres, DougScaledPres$POP)

DougScaledPres$VARS <- "VAR_inland"
DougScaledPres$VARS[DougScaledPres$POP %in% c("DNA_1", "DNA_2")] <- "VAR_coast"
DougScaledPres$VARS[DougScaledPres$POP %in% c("DNA_6")] <- "VAR_mexican"

#==================================#
# Load the varying coefficient model #
#==================================#

load(file="Rdata/fgamSVC.Rdata")
load(file="Rdata/fgamSVC_predTerms.Rdata")
load(file="Rdata/fgamSVC_clusters.Rdata")
summary(fgamSVC)
AIC(fgamSVC)

# predicted occurance probabilities
allPreds <- predict(fgamSVC, newdata=DougScaled, type="response")
# model coefficients at presence observations
predTerms <- predict(fgamSVC, newdata=DougScaled[DougScaled$PRES==1,], type="terms")
save(predTerms, file="Rdata/fgamSVC_predTerms.Rdata")

# Overlap of clusters with varieties.
useCluster <- clusters[[6]]$clustering
cm <- confusion_matrix(as.numeric(as.factor(DougScaledPres$VARS)), useCluster)
xtable(cm, digits=0)
rand.index(as.numeric(as.factor(DougScaledPres$VARS)), useCluster)

#========================#
# Clustering & Diagnosis #
#========================#

coefficient_matrix <- predTerms[,c(2,7)]
coefficient_matrix_inland <- coefficient_matrix[DougScaledPres$VARS == "VAR_inland",]
coefficient_matrix_coast <- coefficient_matrix[DougScaledPres$VARS == "VAR_coast",]
coefficient_matrix_mexican <- coefficient_matrix[DougScaledPres$VARS == "VAR_mexican",]


clusters_inland <- list()
for (i in c(2:15)){
  clusters_inland[[i]] <- pam(coefficient_matrix_inland, k=i, medoids = 'random', variant = 'faster')
}
save(clusters_inland, file="Rdata/fgamSVC_clusters_inland.Rdata")

clusters_coast <- list()
for (i in c(2:15)){
  clusters_coast[[i]] <- pam(coefficient_matrix_coast, k=i, medoids = 'random', variant = 'faster')
}
save(clusters_coast, file="Rdata/fgamSVC_clusters_coast.Rdata")

clusters_mexican <- list()
for (i in c(1:15)){
  clusters_mexican[[i]] <- pam(coefficient_matrix_mexican, k=i, medoids = 'random', variant = 'faster')
}
save(clusters_mexican, file="Rdata/fgamSVC_clusters_mexican.Rdata")


# Optimal cluster number with GAP statistics
gaps <- clusGap(coefficient_matrix_inland, FUN = kmeans, K.max = 15, B=100)
par(mar=c(5,5,4,1))
plot(gaps, las=1, xpd=T, main='GAP statistics', xlab="Number of Clusters")

# Optimal cluster number with GAP statistics
gaps <- clusGap(coefficient_matrix_coast, FUN = kmeans, K.max = 15, B=100)
par(mar=c(5,5,4,1))
plot(gaps, las=1, xpd=T, main='GAP statistics', xlab="Number of Clusters")

# Optimal cluster number with GAP statistics
gaps <- clusGap(coefficient_matrix_mexican, FUN = kmeans, K.max = 15, B=100)
par(mar=c(5,5,4,1))
plot(gaps, las=1, xpd=T, main='GAP statistics', xlab="Number of Clusters")


# Rand index for various cluster numbers
ri <- numeric(length = length(c(2:15)))
for (i in c(2:15)){
  useCluster = clusters[[i]]$cluster
  ri[i-1] <- rand.index(useCluster, as.integer(as.factor(DougPres$POP)))
}
  
print("Maximum Rand Idex reached with cluster number:")
print(which.max(ri))
print("Rand Index:")
print(ri[10])

cfm <- confusion_matrix(as.integer(as.factor(DougPres$POP)), clusters[[6]]$cluster)
cfm

cluster <- clusters_coast[[6]]
useCluster <- cluster$clustering

fviz_cluster(cluster, data=coefficient_matrix_coast, ellipse.type = 'convex')

#==============#
# Overview fig #
#==============#

load(file="Rdata/fgamSVC_clusters_mexican.Rdata")
load(file="Rdata/fgamSVC_clusters_inland.Rdata")
load(file="Rdata/fgamSVC_clusters_coast.Rdata")

x <- c("maps", "ggplot2")
lapply(x, library, character.only=TRUE)
  
## PLOT DNA types ##

useCluster_i <- clusters_inland[[5]]$clustering
p1 = Ecotypes(DougScaledPres[DougScaledPres$VARS == "VAR_inland",], useCluster_i)
  
## PLOT Cluster ##
useCluster_c <- clusters_coast[[6]]$clustering
p2 = Ecotypes(DougScaledPres[DougScaledPres$VARS == "VAR_coast",], useCluster_c)

## PLOT Cluster ##
useCluster_m <- clusters_mexican[[1]]$clustering
p3 = Ecotypes(DougScaledPres[DougScaledPres$VARS == "VAR_mexican",], useCluster_m)

pdf(paste0("plots/varieties_cluster.pdf"), width = 9, height = 9)
ggpubr::ggarrange(p2, p1, p3, ncol = 2, nrow = 2, hjust=-4, labels= c("a", "b", "c", ""))
dev.off()    

