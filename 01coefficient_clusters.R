library(mgcv)
library(fossil)
library(cluster)
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


DougSample <- DougPres[sample(1:nrow(DougPres), size=500),]
  
# PCA 
par(mfrow=c(1,4))
pca.Doug <- prcomp((Douglas[,c(8:13)]))
pca.Dougs <- prcomp((DougScaled[,c(8:13)]))
pca.preds <- prcomp(predTerms)
pca.predss <- prcomp(scale(predTerms))
plot(pca.Doug)
plot(pca.Dougs)
plot(pca.preds)
plot(pca.predss)
pca.preds$rotation
pca.predss$rotation
pca.Doug$rotation


#==================================#
# Load a varying coefficient model #
#==================================#

load(file="Rdata/fgamSVC.Rdata")
load(file="Rdata/fgamSVC_predTerms.Rdata")
load(file="Rdata/fgamSVC_clusters.Rdata")
summary(fgamSVC)
AIC(fgamSVC)
xtable(summary(fgamSVC)$s.table, digits = 2)

mods <- c("ref", "basic", "non-standardized", "scaled-location", "neutral")
for (mod in mods){
  if (mod == "ref"){
    load("Rdata/fgamtrs.Rdata")
    p = "r"
   }else if (mod=="basic"){
    load("Rdata/fgamSVCtrs2.Rdata")
    p = 2
  }else if (mod=="non-standardized"){
    load("Rdata/fgamSVCtrs3.Rdata")
    p = 3
   }else if (mod=="scaled-location"){
    load("Rdata/fgamSVCtrs4.Rdata")
    p = 4
  }else if (mod=="neutral"){
    load("Rdata/fgamSVCtrsN.Rdata")
    p = "N"
  }
}

# predicted occurance probabilities
allPreds <- predict(fgamSVC, newdata=DougScaled, type="response")
# model coefficients at presence observations
predTerms <- predict(fgamSVC, newdata=DougScaled[DougScaled$PRES==1,], type="terms")
save(predTerms, file="Rdata/fgamSVC_predTerms.Rdata")


#========================#
# Clustering & Diagnosis #
#========================#

coefficient_matrix <- predTerms[,c(2,7)]

clusters <- list()
for (i in c(2:15)){
  clusters[[i]] <- pam(scale(coefficient_matrix), k=i, medoids = 'random', variant = 'faster')
}
save(clusters, file="Rdata/fgamSVC_clusters_scaled.Rdata")

# Variety clusters?
cluster <- clusters[[10]]
useCluster <- cluster$clustering
p4 = Ecotypes(DougScaledPres, useCluster)
p4 

# Optimal cluster number with GAP statistics
gaps <- clusGap(coefficient_matrix, FUN = kmeans, K.max = 15, B=100)
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
print(ri[9])

cfm <- confusion_matrix(as.integer(as.factor(DougPres$POP)), clusters[[6]]$cluster)
xtable(cfm, digits = 0)
inpercentageRows(cfm, rowSums(cfm))

cluster <- clusters[[6]]
useCluster <- cluster$clustering

fviz_cluster(cluster, data=coefficient_matrix, ellipse.type = 'convex',
             ggtheme = theme_minimal(), labelsize = 14)

#==============#
# Overview fig #
#==============#

x <- c("maps", "ggplot2")
lapply(x, library, character.only=TRUE)
  
## OCCURRANCE PROB #
p1 = occurance_predictions(DougScaled, allPreds)

# Observed OCCURRANCES
p2 = occurances(DougScaled)


## PLOT DNA types ##
DNAdata = get_DNA_data()
p3 = DNA_types(DougScaledPres, DNAdata)
  
## PLOT Cluster ##
p4 = Ecotypes(DougScaledPres, useCluster)

pdf(paste0("plots/Overview.pdf"), width = 9, height = 9)
ggpubr::ggarrange(p2, p1, p3, p4,ncol = 2, nrow = 2, hjust=-4, labels= c("a", "b", "c", "d"))
dev.off()    


#=========================#
# High-resoluted clusters #
#=========================#

cluster <- clusters[[10]]
useCluster <- cluster$clustering
p4 = Ecotypes(DougScaledPres, useCluster)
p4  

fviz_cluster(cluster, data=coefficient_matrix, ellipse.type = 'convex')
