require(ggplot2)
require(classInt)
require(mgcv)
require(clusteval)

source("utils.R")
source("plots.R")

Douglas <- get_Douglas_data()

DougPres <- Douglas[Douglas$PRES==1,]

# data set complete - scaled
DougScaled <- Douglas
DougScaled[,c(10:41)] <- scale(Douglas[,c(10:41)])

# data set presences - scaled
DougScaledPres <- DougScaled[DougScaled$PRES==1,]

load("Rdata/fgamSVCtsr2.Rdata")
load("Rdata/fgamSVCtsrN.Rdata")
load("Rdata/predTerms2.Rdata")
load("Rdata/states.Rdata")
load("Rdata/states_ca.Rdata")
load("Rdata/states_mex.Rdata")

countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)

predTermsN = predict(fgamSVCtsrN, newdata = DougScaledPres, type="terms")
#save(predTermsN, file="Rdata/predTermsN.Rdata")

summary(predTermsN)
summary(scale(predTerms))
summary(scale(predTermsN))

# coefficient correction
coeffs_mod = (sqrt((scale(predTerms)-scale(predTermsN))^2))
coefficient_maps(coeffs_mod)

clustered.predictions = kmeans_clustering(coeffs_mod, n_clusters=15)
save(clustered.predictions, file="Rdata/clusters_EN.Rdata")

useCluster = clustered.predictions[[2]]$cluster

# Load the data frame that contains the computed clusterwise responses.
load("Rdata/effectplots_EN.Rdata")

pdf("figures/Effects_EN.pdf", height = 12, width = 7)
effectplots(df, "effectplots_EN")
dev.off()

DougScaledPres = DougScaled[DougScaled$PRES==1,]

# Plot Ecotypes for choice of two clusters.
pdf("figures/Ecotypes_EN.pdf", height=7, width = 7)
Ecotypes(DougScaledPres, clustered.predictions[[2]]$cluster)
dev.off()

confusion_matrix(DougScaledPres$POP, useCluster)

## Cluster similarity: Ecotypes and neutral model ##
cluster_similarity(DougScaledPres$POP, useCluster, similarity = c("rand")) #0.7062538
load("Rdata/clustersN.Rdata")
clustersN = clustered.predictionsN[[6]]$cluster
cluster_similarity(clustersN, useCluster, similarity = c("rand")) #0.6917944
load("Rdata/clusters2.Rdata")
clustersE = clusters2[[2]]$cluster
cluster_similarity(clustersE, useCluster, similarity = c("rand")) # 0.7042829