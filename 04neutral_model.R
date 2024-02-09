require(ggplot2)
require(mgcv)
require(maps)
require(raster)

source("utils.R")
source("plots.R")

load("Rdata/states.Rdata")
load("Rdata/states_mex.Rdata")
load("Rdata/states_ca.Rdata")

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


#===============#
# Neutral Model #
#===============#

# Create dataset with randomized presences and absences.
DouglasScaledN = DougScaled
DouglasScaledN$PRES = sample(DouglasScaledN$PRES)

# fit svcm to neutral data
fgamSVCneutral <- gam(PRES ~ s(y, x, k=100) + s(y, x, by=TD, k=100) + s(y, x, by=TD2, k=100) + s(y, x, by= PPT_sm, k=100) + s(y, x, by=PPT_sm2, k=100) + s(y, x, by=MWMT, k=100) + s(y, x, by=MWMT2, k=100), method="ML", family=binomial, control=gam.control(maxit=1000), data=DouglasScaledN)
save(fgamSVCneutral, file="Rdata/fgamSVCneutral.Rdata")
summary(fgamSVCneutral) # very bad model.

xtable(summary(fgamSVCneutral)$s.table, digits = 2)
AIC(fgamSVCneutral)
##k-means clustering##

# predict with neutral model to only true presences.
predTermsN = predict(fgamSVCneutral, newdata = DougScaledPres, type="terms")

# cluster neutral model coefficients
coefficient_matrix <- predTermsN[,c(2,7)]

clusters <- list()
for (i in c(2:15)){
  clusters[[i]] <- pam(scale(coefficient_matrix), k=i, medoids = 'random', variant = 'faster')
}
save(clusters, file="Rdata/fgamSVC_clusters_neutral.Rdata")

# Plot clusters#
load("Rdata/fgamSVC_clusters_neutral.Rdata")
useCluster <- clusters[[10]]$cluster

p4 <- Ecotypes(DougScaledPres, useCluster)

pdf("plots/neutralClusters10.pdf")
p4
dev.off()

# Plot model coefficients #
#coefficient_maps()

############################
# Cluster similarities     #
############################
useCluster <- clusters[[6]]$cluster
cm <- confusion_matrix(useCluster, as.numeric(as.factor(DougScaledPres$POP)))
xtable(cm, digits = 0)
rand.index(useCluster, as.numeric(as.factor(DougScaledPres$POP)))

useCluster <- clusters[[10]]$cluster
rand.index(useCluster, as.numeric(as.factor(DougScaledPres$POP)))

## Confusion of neutral clusters and ecotypic clusters
load("Rdata/fgamSVC_clusters_neutral.Rdata")
nclust <- clusters[[6]]$cluster
load("Rdata/fgamSVC_clusters.Rdata")
eclust <- clusters[[6]]$cluster

cn <- confusion_matrix(nclust, eclust)
rand.index(nclust, eclust)

# confusion of "neutral" ecotypes with ecotypes
inpercentageRows(cn, sums=rowSums((cn)))
inpercentageCols(cn, sums=colSums((cn)))


