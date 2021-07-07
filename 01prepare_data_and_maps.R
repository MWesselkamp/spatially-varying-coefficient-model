require(mgcv)

setwd("~/Documents/Projects/Spatially-varying-coefficient-model")

source("utils.R")
source("plots.R")

## Prepare data ##
##################

Douglas <- get_Douglas_data()

# only presences - unscaled
DougPres <- Douglas[Douglas$PRES==1,]

# data set complete - scaled
DougScaled <- Douglas
DougScaled[,c(10:22)] <- scale(Douglas[,c(10:22)])

# data set presences - scaled
DougScaledPres <- DougScaled[DougScaled$PRES==1,]

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


## Spatially varying coefficient model ##
#########################################

fgamSVCtsr <- gam(PRES ~ s(y, x, k=100) + s(y, x, by=TD, k=100) + 
                    s(y, x, by=TD2, k=100) + s(y, x, by= PPT_sm, k=100) + 
                    s(y, x, by=PPT_sm2, k=100) + s(y, x, by=MWMT, k=100) + 
                    s(y, x, by=MWMT2, k=100), 
                  method="ML", family=binomial, control=gam.control(maxit=1000), 
                  data=Douglas)

save(fgamSVCtrs, file="Rdata/fgamSVCtsr3.Rdata")
summary(fgamSVCtsr)

# grs = c("Lat", "Long"): load("Rdata/fgamSVCtsr.Rdata")
# grs = c("x", "y"):
load("Rdata/fgamSVCtsr2.Rdata")

#TABLE OF MODEL SUMMARY
anova(fgamSVCtrs2)
xtable(summary(fgamSVCtrs2)$s.table, digits = 2)

# predicted occurance probabilities
allPreds <- predict(fgamSVCtsr, newdata=DougScaled, type="response")
save(allPreds, file="Rdata/allPreds2.Rdata")

# model coefficients at presence observations
predTerms <- predict(fgamSVCtsr, newdata=DougScaled[DougScaled$PRES==1,], type="terms")
save(predTerms, file="Rdata/predTerms2.Rdata")

load("Rdata/predTerms2.Rdata")
summary(predTerms)
load("Rdata/predTermsN.Rdata")
summary(predTermsN)

predTermsF <- sqrt((predTerms-predTermsN)^2)
#===================#
# kmeans clustering #
#===================#

clustered.predictions = kmeans_clustering(predTerms, 15)

which.max(sapply(clustered.predictions, function(x) x$tot.withinss/x$totss))

useCluster = clustered.predictions[[6]]$cluster
xtable::xtable(confusion_matrix(DougPres$POP, useCluster), digits=0)

#==============#
# Overview fig #
#==============#

# LOAD TO REPRODUCE FIGURES
load("Rdata/clusters2.Rdata")
useCluster = clusters2[[3]]$cluster
load("Rdata/allPreds2.Rdata")

x <- c("maps", "ggplot2")
lapply(x, library, character.only=TRUE)


## PLOT OCCURRANCE PROB #

p1 = occurance_predictions(DougScaled, allPreds)

# Observed occurances

p2 = occurances(DougScaled)

## PLOT DNA types ##

DNAdata = get_DNA_data()

p3 = DNA_types(DougScaledPres, DNAdata)

## PLOT Cluster ##

p4 = Ecotypes(DougScaledPres, useCluster, n_colors = 6)

 
pdf("figures/Overview.pdf", width = 9, height = 9)
ggpubr::ggarrange(p2, p1, p3, p4,ncol = 2, nrow = 2, hjust=-4, labels= c("a", "b", "c", "d"))
dev.off()    

 