require(dismo)
require(gbm) # boosted regression trees
require(mda) # discriminant function analysis


source("plots.R")
source("utils.R")

#########################################
# Validation with genotypes (Wei et al.)#
#########################################

DNAdata <- get_DNA_data()

load("Rdata/clusters.Rdata")
Ecoclusters_all = clusters[[6]]$cluster
Ecotypes(DougScaledPres, Ecoclusters_all)

DNA_types(DougScaledPres, DNAdata)

## Predict with fitted model to populations, apply k-means and compare genotypes and ecotypes.

load("Rdata/fgamSVC.Rdata")

DNAdata_scaled <- scale(DNAdata[,12:17])
preds = predict(fgamSVC, DNAdata_scaled, type="response")
predTermsDNA = predict(fgamSVC, DNAdata_scaled, type="terms")
Ecoclusters_Wei = kmeans(predTermsDNA, centers = 6, nstart = 25)

## similarity of ecotypes and genotypes.
#confusion_matrix(DNAdata$Genotype, Ecoclusters_Wei$cluster)
cluster_similarity(Ecoclusters_Wei$cluster, DNAdata$Genotype, "rand")
cluster_similarity(Ecoclusters_Wei$cluster, clusters[[6]]$cluster, "rand")

##################################
# Discriminant Function Analysis #
##################################

Ecotypes = clusters2[[6]]$cluster
DouglasScaledPres = cbind(DouglasScaledPres, Ecotypes)

fda1 = fda(Ecotypes ~ TD + TD2 + PPT_sm + PPT_sm2 + MWMT + MWMT2, data=DouglasScaledPres)
predictions = predict(fda1, DouglasScaledPres)
mean(predictions == Ecotypes)

fda1$percent.explained
