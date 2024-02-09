require(ggplot2)
require(mgcv)
require(maps)
require(raster)
require(rgdal)
require(randomForest)
require(dismo)
require(gbm) # boosted regression trees
require(clusteval) # cluster similarity
require(mda) # discriminant function analysis
require(corrplot) # correlation plot
require(Hmisc) # varclus plot

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

confusion_matrix(nclust, eclust)
rand.index(nclust, eclust)


## Cluster similarity: Ecotypes/neutral model and Rehfeldt genotypes ##
confusion_matrix(as.numeric(DNAdata$Genotype), eclust)
inpercentageRows(cn, sums=rowSums((cn)))
inpercentageCols(cn, sums=colSums((cn)))

# confusion of "true" ecotypes with genoptypes
cm1 = confusion_matrix(as.numeric(DouglasScaledPres$POP)-1, eclust)
inpercentageRows(cm1, rowSums(cm1))
# confusion of "neutral" ecotypes with genoptypes
cm2 = confusion_matrix(as.numeric(DouglasScaledPres$POP)-1, nclust)
inpercentageRows(cm2, rowSums(cm2))


################################
# Determine number of clusters #
################################

plot_wss(clusters2)

#########################################
# Validation with genotypes (Wei et al.)#
#########################################

DNAdata <- get_DNA_data()

load("Rdata/clusters2.Rdata")
Ecoclusters_all = clusters2[[6]]$cluster

RColorBrewer::brewer.pal(6, "Set1")
greys = grey.colors(6, start=0.2, end=1.0)
col.6 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")
col.d = c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20", "#bd0026") # same colours as in Paper overview
col.12 = c(col.6, greys)

pDNA = ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data = DougScaledPres, aes(x=Long, y=Lat, col=factor(POP)), size=0.8, alpha=1) + 
  geom_point(data=DNAdata, aes(x=Longitude, y=Latitude, col="white", fill=factor(Genotype)), size=1.2, shape=21) +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_color_manual(values = c(greys, "black"), guide = FALSE) +
  scale_fill_manual(values = c(col.d), name = "Genotypes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())  +
  guides(fill = guide_legend(override.aes = list(size=2)))



## Predict with fitted model to populations, apply k-means and compare genotypes and ecotypes.

load("Rdata/fgamSVCtsr2.Rdata")

DNAdata_scaled <- scale(DNAdata[,12:17])
preds = predict(fgamSVCtrs2, DNAdata_scaled, type="response")
predTermsDNA = predict(fgamSVCtrs2, DNAdata_scaled, type="terms")
Ecoclusters_Wei = kmeans(predTermsDNA, centers = 6, nstart = 25)

## similarity of ecotypes and genotypes.
#confusion_matrix(DNAdata$Genotype, Ecoclusters_Wei$cluster)
cluster_similarity(Ecoclusters_Wei$cluster, DNAdata$Genotype, "rand")
cluster_similarity(Ecoclusters_Wei$cluster, clusters2[[6]]$cluster, "rand")


## Show ecotypic classifications for both, the large dataset (ecotypes in grey colors) and the genotyped populations (ecotypes in rainbow colors).

RColorBrewer::brewer.pal(6, "Dark2")
col.6.2 = c("#F3DF6C", "#CEAB07", "#798E87","#C93312", "#CCC591","#C27D38")
col.l = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF") # same colours as in effect plots.

pECO = ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data = DouglasScaledPres, aes(x=Long, y=Lat, col=factor(Ecoclusters_all)), size=0.8, alpha=1) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) +
  scale_color_manual(values = c(col.l), name="Ecotypes") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=2)))


## Combine with Raw data and model predicitons to new overview plot ##

load("~/ScAdditional/PaperSVCM/Rdata/allPreds2.Rdata")

p1 = ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data=DouglasScaled, aes(x=Long, y=Lat, col=allPreds), size = 0.7)  + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="SVCM \nPredictions", low="#003366", high="#99CCFF", breaks = c(0.25, 0.5, 0.75), labels = c("0.25", "", "0.75")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) 

p2 = ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data=DouglasScaled, aes(x=Long, y=Lat, col=factor(DouglasScaled$PRES, levels = c(1,0))), size = 0.7) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold"), legend.text = element_text(size=11),
        legend.key=element_rect(fill = "white", size=0.5),
        legend.key.size = unit(5,"points")) + 
  guides(color = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values=c("#99CCFF", "#003366"), labels=c("present", "absent"),
                     name="Observation \nStatus")


pdf("figures/OverviewR.pdf")
ggpubr::ggarrange(p2, p1, pDNA, pECO,ncol = 2, nrow = 2, hjust=-4)
dev.off()  


##################################
# Discriminant Function Analysis #
##################################

Ecotypes = clusters2[[6]]$cluster
DouglasScaledPres = cbind(DouglasScaledPres, Ecotypes)

fda1 = fda(Ecotypes ~ TD + TD2 + PPT_sm + PPT_sm2 + MWMT + MWMT2, data=DouglasScaledPres)
predictions = predict(fda1, DouglasScaledPres)
mean(predictions == Ecotypes)

fda1$percent.explained
