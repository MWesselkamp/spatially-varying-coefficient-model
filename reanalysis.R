require(ggplot2)
require(mgcv)
require(maps)
require(raster)
require(rgdal)
require(randomForest)
require(dismo)
require(gbm)
require(clusteval)

#=================#
# stationary GAM  #
#=================#

Douglas <- read.csv("data/DF Plot Data (Norm_6190 Climate).csv")
Douglas$PPT_sm[which(Douglas$PPT_sm == -1)] <- NA # remove incorrect value 
Douglas <- na.omit(Douglas)
row.names(Douglas) <- NULL # reset the rownames to index
Douglas$TD2 <- Douglas$TD^2
Douglas$MWMT2 <- Douglas$MWMT^2
Douglas$PPT_sm2 <- Douglas$PPT_sm^2
Douglas$Tave_wt2 <- Douglas$Tave_wt^2
Douglas$Elev2 <- Douglas$Elev^2

DouglasScaled = Douglas
DouglasScaled[,(10:22)] = scale(DouglasScaled[,(10:22)])
DouglasScaledPres = DouglasScaled[which(DouglasScaled$PRES==1),]

# Take a subset eof the data for computational reasons
DouglasSample = Douglas[sample(nrow(Douglas), 10000),]

# Plot the spatial distribution of subsetted datapoints
x <- c("maps", "ggplot2")
lapply(x, library, character.only=TRUE)
countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)

ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="grey70") +
  geom_point(data=DouglasSample, aes(x=Long, y=Lat, col=factor(DouglasSample$PRES, levels = c(1,0))), size = 0.05) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=11, face="bold"), legend.text = element_text(size=11),
        legend.key=element_rect(fill = "white", size=0.5),
        legend.key.size = unit(5,"points"),
        plot.margin = margin(0.1,0.1,0.1,0.1, "cm")) + 
  guides(color = guide_legend(override.aes = list(size=3))) +
  scale_color_manual(values=c("#99CCFF", "#003366"), labels=c("present", "absent"),
                     name="Observation \nStatus")

# analyse the correlation between variables as a correlation plot and cluster
cormat = cor(as.matrix(DouglasSample[,10:36]))

require(corrplot)
corrplot(cormat, number.digits = 1)

require(Hmisc)
DouglasSample = DouglasSample[,-c(1:9,16,12,22,26,29,31)]
DouglasSample = DouglasSample[,-c(9:15)]
DouglasSample = DouglasSample[,-c(10:12, 14)]
DouglasSample = DouglasSample[,-c(11:13)]
DouglasSample = DouglasSample[,-c(3,4)]
DouglasSample = DouglasSample[,1:6]
DouglasSample = DouglasSample[,-1]

v <- as.formula(paste("~", names(DouglasSample), collapse="+"))
par(mfrow=c(1,1))
plot(varclus(v, data=DouglasSample))
abline(h=0.5)

vmat = varclus(v, data = DouglasSample)


## Rescale the smaller data set
DouglasSampleScaled = cbind(DouglasSample[,c(1:9)], DouglasSample$MWMT, DouglasSample$MWMT2, DouglasSample$PPT_sm, DouglasSample$PPT_sm2, DouglasSample$TD, DouglasSample$TD2, DouglasSample$PPT_wt, DouglasSample$MDMP)
names(DouglasSampleScaled)[10:17] = c("MWMT","MWMT2", "PPT_sm", "PPT_sm2", "TD", "TD2", "PPT_wt", "MDMP")

DouglasSampleScaled$PPT_wt2 = DouglasSampleScaled$PPT_wt^2
DouglasSampleScaled$MDMP2 = DouglasSampleScaled$MDMP^2


DouglasSampleScaled[,c(10:19)] = scale(DouglasSampleScaled[,c(10:19)])


## Run stationary GAMs with the leftover non-correlated variables ##
require(mgcv)

## Run simplest GAM
fm = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)

save(fm, file="Rdata/fgamStatSimple.Rdata")
summary(fm)
plot(fm)

# Use all five variables and quadrativ terms

fm = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100) + s(PPT_wt, k=100) + s(PPT_wt2, k=100) + s(MDMP, k=100) + s(MDMP2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)

save(fm, file="Rdata/fgamStatSimple2.Rdata")
summary(fm)
plot(fm)


## Variable importance ##
# only with non-correlated variables

## RandomForest
rf = randomForest(as.factor(PRES) ~ TD + MWMT + PPT_wt  + PPT_sm + MDMP, data=DouglasSample)
varImpPlot(rf)
save(rf, file="Rdata/randomForest.Rdata")

## Boosted regression trees
# same variables. Test for interactions (reviewer recommendations)

brt1 = gbm.step(data=DouglasSample, gbm.x = c(11,15,17,18,19), gbm.y = 7, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)
save(brt1, file="Rdata/regressiontrees.Rdata")

summary(brt1)
relative.influence(brt1)

## Another simple stationary gam

fm = gam(PRES ~ s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100) + s(PPT_wt, k=100) + s(PPT_wt2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)

save(fm, file="Rdata/fgamStatSimple3.Rdata")
summary(fm)
plot(fm)


## Spatially variable coefficient model with five predictors.



#=============#
# Neutral Map #
#=============#

# Create dataset with randomized presences and absences.
DouglasScaledN = DouglasScaled
DouglasScaledN$PRES = sample(DouglasScaledN$PRES)

# fit svcm to neutral data
fgamSVCtsrN <- gam(PRES ~ s(y, x, k=100) + s(y, x, by=TD, k=100) + s(y, x, by=TD2, k=100) + s(y, x, by= PPT_sm, k=100) + s(y, x, by=PPT_sm2, k=100) + s(y, x, by=MWMT, k=100) + s(y, x, by=MWMT2, k=100), method="ML", family=binomial, control=gam.control(maxit=1000), data=DouglasScaledN)

save(fgamSVCtsrN, file="Rdata/fgamSVCtsrN.Rdata")

summary(fgamSVCtsrN) # very bad model.

##k-means clustering##

# predict with neutral model to only true presences.
predTermsN = predict(fgamSVCtsrN, newdata = DougScaledPres, type="terms")

# cluster neutral model coefficients
clustered.predictionsN <- vector(mode = "list", length = 15)
for (i in 2:16){
  clustered.predictionsN[[i]] <- kmeans(predTermsN, centers=i, nstart=25)
}
save(clustered.predictionsN, file="Rdata/clusteredpredictionsN.Rdata")


# Plot clusters#
countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 180)

load("states.Rdata")
load("states_mex.Rdata")
load("states_ca.Rdata")

useCluster <- clustered.predictionsN[[12]]$cluster
col.6 <- c("#F3DF6C", "#CEAB07", "#798E87","#C93312", "#CCC591","#C27D38")
col.13 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')

p4 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="grey70") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey85") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey85")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey85") +
  geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=factor(useCluster)), size=0.01) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=11, face="bold"), legend.text = element_text(size=11),
        legend.key=element_rect(fill = "white", size=0.5),
        plot.margin = margin(0.1,0.1,0.1,0.1, "cm")) + 
  guides(color = guide_legend(override.aes = list(size=3))) + 
  scale_color_manual(name = "Cluster", values = col.13)

pdf("figures/p4_12_neutral.pdf", width=6)
p4
dev.off()

# Plot model coefficients #

require(classInt)
quantiles.1 <- classIntervals(scale(predTermsN)[,2], 25, style="quantile")
legend_labels.1 <- round(quantiles.1$brks[seq(1,26,5)], 1)

predTerms_cut.1 <- cut(scale(predTermsN)[,2], breaks = quantiles.1$brks, include.lowest = TRUE)
#cols <- colorRampPalette(c("#003366", "#99CCFF"))(30)
preds <- (predTermsN)[,1]

pTD <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="grey55") + 
  #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey68") + 
  #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey68")  + 
  #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey68") +
  geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.1)), size = 0.05)  + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="TD ",low="#106607", high="#75F567", breaks= c(0,5,10,15,20,25),labels=as.character(legend_labels.1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"),legend.text = element_text(size=8),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

pdf("figures/coefs_TD.pdf", width=6)
pTD
dev.off()


############################
# Confusion with genotypes #
############################

confusion_matrix = function(genotype, ecotype){
  mat = matrix(0, nrow=length(unique(ecotype)), ncol=length(unique(genotype)))
  for (i in 1:length(ecotype)){
    mat[ecotype[i], genotype[i]] = mat[ecotype[i], genotype[i]]+1
  }
  return(mat)
}

## Confusion of neutral clusters and ecotypic clusters
load("clusters2.Rdata")
nclust <- clustered.predictionsN[[12]]$cluster
eclust <- clusters2[[13]]$cluster

cn = confusion_matrix(nclust, eclust)

## Cluster similarity: Ecotypes and neutral model ##
require(clusteval)
cluster_similarity(nclust, eclust, similarity = c("jaccard"))
# Rand's: [14 clusters] 0.842104, [12 Clusters] 0.8208949, [6 clusters] 0.6966142
# Jaccard: [14 clusters] 0.1456683, [12 clusters] 0.1583153, [6 clusteres] 0.2151505

## Cluster similartiy: Ecotypes/neutral model and genotypes ##
cluster_similarity(as.numeric(DouglasScaledPres$POP)-1, eclust, similarity=c("rand"))
#[1] 0.7173828
cluster_similarity(as.numeric(DouglasScaledPres$POP)-1, nclust, similarity=c("rand"))
#[1] 0.7120244


## Confusion matrix: 12 clusters in 6 genotypes ## 
cm12 = confusion_matrix(genotype = as.numeric(DouglasScaledPres$POP)-1, ecotype = eclust)
cm12
cluster_similarity(as.numeric(DouglasScaledPres$POP)-1, eclust, similarity=c("rand"))
# Rand's: [14 clusters] 0.772675, [12 clusters] 0.7888955, [6 clusters] 0.7173828


# percentage of ecotypic presences assigned to neutral cluster

inpercentageRows = function(m, sums){
  for(i in 1:nrow(m)){
    m[i,] = round((m[i,]/sums[i])*100, 2)
  }
  return(m)
}

inpercentageRows(cn, sums=rowSums((cn)))

# percentage of neutral clusters assigned to ecotypes

inpercentageCols = function(m, sums){
  for(i in 1:ncol(m)){
    m[,i] = round((m[,i]/sums[i])*100, 2)
  }
  return(m)
}

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

wss = numeric(14)
for (i in 1:14){
  wss[i] = clusters2[[i+1]]$tot.withinss  
}
plot(1:14, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


#########################################
# Validation with genotypes (Wei et al.)#
#########################################

DNAdata = read.csv("~/ScAdditional/PaperSVCM/data/Doug-Fir DNA data (Wei et al, 2011).csv")


col.6 <- c("#F3DF6C", "#CEAB07", "#798E87","#C93312", "#CCC591","#C27D38")
pDNA <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="grey70") + 
  #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey85") + 
  #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey85")  + 
  #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey85") +
  geom_point(data=DNAdata, aes(x=Longitude, y=Latitude, col=factor(Genotype)), size=1) +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_manual(values = col.6, name = "Genotype")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"),legend.text = element_text(size=8),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))


## Extract BIOCLIM variables for 44 genotypes populations.

bioclim = getData("worldclim", var="bio", res=10)
mypoints = data.frame(long=DNAdata$Longitude, lat=DNAdata$Latitude)
myvars = as.data.frame(extract(bioclim, mypoints))
myvars_small = data.frame(TD = myvars$bio7/10, MWMT = myvars$bio10/10, PPT_sm = myvars$bio18)

## prepare data for SVCM
DNAdata = cbind(DNAdata, myvars_small)
DNAdata$TD2 = DNAdata$TD^2
DNAdata$MWMT2 = DNAdata$MWMT^2
DNAdata$PPT_sm2 = DNAdata$PPT_sm^2

DNAdata_scaled = DNAdata
DNAdata_scaled[,(12:17)] = scale(DNAdata_scaled[,(12:17)])

## convert latitude and longitude to UTM coordinates.
xy = data.frame(x = DNAdata_scaled[,5], y = DNAdata_scaled[,4])
coordinates(xy) = c("x", "y")
proj4string(xy) = CRS("+proj=longlat +datum=WGS84")
res = spTransform(xy, "+proj=utm +datum=WGS84")
DNAdata_scaled = cbind(DNAdata_scaled, res@coords)

## Predict with fitted model to populations, apply k-means and compare genotypes and ecotypes.

load("~/Sc_Master/SVCM/Rdata/fgamSVCtrs2.Rdata")
load("~/Sc_Master/SVCM/Rdata/clusters2.Rdata")

preds = predict(fgamSVCtrs2, DNAdata_scaled, type="response")
predTermsDNA = predict(fgamSVCtrs2, DNAdata_scaled, type="terms")

Ecoclusters_Wei = kmeans(predTermsDNA, centers = 6, nstart = 25)

## similarity of ecotypes and genotypes.
confusion_matrix(DNAdata$Genotype, Ecoclusters_Wei$cluster)
cluster_similarity(DNAdata$Genotype, Ecoclusters_Wei$cluster, "rand")


## Show ecotypic classifications for both, the large dataset (ecotypes in grey colors) and the genotyped populations (ecotypes in rainbow colors).
Ecoclusters_all = clusters2[[6]]$cluster+70

greys = grey.colors(6, start=0.2, end=0.8)
rainbows = rainbow(n=6)
col.12 <- c(rainbows, greys)
pDNA <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey85") + 
  #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey85")  + 
  #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey85") +
  geom_point(data = DouglasScaledPres, aes(x=Long, y=Lat, col=factor(Ecoclusters_all)), size=1.5, alpha=0.5) + 
  geom_point(data=DNAdata_scaled, aes(x=Longitude, y=Latitude, col=factor(Ecoclusters_Wei$cluster)), size=1) +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_manual(values = col.12, name = "Genotype")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none", 
        legend.title = element_text(size=8, face="bold"),legend.text = element_text(size=8),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
pDNA



