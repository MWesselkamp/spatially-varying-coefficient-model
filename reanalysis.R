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
 
#=========================================#
# Variable selection with stationary GAM  #
#=========================================#

Douglas <- read.csv("data/DF Plot Data (Norm_6190 Climate).csv")
Douglas$PPT_sm[which(Douglas$PPT_sm == -1)] <- NA # remove incorrect value 
Douglas <- na.omit(Douglas)
row.names(Douglas) <- NULL # reset the rownames to index
Douglas$TD2 <- Douglas$TD^2
Douglas$MWMT2 <- Douglas$MWMT^2
Douglas$PPT_sm2 <- Douglas$PPT_sm^2
Douglas$PPT_wt2 <- Douglas$PPT_wt^2
Douglas$MDMP2 <- Douglas$MDMP^2

DouglasScaled = Douglas
DouglasScaled[,c(10:22, 37:43)] = scale(DouglasScaled[,c(10:22,37:43)])
DouglasScaledPres = DouglasScaled[which(DouglasScaled$PRES==1),]

# Take a subset eof the data for computational reasons
DouglasSample = Douglas[sample(nrow(Douglas), 10000),]

# Plot the spatial distribution of subsetted datapoints

countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)
load("Rdata/states.Rdata")
load("Rdata/states_mex.Rdata")
load("Rdata/states_ca.Rdata")

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
corrplot(cormat, number.digits = 1)

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

# Run simplest GAM
fm = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)

save(fm, file="Rdata/fgamStatSimple.Rdata")
summary(fm)
plot(fm)

# Use all five variables and quadrativ terms

fm = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100) + s(PPT_wt, k=100) + s(PPT_wt2, k=100) + s(MDMP, k=100) + s(MDMP2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasScaled)

save(fm, file="Rdata/fgamStatSimple2.Rdata")

load("Rdata/reanalysis/fgamStatSimple2.Rdata")
summary(fm)
plot(fm)

fm1 = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100) + s(PPT_wt, k=100) + s(PPT_wt2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)

fm2 = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)


## Variable importance ##
# only with non-correlated variables

## RandomForest
rf = randomForest(as.factor(PRES) ~ TD + MWMT + PPT_wt  + PPT_sm + MDMP, data=DouglasSample)
randomForest::varImpPlot(rf)
save(rf, file="Rdata/randomForest.Rdata")

## Boosted regression trees
# same variables. Test for interactions (reviewer recommendations)

brt1 = gbm.step(data=DouglasSample, gbm.x = c(11,15,17,18,19), gbm.y = 7, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)
save(brt1, file="Rdata/regressiontrees.Rdata")

summary(brt1)
sort(relative.influence(brt1))

## Another simple stationary gam

fm = gam(PRES ~ s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100) + s(PPT_wt, k=100) + s(PPT_wt2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)

save(fm, file="Rdata/fgamStatSimple3.Rdata")
summary(fm)
plot(fm)


## Spatially variable coefficient model with five predictors?


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
load("Rdata/reanalysis/clusteredpredictionsN.Rdata")
useCluster <- clustered.predictionsN[[12]]$cluster
col.6 = wesanderson::wes_palette("Cavalcanti1", 6, type = "continuous") 
col.12 = wesanderson::wes_palette("Cavalcanti1", 12, type = "continuous") 

p4 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=factor(useCluster)), size=1.0) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_color_manual(values = c(col.12), name = "Ecotypes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())  +
  guides(color = guide_legend(override.aes = list(size=2)))

pdf("figures/neutralClusters12.pdf")
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
# Cluster similarities     #
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
eclust <- clusters2[[6]]$cluster

confusion_matrix(nclust, eclust)

## Cluster similarity: Ecotypes and neutral model ##
cluster_similarity(nclust, eclust, similarity = c("rand"))
# Rand's: [12 Clusters] 0.8208949, [6 clusters] 0.6966142

## Cluster similartiy: Ecotypes/neutral model and Rehfeldt genotypes ##
cluster_similarity(as.numeric(DougScaledPres$POP)-1, eclust, similarity=c("rand"))
#[6 clusters] 0.7173828 [12 clusters] 0.7888955
cluster_similarity(as.numeric(DougScaledPres$POP)-1, nclust, similarity=c("rand"))
#[6 clusters] 0.7120244 [12 clusters] 0.7493401


## Cluster similarity: Ecotypes/neutral model and Rehfeldt genotypes ##
confusion_matrix(as.numeric(DNAdata$Genotype), eclust)


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
  geom_point(data = DouglasScaledPres, aes(x=Long, y=Lat, col=factor(POP)), size=0.8, alpha=1) + 
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


## Extract BIOCLIM variables for 44 genotypes populations.

bioclim = getData("worldclim", var="bio", lon=DNAdata$Longitude, lat=DNAdata$Latitude,res=0.5)
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
DNAdata = cbind(DNAdata, res@coords)

save(DNAdata, file="Rdata/DNAdata.Rdata")
save(DNAdata_scaled, file="Rdata/DNAdata_scaled.Rdata")

## Predict with fitted model to populations, apply k-means and compare genotypes and ecotypes.

load("Rdata/fgamSVCtsr2.Rdata")

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

#====================#
# Graphical Abstract #
#====================#

useCluster <- clusters2[[6]]$cluster
col.l <-  c("#F3DF6C", "#CEAB07", "#798E87","#C93312", "#CCC591","#C27D38") # coffee colors

p4 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=factor(useCluster)), size=0.8) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_color_manual(values = c(col.l), name = "Ecotypes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())  +
  guides(color = guide_legend(override.aes = list(size=2)))

pdf("figures/graphicalAbstract.pdf")
p4
dev.off()
