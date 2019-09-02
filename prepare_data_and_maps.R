setwd("~/ScAdditional/PaperSVCM")

## Prepare data ##
##################

Douglas <- read.csv("data/DF Plot Data (Norm_6190 Climate).csv")
Douglas$PPT_sm[which(Douglas$PPT_sm == -1)] <- NA # remove incorrect value 
Douglas <- na.omit(Douglas)
row.names(Douglas) <- NULL # reset the rownames to index

Douglas <- Douglas[,c(1:5,7, 9, 11, 15, 17)]
Douglas$TD2 <- Douglas$TD^2
Douglas$MWMT2 <- Douglas$MWMT^2
Douglas$PPT_sm2 <- Douglas$PPT_sm^2

# only presences - unscaled
DougPres <- Douglas[Douglas$PRES==1,]

# data set complete - scaled
DougScaled <- Douglas
DougScaled[,c(8:13)] <- scale(Douglas[,c(8:13)])

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
require(mgcv)

cores <- detectCores() -1
c1 <- makeCluster(cores)
registerDoSNOW(c1)

fgamSVCtrs <- gam(PRES ~ s(Lat, Long, k=100) + s(Lat, Long, by=TD, k=100) + s(Lat, Long, by=TD2, k=100) + s(Lat, Long, by= PPT_sm, k=100) + s(Lat, Long, by=PPT_sm2, k=100) + s(Lat, Long, by=MWMT, k=100) + s(Lat, Long, by=MWMT2, k=100), method="ML", family=binomial, control=gam.control(maxit=1000), data=DougScaled)

save(fgamSVCtrs, file="Rdata/fgamSVCtrs.Rdata")
stopCluster(c1)
summary(fgamSVCtsr)
fgamSVCtsr$boundary

load("Rdata/fgamSVCtsr.Rdata")


## Spatially varying coefficient model with x/y ##
##################################################

fgamSVCtrs2 <- gam(PRES ~ s(y, x, k=100) + s(y, x, by=TD, k=100) + s(y, x, by=TD2, k=100) + s(y, x, by= PPT_sm, k=100) + s(y, x, by=PPT_sm2, k=100) + s(y, x, by=MWMT, k=100) + s(y, x, by=MWMT2, k=100), method="ML", family=binomial, control=gam.control(maxit=1000), data=DougScaled)

save(fgamSVCtrs2, file="Rdata/fgamSVCtsr2.Rdata")

load("Rdata/fgamSVCtsr2.Rdata")

fgamSVCtsr <- fgamSVCtrs2
rm(fgamSVCtrs2)

#TABLE OF MODEL SUMMARY
anova(fgamSVCtsr)
xtable(summary(fgamSVCtsr)$s.table, digits = 2)


# all model predictions #
allPreds <- predict(fgamSVCtsr, newdata=DougScaled, type="response")
save(allPreds, file="Rdata/allPreds2.Rdata")

# model coefficients at presence observations

predTerms <- predict(fgamSVCtsr, newdata=DougScaled[DougScaled$PRES==1,], type="terms")
save(predTerms, file="Rdata/predTerms2.Rdata")

load("predTerms2.Rdata")
summary(predTerms)

## kmeans clustering ###

clustered.predictions <- vector(mode = "list", length = 30)
for (i in 2:length(clustered.predictions)+1){
  clustered.predictions[[i]] <- kmeans(predTerms, centers=i, nstart=25)
}
wss = numeric(length(clustered.predictions))
for (i in 3:31){
  wss[i] = clustered.predictions[[i]]$tot.withinss/clustered.predictions[[i]]$totss
}
plot(3:31, wss[3:31], type="b", ylim=c(0,0.5), xlab="Number of Clusters",
     ylab="Total within groups sum of squares")


clustered.predictions.scaled <- vector(mode = "list", length = 15)
for (i in 2:16){
  clustered.predictions.scaled[[i]] <- kmeans(scale(predTerms), centers=i, nstart=25)
}
save(clustered.predictions.scaled, file="clusters2scaled.Rdata")

which.max(sapply(clustered.predictions, function(x) x$tot.withinss/x$totss))

# LOAD TO REPRODUCE FIGURES
load("Rdata/clusters2.Rdata")
clustered.predictions <- clusters2
load("Rdata/allPreds2.Rdata")

### SOME MAPS ###
#################

x <- c("maps", "ggplot2")
lapply(x, library, character.only=TRUE)

countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)


#---- get states borders for USA, Canada and Mexico
provinces_ca <- map("worldHires", "Canada")
states <- raster::getData("GADM", country="USA", level=1)
save(states, file="output/states.Rdata")

ca <- raster::getData("GADM", country="Canada", level=1)
provinces <- c("British Columbia", "Alberta", "Saskatchewan")
ca_smaller <- ca[ca$NAME_1 %in% provinces,] # only use relevant canadian states to save memory

mex <- raster::getData("GADM", country="Mexico", level=1)

ca_smaller@data$id <- rownames(ca_smaller@data) # define an explicit relationship between data and polygons associated with the data
ca_provinces <- fortify(ca_smaller, region = "id") # convert class of spdf into df
ca_provincesDF <- merge(ca_provinces, ca_smaller@data, by="id") # merge if some important information is missing
save(ca_provinces, file="output/ca_provinces.Rdata")

mex@data$id <- rownames(mex@data)
mex_states <- fortify(mex, region="id") 
save(mex_states, file="output/mex_states.Rdata")

## PLOT OCCURRANCE PROB ##
require(ggplot2)
load("Rdata/states.Rdata")
load("Rdata/states_ca.Rdata")
load("Rdata/states_mex.Rdata")


p1 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="grey70") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey85") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey85")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey85") +
  geom_point(data=DougScaled, aes(x=Long, y=Lat, col=allPreds), size = 0.05)  + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="SVCM \nPredictions", low="#003366", high="#99CCFF") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=11, face="bold"),legend.text = element_text(size=9),
        plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))

pdf("figures/p1.pdf", width=6)
p1
dev.off() 

## PLOT RAW DATA ##

p2 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="grey70") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey85") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey85")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey85") +
  geom_point(data=DougScaled, aes(x=Long, y=Lat, col=factor(DougScaled$PRES, levels = c(1,0))), size = 0.05) + 
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

pdf("figures/p2.pdf", width=6)
p2
dev.off() 

## PLOT DNA types ##
p3 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="grey70") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey85") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey85")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey85") +
  geom_point(data=DougPres, aes(x=Long, y=Lat, col=as.factor(DougPres$POP)), size=0.01)+
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=11, face="bold"), legend.text = element_text(size=11),
        legend.key=element_rect(fill = "white", size=0.5),
        plot.margin = margin(0.1,0.1,0.1,0.1, "cm")) + 
  guides(color = guide_legend(override.aes = list(size=3))) +
  scale_color_manual(values=c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20", "#bd0026"),
                     name="DNA types")

pdf("figures/p3.pdf", width=6)
p3
dev.off()

## PLOT Cluster ##
useCluster <- clustered.predictions[[6]]$cluster
col.l <- c("#F3DF6C", "#CEAB07", "#798E87","#C93312", "#CCC591","#C27D38")

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
  ggsci::scale_color_d3(name="Ecotypes")

pdf("figures/p4.pdf", width=6)
p4
dev.off()

# Multipanel plot
 
pdf("figures/Overview22.pdf")
ggpubr::ggarrange(p2, p1, p3, p4,ncol = 2, nrow = 2, hjust=-4)
dev.off()    

### PLOT WITH SCALED COEFFICIENT MATRIX

useCluster <- clustered.predictions.scaled[[6]]$cluster

p4.s <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="grey55") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey68") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey68")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey68") +
  geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=factor(useCluster)), size=0.01) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"), 
        legend.key=element_rect(fill = "white", size=0.5),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  guides(color = guide_legend(override.aes = list(size=3))) +
  scale_color_manual(name="Ecotype", values = col.l)


pdf("figures/ecotypes6-scaled.pdf")
p4.s
dev.off()  

### Same plot with 13 clusters
load("Rdata/clusters2.Rdata")
useCluster <- clusters2[[12]]$cluster
col.13 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')


p12 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="grey55") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey68") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey68")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey68") +
  geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=factor(useCluster)), size=0.01) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"), 
        legend.key=element_rect(fill = "white", size=0.5),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  guides(color = guide_legend(override.aes = list(size=3))) + 
  scale_color_manual(name = "Ecotypes", values = col.13)
p12

pdf("figures/ecotypes12.pdf")
p12
dev.off()
