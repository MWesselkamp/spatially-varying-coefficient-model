require(mgcv)

setwd("~/ScAdditional/PaperSVCM")

source("utils.R")


## Prepare data ##
##################

Douglas <- get_Douglas_data()

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
fit_SVCM <- function(data = DougScaled, grs = c("x", "y")){
  
  # grs either one of  c("x", "y") or  c("Lat", "Long")
  
  Lat = data[grs[1]]
  Long = data[grs[2]]
  
  cores <- detectCores() -1
  c1 <- makeCluster(cores)
  registerDoSNOW(c1)

  fgamSVCtsr <- gam(PRES ~ s(Lat, Long, k=100) + s(Lat, Long, by=TD, k=100) + s(Lat, Long, by=TD2, k=100) + s(Lat, Long, by= PPT_sm, k=100) + s(Lat, Long, by=PPT_sm2, k=100) + s(Lat, Long, by=MWMT, k=100) + s(Lat, Long, by=MWMT2, k=100), method="ML", family=binomial, control=gam.control(maxit=1000), data=DougScaled)

  save(fgamSVCtrs, file="Rdata/fgamSVCtsr.Rdata")
  stopCluster(c1)
  summary(fgamSVCtsr)

  }
# grs = c("Lat", "Long")
load("Rdata/fgamSVCtsr.Rdata")
# grs = c("x", "y") 
load("Rdata/fgamSVCtsr2.Rdata")

fgamSVCtsr <- fgamSVCtrs2
rm(fgamSVCtrs2)

####################################################

#TABLE OF MODEL SUMMARY
anova(fgamSVCtsr)
xtable(summary(fgamSVCtsr)$s.table, digits = 2)

# all model predictions #
allPreds <- predict(fgamSVCtsr, newdata=DougScaled, type="response")
save(allPreds, file="Rdata/allPreds2.Rdata")

# model coefficients at presence observations
predTerms <- predict(fgamSVCtsr, newdata=DougScaled[DougScaled$PRES==1,], type="terms")
save(predTerms, file="Rdata/predTerms2.Rdata")

load("Rdata/predTerms2.Rdata")
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
save(clustered.predictions.scaled, file="Rdata/clusters2scaled.Rdata")

which.max(sapply(clustered.predictions, function(x) x$tot.withinss/x$totss))

### SOME MAPS ###
#################

# LOAD TO REPRODUCE FIGURES
load("Rdata/clusters2.Rdata")
Ecoclusters_all = clusters2[[6]]$cluster
load("Rdata/allPreds2.Rdata")


load("Rdata/states.Rdata")
load("Rdata/states_ca.Rdata")
load("Rdata/states_mex.Rdata")

x <- c("maps", "ggplot2")
lapply(x, library, character.only=TRUE)

countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)

## PLOT OCCURRANCE PROB #

p1 = ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data=DougScaled, aes(x=Long, y=Lat, col=allPreds), size = 0.7)  + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="SVCM \nPredictions", low="#003366", high="#99CCFF", breaks = c(0.25, 0.5, 0.75), labels = c("0.25", "", "0.75")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) 

# SVCM Predictions

p2 = ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data=DougScaled, aes(x=Long, y=Lat, col=factor(PRES, levels = c(1,0))), size = 0.7) + 
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

## PLOT DNA types ##
DNAdata = read.csv("data/Doug-Fir DNA data (Wei et al, 2011).csv")

RColorBrewer::brewer.pal(6, "Set1")
greys = grey.colors(6, start=0.2, end=1.0)
col.6 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")
col.d = c("#ffffb2","#fed976","#feb24c","#f48330","#f03b20", "#bd0026") # same colours as in Paper overview
col.12 = c(col.6, greys)

p3 = ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data = DougScaledPres, aes(x=Long, y=Lat, col=factor(POP)), size=0.8, alpha=1) + 
  geom_point(data=DNAdata, aes(x=Longitude, y=Latitude, color=factor(Genotype), fill = factor(Genotype)), size=1.5, shape=21) +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_color_manual(values = c(col.d, greys), guide = FALSE) +
  scale_fill_manual(values = c(col.d), name = "Genotypes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())  +
  guides(fill = guide_legend(override.aes = list(size=2)))

## PLOT Cluster ##

RColorBrewer::brewer.pal(6, "Dark2")
col.6.2 = c("#F3DF6C", "#CEAB07", "#798E87","#C93312", "#CCC591","#C27D38")
col.l = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF") # same colours as in effect plots.

p4 = ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data = DougScaledPres, aes(x=Long, y=Lat, col=factor(Ecoclusters_all)), size=0.8, alpha=1) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) +
  scale_color_manual(values = c(col.l), name="Ecotypes") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=2)))

 
pdf("figures/Overview.pdf", width = 9, height = 9)
ggpubr::ggarrange(p2, p1, p3, p4,ncol = 2, nrow = 2, hjust=-4, labels= c("a", "b", "c", "d"))
dev.off()    


### PLOT WITH SCALED COEFFICIENT MATRIX
load("Rdata/clusters2scaled.Rdata")
useCluster <- clustered.predictions.scaled[[6]]$cluster
col.l <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF")

p4.s <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") +
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=factor(useCluster)), size=1.0) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_color_manual(values = c(col.l), name = "Ecotypes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())  +
  guides(color = guide_legend(override.aes = list(size=2)))


pdf("figures/ecotypes6-scaled.pdf")
p4.s
dev.off()  

### Same plot with 12 clusters
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

xtable::xtable(confusion_matrix(DougPres$POP, useCluster), digits=0)
