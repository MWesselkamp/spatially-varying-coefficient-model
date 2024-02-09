require(ggplot2)
require(classInt)
require(mgcv)

source("utils.R")
source("plots.R")

Douglas <- get_Douglas_data()

DougPres <- Douglas[Douglas$PRES==1,]

# data set complete - scaled
DougScaled <- Douglas
DougScaled[,c(10:19)] <- scale(Douglas[,c(10:19)])

# data set presences - scaled
DougScaledPres <- DougScaled[DougScaled$PRES==1,]

load("Rdata/fgamSVC.Rdata")
load("Rdata/fgamSVC_predTerms.Rdata")
load("Rdata/fgamSVC_clusters.Rdata")

countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)

load("Rdata/states.Rdata")
load("Rdata/states_ca.Rdata")
load("Rdata/states_mex.Rdata")

### Coefficient MAPS ##
#######################

pdf(file = "plots/coeffient_maps.pdf", width = 10, height = 10)
coefficient_maps(predTerms, dat=DougPres)
dev.off()

## Remove ouliers before mapping ## 
###################################
predTerms.e <- cbind(DougScaledPres[,c(4:5)], predTerms)
predTerms.r <- apply(as.data.frame(predTerms.e[,4]), 2, quantile, c(0.99,0.01))
preds.r.in <- as.data.frame(predTerms.e[which((predTerms.e[,4] < predTerms.r[1]) & (predTerms.e[,4] > predTerms.r[2])),])
preds.r.out <- as.data.frame(predTerms.e[which(!((predTerms.e[,4] < predTerms.r[1]) & (predTerms.e[,4] > predTerms.r[2]))),])

pTDr <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="grey55") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey68") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey68")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey68") +
  geom_point(data=preds.r.in, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.in[,4]))), size = 0.05)  + 
  geom_point(data=preds.r.out, aes(x=Long, y=Lat, col=preds.r.out[,4]), size = 0.05, colour="white") +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="TD2 ",low="#106607", high="#75F567")+ # values = cols)+ #,  breaks = quantiles$brks) + #low="#003366", high="#99CCFF",
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"),legend.text = element_text(size=8),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

predTerms.r2 <- apply(as.data.frame(predTerms.e[,5]), 2, quantile, c(0.99,0.01))
preds.r.in2 <- as.data.frame(predTerms.e[which((predTerms.e[,5] < predTerms.r2[1]) & (predTerms.e[,5] > predTerms.r2[2])),])
preds.r.out2 <- as.data.frame(predTerms.e[which(!((predTerms.e[,5] < predTerms.r2[1]) & (predTerms.e[,5] > predTerms.r2[2]))),])

pTDr2 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="grey55") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey68") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey68")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey68") +
  geom_point(data=preds.r.in2, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.in2[,5]))), size = 0.05)  + 
  geom_point(data=preds.r.out2, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.out2[,4]))), size = 0.05, colour="white") +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="TD2 ",low="#106607", high="#75F567")+ # values = cols)+ #,  breaks = quantiles$brks) + #low="#003366", high="#99CCFF",
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"),legend.text = element_text(size=8),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))


predTerms.r3 <- apply(as.data.frame(predTerms.e[,6]), 2, quantile, c(0.99,0.01))
preds.r.in3 <- as.data.frame(predTerms.e[which((predTerms.e[,6] < predTerms.r3[1]) & (predTerms.e[,6] > predTerms.r3[2])),])
preds.r.out3 <- as.data.frame(predTerms.e[which(!((predTerms.e[,6] < predTerms.r3[1]) & (predTerms.e[,6] > predTerms.r3[2]))),])

pPPTr <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="grey55") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey68") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey68")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey68") +
  geom_point(data=preds.r.in3, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.in3[,6]))), size = 0.05)  + 
  geom_point(data=preds.r.out3, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.out3[,6]))), size = 0.05, colour="white") +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="PPTsm ",low="#003366", high="#99CCFF")+ # values = cols)+ #,  breaks = quantiles$brks) + #low="#003366", high="#99CCFF",
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"),legend.text = element_text(size=8),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

predTerms.r4 <- apply(as.data.frame(predTerms.e[,7]), 2, quantile, c(0.99,0.01))
preds.r.in4 <- as.data.frame(predTerms.e[which((predTerms.e[,7] < predTerms.r4[1]) & (predTerms.e[,7] > predTerms.r4[2])),])
preds.r.out4 <- as.data.frame(predTerms.e[which(!((predTerms.e[,7] < predTerms.r4[1]) & (predTerms.e[,7] > predTerms.r4[2]))),])

pPPTr2 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="grey55") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey68") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey68")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey68") +
  geom_point(data=preds.r.in4, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.in4[,7]))), size = 0.05)  + 
  geom_point(data=preds.r.out4, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.out4[,7]))), size = 0.05, colour="white") +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="PPTsm2 ",low="#003366", high="#99CCFF")+ # values = cols)+ #,  breaks = quantiles$brks) + #low="#003366", high="#99CCFF",
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"),legend.text = element_text(size=8),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

predTerms.r5 <- apply(as.data.frame(predTerms.e[,8]), 2, quantile, c(0.99,0.01))
preds.r.in5 <- as.data.frame(predTerms.e[which((predTerms.e[,8] < predTerms.r5[1]) & (predTerms.e[,8] > predTerms.r5[2])),])
preds.r.out5 <- as.data.frame(predTerms.e[which(!((predTerms.e[,8] < predTerms.r5[1]) & (predTerms.e[,8] > predTerms.r5[2]))),])

pMWMTr <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="grey55") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey68") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey68")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey68") +
  geom_point(data=preds.r.in5, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.in5[,8]))), size = 0.05)  + 
  geom_point(data=preds.r.out5, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.out5[,8]))), size = 0.05, colour="white") +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="MWMT ",low="#F53A14", high="#FDB6A7")+ # values = cols)+ #,  breaks = quantiles$brks) + #low="#003366", high="#99CCFF",
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"),legend.text = element_text(size=8),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

predTerms.r6 <- apply(as.data.frame(predTerms.e[,9]), 2, quantile, c(0.99,0.01))
preds.r.in6 <- as.data.frame(predTerms.e[which((predTerms.e[,9] < predTerms.r6[1]) & (predTerms.e[,9] > predTerms.r6[2])),])
preds.r.out6 <- as.data.frame(predTerms.e[which(!((predTerms.e[,9] < predTerms.r6[1]) & (predTerms.e[,9] > predTerms.r6[2]))),])

pMWMTr2 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="grey55") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey68") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey68")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey68") +
  geom_point(data=preds.r.in6, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.in6[,9]))), size = 0.05)  + 
  geom_point(data=preds.r.out6, aes(x=Long, y=Lat, col=as.numeric(scale(preds.r.out6[,9]))), size = 0.05, colour="white") +
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_colour_gradient(name="MWMT ",low="#F53A14", high="#FDB6A7")+ # values = cols)+ #,  breaks = quantiles$brks) + #low="#003366", high="#99CCFF",
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=8, face="bold"),legend.text = element_text(size=8),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))


pdf("figures/coefsInSpace_rescaled_noOuties.pdf")
ggpubr::ggarrange(pTDr, pPPTr, pMWMTr, pTDr2, pPPTr2, pMWMTr2, ncol = 3, nrow = 2)
dev.off()