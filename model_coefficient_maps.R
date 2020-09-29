require(ggplot2)
require(classInt)

source("utils.R")

Douglas <- get_Douglas_data()

DougPres <- Douglas[Douglas$PRES==1,]

# data set complete - scaled
DougScaled <- Douglas
DougScaled[,c(8:13)] <- scale(Douglas[,c(8:13)])

# data set presences - scaled
DougScaledPres <- DougScaled[DougScaled$PRES==1,]

load("Rdata/fgamSVCtsr2.Rdata")
load("Rdata/predTerms2.Rdata")

countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)

load("Rdata/states.Rdata")
load("Rdata/states_ca.Rdata")
load("Rdata/states_mex.Rdata")

### Coefficient MAPS ##
#######################

coefficient_maps <- function(coefs){
  
  # in order to show variation and capture outliers, use quantile color classes.
  quantiles.1 <- classIntervals(coefs[,2], 25, style="quantile")
  legend_labels.1 <- round(quantiles.1$brks[seq(1,26,length.out = 5)], 0)
  predTerms_cut.1 <- cut(coefs[,2], breaks = quantiles.1$brks, include.lowest = TRUE)

  pTD <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") + 
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.1)), size = 0.05)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="TD ",low="#106607", high="#75F567", breaks= c(0,6.25,12.5,18.75,25),labels=c(as.character(legend_labels.1))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=13, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) 

  ## TD2 
  quantiles.2 <- classIntervals(coefs[,3], 25, style="quantile")
  legend_labels.2 <- round(quantiles.2$brks[seq(1,26,length.out=5)], 0)
  predTerms_cut.2 <- cut(coefs[,3], breaks = quantiles.2$brks, include.lowest = TRUE)

  pTD2 <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") + 
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.2)), size = 0.05)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="TD2 ",low="#106607", high="#75F567", breaks= c(0,6.25,12.5,18.75,25),labels=as.character(legend_labels.2))+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=13, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) 


  ## PPTsm
  quantiles.3 <- classIntervals(coefs[,4], 25, style="quantile")
  legend_labels.3 <- round(quantiles.3$brks[seq(1,26,length.out=5)], 0)
  predTerms_cut.3 <- cut(coefs[,4], breaks = quantiles.3$brks, include.lowest = TRUE)

  pPPT <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") +
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.3)), size = 0.05)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="PPTsm ",low="#003366", high="#99CCFF", breaks= c(0,6.25,12.5,18.75,25),labels=as.character(legend_labels.3))+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=13, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) 


  ## PPTsm2
  quantiles.4 <- classIntervals(coefs[,5], 25, style="quantile")
  legend_labels.4 <- round(quantiles.4$brks[seq(1,26,length.out=5)], 0)
  predTerms_cut.4 <- cut(coefs[,5], breaks = quantiles.4$brks, include.lowest = TRUE)

  pPPT2 <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") + 
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.4)), size = 0.05)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="PPTsm2 ",low="#003366", high="#99CCFF", breaks= c(0,6.25,12.5,18.75,25),labels=as.character(legend_labels.4))+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=13, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) 

## MTWM

  quantiles.5 <- classIntervals(coefs[,6], 25, style="quantile")
  legend_labels.5 <- round(quantiles.5$brks[seq(1,26,length.out = 5)], 0)
  predTerms_cut.5 <- cut(coefs[,6], breaks = quantiles.5$brks, include.lowest = TRUE)

  pMTWM <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") + 
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.5)), size = 0.05)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="MTWM ",low="#F53A14", high="#FDB6A7", breaks= c(0,6.25,12.5,18.75,25),labels=as.character(legend_labels.5))+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=13, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) 

  quantiles.6 <- classIntervals(coefs[,7], 25, style="quantile")
  legend_labels.6 <- round(quantiles.6$brks[seq(1,26,length.out = 5)], 0)
  predTerms_cut.6 <- cut(coefs[,7], breaks = quantiles.6$brks, include.lowest = TRUE)

  pMTWM2 <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") +
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.6)), size = 0.05)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="MTWM2 ",low="#F53A14", high="#FDB6A7", breaks= c(0,6.25,12.5,18.75,25),labels=as.character(legend_labels.6))+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=13, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())

  quantiles.7 <- classIntervals(coefs[,1], 25, style="quantile")
  legend_labels.7 <- round(quantiles.7$brks[seq(1,26,length.out = 5)], 0)
  predTerms_cut.7 <- cut(coefs[,1], breaks = quantiles.7$brks, include.lowest = TRUE)

  pIntercept <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") +
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.7)), size = 0.05)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="Intercept  ",low="#f7e1c2", high="#784a09", breaks= c(0,6.25,12.5,18.75,25),labels=as.character(legend_labels.7))+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=13, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())

  ggpubr::ggarrange(pMTWM, pTD, pPPT, pMTWM2, pTD2, pPPT2, pIntercept, 
                  ncol = 3, nrow = 3, hjust=-4)
}

coefficient_maps(predTerms)


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