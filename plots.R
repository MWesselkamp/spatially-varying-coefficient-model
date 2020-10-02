#================#
# Plot functions #
#================#
library(ggplot2)
library(maps)


countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)

load("Rdata/states.Rdata")
load("Rdata/states_ca.Rdata")
load("Rdata/states_mex.Rdata")

occurance_predictions <- function(dat = DougScaled, predictions = allPreds){
  
  p1 = ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=dat, aes(x=Long, y=Lat, col=predictions), size = 0.7)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="SVCM \nPredictions", low="#003366", high="#99CCFF", breaks = c(0.25, 0.5, 0.75), labels = c("0.25", "", "0.75")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour="grey45", fill=NA, size=1), 
          panel.background = element_blank(), axis.title = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
          legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) 
  
  return(p1)
}

occurances <- function(dat= DougScaled){
  
  p2 = ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=dat, aes(x=Long, y=Lat, col=factor(PRES, levels = c(1,0))), size = 0.7) + 
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
  
  return(p2)
}

DNA_types = function(dat = DougScaledPres, dna_dat=DNAdata){
  
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
    geom_point(data = dat, aes(x=Long, y=Lat, col=factor(POP)), size=0.8, alpha=1) + 
    geom_point(data=dna_dat, aes(x=Longitude, y=Latitude, color=factor(Genotype), fill = factor(Genotype)), size=1.5, shape=21) +
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_color_manual(values = c(col.d, greys), guide = FALSE) +
    scale_fill_manual(values = c(col.d), name = "Genotypes")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour="grey45", fill=NA, size=1), 
          panel.background = element_blank(), axis.title = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
          legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())  +
    guides(fill = guide_legend(override.aes = list(size=2)))
  
  return(p3)
  
}

Ecotypes = function(dat=DougScaledPres, cluster=useCluster, n_colors=6){
  
  RColorBrewer::brewer.pal(6, "Dark2")
  #col.6.2 = c("#F3DF6C", "#CEAB07", "#798E87","#C93312", "#CCC591","#C27D38")
  if (n_colors==6){
    cols = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF")
    }else{
    cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')
    }

  p4 = ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
    geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data = DougScaledPres, aes(x=Long, y=Lat, col=factor(cluster)), size=0.8, alpha=1) + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) +
    scale_color_manual(values = cols, name="Ecotypes") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour="grey45", fill=NA, size=1), 
          panel.background = element_blank(), axis.title = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
          legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) +
    guides(color = guide_legend(override.aes = list(size=2)))
 
  return(p4) 
}

coefficient_maps <- function(coefs, round_to = 0){
  
  # in order to show variation and capture outliers, use quantile color classes.
  quantiles.1 <- classIntervals(coefs[,2], 25, style="quantile")
  legend_labels.1 <- round(quantiles.1$brks[seq(1,26,length.out = 5)], round_to)
  predTerms_cut.1 <- cut(coefs[,2], breaks = quantiles.1$brks, include.lowest = TRUE)
  
  pTD <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") + 
    #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
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
  legend_labels.2 <- round(quantiles.2$brks[seq(1,26,length.out=5)], round_to)
  predTerms_cut.2 <- cut(coefs[,3], breaks = quantiles.2$brks, include.lowest = TRUE)
  
  pTD2 <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") + 
    #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
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
  legend_labels.3 <- round(quantiles.3$brks[seq(1,26,length.out=5)], round_to)
  predTerms_cut.3 <- cut(coefs[,4], breaks = quantiles.3$brks, include.lowest = TRUE)
  
  pPPT <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") +
    #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
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
  legend_labels.4 <- round(quantiles.4$brks[seq(1,26,length.out=5)], round_to)
  predTerms_cut.4 <- cut(coefs[,5], breaks = quantiles.4$brks, include.lowest = TRUE)
  
  pPPT2 <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") + 
    #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
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
  legend_labels.5 <- round(quantiles.5$brks[seq(1,26,length.out = 5)], round_to)
  predTerms_cut.5 <- cut(coefs[,6], breaks = quantiles.5$brks, include.lowest = TRUE)
  
  pMTWM <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") + 
    #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.5)), size = 0.05)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="MTWM ",low="#F53A14", high="#FDB6A7", breaks= c(0,6.25,12.5,18.75,25),labels=as.character(legend_labels.5))+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour="grey45", fill=NA, size=1), 
          panel.background = element_blank(), axis.title = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
          legend.title = element_text(size=13, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank()) 
  
  quantiles.6 <- classIntervals(coefs[,7], 25, style="quantile")
  legend_labels.6 <- round(quantiles.6$brks[seq(1,26,length.out = 5)], round_to)
  predTerms_cut.6 <- cut(coefs[,7], breaks = quantiles.6$brks, include.lowest = TRUE)
  
  pMTWM2 <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") +
    #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
    geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=as.numeric(predTerms_cut.6)), size = 0.05)  + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
    scale_colour_gradient(name="MTWM2 ",low="#F53A14", high="#FDB6A7", breaks= c(0,6.25,12.5,18.75,25),labels=as.character(legend_labels.6))+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour="grey45", fill=NA, size=1), 
          panel.background = element_blank(), axis.title = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
          legend.title = element_text(size=13, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())
  
  quantiles.7 <- classIntervals(coefs[,1], 25, style="quantile")
  legend_labels.7 <- round(quantiles.7$brks[seq(1,26,length.out = 5)], round_to)
  predTerms_cut.7 <- cut(coefs[,1], breaks = quantiles.7$brks, include.lowest = TRUE)
  
  pIntercept <- ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), colour="grey47", fill="black") +
    #geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
    #geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
    #geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
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


effectplots = function(df,var, save = TRUE){
  
  load(paste0("Rdata/", var, ".Rdata"))
  
  df$ID_f <- factor(df$ID, levels = c("MTWM","TD","PPT_sm"), labels = c("MTWM (°C)", "TD (°C)", "PPT_sm (mm)"))
  df$Cluster <- gsub(df$Cluster, pattern = "Cluster", replacement = "Ecotype")
  
  #if (save) pdf(paste0("figures/Effects", var, ".pdf"), height = 12, width = 7)
    
  
  p = ggplot(df) + ylim(0,1) + geom_path(aes(x=newseq, y=median, col=as.factor(Cluster))) + labs(y="Occurrence Probability")+
    geom_ribbon(aes(ymin = third, ymax = second, x=newseq), fill="grey50" , alpha=0.3) + 
    theme(aspect.ratio = 1, axis.title.x = element_blank(), strip.text.x = element_text(face="bold"),
          strip.text.y = element_text(face="bold"),panel.background = element_blank(), 
          legend.key = element_blank(), legend.position = "none") + 
    facet_grid(Cluster~ID_f, scales="free_x") + 
    ggsci::scale_color_d3(name="Cluster")
  
  #if (save) dev.off()
    
  return(p)
  
}
