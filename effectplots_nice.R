library(foreach)
library(doParallel)
library(mgcv)

newTD <- seq(min(Douglas$TD), max(Douglas$TD), len=100)
newTD2 <- newTD^2
newTDscaled <- (newTD - mean(Douglas$TD))/sd(Douglas$TD)
newTD2scaled <- (newTD2 - mean(Douglas$TD^2))/sd(Douglas$TD^2)

newPPTsm <- seq(min(Douglas$PPT_sm), max(Douglas$PPT_sm), len=100)
newPPTsm2 <- newPPTsm^2
newPPTsmscaled <- (newPPTsm - mean(Douglas$PPT_sm))/sd(Douglas$PPT_sm)
newPPTsm2scaled <- (newPPTsm2 - mean(Douglas$PPT_sm^2))/sd(Douglas$PPT_sm^2)


newMWMT <- seq(min(Douglas$MWMT), max(Douglas$MWMT), len=100)
newMWMT2 <- newMWMT^2
newMWMTscaled <- (newMWMT - mean(Douglas$MWMT))/sd(Douglas$MWMT)
newMWMT2scaled <- (newMWMT2 - mean(Douglas$MWMT^2))/sd(Douglas$MWMT^2)

load("clusters2.Rdata")
# generate new data frame at mean location of clust_members, marginalized to predictor TD
# other predictors are fixed to the mean value of cluster locations.
#useCluster <- clustered.predictions[[6]]
#useCluster <- useCluster$cluster
useCluster <- clusters2[[6]]$cluster
df.com <- vector(mode = "list", length=6)

for (j in 1:6){
  
  clustmembers <- as.numeric(attr(useCluster[which(useCluster==j)], "names"))
  cluster <- paste("Cluster", j, sep=" ")
  
  #c1 <- detectCores()-1
  #registerDoParallel(c1)
  effTD1 <- vector(mode="list", length = length(clustmembers) )

  for (i in 1:length(clustmembers)){
  
    newTDdata <- data.frame(y=DougScaled$y[clustmembers[i]], x=DougScaled$x[clustmembers[i]], TD=newTDscaled, TD2=newTD2scaled, PPT_sm=mean(DougScaled$PPT_sm[clustmembers]), PPT_sm2=mean(DougScaled$PPT_sm2[clustmembers]), MWMT=mean(DougScaled$MWMT[clustmembers]), MWMT2=mean(DougScaled$MWMT2[clustmembers]))
    effTD1[[i]] <- predict.gam(fgamSVCtrs2, newdata = newTDdata, type="response")
  
  }

  #doParallel::stopImplicitCluster()
  df <- do.call(rbind, effTD1)

  df <- as.data.frame(t(apply(df, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles

  colnames(df) <- c("first", "second", "median", "third", "fourth")
  df$ID <- "TD"
  df$Cluster <- cluster
  df$newseq <- newTD

  effPPTsm <- vector(mode="list", length = length(clustmembers) )
  
  for (i in 1:length(clustmembers)){
    
    newPPTsmdata <- data.frame(y=DougScaled$y[clustmembers[i]], x=DougScaled$x[clustmembers[i]], TD=mean(DougScaled$TD[clustmembers]), TD2=mean(DougScaled$TD2[clustmembers]), PPT_sm=newPPTsmscaled, PPT_sm2=newPPTsm2scaled, MWMT=mean(DougScaled$MWMT[clustmembers]), MWMT2=mean(DougScaled$MWMT2[clustmembers]))
    effPPTsm[[i]] <- predict.gam(fgamSVCtrs2, newdata = newPPTsmdata, type="response")
    
  }
  
  #doParallel::stopImplicitCluster()
  df2 <- do.call(rbind, effPPTsm)
  
  df2 <- as.data.frame(t(apply(df2, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles
  
  colnames(df2) <- c("first", "second", "median", "third", "fourth")
  df2$ID <- "PPT_sm"
  df2$Cluster <- cluster
  df2$newseq <- newPPTsm

  #c1 <- detectCores()-1
  #registerDoParallel(c1)
  effMWMT <- vector(mode="list", length = length(clustmembers) )

  for (i in 1:length(clustmembers)){
    newMWMTdata <- data.frame(y=DougScaled$y[clustmembers[i]], x=DougScaled$x[clustmembers[i]], TD=mean(DougScaled$TD[clustmembers]), TD2=mean(DougScaled$TD2[clustmembers]), PPT_sm=mean(DougScaled$PPT_sm[clustmembers]), PPT_sm2=mean(DougScaled$PPT_sm2[clustmembers]), MWMT=newMWMTscaled, MWMT2=newMWMT2scaled)
    effMWMT[[i]] <- predict.gam(fgamSVCtrs2, newdata = newMWMTdata, type="response")
  
  }

  #doParallel::stopImplicitCluster()
  df3 <- do.call(rbind, effMWMT)
  df3 <- as.data.frame(t(apply(df3, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles
  
  colnames(df3) <- c("first", "second", "median", "third", "fourth")
  df3$ID <- "MWMT"
  df3$Cluster <- cluster
  df3$newseq <- newMWMT

  df.com[[j]] <- rbind(df, df2, df3)

}

df.all <- do.call(rbind, df.com)
save(df.all, file="Rdata/effectplot_predictions_df2.Rdata")

load("Rdata/effectplot_predictions_df2.Rdata")

df.all$ID_f <- factor(df.all$ID, levels = c("MWMT","TD","PPT_sm"), labels = c("MTWM (°C)", "TD (°C)", "PPT_sm (mm)"))
df.all$Cluster <- gsub(df.all$Cluster, pattern = "Cluster", replacement = "Ecotype")

pdf("figures/Effects.pdf", height = 12, width = 7)

ggplot(df.all) + ylim(0,1) + geom_path(aes(x=newseq, y=median, col=as.factor(Cluster))) + labs(y="Occurrence Probability")+
  geom_ribbon(aes(ymin = third, ymax = second, x=newseq), fill="grey50" , alpha=0.3) + 
  theme(aspect.ratio = 1, axis.title.x = element_blank(), strip.text.x = element_text(face="bold"),
        strip.text.y = element_text(face="bold"),panel.background = element_blank(), 
        legend.key = element_blank(), legend.position = "none") + 
  facet_grid(Cluster~ID_f, scales="free_x") + 
  ggsci::scale_color_d3(name="Cluster")

dev.off()
