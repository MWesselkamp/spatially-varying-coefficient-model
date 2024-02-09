#library(doSNOW)
library(mgcv)

load("Rdata/fgamSVC.Rdata")
load("Rdata/fgamSVC_predTerms.Rdata")
load("Rdata/fgamSVC_clusters.Rdata")
load("Rdata/DougScaled.Rdata")
load("Rdata/Douglas.Rdata")

source("utils.R")

newTD <- seq(min(Douglas$TD), max(Douglas$TD), len=100)
newTD2 <- newTD^2

newPPTsm <- seq(min(Douglas$PPT_sm), max(Douglas$PPT_sm), len=100)
newPPTsm2 <- newPPTsm^2

newMWMT <- seq(min(Douglas$MWMT), max(Douglas$MWMT), len=100)
newMWMT2 <- newMWMT^2
  
newTDscaled <- (newTD - mean(Douglas$TD))/sd(Douglas$TD)
newTD2scaled <- (newTD2 - mean(Douglas$TD^2))/sd(Douglas$TD^2)
    
newPPTsmscaled <- (newPPTsm - mean(Douglas$PPT_sm))/sd(Douglas$PPT_sm)
newPPTsm2scaled <- (newPPTsm2 - mean(Douglas$PPT_sm^2))/sd(Douglas$PPT_sm^2)
    
newMWMTscaled <- (newMWMT - mean(Douglas$MWMT))/sd(Douglas$MWMT)
newMWMT2scaled <- (newMWMT2 - mean(Douglas$MWMT^2))/sd(Douglas$MWMT^2)
    
DougScaled$x <- (DougScaled$x - min(DougScaled$x))/1000
DougScaled$y <- (DougScaled$y - min(DougScaled$y))/1000

compute_response_curves = function(clus){
  
  clustmembers <- as.numeric(attr(useCluster[which(useCluster==clus)], "names"))
  cluster <- paste("Cluster", clus, sep=" ")
  
  effTD1 <- vector(mode="list", length = length(clustmembers) )

  for (i in 1:length(clustmembers)){
  
    newTDdata <- data.frame(y=DougScaled$y[clustmembers[i]], x=DougScaled$x[clustmembers[i]], TD=newTDscaled, TD2=newTD2scaled, PPT_sm=mean(DougScaled$PPT_sm[clustmembers]), PPT_sm2=mean(DougScaled$PPT_sm2[clustmembers]), MWMT=mean(DougScaled$MWMT[clustmembers]), MWMT2=mean(DougScaled$MWMT2[clustmembers]))
    effTD1[[i]] <- predict.gam(fgamSVC, newdata = newTDdata, type="response")
  
  }

  df <- do.call(rbind, effTD1)

  df <- as.data.frame(t(apply(df, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles

  colnames(df) <- c("first", "second", "median", "third", "fourth")
  df$ID <- "TD"
  df$Cluster <- cluster
  df$newseq <- newTD

  effPPTsm <- vector(mode="list", length = length(clustmembers) )
  
  for (i in 1:length(clustmembers)){
    
    newPPTsmdata <- data.frame(y=DougScaled$y[clustmembers[i]], x=DougScaled$x[clustmembers[i]], TD=mean(DougScaled$TD[clustmembers]), TD2=mean(DougScaled$TD2[clustmembers]), PPT_sm=newPPTsmscaled, PPT_sm2=newPPTsm2scaled, MWMT=mean(DougScaled$MWMT[clustmembers]), MWMT2=mean(DougScaled$MWMT2[clustmembers]))
    effPPTsm[[i]] <- predict.gam(fgamSVC, newdata = newPPTsmdata, type="response")
    
  }
  
  df2 <- do.call(rbind, effPPTsm)
  
  df2 <- as.data.frame(t(apply(df2, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles
  
  colnames(df2) <- c("first", "second", "median", "third", "fourth")
  df2$ID <- "PPT_sm"
  df2$Cluster <- cluster
  df2$newseq <- newPPTsm

  effMWMT <- vector(mode="list", length = length(clustmembers) )

  for (i in 1:length(clustmembers)){
    newMWMTdata <- data.frame(y=DougScaled$y[clustmembers[i]], x=DougScaled$x[clustmembers[i]], TD=mean(DougScaled$TD[clustmembers]), TD2=mean(DougScaled$TD2[clustmembers]), PPT_sm=mean(DougScaled$PPT_sm[clustmembers]), PPT_sm2=mean(DougScaled$PPT_sm2[clustmembers]), MWMT=newMWMTscaled, MWMT2=newMWMT2scaled)
    effMWMT[[i]] <- predict.gam(fgamSVC, newdata = newMWMTdata, type="response")
  
  }

  df3 <- do.call(rbind, effMWMT)
  df3 <- as.data.frame(t(apply(df3, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles
  
  colnames(df3) <- c("first", "second", "median", "third", "fourth")
  df3$ID <- "MTWM"
  df3$Cluster <- cluster
  df3$newseq <- newMWMT

  df_responses[[clus]] <- rbind(df, df2, df3)
  
}

nclusters = 10
useCluster <- clusters[[nclusters]]$clustering

df_responses <- vector(mode = "list", length=nclusters)

for (j in 1:nclusters){
  df_responses[[j]] = compute_response_curves(clus = j)
}

df = do.call(rbind, df_responses)
summary(df)

save(df, file=paste0("Rdata/effectplots", nclusters, ".Rdata"))

#cores = 12
#c1 <- makeCluster(cores)
#registerDoSNOW(c1)
#res = foreach(j = 1:6, .combine=rbind) %dopar% compute_response_curves(j)
#save(res, file="Rdata/effectplot_predictions_df2.Rdata")
#stopCluster(c1)


source("plots.R")
pdf(file=paste0("plots/SVCgam_effectplots", nclusters, ".pdf"), width = 8, height = 12)
effectplots(df, nclusters = nclusters)
dev.off() 

