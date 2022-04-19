library(mgcv)

mod = "basic"

if (mod == "ref"){
  load("Rdata/clustersr.Rdata")
  load("Rdata/fgamtrs.Rdata")
  p=""
}else if (mod=="basic"){
  load("Rdata/clusters2.Rdata")
  load("Rdata/fgamSVCtrs2.Rdata")
  p=2
}else if (mod=="non-standardized"){
  load("Rdata/clusters3.Rdata")
  load("Rdata/fgamSVCtrs3.Rdata")
  p=3
}else if (mod=="scaled-location"){
  load("Rdata/clusters4.Rdata")
  load("Rdata/fgamSVCtrs4.Rdata")
  p=4
}else if (mod=="neutral"){
  load("Rdata/clustersN.Rdata")
  load("Rdata/fgamSVCtrsN.Rdata")
  p="N"
}

load("Rdata/DougScaled.Rdata")
load("Rdata/Douglas.Rdata")

source("utils.R")

newTD <- seq(min(Douglas$TD), max(Douglas$TD), len=100)
newTD2 <- newTD^2

newPPTsm <- seq(min(Douglas$PPT_sm), max(Douglas$PPT_sm), len=100)
newPPTsm2 <- newPPTsm^2

newMWMT <- seq(min(Douglas$MWMT), max(Douglas$MWMT), len=100)
newMWMT2 <- newMWMT^2

useCluster <- clusters[[6]]$cluster

if (mod=="basic"){
  
  newTDscaled <- (newTD - mean(Douglas$TD))/sd(Douglas$TD)
  newTD2scaled <- (newTD2 - mean(Douglas$TD^2))/sd(Douglas$TD^2)
  
  newPPTsmscaled <- (newPPTsm - mean(Douglas$PPT_sm))/sd(Douglas$PPT_sm)
  newPPTsm2scaled <- (newPPTsm2 - mean(Douglas$PPT_sm^2))/sd(Douglas$PPT_sm^2)
  
  newMWMTscaled <- (newMWMT - mean(Douglas$MWMT))/sd(Douglas$MWMT)
  newMWMT2scaled <- (newMWMT2 - mean(Douglas$MWMT^2))/sd(Douglas$MWMT^2)
  
} else if (mod=="scaled-location"){
  
    newTDscaled <- (newTD - mean(Douglas$TD))/sd(Douglas$TD)
    newTD2scaled <- (newTD2 - mean(Douglas$TD^2))/sd(Douglas$TD^2)
    
    newPPTsmscaled <- (newPPTsm - mean(Douglas$PPT_sm))/sd(Douglas$PPT_sm)
    newPPTsm2scaled <- (newPPTsm2 - mean(Douglas$PPT_sm^2))/sd(Douglas$PPT_sm^2)
    
    newMWMTscaled <- (newMWMT - mean(Douglas$MWMT))/sd(Douglas$MWMT)
    newMWMT2scaled <- (newMWMT2 - mean(Douglas$MWMT^2))/sd(Douglas$MWMT^2)
    
    DougScaled$x <- (DougScaled$x - min(DougScaled$x))
    DougScaled$y <- (DougScaled$y - min(DougScaled$y))

} else if (mod=="non-standardized"){
  
  newTDscaled <- newTD
  newTD2scaled <- newTD2
  
  newPPTsmscaled <- newPPTsm
  newPPTsm2scaled <- newPPTsm2
  
  newMWMTscaled <- newMWMT
  newMWMT2scaled <- newMWMT2
  
  DougScaled <- Douglas
}


compute_response_curves = function(j){
  
  clustmembers <- as.numeric(attr(useCluster[which(useCluster==j)], "names"))
  cluster <- paste("Cluster", j, sep=" ")
  
  effTD1 <- vector(mode="list", length = length(clustmembers) )

  for (i in 1:length(clustmembers)){
  
    newTDdata <- data.frame(y=DougScaled$y[clustmembers[i]], x=DougScaled$x[clustmembers[i]], TD=newTDscaled, TD2=newTD2scaled, PPT_sm=mean(DougScaled$PPT_sm[clustmembers]), PPT_sm2=mean(DougScaled$PPT_sm2[clustmembers]), MWMT=mean(DougScaled$MWMT[clustmembers]), MWMT2=mean(DougScaled$MWMT2[clustmembers]))
    effTD1[[i]] <- predict.gam(fgamSVCtrs, newdata = newTDdata, type="response")
  
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
    effPPTsm[[i]] <- predict.gam(fgamSVCtrs, newdata = newPPTsmdata, type="response")
    
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
    effMWMT[[i]] <- predict.gam(fgamSVCtrs, newdata = newMWMTdata, type="response")
  
  }

  df3 <- do.call(rbind, effMWMT)
  df3 <- as.data.frame(t(apply(df3, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles
  
  colnames(df3) <- c("first", "second", "median", "third", "fourth")
  df3$ID <- "MTWM"
  df3$Cluster <- cluster
  df3$newseq <- newMWMT

  df.com[[j]] <- rbind(df, df2, df3)
  
}

useCluster <- clustered.predictions[[6]]$cluster
df.com <- vector(mode = "list", length=6)

for (i in 1:6){
 df.com[[i]] = compute_response_curves(i)
}

df = do.call(rbind, df.com)
summary(df)

save(df, file=paste0("Rdata/effectplots", p, ".Rdata"))

#cores = 12
#c1 <- makeCluster(cores)
#registerDoSNOW(c1)
#res = foreach(j = 1:6, .combine=rbind) %dopar% compute_response_curves(j)
#save(res, file="Rdata/effectplot_predictions_df2.Rdata")
#stopCluster(c1)
