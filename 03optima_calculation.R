load("Rdata/predTerms2.Rdata")
load("Rdata/clusters2.Rdata")
load("Rdata/DougScaled.Rdata")

useCluster <- clusters2[[6]]$cluster
DougScaledPres <- DougScaled[which(DougScaled$PRES==1),]

tail(predTerms)
head(predTerms)


### CALCULATE OPTIMA ###
########################

optima <- function(preds){
  TDopt = vector(mode = "list", length=nrow(preds))
  PPTsmopt = vector(mode = "list", length=nrow(preds))
  MWMTopt = vector(mode = "list", length=nrow(preds))
  
  for (i in 1:nrow(preds)){
    
    pred <- preds[i,]
    
    if (pred[3] < 0){
      TDmax <- -(pred[2]/pred[3]*2)
      TDopt[i] <- TDmax
    } else {
      TDopt[i] <- NA
    }
    
    if (pred[5] < 0){
      PPTsmmax <- -(pred[4]/pred[5]*2)
      PPTsmopt[i] <- PPTsmmax
    } else {
      PPTsmopt[i] <-  NA
    }
    
    if (pred[7] < 0){
      MWMTmax <- -(pred[6]/pred[7]*2)
      MWMTopt[i] <- MWMTmax
    } else {
      MWMTopt[i] <- NA
    }
  }
  optima <- cbind(TDopt, PPTsmopt, MWMTopt)
  return(optima)
}
save(optima, file="functions/optima.Rdata")
optimas <- as.data.frame(optima(predTerms))

plot(optimas$TDopt, optimas$PPTsmopt)

##rescale optimas

optims_r <- sapply(optimas[,1], function(x) {x*sd(Douglas$TD) + mean(Douglas$TD)})
optims_r2 <- sapply(optimas[,2], function(x) {x*sd(Douglas$PPT_sm) + mean(Douglas$PPT_sm)})
optims_r3 <- sapply(optimas[,3], function(x) {x*sd(Douglas$MWMT) + mean(Douglas$MWMT)})

optims_r <- cbind(optims_r, optims_r2, optims_r3)

summary(as.data.frame(optims_r))

## Check optima in environmental space
clusters <- as.data.frame(cbind(optims_r, useCluster))
head(clusters)

ggplot(clusters) + geom_point(aes(x=optims_r, y=optims_r2, color=as.factor(useCluster))) + xlim(-800,800) + ylim(-400,400)
ggplot(clusters) + geom_point(aes(x=optims_r2, y=optims_r3, color=as.factor(useCluster))) + xlim(-800,800) + ylim(-400,400)

### Minimum and Maximum Predicitons ###
#######################################

allmembers <- as.numeric(attr(useCluster, "names"))


## TD ##
########

TDmin <- min(Douglas$TD) # pull minimum value for TD from data set
# design decision: use complete environemnt (douglas$TD) or habitat environment (DougPres$TD)
TDmin <- which((Douglas$TD)==TDmin) # sample index for observation with this value to receive the according scaled value
pminTD <- vector(mode = "list", length = length(allmembers))
for (i in 1:length(allmembers)){
  newTDdata <- data.frame(y=DougScaledPres$y[i], x=DougScaledPres$x[i], TD=DougScaled$TD[TDmin], TD2=DougScaled$TD2[TDmin], PPT_sm=0, PPT_sm2=0, MWMT=0, MWMT2=0)
  # new data set for each observation with minimum values from the scaled data set
  pminTD[i] <- predict(fgamSVCtrs2, newdata = newTDdata, type="response")
  # minimum prediction for each observation
}

TDmax <- max(Douglas$TD)
TDmax <- which((Douglas$TD)==TDmax)
pmaxTD <- vector(mode = "list", length = length(allmembers))
for (i in 1:length(allmembers)){
  newTDdata <- data.frame(y=DougScaledPres$y[i], x=DougScaledPres$x[i], TD=DougScaled$TD[TDmax], TD2=DougScaled$TD2[TDmax], PPT_sm=0, PPT_sm2=0, MWMT=0, MWMT2=0)
  pmaxTD[i] <- predict(fgamSVCtrs2, newdata = newTDdata, type="response")
}


## PPT_sm ##
############

PPTsmmin <- min((na.exclude(Douglas)$PPT_sm))
PPTsmmin <- sample(which((Douglas$PPT_sm)==PPTsmmin),1) # 1 Nur ein PPTsmmin verwenden
pminPPTsm <- vector(mode = "list", length = length(allmembers))
for (i in 1:length(allmembers)){
  newPPTsmdata <- data.frame(y=DougScaledPres$y[i], x=DougScaledPres$x[i], PPT_sm=DougScaled$PPT_sm[PPTsmmin], PPT_sm2=DougScaled$PPT_sm2[PPTsmmin], TD=0, TD2=0, MWMT=0, MWMT2=0)
  pminPPTsm[i] <- predict(fgamSVCtrs2, newdata = newPPTsmdata, type="response")
}

PPTsmmax <- max(na.exclude(Douglas)$PPT_sm)
PPTsmmax <- which((Douglas$PPT_sm)==PPTsmmax)
pmaxPPTsm <- numeric()
for (i in 1:length(allmembers)){
  newPPTsmdata <- data.frame(y=DougScaledPres$y[i], x=DougScaledPres$x[i], PPT_sm=DougScaled$PPT_sm[PPTsmmax], PPT_sm2=DougScaled$PPT_sm2[PPTsmmax], TD=0, TD2=0, MWMT=0, MWMT2=0)
  pmaxPPTsm[i] <- predict(fgamSVCtrs2, newdata = newPPTsmdata, type="response")
}


## MWMT ####
############


MWMTmin <- min(Douglas$MWMT)
MWMTmin <- which(Douglas$MWMT==MWMTmin)
pminMWMT <- vector(mode = "list", length = length(allmembers))
for (i in 1:length(allmembers)){
  newMWMTdata <- data.frame(y=DougScaledPres$y[i], x=DougScaledPres$x[i], MWMT=DougScaled$MWMT[MWMTmin], MWMT2=DougScaled$MWMT2[MWMTmin], PPT_sm=0, PPT_sm2=0, TD=0, TD2=0)
  pminMWMT[i] <- predict(fgamSVCtrs2, newdata = newMWMTdata, type="response")
}

MWMTmax <- max(Douglas$MWMT)
MWMTmax <- sample(which(Douglas$MWMT==MWMTmax),1)
pmaxMWMT <- numeric()
for (i in 1:length(allmembers)){
  newMWMTdata <- data.frame(y=DougScaledPres$y[i], x=DougScaledPres$x[i], MWMT=DougScaled$MWMT[MWMTmax], MWMT2=DougScaled$MWMT2[MWMTmax], PPT_sm=0, PPT_sm2=0, TD=0, TD2=0)
  pmaxMWMT[i] <- predict(fgamSVCtrs2, newdata = newMWMTdata, type="response")
}


pTD <- cbind(pminTD, pmaxTD)
pPPTsm <- cbind(pminPPTsm, pmaxPPTsm)
pMWMT <- cbind(pminMWMT, pmaxMWMT)
pminmax <- cbind(pTD, pPPTsm)
pminmax <- cbind(pminmax, pMWMT)
pminmax <- as.data.frame(pminmax)

pminmax <- as.data.frame(t(apply(pminmax, 1, unlist)))
save(pminmax, file="Rdata/pminmax2.Rdata")

TDminmax <- vector(mode = "list", length = length(nrow(pminmax)))
for (i in 1:nrow(pminmax)){
  if (pminmax$pminTD[i] > pminmax$pmaxTD[i]){
    TDminmax[i] <- min(Douglas$TD)
  }else{
    TDminmax[i] <- max(Douglas$TD)
  }
}
PPTminmax <- vector(mode = "list", length = length(nrow(pminmax)))
for (i in 1:nrow(pminmax)){
  if (pminmax$pminPPTsm[i] > pminmax$pmaxPPTsm[i]){
    PPTminmax[i] <- min(Douglas$PPT_sm)
  }else{
    PPTminmax[i] <- max(Douglas$PPT_sm)
  }
}
MWMTminmax <- vector(mode = "list", length = length(nrow(pminmax)))
for (i in 1:nrow(pminmax)){
  if (pminmax$pminMWMT[i] > pminmax$pmaxMWMT[i]){
    MWMTminmax[i] <- min(Douglas$MWMT)
  } else  {
    MWMTminmax[i] <- max(Douglas$MWMT)
  }
}

MinMax <- as.data.frame(cbind(as.numeric(TDminmax), as.numeric(PPTminmax), as.numeric(MWMTminmax)))
save(MinMax, file="Rdata/MinMax2.Rdata")

### REPLACE MISSING OPTIMA BY MINIMUM MAXIMUM PROBS ###
#######################################################

optima_2 <- function(optims, minmax, min, max){
  for (i in 1:length(optims)){
    if (is.na(optims[i]) | optims[i] < min | optims[i] > max){
      optims[i] <- minmax[i]
    } 
  }
  return(optims)
}

TDopts <- optima_2(optims = optims_r[,1], minmax=MinMax[,1], min = min(Douglas$TD), max= max(Douglas$TD))
PPTopts <- optima_2(optims = optims_r[,2], minmax=MinMax[,2], min = min(Douglas$PPT_sm), max= max(Douglas$PPT_sm))
MWMTopts <- optima_2(optims = optims_r[,3], minmax=MinMax[,3], min = min(Douglas$MWMT), max= max(Douglas$MWMT))

Opts <- data.frame("TDopts"= TDopts, "PPTsmopts"=PPTopts, "MWMTopts" = MWMTopts)
summary(Opts)
save(Opts, file="Rdata/Opts2.Rdata")

### CALCULATE CLOSEST POINT IN ENVIRONMENTAL SPACE ###
######################################################
require("Biobase")
names(Opts) <- c("TD", "PPT_sm", "MWMT")
gri <- data.frame("TD" = Douglas$TD, "PPT_sm"= Douglas$PPT_sm, "MWMT"=Douglas$MWMT)
gri <- gri[sample(nrow(gri), 10000),]
 
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
euc.dist2 <- function(x1, x2) dist(rbind(x1,x2), method = "euclidean")

microbenchmark(euc.dist(gri[1,], Opts[1,]), euc.dist2(gri[1,], Opts[1,]), times=100)
# use the dist function

start.time <- Sys.time()

Proj_Opts <- matrix(NA, ncol=3, nrow=nrow(Opts))

for (j in 1:nrow(Opts)){
  opts <- as.numeric(scale(Opts)[j,])
  dists <- apply(scale(gri), 1, function(x) euc.dist2(opts, x) )
  Proj_Opts[j,] <- as.numeric(gri[which.min(dists),])
}

save(Proj_Opts, file="Rdata/Proj_Opts2.Rdata")
end.time <- Sys.time()

start.time - end.time

summary(Proj_Opts)
