 ### Minimum and Maximum Predicitons ###
#######################################
source("utils.R")

load("Rdata/clusters2.Rdata")
useCluster = clusters2[[6]]$cluster
allmembers <- as.numeric(attr(useCluster, "names"))

Douglas = get_Douglas_data()

load("Rdata/DougScaled.Rdata")
DougScaledPres = DougScaled[DougScaled$PRES==1,]

## TD ##
########

TDmin <- min(Douglas$TD) # pull minimum value for TD from data set
# design decision: use complete environemnt (douglas$TD) or habitat environment (DougPres$TD)
TDmin <- which((Douglas$TD)==TDmin) # sample index for observation with this value to receive the according scaled value
pminTD <- vector(mode = "list", length = length(allmembers))
for (i in 1:length(allmembers)){
  newTDdata <- data.frame(Lat=DougScaledPres$Lat[i], Long=DougScaledPres$Long[i], TD=DougScaled$TD[TDmin,], TD2=DougScaled$TD2[TDmin,], PPT_sm=0, PPT_sm2=0, MWMT=0, MWMT2=0)
  # new data set for each observation with minimum values from the scaled data set
  pminTD[i] <- predict(fgamSVCtsr, newdata = newTDdata, type="response")
 # minimum prediction for each observation
}

TDmax <- max(Douglas$TD)
TDmax <- which((Douglas$TD)==TDmax)
pmaxTD <- vector(mode = "list", length = length(allmembers))
for (i in 1:length(allmembers)){
  newTDdata <- data.frame(Lat=DougScaledPres$Lat[i], Long=DougScaledPres$Long[i], TD=DougScaled$TD[TDmax,], TD2=DougScaled$TD2[TDmax,], PPT_sm=0, PPT_sm2=0, MWMT=0, MWMT2=0)
  pmaxTD[i] <- predict(fgamSVCtsr, newdata = newTDdata, type="response")
}


## PPT_sm ##
############

PPTsmmin <- min((na.exclude(Douglas)$PPT_sm))
PPTsmmin <- sample(which((Douglas$PPT_sm)==PPTsmmin),1) # 1 Nur ein PPTsmmin verwenden
pminPPTsm <- vector(mode = "list", length = length(allmembers))
for (i in 1:length(allmembers)){
  newPPTsmdata <- data.frame(Lat=DougScaledPres$Lat[i], Long=DougScaledPres$Long[i], PPT_sm=DougScaled$PPT_sm[PPTsmmin,], PPT_sm2=DougScaled$PPT_sm2[PPTsmmin,], TD=0, TD2=0, MWMT=0, MWMT2=0)
  pminPPTsm[i] <- predict(fgamSVCtsr, newdata = newPPTsmdata, type="response")
}

PPTsmmax <- max(na.exclude(Douglas)$PPT_sm)
PPTsmmax <- which((Douglas$PPT_sm)==PPTsmmax)
pmaxPPTsm <- numeric()
for (i in 1:length(allmembers)){
  newPPTsmdata <- data.frame(Lat=DougScaledPres$Lat[i], Long=DougScaledPres$Long[i], PPT_sm=DougScaled$PPT_sm[PPTsmmax,], PPT_sm2=DougScaled$PPT_sm2[PPTsmmax,], TD=0, TD2=0, MWMT=0, MWMT2=0)
  pmaxPPTsm[i] <- predict(fgamSVCtsr, newdata = newPPTsmdata, type="response")
}


## MWMT ####
############


MWMTmin <- min(Douglas$MWMT)
MWMTmin <- which(Douglas$MWMT==MWMTmin)
pminMWMT <- vector(mode = "list", length = length(allmembers))
for (i in 1:length(allmembers)){
  newMWMTdata <- data.frame(Lat=DougScaledPres$Lat[i], Long=DougScaledPres$Long[i], MWMT=DougScaled$MWMT[MWMTmin,], MWMT2=DougScaled$MWMT2[MWMTmin,], PPT_sm=0, PPT_sm2=0, TD=0, TD2=0)
  pminMWMT[i] <- predict(fgamSVCtsr, newdata = newMWMTdata, type="response")
}

MWMTmax <- max(Douglas$MWMT)
MWMTmax <- sample(which(Douglas$MWMT==MWMTmax),1)
pmaxMWMT <- numeric()
for (i in 1:length(allmembers)){
  newMWMTdata <- data.frame(Lat=DougScaledPres$Lat[i], Long=DougScaledPres$Long[i], MWMT=DougScaled$MWMT[MWMTmax,], MWMT2=DougScaled$MWMT2[MWMTmax,], PPT_sm=0, PPT_sm2=0, TD=0, TD2=0)
  pmaxMWMT[i] <- predict(fgamSVCtsr, newdata = newMWMTdata, type="response")
}


pTD <- cbind(pminTD, pmaxTD)
pPPTsm <- cbind(pminPPTsm, pmaxPPTsm)
pMWMT <- cbind(pminMWMT, pmaxMWMT)
pminmax <- cbind(pTD, pPPTsm)
pminmax <- cbind(pminmax, pMWMT)
pminmax <- as.data.frame(pminmax)
head(pminmax)
save(pminmax, file="Rdata/newcalculations/pminmax.Rdata")

TDminmax <- vector(mode = "list", length = length(nrow(pminmax)))
for (i in 1:nrow(pminmax)){
  if (pminmax[,1] > pminmax[,2]){
    TDminmax[i] <- min(Douglas$TD)
  } else  {
    TDminmax[i] <- max(Douglas$TD)
  }
}
PPTminmax <- vector(mode = "list", length = length(nrow(pminmax)))
for (i in 1:nrow(pminmax)){
  if (pminmax[,3] > pminmax[,4]){
    PPTminmax[i] <- min(Douglas$PPT_sm)
  } else  {
    PPTminmax[i] <- max(Douglas$PPT_sm)
  }
}
MWMTminmax <- vector(mode = "list", length = length(nrow(pminmax)))
for (i in 1:nrow(pminmax)){
  if (pminmax[,3] > pminmax[,4]){
    MWMTminmax[i] <- min(Douglas$MWMT)
  } else  {
    MWMTminmax[i] <- max(Douglas$MWMT)
  }
}

MinMax <- cbind(TDminmax, PPTminmax, MWMTminmax)
save(MinMax, file="Rdata/newcalculations/MinMax.Rdata")

