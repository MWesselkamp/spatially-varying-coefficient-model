setwd("~/8787")
require(mgcv)

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

DougScaled <- Douglas
DougScaled[,c(8:13)] <- scale(Douglas[,c(8:13)])


## Spatially varying coefficient model with x/y ##
##################################################

fgamSVCtrs2 <- gam(PRES ~ s(y, x, k=100) + s(y, x, by=TD, k=100) + s(y, x, by=TD2, k=100) + s(y, x, by= PPT_sm, k=100) + s(y, x, by=PPT_sm2, k=100) + s(y, x, by=MWMT, k=100) + s(y, x, by=MWMT2, k=100), method="ML", family=binomial, control=gam.control(maxit=1000), data=DougScaled)

save(fgamSVCtrs2, file="Rdata/fgamSVCtrs2.Rdata")

