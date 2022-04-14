require(mgcv)
require(parallel)

## Prepare data ##
##################

Douglas <- read.csv("~/svcm/DF_Plot_Data_Norm_6190.csv")
Douglas$PPT_sm[which(Douglas$PPT_sm == -1)] <- NA # remove incorrect value 
Douglas <- na.omit(Douglas)
row.names(Douglas) <- NULL # reset the rownames to index

Douglas <- Douglas[,c(1:5,7, 9, 11, 15, 17)]
Douglas$TD2 <- Douglas$TD^2
Douglas$MWMT2 <- Douglas$MWMT^2
Douglas$PPT_sm2 <- Douglas$PPT_sm^2

DougScaled <- Douglas
DougScaled[,c(8:13)] <- scale(Douglas[,c(8:13)])
DougScaled$x <- (DougScaled$x - min(DougScaled$x))
DougScaled$y <- (DougScaled$y - min(DougScaled$y))

nc <- detectCores()
print(nc)
if (detectCores()>1) { ## no point otherwise
  cl <- makeCluster(nc-1)
  ## could also use makeForkCluster, but read warnings first!
} else cl <- NULL


## GAM without spatially varying coefficients ##
################################################

system.time(fgamtrs <- gam(PRES ~ s(TD, bs = "ts", k=100) +
                    s(TD2, bs ="ts", k=100) + s(PPT_sm, bs = "ts", k=100) +
                    s(PPT_sm2, bs = "ts", k=100) + s(MWMT, bs = "ts", k=100) +
                    s(MWMT2, bs = "ts", k=100),
                  method="ML", family=binomial, control=gam.control(maxit=1000, nthreads=nc),
                  data=DougScaled))

save(fgamtrs, file="~/svcm/fgamtrs.Rdata")

## Spatially varying coefficient model with x/y ##
##################################################

system.time(fgamSVCtrs <- gam(PRES ~ s(y, x, bs = "ts", k=100) + s(y, x, by=TD, bs = "ts", k=100) + 
                    s(y, x, by=TD2, bs ="ts", k=100) + s(y, x, by= PPT_sm, bs = "ts", k=100) + 
                    s(y, x, by=PPT_sm2, bs = "ts", k=100) + s(y, x, by=MWMT, bs = "ts", k=100) + 
                    s(y, x, by=MWMT2, bs = "ts", k=100), 
                  method="ML", family=binomial, control=gam.control(maxit=1000, nthreads=nc), 
                  data=DougScaled))

#save(fgamSVCtrs, file="~/svcm/fgamSVCtrs2.Rdata")

## Use not-standardized climate predictors ##
#############################################


system.time(fgamSVCtrs <- gam(PRES ~ s(y, x, bs = "ts", k=100) + s(y, x, by=TD, bs = "ts", k=100) +
                    s(y, x, by=TD2, bs ="ts", k=100) + s(y, x, by= PPT_sm, bs = "ts", k=100) +
                    s(y, x, by=PPT_sm2, bs = "ts", k=100) + s(y, x, by=MWMT, bs = "ts", k=100) +
                    s(y, x, by=MWMT2, bs = "ts", k=100),
                  method="ML", family=binomial, control=gam.control(maxit=1000, nthreads=nc),
                  data=Douglas))

save(fgamSVCtrs, file="~/svcm/fgamSVCtrs3.Rdata")


## Move x and y to 0 at the coordinate system ##
################################################

DougScaled$x <- (DougScaled$x - mean(DougScaled$x))
DougScaled$y <- (DougScaled$y - mean(DougScaled$y))


system.time(fgamSVCtrs <- gam(PRES ~ s(y, x, bs = "ts", k=100) + s(y, x, by=TD, bs = "ts", k=100) +
                    s(y, x, by=TD2, bs ="ts", k=100) + s(y, x, by= PPT_sm, bs = "ts", k=100) +
                    s(y, x, by=PPT_sm2, bs = "ts", k=100) + s(y, x, by=MWMT, bs = "ts", k=100) +
                    s(y, x, by=MWMT2, bs = "ts", k=100),
                  method="ML", family=binomial, control=gam.control(maxit=1000, nthreads=nc),
                  data=DougScaled))

save(fgamSVCtrs, file="~/svcm/fgamSVCtrs4.Rdata")