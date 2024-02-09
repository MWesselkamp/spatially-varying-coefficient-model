require(mgcv)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Projects/Spatially-varying-coefficient-model")

source('utils.R')
## Prepare data ##
##################

Douglas <- get_Douglas_data()

# only presences - unscaled
DougPres <- Douglas[Douglas$PRES==1,]

# data set complete - scaled
DougScaled <- Douglas
DougScaled[,c(10:19)] <- scale(Douglas[,c(10:19)])


## Spatially varying coefficient model with x/y ##
##################################################

fgamSVC <- gam(PRES ~ s(y, x, bs = "gp", k=100) + s(y, x, by=TD,bs = "gp", k=100) + s(y, x, by=TD2,bs = "gp", k=100) + s(y, x, by= PPT_sm,bs = "gp", k=100) + s(y, x, by=PPT_sm2, bs = "gp",k=100) + s(y, x, by=MWMT,bs = "gp", k=100) + s(y, x, by=MWMT2,bs = "gp", k=100), method="ML", family=binomial, control=gam.control(maxit=1000), data=DougScaled)
save(fgamSVC, file="Rdata/fgamSVC.Rdata")

summary(fgamSVC)
