require(ggplot2)
require(mgcv)
require(maps)
require(raster)
require(rgdal)
require(randomForest)
require(dismo)
require(gbm) # boosted regression trees
require(clusteval) # cluster similarity
require(mda) # discriminant function analysis
require(corrplot) # correlation plot
require(Hmisc) # varclus plot

source("utils.R")

load("Rdata/states.Rdata")
load("Rdata/states_mex.Rdata")
load("Rdata/states_ca.Rdata")

#=========================================#
# Variable selection with stationary GAM  #
#=========================================#

Douglas <- get_DouglasFir_data()

DouglasScaled = Douglas
DouglasScaled[,c(10:22, 37:43)] = scale(DouglasScaled[,c(10:22,37:43)])
DouglasScaledPres = DouglasScaled[which(DouglasScaled$PRES==1),]

# Take a subset eof the data for computational reasons
DouglasSample = Douglas[sample(nrow(Douglas), 10000),]

# Plot the spatial distribution of subsetted datapoints

countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)

ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="grey70") +
  geom_point(data=DouglasSample, aes(x=Long, y=Lat, col=factor(DouglasSample$PRES, levels = c(1,0))), size = 0.05) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", 
        legend.title = element_text(size=11, face="bold"), legend.text = element_text(size=11),
        legend.key=element_rect(fill = "white", size=0.5),
        legend.key.size = unit(5,"points"),
        plot.margin = margin(0.1,0.1,0.1,0.1, "cm")) + 
  guides(color = guide_legend(override.aes = list(size=3))) +
  scale_color_manual(values=c("#99CCFF", "#003366"), labels=c("present", "absent"),
                     name="Observation \nStatus")

# analyse the correlation between variables as a correlation plot and cluster
cormat = cor(as.matrix(DouglasSample[,10:36]))
corrplot(cormat, number.digits = 1)

DouglasSample = DouglasSample[,-c(1:9,16,12,22,26,29,31)]
DouglasSample = DouglasSample[,-c(9:15)]
DouglasSample = DouglasSample[,-c(10:12, 14)]
DouglasSample = DouglasSample[,-c(11:13)]
DouglasSample = DouglasSample[,-c(3,4)]
DouglasSample = DouglasSample[,1:6]
DouglasSample = DouglasSample[,-1]

v <- as.formula(paste("~", names(DouglasSample), collapse="+"))
par(mfrow=c(1,1))
plot(varclus(v, data=DouglasSample))
abline(h=0.5)

vmat = varclus(v, data = DouglasSample)


## Rescale the smaller data set
DouglasSampleScaled = cbind(DouglasSample[,c(1:9)], DouglasSample$MWMT, DouglasSample$MWMT2, DouglasSample$PPT_sm, DouglasSample$PPT_sm2, DouglasSample$TD, DouglasSample$TD2, DouglasSample$PPT_wt, DouglasSample$MDMP)
names(DouglasSampleScaled)[10:17] = c("MWMT","MWMT2", "PPT_sm", "PPT_sm2", "TD", "TD2", "PPT_wt", "MDMP")

DouglasSampleScaled$PPT_wt2 = DouglasSampleScaled$PPT_wt^2
DouglasSampleScaled$MDMP2 = DouglasSampleScaled$MDMP^2

DouglasSampleScaled[,c(10:19)] = scale(DouglasSampleScaled[,c(10:19)])


## Run stationary GAMs with the leftover non-correlated variables ##

# Run simplest GAM
fm = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)

save(fm, file="Rdata/fgamStatSimple.Rdata")
summary(fm)
plot(fm)

# Use all five variables and quadrativ terms

fm = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100) + s(PPT_wt, k=100) + s(PPT_wt2, k=100) + s(MDMP, k=100) + s(MDMP2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasScaled)

save(fm, file="Rdata/fgamStatSimple2.Rdata")

load("Rdata/reanalysis/fgamStatSimple2.Rdata")
summary(fm)
plot(fm)

fm1 = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100) + s(PPT_wt, k=100) + s(PPT_wt2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)

fm2 = gam(PRES ~ s(TD, k=100) + s(TD2, k=100) + s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)


## Variable importance ##
# only with non-correlated variables

## RandomForest
rf = randomForest(as.factor(PRES) ~ TD + MWMT + PPT_wt  + PPT_sm + MDMP, data=DouglasSample)
randomForest::varImpPlot(rf)
save(rf, file="Rdata/randomForest.Rdata")

## Boosted regression trees
# same variables. Test for interactions (reviewer recommendations)

brt1 = gbm.step(data=DouglasSample, gbm.x = c(11,15,17,18,19), gbm.y = 7, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)
save(brt1, file="Rdata/regressiontrees.Rdata")

summary(brt1)
sort(relative.influence(brt1))

## Another simple stationary gam

fm = gam(PRES ~ s(MWMT, k=100) + s(MWMT2, k=100) + s(PPT_sm, k=100) + s(PPT_sm2, k=100) + s(PPT_wt, k=100) + s(PPT_wt2, k=100), control = gam.control(maxit = 1000), family = binomial, method="ML", data=DouglasSampleScaled)

save(fm, file="Rdata/fgamStatSimple3.Rdata")
summary(fm)
plot(fm)


## Spatially variable coefficient model with five predictors?


#====================#
# Graphical Abstract #
#====================#

useCluster <- clusters2[[6]]$cluster
col.l <-  c("#F3DF6C", "#CEAB07", "#798E87","#C93312", "#CCC591","#C27D38") # coffee colors

p4 <- ggplot() + 
  geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
  geom_path(data= states, aes(x=long, y=lat, group=group), col="grey25") + 
  geom_path(data=ca_provinces,aes(x=long, y=lat, group=group), col="grey25")  + 
  geom_path(data=mex_states, aes(x=long, y=lat, group=group), col="grey25") +
  geom_point(data=DougScaledPres, aes(x=Long, y=Lat, col=factor(useCluster)), size=0.8) + 
  coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) + 
  scale_color_manual(values = c(col.l), name = "Ecotypes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),legend.text = element_text(size=11), legend.key.size = unit(0.5, units = "cm"), legend.key = element_blank())  +
  guides(color = guide_legend(override.aes = list(size=2)))

pdf("figures/graphicalAbstract.pdf")
p4
dev.off()
