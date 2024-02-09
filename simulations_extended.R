library(dplyr)
library(mgcv) # gam
library(fossil) # rand.index
library(cluster) # GAP statistics
library(factoextra) # visualize clusters

set.seed(42)

Ecotypes = function(dat, cluster){
  
  cols = c('#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#e6194b', '#3cb44b',    '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')
  
  p4 = ggplot() + 
    geom_polygon(data=northA, aes(x=long, y=lat, group=group), fill="black") + 
    geom_point(data = dat, aes(x=Long, y=Lat, col=factor(cluster)), size=0.8, alpha=1) + 
    coord_fixed(ratio=1, xlim = c(-130, -95), ylim=c(15, 55)) +
    scale_color_manual(values = cols, name="Ecotypes") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour="grey45", fill=NA, size=1), 
          panel.background = element_blank(), axis.title = element_blank(), 
          axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
          legend.title = element_text(size=16, face="bold"),legend.text = element_text(size=14, face="bold"), legend.key.size = unit(0.8, units = "cm"), legend.key = element_blank()) +
    guides(color = guide_legend(override.aes = list(size=4)))
  
  return(p4) 
}

linkLinearResponse <- function(linearResponse, n){
  family = binomial()
  linkResponse = family$linkinv(linearResponse)
  observedResponse = rbinom(n = n, 1, prob = linkResponse)
  hist(observedResponse)
  return(observedResponse)
}

countries <- map_data("world")
northA <- subset(countries, region %in% c("USA", "Mexico", "Canada") & long < 170)

source("utils.R")

Douglas <- get_Douglas_data()
DougScaled <- Douglas
DougScaled[,c(10:19)] <- scale(Douglas[,c(10:19)])
DougScaledPres <- DougScaled[DougScaled$PRES==1,]

# number of samples from the coordinates
n_samples <- 30000

# Sample from data set
doug_rs <- Douglas %>% 
  sample_n(n_samples) %>% 
  select(ID, y, x, Lat, Long, PRES, POP, MWMT, TD, PPT_sm) %>% 
  mutate(POP_num = (as.numeric(as.factor(POP)))) %>% 
  mutate(MWMT2 = MWMT^2,TD2 = TD^2,PPT_sm2 = PPT_sm^2) %>% 
  mutate(y = (y-min(y))/1000, x = (x-min(x))/1000) # center UTM coords and transform to kilometres, https://www.sciencedirect.com/science/article/pii/S2211675320300646?via%3Dihub#sec6

doug_rs_scaled <- doug_rs
doug_rs_scaled <- doug_rs_scaled[-sample(which(doug_rs_scaled$PRES==0), size=7000),]
doug_rs_scaled[,c(8:10,12:14)] <- scale(doug_rs_scaled[,c(8:10,12:14)])
# downsample zeros in the sample.
hist(doug_rs_scaled$PRES)
sum(doug_rs_scaled$PRES)

corrplot::corrplot(cor(doug_rs_scaled[,c(8:10)]))

#===================================#
# Fit six different ecotypic models #
#===================================#

# predict to actual presences
doug_rs_scaled_pres <- doug_rs_scaled[doug_rs_scaled$PRES==1,]

linearResponse1 <- -14.00*doug_rs_scaled$TD + 17.18*doug_rs_scaled$TD2 -1.38*doug_rs_scaled$PPT_sm + 3.35*doug_rs_scaled$PPT_sm2 - 36.99*doug_rs_scaled$MWMT + 47.65*doug_rs_scaled$MWMT2 
observedResponse1 <- linkLinearResponse(linearResponse1, length(linearResponse1))
fm1 <- gam(observedResponse1 ~ s(x,y, bs = "gp") + s(x,y, by = MWMT, bs = "gp") + s(x,y, by = TD, bs = "gp") + s(x,y, by = PPT_sm, bs = "gp") + s(x,y, by = MWMT2, bs = "gp") + s(x,y, by = TD2, bs = "gp") + s(x,y, by = PPT_sm2, bs = "gp"), data = doug_rs_scaled,method="ML", family = binomial)
save(fm1, file = "Rdata/simulations_fm1.Rdata")
summary(fm1)
# without effect of location: R-sq.(adj) =  0.666   Deviance explained = 61.4%
# if we add a small effect of location: R-sq.(adj) =  0.697   Deviance explained = 65.1%
preds1 <- predict(fm1, newdata = doug_rs_scaled, type="response")
hist(preds1)

linearResponse2 <- -16.49*doug_rs_scaled$TD + 19.55*doug_rs_scaled$TD2 -3.53*doug_rs_scaled$PPT_sm + 6.94 *doug_rs_scaled$PPT_sm2 -16.06*doug_rs_scaled$MWMT + 23.02*doug_rs_scaled$MWMT2 
observedResponse2 <- linkLinearResponse(linearResponse2, length(linearResponse2))
fm2 <- gam(observedResponse2 ~ s(x,y, bs = "gp") + s(x,y, by = MWMT, bs = "gp") + s(x,y, by = TD, bs = "gp") + s(x,y, by = PPT_sm, bs = "gp") + s(x,y, by = MWMT2, bs = "gp") + s(x,y, by = TD2, bs = "gp") + s(x,y, by = PPT_sm2, bs = "gp"), data = doug_rs_scaled,method="ML", family = binomial)
save(fm2, file = "Rdata/simulations_fm2.Rdata")
summary(fm2)
preds2 <- predict(fm2, newdata = doug_rs_scaled, type="response")
hist(preds2)

linearResponse3 <- -14.21*doug_rs_scaled$TD + 19.52*doug_rs_scaled$TD2 -14.91*doug_rs_scaled$PPT_sm +
  29.95*doug_rs_scaled$PPT_sm2 + 1.77*doug_rs_scaled$MWMT + 1.68*doug_rs_scaled$MWMT2
observedResponse3 <- linkLinearResponse(linearResponse3, length(linearResponse3))
fm3 <- gam(observedResponse3 ~ s(x,y, bs = "gp") + s(x,y, by = MWMT, bs = "gp") + s(x,y, by = TD, bs = "gp") + s(x,y, by = PPT_sm, bs = "gp") + s(x,y, by = MWMT2, bs = "gp") + s(x,y, by = TD2, bs = "gp") + s(x,y, by = PPT_sm2, bs = "gp"), data = doug_rs_scaled,method="ML", family = binomial)
save(fm3, file = "Rdata/simulations_fm3.Rdata")
summary(fm3)
preds3 <- predict(fm3, newdata = doug_rs_scaled, type="response")
hist(preds3)

linearResponse4 <- 0.13*doug_rs_scaled$TD + 2.24*doug_rs_scaled$TD2 -2.99*doug_rs_scaled$PPT_sm +
  6.19*doug_rs_scaled$PPT_sm2 - 4.69*doug_rs_scaled$MWMT + 8.19*doug_rs_scaled$MWMT2
observedResponse4 <- linkLinearResponse(linearResponse4, length(linearResponse4))
fm4 <- gam(observedResponse4 ~ s(x,y, bs = "gp") + s(x,y, by = MWMT, bs = "gp") + s(x,y, by = TD, bs = "gp") + s(x,y, by = PPT_sm, bs = "gp") + s(x,y, by = MWMT2, bs = "gp") + s(x,y, by = TD2, bs = "gp") + s(x,y, by = PPT_sm2, bs = "gp"), data = doug_rs_scaled, method="ML",family = binomial)
save(fm4, file = "Rdata/simulations_fm4.Rdata")
summary(fm4)
preds4 <- predict(fm4, newdata = doug_rs_scaled, type="response")
hist(preds4)

linearResponse5 <- -1.07*doug_rs_scaled$TD + 3.58*doug_rs_scaled$TD2 -2.42*doug_rs_scaled$PPT_sm + 4.97*doug_rs_scaled$PPT_sm2 - 14.38*doug_rs_scaled$MWMT + 19.33*doug_rs_scaled$MWMT2
observedResponse5 <- linkLinearResponse(linearResponse5, length(linearResponse5))
fm5 <- gam(observedResponse5 ~  s(x,y, bs = "gp") + s(x,y, by = MWMT, bs = "gp") + s(x,y, by = TD, bs = "gp") + s(x,y, by = PPT_sm, bs = "gp") + s(x,y, by = MWMT2, bs = "gp") + s(x,y, by = TD2, bs = "gp") + s(x,y, by = PPT_sm2, bs = "gp"), data = doug_rs_scaled,method="ML", family = binomial)
save(fm5, file = "Rdata/simulations_fm5.Rdata")
summary(fm5)
preds5 <- predict(fm5, newdata = doug_rs_scaled, type="response")
hist(preds5)

linearResponse6 <- -1.91*doug_rs_scaled$TD + 4.58*doug_rs_scaled$TD2 -2.45*doug_rs_scaled$PPT_sm + 5.06*doug_rs_scaled$PPT_sm2 - 24.17*doug_rs_scaled$MWMT + 30.60*doug_rs_scaled$MWMT2
observedResponse6 <- linkLinearResponse(linearResponse6, length(linearResponse6))
fm6 <- gam(observedResponse6 ~ s(x,y, bs = "gp") + s(x,y, by = MWMT, bs = "gp") + s(x,y, by = TD, bs = "gp") + s(x,y, by = PPT_sm, bs = "gp") + s(x,y, by = MWMT2, bs = "gp") + s(x,y, by = TD2, bs = "gp") + s(x,y, by = PPT_sm2, bs = "gp"), data = doug_rs_scaled,method="ML", family = binomial)
save(fm6, file = "Rdata/simulations_fm6.Rdata")
summary(fm6)
preds6 <- predict(fm6, newdata = doug_rs_scaled, type="response")
hist(preds6)


pmat <- data.frame(ID = doug_rs_scaled$ID, x = doug_rs_scaled$x, y = doug_rs_scaled$y, 'A' = 0, 'B' = 0, 'C' = 0, 'D' = 0, 'E' = 0, 'F' = 0, 'N' = 0)
pmat$A <- preds1 # [doug_rs_scaled$PRES==1,]
pmat$B <- preds2 # [doug_rs_scaled$PRES==1,]
pmat$C <- preds3 # usw..
pmat$D <- preds4
pmat$E <- preds5
pmat$F <- preds6

pmat_pres <- pmat[doug_rs_scaled$PRES==1,]
summary(pmat_pres)
ecotypes <- apply(pmat_pres[,c(4:9)], 1, function(x) sample(LETTERS[1:6], size = 1,prob = x))
#ecotypes <- apply(pmat[,c(4:9)], 1, function(x) sample(LETTERS[1:6], size = 1,prob = x))

Ecotypes(doug_rs_scaled_pres, cluster = ecotypes)
#Ecotypes(doug_rs_scaled, cluster = ecotypes)
#doug_rs_scaled$POP_sim <- ecotypes
#=======================================#
# Create SVCM model and extract cluster #
#=======================================#

doug_rs_scaled$POP_sim <- 'N'
doug_rs_scaled$POP_sim[doug_rs_scaled$PRES==1] <- ecotypes
doug_rs_scaled_pres <- doug_rs_scaled[doug_rs_scaled$PRES==1,]
barplot(table(doug_rs_scaled$POP_sim), xlab = "Ecotype", ylab = "Frequency")
Ecotypes(doug_rs_scaled_pres, cluster = doug_rs_scaled_pres$POP_sim)

# Add simulated presences to data: PRES_sim
# add additional stochasticity?
prob_means <- apply(pmat[,c(4:9)], 2, mean)
PRES_sim <- numeric(length = length(ecotypes))
for (i in 1:length(ecotypes)){
  PRES_sim[i] <- rbinom(1, 1, prob = pmat_pres[i, ecotypes[i]]) #/prob_means[ecotypes[i]]
}
sum(PRES_sim)

# or maybe not?
doug_rs_scaled$PRES_sim <- 0 #PRES_sim# 0
doug_rs_scaled$PRES_sim[doug_rs_scaled$POP_sim !='N'] <- PRES_sim
doug_rs_scaled_pres <- doug_rs_scaled[doug_rs_scaled$PRES_sim==1,]
Ecotypes(doug_rs_scaled_pres, cluster = doug_rs_scaled_pres$POP_sim)
# maybe remove the random zeros
doug_rs_scaled_sim <- doug_rs_scaled[-which((doug_rs_scaled$POP_sim !='N') & (doug_rs_scaled$PRES_sim==0)),]

# use gaussian processes: https://agile-giss.copernicus.org/articles/3/31/2022/agile-giss-3-31-2022.pdf
fmSVCM <- mgcv::gam(PRES_sim ~ s(x,y, bs = "gp") + s(x,y, by = MWMT, bs = "gp") + s(x,y, by = TD, bs = "gp") + s(x,y, by = PPT_sm, bs = "gp") + s(x,y, by = MWMT2, bs = "gp") + s(x,y, by = TD2, bs = "gp") + s(x,y, by = PPT_sm2, bs = "gp"), data = doug_rs_scaled_sim, method="ML", family = binomial)
summary(fmSVCM) 
length(fmSVCM$smooth)
# If we add stochasticity and keep zeros (PRES_sim): R-sq.(adj) =  0.422   Deviance explained = 45.5%
# If we add stochasticity and remove zeros (PRES_sim =! 0): R-sq.(adj) =  0.697   Deviance explained = 68.5%
# If we use all raw presences (1): # R-sq.(adj) =   0.72   Deviance explained =   69%, n = 8000

# Extract spatially varying coefficients where Douglas Fir is present.
doug_rs_scaled_pres <- doug_rs_scaled[doug_rs_scaled$PRES_sim==1,]
predTerm <- predict(fmSVCM, newdata = doug_rs_scaled_pres, type="terms")

gaps <- clusGap(predTerm, FUN = kmeans, K.max = 20, B=100)
par(mar=c(5,5,4,1))
plot(gaps, las=1, xpd=T)

cluster <-  pam(predTerm, k=6, medoids = 'random', variant = 'faster')
fviz_cluster(cluster, data=predTerm, ellipse.type = 'convex')
useCluster <- cluster$clustering

rand.index(doug_rs_scaled_pres$POP_sim, useCluster)
table(doug_rs_scaled_pres$POP_sim, useCluster)

Ecotypes(doug_rs_scaled_pres, cluster = as.factor(useCluster))
coefficient_maps(predTerm, dat = doug_rs_scaled_pres)

#======================#
# Plot response curves #
#======================#

compute_response_curves = function(j){
  
  newTD <- seq(min(doug_rs$TD), max(doug_rs$TD), len=100)
  newTD2 <- newTD^2
  newPPTsm <- seq(min(doug_rs$PPT_sm), max(doug_rs$PPT_sm), len=100)
  newPPTsm2 <- newPPTsm^2
  newMWMT <- seq(min(doug_rs$MWMT), max(doug_rs$MWMT), len=100)
  newMWMT2 <- newMWMT^2
  
  newTDscaled <- (newTD - mean(doug_rs$TD))/sd(doug_rs$TD)
  newTD2scaled <- (newTD2 - mean(doug_rs$TD2))/sd(doug_rs$TD2)
  newPPTsmscaled <- (newPPTsm - mean(doug_rs$PPT_sm))/sd(doug_rs$PPT_sm)
  newPPTsm2scaled <- (newPPTsm2 - mean(doug_rs$PPT_sm2))/sd(doug_rs$PPT_sm2)
  newMWMTscaled <- (newMWMT - mean(doug_rs$MWMT))/sd(doug_rs$MWMT)
  newMWMT2scaled <- (newMWMT2 - mean(doug_rs$MWMT2))/sd(doug_rs$MWMT2)
  
  clustmembers <- which(useCluster==j)
  cluster <- paste("Cluster", j, sep=" ")
  
  effTD1 <- vector(mode="list", length = length(clustmembers) )
  
  for (i in 1:length(clustmembers)){
    
    newTDdata <- data.frame(y=doug_rs_scaled$y[clustmembers[i]], x=doug_rs_scaled$x[clustmembers[i]], TD=newTDscaled, TD2=newTD2scaled, PPT_sm=mean(doug_rs_scaled$PPT_sm[clustmembers]), PPT_sm2=mean(doug_rs_scaled$PPT_sm2[clustmembers]), MWMT=mean(doug_rs_scaled$MWMT[clustmembers]), MWMT2=mean(doug_rs_scaled$MWMT2[clustmembers]))
    effTD1[[i]] <- predict.gam(fmSVCM, newdata = newTDdata, type="response")
    
  }
  
  effects <- do.call(rbind, effTD1)
  
  df <- as.data.frame(t(apply(effects, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles
  colnames(df) <- c("first", "second", "median", "third", "fourth")
  df$mean <- apply(effects, 2, mean)
  df$ID <- "TD"
  df$Cluster <- cluster
  df$newseq <- newTD
  
  effPPTsm <- vector(mode="list", length = length(clustmembers) )
  
  for (i in 1:length(clustmembers)){
    
    newPPTsmdata <- data.frame(y=doug_rs_scaled$y[clustmembers[i]], x=doug_rs_scaled$x[clustmembers[i]], TD=mean(doug_rs_scaled$TD[clustmembers]), TD2=mean(doug_rs_scaled$TD2[clustmembers]), PPT_sm=newPPTsmscaled, PPT_sm2=newPPTsm2scaled, MWMT=mean(doug_rs_scaled$MWMT[clustmembers]), MWMT2=mean(doug_rs_scaled$MWMT2[clustmembers]))
    effPPTsm[[i]] <- predict.gam(fmSVCM, newdata = newPPTsmdata, type="response")
    
  }
  
  effects <- do.call(rbind, effPPTsm)
  
  df2 <- as.data.frame(t(apply(effects, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles
  colnames(df2) <- c("first", "second", "median", "third", "fourth")
  df2$mean <- apply(effects, 2, mean)
  df2$ID <- "PPT_sm"
  df2$Cluster <- cluster
  df2$newseq <- newPPTsm
  
  effMWMT <- vector(mode="list", length = length(clustmembers) )
  
  for (i in 1:length(clustmembers)){
    newMWMTdata <- data.frame(y=doug_rs_scaled$y[clustmembers[i]], x=doug_rs_scaled$x[clustmembers[i]], TD=mean(doug_rs_scaled$TD[clustmembers]), TD2=mean(doug_rs_scaled$TD2[clustmembers]), PPT_sm=mean(doug_rs_scaled$PPT_sm[clustmembers]), PPT_sm2=mean(doug_rs_scaled$PPT_sm2[clustmembers]), MWMT=newMWMTscaled, MWMT2=newMWMT2scaled)
    effMWMT[[i]] <- predict.gam(fmSVCM, newdata = newMWMTdata, type="response")
    
  }
  
  effects <- do.call(rbind, effMWMT)
  df3 <- as.data.frame(t(apply(effects, 2, quantile, c(0.95, 0.75, 0.5, 0.25, 0.05)))) # generate quantiles
  colnames(df3) <- c("first", "second", "median", "third", "fourth")
  df3$mean <- apply(effects, 2, mean)
  df3$ID <- "MTWM"
  df3$Cluster <- cluster
  df3$newseq <- newMWMT
  
  df_reponses[[j]] <- rbind(df, df2, df3)
  
}


effectplots = function(df){

  
  df$ID_f <- factor(df$ID, levels = c("MTWM","TD","PPT_sm"), labels = c("MTWM (°C)", "TD (°C)", "PPT_sm (mm)"))
  df$Cluster <- gsub(df$Cluster, pattern = "Cluster", replacement = "Ecotype")

  cols = c('#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#e6194b', '#3cb44b',    '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')
  
  p = ggplot(df) + ylim(0,1) + 
    geom_path(aes(x=newseq, y=mean, col=as.factor(Cluster))) + labs(y="Occurrence Probability")+
    geom_ribbon(aes(ymin = first, ymax = fourth, x=newseq), fill="grey50" , alpha=0.3) + 
    theme(aspect.ratio = 1, axis.title.x = element_blank(), strip.text.x = element_text(face="bold"),
          strip.text.y = element_text(face="bold"),panel.background = element_blank(), 
          legend.key = element_blank(), legend.position = "none") + 
    facet_grid(Cluster~ID_f, scales="free_x") + 
    scale_color_manual(values = cols, name="Ecotypes")

  return(p)
  
}

df_reponses <- vector(mode = "list", length=6)

for (i in 1:6){
  df_reponses[[i]] = compute_response_curves(i)
}

df = do.call(rbind, df_reponses)
summary(df)

effectplots(df)
