require(mgcv)

#setwd("~/Documents/Projects/Spatially-varying-coefficient-model")


## Prepare data ##
##################
 
Douglas <- read.csv("~/Spatially-varying-coefficient-model/data/DF_Plot_Data_Norm_6190.csv")
Douglas$PPT_sm[which(Douglas$PPT_sm == -1)] <- NA # remove incorrect value
Douglas <- na.omit(Douglas) 
row.names(Douglas) <- NULL # reset the rownames to index
Douglas$TD2 <- Douglas$TD^2
Douglas$MWMT2 <- Douglas$MWMT^2
Douglas$PPT_sm2 <- Douglas$PPT_sm^2
Douglas$PPT_wt2 <- Douglas$PPT_wt^2
Douglas$MDMP2 <- Douglas$MDMP^2

# only presences - unscaled
DougPres <- Douglas[Douglas$PRES==1,]

# data set complete - scaled
DougScaled <- Douglas
DougScaled[,c(10:22)] <- scale(Douglas[,c(10:22)])
summary(DougScaled)

DougScaled4 <- DougScaled
DougScaled4$x <- DougScaled4$x-min(DougScaled4$x)
DougScaled4$y <- DougScaled4$y-min(DougScaled4$y)

# data set presences - scaled
DougScaledPres <- DougScaled[DougScaled$PRES==1,]

DougSample <- DougPres[sample(1:nrow(DougPres), size=500),]
 
#==================================#
# Load a varying coefficient model #
#==================================#

mods = c("ref", "basic","non-standardized", "scaled-location")
  
for (mod in mods){
    if (mod == "ref"){
       load("Rdata/fgamtrs.Rdata")
       p = "r"
       
       # predicted occurance probabilities
        allPreds <- predict(fgamtrs, newdata=DougScaled, type="response")
  	save(allPreds, file=paste0("Rdata/allPreds", p, ".Rdata"))
  	
  	# model coefficients at presence observations
  	predTerms <- predict(fgamtrs, newdata=DougScaled[DougScaled$PRES==1,], type="terms")
  	save(predTerms, file=paste0("Rdata/predTerms", p, ".Rdata"))
  	
    }else if (mod=="basic"){
    	load("Rdata/fgamSVCtrs2.Rdata")
 	p = 2
  	
  	# predicted occurance probabilities
  	allPreds <- predict(fgamSVCtrs, newdata=DougScaled, type="response")
  	save(allPreds, file=paste0("Rdata/allPreds", p, ".Rdata"))
  	# model coefficients at presence observations
  	predTerms <- predict(fgamSVCtrs, newdata=DougScaled[DougScaled$PRES==1,], type="terms")
  	save(predTerms, file=paste0("Rdata/predTerms", p, ".Rdata"))
  	
    }else if (mod=="non-standardized"){
        load("Rdata/fgamSVCtrs3.Rdata")
 	p = 3
 	 
 	# predicted occurance probabilities
  	allPreds <- predict(fgamSVCtrs, newdata=Douglas, type="response")
  	save(allPreds, file=paste0("Rdata/allPreds", p, ".Rdata"))
  	# model coefficients at presence observations
  	predTerms <- predict(fgamSVCtrs, newdata=Douglas[Douglas$PRES==1,], type="terms")
  	save(predTerms, file=paste0("Rdata/predTerms", p, ".Rdata"))
  	
    }else if (mod=="scaled-location"){
        load("Rdata/fgamSVCtrs4.Rdata")
  	p = 4
  	
 	# predicted occurance probabilities
  	allPreds <- predict(fgamSVCtrs, newdata=DougScaled, type="response")
  	save(allPreds, file=paste0("Rdata/allPreds", p, ".Rdata"))
  	# model coefficients at presence observations
  	predTerms <- predict(fgamSVCtrs, newdata=DougScaled[DougScaled$PRES==1,], type="terms")
  	save(predTerms, file=paste0("Rdata/predTerms", p, ".Rdata"))
  	
    }else if (mod=="neutral"){
        load("Rdata/fgamSVCtrsN.Rdata")
        p = "N"
  	# predicted occurance probabilities
  	allPreds <- predict(fgamSVCtrs, newdata=DougScaled, type="response")
  	save(allPreds, file=paste0("Rdata/allPreds", p, ".Rdata"))
  	
  	# model coefficients at presence observations
  	predTerms <- predict(fgamSVCtrs, newdata=DougScaled[DougScaled$PRES==1,], type="terms")
  	save(predTerms, file=paste0("Rdata/predTerms", p, ".Rdata"))
    }

    #===================#
    # kmeans clustering #
    #===================#
    
    if (mod == "ref"){
       load("Rdata/predTerms.Rdata")
    }else if (mod=="basic"){
    	load("Rdata/predTerms2.Rdata")
    }else if (mod=="non-standardized"){
        load("Rdata/predTerms3.Rdata")
    }else if (mod=="scaled-location"){
        load("Rdata/predTerms4.Rdata")
    }else if (mod=="neutral"){
        load("Rdata/predTermsN.Rdata")
    }
    
    
    # Cluster without intercept.
    clusters <- vector(mode = "list", length = 15)
    for (i in 1:length(clusters)){
        clusters[[i]] <- kmeans(predTerms[,2:7], centers=i, nstart=25)
    }

    save(clusters, file=paste0("Rdata/clusters", p, ".Rdata"))

}