#===================#
# Utility functions #
#===================#
library(maps)
library(raster)
library(rgdal)

get_Douglas_data <- function(){
  
  Douglas <- read.csv("data/DF Plot Data (Norm_6190 Climate).csv")
  Douglas$PPT_sm[which(Douglas$PPT_sm == -1)] <- NA # remove incorrect value 
  Douglas <- na.omit(Douglas)
  row.names(Douglas) <- NULL # reset the rownames to index
  Douglas$TD2 <- Douglas$TD^2
  Douglas$MWMT2 <- Douglas$MWMT^2
  Douglas$PPT_sm2 <- Douglas$PPT_sm^2
  Douglas$PPT_wt2 <- Douglas$PPT_wt^2
  Douglas$MDMP2 <- Douglas$MDMP^2
  
  return(Douglas)
}

get_state_borders <- function(){
  
  #---- get states borders for USA, Canada and Mexico
  provinces_ca <- map("worldHires", "Canada")
  states <- raster::getData("GADM", country="USA", level=1)
  save(states, file="output/states.Rdata")
  
  ca <- raster::getData("GADM", country="Canada", level=1)
  provinces <- c("British Columbia", "Alberta", "Saskatchewan")
  ca_smaller <- ca[ca$NAME_1 %in% provinces,] # only use relevant canadian states to save memory
  
  mex <- raster::getData("GADM", country="Mexico", level=1)
  
  ca_smaller@data$id <- rownames(ca_smaller@data) # define an explicit relationship between data and polygons associated with the data
  ca_provinces <- fortify(ca_smaller, region = "id") # convert class of spdf into df
  ca_provincesDF <- merge(ca_provinces, ca_smaller@data, by="id") # merge if some important information is missing
  save(ca_provinces, file="output/ca_provinces.Rdata")
  
  mex@data$id <- rownames(mex@data)
  mex_states <- fortify(mex, region="id") 
  save(mex_states, file="output/mex_states.Rdata")
}

get_DNA_data <- function(){
  
  DNAdata = read.csv("data/Doug-Fir DNA data (Wei et al, 2011).csv")
  
  bioclim = getData("worldclim", var="bio", res = 5)
  mypoints = data.frame(long=DNAdata$Longitude, lat=DNAdata$Latitude)
  myvars = as.data.frame(extract(bioclim, mypoints))
  myvars_small = data.frame(TD = myvars$bio7/10, MWMT = myvars$bio10/10, PPT_sm = myvars$bio18)
  
  ## prepare data for SVCM
  DNAdata = cbind(DNAdata, myvars_small)
  DNAdata$TD2 = DNAdata$TD^2
  DNAdata$MWMT2 = DNAdata$MWMT^2
  DNAdata$PPT_sm2 = DNAdata$PPT_sm^2
  
  ## convert latitude and longitude to UTM coordinates.
  xy = data.frame(x = DNAdata_scaled[,5], y = DNAdata_scaled[,4])
  coordinates(xy) = c("x", "y")
  proj4string(xy) = CRS("+proj=longlat +datum=WGS84")
  res = spTransform(xy, "+proj=utm +datum=WGS84")

  DNAdata = cbind(DNAdata, res@coords)
  
  return(DNAdata)
}


kmeans_clustering = function(coeffs, n_clusters){
  
  clustered.predictions <- vector(mode = "list", length = n_clusters)
  
  for (i in 1:length(clustered.predictions)){
    clustered.predictions[[i]] <- kmeans(coeffs, centers=i, nstart=25)
  }
  
  wss = numeric(length(clustered.predictions))
  for (i in 2:length(clustered.predictions)){
    wss[i] = clustered.predictions[[i]]$tot.withinss/clustered.predictions[[i]]$totss
  }
  plot(2:length(clustered.predictions)+1, wss[2:length(clustered.predictions)+1], type="b", ylim=c(0,1), xlab="Number of Clusters", ylab="Total within groups sum of squares")
  return(clustered.predictions)
  
}


confusion_matrix = function(genotype, ecotype){
  
  if(is.factor(genotype)){
    genotype = as.integer(genotype)-1
  }
  mat = matrix(0, nrow=length(unique(ecotype)), ncol=length(unique(genotype)))
  
  for (i in 1:length(ecotype)){
    mat[ecotype[i], genotype[i]] = mat[ecotype[i], genotype[i]]+1
  }
  
  return(mat)
}

# percentage of ecotypic presences assigned to neutral cluster

inpercentageRows = function(m, sums){
  for(i in 1:nrow(m)){
    m[i,] = round((m[i,]/sums[i])*100, 2)
  }
  return(m)
}

# percentage of neutral clusters assigned to ecotypes

inpercentageCols = function(m, sums){
  for(i in 1:ncol(m)){
    m[,i] = round((m[,i]/sums[i])*100, 2)
  }
  return(m)
}

plot_wss <- function(cluster_list){
  
  wss = numeric(14)
  for (i in 1:14){
    wss[i] = cluster_list[[i+1]]$tot.withinss  
  }
  plot(1:14, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  
}