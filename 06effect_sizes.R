library(xtable)
library(dplyr)
library(tidyr)

#===================#
# Coefficient sizes #
#===================#

load(file="Rdata/fgamSVC.Rdata")
load(file="Rdata/fgamSVC_predTerms.Rdata")
load(file="Rdata/fgamSVC_clusters.Rdata")

useCluster <- clusters[[6]]$clustering
mean_responses <- vector(mode="list")
sd_responses <- vector(mode="list")
for (clus in unique(useCluster)){
  mean_responses[[clus]] <- apply(predTerms[useCluster==clus,], 2, mean)
  sd_responses[[clus]]  <- apply(predTerms[useCluster==clus,], 2, sd)
}

mean_responses <- do.call(rbind,mean_responses)
sd_responses <- do.call(rbind,sd_responses)

mean_responses
sd_responses

xtable(mean_responses, digits = 2)
xtable(sd_responses, digits = 2)

#==========================#
# Estimated response sizes #
#==========================#

estimated_responses <- predict(fgamSVC, newdata = DougScaledPres, type="response")

mean_responses <- vector(mode="list")
sd_responses <- vector(mode="list")
for (clus in unique(useCluster)){
  mean_responses[[clus]] <- mean(estimated_responses[useCluster==clus])
  sd_responses[[clus]]  <- sd(estimated_responses[useCluster==clus])
}

mean_responses <- do.call(rbind,mean_responses)
sd_responses <- do.call(rbind,sd_responses)

mean_responses
sd_responses

xtable(mean_responses, digits = 2)
xtable(sd_responses, digits = 2)

#============================#
# Conditional response sizes #
#============================#

load(file="Rdata/effectplots.Rdata")

responses <- df %>% 
  group_by(ID, Cluster) %>% 
  summarise(median = max(median)) %>% 
  pivot_wider(names_from = Cluster, values_from = median)

xtable(responses, digits = 2)

