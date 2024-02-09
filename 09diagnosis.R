#===========#
# Load Data #
#===========#

Douglas <- get_Douglas_data()
# only presences - unscaled
DougPres <- Douglas[Douglas$PRES==1,]
# data set complete - scaled
DougScaled <- Douglas
DougScaled[,c(10:22)] <- scale(Douglas[,c(10:22)])
# data set presences - scaled
DougScaledPres <- DougScaled[DougScaled$PRES==1,]


load("Rdata/clusters2.Rdata")
useCluster = clusters2[[6]]$cluster
DougPres$Ecotype <- as.factor(useCluster)


RColorBrewer::brewer.pal(6, "Dark2")
#col.6.2 = c("#F3DF6C", "#CEAB07", "#798E87","#C93312", "#CCC591","#C27D38")

cols = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF")
cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')


ggplot(DougPres, aes(x=MWMT, fill = Ecotype, alpha=0.6)) + 
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = cols) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), 
        legend.position = "bottom",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=11), 
        legend.key.size = unit(0.5, units = "cm"), 
        legend.key = element_blank(),
        aspect.ratio = 1)  


ggplot(DougPres) + 
  geom_point(aes(x=MWMT, y = TD, colour = Ecotype), alpha=0.6) +
  scale_color_manual(values = cols) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="grey45", fill=NA, size=1), 
        panel.background = element_blank(), 
        legend.position = "bottom",
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.5, units = "cm"), 
        legend.key = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size=14),
        axis.text  = element_text(size=12)) 
