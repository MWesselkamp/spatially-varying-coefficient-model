DF_IR <- read.csv("~/Documents/Projects/Spatially-varying-coefficient-model/data/Doug_fir_MSc_Database.csv")
unique(DF_IR$PROV_GROUP)
length(unique(DF_IR$PROV_GROUP))
length(unique(DF_IR$PROV_NAME))

unique(DF_IR$PROV_LATITUDE)
unique(DF_IR$PROV_IUFRO_ID)
unique(DF_IR$SEEDZONE)
unique(DF_IR$PROV_ECOREG1)
unique(DF_IR$PROV_ECOREG4)

length(unique(DF_IR$PROV_MAP))
length(unique(DF_IR$PROV_MAT))
length(unique(DF_IR$PROV_Tmin12))

summary(DF_IR)

DFIR_small <- DF_IR[DF_IR$SITE_COUNTRY]