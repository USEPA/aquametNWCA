## ----mets.1--------------------------------------------------------------
library(aquametNWCA)

head(VascPlantEx)
exPlant <- prepareData(VascPlantEx, sampID='UID')

## ----mets.2--------------------------------------------------------------
mets <- calcVMMImets(exPlant$byUIDspp, sampID='UID')

head(mets)

# Now cast this data frame into wide format
library(reshape2)
mets.wide <- dcast(mets,UID~PARAMETER,value.var='RESULT')
head(mets.wide)


## ----mmi.1---------------------------------------------------------------
sites <- data.frame(UID=seq(1,10), NWCA_ECO4=c('CPL','CPL','IPL','CPL','CPL','IPL','CPL','CPL','CPL','CPL')
    ,NWCA_WET_GRP=c('PRLW','PRLW','PRLH','PRLW','PRLH','PRLW','EH','EH','PRLH','PRLH')
    ,ECO_X_WETGRP=c('CPL-PRLW','CPL-PRLW','IPL-PRLH','CPL-PRLW','CPL-PRLH','IPL-PRLW','ALL-EH'
                  ,'ALL-EH','CPL-PRLH','CPL-PRLH'),stringsAsFactors=F)

mets.wide.1 <- merge(sites, mets.wide, by='UID')

vmmi <- calcVMMI_fromMets(mets.wide.1, sampID='UID')
head(vmmi)


## ----mmi.2---------------------------------------------------------------
vmmi.alt <- calcVMMI_fromMets(mets.wide, sampID='UID')

head(vmmi.alt)

