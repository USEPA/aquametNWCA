## ----mets.1--------------------------------------------------------------
library(aquametNWCA)

head(VascPlantEx)
exPlant <- nwcaVegData(VascPlantEx, sampID='UID', cValReg='STATE')

## ----mets.2--------------------------------------------------------------
mets <- calcVMMImets(exPlant$byUIDspp, sampID='UID')

head(mets) # Already in wide format



## ----mmi.1---------------------------------------------------------------
sites <- data.frame(UID=seq(1,10), NWCA_ECO4=c('CPL','CPL','IPL','CPL','CPL','IPL','CPL','CPL','CPL','CPL')
    ,NWCA_WET_GRP=c('PRLW','PRLW','PRLH','PRLW','PRLH','PRLW','EH','EH','PRLH','PRLH')
    ,ECO_X_WETGRP=c('CPL-PRLW','CPL-PRLW','IPL-PRLH','CPL-PRLW','CPL-PRLH','IPL-PRLW','ALL-EH'
                  ,'ALL-EH','CPL-PRLH','CPL-PRLH'),stringsAsFactors=F)

mets.1 <- merge(sites, mets, by='UID')

vmmi <- calcVMMI_fromMets(mets.1, sampID='UID')
head(vmmi)


## ----mmi.2---------------------------------------------------------------
vmmi.alt <- calcVMMI_fromMets(mets, sampID='UID')

head(vmmi.alt)

