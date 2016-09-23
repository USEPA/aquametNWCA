## ----allVasc.1-----------------------------------------------------------
library(aquametNWCA2)

head(VascPlantEx)

## ----allVasc.2-----------------------------------------------------------
outdf <- calcVascPlantMets(dfIn=VascPlantEx,taxaIn=taxaNWCA,taxaCC=ccNatNWCA,taxaWIS=wisNWCA)

names(outdf)

## ----subVasc.1-----------------------------------------------------------
sumdf <- prepareData(indf=VascPlantEx, inTaxa=taxaNWCA, inNatCC=ccNatNWCA, inWIS=wisNWCA)
str(sumdf)

## ----subVasc.2-----------------------------------------------------------
outRich <- calcRichness(sumdf$byUIDspp,sumdf$byPlotspp,sumdf$byUIDgen,sumdf$byPlotgen,sumdf$byUIDfam,sumdf$byPlotfam)
head(outRich)
unique(outRich$PARAMETER)

## ----subVasc.3-----------------------------------------------------------
outCat <- calcCategory(sumdf$byUIDspp)
head(outCat)
unique(outCat$PARAMETER) 

## ----subVasc.4-----------------------------------------------------------
outBC <- calcBCmets(sumdf$byPlotspp)
head(outBC)
unique(outBC$PARAMETER)

## ----subVasc.5-----------------------------------------------------------
allVasc <- rbind(outRich, outCat, outBC)

library(reshape2)
outVasc <- dcast(allVasc,UID~PARAMETER,value.var='RESULT')
head(outVasc)

## ----allTree.1-----------------------------------------------------------
head(TreesEx)

outTree <- calcTreeMets(TreesEx)
head(outTree)

## ----subTree.1-----------------------------------------------------------
outSnag <- calcSnagMets(TreesEx)
head(outSnag)
unique(outSnag$PARAMETER)

## ----allVG---------------------------------------------------------------
head(Vtype_GrCovEx)

outVG <- calcVtype_GcovMets(Vtype_GrCovEx)

head(outVG)

## ----subVG---------------------------------------------------------------
library(plyr)
nplots <- ddply(Vtype_GrCovEx,c('UID'),summarize,NPLOTS=length(unique(PLOT)))

outST <- calcSandTMets(Vtype_GrCovEx,nplots)
head(outST)
unique(outST$PARAMETER)

outVstrat <- calcVascStratMets(Vtype_GrCovEx,nplots)
head(outVstrat)
unique(outVstrat$PARAMETER)

# Now we can combine these two outputs and cast them into wide format
outAllVG <- rbind(outST,outVstrat)
outAllVG.wide <- dcast(outAllVG,UID~PARAMETER,value.var='RESULT')
head(outAllVG.wide)

