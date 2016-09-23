## ----allVasc.1-----------------------------------------------------------
library(aquametNWCA2)

head(VascPlantEx)

## ----allVasc.2-----------------------------------------------------------
outdf <- calcVascPlantMets(dfIn=VascPlantEx,taxaIn=taxaNWCA,taxaCC=ccNatNWCA,taxaWIS=wisNWCA,sampID='UID')

names(outdf)

## ----subVasc.1-----------------------------------------------------------
sumdf <- prepareData(indf=VascPlantEx, sampID='UID', inTaxa=taxaNWCA, inNatCC=ccNatNWCA, inWIS=wisNWCA)
str(sumdf)

## ----subVasc.2-----------------------------------------------------------
outRich <- calcRichness(sumdf$byUIDspp,sumdf$byPlotspp,sumdf$byUIDgen,sumdf$byPlotgen,sumdf$byUIDfam,sumdf$byPlotfam
                        ,sampID='UID')
head(outRich)
unique(outRich$PARAMETER)

## ----subVasc.3-----------------------------------------------------------
outCat <- calcCategory(sumdf$byUIDspp, sampID='UID')
head(outCat)
unique(outCat$PARAMETER) 

## ----subVasc.4-----------------------------------------------------------
outBC <- calcBCmets(sumdf$byPlotspp, sampID='UID')
head(outBC)
unique(outBC$PARAMETER)

## ----subVasc.5-----------------------------------------------------------
allVasc <- rbind(outRich, outCat, outBC)

library(reshape2)
outVasc <- dcast(allVasc,UID~PARAMETER,value.var='RESULT')
head(outVasc)

## ----allTree.1-----------------------------------------------------------
head(TreesEx)

outTree <- calcTreeMets(TreesEx, sampID='UID')
head(outTree)

## ----subTree.1-----------------------------------------------------------
outSnag <- calcSnagMets(TreesEx, sampID='UID')
head(outSnag)
unique(outSnag$PARAMETER)

## ----allVG---------------------------------------------------------------
head(Vtype_GrCovEx)

outVG <- calcVtype_GcovMets(Vtype_GrCovEx, sampID='UID')

head(outVG)

## ----subVG---------------------------------------------------------------
library(plyr)
nplots <- ddply(Vtype_GrCovEx,c('UID'),summarize,NPLOTS=length(unique(PLOT)))

outST <- calcSandTMets(Vtype_GrCovEx, nplots, sampID='UID')
head(outST)
unique(outST$PARAMETER)

outVstrat <- calcVascStratMets(Vtype_GrCovEx, nplots, sampID='UID')
head(outVstrat)
unique(outVstrat$PARAMETER)

# Now we can combine these two outputs and cast them into wide format
outAllVG <- rbind(outST,outVstrat)
outAllVG.wide <- dcast(outAllVG,UID~PARAMETER,value.var='RESULT')
head(outAllVG.wide)

