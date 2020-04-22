# NWCA_TestVegType-GrndCover_MetricCode.r
# Purpose: Test veg type and ground cover metric code to make sure changes in code do not change 
# initial metric calculations beyond the test data
#
# Created 4/22/2020 by Karen Blocksom
#############################################################################

require(gtools) #smartbind
require(plyr) # mutate and ddply
library(aquametNWCA)
library(dplyr)

vegtype <- read.delim("L:/priv/CORFiles/IM-TH007/data/im/nwca/raw/tabfiles/vegtype.tab",header=TRUE,sep='\t',stringsAsFactors=F) %>%
  subset(select=c('UID','PLOT','PARAMETER','RESULT'))

vegPlots <- filter(vegtype, PARAMETER!='SANDT_CLASS') %>%
  ddply(c('UID'), summarise, NPLOTS=length(unique(PLOT)))

vtypeMets <- calcVtype_GcovMets(vegtype, vegPlots, sampID='UID')

vtypeMets.long <- tidyr::pivot_longer(vtypeMets, cols=names(vtypeMets)[names(vtypeMets)!='UID'], 
                                      names_to='PARAMETER', values_to='RESULT')

curMets <- read.delim("L:/priv/CORFiles/IM-TH007/data/im/nwca/data/tabfiles/vegMetrics.tab", sep='\t',
                      stringsAsFactors=F)

curMets.long <- select(curMets, -PUBLICATION_DATE,-DATE_COL,-SITE_ID,-SITE_USE,-STATE,-VISIT_NO) %>%
  reshape2::melt(id.var=c('UID'), variable.name='PARAMETER', value.name='RESULT')

matchVtype <- merge(curMets.long, vtypeMets.long, by=c('UID','PARAMETER'))

diffVtype <- filter(matchVtype, RESULT.x!=RESULT.y) # No differences in values
