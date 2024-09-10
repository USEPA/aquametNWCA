# NWCA_TestTreeMetricCode.r
# Purpose: Test tree metric code to make sure changes in code do not change 
# initial metric calculations beyond the test data
#
# Created 3/11/2020 by Karen Blocksom
#############################################################################

require(gtools) #smartbind
require(plyr) # mutate and ddply
require(car) # Recode
require(reshape2) # melt and dcast
library(aquametNWCA)
library(dplyr)

tree <- read.delim("L:/priv/CORFiles/IM-TH007/data/im/nwca/raw/tabfiles/tree.tab",header=TRUE,sep='\t',stringsAsFactors=F) %>%
  subset(select=c('UID','PLOT','PAGE','LINE','PARAMETER','RESULT'))

treePlots <- ddply(tree, c('UID'), summarise, NPLOTS=length(unique(PLOT)))

treeMets <- calcTreeMets(tree, treePlots, sampID='UID')

treeMets.long <- reshape2::melt(treeMets, id.vars=c('UID'), variable.name='PARAMETER', value.name='RESULT')

curMets <- read.delim("L:/priv/CORFiles/IM-TH007/data/im/nwca/data/tabfiles/vegMetrics.tab", sep='\t',
                          stringsAsFactors=F)

curMets.long <- select(curMets, -PUBLICATION_DATE,-DATE_COL,-SITE_ID,-SITE_USE,-STATE,-VISIT_NO) %>%
  reshape2::melt(id.var=c('UID'), variable.name='PARAMETER', value.name='RESULT')

matchMets <- merge(curMets.long, treeMets.long, by=c('UID','PARAMETER'), all.y=TRUE)

msgCur <- filter(matchMets, is.na(RESULT.x)) # These appear to be 50 Alaska sites, which we did not calculate metrics for
table(msgCur$UID)
table(msgCur$PARAMETER) # All tree metrics are missing in the metrics file

diffMets <- filter(matchMets, RESULT.x!=RESULT.y) # Zero mismatches


