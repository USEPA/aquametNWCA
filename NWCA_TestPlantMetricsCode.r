# NWCA_TestPlantMetricCode.r
# Purpose: Test plant metric code to make sure changes in code do not change 
# initial metric calculations beyond the test data
#
# Created 4/6/2020 by Karen Blocksom
#############################################################################

require(gtools) #smartbind
require(plyr) # mutate and ddply
require(car) # Recode
require(reshape2) # melt and dcast
library(aquametNWCA)
library(dplyr)

plants <- read.delim("L:/priv/CORFiles/IM-TH007/data/im/nwca/raw/tabfiles/plant_wide.tab",header=TRUE,sep='\t',stringsAsFactors=F) %>%
  subset(!is.na(COVER), select=c('UID','STATE','PLOT','PAGE','LINE','SPECIES_NAME_ID','SPECIES','COVER')) %>%
  mutate(SPECIES=gsub("\u00D7","x", SPECIES)) # Replace the times symbol x with a lower case x

sites <- read.delim("L:/priv/CORFiles/IM-TH007/data/im/nwca/data/tabfiles/siteinfo_8-13-2015.tab",header=TRUE,sep='\t',stringsAsFactors=F) %>%
  select(UID, USAC_REGION)

plants.1 <- merge(sites, plants, by='UID') %>%
  plyr::rename(c('SPECIES'='USDA_NAME'))

plantPlots <- ddply(plants, c('UID'), summarise, NPLOTS=length(unique(PLOT)))

plantMets <- calcVascPlantMets(plants.1, sampID='UID')

plantMets.long <- reshape2::melt(plantMets, id.vars=c('UID'), variable.name='PARAMETER', value.name='RESULT') %>%
  mutate(PARAMETER=as.character(PARAMETER))

curMets <- read.delim("L:/priv/CORFiles/IM-TH007/data/im/nwca/data/tabfiles/vegmetrics.tab",header=TRUE,sep='\t',stringsAsFactors=F)

curMets.long <- select(curMets, -PUBLICATION_DATE, -DATE_COL, -SITE_ID, -SITE_USE, -VISIT_NO, -STATE) %>%
  reshape2::melt(id.vars=c('UID'), variable.name='PARAMETER', value.name='RESULT') %>%
  mutate(PARAMETER=as.character(PARAMETER))

matchMets <- merge(curMets.long, plantMets.long, by=c('UID','PARAMETER')) %>%
  mutate(RESULT.x=as.numeric(RESULT.x))

diffMets <- dplyr::filter(matchMets, abs(RESULT.x-RESULT.y)>0.001)
