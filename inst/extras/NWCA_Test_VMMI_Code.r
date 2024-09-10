# NWCA_Test_VMMI_Code.r

require(gtools) #smartbind
require(plyr) # mutate and ddply
require(car) # Recode
require(reshape2) # melt and dcast
library(aquametNWCA)
library(dplyr)
library(Hmisc)
library(tidyr)

vmets <- read.delim("O:/priv/cphea/pesd/cor/CORFiles/IM-TH007/data/im/nwca16/data/tabfiles/nwca16_Vegmetrics.tab",header=TRUE,sep='\t',stringsAsFactors=F) 

sites <- read.delim("O:/priv/cphea/pesd/cor/CORFiles/IM-TH007/data/im/nwca16/data/tabfiles/nwca16_siteinfo.tab",header=TRUE,sep='\t',stringsAsFactors=F)

vmets.1 <- merge(sites[,c('UID','WETCLS_GRP','RPT_UNIT')], vmets, by='UID')

vmets.eh_ew <- subset(vmets.1, WETCLS_GRP %in% c('EH','EW'))


test.eh_ew <- calcVMMI_2016_fromMets(vmets.eh_ew, sampID='UID', wetcls_grp='WETCLS_GRP')

test.eh_ew.long <- mutate_all(test.eh_ew, as.character) %>%
  pivot_longer(cols=N_ANNUAL_SC16:VEGCOND_2016, names_to='PARAMETER', values_to='RESULT',
               values_drop_na = TRUE)

curMMI <- read.delim("O:/priv/cphea/pesd/cor/CORFiles/IM-TH007/data/im/nwca16/data/tabfiles/nwca16_VMMI_condition.tab",header=TRUE,sep='\t',stringsAsFactors=F) 

curMMI.num <- pivot_longer(curMMI, cols = VMMI_2016:XRCOV_GRAMINOID_SC16, 
                            names_to='PARAMETER', values_to='RESULT', values_drop_na = TRUE) 

curMMI.cond <- pivot_longer(curMMI, cols = VEGCOND_2016, 
                            names_to='PARAMETER', values_to='RESULT', values_drop_na = TRUE) 

matchMMI.num <- merge(curMMI.num, test.eh_ew.long, by = c('UID', 'PARAMETER')) %>%
  mutate(RESULT.y = as.numeric(RESULT.y)) %>%
  mutate(RESULT.x = ifelse(PARAMETER=='VMMI_2016', round(RESULT.x, 1), round(RESULT.x, 2)))

filter(matchMMI.num, abs(RESULT.x-RESULT.y)>0.001)

matchMMI.cond <- merge(curMMI.cond, test.eh_ew.long, by = c('UID', 'PARAMETER')) 

filter(matchMMI.cond, RESULT.x!=RESULT.y)
