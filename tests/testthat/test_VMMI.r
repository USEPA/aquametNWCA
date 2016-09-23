library(aquametNWCA2)
library(testthat)

context("VMMI and condition values correct")

test_that("VMMI and condition using ECO_X_WETGRP correct",
          {testIn <- merge(testSites,testMets,by='UID') %>%
            dplyr::select(UID,ECO_X_WETGRP,PARAMETER,RESULT) %>%
            dplyr::filter(PARAMETER %in% c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')) %>%
            reshape2::dcast(UID+ECO_X_WETGRP~PARAMETER,value.var='RESULT')
          testOut <- calcVMMI_fromMets(testIn)
          testOut.long <- reshape2::melt(testOut,id.vars=c('UID','ECO_X_WETGRP'),variable.name='PARAMETER'
                                         ,value.name='RESULT') %>%
            plyr::mutate(PARAMETER=as.character(PARAMETER))
          compOut <- merge(testMets,testOut.long,by=c('UID','PARAMETER'))
          expect_true(nrow(compOut)==60)
          compOut.char <- subset(compOut,PARAMETER=='VEGCOND')
          expect_equal(compOut.char$RESULT.x,compOut.char$RESULT.y)
          compOut.num <- subset(compOut,PARAMETER!='VEGCOND')
          compOut.num <- dplyr::mutate(compOut.num,RESULT.x=as.numeric(RESULT.x),RESULT.y=as.numeric(RESULT.y))
          expect_equal(compOut.num$RESULT.x,compOut.num$RESULT.y,tolerance=0.001)    
          }
)

test_that("VMMI and condition using NWCA_ECO4 and NWCA_WET_GRP correct",
          {testIn <- merge(testSites,testMets,by='UID') %>%
             dplyr::select(UID,NWCA_ECO4,NWCA_WET_GRP,PARAMETER,RESULT) %>%
             dplyr::filter(PARAMETER %in% c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')) %>%
             reshape2::dcast(UID+NWCA_ECO4+NWCA_WET_GRP~PARAMETER,value.var='RESULT')
           testOut <- calcVMMI_fromMets(testIn)
           testOut.long <- reshape2::melt(testOut,id.vars=c('UID','ECO_X_WETGRP'),variable.name='PARAMETER'
                                          ,value.name='RESULT') %>%
             plyr::mutate(PARAMETER=as.character(PARAMETER))
           compOut <- merge(testMets,testOut.long,by=c('UID','PARAMETER'))
           expect_true(nrow(compOut)==60)
           compOut.char <- subset(compOut,PARAMETER=='VEGCOND')
           expect_equal(compOut.char$RESULT.x,compOut.char$RESULT.y)
           compOut.num <- subset(compOut,PARAMETER!='VEGCOND')
           compOut.num <- dplyr::mutate(compOut.num,RESULT.x=as.numeric(RESULT.x),RESULT.y=as.numeric(RESULT.y))
           expect_equal(compOut.num$RESULT.x,compOut.num$RESULT.y,tolerance=0.001)                
          })

test_that("VMMI lacking ecoregion and wetland type correct",
          {testIn <- merge(testSites,testMets,by='UID') %>%
             dplyr::select(UID,PARAMETER,RESULT) %>%
             dplyr::filter(PARAMETER %in% c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')) %>%
             reshape2::dcast(UID~PARAMETER,value.var='RESULT')
           testOut <- calcVMMI_fromMets(testIn)
           testOut.long <- reshape2::melt(testOut,id.vars=c('UID'),variable.name='PARAMETER'
                                          ,value.name='RESULT') %>%
             plyr::mutate(PARAMETER=as.character(PARAMETER))
           compOut <- merge(testMets,testOut.long,by=c('UID','PARAMETER'))
           expect_true(nrow(compOut)==50)
           compOut <- subset(compOut,PARAMETER!='VEGCOND')
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x),RESULT.y=as.numeric(RESULT.y))
           expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.001)                
          })