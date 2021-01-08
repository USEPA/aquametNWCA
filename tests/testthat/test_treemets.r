library(aquametNWCA)
library(testthat)

context("Tree metric functions")

nplots <- plyr::ddply(testTree,c('UID'),dplyr::summarise
                                ,NPLOTS=length(unique(PLOT)))

test_that("Snag metric values correct",
          {
           snagOut <- calcSnagMets(testTree,nplots,sampID='UID')
           compOut <- merge(testMets,snagOut,by=c('UID','PARAMETER'))
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
           expect_true(nrow(snagOut)==140)
           expect_true(nrow(compOut)==nrow(snagOut))
           expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)            
          })

test_that("Tree count metric values correct",
          {
           tcntOut <- calcTreeCntMets(testTree,nplots,sampID='UID')
           compOut <- merge(testMets,tcntOut,by=c('UID','PARAMETER'))
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
           expect_true(nrow(tcntOut)==220)
           expect_true(nrow(compOut)==nrow(tcntOut))
           expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)                      
          })

test_that("Tree cover metric values correct",
          {
           tcvrOut <- calcTreeCoverMets(testTree,nplots,sampID='UID')
           compOut <- merge(testMets,tcvrOut,by=c('UID','PARAMETER'))
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
           expect_true(nrow(tcvrOut)==400)
           expect_true(nrow(compOut)==nrow(tcvrOut))
           expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)            
           
          })

test_that("All tree metric values correct",
          {
            treeOut <- calcTreeMets(testTree,nplots,sampID='UID')
            treeOut.long <- reshape2::melt(treeOut,id.vars=c('UID'),variable.name='PARAMETER',value.name='RESULT')
            compOut <- merge(testMets,treeOut.long,by=c('UID','PARAMETER'))
            compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
            expect_true(nrow(treeOut.long)==760)
            expect_true(nrow(compOut)==nrow(treeOut.long))
            expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)                      
            
          })
