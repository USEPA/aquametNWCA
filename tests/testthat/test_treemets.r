library(aquametNWCA)
library(testthat)

context("Tree metric functions")


test_that("Snag metric values correct",
          {snagOut <- calcSnagMets(testTree,sampID='UID')
           compOut <- merge(testMets,snagOut,by=c('UID','PARAMETER'))
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
           expect_true(nrow(snagOut)==140)
           expect_true(nrow(compOut)==nrow(snagOut))
           expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.001)            
          })

test_that("Tree count metric values correct",
          {tcntOut <- calcTreeCntMets(testTree,sampID='UID')
           compOut <- merge(testMets,tcntOut,by=c('UID','PARAMETER'))
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
           expect_true(nrow(tcntOut)==160)
           expect_true(nrow(compOut)==nrow(tcntOut))
           expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.001)                      
          })

test_that("Tree cover metric values correct",
          {tcvrOut <- calcTreeCoverMets(testTree,sampID='UID')
           compOut <- merge(testMets,tcvrOut,by=c('UID','PARAMETER'))
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
           expect_true(nrow(tcvrOut)==400)
           expect_true(nrow(compOut)==nrow(tcvrOut))
           expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.001)            
           
          })

test_that("All tree metric values correct",
          {
            treeOut <- calcTreeMets(testTree,sampID='UID')
            treeOut.long <- reshape2::melt(treeOut,id.vars=c('UID'),variable.name='PARAMETER',value.name='RESULT')
            compOut <- merge(testMets,treeOut.long,by=c('UID','PARAMETER'))
            compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
            expect_true(nrow(treeOut.long)==700)
            expect_true(nrow(compOut)==nrow(treeOut.long))
            expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.001)                      
            
          })