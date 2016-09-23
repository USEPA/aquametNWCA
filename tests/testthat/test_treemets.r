library(aquametNWCA2)
library(testthat)

context("Tree metric functions")


test_that("Snag metric values correct",
          {snagOut <- calcSnagMets(testTree)
           compOut <- merge(testMets,snagOut,by=c('UID','PARAMETER'))
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
           expect_true(nrow(snagOut)==140)
           expect_true(nrow(compOut)==nrow(snagOut))
           expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.001)            
          })

test_that("Tree count metric values correct",
          {tcntOut <- calcTreeCntMets(testTree)
           compOut <- merge(testMets,tcntOut,by=c('UID','PARAMETER'))
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
           expect_true(nrow(tcntOut)==160)
           expect_true(nrow(compOut)==nrow(tcntOut))
           expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.001)                      
          })

test_that("Tree cover metric values correct",
          {tcvrOut <- calcTreeCoverMets(testTree)
           compOut <- merge(testMets,tcvrOut,by=c('UID','PARAMETER'))
           compOut <- dplyr::mutate(compOut,RESULT.x=as.numeric(RESULT.x))
           expect_true(nrow(tcvrOut)==400)
           expect_true(nrow(compOut)==nrow(tcvrOut))
           expect_equal(compOut$RESULT.x,compOut$RESULT.y,tolerance=0.001)            
           
          })