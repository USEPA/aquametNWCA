library(aquametNWCA)
library(testthat)

context("Tree metric functions")

nplots <- aggregate(x = list(NPLOTS = testTree$PLOT), 
                    by = testTree[c("UID")],
                    FUN = function(x){length(unique(x))})

test_that("Snag metric values correct",
          {
           snagOut <- calcSnagMets(testTree,nplots,sampID='UID')
           compOut <- merge(testMets,snagOut,by=c('UID','PARAMETER'))
           compOut$RESULT.x <- with(compOut,as.numeric(RESULT.x))
           expect_true(nrow(snagOut)==140)
           expect_true(nrow(compOut)==nrow(snagOut))
           expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)            
          })

test_that("Tree count metric values correct",
          {
           tcntOut <- calcTreeCntMets(testTree,nplots,sampID='UID')
           compOut <- merge(testMets,tcntOut,by=c('UID','PARAMETER'))
           compOut$RESULT.x <- with(compOut,as.numeric(RESULT.x))
           expect_true(nrow(tcntOut)==220)
           expect_true(nrow(compOut)==nrow(tcntOut))
           expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)                      
          })

test_that("Tree cover metric values correct",
          {
           tcvrOut <- calcTreeCoverMets(testTree,nplots,sampID='UID')
           compOut <- merge(testMets,tcvrOut,by=c('UID','PARAMETER'))
           compOut$RESULT.x <- with(compOut,as.numeric(RESULT.x))
           expect_true(nrow(tcvrOut)==400)
           expect_true(nrow(compOut)==nrow(tcvrOut))
           expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)            
           
          })

test_that("All tree metric values correct",
          {
            treeOut <- calcTreeMets(testTree,nplots,sampID='UID')
            varNames <- names(treeOut)[names(treeOut)!='UID']
            treeOut.long <- reshape(treeOut, idvar = 'UID', direction = 'long',
                                    times = varNames, varying = varNames,
                                    timevar = 'PARAMETER', v.names = 'RESULT')
            compOut <- merge(testMets,treeOut.long,by=c('UID','PARAMETER'))
            compOut$RESULT.x <- with(compOut,as.numeric(RESULT.x))
            expect_true(nrow(treeOut.long)==760)
            expect_true(nrow(compOut)==nrow(treeOut.long))
            expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)                      
            
          })

