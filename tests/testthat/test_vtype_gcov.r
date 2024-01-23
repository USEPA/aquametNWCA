library(aquametNWCA)
library(testthat)

context("Vegetation types and ground cover metric functions")

testVCGC.sub <- subset(testVCGC, PARAMETER != 'SANDT_CLASS')
nplots <- aggregate(x = list(NPLOTS = testVCGC.sub$PLOT), 
                    by = testVCGC.sub[c("UID")],
                    FUN = function(x){length(unique(x))})


test_that("S & T metric values correct",
        {
                 sandtOut <- calcSandTMets(testVCGC,nplots)
                 compOut <- merge(testMets,sandtOut,by=c('UID','PARAMETER'))
                 expect_true(nrow(sandtOut)==50)
                 expect_true(nrow(compOut)==nrow(sandtOut))
                 compOut.char <- subset(compOut,PARAMETER=='DOM_SANDT')
                 expect_equal(compOut.char$RESULT.x,compOut.char$RESULT.y)
                 compOut.num <- subset(compOut,PARAMETER!='DOM_SANDT')
                 expect_equal(compOut.num$RESULT.x,compOut.num$RESULT.y,tolerance=0.001)            
        })

test_that("Vascular stratum metric values correct",
        {
         vstratOut <- calcVascStratMets(testVCGC,nplots)
         compOut <- merge(testMets,vstratOut,by=c('UID','PARAMETER'))
         expect_true(nrow(vstratOut)==430)
         expect_true(nrow(compOut)==nrow(vstratOut))
         expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.01)            
        })

test_that("Non-vascular metric values correct",
        {
         nvascOut <- calcNonvascMets(testVCGC,nplots)
         compOut <- merge(testMets,nvascOut,by=c('UID','PARAMETER'))
         expect_true(nrow(nvascOut)==170)
         expect_true(nrow(compOut)==nrow(nvascOut))
         expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)            
        })

test_that("Water cover metric values correct",
        {
         wcovOut <- calcWcovMets(testVCGC,nplots)
         compOut <- merge(testMets,wcovOut,by=c('UID','PARAMETER'))
         expect_true(nrow(wcovOut)==180)
         expect_true(nrow(compOut)==nrow(wcovOut))
         expect_equal(as.numeric(compOut$RESULT.x),as.numeric(compOut$RESULT.y),tolerance=0.001)  
        })

test_that("Bare ground and litter metric values correct",
        {
         bglitOut <- calcBareGround_LitterMets(testVCGC,nplots)
         compOut <- merge(testMets,bglitOut,by=c('UID','PARAMETER'))
         expect_true(nrow(bglitOut)==250)
         expect_true(nrow(compOut)==nrow(bglitOut))
         compOut.char <- subset(compOut,PARAMETER=='LITTER_TYPE')
         expect_equal(compOut.char$RESULT.x,compOut.char$RESULT.y)
         compOut.num <- subset(compOut,PARAMETER!='LITTER_TYPE')
         expect_equal(as.numeric(compOut.num$RESULT.x),as.numeric(compOut.num$RESULT.y),tolerance=0.01)            
        })

test_that("Veg strata and ground cover metric values correct",
          {
                  vcgcOut <- calcVtype_GcovMets(testVCGC, nplots, sampID='UID', survyear='2011')
                  varNames <- names(vcgcOut)[names(vcgcOut)!='UID']
                  vcgcOut.long <- reshape(vcgcOut, idvar = 'UID', direction = 'long',
                                         times = varNames, varying = varNames,
                                         timevar = 'PARAMETER', v.names = 'RESULT')
                  compOut <- merge(testMets,vcgcOut.long,by=c('UID','PARAMETER'))
                  expect_true(nrow(vcgcOut.long)==1080)
                  expect_true(nrow(compOut)==nrow(vcgcOut.long))
                  compOut.char <- subset(compOut,PARAMETER %in% c('DOM_SANDT','LITTER_TYPE'))
                  expect_equal(compOut.char$RESULT.x,compOut.char$RESULT.y)
                  compOut.num <- subset(compOut,!PARAMETER %in% c('DOM_SANDT','LITTER_TYPE'))
                  expect_equal(as.numeric(compOut.num$RESULT.x),as.numeric(compOut.num$RESULT.y),
                               tolerance=0.01)            
  })
