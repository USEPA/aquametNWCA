library(aquametNWCA)
library(testthat)

context("Vascular plant metric functions")

# Create an input dataset and an expected output dataset

test_that("Data frames created with correct structure", {
  testDFs <- nwcaVegInput(
    sampID = "UID", tvar = "USDA_NAME", testVasc, taxon_name = "USDA_NAME",
    taxaNWCA, state = "STATE", coeReg = "USAC_REGION", cValReg = "STATE"
  )
  expect_true(length(testDFs) == 2)
  expect_true(class(testDFs) == "list")
  expect_equal(names(testDFs), c("byUID", "byPlot"))
  expect_equivalent(names(testDFs$byUID), c("UID", "STATE", "USAC_REGION", "TAXON", "NUM", "XABCOV", "NPLOTS", "DISTINCT"))
  expect_equivalent(names(testDFs$byPlot), c("UID", "STATE", "USAC_REGION", "TAXON", "PLOT", "COVER", "DISTINCT"))
})


dfTest <- prepareData(testVasc,
  sampID = "UID", inTaxa = taxaNWCA,
  inNat = ccNatNWCA,
  inCVal = ccNatNWCA, inWIS = wisNWCA,
  taxon_name = "USDA_NAME", state = "STATE",
  coeReg = "USAC_REGION", cValReg = "STATE"
)

test_that("Richness metric values are correct", {
  richOut <- calcRichness(dfTest[[1]], dfTest[[2]], dfTest[[3]], dfTest[[4]], dfTest[[5]], dfTest[[6]])
  compOut <- merge(testMets, richOut, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(compOut) == nrow(richOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.001)
})


# Create dataset for other metric calculation
# Use sXRCOV and sRFREQ instead so to avoid issues with calculations
testForCalc <- dfTest$byUIDspp


test_that("Category metric function values correct", {
  catOut <- calcCategory(testForCalc)
  compOut <- merge(testMets, catOut, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(catOut) == 720)
  expect_true(nrow(compOut) == nrow(catOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.01)
})

test_that("Duration metric function values correct", {
  durOut <- calcDuration(testForCalc)
  compOut <- merge(testMets, durOut, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(durOut) == 480)
  expect_true(nrow(compOut) == nrow(durOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.001)
})

test_that("Growth habit metric function values correct", {
  grhOut <- calcGrowthHabit(testForCalc)
  compOut <- merge(testMets, grhOut, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(grhOut) == 1160)
  expect_true(nrow(compOut) == nrow(grhOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.01)
})

test_that("WIS metric function values correct", {
  wisOut <- calcWIS(testForCalc)
  compOut <- merge(testMets, wisOut, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(wisOut) == 430)
  expect_true(nrow(compOut) == nrow(wisOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.001)
})

test_that("CC metric function values correct", {
  ccOut <- calcCC(testForCalc)
  compOut <- merge(testMets, ccOut, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(ccOut) == 320)
  expect_true(nrow(compOut) == nrow(ccOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.001)
})

test_that("Native metric function values correct", {
  natOut <- calcNative(testForCalc)
  compOut <- merge(testMets, natOut, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(natOut) == 300)
  expect_true(nrow(compOut) == nrow(natOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.01)
})

test_that("Diversity metric function values correct", {
  divOut <- calcDiversity(testForCalc)
  compOut <- merge(testMets, divOut, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(divOut) == 120)
  expect_true(nrow(compOut) == nrow(divOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.001)
})

test_that("Mean Bray-Curtis metric values correct", {
  bcOut <- calcBCmets(dfTest[[2]])
  compOut <- merge(testMets, bcOut, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(bcOut) == 20)
  expect_true(nrow(compOut) == nrow(bcOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.001)
})

test_that("VMMI metric calculations", {
  metOut <- calcVMMImets(testForCalc)
  varNames <- names(metOut)[names(metOut) != "UID"]
  metOut.long <- reshape(metOut,
    idvar = "UID", direction = "long",
    times = varNames, varying = varNames,
    timevar = "PARAMETER", v.names = "RESULT"
  )
  compOut <- merge(testMets, metOut.long, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(metOut.long) == 40)
  expect_true(nrow(compOut) == nrow(metOut.long))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.001)
})

test_that("All vascular plant metrics correct", {
  metOut <- calcVascPlantMets(testVasc,
    taxon_name = "USDA_NAME", taxaIn = taxaNWCA,
    taxaNat = ccNatNWCA, taxaCC = ccNatNWCA, taxaWIS = wisNWCA,
    sampID = "UID",
    state = "STATE", coeReg = "USAC_REGION", cValReg = "STATE"
  )
  varNames <- names(metOut)[names(metOut) != "UID"]
  metOut.long <- reshape(metOut,
    idvar = "UID", direction = "long",
    times = varNames, varying = varNames,
    timevar = "PARAMETER", v.names = "RESULT"
  )
  compOut <- merge(testMets, metOut.long, by = c("UID", "PARAMETER"))
  compOut$RESULT.x <- with(compOut, as.numeric(RESULT.x))
  expect_true(nrow(metOut.long) == 3920)
  expect_true(nrow(metOut.long) == nrow(compOut))
  expect_equal(compOut$RESULT.x, compOut$RESULT.y, tolerance = 0.001)
})
