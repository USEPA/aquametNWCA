library(aquametNWCA)
library(testthat)

context("VMMI and condition values correct")

test_that("VMMI and condition using ECO_X_WETGRP correct", {
  testIn <- merge(testSites, testMets, by = "UID")
  testIn <- subset(testIn, PARAMETER %in% c(
    "FQAI_ALL", "N_TOL", "RIMP_NATSPP",
    "XRCOV_MONOCOTS_NAT"
  ),
  select = c("UID", "ECO_X_WETGRP", "PARAMETER", "RESULT")
  )

  # dplyr::select(UID,ECO_X_WETGRP,PARAMETER,RESULT) %>%
  # dplyr::filter(PARAMETER %in% c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')) %>%
  testIn.wide <- reshape(testIn,
    idvar = c("UID", "ECO_X_WETGRP"),
    direction = "wide", timevar = "PARAMETER",
    v.names = "RESULT"
  )
  names(testIn.wide) <- gsub("RESULT\\.", "", names(testIn.wide))
  # reshape2::dcast(UID+ECO_X_WETGRP~PARAMETER,value.var='RESULT')
  testOut <- calcVMMI_fromMets(testIn.wide)
  varNames <- names(testOut)[!names(testOut) %in% c("UID", "ECO_X_WETGRP")]
  testOut.long <- reshape(testOut,
    idvar = c("UID", "ECO_X_WETGRP"), direction = "long",
    times = varNames, varying = varNames,
    timevar = "PARAMETER", v.names = "RESULT"
  )
  testOut.long$PARAMETER <- as.character(testOut.long$PARAMETER)
  # testOut.long <- reshape2::melt(testOut,id.vars=c('UID','ECO_X_WETGRP'),variable.name='PARAMETER'
  #                                ,value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(testMets, testOut.long, by = c("UID", "PARAMETER"))
  expect_true(nrow(compOut) == 60)
  compOut.char <- subset(compOut, PARAMETER == "VEGCOND")
  expect_equal(compOut.char$RESULT.x, compOut.char$RESULT.y)
  compOut.num <- subset(compOut, PARAMETER != "VEGCOND")
  expect_equal(as.numeric(compOut.num$RESULT.x), as.numeric(compOut.num$RESULT.y), tolerance = 0.001)
})

test_that("VMMI and condition using NWCA_ECO4 and NWCA_WET_GRP correct", {
  testIn <- merge(testSites, testMets, by = "UID")
  testIn <- subset(testIn, PARAMETER %in% c(
    "FQAI_ALL", "N_TOL", "RIMP_NATSPP",
    "XRCOV_MONOCOTS_NAT"
  ),
  select = c("UID", "NWCA_ECO4", "NWCA_WET_GRP", "PARAMETER", "RESULT")
  )

  # testIn <- merge(testSites,testMets,by='UID') %>%
  #  dplyr::select(UID,NWCA_ECO4,NWCA_WET_GRP,PARAMETER,RESULT) %>%
  #  dplyr::filter(PARAMETER %in% c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')) %>%
  #  reshape2::dcast(UID+NWCA_ECO4+NWCA_WET_GRP~PARAMETER,value.var='RESULT')
  testIn.wide <- reshape(testIn,
    idvar = c("UID", "NWCA_ECO4", "NWCA_WET_GRP"),
    direction = "wide", timevar = "PARAMETER",
    v.names = "RESULT"
  )
  names(testIn.wide) <- gsub("RESULT\\.", "", names(testIn.wide))

  testOut <- calcVMMI_fromMets(testIn.wide)
  varNames <- names(testOut)[!names(testOut) %in% c("UID", "ECO_X_WETGRP")]
  testOut.long <- reshape(testOut,
    idvar = c("UID", "ECO_X_WETGRP"), direction = "long",
    times = varNames, varying = varNames,
    timevar = "PARAMETER", v.names = "RESULT"
  )
  testOut.long$PARAMETER <- as.character(testOut.long$PARAMETER)
  # testOut.long <- reshape2::melt(testOut,id.vars=c('UID','ECO_X_WETGRP'),variable.name='PARAMETER'
  #                                ,value.name='RESULT') %>%
  #   plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(testMets, testOut.long, by = c("UID", "PARAMETER"))
  expect_true(nrow(compOut) == 60)
  compOut.char <- subset(compOut, PARAMETER == "VEGCOND")
  expect_equal(compOut.char$RESULT.x, compOut.char$RESULT.y)
  compOut.num <- subset(compOut, PARAMETER != "VEGCOND")
  expect_equal(as.numeric(compOut.num$RESULT.x), as.numeric(compOut.num$RESULT.y), tolerance = 0.001)
})

test_that("VMMI lacking ecoregion and wetland type correct", {
  testIn <- merge(testSites, testMets, by = "UID")
  testIn <- subset(testIn, PARAMETER %in% c(
    "FQAI_ALL", "N_TOL", "RIMP_NATSPP",
    "XRCOV_MONOCOTS_NAT"
  ),
  select = c("UID", "PARAMETER", "RESULT")
  )
  # testIn <- merge(testSites,testMets,by='UID') %>%
  #  dplyr::select(UID,PARAMETER,RESULT) %>%
  #  dplyr::filter(PARAMETER %in% c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')) %>%
  #  reshape2::dcast(UID~PARAMETER,value.var='RESULT')
  testIn.wide <- reshape(testIn,
    idvar = c("UID"),
    direction = "wide", timevar = "PARAMETER",
    v.names = "RESULT"
  )
  names(testIn.wide) <- gsub("RESULT\\.", "", names(testIn.wide))

  testOut <- calcVMMI_fromMets(testIn.wide)

  # testOut <- calcVMMI_fromMets(testIn)
  varNames <- names(testOut)[!names(testOut) %in% c("UID", "ECO_X_WETGRP")]
  testOut.long <- reshape(testOut,
    idvar = c("UID"), direction = "long",
    times = varNames, varying = varNames,
    timevar = "PARAMETER", v.names = "RESULT"
  )
  testOut.long$PARAMETER <- as.character(testOut.long$PARAMETER)

  # testOut.long <- reshape2::melt(testOut,id.vars=c('UID'),variable.name='PARAMETER'
  #                                ,value.name='RESULT') %>%
  # plyr::mutate(PARAMETER=as.character(PARAMETER))
  compOut <- merge(testMets, testOut.long, by = c("UID", "PARAMETER"))
  expect_true(nrow(compOut) == 50)
  compOut <- subset(compOut, PARAMETER != "VEGCOND")
  expect_equal(as.numeric(compOut$RESULT.x), as.numeric(compOut$RESULT.y), tolerance = 0.001)
})
