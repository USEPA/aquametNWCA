#' @export
#'
#' @title Calculate bare ground and litter metrics
#'
#' @description Calculate bare ground and litter metrics based on
#' ground cover data. This function is called by calcVtype_GcovMets().
#'
#' @param dataIn A data frame containing the following variables:
#' \itemize{
#' \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#' \item PLOT (Sample plot from which data were collected)
#'
#' \item PARAMETER (specific measurement type)
#'
#' \item RESULT (measured value)
#' }
#' The following parameters are used in
#' calculating ground cover metrics: 'LITTER_THATCH', 'LITTER_FORB',
#' 'LITTER_CONIFER', 'LITTER_DECID', 'LITTER_BROADLEAF',
#' 'LITTER_DEPTH_SW', 'LITTER_DEPTH_NE', 'DEPTH_SW' (2016),
#' 'DEPTH_NE' (2016), TOTAL_LITTER', 'WD_FINE',
#' 'WD_COARSE', 'EXPOSED_SOIL', 'EXPOSED_GRAVEL', 'EXPOSED_ROCK',
#'  'PREDOMINANT_LITTER' (2016).
#'
#' Additional parameters or variables are ignored.
#' @param nPlot A data frame with the
#' number of plots sampled associated with each sample,
#' including \emph{sampID} variables and NPLOTS.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @param survyear A string with the survey year. The default is '2011'.
#' Starting in 2016, parameter structure changed such that, instead of
#' individual LITTER_XXXX parameters to indicate the preominant litter
#' type, PREDOMINANT_LITTER was used. In addition, DEPTH_NE and DEPTH_SW
#' replaced LITTER_DEPTH_NE and LITTER_DEPTH_SW.
#'
#' @details If any of the parameters are missing, they are assumed to be
#' zeros (if numeric), and metric values associated with any metrics that
#' cannot be calculated due to missing parameters are set to a
#' standardized value.
#'
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. Data frame contains
#'   \emph{sampID} variables, PARAMETER, RESULT, where values of PARAMETER
#'   consist of the metric name concatenated with trait value, with valid values
#'   of: FREQ_BAREGD, FREQ_EXPOSED_GRAVEL, FREQ_EXPOSED_ROCK, FREQ_EXPOSED_SOIL,
#'   FREQ_LITTER, FREQ_WD_COARSE, FREQ_WD_FINE, IMP_BAREGD, IMP_EXPOSED_GRAVEL,
#'   IMP_EXPOSED_ROCK, IMP_EXPOSED_SOIL, IMP_LITTER, IMP_WD_COARSE, IMP_WD_FINE,
#'   XCOV_BAREGD, XCOV_EXPOSED_GRAVEL, XCOV_EXPOSED_ROCK, XCOV_EXPOSED_SOIL,
#'   XCOV_LITTER, XCOV_WD_COARSE, XCOV_WD_FINE. A list of metric descriptions is
#'   provided in the document named \href{https://github.com/USEPA/aquametNWCA/blob/main/inst/VType_GrdCvr_Metric_Descriptions.pdf}{VegTypes_GrdCover_Metric_Descriptions.pdf}
#'   
#'   included in the help directory for the package.
#'
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#'
#' @examples
#' head(Vtype_GrCovEx)
#' # Create data frame with number of plots sampled for each UID
#' nplots <- data.frame(UID = seq(1:10), NPLOTS = rep(5, 10), stringsAsFactors = F)
#' # alternative approach to creating this data frame
#' nplots <- aggregate(
#'   x = list(NPLOTS = Vtype_GrCovEx$PLOT),
#'   by = Vtype_GrCovEx[c("UID")],
#'   FUN = function(x) {
#'     length(unique(x))
#'   }
#' )
#'
#' bgEx <- calcBareGround_LitterMets(Vtype_GrCovEx, nplots)
#'
#' head(bgEx)
#' unique(bgEx$PARAMETER)
calcBareGround_LitterMets <- function(dataIn, nPlot, sampID = "UID", survyear = "2011") {
  ## Now merge back with input df
  dataIn1 <- merge(dataIn, nPlot, by = sampID)

  # Create vector of all samples in dataset
  for (i in 1:length(sampID)) {
    if (i == 1) {
      dataIn1$SAMPID <- dataIn1[, sampID[i]]
    } else {
      dataIn1$SAMPID <- paste(dataIn1$SAMPID, dataIn1[, sampID[i]], sep = ".")
    }
  }
  samples <- unique(subset(dataIn1, select = c(sampID, "SAMPID")))

  allUIDs <- data.frame(SAMPID = unique(dataIn1$SAMPID), stringsAsFactors = F)

  # LITTER TYPES
  # Need to calculate the number of quadrats sampled using the NE and SW parameters
  # Only one set of parameters or the other will be in any dataset
  litter.sub <- subset(dataIn1, PARAMETER %in% c("LITTER_DEPTH_NE", "LITTER_DEPTH_SW", "DEPTH_NE", "DEPTH_SW"))
  # Calculate the number of quadrats (SW and NE)
  numQuads <- aggregate(
    x = list(NQUADS = litter.sub$RESULT),
    by = litter.sub[c("SAMPID")], FUN = length
  )
  # Subset data for litter and bare ground parameters
  litter1 <- subset(dataIn1, PARAMETER %in% c(
    "LITTER_THATCH", "LITTER_FORB", "LITTER_CONIFER",
    "LITTER_DECID", "LITTER_BROADLEAF",
    "LITTER_DEPTH_SW", "LITTER_DEPTH_NE",
    "TOTAL_LITTER", "WD_FINE", "WD_COARSE",
    "EXPOSED_SOIL", "EXPOSED_GRAVEL",
    "EXPOSED_ROCK", "DEPTH_SW", "DEPTH_NE",
    "PREDOMINANT_LITTER"
  ))
  # Merge litter subset with number of quadrats
  litter2 <- merge(litter1, numQuads, by = "SAMPID")
  # Subset to needed variables
  litter2 <- subset(litter2, select = c("SAMPID", "PLOT", "PARAMETER", "RESULT", "NPLOTS", "NQUADS"))
  # If survey year 2011, differs in parameter names and format of data from later years
  if (survyear == "2011") {
    littype <- subset(litter2, PARAMETER %in% c(
      "LITTER_THATCH", "LITTER_FORB", "LITTER_CONIFER", "LITTER_DECID",
      "LITTER_BROADLEAF", "LITTER_NONE"
    ) & !is.na(RESULT))
    # Count number of litter types
    rr1 <- aggregate(
      x = list(N_LITTER_TYPE = littype$PARAMETER),
      by = littype[c("SAMPID")],
      FUN = function(x) {
        length(unique(x))
      }
    )
    # Merge this count back to litter data
    rr1 <- merge(littype, rr1, by = "SAMPID")
    # Count plots by litter type
    rr2 <- aggregate(
      x = list(NUM = rr1$PLOT), by = rr1[c("SAMPID", "PARAMETER", "N_LITTER_TYPE")],
      FUN = length
    )
    # Find dominant litter type
    rr3 <- aggregate(x = list(MAXN = rr2$NUM), by = rr2[c("SAMPID")], FUN = max)
    rr3 <- merge(rr2, rr3, by = "SAMPID")
    rr4 <- subset(rr3, MAXN == NUM)
    # Identify max litter types and combine if more than one with same frequency
    rr5 <- aggregate(
      x = list(LITTER_TYPE = rr4$PARAMETER),
      by = rr4[c("SAMPID", "N_LITTER_TYPE")],
      FUN = function(x) {
        paste(x, collapse = "_")
      }
    )
    rr5$LITTER_TYPE <- with(rr5, gsub("LITTER_", "", LITTER_TYPE))
  } else {
    # Subset to predominant litter parameter where not missing
    littype <- subset(litter2, PARAMETER == "PREDOMINANT_LITTER" & !is.na(RESULT))
    # Count the total number of litter types
    rr1 <- aggregate(
      x = list(N_LITTER_TYPE = littype$RESULT),
      by = littype[c("SAMPID")],
      FUN = function(x) {
        length(unique(x))
      }
    )
    # Merge back with data
    rr1 <- merge(littype, rr1, by = "SAMPID")
    # Count number of each litter type
    rr2 <- aggregate(
      x = list(NUM = rr1$PLOT), by = rr1[c("SAMPID", "RESULT", "N_LITTER_TYPE")],
      FUN = length
    )
    # Find max litter type
    rr3 <- aggregate(x = list(MAXN = rr2$NUM), by = rr2[c("SAMPID")], FUN = max)
    rr3 <- merge(rr2, rr3, by = "SAMPID")
    rr4 <- subset(rr3, MAXN == NUM)
    # Combine litter types if more than one
    rr5 <- aggregate(
      x = list(LITTER_TYPE = rr4$RESULT),
      by = rr4[c("SAMPID", "N_LITTER_TYPE")],
      FUN = function(x) {
        paste(x, collapse = "_")
      }
    )
  }
  # To determine median depth, we must account for any quadrats without depth recorded
  litdep <- subset(litter2, PARAMETER %in% c("LITTER_DEPTH_SW", "LITTER_DEPTH_NE", "DEPTH_NE", "DEPTH_SW"))
  # Sum litter depth
  ss1.sum <- aggregate(
    x = list(XDEPTH_LITTER = litdep$RESULT),
    by = litdep[c("SAMPID", "NQUADS")],
    FUN = function(x) {
      sum(as.numeric(x))
    }
  )
  # Calculate number of depths
  ss1.length <- aggregate(
    x = list(NSAMP = litdep$RESULT),
    by = litdep[c("SAMPID", "NQUADS")],
    FUN = length
  )
  # Merge data frames
  ss1 <- merge(ss1.sum, ss1.length, by = c("SAMPID", "NQUADS"))
  # Use the merged values to calculate the mean litter depth
  ss1$XDEPTH_LITTER <- with(ss1, round(XDEPTH_LITTER / NQUADS, 2))
  # calculate number of samples with litter depth info
  litdep1 <- aggregate(
    x = list(NSAMP = litdep$RESULT),
    by = litdep[c("SAMPID", "NPLOTS", "NQUADS")],
    FUN = length
  )
  # Merge back with data
  litdep1 <- merge(litdep, litdep1, by = c("SAMPID", "NPLOTS", "NQUADS"))
  # Calculate median litter depth
  tt <- aggregate(
    x = list(MEDDEPTH_LITTER = litdep1$RESULT),
    by = litdep1[c("SAMPID")],
    FUN = function(x) {
      median(as.numeric(x))
    }
  )
  # Melt litter type metrics data frame
  varNames <- names(rr5)[!names(rr5) %in% "SAMPID"]
  rr5.long <- reshape(rr5,
    idvar = "SAMPID", direction = "long",
    times = varNames, varying = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )
  # Melt litter depth metrics, select subset of variables
  ss1.long <- reshape(ss1,
    idvar = "SAMPID", direction = "long",
    times = "XDEPTH_LITTER", varying = "XDEPTH_LITTER",
    timevar = "METRIC", v.names = "RESULT"
  )
  ss1.long <- subset(ss1.long, select = c("SAMPID", "METRIC", "RESULT"))

  # Melt median litter depth data frame
  varNames <- names(tt)[!names(tt) %in% "SAMPID"]
  tt.long <- reshape(tt,
    idvar = "SAMPID", direction = "long",
    times = varNames, varying = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )
  # Combine datasets into one data frame
  loutdf <- rbind(rr5.long, ss1.long, tt.long)
  loutdf$METRIC <- with(loutdf, as.character(METRIC))

  # Cast output data frame wide and drop prefix from variable names added by reshape()
  loutdf1 <- reshape(loutdf,
    idvar = "SAMPID", direction = "wide",
    timevar = "METRIC", v.names = "RESULT"
  )
  names(loutdf1) <- gsub("RESULT\\.", "", names(loutdf1))

  # Create empty data frame with all expected metrics
  empty_lit <- data.frame(t(rep(NA, 4)), stringsAsFactors = FALSE)
  names(empty_lit) <- c("N_LITTER_TYPE", "LITTER_TYPE", "XDEPTH_LITTER", "MEDDEPTH_LITTER")

  # Merge empty data frame with output data frame, drop row with missing SAMPID
  loutdf2 <- subset(merge(loutdf1, empty_lit, all = TRUE), !is.na(SAMPID))
  # Merge full sample list with output data frame
  loutdf3 <- merge(allUIDs, loutdf2, by = "SAMPID", all.x = T)
  # Melt data frame
  varNames <- names(loutdf3)[!names(loutdf3) %in% "SAMPID"]
  loutdf4 <- reshape(loutdf3,
    idvar = "SAMPID", direction = "long",
    times = varNames, varying = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )
  loutdf4$METRIC <- with(loutdf4, as.character(METRIC))

  print("Done with litter metrics")
  # Cast data frame wide and drop variable name prefix added by reshape()
  litterOut <- reshape(loutdf4,
    idvar = "SAMPID", direction = "wide",
    timevar = "METRIC", v.names = "RESULT"
  )
  names(litterOut) <- gsub("RESULT\\.", "", names(litterOut))

  # BARE GROUND
  # Subset data to only parameters related to woody debris and bare ground, total litter
  # Keep only positive RESULT values
  bgrd <- subset(dataIn1, PARAMETER %in% c(
    "EXPOSED_SOIL", "EXPOSED_GRAVEL", "EXPOSED_ROCK",
    "WD_FINE", "WD_COARSE", "TOTAL_LITTER"
  ) &
    RESULT != 0 & !is.na(RESULT), select = c("SAMPID", "PLOT", "PARAMETER", "RESULT", "NPLOTS"))
  # Need to create values for new PARAMETER='BAREGD' based on occurrence of either EXPOSED_SOIL,
  # EXPOSED_GRAVEL, or EXPOSED_ROCK at site
  bgrd.sub <- subset(bgrd, PARAMETER %in% c("EXPOSED_SOIL", "EXPOSED_GRAVEL", "EXPOSED_ROCK"))
  # Sum bare ground cover
  bgrd1 <- aggregate(
    x = list(RESULT = bgrd.sub$RESULT),
    by = bgrd.sub[c("SAMPID", "PLOT", "NPLOTS")],
    FUN = function(x) {
      as.character(sum(as.numeric(x)))
    }
  )
  bgrd1$PARAMETER <- "BAREGD"
  # Add these data back into full dataset
  bgrdIn <- rbind(bgrd, bgrd1)

  # Need to fill in zeros if plot sampled and variable is zero
  # Cast data frame, then drop prefix added to variable names by reshape()
  bgrdIn1.wide <- reshape(bgrdIn,
    idvar = c("SAMPID", "PLOT", "NPLOTS"), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(bgrdIn1.wide) <- gsub("RESULT\\.", "", names(bgrdIn1.wide))
  # Melt data frame, fill in missing RESULT values with 0
  varNames <- names(bgrdIn1.wide)[!names(bgrdIn1.wide) %in% c("SAMPID", "PLOT", "NPLOTS")]
  bgrdIn1 <- reshape(bgrdIn1.wide,
    idvar = c("SAMPID", "PLOT", "NPLOTS"), direction = "long",
    times = varNames, varying = varNames,
    timevar = "PARAMETER", v.names = "RESULT"
  )
  bgrdIn1$RESULT <- with(bgrdIn1, ifelse(is.na(RESULT), 0, as.numeric(RESULT)))
  # Keep only values from above where RESULT is >0
  bgrdIn1.sub <- subset(bgrdIn1, RESULT != 0)
  # Count the number of occurrences of each type of bare ground
  bgindf1.freq <- aggregate(
    x = list(FREQ = bgrdIn1.sub$RESULT),
    by = bgrdIn1.sub[c("SAMPID", "PARAMETER", "NPLOTS")],
    FUN = length
  )
  # Sum the cover of each type of bare ground
  bgindf1.sum <- aggregate(
    x = list(XCOV = bgrdIn1.sub$RESULT),
    by = bgrdIn1.sub[c("SAMPID", "PARAMETER", "NPLOTS")],
    FUN = function(x) {
      sum(as.numeric(x))
    }
  )
  # Merge these calculated values, then calculate the metrics based on them
  bgindf1 <- merge(bgindf1.sum, bgindf1.freq, by = c("SAMPID", "PARAMETER", "NPLOTS"), all = TRUE)
  bgindf1$FREQ <- with(bgindf1, round(FREQ / NPLOTS * 100, 2))
  bgindf1$XCOV <- with(bgindf1, round(XCOV / NPLOTS, 2))
  bgindf1$IMP <- with(bgindf1, round((FREQ + XCOV) / 2, 2))

  # Melt output data frame
  bgoutdf <- reshape(bgindf1,
    idvar = c("SAMPID", "PARAMETER"), direction = "long",
    times = c("XCOV", "FREQ", "IMP"), varying = c("XCOV", "FREQ", "IMP"),
    timevar = "variable", v.names = "RESULT"
  )
  # Create parameter name based on combinations of variable and PARAMETER
  bgoutdf$PARAMETER <- with(bgoutdf, ifelse(PARAMETER == "TOTAL_LITTER",
    paste(variable, "LITTER", sep = "_"),
    paste(variable, PARAMETER, sep = "_")
  ))
  bgoutdf <- bgoutdf[, c("SAMPID", "PARAMETER", "RESULT")]
  # Cast data frame and drop prefix to variable name added by reshape()
  bgoutdf1.wide <- reshape(bgoutdf,
    idvar = "SAMPID", direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(bgoutdf1.wide) <- gsub("RESULT\\.", "", names(bgoutdf1.wide))
  # Melt data frame
  varNames <- names(bgoutdf1.wide)[!names(bgoutdf1.wide) %in% c("SAMPID")]
  bgoutdf1 <- reshape(bgoutdf1.wide,
    idvar = c("SAMPID"), direction = "long",
    times = varNames, varying = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )
  # Fill in zeroes where missing
  bgoutdf1$RESULT <- with(bgoutdf1, ifelse(is.na(RESULT), 0, RESULT))
  bgoutdf1$METRIC <- with(bgoutdf1, as.character(METRIC))
  # Cast data wide and drop variable name prefix added by reshape()
  bgoutdf2 <- reshape(bgoutdf1,
    idvar = "SAMPID", direction = "wide",
    timevar = "METRIC", v.names = "RESULT"
  )
  names(bgoutdf2) <- gsub("RESULT\\.", "", names(bgoutdf2))
  # Create empty data frame
  empty_bg <- data.frame(t(rep(NA, 21)), stringsAsFactors = FALSE)
  names(empty_bg) <- c(
    "FREQ_BAREGD", "FREQ_EXPOSED_GRAVEL", "FREQ_EXPOSED_ROCK",
    "FREQ_EXPOSED_SOIL", "FREQ_LITTER", "FREQ_WD_COARSE",
    "FREQ_WD_FINE", "IMP_BAREGD", "IMP_EXPOSED_GRAVEL",
    "IMP_EXPOSED_ROCK", "IMP_EXPOSED_SOIL", "IMP_LITTER",
    "IMP_WD_COARSE", "IMP_WD_FINE", "XCOV_BAREGD",
    "XCOV_EXPOSED_GRAVEL", "XCOV_EXPOSED_ROCK", "XCOV_EXPOSED_SOIL",
    "XCOV_LITTER", "XCOV_WD_COARSE", "XCOV_WD_FINE"
  )
  # Merge empty data frame with output data frame and drop row with missing SAMPID
  bgoutdf3 <- subset(merge(bgoutdf2, empty_bg, all = TRUE), !is.na(SAMPID))
  # Merge full sample list with output data frame
  bgoutdf4 <- merge(allUIDs, bgoutdf3, by = "SAMPID", all.x = T)
  # Melt resulting data frame and fill in 0 for RESULT where missing
  varNames <- names(bgoutdf4)[!names(bgoutdf4) %in% c("SAMPID")]
  bgoutdf5 <- reshape(bgoutdf4,
    idvar = "SAMPID", direction = "long",
    times = varNames, varying = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )
  bgoutdf5$METRIC <- with(bgoutdf5, as.character(METRIC))
  bgoutdf5$RESULT <- with(bgoutdf5, ifelse(is.na(RESULT), 0, RESULT))

  print("Done with bare ground metrics")
  # Cast data frame wide and drop variable name prefix added by reshape()
  bgrdOut <- reshape(bgoutdf5,
    idvar = "SAMPID", direction = "wide",
    timevar = "METRIC", v.names = "RESULT"
  )
  names(bgrdOut) <- gsub("RESULT\\.", "", names(bgrdOut))

  # Now combine bare ground and litter metrics to perform a final check
  combOut <- merge(litterOut, bgrdOut, by = "SAMPID", all = T)

  # Melt data frame and fill in 0 where missing, except if metric is LITTER_TYPE
  varNames <- names(combOut)[!names(combOut) %in% c("SAMPID")]
  combOut.1 <- reshape(combOut,
    idvar = "SAMPID", direction = "long",
    times = varNames, varying = varNames,
    timevar = "variable", v.names = "value"
  )
  combOut.1$value <- with(combOut.1, ifelse(is.na(value) & variable != "LITTER_TYPE", 0, value))
  # Need to determine whether missing LITTER_TYPE value should be ABSENT or NONE
  # Cast output data frame and drop variable name prefix added by reshape()
  combOut.wide <- reshape(combOut.1,
    idvar = "SAMPID", direction = "wide",
    timevar = "variable", v.names = "value"
  )
  names(combOut.wide) <- gsub("value\\.", "", names(combOut.wide))
  # If LITTER_TYPE is missing and mean litter depth, mean litter cover, and litter frequency are all 0,
  # set LITTER_TYPE to 'ABSENT', otherwise set missing LITTER_TYPE to 'NONE'
  combOut.wide$LITTER_TYPE <- with(combOut.wide, ifelse(is.na(LITTER_TYPE) & XDEPTH_LITTER == 0 & XCOV_LITTER == 0 & FREQ_LITTER == 0, "ABSENT", ifelse(is.na(LITTER_TYPE), "NONE", LITTER_TYPE)))
  # Melt output data frame
  varNames <- names(combOut.wide)[!names(combOut.wide) %in% c("SAMPID")]
  combOut.long <- reshape(combOut.wide,
    idvar = "SAMPID", direction = "long",
    times = varNames, varying = varNames,
    timevar = "PARAMETER", v.names = "RESULT"
  )
  # Merge data with samples to bring back in sampID variables, then select only needed variables
  finOut <- merge(samples, combOut.long, by = "SAMPID")
  finOut <- finOut[c(sampID, "PARAMETER", "RESULT")]

  return(finOut)
}
