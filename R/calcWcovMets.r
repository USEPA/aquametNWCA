#' @export
#'
#' @title Calculate water cover metrics
#'
#' @description This function calculates metrics water
#' cover and depth based on ground cover data. This function is
#' called by calcVtype_GcovMets().
#'
#' @param dataIn A data frame containing the following variables:
#' \itemize{
#' \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#' \item PLOT: Sample plot from which data
#' were collected
#'
#' \item PARAMETER: specific measurement type
#'
#' \item RESULT: measured value
#' }
#' The following parameters are used in
#' calculating water cover metrics: 'TIME', 'MINIMUM_DEPTH',
#' 'MAXIMUM_DEPTH', 'PREDOMINANT_DEPTH', 'TOTAL_WATER',
#' 'WATER_NOVEG', 'WATER_AQVEG', 'WATER_EMERGVEG'.
#' Additional parameters or variables are ignored.
#' WATER_NOVEG, WATER_AQVEG, MINIMUM_DEPTH, and
#' MAXIMUM_DEPTH not measured in 2016.
#' @param nPlot A data frame with the
#' number of plots sampled associated with each
#' sample with \emph{sampID} variables and NPLOTS
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples,
#' 'UID' by default.
#' @param survyear A string with the survey year. The default is '2011'.
#' Starting in 2016, only TIME, TOTAL_WATER, and PREDOMINANT_DEPTH is
#' measured.
#'
#' @details If any of the parameters are missing, they are assumed
#' to be zeros (if numeric), and metric values associated with any
#' metrics that cannot be calculated due to missing parameters are set
#' to a standardized value.
#'
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. Data frame contains
#'   \emph{sampID} variables, PARAMETER, RESULT, where values of PARAMETER
#'   consist of the metric name concatenated with trait value, with valid values
#'   of: MIN_H2O_DEPTH, MAX_H2O_DEPTH, XH2O_DEPTH_AA, MIN_COV_H2O, MAX_COV_H2O,
#'   FREQ_H2O, FREQ_H2O_AQVEG, FREQ_H2O_EMERGVEG, FREQ_H2O_NOVEG, XCOV_H2O,
#'   XCOV_H2O_AQVEG, XCOV_H2O_EMERGVEG, XCOV_H2O_NOVEG, IMP_H2O, IMP_H2O_AQVEG,
#'   IMP_H2O_EMERGVEG, IMP_H2O_NOVEG, XH2O_DEPTH. A list of metric descriptions
#'   is provided in the document named \href{https://github.com/USEPA/aquametNWCA/blob/main/inst/VType_GrdCvr_Metric_Descriptions.pdf}{VegTypes_GrdCover_Metric_Descriptions.pdf}
#'   included in the help directory for the package.
#'
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#'
#' @examples
#' head(Vtype_GrCovEx)
#' # Create data frame with number of plots sampled for each sampID
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
#' wcovEx <- calcWcovMets(Vtype_GrCovEx, nplots)
#'
#' head(wcovEx)
#' unique(wcovEx$METRIC)
calcWcovMets <- function(dataIn, nPlot, sampID = "UID", survyear = "2011") {
  # Merge nPlot with dataIn by sampID variables
  dataIn1 <- merge(dataIn, nPlot, by = sampID)

  # Create SAMPID variable based on combination of variables in sampID
  for (i in 1:length(sampID)) {
    if (i == 1) {
      dataIn1$SAMPID <- dataIn1[, sampID[i]]
    } else {
      dataIn1$SAMPID <- paste(dataIn1$SAMPID, dataIn1[, sampID[i]], sep = ".")
    }
  }
  # select unique set of samples along with SAMPID
  samples <- unique(subset(dataIn1, select = c(sampID, "SAMPID")))
  # Create vector of all samples in dataset
  allUIDs <- data.frame(SAMPID = unique(dataIn1$SAMPID), stringsAsFactors = F)

  # WATER COVER AND DEPTH
  # Subset data to just water cover parameters
  wdep <- subset(dataIn1, PARAMETER %in% c(
    "TIME", "MINIMUM_DEPTH", "MAXIMUM_DEPTH", "PREDOMINANT_DEPTH", "TOTAL_WATER",
    "WATER_NOVEG", "WATER_AQVEG", "WATER_EMERGVEG"
  ))
  # Ensure RESULT is numeric
  wdep$RESULT <- with(wdep, as.numeric(RESULT))

  # Subset variables
  wdep <- subset(wdep, select = c("SAMPID", "PLOT", "NPLOTS", "PARAMETER", "RESULT"))

  # Cast data wide and drop prefix from reshape()
  wdep1 <- reshape(wdep,
    idvar = c("SAMPID", "PLOT", "NPLOTS"), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(wdep1) <- gsub("RESULT\\.", "", names(wdep1))

  # Subset to remove missing and blank RESULT values for total water and aquatic vegetation
  # These are not calculated variables, just direct from the field form
  wdep2 <- subset(wdep, RESULT %nin% c("") & !is.na(RESULT) & PARAMETER %in% c("TOTAL_WATER", "WATER_NOVEG", "WATER_AQVEG", "WATER_EMERGVEG"))
  # Cast wide and drop prefix added by reshape()
  wdep2a.wide <- reshape(wdep2,
    idvar = c("SAMPID", "PLOT", "NPLOTS"), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(wdep2a.wide) <- gsub("RESULT\\.", "", names(wdep2a.wide))
  # Melt data frame and fill in missing RESULT values with 0
  varNames <- names(wdep2a.wide)[!names(wdep2a.wide) %in% c("SAMPID", "PLOT", "NPLOTS")]
  wdep2a <- reshape(wdep2a.wide,
    idvar = c("SAMPID", "PLOT", "NPLOTS"), direction = "long",
    times = varNames, varying = varNames,
    timevar = "PARAMETER", v.names = "RESULT"
  )
  wdep2a$RESULT <- with(wdep2a, ifelse(is.na(RESULT), 0, as.numeric(RESULT)))

  # Variables measured in 2011 and 2016 (and later) changed, so some metrics only apply to 2011
  if (survyear == "2011") {
    # calculate minimum depth and water cover
    wat1.min <- aggregate(
      x = list(MIN_H2O_DEPTH = wdep1$MINIMUM_DEPTH, MIN_COV_H2O = wdep1$TOTAL_WATER),
      by = wdep1[c("SAMPID", "NPLOTS")],
      FUN = function(x) {
        ifelse(any(!is.na(x)), min(x, na.rm = TRUE), NA)
      }
    )
    # Calculate maximum depth and water cover
    wat1.max <- aggregate(
      x = list(MAX_H2O_DEPTH = wdep1$MAXIMUM_DEPTH, MAX_COV_H2O = wdep1$TOTAL_WATER),
      by = wdep1[c("SAMPID", "NPLOTS")],
      FUN = function(x) {
        ifelse(any(!is.na(x)), max(x, na.rm = TRUE), NA)
      }
    )
    # Calculate first step in mean water depth
    wat1.sum <- aggregate(
      x = list(XH2O_DEPTH_AA = wdep1$PREDOMINANT_DEPTH),
      by = wdep1[c("SAMPID", "NPLOTS")],
      FUN = function(x) {
        ifelse(any(!is.na(x)),
          sum(x, na.rm = TRUE), NA
        )
      }
    )
    # Complete mean water depth calculation
    wat1.sum$XH2O_DEPTH_AA <- with(wat1.sum, round(XH2O_DEPTH_AA / NPLOTS, 2))

    # Merge metric data frames
    wat1 <- merge(wat1.min, wat1.max, by = c("SAMPID", "NPLOTS"))
    wat1 <- merge(wat1, wat1.sum, by = c("SAMPID", "NPLOTS"))
  } else { # If not 2011, more limited data collected for water
    # Minimum water cover
    wat1.min <- aggregate(
      x = list(MIN_COV_H2O = wdep1$TOTAL_WATER),
      by = wdep1[c("SAMPID", "NPLOTS")],
      FUN = function(x) {
        ifelse(any(!is.na(x)), min(x, na.rm = TRUE), NA)
      }
    )
    # Maximum water cover
    wat1.max <- aggregate(
      x = list(MAX_COV_H2O = wdep1$TOTAL_WATER),
      by = wdep1[c("SAMPID", "NPLOTS")],
      FUN = function(x) {
        ifelse(any(!is.na(x)), max(x, na.rm = TRUE), NA)
      }
    )
    # Sum depths for calculating mean water depth
    wat1.sum <- aggregate(
      x = list(XH2O_DEPTH_AA = wdep1$PREDOMINANT_DEPTH),
      by = wdep1[c("SAMPID", "NPLOTS")],
      FUN = function(x) {
        ifelse(any(!is.na(x)),
          sum(x, na.rm = TRUE), NA
        )
      }
    )
    # Finish calculating mean water depth
    wat1.sum$XH2O_DEPTH_AA <- with(wat1.sum, round(XH2O_DEPTH_AA / NPLOTS, 2))
    # Merge metric data frames
    wat1 <- merge(wat1.min, wat1.max, by = c("SAMPID", "NPLOTS"))
    wat1 <- merge(wat1, wat1.sum, by = c("SAMPID", "NPLOTS"))
  }
  # Melt output metrics
  varNames <- names(wat1)[!names(wat1) %in% c("SAMPID", "NPLOTS")]
  wat1a <- reshape(wat1,
    idvar = c("SAMPID"), direction = "long",
    times = varNames, varying = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )

  ## Fix Inf values to 0s
  wat1a$RESULT <- with(wat1a, ifelse(RESULT %in% c("Inf", "-Inf"), 0,
    ifelse(RESULT == "Inf_-Inf", "", RESULT)
  ))
  wat1a$METRIC <- with(wat1a, as.character(METRIC))
  wat1a$NPLOTS <- NULL # drop variable

  # Keep only positive values from water parameters
  wdep2a.sub <- subset(wdep2a, as.numeric(RESULT) > 0)
  # Sum water cover by type of vegetation
  wat2.sum <- aggregate(
    x = list(XCOV_H2O = wdep2a.sub$RESULT),
    by = wdep2a.sub[c("SAMPID", "PARAMETER", "NPLOTS")],
    FUN = function(x) {
      sum(as.numeric(x))
    }
  )
  # Count number of occurrences of each parameter
  wat2.freq <- aggregate(
    x = list(FREQ_H2O = wdep2a.sub$RESULT),
    by = wdep2a.sub[c("SAMPID", "PARAMETER", "NPLOTS")],
    FUN = function(x) {
      length(as.numeric(x))
    }
  )
  # Merge data frames
  wat2 <- merge(wat2.sum, wat2.freq, by = c("SAMPID", "PARAMETER", "NPLOTS"), all = TRUE)
  # Calculate water cover metrics by parameter
  wat2$FREQ_H2O <- with(wat2, round(FREQ_H2O / NPLOTS * 100, 2))
  wat2$XCOV_H2O <- with(wat2, round(XCOV_H2O / NPLOTS, 2))
  wat2$IMP_H2O <- with(wat2, round((FREQ_H2O + XCOV_H2O) / 2, 2))
  wat2$NPLOTS <- NULL # drop NPLOTS
  # Melt output data frame
  wat2a <- reshape(wat2,
    idvar = c("SAMPID", "PARAMETER"), direction = "long",
    times = c("FREQ_H2O", "XCOV_H2O", "IMP_H2O"),
    varying = c("FREQ_H2O", "XCOV_H2O", "IMP_H2O"),
    timevar = "METRIC", v.names = "RESULT"
  )
  # Update METRIC by combining PARAMETER and METRIC values, except for TOTAL_WATER parameter
  wat2a$METRIC <- with(wat2a, ifelse(PARAMETER == "TOTAL_WATER", as.character(METRIC),
    as.character(paste(METRIC, substring(PARAMETER, 7), sep = "_"))
  ))
  # Drop PARAMETER
  wat2a$PARAMETER <- NULL

  ## Now pull only PREDOMINANT_DEPTH that is not missing and >0 to calculate XH2O_DEPTH
  wdep1.pos <- subset(wdep1, PREDOMINANT_DEPTH > 0 & !is.na(PREDOMINANT_DEPTH))
  # Calculate mean water depth based on PREDOMINANT_DEPTH
  # Differs from XH2O_DEPTH_AA because denominator is number of plots sampled in _AA and just
  # number of values measured >0 for this version
  wat3 <- aggregate(
    x = list(RESULT = wdep1.pos$PREDOMINANT_DEPTH),
    by = wdep1.pos[c("SAMPID")],
    FUN = function(x) {
      mean(x, na.rm = TRUE)
    }
  )
  wat3$METRIC <- "XH2O_DEPTH"
  # Combine with previously calculated metrics
  watMet <- rbind(wat1a, wat2a, wat3)
  # Cast metric data frame wide and drop prefix added by reshape() to variable names
  watMet1 <- reshape(watMet,
    idvar = "SAMPID", direction = "wide",
    timevar = "METRIC", v.names = "RESULT"
  )
  names(watMet1) <- gsub("RESULT\\.", "", names(watMet1))

  # Create empty data frame based on year of survey - more metrics in 2011 than 2016.
  if (survyear == "2011") {
    empty_wat <- data.frame(t(rep(NA, 18)), stringsAsFactors = F)
    names(empty_wat) <- c(
      "MIN_H2O_DEPTH", "MAX_H2O_DEPTH", "XH2O_DEPTH_AA", "MIN_COV_H2O",
      "MAX_COV_H2O", "FREQ_H2O", "FREQ_H2O_AQVEG", "FREQ_H2O_EMERGVEG",
      "FREQ_H2O_NOVEG", "XCOV_H2O", "XCOV_H2O_AQVEG", "XCOV_H2O_EMERGVEG",
      "XCOV_H2O_NOVEG", "IMP_H2O", "IMP_H2O_AQVEG", "IMP_H2O_EMERGVEG",
      "IMP_H2O_NOVEG", "XH2O_DEPTH"
    )
  } else {
    empty_wat <- data.frame(t(rep(NA, 7)), stringsAsFactors = F)
    names(empty_wat) <- c(
      "XH2O_DEPTH_AA", "MIN_COV_H2O", "MAX_COV_H2O", "FREQ_H2O",
      "XCOV_H2O", "IMP_H2O", "XH2O_DEPTH"
    )
  }
  # Merge empty data frame with metric data frame, remove row with missing SAMPID
  watMet2 <- subset(merge(watMet1, empty_wat, all = TRUE), !is.na(SAMPID))
  # Merge with full list of samples to ensure all are included in output file
  watMet3 <- merge(allUIDs, watMet2, by = "SAMPID", all.x = T)
  # Melt data frame and substitute missing RESULT values with 0
  varNames <- names(watMet3)[!names(watMet3) %in% c("SAMPID", "NPLOTS")]
  watMet4 <- reshape(watMet3,
    idvar = "SAMPID", direction = "long",
    times = varNames, varying = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )

  watMet4$METRIC <- with(watMet4, as.character(METRIC))
  watMet4$RESULT <- with(watMet4, ifelse(is.na(RESULT), "0", RESULT))

  print("Done with water depth metrics")
  # Merge samples with data frame to get sampID variables back, then drop SAMPID
  # Also rename METRIC to PARAMETER
  watMetOut <- merge(samples, watMet4, by = "SAMPID")
  watMetOut$PARAMETER <- watMetOut$METRIC
  watMetOut <- watMetOut[, c(sampID, "PARAMETER", "RESULT")]

  return(watMetOut)
}
