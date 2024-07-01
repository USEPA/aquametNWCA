#' @export
#'
#' @title Calculate snag metrics
#'
#' @description This function calculates all of the snag
#' metrics from tree data, and is called by function
#' calcTreeMetrics().
#'
#' @param treeIn A data frame containing the following variables:
#' \itemize{
#' \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#' \item PLOT: Sample plot from which data were collected
#'
#' \item PARAMETER: specific measurement type
#'
#' \item RESULT: measured value
#' }
#' The following parameters are used in
#' calculating tree metrics: 'XXTHIN_SNAG', 'XTHIN_SNAG', 'THIN_SNAG',
#' 'JR_SNAG', 'THICK_SNAG', 'XTHICK_SNAG', 'XXTHICK_SNAG'.
#' Additional parameters or variables are ignored.
#'
#' @param nPlot A data frame with the
#' number of plots sampled associated with
#' each sample, with \emph{sampID} variables and NPLOTS.
#'
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default.
#'
#' @details If any of the parameters are missing, they are assumed
#' to be zeros, and metric values associated with any metrics that
#' cannot be calculated due to missing parameters are set to 0.
#'
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. The data frame contains the
#'   \emph{sampID} variables, PARAMETER, RESULT, where PARAMETER values include:
#'   TOTN_SNAGS, XN_SNAGS, TOTN_JR_SNAG, TOTN_THICK_SNAG, TOTN_THIN_SNAG,
#'   TOTN_XTHICK_SNAG, TOTN_XTHIN_SNAG, TOTN_XXTHIN_SNAG, XN_JR_SNAG,
#'   XN_THICK_SNAG, XN_THIN_SNAG, XN_XTHICK_SNAG, XN_XTHIN_SNAG, and
#'   XN_XXTHIN_SNAG. A list of metric descriptions is provided in the document
#'   named Tree_Metric_Descriptions.pdf included in the help directory for the
#'   package.
#'
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#'
#' @examples
#' head(TreesEx)
#'
#' nplots <- aggregate(
#'   x = list(NPLOTS = TreesEx$PLOT), by = TreesEx[c("UID")],
#'   FUN = function(x) {
#'     length(unique(x))
#'   }
#' )
#'
#' snagEx <- calcSnagMets(TreesEx, nPlot = nplots, sampID = "UID")
#'
#' head(snagEx)
#' unique(snagEx$PARAMETER)
calcSnagMets <- function(treeIn, nPlot, sampID = "UID") {
  # Merge nPlot data frame and treeIn data frame by variables in sampID
  treeIn <- merge(treeIn, nPlot, by = sampID, all.y = TRUE)

  # Create SAMPID variable based on those listed in sampID argument
  for (i in 1:length(sampID)) {
    if (i == 1) {
      treeIn$SAMPID <- treeIn[, sampID[i]]
    } else {
      treeIn$SAMPID <- paste(treeIn$SAMPID, treeIn[, sampID[i]], sep = ".")
    }
  }
  # Create data frame of unique samples in dataset, with sampID variables and SAMPID
  samples <- unique(subset(treeIn, select = c(sampID, "SAMPID")))
  # Create data frame of just unique SAMPID variables
  allUIDs <- data.frame(SAMPID = unique(treeIn$SAMPID), stringsAsFactors = F)

  # Subset data to just those parameters for snag counts
  snags <- subset(treeIn, PARAMETER %in% c(
    "XXTHIN_SNAG", "XTHIN_SNAG", "THIN_SNAG", "JR_SNAG",
    "THICK_SNAG", "XTHICK_SNAG", "XXTHICK_SNAG"
  ))

  # If the resulting subset has at least 1 row, proceed with calculations
  if (nrow(snags) > 0) {
    # Sum RESULT values by PARAMETER
    snagsOut <- aggregate(
      x = list(TOTN = snags$RESULT), by = snags[c("SAMPID", "PARAMETER", "NPLOTS")],
      FUN = function(x) {
        sum(as.numeric(x))
      }
    )

    # Calculate XN using results above and NPLOTS from nPlot input
    snagsOut$XN <- with(snagsOut, round((TOTN / NPLOTS), 2))
    # Melt data frame into long format
    snagsOut1 <- reshape(snagsOut,
      idvar = c("SAMPID", "PARAMETER"), direction = "long",
      varying = c("TOTN", "XN"), times = c("TOTN", "XN"),
      timevar = "METRIC", v.names = "RESULT"
    )
    # Combine METRIC and PARAMETER to create final METRIC name, then subset to just relevant variables
    snagsOut1$METRIC <- with(snagsOut1, paste(as.character(METRIC), PARAMETER, sep = "_"))
    snagsOut1 <- subset(snagsOut1, select = c("SAMPID", "METRIC", "RESULT"))
    # Calculate summed TOTN and XN as overall snags metrics, then round XN to 2 digits
    totsnags <- aggregate(
      x = list(TOTN_SNAGS = snagsOut$TOTN, XN_SNAGS = snagsOut$XN),
      by = snagsOut[c("SAMPID")], FUN = function(x) {
        sum(as.numeric(x))
      }
    )

    totsnags$XN_SNAGS <- with(totsnags, round(XN_SNAGS, 2))
    # Melt data frame into long format
    totsnags1 <- reshape(totsnags,
      idvar = "SAMPID", direction = "long",
      varying = c("TOTN_SNAGS", "XN_SNAGS"), times = c("TOTN_SNAGS", "XN_SNAGS"),
      timevar = "METRIC", v.names = "RESULT"
    )
    # Combine totals with other metrics
    allSnagsOut <- rbind(totsnags1, snagsOut1)
    # Cast resulting data frame into wide format
    allSnagsOut1 <- reshape(allSnagsOut,
      idvar = "SAMPID", direction = "wide",
      timevar = "METRIC", v.names = "RESULT"
    )
    # Remove prefix from variable names
    names(allSnagsOut1) <- gsub("RESULT\\.", "", names(allSnagsOut1))
    # Create empty data frame to ensure that all metrics are represented in output data frame
    empty_snags <- data.frame(t(rep(NA, 14)), stringsAsFactors = F)

    names(empty_snags) <- c(
      "TOTN_SNAGS", "XN_SNAGS", "TOTN_JR_SNAG", "TOTN_THICK_SNAG",
      "TOTN_THIN_SNAG", "TOTN_XTHICK_SNAG", "TOTN_XTHIN_SNAG",
      "TOTN_XXTHIN_SNAG", "XN_JR_SNAG", "XN_THICK_SNAG", "XN_THIN_SNAG",
      "XN_XTHICK_SNAG", "XN_XTHIN_SNAG", "XN_XXTHIN_SNAG"
    )
    # Merge empty data frame with output df and remove the empty row - this ensures all metrics
    # are included.
    allSnagsOut2 <- subset(merge(allSnagsOut1, empty_snags, all = TRUE), !is.na(SAMPID))

    # Merge with the all UIDs and fill in missing values with zeroes
    allSnagsOut3 <- merge(allUIDs, allSnagsOut2, by = "SAMPID", all.x = T)
    # Melt merge data frame
    varNames <- names(allSnagsOut3)[!names(allSnagsOut3) %in% c("SAMPID")]
    allSnagsOut4 <- reshape(allSnagsOut3,
      idvar = "SAMPID", direction = "long",
      varying = varNames, times = varNames,
      timevar = "METRIC", v.names = "RESULT"
    )
    allSnagsOut4$METRIC <- with(allSnagsOut4, as.character(METRIC))
    # Replace missing values with 0
    allSnagsOut4$RESULT <- with(allSnagsOut4, ifelse(is.na(RESULT), 0, RESULT))
  } else { # Otherwise create empty data frame
    # Create empty data frame as above
    empty_snags <- data.frame(t(rep(NA, 14)), stringsAsFactors = F)

    names(empty_snags) <- c(
      "TOTN_SNAGS", "XN_SNAGS", "TOTN_JR_SNAG", "TOTN_THICK_SNAG",
      "TOTN_THIN_SNAG", "TOTN_XTHICK_SNAG", "TOTN_XTHIN_SNAG",
      "TOTN_XXTHIN_SNAG", "XN_JR_SNAG", "XN_THICK_SNAG", "XN_THIN_SNAG",
      "XN_XTHICK_SNAG", "XN_XTHIN_SNAG", "XN_XXTHIN_SNAG"
    )
    # Merge list of SAMPIDs with empty data frame
    allSnagsOut <- merge(data.frame(SAMPID = rep(unique(treeIn$SAMPID)), stringsAsFactors = F), empty_snags, all = TRUE)
    # Remove the row without SAMPID
    allSnagsOut1 <- subset(allSnagsOut, !is.na(SAMPID))
    # Melt resulting data frame
    varNames <- names(allSnagsOut1)[!names(allSnagsOut3) %in% c("SAMPID")]
    allSnagsOut2 <- reshape(allSnagsOut1,
      idvar = "SAMPID", direction = "long",
      varying = varNames, times = varNames,
      timevar = "METRIC", v.names = "RESULT"
    )

    allSnagsOut2$METRIC <- with(allSnagsOut2, as.character(METRIC))
    # Set RESULT to 0 for all sites
    allSnagsOut2$RESULT <- 0

    allSnagsOut4 <- allSnagsOut2
  }
  print("Done with snag metrics")
  # Merge samples data frame to get back to the original sampID variables, rename METRIC to PARAMETER,
  # then drop METRIC and SAMPID
  treeOut <- merge(samples, allSnagsOut4, by = "SAMPID")
  treeOut$PARAMETER <- treeOut$METRIC
  treeOut$METRIC <- NULL
  treeOut$SAMPID <- NULL

  return(treeOut)
}
