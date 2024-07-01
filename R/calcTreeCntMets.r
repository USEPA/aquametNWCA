# Tree count metrics function
#' @export
#'
#' @title Calculate tree count metrics
#'
#' @description This function calculates tree count metrics
#' using tree data, and is called by function calcTreeMetrics()
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
#' The following parameters are used in calculating tree metrics:
#' 'XXTHIN_TREE', 'XTHIN_TREE', 'THIN_TREE', 'JR_TREE', 'THICK_TREE',
#' 'XTHICK_TREE', 'XXTHICK_TREE'. Additional parameters or variables
#' are ignored.
#' @param nPlot A data frame with the
#' number of plots sampled associated with
#' each sample, with \emph{sampID} variables and NPLOTS.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#'
#' @details If any of the parameters are missing, they are assumed to be
#' zeros, and metric values associated with any metrics that cannot be
#' calculated due to missing parameters are set to 0.
#'
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. The data frame contains the
#'   \emph{sampID} variables, PARAMETER, RESULT, where PARAMETER values include:
#'   TOTN_TREES, XN_TREES, TOTN_JR_TREE, TOTN_THICK_TREE, TOTN_THIN_TREE,
#'   TOTN_XTHICK_TREE, TOTN_XTHIN_TREE, TOTN_XXTHICK_TREE, TOTN_XXTHIN_TREE,
#'   XN_JR_TREE, XN_THICK_TREE, XN_THIN_TREE, XN_XTHICK_TREE, XN_XTHIN_TREE,
#'   XN_XXTHICK_TREE, XN_XXTHIN_TREE. A list of metric descriptions is provided
#'   in the document named Tree_Metric_Descriptions.pdf included in the help
#'   directory for the package.
#'
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#'
#' @examples
#' head(TreesEx)
#' nplots <- aggregate(
#'   x = list(NPLOTS = TreesEx$PLOT), by = TreesEx[c("UID")],
#'   FUN = function(x) {
#'     length(unique(x))
#'   }
#' )
#'
#' tcntEx <- calcTreeCntMets(TreesEx, nPlot = nplots, sampID = "UID")
#'
#' head(tcntEx)
#' unique(tcntEx$PARAMETER)
calcTreeCntMets <- function(treeIn, nPlot, sampID = "UID") {
  # Merge nPlot with treeIn data frame by sampID variables
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
  # Subset treeIn to keep only tree count parameters
  tcnts <- subset(treeIn, PARAMETER %in% c(
    "XXTHIN_TREE", "XTHIN_TREE", "THIN_TREE", "JR_TREE",
    "THICK_TREE", "XTHICK_TREE", "XXTHICK_TREE"
  ))

  # If resulting dataset has at least one row, perform calculations below
  if (nrow(tcnts) > 0) {
    # Sum counts by PARAMETER
    tcntsOut <- aggregate(
      x = list(TOTN = tcnts$RESULT), by = tcnts[c("SAMPID", "PARAMETER", "NPLOTS")],
      FUN = function(x) {
        sum(as.numeric(x))
      }
    )

    # Use above results to calculate XN and round result to 2 digits
    tcntsOut$XN <- with(tcntsOut, round(TOTN / NPLOTS, 2))
    # Melt data frame into long format
    tcntsOut1 <- reshape(tcntsOut,
      idvar = c("SAMPID", "PARAMETER"), direction = "long",
      varying = c("TOTN", "XN"), times = c("TOTN", "XN"),
      timevar = "METRIC", v.names = "RESULT"
    )
    # Update value of METRIC by combining METRIC and PARAMETER
    tcntsOut1$METRIC <- with(tcntsOut1, paste(as.character(METRIC), PARAMETER, sep = "_"))
    # Subset to relevant variables
    tcntsOut1 <- subset(tcntsOut1, select = c("SAMPID", "METRIC", "RESULT"))
    # Calculate totals across all sizes for both TOTN and XN, then round to 2 digits
    tottrees <- aggregate(
      x = list(TOTN_TREES = tcntsOut$TOTN, XN_TREES = tcntsOut$XN),
      by = tcntsOut[c("SAMPID")], FUN = sum
    )
    tottrees$XN_TREES <- with(tottrees, round(XN_TREES, 2))
    # Melt data frame in order to merge with parameter-based metrics above
    tottrees1 <- reshape(tottrees,
      idvar = "SAMPID", direction = "long",
      varying = c("TOTN_TREES", "XN_TREES"), times = c("TOTN_TREES", "XN_TREES"),
      timevar = "METRIC", v.names = "RESULT"
    )
    # Combine totals with other metrics
    allTreesOut <- rbind(tottrees1, tcntsOut1)
    # Cast data frame wide and remove prefix in resulting data frame
    allTreesOut1 <- reshape(allTreesOut,
      idvar = "SAMPID", direction = "wide",
      timevar = "METRIC", v.names = "RESULT"
    )
    names(allTreesOut1) <- gsub("RESULT\\.", "", names(allTreesOut1))
    # Create empty data frame with names of all metrics that should be calculated
    empty_trees <- data.frame(t(rep(NA, 16)), stringsAsFactors = F)
    names(empty_trees) <- c(
      "TOTN_TREES", "XN_TREES", "TOTN_JR_TREE", "TOTN_THICK_TREE", "TOTN_THIN_TREE", "TOTN_XTHICK_TREE", "TOTN_XTHIN_TREE",
      "TOTN_XXTHICK_TREE", "TOTN_XXTHIN_TREE", "XN_JR_TREE", "XN_THICK_TREE", "XN_THIN_TREE", "XN_XTHICK_TREE", "XN_XTHIN_TREE",
      "XN_XXTHICK_TREE", "XN_XXTHIN_TREE"
    )

    # Merge empty data frame with metrics, then drop row with missing SAMPID
    allTreesOut2 <- subset(merge(allTreesOut1, empty_trees, all = TRUE), !is.na(SAMPID))
    # Merge resulting data frame with list of allUIDs to make sure all samples are represented
    allTreesOut3 <- merge(allUIDs, allTreesOut2, by = "SAMPID", all.x = T)

    # Melt resulting data frame, then set missing RESULT values to 0
    varNames <- names(allTreesOut3)[!names(allTreesOut3) %in% c("SAMPID")]
    allTreesOut4 <- reshape(allTreesOut3,
      idvar = "SAMPID", direction = "long",
      varying = varNames, times = varNames,
      timevar = "METRIC", v.names = "RESULT"
    )
    allTreesOut4$METRIC <- with(allTreesOut4, as.character(METRIC))
    allTreesOut4$RESULT <- with(allTreesOut4, ifelse(is.na(RESULT), 0, RESULT))
    # Cast data frame wide and remove prefix added by reshape() function
    allTreesOut4.wide <- reshape(allTreesOut4,
      idvar = "SAMPID", direction = "wide",
      timevar = "METRIC", v.names = "RESULT"
    )
    names(allTreesOut4.wide) <- gsub("RESULT\\.", "", names(allTreesOut4.wide))
    # Calculate additional TOTN and XN metrics based on combinations of previously calculated metrics
    allTreesOut4.wide$TOTN_SMALL <- with(allTreesOut4.wide, TOTN_XXTHIN_TREE + TOTN_XTHIN_TREE)
    allTreesOut4.wide$TOTN_MID <- with(allTreesOut4.wide, TOTN_THIN_TREE + TOTN_JR_TREE)
    allTreesOut4.wide$TOTN_LARGE <- with(allTreesOut4.wide, TOTN_THICK_TREE + TOTN_XTHICK_TREE + TOTN_XXTHICK_TREE)
    allTreesOut4.wide$XN_SMALL <- with(allTreesOut4.wide, XN_XXTHIN_TREE + XN_XTHIN_TREE)
    allTreesOut4.wide$XN_MID <- with(allTreesOut4.wide, XN_THIN_TREE + XN_JR_TREE)
    allTreesOut4.wide$XN_LARGE <- with(allTreesOut4.wide, XN_THICK_TREE + XN_XTHICK_TREE + XN_XXTHICK_TREE)
    # Melt metrics into long format
    varNames <- names(allTreesOut4.wide)[!names(allTreesOut4.wide) %in% c("SAMPID")]
    allTreesOut5 <- reshape(allTreesOut4.wide,
      idvar = "SAMPID", direction = "long",
      varying = varNames, times = varNames,
      timevar = "METRIC", v.names = "RESULT"
    )
  } else { # If no tree counts in input dataset, create empty data frame, then add variable names
    empty_trees <- data.frame(t(rep(NA, 16)), stringsAsFactors = F)
    names(empty_trees) <- c(
      "TOTN_TREES", "XN_TREES", "TOTN_JR_TREE", "TOTN_THICK_TREE",
      "TOTN_THIN_TREE", "TOTN_XTHICK_TREE", "TOTN_XTHIN_TREE",
      "TOTN_XXTHICK_TREE", "TOTN_XXTHIN_TREE", "XN_JR_TREE", "XN_THICK_TREE",
      "XN_THIN_TREE", "XN_XTHICK_TREE", "XN_XTHIN_TREE",
      "XN_XXTHICK_TREE", "XN_XXTHIN_TREE", "TOTN_SMALL", "TOTN_MID", "TOTN_LARGE",
      "XN_SMALL", "XN_MID", "XN_LARGE"
    )
    # Merge empty data frame with list of samples
    allTreesOut <- merge(data.frame(SAMPID = rep(unique(treeIn$SAMPID)), stringsAsFactors = F),
      empty_trees,
      all = TRUE
    )
    # Subset to remove row with missing SAMPID
    allTreesOut1 <- subset(allTreesOut, !is.na(SAMPID))
    # Melt output data frame and set all RESULT values to 0
    varNames <- names(allTreesOut1)[!names(allTreesOut1) %in% c("SAMPID")]
    allTreesOut2 <- reshape(allTreesOut2,
      idvar = "SAMPID", direction = "long",
      varying = varNames, times = varNames,
      timevar = "METRIC", v.names = "RESULT"
    )

    allTreesOut2$METRIC <- with(allTreesOut2, as.character(METRIC))
    allTreesOut2$RESULT <- 0

    allTreesOut5 <- allTreesOut2
  }

  print("Done with tree count metrics")
  # Merge with samples data frame to get back sampID variables, rename METRIC to PARAMETER,
  # then drop METRIC and SAMPID variables.
  treeOut <- merge(samples, allTreesOut5, by = "SAMPID")
  treeOut$PARAMETER <- treeOut$METRIC
  treeOut$METRIC <- NULL
  treeOut$SAMPID <- NULL

  return(treeOut)
}
