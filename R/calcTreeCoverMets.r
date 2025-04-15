## Tree species and cover metric function
#' @export
#'
#' @title Calculate tree cover metrics
#'
#' @description This function calculates tree cover metrics
#' using tree data, and is called by function calcTreeMetrics().
#'
#' @param treeIn A data frame containing the following variables:
#' \itemize{
#' \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#' \item PAGE: Page from field form
#'
#' \item LINE: Line number from field form
#'
#' \item PLOT: Sample plot from which data were collected
#'
#' \item PARAMETER: specific measurement type
#'
#' \item RESULT: measured value
#' }
#' The following parameters are used in
#' calculating tree metrics: 'VSMALL_TREE', 'SMALL_TREE', 'LMED_TREE'
#' , 'HMED_TREE', 'TALL_TREE', 'VTALL_TREE'.  Additional parameters
#' or variables are ignored.
#' @param nPlot A data frame with the
#' number of plots sampled associated with
#' each sample, with \emph{sampID} variables and NPLOTS.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#'
#' @details If any of the parameters are missing, they are assumed
#' to be zeros, and metric values associated with any metrics that
#' cannot be calculated due to missing parameters are set to 0.
#'
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. The data frame contains the
#'   \emph{sampID} variables, PARAMETER, RESULT, where PARAMETER values include:
#'   N_TREESPP, N_TALL_TREE, N_HMED_TREE, N_LMED_TREE, N_SMALL_TREE,
#'   N_VSMALL_TREE, N_VTALL_TREE, N_TREE_UPPER, N_TREE_MID, N_TREE_GROUND,
#'   PCTN_TREE_UPPER, PCTN_TREE_MID, PCTN_TREE_GROUND. A list of metric
#'   descriptions is provided in the document named \href{https://github.com/USEPA/aquametNWCA/blob/main/inst/Tree_Metric_Descriptions.pdf}{Tree_Metric_Descriptions.pdf}
#'   included in the help directory for the package.
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
#' tcvrEx <- calcTreeCoverMets(TreesEx, nPlot = nplots, sampID = "UID")
#'
#' head(tcvrEx)
#' unique(tcvrEx$PARAMETER)
calcTreeCoverMets <- function(treeIn, nPlot, sampID = "UID") {
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

  ##### TREE SPECIES METRICS ##############################
  # Subset treeIn to keep only tree height parameters
  tcvr <- subset(treeIn, PARAMETER %in% c(
    "VSMALL_TREE", "SMALL_TREE", "LMED_TREE", "HMED_TREE",
    "TALL_TREE", "VTALL_TREE"
  ))
  # If resulting dataset has at least one row, perform calculations below
  if (nrow(tcvr) > 0) {
    # Subset input data to tree species data and cast data frame wide and
    # remove prefix added by reshape()
    tspp <- reshape(
      subset(treeIn, PARAMETER == "TREE_SPECIES",
        select = c("SAMPID", "PAGE", "LINE", "PLOT", "PARAMETER", "RESULT")
      ),
      idvar = c("SAMPID", "PAGE", "LINE", "PLOT"), direction = "wide",
      timevar = "PARAMETER", v.names = "RESULT"
    )
    names(tspp) <- gsub("RESULT\\.", "", names(tspp))

    # Merge tree species and tree cover data
    tcvr1 <- merge(tspp, tcvr, by = c("SAMPID", "PLOT", "PAGE", "LINE"), all = TRUE)
    # Create alternate parameter values by combining original metric values
    tcvr1$PARAM_ALT <- NA
    tcvr1$PARAM_ALT[tcvr1$PARAMETER %in% c("VSMALL_TREE", "SMALL_TREE")] <- "TREE_GROUND"
    tcvr1$PARAM_ALT[tcvr1$PARAMETER %in% c("LMED_TREE", "HMED_TREE")] <- "TREE_MID"
    tcvr1$PARAM_ALT[tcvr1$PARAMETER %in% c("TALL_TREE", "VTALL_TREE")] <- "TREE_UPPER"
    # Calculate number of species
    ntrspp <- aggregate(
      x = list(N_TREESPP = tcvr1$TREE_SPECIES),
      by = tcvr1[c("SAMPID")], FUN = function(x) {
        length(unique(x))
      }
    )
    # Merge count with cover and species data from above
    totspp <- merge(tcvr1, ntrspp, by = c("SAMPID"))
    # Melt data frame into long format for metric calculation
    totspp1 <- reshape(ntrspp,
      idvar = "SAMPID", direction = "long",
      varying = "N_TREESPP", times = "N_TREESPP",
      timevar = "METRIC", v.names = "RESULT"
    )

    # Subset to values above 0 and not missing
    totspp.pos <- subset(totspp, RESULT != "0" & !is.na(RESULT),
      select = c("SAMPID", "PARAMETER", "PARAM_ALT", "TREE_SPECIES")
    )
    # Calculate number of unique species by PARAMETER value (N metrics)
    tspp1a <- aggregate(
      x = list(N = totspp.pos$TREE_SPECIES),
      by = totspp.pos[c("SAMPID", "PARAMETER")],
      FUN = function(x) {
        length(unique(x))
      }
    )
    # Drop for any cases where PARAMETER is missing
    tspp1a <- subset(tspp1a, !is.na(PARAMETER))
    # Calculate number of unique species by PARAM_ALT value
    tspp1b <- aggregate(
      x = list(N = totspp.pos$TREE_SPECIES),
      by = totspp.pos[c("SAMPID", "PARAM_ALT")],
      FUN = function(x) {
        length(unique(x))
      }
    )
    # Drop for any cases where PARAM_ALT missing
    tspp1b <- subset(tspp1b, !is.na(PARAM_ALT))
    # Merge above data frame with data frame with total number of species
    tspp1c <- merge(tspp1b, ntrspp, by = "SAMPID")
    tspp1c <- subset(tspp1c, !is.na(PARAM_ALT))
    # Calculate PCTN (% species) and round to 2 digits
    tspp1c$PCTN <- with(tspp1c, round((N / N_TREESPP) * 100, 2))
    # Drop total number of tree species metric
    tspp1c$N_TREESPP <- NULL

    # Melt first and third data frames created above, updating metric name by
    # combining with PARAMETER value
    tspp2a <- reshape(tspp1a,
      idvar = c("SAMPID", "PARAMETER"), direction = "long",
      varying = "N", times = "N",
      timevar = "METRIC", v.names = "RESULT"
    )
    tspp2a$METRIC <- with(tspp2a, paste(METRIC, PARAMETER, sep = "_"))

    tspp2c <- reshape(tspp1c,
      idvar = c("SAMPID", "PARAM_ALT"), direction = "long",
      varying = c("PCTN", "N"), times = c("PCTN", "N"),
      timevar = "METRIC", v.names = "RESULT"
    )
    tspp2c$METRIC <- with(tspp2c, paste(METRIC, PARAM_ALT, sep = "_"))
    # Combine total number of species data frame and two sets of metrics above
    tsppOut <- rbind(totspp1, subset(tspp2a, select = -PARAMETER), subset(tspp2c, select = -PARAM_ALT))
    # Cast resulting data frame wide and remove prefix added by reshape
    tsppOut1 <- reshape(tsppOut,
      idvar = c("SAMPID"), direction = "wide",
      timevar = "METRIC", v.names = "RESULT"
    )
    names(tsppOut1) <- gsub("RESULT\\.", "", names(tsppOut1))

    # Create empty data frame with full set of expected metrics
    empty_tspp <- data.frame(t(rep(NA, 13)))
    names(empty_tspp) <- c(
      "N_TREESPP", "N_TALL_TREE", "N_HMED_TREE", "N_LMED_TREE", "N_SMALL_TREE",
      "N_VSMALL_TREE", "N_VTALL_TREE", "N_TREE_UPPER",
      "N_TREE_MID", "N_TREE_GROUND", "PCTN_TREE_UPPER", "PCTN_TREE_MID",
      "PCTN_TREE_GROUND"
    )
    # Merge output data frame with empty data frame and remove row missing SAMPID
    tsppOut2 <- subset(merge(tsppOut1, empty_tspp, all = TRUE), !is.na(SAMPID))
    # Merge full set of UIDs with resulting data frame to ensure all samples included
    tsppOut3 <- merge(allUIDs, tsppOut2, by = "SAMPID", all.x = T)
    # Melt resulting data frame, then replace missing RESULT values with 0
    varNames <- names(tsppOut3)[!names(tsppOut3) %in% c("SAMPID")]
    tsppOut4 <- reshape(tsppOut3,
      idvar = "SAMPID", direction = "long",
      varying = varNames, times = varNames,
      timevar = "METRIC", v.names = "RESULT"
    )
    tsppOut4$METRIC <- with(tsppOut4, as.character(METRIC))
    tsppOut4$RESULT <- with(tsppOut4, ifelse(is.na(RESULT), 0, RESULT))
  } else {
    # If no values for these parameters in input data frame, create empty data frame
    empty_tspp <- data.frame(t(rep(NA, 13)), stringsAsFactors = F)
    names(empty_tspp) <- c(
      "N_TREESPP", "N_TALL_TREE", "N_HMED_TREE", "N_LMED_TREE", "N_SMALL_TREE",
      "N_VSMALL_TREE", "N_VTALL_TREE", "N_TREE_UPPER",
      "N_TREE_MID", "N_TREE_GROUND", "PCTN_TREE_UPPER", "PCTN_TREE_MID",
      "PCTN_TREE_GROUND"
    )
    # Merge data frame with list of all UIDs in input data frame, dropping row with missing SAMPID
    tsppOut <- merge(data.frame(SAMPID = rep(unique(treeIn$SAMPID)), stringsAsFactors = F), empty_tspp, all = TRUE)
    tsppOut1 <- subset(tsppOut, !is.na(SAMPID))
    # Melt data frame and set all RESULT values to 0
    varNames <- names(tsppOut1)[!names(tsppOut1) %in% c("SAMPID")]
    tsppOut2 <- reshape(tsppOut1,
      idvar = "SAMPID", direction = "long",
      varying = varNames, times = varNames,
      timevar = "METRIC", v.names = "RESULT"
    )
    tsppOut2$METRIC <- with(tsppOut2, as.character(METRIC))
    tsppOut2$RESULT <- 0
    tsppOut4 <- tsppOut2
  }
  print("Done with tree species metrics")

  ## TREE COVER METRICS ###########
  ## Sum by species within plot
  if (nrow(tcvr) > 0) {
    # Sum cover by PARAMETER
    tcvr2a <- aggregate(
      x = list(COV = tcvr1$RESULT),
      by = tcvr1[c("SAMPID", "PLOT", "NPLOTS", "PARAMETER", "TREE_SPECIES")],
      FUN = function(x) {
        sum(as.numeric(x))
      }
    )
    # Cap values at 100 percent and keep only non-zero values
    tcvr2a$COV <- ifelse(tcvr2a$COV > 100, 100, tcvr2a$COV)
    tcvr2a <- subset(tcvr2a, COV != 0)
    # Sum cover by PARAM_ALT
    tcvr2b <- aggregate(
      x = list(COV = tcvr1$RESULT),
      by = tcvr1[c("SAMPID", "PLOT", "NPLOTS", "PARAM_ALT", "TREE_SPECIES")],
      FUN = function(x) {
        sum(as.numeric(x))
      }
    )
    # Cap values at 100 percent and keep only non-zero values
    tcvr2b$COV <- ifelse(tcvr2b$COV > 100, 100, tcvr2b$COV)
    tcvr2b <- subset(tcvr2b, COV != 0)
    # Take unique values of NPLOTS
    tcvr3a.uniq <- aggregate(
      x = list(uniqN = tcvr2a$NPLOTS),
      by = tcvr2a[c("SAMPID", "PARAMETER")], FUN = unique
    )
    # Calculate number of unique plots by parameter
    tcvr3a.length <- aggregate(
      x = list(uniqPlot = tcvr2a$PLOT),
      by = tcvr2a[c("SAMPID", "PARAMETER")],
      FUN = function(x) {
        length(unique(x))
      }
    )
    # Sum Cover by PARAMETER
    tcvr3a.sum <- aggregate(
      x = list(sumcov = tcvr2a$COV),
      by = tcvr2a[c("SAMPID", "PARAMETER")],
      FUN = function(x) {
        sum(as.numeric(x))
      }
    )
    # Merge length with unique plots and with cover data frame
    tcvr3a <- merge(tcvr3a.length, tcvr3a.uniq, by = c("SAMPID", "PARAMETER"))
    tcvr3a <- merge(tcvr3a, tcvr3a.sum, by = c("SAMPID", "PARAMETER"))
    # Now complete calculations of FREQ, XCOV, and IMP for each sample and parameter
    tcvr3a$FREQ <- with(tcvr3a, round((uniqPlot / uniqN) * 100, 2))
    tcvr3a$XCOV <- with(tcvr3a, round(sumcov / uniqN, 2))
    tcvr3a$IMP <- with(tcvr3a, round((FREQ + XCOV) / 2, 2))
    # Subset to select relevant variables
    tcvr3a <- subset(tcvr3a, select = c("SAMPID", "PARAMETER", "FREQ", "XCOV", "IMP"))
    # Again get unique NPLOTS values by sample
    tcvr3b.uniq <- aggregate(
      x = list(uniqN = tcvr2b$NPLOTS),
      by = tcvr2b[c("SAMPID", "PARAM_ALT")], FUN = unique
    )
    # Calculate number of unique plots by PARAM_ALT
    tcvr3b.length <- aggregate(
      x = list(uniqPlot = tcvr2b$PLOT),
      by = tcvr2b[c("SAMPID", "PARAM_ALT")],
      FUN = function(x) {
        length(unique(x))
      }
    )
    # Sum Cover by PARAM_ALT
    tcvr3b.sum <- aggregate(
      x = list(sumcov = tcvr2b$COV),
      by = tcvr2b[c("SAMPID", "PARAM_ALT")],
      FUN = function(x) {
        sum(as.numeric(x))
      }
    )
    # Merge length with unique plots and with cover data frame
    tcvr3b <- merge(tcvr3b.length, tcvr3b.uniq, by = c("SAMPID", "PARAM_ALT"))
    tcvr3b <- merge(tcvr3b, tcvr3b.sum, by = c("SAMPID", "PARAM_ALT"))
    # Now complete calculations of FREQ, XCOV, and IMP for each sample and PARAM_ALT
    tcvr3b$FREQ <- with(tcvr3b, round((uniqPlot / uniqN) * 100, 2))
    tcvr3b$XCOV <- with(tcvr3b, round(sumcov / uniqN, 2))
    tcvr3b$IMP <- with(tcvr3b, round((FREQ + XCOV) / 2, 2))
    # Subset to keep relevant variables
    tcvr3b <- subset(tcvr3b, select = c("SAMPID", "PARAM_ALT", "FREQ", "XCOV", "IMP"))
    # Melt data frame by PARAMETER and update METRIC value by combining METRIC and PARAMETER values
    varNames.a <- names(tcvr3a)[!names(tcvr3a) %in% c("SAMPID", "PARAMETER")]
    tcvr4a <- reshape(tcvr3a,
      idvar = c("SAMPID", "PARAMETER"), direction = "long",
      varying = varNames.a, times = varNames.a,
      timevar = "METRIC", v.names = "RESULT"
    )
    tcvr4a$METRIC <- with(tcvr4a, paste(METRIC, PARAMETER, sep = "_"))
    # Melt data frame by PARAM_ALT and update METRIC value by combining METRIC and PARAMETER values
    varNames.b <- names(tcvr3b)[!names(tcvr3b) %in% c("SAMPID", "PARAM_ALT")]
    tcvr4b <- reshape(tcvr3b,
      idvar = c("SAMPID", "PARAM_ALT"), direction = "long",
      varying = varNames.b, times = varNames.b,
      timevar = "METRIC", v.names = "RESULT"
    )
    tcvr4b$METRIC <- with(tcvr4b, paste(METRIC, PARAM_ALT, sep = "_"))
    # Combine resulting data frames after dropping PARAMETER and PARAM_ALT variables
    tcvrOut <- rbind(subset(tcvr4a, select = -PARAMETER), subset(tcvr4b, select = -PARAM_ALT))
    # Cast resulting data frame and drop prefix added by reshape()
    tcvrOut1 <- reshape(tcvrOut,
      idvar = "SAMPID", direction = "wide",
      timevar = "METRIC", v.names = "RESULT"
    )
    names(tcvrOut1) <- gsub("RESULT\\.", "", names(tcvrOut1))
    # Create empty data frame
    empty_tcvr <- data.frame(t(rep(NA, 27)), stringsAsFactors = F)
    names(empty_tcvr) <- c(
      "FREQ_TALL_TREE", "FREQ_HMED_TREE", "FREQ_LMED_TREE", "FREQ_SMALL_TREE",
      "FREQ_VSMALL_TREE", "FREQ_VTALL_TREE", "XCOV_TALL_TREE",
      "XCOV_HMED_TREE", "XCOV_LMED_TREE", "XCOV_SMALL_TREE", "XCOV_VSMALL_TREE",
      "XCOV_VTALL_TREE", "IMP_TALL_TREE", "IMP_HMED_TREE",
      "IMP_LMED_TREE", "IMP_SMALL_TREE", "IMP_VSMALL_TREE", "IMP_VTALL_TREE",
      "FREQ_TREE_UPPER", "FREQ_TREE_MID", "FREQ_TREE_GROUND",
      "XCOV_TREE_UPPER", "XCOV_TREE_MID", "XCOV_TREE_GROUND",
      "IMP_TREE_UPPER", "IMP_TREE_MID", "IMP_TREE_GROUND"
    )
    # Merge with output data frame and drop row with missing SAMPID
    tcvrOut2 <- subset(merge(tcvrOut1, empty_tcvr, all = TRUE), !is.na(SAMPID))
    # Merge with allUIDs to make sure all samples represented
    tcvrOut3 <- merge(allUIDs, tcvrOut2, by = "SAMPID", all.x = T)
    # Melt data frame and set missing RESULT values to 0
    varNames <- names(tcvrOut3)[!names(tcvrOut3) %in% c("SAMPID")]
    tcvrOut4 <- reshape(tcvrOut3,
      idvar = "SAMPID", direction = "long",
      varying = varNames, times = varNames,
      timevar = "METRIC", v.names = "RESULT"
    )

    tcvrOut4$METRIC <- with(tcvrOut4, as.character(METRIC))
    tcvrOut4$RESULT <- with(tcvrOut4, ifelse(is.na(RESULT), 0, RESULT))
  } else { # If no data for this subset, create empty data frame
    empty_tcvr <- data.frame(t(rep(NA, 27)), stringsAsFactors = F)
    names(empty_tcvr) <- c(
      "FREQ_TALL_TREE", "FREQ_HMED_TREE", "FREQ_LMED_TREE", "FREQ_SMALL_TREE",
      "FREQ_VSMALL_TREE", "FREQ_VTALL_TREE", "XCOV_TALL_TREE",
      "XCOV_HMED_TREE", "XCOV_LMED_TREE", "XCOV_SMALL_TREE", "XCOV_VSMALL_TREE",
      "XCOV_VTALL_TREE", "IMP_TALL_TREE", "IMP_HMED_TREE",
      "IMP_LMED_TREE", "IMP_SMALL_TREE", "IMP_VSMALL_TREE", "IMP_VTALL_TREE",
      "FREQ_TREE_UPPER", "FREQ_TREE_MID", "FREQ_TREE_GROUND",
      "XCOV_TREE_UPPER", "XCOV_TREE_MID", "XCOV_TREE_GROUND",
      "IMP_TREE_UPPER", "IMP_TREE_MID", "IMP_TREE_GROUND"
    )
    # Merge empty data frame with SAMPID list, drop missing SAMPID
    tcvrOut <- merge(data.frame(SAMPID = rep(unique(treeIn$SAMPID)), stringsAsFactors = F),
      empty_tcvr,
      all = TRUE
    )
    tcvrOut1 <- subset(tcvrOut, !is.na(SAMPID))
    # Melt data frame and set RESULT to 0
    varNames <- names(tcvrOut1)[!names(tcvrOut1) %in% c("SAMPID")]
    tcvrOut2 <- reshape(tcvrOut1,
      idvar = "SAMPID", direction = "long",
      varying = varNames, times = varNames,
      timevar = "METRIC", v.names = "RESULT"
    )

    tcvrOut2$METRIC <- with(tcvrOut2, as.character(METRIC))
    tcvrOut2$RESULT <- with(tcvrOut2, ifelse(is.na(RESULT), 0, RESULT))

    tcvrOut4 <- tcvrOut2
  }
  print("Done with tree cover metrics")
  # Combine cover and species output data frames
  treeOut <- rbind(tsppOut4, tcvrOut4)
  treeOut$PARAMETER <- treeOut$METRIC
  treeOut$METRIC <- NULL
  # Merge with samples to get back original sampID variables, then drop SAMPID
  treeOut.1 <- merge(samples, treeOut, by = "SAMPID")
  treeOut.1 <- subset(treeOut.1, select = -SAMPID)

  return(treeOut.1)
}
