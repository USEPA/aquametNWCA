#### VASCULAR STRATA
#' @export
#'
#' @title Calculate vascular strata metrics
#'
#' @description Calculate vascular strata metrics using
#' vegetation type data. This function is called by
#' calcVtype_GcovMets().
#'
#' @param dataIn A data frame containing the following
#' variables:
#' \itemize{
#' \item sampID: Variable(s) identified in the \emph{sampID} argument
#'
#' \item PLOT: Sample plot from which data were collected
#'
#' \item PARAMETER: specific measurement type
#'
#' \item RESULT: measured value
#' }
#' The following parameters are used in
#' calculating vascular strata metrics: 'SUBMERGED_AQ', 'FLOATING_AQ',
#' 'LIANAS', 'VTALL_VEG', 'TALL_VEG', 'HMED_VEG', 'MED_VEG',
#' 'SMALL_VEG', 'VSMALL_VEG'.
#' Additional parameters or variables are ignored.
#' @param nPlot A data frame with the
#' number of plots sampled associated with each
#' sample with \emph{sampID} variables and NPLOTS
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples,
#' 'UID' by default.
#'
#' @details If any of the parameters are missing, they are
#' assumed to be zeros (if numeric), and metric values associated
#' with any metrics that cannot be calculated due to missing
#' parameters are set to a standardized value.
#'
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. Data frame contains
#'   \emph{sampID} variables, PARAMETER, RESULT, where values of PARAMETER
#'   consist of the metric name concatenated with trait value, with valid values
#'   of: FREQ_FLOATING_AQ, FREQ_HMED_VEG, FREQ_LIANAS, FREQ_MED_VEG,
#'   FREQ_SMALL_VEG, FREQ_SUBMERGED_AQ, FREQ_TALL_VEG, FREQ_VSMALL_VEG,
#'   FREQ_VTALL_VEG, IMP_FLOATING_AQ, IMP_HMED_VEG, IMP_LIANAS, IMP_MED_VEG,
#'   IMP_SMALL_VEG, IMP_SUBMERGED_AQ, IMP_TALL_VEG, IMP_VSMALL_VEG,
#'   IMP_VTALL_VEG, XCOV_FLOATING_AQ, XCOV_HMED_VEG, XCOV_LIANAS, XCOV_MED_VEG,
#'   XCOV_SMALL_VEG, XCOV_SUBMERGED_AQ, XCOV_TALL_VEG, XCOV_VSMALL_VEG,
#'   XCOV_VTALL_VEG, XRCOV_FLOATING_AQ, XRCOV_HMED_VEG, XRCOV_LIANAS,
#'   XRCOV_MED_VEG, XRCOV_SMALL_VEG, XRCOV_SUBMERGED_AQ, XRCOV_TALL_VEG,
#'   XRCOV_VSMALL_VEG, XRCOV_VTALL_VEG, H_VASC_STRATA, J_VASC_STRATA,
#'   D_VASC_STRATA. A list of metric descriptions is provided in the document
#'   \href{https://github.com/USEPA/aquametNWCA/blob/main/inst/VType_GrdCvr_Metric_Descriptions.pdf}{VegTypes_GrdCover_Metric_Descriptions.pdf} included in the help directory
#'   for the package.
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
#' stratEx <- calcVascStratMets(Vtype_GrCovEx, nplots)
#'
#' head(stratEx)
#' unique(stratEx$METRIC)
calcVascStratMets <- function(dataIn, nPlot, sampID = "UID") {
  # Merge dataIn with nPlot
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
  # Subset to just parameters for vascular strata and drop any with 0 values or missing RESULT
  vstrat <- subset(dataIn1, PARAMETER %in% c(
    "SUBMERGED_AQ", "FLOATING_AQ", "LIANAS", "VTALL_VEG",
    "TALL_VEG", "HMED_VEG", "MED_VEG",
    "SMALL_VEG", "VSMALL_VEG"
  ) & !is.na(RESULT) & RESULT != "0")
  # Sum cover for total cover of vascular strata
  vstrat.sum <- aggregate(
    x = list(XTOTCOV_VASC_STRATA = vstrat$RESULT),
    by = vstrat[c("SAMPID", "NPLOTS")],
    FUN = function(x) {
      sum(as.numeric(x))
    }
  )
  # Complete calculation of mean cover by dividing by number of plots
  vstrat.sum$XTOTCOV_VASC_STRATA <- with(vstrat.sum, XTOTCOV_VASC_STRATA / NPLOTS)
  # Calculate number of vascular strata and sampled plots
  vstrat.cnt <- aggregate(
    x = list(N_VASC_STRATA = vstrat$PARAMETER, PLOTSAMP = vstrat$PLOT),
    by = vstrat[c("SAMPID")], FUN = function(x) {
      length(unique(x))
    }
  )
  # Merge back with previous data frames
  vstrat <- merge(vstrat, vstrat.sum, by = c("SAMPID", "NPLOTS"))
  vstrat <- merge(vstrat, vstrat.cnt, by = "SAMPID")
  # Number of vascular strata by plot for later calculation
  vstratPlot.n <- aggregate(
    x = list(N_VSTRATA = vstrat$PARAMETER),
    by = vstrat[c(
      "SAMPID", "PLOT", "N_VASC_STRATA", "NPLOTS",
      "XTOTCOV_VASC_STRATA", "PLOTSAMP"
    )],
    FUN = function(x) {
      length(unique(x))
    }
  )
  # Sum cover values by plot
  vstratPlot.sum <- aggregate(
    x = list(SUM_PLOT = vstrat$RESULT),
    by = vstrat[c(
      "SAMPID", "PLOT", "N_VASC_STRATA", "NPLOTS",
      "XTOTCOV_VASC_STRATA", "PLOTSAMP"
    )],
    FUN = function(x) {
      sum(as.numeric(x))
    }
  )
  # Merge these two calculations together
  vstratPlot <- merge(vstratPlot.n, vstratPlot.sum, by = c(
    "SAMPID", "PLOT", "N_VASC_STRATA",
    "NPLOTS", "XTOTCOV_VASC_STRATA", "PLOTSAMP"
  ))
  # Calculate mean number of vascular strata
  vstrat1.sum <- aggregate(
    x = list(XN_VASC_STRATA = vstratPlot$N_VSTRATA),
    by = vstratPlot[c("SAMPID", "N_VASC_STRATA", "XTOTCOV_VASC_STRATA", "NPLOTS", "PLOTSAMP")],
    FUN = sum
  )
  # Calculate minimum number of strata
  vstrat1.min <- aggregate(
    x = list(minN = vstratPlot$N_VSTRATA),
    by = vstratPlot[c("SAMPID", "N_VASC_STRATA", "XTOTCOV_VASC_STRATA", "NPLOTS", "PLOTSAMP")],
    FUN = min
  )
  # Calculate maximum number of strata
  vstrat1.max <- aggregate(
    x = list(maxN = vstratPlot$N_VSTRATA),
    by = vstratPlot[c("SAMPID", "N_VASC_STRATA", "XTOTCOV_VASC_STRATA", "NPLOTS", "PLOTSAMP")],
    FUN = max
  )
  # Merge data sets together
  vstrat1 <- merge(vstrat1.sum, vstrat1.min,
    by = c("SAMPID", "N_VASC_STRATA", "XTOTCOV_VASC_STRATA", "NPLOTS", "PLOTSAMP")
  )
  vstrat1 <- merge(vstrat1, vstrat1.max, c("SAMPID", "N_VASC_STRATA", "XTOTCOV_VASC_STRATA", "NPLOTS", "PLOTSAMP"))
  # Final calculations for mean and range of number of vascular strata
  vstrat1$XN_VASC_STRATA <- with(vstrat1, XN_VASC_STRATA / NPLOTS)
  vstrat1$RG_VASC_STRATA <- with(vstrat1, ifelse(PLOTSAMP == NPLOTS, maxN - minN, maxN - 0))
  # Drop unnecessary variables
  vstrat1$minN <- NULL
  vstrat1$maxN <- NULL
  vstrat1$PLOTSAMP <- NULL
  vstrat1$NPLOTS <- NULL
  # Create empty data frame
  empty_vstrat <- data.frame(t(rep(NA, 4)), stringsAsFactors = FALSE)
  names(empty_vstrat) <- c("N_VASC_STRATA", "XTOTCOV_VASC_STRATA", "XN_VASC_STRATA", "RG_VASC_STRATA")
  # Merge empty data frame output data frame from above, drop row with missing SAMPID
  vstrat2 <- subset(merge(vstrat1, empty_vstrat, all = TRUE), !is.na(SAMPID))
  # Merge full set of samples with output to ensure all samples represented
  vstrat3 <- merge(allUIDs, vstrat2, by = "SAMPID", all.x = T)
  # Melt data frame, then replace missing RESULT with 0
  varNames <- names(vstrat3)[!names(vstrat3) %in% c("SAMPID")]
  vstratMet <- reshape(vstrat3,
    idvar = "SAMPID", direction = "long",
    varying = varNames, times = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )

  vstratMet$METRIC <- with(vstratMet, as.character(METRIC))
  vstratMet$RESULT <- with(vstratMet, ifelse(is.na(RESULT), 0, RESULT))

  ## Calculate frequency, mean cover, relative mean cover, and relative importance by vascular stratum
  vstrat.pos <- subset(vstrat, RESULT != 0)
  # Calculate precursor for FREQ calculation by parameter
  indf1.length <- aggregate(
    x = list(FREQ = vstrat.pos$PLOT),
    by = vstrat.pos[c("SAMPID", "PARAMETER", "NPLOTS", "XTOTCOV_VASC_STRATA")],
    FUN = length
  )
  # Sum cover by parameter
  indf1.sum <- aggregate(
    x = list(XCOV = vstrat.pos$RESULT),
    by = vstrat.pos[c("SAMPID", "PARAMETER", "NPLOTS", "XTOTCOV_VASC_STRATA")],
    FUN = function(x) {
      sum(as.numeric(x))
    }
  )
  # Merge calculated values together, then calculate additional metrics
  indf1 <- merge(indf1.length, indf1.sum, by = c("SAMPID", "PARAMETER", "NPLOTS", "XTOTCOV_VASC_STRATA"), all = TRUE)
  indf1$FREQ <- with(indf1, round(FREQ / NPLOTS * 100, 2))
  indf1$XCOV <- with(indf1, round(XCOV / NPLOTS, 2))
  indf1$XRCOV <- with(indf1, round(XCOV / XTOTCOV_VASC_STRATA * 100, 2))
  indf1$IMP <- with(indf1, round((FREQ + XCOV) / 2, 2))

  # Now drop variables
  indf1$XTOTCOV_VASC_STRATA <- NULL
  indf1$NPLOTS <- NULL
  # Melt data frame
  outdf <- reshape(indf1,
    idvar = c("SAMPID", "PARAMETER"), direction = "long",
    varying = c("FREQ", "XCOV", "XRCOV", "IMP"), times = c("FREQ", "XCOV", "XRCOV", "IMP"),
    timevar = "variable", v.names = "RESULT"
  )
  # Update parameter value by combining variable from above with PARAMETER
  outdf$PARAMETER <- with(outdf, paste(variable, PARAMETER, sep = "_"))
  # Drop variable
  outdf$variable <- NULL
  # Cast output data frame to wide format and drop prefix added to variable names by reshape()
  outdf.wide <- reshape(outdf,
    idvar = c("SAMPID"), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(outdf.wide) <- gsub("RESULT\\.", "", names(outdf.wide))
  # Melt again, now that we have all metrics represented for all samples
  varNames <- names(outdf.wide)[!names(outdf.wide) %in% c("SAMPID")]
  outdf1 <- reshape(outdf.wide,
    idvar = c("SAMPID"), direction = "long",
    varying = varNames, times = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )
  # Set RESULT to 0 where missing
  outdf1$RESULT[is.na(outdf1$RESULT)] <- 0

  # Calculate diversity indices based on vascular strata
  # Start with Shannon diversity (H)
  div1.h <- aggregate(
    x = list(H_VASC_STRATA = indf1$XRCOV),
    by = indf1[c("SAMPID")], FUN = function(x) {
      round(-1 * sum((x / 100) * log(x / 100)), 4)
    }
  )
  # Now Pielou's evenness - initial step in calculation
  div1.j <- aggregate(
    x = list(J_VASC_STRATA = indf1$PARAMETER),
    by = indf1[c("SAMPID")], FUN = length
  )
  # Simpson Diversity index
  div1.d <- aggregate(
    x = list(D_VASC_STRATA = indf1$XRCOV),
    by = indf1[c("SAMPID")],
    FUN = function(x) {
      round(1 - sum((x / 100)^2), 4)
    }
  )
  # Merge data frames
  div1 <- merge(div1.h, div1.j, by = "SAMPID", all = TRUE)
  div1 <- merge(div1, div1.d, by = "SAMPID", all = TRUE)
  # Finish calculation of J, then fill in 0 where missing
  div1$J_VASC_STRATA <- with(div1, round(H_VASC_STRATA / log(J_VASC_STRATA), 4))
  div1$J_VASC_STRATA <- with(div1, ifelse(!is.na(J_VASC_STRATA), J_VASC_STRATA, 0))
  # Melt data frame
  div1.long <- reshape(div1,
    idvar = "SAMPID", direction = "long",
    varying = c("H_VASC_STRATA", "J_VASC_STRATA", "D_VASC_STRATA"),
    times = c("H_VASC_STRATA", "J_VASC_STRATA", "D_VASC_STRATA"),
    timevar = "METRIC", v.names = "RESULT"
  )
  # Combine with previously calculated metrics above
  outdf2 <- rbind(outdf1, div1.long)
  # Create empty data frame
  empty_vtype <- data.frame(t(rep(NA, 39)), stringsAsFactors = FALSE)
  names(empty_vtype) <- c(
    "FREQ_FLOATING_AQ", "FREQ_HMED_VEG", "FREQ_LIANAS", "FREQ_MED_VEG",
    "FREQ_SMALL_VEG", "FREQ_SUBMERGED_AQ",
    "FREQ_TALL_VEG", "FREQ_VSMALL_VEG", "FREQ_VTALL_VEG",
    "IMP_FLOATING_AQ", "IMP_HMED_VEG", "IMP_LIANAS",
    "IMP_MED_VEG", "IMP_SMALL_VEG", "IMP_SUBMERGED_AQ",
    "IMP_TALL_VEG", "IMP_VSMALL_VEG", "IMP_VTALL_VEG",
    "XCOV_FLOATING_AQ", "XCOV_HMED_VEG", "XCOV_LIANAS",
    "XCOV_MED_VEG", "XCOV_SMALL_VEG", "XCOV_SUBMERGED_AQ",
    "XCOV_TALL_VEG", "XCOV_VSMALL_VEG", "XCOV_VTALL_VEG",
    "XRCOV_FLOATING_AQ", "XRCOV_HMED_VEG", "XRCOV_LIANAS",
    "XRCOV_MED_VEG", "XRCOV_SMALL_VEG", "XRCOV_SUBMERGED_AQ",
    "XRCOV_TALL_VEG", "XRCOV_VSMALL_VEG", "XRCOV_VTALL_VEG",
    "H_VASC_STRATA", "J_VASC_STRATA", "D_VASC_STRATA"
  )
  # Cast output data frame from above and drop prefix added by reshape()
  outdf3 <- reshape(outdf2,
    idvar = "SAMPID", direction = "wide",
    timevar = "METRIC", v.names = "RESULT"
  )

  names(outdf3) <- gsub("RESULT\\.", "", names(outdf3))
  # Merge with empty data frame and drop row with missing SAMPID
  outdf4 <- subset(merge(outdf3, empty_vtype, all = TRUE), !is.na(SAMPID))
  # Merge with full list of samples to make sure all are represented
  outdf5 <- merge(allUIDs, outdf4, by = "SAMPID", all.x = T)
  # Melt output data frame and fill in 0 for appropriate metrics where missing values
  varNames <- names(outdf5)[!names(outdf5) %in% c("SAMPID")]
  outdf6 <- reshape(outdf5,
    idvar = "SAMPID", direction = "long",
    varying = varNames, times = varNames,
    timevar = "METRIC", v.names = "RESULT"
  )
  outdf6$METRIC <- with(outdf6, as.character(METRIC))
  outdf6$RESULT <- with(outdf6, ifelse(METRIC %nin% c("D_VASC_STRATA", "H_VASC_STRATA") & is.na(RESULT), 0, RESULT))

  print("Done with vascular strata metrics")

  # Now combine output data frames, rename METRIC to PARAMETER and drop METRIC
  vtOut <- rbind(outdf6, vstratMet)
  vtOut$PARAMETER <- vtOut$METRIC
  vtOut$METRIC <- NULL
  # Merge with samples to get back to sampID variables
  vtOut.1 <- merge(samples, vtOut, by = "SAMPID")
  # Select subset of variables
  vtOut.1 <- vtOut.1[, c(sampID, "PARAMETER", "RESULT")]

  return(vtOut.1)
}
