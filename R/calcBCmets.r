#' @export
#'
#' @title Calculate Bray-Curtis metrics
#'
#' @description This function calculates the mean Bray-Curtis
#' distances among plots for all species, and includes a version using only
#' native species if the variable NWCA_NATSTAT is included in the input
#' data frame. This variable is found in the ccNatNWCA dataset.
#'
#' @param vascIn Data frame containing cover data summarized by
#' \emph{sampID} variables, PLOT, and
#' DISTINCT at the species level. Must also contain at least one of the
#' following: TAXON (taxon name) or SPECIES_NAME_ID (numeric code
#' for taxon). If NWCA_NATSTAT is included, a native species version is
#' also calculated.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#'
#' @details This function calculates metrics based on all species, and on
#' native species if the appropriate variable is included in the
#' input data frame.
#'
#' @return Data frame containing \emph{sampID} variables, PARAMETER, and RESULT,
#'   with one row of results per parameter and UID. The values for PARAMETER
#'   consist of the metric name concatenated with taxonomic level (represented
#'   as TAXLEVEL below):
#' \itemize{
#' \item XBCDIST_SPP: Mean Bray-Curtis distance among plots in sample
#'
#' \item XBCDIST_NATSPP: Mean Bray-Curtis distance among plots in sample
#' based only on native species.
#' }

calcBCmets <- function(vascIn, sampID = "UID") {
  vascIn <- as.data.frame(vascIn) # Do this in case read in as a tibble or data.table, which might cause problems
  # Need to account for cases where there is no SPECIES_NAME_ID by creating one just for this
  # calculation
  if ("SPECIES_NAME_ID" %nin% names(vascIn)) {
    uniqNames <- data.frame(TAXON = unique(vascIn$TAXON), stringsAsFactors = F)
    uniqNames$SPECIES_NAME_ID <- seq(from = 1, to = nrow(uniqNames))
    vascIn <- merge(vascIn, uniqNames, by = "TAXON")
  }
  vascIn$COVER <- as.numeric(vascIn$COVER)

  forDist <- aggregate(
    x = list(COVER = vascIn$COVER),
    by = vascIn[, c(sampID, "PLOT", "SPECIES_NAME_ID")],
    FUN = sum
  )
  # This df needs to be in wide format
  forDist$SPECIES <- with(forDist, paste("s", SPECIES_NAME_ID, sep = ""))

  meanBC <- int.calcXBC(forDist, sampID)

  xbcOut <- meanBC

  if ("NWCA_NATSTAT" %in% names(vascIn)) {
    forDist.nat <- aggregate(
      x = list(COVER = vascIn$COVER),
      by = vascIn[, c(sampID, "PLOT", "SPECIES_NAME_ID", "NWCA_NATSTAT")],
      FUN = sum
    )

    forDist.nat$SPECIES <- with(forDist.nat, paste("s", SPECIES_NAME_ID, sep = ""))
    forDist.nat <- subset(forDist.nat, NWCA_NATSTAT == "NAT")

    meanBC_nat <- int.calcXBC(forDist.nat, sampID)
    meanBC_nat$XBCDIST_NATSPP <- meanBC_nat$XBCDIST_SPP
    meanBC_nat$XBCDIST_SPP <- NULL

    xbcOut <- merge(xbcOut, meanBC_nat, by = sampID, all.x = T)
  }

  xbcOut.1 <- reshape(xbcOut,
    idvar = sampID, direction = "long",
    varying = names(xbcOut)[!names(xbcOut) %in% c(sampID)],
    timevar = "PARAMETER", v.names = "RESULT",
    times = names(xbcOut)[!names(xbcOut) %in% c(sampID)]
  )

  xbcOut.1$RESULT <- with(xbcOut.1, ifelse(is.na(RESULT), 0, RESULT))
  xbcOut.1$PARAMETER <- as.character(xbcOut.1$PARAMETER)

  return(xbcOut.1)
}
