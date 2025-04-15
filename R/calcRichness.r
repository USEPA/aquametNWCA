#' @export
#'
#' @title Calculate richness metrics
#'
#' @description This function calculates richness metrics using plot- and
#' UID-based datasets at the species, genus, and family levels.
#'
#' @param byUIDspp Data frame containing species data summarized to the UID
#' and TAXON level, with the following variables:
#' \itemize{
#'  \item sampID: Variables identified by \emph{sampID} argument
#'
#' \item TAXON: Taxon name
#'
#' \item COVER: Total percent cover for taxon
#'
#' \item DISTINCT: Distinctness (0/1) of taxon within sample, where
#' taxa above species level are not distinct (0) if taxa
#' at a lower level are also included in sample.
#' }
#' @param byPlotspp Data frame containing species data summarized
#' to the UID, PLOT, and TAXON level, with the following variables:
#' \itemize{
#'  \item sampID: Variables identified in sampID argument
#'
#' \item PLOT: Plot from which data were collected
#'
#' \item TAXON: Taxon name
#'
#' \item COVER: Total percent cover for taxon
#'
#' \item DISTINCT: Distinctness (0/1) of taxon within sample, where
#' taxa above species level are not distinct (0) if taxa
#' at a lower level are also included in sample.
#' }
#' @param byUIDgen Data frame containing genus-level data summarized to the
#'   sampID variables and TAXON level and containing the same variables as
#'   byUIDspp (described above).
#' @param byPlotgen Data frame containing genus-level data summarized to the
#'   sampID variables, PLOT, and TAXON level and containing the same variables
#'   as byPlotspp (described above).
#' @param byUIDfam  Data frame containing family-level data summarized to the
#'   sampID variables and TAXON level and containing the same variables as
#'   byUIDspp (described above).
#' @param byPlotfam Data frame containing family-level data summarized to the
#'   sampID variables, PLOT, and TAXON level and containing the same variables
#'   as byPlotspp (described above).
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#'
#' @details The prepareData() function creates a list object with
#' all of the necessary input data frames. For each taxonomic level,
#' the function createDFs() creates a list with a data frame summarized
#' by UID and one by sampID variables and PLOT.
#'
#' @return   Data frame containing \emph{sampID} variables, PARAMETER, and
#'   RESULT, with one row of results per parameter and \emph{sampID}. The values
#'   for PARAMETER consist of the metric name concatenated with taxonomic level
#'   (represented as SPP, GEN, and FAM below):
#' \itemize{
#' \item TOTN_SPP, TOTN_GEN, TOTN_FAM: Number of unique taxa in sample (UID)
#'
#' \item XN_SPP, XN_GEN, XN_FAM: Mean number of taxa per plot
#'
#' \item MEDN_SPP, MEDN_GEN, MEDN_FAM: Median number of taxa per plot
#'
#' \item SDN_SPP, SDN_GEN, SDN_FAM: Standard deviation of number of taxa per plot
#'
#' \item N_PLOTS: Number of plots sampled for UID
#' }
#'
#' If NWCA_NATSTAT is present in the input data frame, the following
#' metrics are calculated. GRP in the metric names below represents
#' native status values of NAT, ADV, CRYP, INTR, and subsets ALIEN
#' (INTR + ADV) and AC (ALIEN and CRYP):
#' \itemize{
#' \item TOTN_GRP: Number of unique taxa in sample (UID)
#'
#' \item XN_GRP: Mean number of taxa per plot
#'
#' \item MEDN_GRP: Median number of taxa per plot
#'
#' \item SDN_GRP: Standard deviation of number of taxa per plot
#' }
#' A list of metric descriptions is provided in the document named 
#' \href{https://github.com/USEPA/aquametNWCA/blob/main/inst/VascPlant_Metric_Descriptions.pdf}{VascPlant_Metric_Descriptions.pdf}
#' @references US Environmental Protection Agency. 2016. National
#' Wetland Condition Assessment: 2011 Technical Report. EPA-843-R-15-006.
#' US Environmental Protection Agency, Washington, DC.
#'
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#'
#' @examples
#' head(VascPlantEx)
#' exPlant <- prepareData(VascPlantEx,
#'   taxon_name = "USDA_NAME",
#'   inTaxa = taxaNWCA, inNat = ccNatNWCA, inCVal = ccNatNWCA,
#'   inWIS = wisNWCA, cValReg = "STATE"
#' )
#'
#' richEx <- calcRichness(
#'   exPlant$byUIDspp, exPlant$byPlotspp,
#'   exPlant$byUIDgen, exPlant$byPlotgen, exPlant$byUIDfam, exPlant$byPlotfam
#' )
#'
#' head(richEx)
#' unique(richEx$PARAMETER)
calcRichness <- function(byUIDspp, byPlotspp, byUIDgen, byPlotgen,
                         byUIDfam, byPlotfam, sampID = "UID") {
  sppRich <- int.calcRich(byUIDspp, byPlotspp, "SPP", sampID)
  genRich <- int.calcRich(byUIDgen, byPlotgen, "GEN", sampID)
  famRich <- int.calcRich(byUIDfam, byPlotfam, "FAM", sampID)

  richOut <- rbind(sppRich, genRich, famRich)

  if ("NWCA_NATSTAT" %in% names(byUIDspp)) {
    natRich <- int.calcRichNS(byUIDspp, byPlotspp, c("NAT"), "NATSPP", sampID)
    advRich <- int.calcRichNS(byUIDspp, byPlotspp, c("ADV"), "ADVSPP", sampID)
    crypRich <- int.calcRichNS(byUIDspp, byPlotspp, c("CRYP"), "CRYPSPP", sampID)
    intrRich <- int.calcRichNS(byUIDspp, byPlotspp, c("INTR"), "INTRSPP", sampID)
    alienRich <- int.calcRichNS(byUIDspp, byPlotspp, c("ADV", "INTR"), "ALIENSPP", sampID)
    acRich <- int.calcRichNS(byUIDspp, byPlotspp, c("INTR", "ADV", "CRYP"), "AC", sampID)

    # Combine all into a single df
    allNSrich <- rbind(natRich, advRich, crypRich, intrRich, alienRich, acRich)

    richOut <- rbind(richOut, allNSrich)
  }
  # Must fill in missing categories for all sites with zeros

  richOut.wide <- reshape(richOut,
    idvar = c(sampID), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(richOut.wide) <- gsub("RESULT\\.", "", names(richOut.wide))

  outdf <- reshape(richOut.wide,
    idvar = sampID, direction = "long",
    varying = names(richOut.wide)[!names(richOut.wide) %in% c(sampID)],
    timevar = "PARAMETER", v.names = "RESULT",
    times = names(richOut.wide)[!names(richOut.wide) %in% c(sampID)]
  )

  outdf$RESULT <- with(outdf, ifelse(is.na(RESULT), 0, RESULT))
  outdf$PARAMETER <- as.character(outdf$PARAMETER)

  return(outdf)
}
