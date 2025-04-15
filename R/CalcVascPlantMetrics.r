#' @export
#'
#' @title Calculate NWCA 2011 vascular plant metrics based on form V-2
#'
#' @description This function calculates the NWCA 2011 vascular plant metrics.
#'   It assumes input data are organized so that each row represents the cover
#'   for a single taxon within a single plot at a site. It also assumes a single
#'   column is used to identify each unique sample.
#'
#' @param vascIn A data frame with the following variables (at a minimum):
#'   \itemize{ \item sampID: Variable(s) identified in \emph{sampID} argument,
#'   'UID' by default
#'
#'   \item PLOT: plot number
#'
#'   \item Variable named in \emph{state}: Two-letter state postal code of site, used to link
#' native status to taxa in native status taxalist (inNat)
#'
#'   \item Variable named in \emph{coeReg}: U.S. Army Corps of Engineers region abbreviation for
#'   sample, to correspond to GEOG_ID in Wetland Indicator Status taxalist (inWIS)
#'
#'   \item Variable named in \emph{cValReg}: NWCA C-value regions: values must match GEOG_ID
#'   in C-value taxalist (inCVal)
#'  \item Variable named in \emph{taxon_name}: Taxon name, must match with taxa data frame
#'
#'   \item COVER: value from 0 to 100 indicating relative cover of taxon within
#'   plot). Data should already be summed by UID, PLOT, and USDA_NAME. }
#' @param taxon_name String containing the name of variable for taxon name in
#' @param taxaIn A data frame containing taxonomy of all taxa found in vascIn,
#'   with the following variables at a minimum:
#'   \itemize{ \item USDA_NAME: USDA
#'   accepted name of taxon in sample
#'
#'   \item FAMILY: family name of taxon
#'
#'   \item GENUS: genus name of taxon
#'
#' \item CATEGORY: USDA PLANTS category
#'
#' \item GROWTH_HABIT: growth habit as designated in USDA PLANTS
#'
#' \item DURATION: duration as designated by USDA PLANTS
#' }
#' If this data frame is not supplied, the data set taxaNWCA included
#' in the package is used.
#' @param taxaNat A data frame containing native
#' status as assigned for NWCA for all taxa in vascIn, with the following
#' variables, at a minimum:
#' \itemize{
#' \item USDA_NAME: USDA accepted name
#'
#' \item GEOG_ID: state to which value applies
#'
#' \item NWCA_NATSTAT: state- and taxon-specific native status
#' }
#' If this data frame is not supplied, the data set ccNWCA included in
#' the package is used.
#' @param taxaCC A data frame containing Coefficient of Conservatism
#' as assigned for NWCA for all taxa in vascIn, with the following
#' variables, at a minimum:
#' \itemize{
#' \item USDA_NAME: USDA accepted name
#'
#' \item GEOG_ID: \emph{cValReg} to which value applies
#'
#' \item NWCA_CC: state- and taxon-specific Coefficient of Conservatism
#' }
#' If this data frame is not supplied, the data set ccNWCA included in
#' the package is used.
#' @param taxaWIS A data frame containing Wetland Indicator Status as assigned
#' by USAC and reconciled by USDA PLANTS, with the following variables, at a
#' minimum:
#' \itemize{
#' \item USDA_NAME: USDA accepted name
#'
#' \item GEOG_ID: USAC region to which value applies
#'
#' \item WIS: Wetland Indicator Status
#' }
#' If this data frame is not supplied, the data set wisNWCA included in the
#' package is used.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#' @param state String containing the name of the state in \emph{vascIn},
#' with default value of 'STATE'
#'
#' @param coeReg String containing the name of the U.S. Army Corps of Engineers
#' region in \emph{vascIn} associated with Wetland Indicator Status,
#' with default value of 'USAC_REGION'
#'
#' @param cValReg String containing the name of the variable in \emph{taxaCC}
#'  which specifies the C-value region.
#' @return Either a character string containing an error message when metric
#' calculation is not successful, or a data frame. The first column(s) of the
#' data frame contain the \emph{sampID} variables, and subsequent columns are
#' named for each metric and
#' contain metric values. 
#' A list of metric descriptions is provided in the document named 
#' \href{https://github.com/USEPA/aquametNWCA/blob/main/inst/VascPlant_Metric_Descriptions.pdf}{VascPlant_Metric_Descriptions.pdf}
#'
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC. \url{https://www.epa.gov/national-aquatic-resource-surveys/national-wetland-condition-assessment-2011-technical-report}
#' @references US Environmental Protection Agency. 2023. National Wetland 
#' Condition Assessment: 2016 Technical Support Document. EPA-841-B-23-001. 
#' \url{https://www.epa.gov/national-aquatic-resource-surveys/national-wetland-condition-assessment-2016-technical-support}
#'
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#'
#' @examples
#' head(VascPlantEx)
#'
#' exOut <- calcVascPlantMets(VascPlantEx,
#'   taxon_name = "USDA_NAME",
#'   taxaIn = taxaNWCA, taxaNat = ccNatNWCA,
#'   taxaCC = ccNatNWCA, taxaWIS = wisNWCA, cValReg = "STATE"
#' )
#'
#' str(exOut)
#' head(exOut)
calcVascPlantMets <- function(vascIn, taxon_name, taxaIn = taxaNWCA, taxaNat = ccNatNWCA, taxaCC = ccNatNWCA,
                              taxaWIS = wisNWCA, sampID = "UID", state = "STATE",
                              coeReg = "USAC_REGION", cValReg = "NWC_CREG") {
  print(cValReg)
  # Organize the vascular plant data in such a way that it is readily useable in the metric
  # calculation functions
  prepDat <- prepareData(vascIn, sampID,
    taxon_name = taxon_name, inTaxa = taxaIn, inNat = taxaNat,
    inCVal = taxaCC, inWIS = taxaWIS, state = state,
    coeReg = coeReg, cValReg = cValReg
  )

  print("Initial datasets prepared for metric calculation.")

  # Calculate richness metrics - requires all of the various outputs from prepareData() above
  richMets <- calcRichness(prepDat$byUIDspp, prepDat$byPlotspp, prepDat$byUIDgen,
    prepDat$byPlotgen, prepDat$byUIDfam,
    prepDat$byPlotfam,
    sampID = sampID
  )

  print("Done calculating richness metrics.")

  # Set side the input that works with all other vascular metric functions
  sppForCalc <- prepDat$byUIDspp

  print("Ready to start calculating trait metrics.")

  # Calculate diversity indices
  divMets <- calcDiversity(sppForCalc, sampID)
  print("Done with diversity indices.")

  # Calculate duration metrics
  durMets <- calcDuration(sppForCalc, sampID)
  print("Done with duration metrics.")

  # Calculate growth habit metrics
  grhMets <- calcGrowthHabit(sppForCalc, sampID)
  print("Done with growth habit metrics.")

  # Calculate category metrics
  catMets <- calcCategory(sppForCalc, sampID)
  print("Done with category metrics.")

  # Calculate Wetland Indicator Status metrics
  wisMets <- calcWIS(sppForCalc, sampID)
  print("Done with WIS metrics.")

  # Calculate C-value metrics
  ccMets <- calcCC(sppForCalc, sampID)
  print("Done with CC metrics.")

  # Calculate native metrics
  natMets <- calcNative(sppForCalc, sampID)
  print("Done with native status metrics.")

  # Calculate mean Bray-Curtis metrics
  xbcMets <- calcBCmets(prepDat$byPlotspp, sampID)
  print("Done with Bray-Curtis distance.")

  # Now create df with just XTOTABCOV by UID
  sppForCalc.in <- unique(subset(sppForCalc, select = c(sampID, "XTOTABCOV")))

  # Melt this data frame in order to add it to other metric data
  varNames <- names(sppForCalc.in)[!names(sppForCalc.in) %in% c(sampID)]
  xtotabcov <- reshape(sppForCalc.in,
    idvar = sampID, direction = "long",
    varying = varNames, times = varNames,
    timevar = "PARAMETER", v.names = "RESULT"
  )

  # COMBINE ALL PLANT METRICS INTO A SINGLE DF
  allOut <- rbind(
    xtotabcov, richMets, divMets, durMets, grhMets, catMets,
    wisMets, ccMets, natMets, xbcMets
  )
  # Cast long data frame into wide format
  finOut <- reshape(allOut,
    idvar = c(sampID), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  # Must remove prefix added to names in reshape()
  names(finOut) <- gsub("RESULT\\.", "", names(finOut))

  return(finOut)
}
