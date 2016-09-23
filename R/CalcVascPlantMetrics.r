#' @export
#' @title Calculate NWCA 2011 vascular plant metrics based on form V-2
#' @description This function calculates the NWCA 2011 vascular plant metrics.
#' It assumes input data are organized so that each row represents the cover
#' for a single taxon within a single plot at a site. It also assumes a single
#' column is used to identify each unique sample.
#' @param dfIn A data frame with the following variables (at a minimum):
#' \itemize{
#' \item sampID - A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#'
#' \item PLOT (numeric variable)
#'
#' \item STATE (state where sample collected)
#'
#' \item USAC_REGION (U.S. Army Corps of Engineers region of site where sample
#' collected)
#'
#' \item USDA_NAME (USDA accepted name of taxon in sample)
#'
#' \item COVER (value from 0 to 100 indicating relative cover of taxon within plot).
#' Data should already be summed by UID, PLOT, and USDA_NAME.
#' }
#' @param taxaIn A data frame containing taxonomy of all taxa found in dfIn,
#' with the following variables at a minimum:
#' \itemize{
#' \item USDA_NAME (USDA accepted name of taxon in sample)
#'
#' \item FAMILY (family name of taxon)
#'
#' \item GENUS (genus name of taxon)
#'
#' \item CATEGORY (USDA PLANTS category)
#'
#' \item GROWTH_HABIT (growth habit as designated in USDA PLANTS)
#'
#' \item DURATION (duration as designated by USDA PLANTS).
#' }
#' If this data frame is not supplied, the data set taxaNWCA included
#' in the package is used.
#' @param taxaCC A data frame containing Coefficient of Conservatism and native
#' status as assigned for NWCA for all taxa in dfIn, with the following
#' variables, at a minimum:
#' \itemize{
#' \item USDA_NAME (USDA accepted name)
#'
#' \item GEOG_ID (state to which value applies)
#'
#' \item NWCA_CC (state- and taxon-specific Coefficient of Conservatism)
#'
#' \item NWCA_NATSTAT (state- and taxon-specific native status)
#' }
#' If this data frame is not supplied, the data set ccNWCA included in
#' the package is used.
#'
#' @param taxaWIS A data frame containing Wetland Indicator Status as assigned
#' by USAC and reconciled by USDA PLANTS, with the following variables, at a
#' minimum:
#' \itemize{
#' \item USDA_NAME (USDA accepted name)
#'
#' \item GEOG_ID (USAC region to which value applies)
#'
#' \item WIS (Wetland Indicator Status)
#' }
#' If this data frame is not supplied, the data set wisNWCA included in the
#' package is used.
#' @return Either a character string containing an error message when metric
#' calculation is not successful, or a data frame. The first column of the
#' data frame is the UID and subsequent columns are named for each metric and
#' contain metric values. A list of metrics is provided in the document named
#' "Vascular_Plant_Metric_Descriptions.pdf" included in the help directory for
#' the package.
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC.
#'
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' @examples
#' head(VascPlantEx)
#'
#' exOut <- calcVascPlantMets(VascPlantEx,taxaNWCA,ccNatNWCA,wisNWCA)
#'
#' str(exOut)
#' head(exOut)
calcVascPlantMets <- function(dfIn,taxaIn=taxaNWCA,taxaCC=ccNatNWCA,taxaWIS=wisNWCA,sampID='UID'){
  prepDat <- prepareData(dfIn,taxaIn,taxaCC,taxaWIS)

  print("Initial datasets prepared for metric calculation.")

  richMets <- calcRichness(prepDat$byUIDspp,prepDat$byPlotspp,prepDat$byUIDgen,prepDat$byPlotgen
                           ,prepDat$byUIDfam,prepDat$byPlotfam)

  print("Done calculating richness metrics.")


  sppForCalc <- prepDat$byUIDspp

  print("Ready to start calculating trait metrics.")

  # Calculate diversity indices
  divMets <- calcDiversity(sppForCalc,sampID)
  print("Done with diversity indices.")

  durMets <- calcDuration(sppForCalc,sampID)
  print("Done with duration metrics.")

  grhMets <- calcGrowthHabit(sppForCalc,sampID)
  print("Done with growth habit metrics.")

  catMets <- calcCategory(sppForCalc,sampID)
  print("Done with category metrics.")

  wisMets <- calcWIS(sppForCalc,sampID)
  print("Done with WIS metrics.")

  ccMets <- calcCC(sppForCalc,sampID)
  print("Done with CC metrics.")

  natMets <- calcNative(sppForCalc,sampID)
  print("Done with native status metrics.")

  xbcMets <- calcBCmets(prepDat$byPlotspp,sampID)
  print("Done with Bray-Curtis distance.")

  # Now create df with just XTOTABCOV by UID
  xtotabcov <- reshape2::melt(unique(subset(sppForCalc,select=c(sampID,'XTOTABCOV'))),id.vars='UID'
                              ,variable.name='PARAMETER',value.name='RESULT')

  ### COMBINE ALL PLANT METRICS INTO A SINGLE DF
  allOut <- rbind(xtotabcov,richMets,divMets,durMets,grhMets,catMets,wisMets,ccMets,natMets,xbcMets)

  formula <- paste(paste(sampID,collapse='+'),'~PARAMETER',sep='')
  finOut <- dcast(allOut,eval(formula),value.var='RESULT')
  return(finOut)
}

