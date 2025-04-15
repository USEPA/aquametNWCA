#' @export
#'
#' @title Calculate NWCA 2011 vegetation types and ground surface
#' attribute metrics based on data collected from field form V-3
#'
#' @description This function calculates the NWCA 2011 vegetation types and
#' ground surface attribute metrics.  It assumes input data are organized
#' so that each row represents the cover for a single taxon within a single
#' plot at a site. It also assumes a single column is used to identify
#' each unique sample. Missing values and blanks should be removed from the
#' input data frame.
#'
#' @param dataIn A data frame containing the following variables:
#' \itemize{
#' \item sampID: Variables identified in the \emph{sampID} argument
#'
#' \item PLOT: Sample plot from which data were collected
#'
#' \item PARAMETER: specific measurement type
#'
#' \item RESULT: measured value
#' }
#' The following parameters are used in
#' calculating tree metrics: 'SANDT_CLASS','SUBMERGED_AQ',
#' 'FLOATING_AQ','LIANAS','VTALL_VEG','TALL_VEG','HMED_VEG',
#' 'MED_VEG','SMALL_VEG','VSMALL_VEG','PEAT_MOSS','BRYOPHYTES',
#' 'LICHENS','ARBOREAL','ALGAE','MACROALGAE','TIME','MINIMUM_DEPTH',
#' 'MAXIMUM_DEPTH','PREDOMINANT_DEPTH','TOTAL_WATER','WATER_NOVEG',
#' 'WATER_AQVEG','WATER_EMERGVEG','LITTER_THATCH','LITTER_FORB',
#' 'LITTER_CONIFER','LITTER_DECID','LITTER_BROADLEAF',
#' 'LITTER_DEPTH_SW','LITTER_DEPTH_NE','TOTAL_LITTER','WD_FINE',
#' 'WD_COARSE','EXPOSED_SOIL','EXPOSED_GRAVEL','EXPOSED_ROCK',
#' 'EXPOSED_SOIL','EXPOSED_GRAVEL','EXPOSED_ROCK','WD_FINE',
#' 'WD_COARSE','TOTAL_LITTER'. Additional parameters or variables
#' are ignored.
#' @param nPlotIn A data frame with the
#' number of plots sampled associated with
#' each sample, with \emph{sampID} variables and NPLOTS.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#' @param survyear A string with the survey year. Some parameters only
#' measured in 2011, applies to water cover metrics and bare ground
#' and litter metrics.
#'
#' @details If any of the parameters are missing, they are assumed to be
#' zeros (if numeric), and metric values associated with any metrics that
#' cannot be calculated due to missing parameters are set to a standardized
#' value.
#'
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. The first column of the
#'   data frame is the \emph{sampID} variables and subsequent columns are named
#'   for each metric and contain metric values. A list of vegetation type and
#'   ground cover metrics is installed with the package. Use
#'   system.file("VegTypes_GrdCover_Metric_Descriptions.pdf", package="aquametNWCA")
#'   to locate the file or use this link: \href{https://github.com/USEPA/aquametNWCA/blob/main/inst/VType_GrdCvr_Metric_Descriptions.pdf}{VegTypes_GrdCover_Metric_Descriptions.pdf}
#'
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC. \url{https://www.epa.gov/national-aquatic-resource-surveys/national-wetland-condition-assessment-2011-technical-report}
#'
#' @references US Environmental Protection Agency. 2023. National Wetland 
#' Condition Assessment: 2016 Technical Support Document. EPA-841-B-23-001. 
#' \url{https://www.epa.gov/national-aquatic-resource-surveys/national-wetland-condition-assessment-2016-technical-support}
#' 
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#'
#' @examples
#' head(Vtype_GrCovEx)
#'
#' subVtype_GrCovEx <- subset(Vtype_GrCovEx, PARAMETER != "SANDT_CLASS")
#'
#' nplots <- aggregate(
#'   x = list(NPLOTS = subVtype_GrCovEx$PLOT),
#'   by = subVtype_GrCovEx[c("UID")],
#'   FUN = function(x) {
#'     length(unique(x))
#'   }
#' )
#'
#' exOut <- calcVtype_GcovMets(Vtype_GrCovEx, nPlotIn = nplots, sampID = "UID")
#'
#' str(exOut)
#' head(exOut)
calcVtype_GcovMets <- function(dataIn, nPlotIn, sampID = "UID", survyear = "2011") {
  # Evaluates whether treeIn contains expected variables, returns error if any are missing
  datNames <- c("UID", "PLOT", "PARAMETER", "RESULT")
  if (any(datNames %nin% names(dataIn))) {
    print(paste("Missing key variables! Should be ", paste(sampID, collapse = ", "),
      "PLOT, PARAMETER, and RESULT.",
      sep = ""
    ))
    return(NULL)
  }

  # Calculate metrics based on S & T class
  sandtOut <- calcSandTMets(dataIn, nPlotIn, sampID)
  # Calculate vascular strata metrics
  vstratOut <- calcVascStratMets(dataIn, nPlotIn, sampID)
  # Calculate non-vascular metrics
  nonvascOut <- calcNonvascMets(dataIn, nPlotIn, sampID)
  # Calculate water cover metrics
  wcovOut <- calcWcovMets(dataIn, nPlotIn, sampID, survyear)
  # Calculate bare ground and litter metrics
  bgLittOut <- calcBareGround_LitterMets(dataIn, nPlotIn, sampID, survyear)

  # Combine all metrics into single data frame
  vgOut <- rbind(sandtOut, vstratOut, nonvascOut, wcovOut, bgLittOut)
  # Cast data frame into wide format
  vgOut.1 <- reshape(vgOut,
    idvar = sampID, direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  # Remove prefix added by reshape() from variable names
  names(vgOut.1) <- gsub("RESULT\\.", "", names(vgOut.1))

  return(vgOut.1)
}
