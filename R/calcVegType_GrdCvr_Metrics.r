#' @export
#' @title Calculate NWCA 2011 vegetation types and ground surface
#' attribute metrics based on data collected from field form V-3
#' @description This function calculates the NWCA 2011 vegetation types and
#' ground surface attribute metrics.  It assumes input data are organized
#' so that each row represents the cover for a single taxon within a single
#' plot at a site. It also assumes a single column is used to identify
#' each unique sample. Missing values and blanks should be removed from the
#' input data frame.
#' @param dataIn A data frame containing the following variables:
#' \itemize{
#' \item sampID - Variables identified in the sampID argument
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
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @details If any of the parameters are missing, they are assumed to be
#' zeros (if numeric), and metric values associated with any metrics that
#' cannot be calculated due to missing parameters are set to a standardized
#' value.
#' @return Either a character string containing an error message when metric
#' calculation is not successful, or a data frame. The first column of the
#' data frame is the sampID variables and subsequent columns are named for each metric and
#' contain metric values. A list of metrics is provided in the document named
#' "VegTypes_GrdCover_Metric_Descriptions.pdf" included in the help directory
#' for the package.
#' @references US Environmental Protection Agency. 2016. National
#' Wetland Condition Assessment: 2011 Technical Report. EPA-843-R-15-006.
#' US Environmental Protection Agency, Washington, DC.
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' @examples
#' head(Vtype_GrCovEx)
#' exOut <- calcVtype_GcovMets(Vtype_GrCovEx)
#'
#' str(exOut)
#' head(exOut)

calcVtype_GcovMets <- function(dataIn,sampID='UID'){

  datNames <- c('UID','PLOT','PARAMETER','RESULT')
  if(any(datNames %nin% names(dataIn))){
    print(paste("Missing key variables! Should be ",paste(sampID,collapse=', '),"PLOT, PARAMETER, and RESULT.",sep=''))
    return(NULL)
  }

  nplots <- ddply(subset(dataIn,PARAMETER!='SANDT_CLASS'),sampID,summarise,NPLOTS=length(unique(PLOT)))

  sandtOut <- calcSandTMets(dataIn,nplots,sampID)
  vstratOut <- calcVascStratMets(dataIn,nplots,sampID)
  nonvascOut <- calcNonvascMets(dataIn,nplots,sampID)
  wcovOut <- calcWcovMets(dataIn,nplots,sampID)
  bgLittOut <- calcBareGround_LitterMets(dataIn,nplots,sampID)


  vgOut <- rbind(sandtOut,vstratOut,nonvascOut,wcovOut,bgLittOut)
  
  formula <- paste(paste(sampID,collapse='+'),'~PARAMETER',sep='')
  vgOut.1 <- reshape2::dcast(vgOut,eval(formula),value.var='RESULT')

  return(vgOut.1)

}
