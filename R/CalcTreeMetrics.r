#' @export
#' 
#' @title Calculate NWCA 2011 tree metrics based on data collected from
#' field form V-4a
#' 
#' @description This function calculates the NWCA 2011 tree metrics.
#' It assumes input data are organized in long format so that each
#' row represents the value for a single parameter value within a single
#' plot at a site.
#' 
#' @param treeIn A data frame containing the following variables:
#' \itemize{
#'    \item sampID: Variable(s) found in the argument \emph{sampID}
#'
#'    \item PAGE: Page number from field form
#'
#'    \item LINE: Line number from field form
#'
#'    \item PLOT: Sample plot from which data were collected
#'
#'    \item PARAMETER: specific measurement type
#'
#'    \item RESULT: measured value
#' }
#' The following parameters are used in
#' calculating tree metrics: 'XXTHIN_SNAG','XTHIN_SNAG','THIN_SNAG',
#' 'JR_SNAG','THICK_SNAG','XTHICK_SNAG','XXTHICK_SNAG','XXTHIN_TREE',
#' 'XTHIN_TREE','THIN_TREE','JR_TREE','THICK_TREE','XTHICK_TREE',
#' 'XXTHICK_TREE','VSMALL_TREE','SMALL_TREE','LMED_TREE','HMED_TREE',
#' 'TALL_TREE','VTALL_TREE'. Additional parameters or variables are
#' ignored.
#' @param nPlotIn A data frame with the 
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
#'   calculation is not successful, or a data frame. The first columns of the
#'   data frame are the sampID variables and subsequent columns are named for
#'   each metric and contain metric values. A list of tree metrics is 
#'   installed with the package. Use 
#'   system.file("Tree_Metric_Descriptions.pdf", package="aquametNWCA") to 
#'   locate the file. 
#'   
#' @references US Environmental Protection Agency. 2016. National
#' Wetland Condition Assessment: 2011 Technical Report. EPA-843-R-15-006.
#' US Environmental Protection Agency, Washington, DC.
#' 
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(TreesEx)
#' 
#' nplots <- plyr::ddply(TreesEx,c('UID'),dplyr::summarise
#' ,NPLOTS=length(unique(PLOT)))
#'
#' exOut <- calcTreeMets(TreesEx, nPlotIn=nplots, sampID='UID')
#'
#' head(exOut)
#' str(exOut)

calcTreeMets <- function(treeIn, nPlotIn, sampID='UID'){

  datNames <- c(sampID,'PAGE','LINE','PLOT','PARAMETER','RESULT')
  if(any(datNames %nin% names(treeIn))){
    print("Missing key variables! Should be UID, PAGE, LINE, PLOT, PARAMETER, and RESULT.")
    return(NULL)
  }

  snagsOut <- calcSnagMets(treeIn, nPlot=nPlotIn, sampID)
  treeCntsOut <- calcTreeCntMets(treeIn, nPlot=nPlotIn, sampID)
  treeCvrOut <- calcTreeCoverMets(treeIn, nPlot=nPlotIn, sampID)

  ## Merge into a single output file in wide format
  treeMetOut <- rbind(snagsOut,treeCntsOut,treeCvrOut)
  
  treeMetOut.1 <- reshape(treeMetOut, idvar = sampID, direction = 'wide',
                          timevar = 'PARAMETER', v.names = 'RESULT')
  names(treeMetOut.1) <- gsub("RESULT\\.", "", names(treeMetOut.1))
 
  return(treeMetOut.1)
}

