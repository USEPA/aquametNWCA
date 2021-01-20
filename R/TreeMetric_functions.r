#' @export
#' 
#' @title Calculate snag metrics
#' 
#' @description This function calculates all of the snag
#' metrics from tree data, and is called by function
#' calcTreeMetrics().
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
#' The following parameters are used in
#' calculating tree metrics: 'XXTHIN_SNAG', 'XTHIN_SNAG', 'THIN_SNAG', 
#' 'JR_SNAG', 'THICK_SNAG', 'XTHICK_SNAG', 'XXTHICK_SNAG'.
#' Additional parameters or variables are ignored.
#' 
#' @param nPlot A data frame with the 
#' number of plots sampled associated with
#' each sample, with \emph{sampID} variables and NPLOTS.
#' 
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default.
#' 
#' @details If any of the parameters are missing, they are assumed
#' to be zeros, and metric values associated with any metrics that
#' cannot be calculated due to missing parameters are set to 0.
#' 
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. The data frame contains the
#'   \emph{sampID} variables, PARAMETER, RESULT, where PARAMETER values include:
#'   TOTN_SNAGS, XN_SNAGS, TOTN_JR_SNAG, TOTN_THICK_SNAG, TOTN_THIN_SNAG,
#'   TOTN_XTHICK_SNAG, TOTN_XTHIN_SNAG, TOTN_XXTHIN_SNAG, XN_JR_SNAG,
#'   XN_THICK_SNAG, XN_THIN_SNAG, XN_XTHICK_SNAG, XN_XTHIN_SNAG, and
#'   XN_XXTHIN_SNAG. A list of metric descriptions is provided in the document
#'   named Tree_Metric_Descriptions.pdf included in the help directory for the
#'   package.
#'   
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(TreesEx)
#' 
#' nplots <- aggregate(x = list(NPLOTS = TreesEx$PLOT), by = TreesEx[c('UID')],
#' FUN = function(x){length(unique(x))})
#' 
#' snagEx <- calcSnagMets(TreesEx, nPlot=nplots, sampID='UID')
#'
#' head(snagEx)
#' unique(snagEx$PARAMETER)
calcSnagMets <- function(treeIn, nPlot, sampID='UID'){
  # Merge nPlot data frame and treeIn data frame by variables in sampID
  treeIn <- merge(treeIn, nPlot, by=sampID, all.y=TRUE)

  # Create SAMPID variable based on those listed in sampID argument
  for(i in 1:length(sampID)){
    if(i==1) treeIn$SAMPID <- treeIn[,sampID[i]]
    else treeIn$SAMPID <- paste(treeIn$SAMPID, treeIn[,sampID[i]], sep='.')
  }
  # Create data frame of unique samples in dataset, with sampID variables and SAMPID
  samples <- unique(subset(treeIn,select=c(sampID,'SAMPID')))
  # Create data frame of just unique SAMPID variables
  allUIDs <- data.frame(SAMPID=unique(treeIn$SAMPID),stringsAsFactors=F)

  # Subset data to just those parameters for snag counts
  snags <- subset(treeIn, PARAMETER %in% c('XXTHIN_SNAG','XTHIN_SNAG','THIN_SNAG','JR_SNAG',
                                           'THICK_SNAG','XTHICK_SNAG','XXTHICK_SNAG'))
  
  # If the resulting subset has at least 1 row, proceed with calculations
  if(nrow(snags)>0){
    # Sum RESULT values by PARAMETER
    snagsOut <- aggregate(x = list(TOTN = snags$RESULT), by = snags[c('SAMPID','PARAMETER','NPLOTS')],
                               FUN = function(x){sum(as.numeric(x))})
    
    # Calculate XN using results above and NPLOTS from nPlot input
    snagsOut$XN <- with(snagsOut, round((TOTN/NPLOTS), 2))
    # Melt data frame into long format
    snagsOut1 <- reshape(snagsOut, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                         varying = c('TOTN','XN'), times = c('TOTN','XN'),
                         timevar = 'METRIC', v.names = 'RESULT')
    # Combine METRIC and PARAMETER to create final METRIC name, then subset to just relevant variables
    snagsOut1$METRIC <- with(snagsOut1, paste(as.character(METRIC), PARAMETER, sep='_'))
    snagsOut1 <- subset(snagsOut1, select = c('SAMPID','METRIC','RESULT'))
    # Calculate summed TOTN and XN as overall snags metrics, then round XN to 2 digits
    totsnags <- aggregate(x = list(TOTN_SNAGS = snagsOut$TOTN, XN_SNAGS = snagsOut$XN),
                          by = snagsOut[c('SAMPID')], FUN = function(x){sum(as.numeric(x))})
    
    totsnags$XN_SNAGS <- with(totsnags, round(XN_SNAGS, 2))
    # Melt data frame into long format
    totsnags1 <- reshape(totsnags, idvar = 'SAMPID', direction = 'long',
                         varying = c('TOTN_SNAGS', 'XN_SNAGS'), times = c('TOTN_SNAGS','XN_SNAGS'),
                         timevar = 'METRIC', v.names = 'RESULT')
    # Combine totals with other metrics
    allSnagsOut <- rbind(totsnags1, snagsOut1)
    # Cast resulting data frame into wide format
    allSnagsOut1 <- reshape(allSnagsOut, idvar = 'SAMPID', direction = 'wide',
                            timevar = 'METRIC', v.names = 'RESULT')
    # Remove prefix from variable names
    names(allSnagsOut1) <- gsub("RESULT\\.", "", names(allSnagsOut1))
    # Create empty data frame to ensure that all metrics are represented in output data frame
    empty_snags <- data.frame(t(rep(NA,14)),stringsAsFactors=F)

    names(empty_snags) <- c("TOTN_SNAGS","XN_SNAGS","TOTN_JR_SNAG","TOTN_THICK_SNAG",
                            "TOTN_THIN_SNAG","TOTN_XTHICK_SNAG","TOTN_XTHIN_SNAG"
                            ,"TOTN_XXTHIN_SNAG","XN_JR_SNAG","XN_THICK_SNAG","XN_THIN_SNAG",
                            "XN_XTHICK_SNAG","XN_XTHIN_SNAG","XN_XXTHIN_SNAG")
    # Merge empty data frame with output df and remove the empty row - this ensures all metrics 
    # are included.
    allSnagsOut2 <- subset(merge(allSnagsOut1, empty_snags, all=TRUE),!is.na(SAMPID))

    # Merge with the all UIDs and fill in missing values with zeroes
    allSnagsOut3 <- merge(allUIDs, allSnagsOut2, by='SAMPID', all.x=T)
    # Melt merge data frame
    varNames <- names(allSnagsOut3)[!names(allSnagsOut3) %in% c('SAMPID')]
    allSnagsOut4 <- reshape(allSnagsOut3, idvar = 'SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')
    allSnagsOut4$METRIC <- with(allSnagsOut4, as.character(METRIC))
    # Replace missing values with 0 
    allSnagsOut4$RESULT <- with(allSnagsOut4, ifelse(is.na(RESULT), 0, RESULT))
    
  }else{ # Otherwise create empty data frame 
    # Create empty data frame as above
    empty_snags <- data.frame(t(rep(NA,14)), stringsAsFactors=F)

    names(empty_snags) <- c("TOTN_SNAGS","XN_SNAGS","TOTN_JR_SNAG","TOTN_THICK_SNAG",
                            "TOTN_THIN_SNAG","TOTN_XTHICK_SNAG","TOTN_XTHIN_SNAG"
                            ,"TOTN_XXTHIN_SNAG","XN_JR_SNAG","XN_THICK_SNAG","XN_THIN_SNAG",
                            "XN_XTHICK_SNAG","XN_XTHIN_SNAG","XN_XXTHIN_SNAG")
    # Merge list of SAMPIDs with empty data frame
    allSnagsOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)),stringsAsFactors=F), empty_snags, all=TRUE)
    # Remove the row without SAMPID
    allSnagsOut1 <- subset(allSnagsOut,!is.na(SAMPID))
    # Melt resulting data frame
    varNames <- names(allSnagsOut1)[!names(allSnagsOut3) %in% c('SAMPID')]
    allSnagsOut2 <- reshape(allSnagsOut1, idvar = 'SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')

    allSnagsOut2$METRIC <- with(allSnagsOut2, as.character(METRIC))
    # Set RESULT to 0 for all sites
    allSnagsOut2$RESULT <- 0

    allSnagsOut4 <- allSnagsOut2
  }
  print("Done with snag metrics")
  # Merge samples data frame to get back to the original sampID variables, rename METRIC to PARAMETER, 
  # then drop METRIC and SAMPID
  treeOut <- merge(samples, allSnagsOut4, by='SAMPID') 
  treeOut$PARAMETER <- treeOut$METRIC
  treeOut$METRIC <- NULL
  treeOut$SAMPID <- NULL

  return(treeOut)
}



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
#' nplots <- aggregate(x = list(NPLOTS = TreesEx$PLOT), by = TreesEx[c('UID')],
#' FUN = function(x){length(unique(x))})
#' 
#' tcntEx <- calcTreeCntMets(TreesEx, nPlot=nplots, sampID='UID')
#'
#' head(tcntEx)
#' unique(tcntEx$PARAMETER)
calcTreeCntMets <- function(treeIn, nPlot, sampID='UID'){
  # Merge nPlot with treeIn data frame by sampID variables
  treeIn <- merge(treeIn, nPlot, by=sampID, all.y=TRUE)
  
  # Create SAMPID variable based on those listed in sampID argument
  for(i in 1:length(sampID)){
    if(i==1) treeIn$SAMPID <- treeIn[,sampID[i]]
    else treeIn$SAMPID <- paste(treeIn$SAMPID, treeIn[,sampID[i]], sep='.')
  }
  # Create data frame of unique samples in dataset, with sampID variables and SAMPID
  samples <- unique(subset(treeIn, select=c(sampID, 'SAMPID')))
  # Create data frame of just unique SAMPID variables
  allUIDs <- data.frame(SAMPID=unique(treeIn$SAMPID), stringsAsFactors=F)
  # Subset treeIn to keep only tree count parameters
  tcnts <- subset(treeIn, PARAMETER %in% c('XXTHIN_TREE','XTHIN_TREE','THIN_TREE','JR_TREE',
                                           'THICK_TREE','XTHICK_TREE','XXTHICK_TREE'))
  
  # If resulting dataset has at least one row, perform calculations below
  if(nrow(tcnts)>0){
    # Sum counts by PARAMETER
    tcntsOut <- aggregate(x = list(TOTN = tcnts$RESULT), by = tcnts[c('SAMPID','PARAMETER','NPLOTS')],
                            FUN = function(x){sum(as.numeric(x))})
    
    # Use above results to calculate XN and round result to 2 digits
    tcntsOut$XN <- with(tcntsOut, round(TOTN/NPLOTS, 2))
    # Melt data frame into long format
    tcntsOut1 <- reshape(tcntsOut, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                         varying = c('TOTN','XN'), times = c('TOTN','XN'),
                         timevar = 'METRIC', v.names = 'RESULT')
    # Update value of METRIC by combining METRIC and PARAMETER
    tcntsOut1$METRIC <- with(tcntsOut1, paste(as.character(METRIC), PARAMETER, sep='_'))
    # Subset to relevant variables
    tcntsOut1 <- subset(tcntsOut1, select = c('SAMPID','METRIC','RESULT'))
    # Calculate totals across all sizes for both TOTN and XN, then round to 2 digits
    tottrees <- aggregate(x = list(TOTN_TREES = tcntsOut$TOTN, XN_TREES = tcntsOut$XN),
                          by = tcntsOut[c('SAMPID')], FUN = sum)
    tottrees$XN_TREES <- with(tottrees, round(XN_TREES, 2))
    # Melt data frame in order to merge with parameter-based metrics above
    tottrees1 <- reshape(tottrees, idvar = 'SAMPID', direction = 'long',
                         varying = c('TOTN_TREES','XN_TREES'), times = c('TOTN_TREES','XN_TREES'),
                         timevar = 'METRIC', v.names = 'RESULT')
    # Combine totals with other metrics
    allTreesOut <- rbind(tottrees1, tcntsOut1)
    # Cast data frame wide and remove prefix in resulting data frame
    allTreesOut1 <- reshape(allTreesOut, idvar = 'SAMPID', direction = 'wide',
                            timevar = 'METRIC', v.names = 'RESULT')
    names(allTreesOut1) <- gsub("RESULT\\.", "", names(allTreesOut1))
    # Create empty data frame with names of all metrics that should be calculated
    empty_trees <- data.frame(t(rep(NA,16)), stringsAsFactors=F)
    names(empty_trees) <- c("TOTN_TREES","XN_TREES","TOTN_JR_TREE","TOTN_THICK_TREE","TOTN_THIN_TREE","TOTN_XTHICK_TREE","TOTN_XTHIN_TREE"
                            ,"TOTN_XXTHICK_TREE","TOTN_XXTHIN_TREE","XN_JR_TREE","XN_THICK_TREE","XN_THIN_TREE","XN_XTHICK_TREE","XN_XTHIN_TREE"
                            ,"XN_XXTHICK_TREE","XN_XXTHIN_TREE")
    
    # Merge empty data frame with metrics, then drop row with missing SAMPID
    allTreesOut2 <- subset(merge(allTreesOut1, empty_trees, all=TRUE), !is.na(SAMPID))
    # Merge resulting data frame with list of allUIDs to make sure all samples are represented
    allTreesOut3 <- merge(allUIDs, allTreesOut2, by='SAMPID', all.x=T)
    
    # Melt resulting data frame, then set missing RESULT values to 0
    varNames <- names(allTreesOut3)[!names(allTreesOut3) %in% c('SAMPID')]
    allTreesOut4 <- reshape(allTreesOut3, idvar = 'SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')
    allTreesOut4$METRIC <- with(allTreesOut4, as.character(METRIC))
    allTreesOut4$RESULT <- with(allTreesOut4, ifelse(is.na(RESULT), 0, RESULT))
    # Cast data frame wide and remove prefix added by reshape() function
    allTreesOut4.wide <- reshape(allTreesOut4, idvar='SAMPID', direction = 'wide',
                                 timevar = 'METRIC', v.names = 'RESULT')
    names(allTreesOut4.wide) <- gsub('RESULT\\.', '', names(allTreesOut4.wide))
    # Calculate additional TOTN and XN metrics based on combinations of previously calculated metrics
    allTreesOut4.wide$TOTN_SMALL <- with(allTreesOut4.wide, TOTN_XXTHIN_TREE+TOTN_XTHIN_TREE)
    allTreesOut4.wide$TOTN_MID <- with(allTreesOut4.wide, TOTN_THIN_TREE+TOTN_JR_TREE)
    allTreesOut4.wide$TOTN_LARGE <- with(allTreesOut4.wide, TOTN_THICK_TREE+TOTN_XTHICK_TREE+TOTN_XXTHICK_TREE)
    allTreesOut4.wide$XN_SMALL <- with(allTreesOut4.wide, XN_XXTHIN_TREE+XN_XTHIN_TREE)
    allTreesOut4.wide$XN_MID <- with(allTreesOut4.wide, XN_THIN_TREE+XN_JR_TREE)
    allTreesOut4.wide$XN_LARGE <- with(allTreesOut4.wide, XN_THICK_TREE+XN_XTHICK_TREE+XN_XXTHICK_TREE)
    # Melt metrics into long format
    varNames <- names(allTreesOut4.wide)[!names(allTreesOut4.wide) %in% c('SAMPID')]
    allTreesOut5 <- reshape(allTreesOut4.wide, idvar='SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')
    
  }else{ # If no tree counts in input dataset, create empty data frame, then add variable names
    empty_trees <- data.frame(t(rep(NA,16)), stringsAsFactors=F)
    names(empty_trees) <- c("TOTN_TREES","XN_TREES","TOTN_JR_TREE","TOTN_THICK_TREE",
                            "TOTN_THIN_TREE","TOTN_XTHICK_TREE","TOTN_XTHIN_TREE",
                            "TOTN_XXTHICK_TREE","TOTN_XXTHIN_TREE","XN_JR_TREE","XN_THICK_TREE",
                            "XN_THIN_TREE","XN_XTHICK_TREE","XN_XTHIN_TREE",
                            "XN_XXTHICK_TREE","XN_XXTHIN_TREE","TOTN_SMALL","TOTN_MID","TOTN_LARGE",
                            "XN_SMALL","XN_MID","XN_LARGE")
    # Merge empty data frame with list of samples
    allTreesOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)), stringsAsFactors=F), 
                         empty_trees, all=TRUE)
    # Subset to remove row with missing SAMPID
    allTreesOut1 <- subset(allTreesOut, !is.na(SAMPID))
    # Melt output data frame and set all RESULT values to 0
    varNames <- names(allTreesOut1)[!names(allTreesOut1) %in% c('SAMPID')]
    allTreesOut2 <- reshape(allTreesOut2, idvar = 'SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')

    allTreesOut2$METRIC <- with(allTreesOut2, as.character(METRIC))
    allTreesOut2$RESULT <- 0

    allTreesOut5 <- allTreesOut2
  }
  
    print("Done with tree count metrics")
  # Merge with samples data frame to get back sampID variables, rename METRIC to PARAMETER, 
  # then drop METRIC and SAMPID variables.
  treeOut <- merge(samples, allTreesOut5, by='SAMPID') 
  treeOut$PARAMETER <- treeOut$METRIC
  treeOut$METRIC <- NULL
  treeOut$SAMPID <- NULL

  return(treeOut)
}



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
#'   descriptions is provided in the document named Tree_Metric_Descriptions.pdf
#'   included in the help directory for the package.
#' 
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(TreesEx)
#' nplots <- aggregate(x = list(NPLOTS = TreesEx$PLOT), by = TreesEx[c('UID')],
#' FUN = function(x){length(unique(x))})
#' 
#' tcvrEx <- calcTreeCoverMets(TreesEx, nPlot=nplots, sampID='UID')
#'
#' head(tcvrEx)
#' unique(tcvrEx$PARAMETER)
calcTreeCoverMets <- function(treeIn, nPlot, sampID='UID'){
  # Merge nPlot with treeIn data frame by sampID variables
  treeIn <- merge(treeIn, nPlot, by=sampID, all.y=TRUE)
  # Create SAMPID variable based on those listed in sampID argument
  for(i in 1:length(sampID)){
    if(i==1) treeIn$SAMPID <- treeIn[,sampID[i]]
    else treeIn$SAMPID <- paste(treeIn$SAMPID,treeIn[,sampID[i]],sep='.')
  }
  # Create data frame of unique samples in dataset, with sampID variables and SAMPID
  samples <- unique(subset(treeIn,select=c(sampID, 'SAMPID')))
  # Create data frame of just unique SAMPID variables
  allUIDs <- data.frame(SAMPID=unique(treeIn$SAMPID), stringsAsFactors=F)
  
  ##### TREE SPECIES METRICS ##############################
  # Subset treeIn to keep only tree height parameters
  tcvr <- subset(treeIn,PARAMETER %in% c('VSMALL_TREE','SMALL_TREE','LMED_TREE','HMED_TREE',
                                         'TALL_TREE','VTALL_TREE'))
  # If resulting dataset has at least one row, perform calculations below
  if(nrow(tcvr)>0){
    # Subset input data to tree species data and cast data frame wide and 
    # remove prefix added by reshape()
    tspp <- reshape(subset(treeIn, PARAMETER=='TREE_SPECIES', 
                           select = c('SAMPID','PAGE','LINE','PLOT','PARAMETER','RESULT')), 
                    idvar = c('SAMPID','PAGE','LINE','PLOT'), direction = 'wide',
                    timevar = 'PARAMETER', v.names = 'RESULT')
    names(tspp) <- gsub('RESULT\\.', '', names(tspp)) 
  
    # Merge tree species and tree cover data
    tcvr1 <- merge(tspp, tcvr,by=c('SAMPID','PLOT','PAGE','LINE'), all=TRUE)
    # Create alternate parameter values by combining original metric values
    tcvr1$PARAM_ALT <- NA
    tcvr1$PARAM_ALT[tcvr1$PARAMETER %in% c('VSMALL_TREE','SMALL_TREE')] <- 'TREE_GROUND'
    tcvr1$PARAM_ALT[tcvr1$PARAMETER %in% c('LMED_TREE','HMED_TREE')] <- 'TREE_MID'
    tcvr1$PARAM_ALT[tcvr1$PARAMETER %in% c('TALL_TREE','VTALL_TREE')] <- 'TREE_UPPER'
    # Calculate number of species 
    ntrspp <- aggregate(x = list(N_TREESPP = tcvr1$TREE_SPECIES), 
                        by = tcvr1[c('SAMPID')], FUN = function(x){length(unique(x))})
    # Merge count with cover and species data from above
    totspp <- merge(tcvr1, ntrspp, by = c('SAMPID'))
    # Melt data frame into long format for metric calculation
    totspp1 <- reshape(ntrspp, idvar = 'SAMPID', direction = 'long',
                       varying = 'N_TREESPP', times = 'N_TREESPP',
                       timevar = 'METRIC', v.names = 'RESULT')
    
    # Subset to values above 0 and not missing
    totspp.pos <- subset(totspp, RESULT!='0' & !is.na(RESULT), 
                         select = c('SAMPID','PARAMETER','PARAM_ALT','TREE_SPECIES'))
    # Calculate number of unique species by PARAMETER value (N metrics)
    tspp1a <- aggregate(x = list(N = totspp.pos$TREE_SPECIES), 
                        by = totspp.pos[c('SAMPID','PARAMETER')],
                        FUN = function(x){length(unique(x))})
    # Drop for any cases where PARAMETER is missing
    tspp1a <- subset(tspp1a, !is.na(PARAMETER))
    # Calculate number of unique species by PARAM_ALT value
    tspp1b <- aggregate(x = list(N = totspp.pos$TREE_SPECIES), 
                        by = totspp.pos[c('SAMPID','PARAM_ALT')],
                        FUN = function(x){length(unique(x))})
    # Drop for any cases where PARAM_ALT missing
    tspp1b <- subset(tspp1b, !is.na(PARAM_ALT))
    # Merge above data frame with data frame with total number of species
    tspp1c <- merge(tspp1b, ntrspp, by='SAMPID')
    tspp1c <- subset(tspp1c, !is.na(PARAM_ALT))
    # Calculate PCTN (% species) and round to 2 digits
    tspp1c$PCTN <- with(tspp1c, round((N/N_TREESPP)*100, 2))
    # Drop total number of tree species metric
    tspp1c$N_TREESPP <- NULL
    
    # Melt first and third data frames created above, updating metric name by 
    # combining with PARAMETER value
    tspp2a <- reshape(tspp1a, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                      varying = 'N', times = 'N',
                      timevar = 'METRIC', v.names = 'RESULT')
    tspp2a$METRIC <- with(tspp2a, paste(METRIC, PARAMETER, sep='_'))

    tspp2c <- reshape(tspp1c, idvar = c('SAMPID','PARAM_ALT'), direction = 'long',
                       varying = c('PCTN','N'), times = c('PCTN','N'),
                       timevar = 'METRIC', v.names = 'RESULT')
    tspp2c$METRIC <- with(tspp2c, paste(METRIC, PARAM_ALT, sep='_'))
    # Combine total number of species data frame and two sets of metrics above
    tsppOut <- rbind(totspp1,subset(tspp2a,select=-PARAMETER),subset(tspp2c,select=-PARAM_ALT))
    # Cast resulting data frame wide and remove prefix added by reshape
    tsppOut1 <- reshape(tsppOut, idvar = c('SAMPID'), direction = 'wide',
                        timevar = 'METRIC', v.names = 'RESULT')
    names(tsppOut1) <- gsub("RESULT\\.", "", names(tsppOut1))
    
    # Create empty data frame with full set of expected metrics
    empty_tspp <- data.frame(t(rep(NA,13)))
    names(empty_tspp) <- c("N_TREESPP","N_TALL_TREE","N_HMED_TREE","N_LMED_TREE","N_SMALL_TREE",
                           "N_VSMALL_TREE","N_VTALL_TREE","N_TREE_UPPER",
                           "N_TREE_MID","N_TREE_GROUND","PCTN_TREE_UPPER","PCTN_TREE_MID",
                           "PCTN_TREE_GROUND")
    # Merge output data frame with empty data frame and remove row missing SAMPID
    tsppOut2 <- subset(merge(tsppOut1, empty_tspp, all=TRUE), !is.na(SAMPID))
    # Merge full set of UIDs with resulting data frame to ensure all samples included
    tsppOut3 <- merge(allUIDs, tsppOut2, by='SAMPID', all.x=T)
    # Melt resulting data frame, then replace missing RESULT values with 0
    varNames <- names(tsppOut3)[!names(tsppOut3) %in% c('SAMPID')]
    tsppOut4 <- reshape(tsppOut3, idvar = 'SAMPID', direction = 'long',
                        varying = varNames, times = varNames,
                        timevar = 'METRIC', v.names = 'RESULT')
    tsppOut4$METRIC <- with(tsppOut4, as.character(METRIC))
    tsppOut4$RESULT <- with(tsppOut4, ifelse(is.na(RESULT), 0, RESULT))
    
  }else{
    # If no values for these parameters in input data frame, create empty data frame
    empty_tspp <- data.frame(t(rep(NA,13)), stringsAsFactors=F)
    names(empty_tspp) <- c("N_TREESPP","N_TALL_TREE","N_HMED_TREE","N_LMED_TREE","N_SMALL_TREE",
                           "N_VSMALL_TREE","N_VTALL_TREE","N_TREE_UPPER"
                           ,"N_TREE_MID","N_TREE_GROUND","PCTN_TREE_UPPER","PCTN_TREE_MID",
                           "PCTN_TREE_GROUND")
    # Merge data frame with list of all UIDs in input data frame, dropping row with missing SAMPID
    tsppOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)), stringsAsFactors=F), empty_tspp, all=TRUE)
    tsppOut1 <- subset(tsppOut, !is.na(SAMPID))
    # Melt data frame and set all RESULT values to 0
    varNames <- names(tsppOut1)[!names(tsppOut1) %in% c('SAMPID')]
    tsppOut2 <- reshape(tsppOut1, idvar = 'SAMPID', direction = 'long',
                        varying = varNames, times = varNames,
                        timevar = 'METRIC', v.names = 'RESULT')
    tsppOut2$METRIC <- with(tsppOut2, as.character(METRIC))
    tsppOut2$RESULT <- 0
    tsppOut4 <- tsppOut2
  }
  print("Done with tree species metrics")

  ## TREE COVER METRICS ###########
  ## Sum by species within plot
  if(nrow(tcvr)>0){
    # Sum cover by PARAMETER
    tcvr2a <- aggregate(x = list(COV = tcvr1$RESULT), 
                        by = tcvr1[c('SAMPID','PLOT','NPLOTS','PARAMETER','TREE_SPECIES')],
                        FUN = function(x){sum(as.numeric(x))})
    # Cap values at 100 percent and keep only non-zero values
    tcvr2a$COV <- ifelse(tcvr2a$COV>100, 100, tcvr2a$COV)
    tcvr2a <- subset(tcvr2a, COV!=0)
    # Sum cover by PARAM_ALT
    tcvr2b <- aggregate(x = list(COV = tcvr1$RESULT), 
                        by = tcvr1[c('SAMPID','PLOT','NPLOTS','PARAM_ALT','TREE_SPECIES')],
                        FUN = function(x){sum(as.numeric(x))})
    # Cap values at 100 percent and keep only non-zero values
    tcvr2b$COV <- ifelse(tcvr2b$COV>100, 100, tcvr2b$COV)
    tcvr2b <- subset(tcvr2b, COV!=0)
    # Take unique values of NPLOTS
    tcvr3a.uniq <- aggregate(x = list(uniqN = tcvr2a$NPLOTS),
                             by = tcvr2a[c('SAMPID','PARAMETER')], FUN = unique)
    # Calculate number of unique plots by parameter
    tcvr3a.length <- aggregate(x = list(uniqPlot = tcvr2a$PLOT),
                               by = tcvr2a[c('SAMPID','PARAMETER')],
                               FUN = function(x){length(unique(x))})
    # Sum Cover by PARAMETER
    tcvr3a.sum <- aggregate(x = list(sumcov = tcvr2a$COV),
                            by = tcvr2a[c('SAMPID','PARAMETER')],
                            FUN = function(x){sum(as.numeric(x))})
    # Merge length with unique plots and with cover data frame
    tcvr3a <- merge(tcvr3a.length, tcvr3a.uniq, by = c('SAMPID','PARAMETER'))
    tcvr3a <- merge(tcvr3a, tcvr3a.sum, by = c('SAMPID','PARAMETER'))
    # Now complete calculations of FREQ, XCOV, and IMP for each sample and parameter
    tcvr3a$FREQ <- with(tcvr3a, round((uniqPlot/uniqN)*100, 2))
    tcvr3a$XCOV <- with(tcvr3a, round(sumcov/uniqN, 2))
    tcvr3a$IMP <- with(tcvr3a, round((FREQ + XCOV)/2, 2))
    # Subset to select relevant variables
    tcvr3a <- subset(tcvr3a, select = c('SAMPID','PARAMETER','FREQ','XCOV','IMP'))
    # Again get unique NPLOTS values by sample
    tcvr3b.uniq <- aggregate(x = list(uniqN = tcvr2b$NPLOTS),
                             by = tcvr2b[c('SAMPID','PARAM_ALT')], FUN = unique)
    # Calculate number of unique plots by PARAM_ALT
    tcvr3b.length <- aggregate(x = list(uniqPlot = tcvr2b$PLOT),
                               by = tcvr2b[c('SAMPID','PARAM_ALT')],
                               FUN = function(x){length(unique(x))})
    # Sum Cover by PARAM_ALT
    tcvr3b.sum <- aggregate(x = list(sumcov = tcvr2b$COV),
                            by = tcvr2b[c('SAMPID','PARAM_ALT')],
                            FUN = function(x){sum(as.numeric(x))})
    # Merge length with unique plots and with cover data frame
    tcvr3b <- merge(tcvr3b.length, tcvr3b.uniq, by = c('SAMPID','PARAM_ALT'))
    tcvr3b <- merge(tcvr3b, tcvr3b.sum, by = c('SAMPID','PARAM_ALT'))
    # Now complete calculations of FREQ, XCOV, and IMP for each sample and PARAM_ALT
    tcvr3b$FREQ <- with(tcvr3b, round((uniqPlot/uniqN)*100, 2))
    tcvr3b$XCOV <- with(tcvr3b, round(sumcov/uniqN, 2))
    tcvr3b$IMP <- with(tcvr3b, round((FREQ + XCOV)/2, 2))
    # Subset to keep relevant variables
    tcvr3b <- subset(tcvr3b, select = c('SAMPID','PARAM_ALT','FREQ','XCOV','IMP'))
    # Melt data frame by PARAMETER and update METRIC value by combining METRIC and PARAMETER values
    varNames.a <- names(tcvr3a)[!names(tcvr3a) %in% c('SAMPID','PARAMETER')]
    tcvr4a <- reshape(tcvr3a, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                      varying = varNames.a, times = varNames.a,
                      timevar = 'METRIC', v.names = 'RESULT')
    tcvr4a$METRIC = with(tcvr4a, paste(METRIC, PARAMETER, sep='_'))
    # Melt data frame by PARAM_ALT and update METRIC value by combining METRIC and PARAMETER values
    varNames.b <- names(tcvr3b)[!names(tcvr3b) %in% c('SAMPID','PARAM_ALT')]
    tcvr4b <- reshape(tcvr3b, idvar = c('SAMPID','PARAM_ALT'), direction = 'long',
                      varying = varNames.b, times = varNames.b,
                      timevar = 'METRIC', v.names = 'RESULT')
    tcvr4b$METRIC = with(tcvr4b, paste(METRIC, PARAM_ALT, sep='_'))
    # Combine resulting data frames after dropping PARAMETER and PARAM_ALT variables
    tcvrOut <- rbind(subset(tcvr4a,select=-PARAMETER),subset(tcvr4b,select=-PARAM_ALT))
    # Cast resulting data frame and drop prefix added by reshape()
    tcvrOut1 <- reshape(tcvrOut, idvar = 'SAMPID', direction = 'wide',
                        timevar = 'METRIC', v.names = 'RESULT')
    names(tcvrOut1) <- gsub("RESULT\\.", '', names(tcvrOut1))
    # Create empty data frame
    empty_tcvr <- data.frame(t(rep(NA,27)), stringsAsFactors=F)
    names(empty_tcvr) <- c("FREQ_TALL_TREE","FREQ_HMED_TREE","FREQ_LMED_TREE","FREQ_SMALL_TREE",
                           "FREQ_VSMALL_TREE","FREQ_VTALL_TREE","XCOV_TALL_TREE",
                           "XCOV_HMED_TREE","XCOV_LMED_TREE","XCOV_SMALL_TREE","XCOV_VSMALL_TREE",
                           "XCOV_VTALL_TREE","IMP_TALL_TREE","IMP_HMED_TREE",
                           "IMP_LMED_TREE","IMP_SMALL_TREE","IMP_VSMALL_TREE","IMP_VTALL_TREE",
                           "FREQ_TREE_UPPER","FREQ_TREE_MID","FREQ_TREE_GROUND",
                           "XCOV_TREE_UPPER","XCOV_TREE_MID","XCOV_TREE_GROUND",
                           "IMP_TREE_UPPER","IMP_TREE_MID","IMP_TREE_GROUND")
    # Merge with output data frame and drop row with missing SAMPID
    tcvrOut2 <- subset(merge(tcvrOut1, empty_tcvr, all=TRUE), !is.na(SAMPID))
    # Merge with allUIDs to make sure all samples represented
    tcvrOut3 <- merge(allUIDs, tcvrOut2, by='SAMPID', all.x=T)
    # Melt data frame and set missing RESULT values to 0
    varNames <- names(tcvrOut3)[!names(tcvrOut3) %in% c('SAMPID')]
    tcvrOut4 <- reshape(tcvrOut3, idvar = 'SAMPID', direction = 'long',
                        varying = varNames, times = varNames,
                        timevar = 'METRIC', v.names = 'RESULT')
    
    tcvrOut4$METRIC <- with(tcvrOut4, as.character(METRIC))
    tcvrOut4$RESULT <- with(tcvrOut4, ifelse(is.na(RESULT), 0, RESULT))

  }else{ # If no data for this subset, create empty data frame
    empty_tcvr <- data.frame(t(rep(NA,27)),stringsAsFactors=F)
    names(empty_tcvr) <- c("FREQ_TALL_TREE","FREQ_HMED_TREE","FREQ_LMED_TREE","FREQ_SMALL_TREE",
                           "FREQ_VSMALL_TREE","FREQ_VTALL_TREE","XCOV_TALL_TREE",
                           "XCOV_HMED_TREE","XCOV_LMED_TREE","XCOV_SMALL_TREE","XCOV_VSMALL_TREE",
                           "XCOV_VTALL_TREE","IMP_TALL_TREE","IMP_HMED_TREE",
                           "IMP_LMED_TREE","IMP_SMALL_TREE","IMP_VSMALL_TREE","IMP_VTALL_TREE",
                           "FREQ_TREE_UPPER","FREQ_TREE_MID","FREQ_TREE_GROUND",
                           "XCOV_TREE_UPPER","XCOV_TREE_MID","XCOV_TREE_GROUND",
                           "IMP_TREE_UPPER","IMP_TREE_MID","IMP_TREE_GROUND")
    # Merge empty data frame with SAMPID list, drop missing SAMPID
    tcvrOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)), stringsAsFactors=F), 
                     empty_tcvr, all=TRUE)
    tcvrOut1 <- subset(tcvrOut, !is.na(SAMPID))
    # Melt data frame and set RESULT to 0
    varNames <- names(tcvrOut1)[!names(tcvrOut1) %in% c('SAMPID')]
    tcvrOut2 <- reshape(tcvrOut1, idvar = 'SAMPID', direction = 'long',
                        varying = varNames, times = varNames,
                        timevar = 'METRIC', v.names = 'RESULT')
    
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
  treeOut.1 <- merge(samples, treeOut, by='SAMPID') 
  treeOut.1 <- subset(treeOut.1, select = -SAMPID)
  
  return(treeOut.1)
}
