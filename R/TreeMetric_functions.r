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

  treeIn <- merge(treeIn,nPlot,by=sampID, all.y=TRUE)
  # Create vector of all samples in dataset
  for(i in 1:length(sampID)){
    if(i==1) treeIn$SAMPID <- treeIn[,sampID[i]]
    else treeIn$SAMPID <- paste(treeIn$SAMPID,treeIn[,sampID[i]],sep='.')
  }
  samples <- unique(subset(treeIn,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(treeIn$SAMPID),stringsAsFactors=F)

  snags <- subset(treeIn,PARAMETER %in% c('XXTHIN_SNAG','XTHIN_SNAG','THIN_SNAG','JR_SNAG','THICK_SNAG','XTHICK_SNAG','XXTHICK_SNAG'))
  if(nrow(snags)>0){
    snagsOut <- aggregate(x = list(TOTN = snags$RESULT), by = snags[c('SAMPID','PARAMETER','NPLOTS')],
                               FUN = function(x){sum(as.numeric(x))})
    
   
    snagsOut$XN <- with(snagsOut, round((TOTN/NPLOTS), 2))
    
    snagsOut1 <- reshape(snagsOut, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                         varying = c('TOTN','XN'), times = c('TOTN','XN'),
                         timevar = 'METRIC', v.names = 'RESULT')
    
    snagsOut1$METRIC <- with(snagsOut1, paste(as.character(METRIC), PARAMETER, sep='_'))
    snagsOut1 <- subset(snagsOut1, select = c('SAMPID','METRIC','RESULT'))

    totsnags <- aggregate(x = list(TOTN_SNAGS = snagsOut$TOTN, XN_SNAGS = snagsOut$XN),
                          by = snagsOut[c('SAMPID')], FUN = function(x){sum(as.numeric(x))})
    totsnags$XN_SNAGS <- with(totsnags, round(XN_SNAGS, 2))
    totsnags1 <- reshape(totsnags, idvar = 'SAMPID', direction = 'long',
                         varying = c('TOTN_SNAGS', 'XN_SNAGS'), times = c('TOTN_SNAGS','XN_SNAGS'),
                         timevar = 'METRIC', v.names = 'RESULT')

    allSnagsOut <- rbind(totsnags1, snagsOut1)
    
    allSnagsOut1 <- reshape(allSnagsOut, idvar = 'SAMPID', direction = 'wide',
                            timevar = 'METRIC', v.names = 'RESULT')
    names(allSnagsOut1) <- gsub("RESULT\\.", "", names(allSnagsOut1))

    empty_snags <- data.frame(t(rep(NA,14)),stringsAsFactors=F)

    names(empty_snags) <- c("TOTN_SNAGS","XN_SNAGS","TOTN_JR_SNAG","TOTN_THICK_SNAG","TOTN_THIN_SNAG","TOTN_XTHICK_SNAG","TOTN_XTHIN_SNAG"
                            ,"TOTN_XXTHIN_SNAG","XN_JR_SNAG","XN_THICK_SNAG","XN_THIN_SNAG","XN_XTHICK_SNAG","XN_XTHIN_SNAG","XN_XXTHIN_SNAG")

    allSnagsOut2 <- subset(merge(allSnagsOut1, empty_snags, all=TRUE),!is.na(SAMPID))

    # Merge with the all UIDs and fill in missing values with zeroes
    allSnagsOut3 <- merge(allUIDs,allSnagsOut2,by='SAMPID',all.x=T)

    varNames <- names(allSnagsOut3)[!names(allSnagsOut3) %in% c('SAMPID')]
    allSnagsOut4 <- reshape(allSnagsOut3, idvar = 'SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')
    allSnagsOut4$METRIC <- with(allSnagsOut4, as.character(METRIC))
    allSnagsOut4$RESULT <- with(allSnagsOut4, ifelse(is.na(RESULT),0,RESULT))
  }else{
    empty_snags <- data.frame(t(rep(NA,14)),stringsAsFactors=F)

    names(empty_snags) <- c("TOTN_SNAGS","XN_SNAGS","TOTN_JR_SNAG","TOTN_THICK_SNAG","TOTN_THIN_SNAG","TOTN_XTHICK_SNAG","TOTN_XTHIN_SNAG"
                            ,"TOTN_XXTHIN_SNAG","XN_JR_SNAG","XN_THICK_SNAG","XN_THIN_SNAG","XN_XTHICK_SNAG","XN_XTHIN_SNAG","XN_XXTHIN_SNAG")

    allSnagsOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)),stringsAsFactors=F), empty_snags, all=TRUE)

    allSnagsOut1 <- subset(allSnagsOut,!is.na(SAMPID))
    
    varNames <- names(allSnagsOut1)[!names(allSnagsOut3) %in% c('SAMPID')]
    allSnagsOut2 <- reshape(allSnagsOut1, idvar = 'SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')

    allSnagsOut2$METRIC <- with(allSnagsOut2, as.character(METRIC))
    allSnagsOut2$RESULT <- 0

    allSnagsOut4 <- allSnagsOut2
  }
  print("Done with snag metrics")

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
  treeIn <- merge(treeIn, nPlot, by=sampID, all.y=TRUE)
  
  for(i in 1:length(sampID)){
    if(i==1) treeIn$SAMPID <- treeIn[,sampID[i]]
    else treeIn$SAMPID <- paste(treeIn$SAMPID,treeIn[,sampID[i]],sep='.')
  }
  samples <- unique(subset(treeIn,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(treeIn$SAMPID),stringsAsFactors=F)
  
  tcnts <- subset(treeIn,PARAMETER %in% c('XXTHIN_TREE','XTHIN_TREE','THIN_TREE','JR_TREE','THICK_TREE','XTHICK_TREE','XXTHICK_TREE'))
  
  
  if(nrow(tcnts)>0){
    tcntsOut <- aggregate(x = list(TOTN = tcnts$RESULT), by = tcnts[c('SAMPID','PARAMETER','NPLOTS')],
                            FUN = function(x){sum(as.numeric(x))})
    

    tcntsOut$XN <- with(tcntsOut, round(TOTN/NPLOTS, 2))

    tcntsOut1 <- reshape(tcntsOut, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                         varying = c('TOTN','XN'), times = c('TOTN','XN'),
                         timevar = 'METRIC', v.names = 'RESULT')
    
    tcntsOut1$METRIC <- with(tcntsOut1, paste(as.character(METRIC), PARAMETER, sep='_'))
    tcntsOut1 <- subset(tcntsOut1, select = c('SAMPID','METRIC','RESULT'))

    tottrees <- aggregate(x = list(TOTN_TREES = tcntsOut$TOTN, XN_TREES = tcntsOut$XN),
                          by = tcntsOut[c('SAMPID')], FUN = sum)
    tottrees$XN_TREES <- with(tottrees, round(XN_TREES, 2))

    tottrees1 <- reshape(tottrees, idvar = 'SAMPID', direction = 'long',
                         varying = c('TOTN_TREES','XN_TREES'), times = c('TOTN_TREES','XN_TREES'),
                         timevar = 'METRIC', v.names = 'RESULT')

    allTreesOut <- rbind(tottrees1, tcntsOut1)
    
    allTreesOut1 <- reshape(allTreesOut, idvar = 'SAMPID', direction = 'wide',
                            timevar = 'METRIC', v.names = 'RESULT')
    names(allTreesOut1) <- gsub("RESULT\\.", "", names(allTreesOut1))

    empty_trees <- data.frame(t(rep(NA,16)),stringsAsFactors=F)
    names(empty_trees) <- c("TOTN_TREES","XN_TREES","TOTN_JR_TREE","TOTN_THICK_TREE","TOTN_THIN_TREE","TOTN_XTHICK_TREE","TOTN_XTHIN_TREE"
                            ,"TOTN_XXTHICK_TREE","TOTN_XXTHIN_TREE","XN_JR_TREE","XN_THICK_TREE","XN_THIN_TREE","XN_XTHICK_TREE","XN_XTHIN_TREE"
                            ,"XN_XXTHICK_TREE","XN_XXTHIN_TREE")

    allTreesOut2 <- subset(merge(allTreesOut1, empty_trees, all=TRUE),!is.na(SAMPID))

    allTreesOut3 <- merge(allUIDs, allTreesOut2, by='SAMPID', all.x=T)

    varNames <- names(allTreesOut3)[!names(allTreesOut3) %in% c('SAMPID')]
    allTreesOut4 <- reshape(allTreesOut3, idvar = 'SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')
    allTreesOut4$METRIC <- with(allTreesOut4, as.character(METRIC))
    allTreesOut4$RESULT <- with(allTreesOut4, ifelse(is.na(RESULT),0,RESULT))

    allTreesOut4.wide <- reshape(allTreesOut4, idvar='SAMPID', direction = 'wide',
                                 timevar = 'METRIC', v.names = 'RESULT')
    names(allTreesOut4.wide) <- gsub('RESULT\\.', '', names(allTreesOut4.wide))
    
    allTreesOut4.wide$TOTN_SMALL <- with(allTreesOut4.wide, TOTN_XXTHIN_TREE+TOTN_XTHIN_TREE)
    allTreesOut4.wide$TOTN_MID <- with(allTreesOut4.wide, TOTN_THIN_TREE+TOTN_JR_TREE)
    allTreesOut4.wide$TOTN_LARGE <- with(allTreesOut4.wide, TOTN_THICK_TREE+TOTN_XTHICK_TREE+TOTN_XXTHICK_TREE)
    allTreesOut4.wide$XN_SMALL <- with(allTreesOut4.wide, XN_XXTHIN_TREE+XN_XTHIN_TREE)
    allTreesOut4.wide$XN_MID <- with(allTreesOut4.wide, XN_THIN_TREE+XN_JR_TREE)
    allTreesOut4.wide$XN_LARGE <- with(allTreesOut4.wide, XN_THICK_TREE+XN_XTHICK_TREE+XN_XXTHICK_TREE)
    
    varNames <- names(allTreesOut4.wide)[!names(allTreesOut4.wide) %in% c('SAMPID')]
    allTreesOut5 <- reshape(allTreesOut4.wide, idvar='SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')
    
  }else{
    empty_trees <- data.frame(t(rep(NA,16)),stringsAsFactors=F)
    names(empty_trees) <- c("TOTN_TREES","XN_TREES","TOTN_JR_TREE","TOTN_THICK_TREE","TOTN_THIN_TREE","TOTN_XTHICK_TREE","TOTN_XTHIN_TREE"
                            ,"TOTN_XXTHICK_TREE","TOTN_XXTHIN_TREE","XN_JR_TREE","XN_THICK_TREE","XN_THIN_TREE","XN_XTHICK_TREE","XN_XTHIN_TREE"
                            ,"XN_XXTHICK_TREE","XN_XXTHIN_TREE","TOTN_SMALL","TOTN_MID","TOTN_LARGE",
                            "XN_SMALL","XN_MID","XN_LARGE")

    allTreesOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)),stringsAsFactors=F), empty_trees, all=TRUE)

    allTreesOut1 <- subset(allTreesOut,!is.na(SAMPID))

    varNames <- names(allTreesOut1)[!names(allTreesOut1) %in% c('SAMPID')]
    allTreesOut2 <- reshape(allTreesOut2, idvar = 'SAMPID', direction = 'long',
                            varying = varNames, times = varNames,
                            timevar = 'METRIC', v.names = 'RESULT')

    allTreesOut2$METRIC <- with(allTreesOut2, as.character(METRIC))
    allTreesOut2$RESULT <- 0

    allTreesOut5 <- allTreesOut2
  }
  
    print("Done with tree count metrics")
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
  treeIn <- merge(treeIn,nPlot,by=sampID, all.y=TRUE)
  
  for(i in 1:length(sampID)){
    if(i==1) treeIn$SAMPID <- treeIn[,sampID[i]]
    else treeIn$SAMPID <- paste(treeIn$SAMPID,treeIn[,sampID[i]],sep='.')
  }
  samples <- unique(subset(treeIn,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(treeIn$SAMPID),stringsAsFactors=F)
  
  ##### TREE SPECIES METRICS ####################################################################################
  tcvr <- subset(treeIn,PARAMETER %in% c('VSMALL_TREE','SMALL_TREE','LMED_TREE','HMED_TREE','TALL_TREE','VTALL_TREE'))

  if(nrow(tcvr)>0){
    tspp <- reshape(subset(treeIn, PARAMETER=='TREE_SPECIES', 
                           select = c('SAMPID','PAGE','LINE','PLOT','PARAMETER','RESULT')), 
                    idvar = c('SAMPID','PAGE','LINE','PLOT'), direction = 'wide',
                    timevar = 'PARAMETER', v.names = 'RESULT')
    names(tspp) <- gsub('RESULT\\.', '', names(tspp)) 
 
    tcvr1 <- merge(tspp,tcvr,by=c('SAMPID','PLOT','PAGE','LINE'),all=TRUE)
    tcvr1$PARAM_ALT <- NA
    tcvr1$PARAM_ALT[tcvr1$PARAMETER %in% c('VSMALL_TREE','SMALL_TREE')] <- 'TREE_GROUND'
    tcvr1$PARAM_ALT[tcvr1$PARAMETER %in% c('LMED_TREE','HMED_TREE')] <- 'TREE_MID'
    tcvr1$PARAM_ALT[tcvr1$PARAMETER %in% c('TALL_TREE','VTALL_TREE')] <- 'TREE_UPPER'

    ntrspp <- aggregate(x = list(N_TREESPP = tcvr1$TREE_SPECIES), 
                        by = tcvr1[c('SAMPID')], FUN = function(x){length(unique(x))})
    totspp <- merge(tcvr1, ntrspp, by = c('SAMPID'))
    totspp1 <- reshape(ntrspp, idvar = 'SAMPID', direction = 'long',
                       varying = 'N_TREESPP', times = 'N_TREESPP',
                       timevar = 'METRIC', v.names = 'RESULT')
    

    totspp.pos <- subset(totspp, RESULT!='0' & !is.na(RESULT), 
                         select = c('SAMPID','PARAMETER','PARAM_ALT','TREE_SPECIES'))
    
    tspp1a <- aggregate(x = list(N = totspp.pos$TREE_SPECIES), 
                        by = totspp.pos[c('SAMPID','PARAMETER')],
                        FUN = function(x){length(unique(x))})
    tspp1a <- subset(tspp1a, !is.na(PARAMETER))

    tspp1b <- aggregate(x = list(N = totspp.pos$TREE_SPECIES), 
                        by = totspp.pos[c('SAMPID','PARAM_ALT')],
                        FUN = function(x){length(unique(x))})
    tspp1b <- subset(tspp1b, !is.na(PARAM_ALT))
    
    tspp1c <- merge(tspp1b, ntrspp, by='SAMPID')
    tspp1c <- subset(tspp1c, !is.na(PARAM_ALT))
    
    tspp1c$PCTN <- with(tspp1c, round((N/N_TREESPP)*100, 2))
    tspp1c$N_TREESPP <- NULL
    

    tspp2a <- reshape(tspp1a, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                      varying = 'N', times = 'N',
                      timevar = 'METRIC', v.names = 'RESULT')
    tspp2a$METRIC <- with(tspp2a, paste(METRIC, PARAMETER, sep='_'))

    tspp2c <- reshape(tspp1c, idvar = c('SAMPID','PARAM_ALT'), direction = 'long',
                       varying = c('PCTN','N'), times = c('PCTN','N'),
                       timevar = 'METRIC', v.names = 'RESULT')
    tspp2c$METRIC <- with(tspp2c, paste(METRIC, PARAM_ALT, sep='_'))

    tsppOut <- rbind(totspp1,subset(tspp2a,select=-PARAMETER),subset(tspp2c,select=-PARAM_ALT))
    tsppOut1 <- reshape(tsppOut, idvar = c('SAMPID'), direction = 'wide',
                        timevar = 'METRIC', v.names = 'RESULT')
    names(tsppOut1) <- gsub("RESULT\\.", "", names(tsppOut1))

    empty_tspp <- data.frame(t(rep(NA,13)))
    names(empty_tspp) <- c("N_TREESPP","N_TALL_TREE","N_HMED_TREE","N_LMED_TREE","N_SMALL_TREE",
                           "N_VSMALL_TREE","N_VTALL_TREE","N_TREE_UPPER",
                           "N_TREE_MID","N_TREE_GROUND","PCTN_TREE_UPPER","PCTN_TREE_MID",
                           "PCTN_TREE_GROUND")

    tsppOut2 <- subset(merge(tsppOut1, empty_tspp, all=TRUE),!is.na(SAMPID))

    tsppOut3 <- merge(allUIDs,tsppOut2,by='SAMPID',all.x=T)

    varNames <- names(tsppOut3)[!names(tsppOut3) %in% c('SAMPID')]
    tsppOut4 <- reshape(tsppOut3, idvar = 'SAMPID', direction = 'long',
                        varying = varNames, times = varNames,
                        timevar = 'METRIC', v.names = 'RESULT')
    tsppOut4$METRIC <- with(tsppOut4, as.character(METRIC))
    tsppOut4$RESULT <- with(tsppOut4, ifelse(is.na(RESULT),0,RESULT))
  }else{
    empty_tspp <- data.frame(t(rep(NA,13)),stringsAsFactors=F)
    names(empty_tspp) <- c("N_TREESPP","N_TALL_TREE","N_HMED_TREE","N_LMED_TREE","N_SMALL_TREE","N_VSMALL_TREE","N_VTALL_TREE","N_TREE_UPPER"
                           ,"N_TREE_MID","N_TREE_GROUND","PCTN_TREE_UPPER","PCTN_TREE_MID","PCTN_TREE_GROUND")

    tsppOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)),stringsAsFactors=F), empty_tspp, all=TRUE)
    tsppOut1 <- subset(tsppOut,!is.na(SAMPID))
    
    varNames <- names(tsppOut1)[!names(tsppOut1) %in% c('SAMPID')]
    tsppOut2 <- reshape(tsppOut1, idvar = 'SAMPID', direction = 'long',
                        varying = varNames, times = varNames,
                        timevar = 'METRIC', v.names = 'RESULT')
    tsppOut2$METRIC <- with(tsppOut2, as.character(METRIC))
    tsppOut2$RESULT <- 0
    tsppOut4 <- tsppOut2
  }
  print("Done with tree species metrics")

  ## TREE COVER METRICS ###########################################################################################
  ## Sum by species within plot
  if(nrow(tcvr)>0){
    tcvr2a <- aggregate(x = list(COV = tcvr1$RESULT), 
                        by = tcvr1[c('SAMPID','PLOT','NPLOTS','PARAMETER','TREE_SPECIES')],
                        FUN = function(x){sum(as.numeric(x))})
    tcvr2a$COV <- ifelse(tcvr2a$COV>100,100,tcvr2a$COV)
    tcvr2a <- subset(tcvr2a, COV!=0)
    
    tcvr2b <- aggregate(x = list(COV = tcvr1$RESULT), 
                        by = tcvr1[c('SAMPID','PLOT','NPLOTS','PARAM_ALT','TREE_SPECIES')],
                        FUN = function(x){sum(as.numeric(x))})
    
    tcvr2b$COV <- ifelse(tcvr2b$COV>100,100,tcvr2b$COV)
    tcvr2b <- subset(tcvr2b, COV!=0)

    tcvr3a.uniq <- aggregate(x = list(uniqN = tcvr2a$NPLOTS),
                             by = tcvr2a[c('SAMPID','PARAMETER')], FUN = unique)
    
    tcvr3a.length <- aggregate(x = list(uniqPlot = tcvr2a$PLOT),
                               by = tcvr2a[c('SAMPID','PARAMETER')],
                               FUN = function(x){length(unique(x))})
    tcvr3a.sum <- aggregate(x = list(sumcov = tcvr2a$COV),
                            by = tcvr2a[c('SAMPID','PARAMETER')],
                            FUN = function(x){sum(as.numeric(x))})
    tcvr3a <- merge(tcvr3a.length, tcvr3a.uniq, by = c('SAMPID','PARAMETER'))
    tcvr3a <- merge(tcvr3a, tcvr3a.sum, by = c('SAMPID','PARAMETER'))
    
    tcvr3a$FREQ <- with(tcvr3a, round((uniqPlot/uniqN)*100, 2))
    tcvr3a$XCOV <- with(tcvr3a, round(sumcov/uniqN, 2))
    tcvr3a$IMP <- with(tcvr3a, round((FREQ + XCOV)/2, 2))
  
    tcvr3a <- subset(tcvr3a, select = c('SAMPID','PARAMETER','FREQ','XCOV','IMP'))
    
    tcvr3b.uniq <- aggregate(x = list(uniqN = tcvr2b$NPLOTS),
                             by = tcvr2b[c('SAMPID','PARAM_ALT')], FUN = unique)
    
    tcvr3b.length <- aggregate(x = list(uniqPlot = tcvr2b$PLOT),
                               by = tcvr2b[c('SAMPID','PARAM_ALT')],
                               FUN = function(x){length(unique(x))})
    tcvr3b.sum <- aggregate(x = list(sumcov = tcvr2b$COV),
                            by = tcvr2b[c('SAMPID','PARAM_ALT')],
                            FUN = function(x){sum(as.numeric(x))})
    tcvr3b <- merge(tcvr3b.length, tcvr3b.uniq, by = c('SAMPID','PARAM_ALT'))
    tcvr3b <- merge(tcvr3b, tcvr3b.sum, by = c('SAMPID','PARAM_ALT'))
    
    tcvr3b$FREQ <- with(tcvr3b, round((uniqPlot/uniqN)*100, 2))
    tcvr3b$XCOV <- with(tcvr3b, round(sumcov/uniqN, 2))
    tcvr3b$IMP <- with(tcvr3b, round((FREQ + XCOV)/2, 2))
    
    tcvr3b <- subset(tcvr3b, select = c('SAMPID','PARAM_ALT','FREQ','XCOV','IMP'))
    
    varNames.a <- names(tcvr3a)[!names(tcvr3a) %in% c('SAMPID','PARAMETER')]
    tcvr4a <- reshape(tcvr3a, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                      varying = varNames.a, times = varNames.a,
                      timevar = 'METRIC', v.names = 'RESULT')
    tcvr4a$METRIC = with(tcvr4a, paste(METRIC, PARAMETER, sep='_'))

    varNames.b <- names(tcvr3b)[!names(tcvr3b) %in% c('SAMPID','PARAM_ALT')]
    tcvr4b <- reshape(tcvr3b, idvar = c('SAMPID','PARAM_ALT'), direction = 'long',
                      varying = varNames.b, times = varNames.b,
                      timevar = 'METRIC', v.names = 'RESULT')
    tcvr4b$METRIC = with(tcvr4b, paste(METRIC, PARAM_ALT, sep='_'))

    tcvrOut <- rbind(subset(tcvr4a,select=-PARAMETER),subset(tcvr4b,select=-PARAM_ALT))
    tcvrOut1 <- reshape(tcvrOut, idvar = 'SAMPID', direction = 'wide',
                        timevar = 'METRIC', v.names = 'RESULT')
    names(tcvrOut1) <- gsub("RESULT\\.", '', names(tcvrOut1))

    empty_tcvr <- data.frame(t(rep(NA,27)),stringsAsFactors=F)
    names(empty_tcvr) <- c("FREQ_TALL_TREE","FREQ_HMED_TREE","FREQ_LMED_TREE","FREQ_SMALL_TREE",
                           "FREQ_VSMALL_TREE","FREQ_VTALL_TREE","XCOV_TALL_TREE",
                           "XCOV_HMED_TREE","XCOV_LMED_TREE","XCOV_SMALL_TREE","XCOV_VSMALL_TREE",
                           "XCOV_VTALL_TREE","IMP_TALL_TREE","IMP_HMED_TREE",
                           "IMP_LMED_TREE","IMP_SMALL_TREE","IMP_VSMALL_TREE","IMP_VTALL_TREE",
                           "FREQ_TREE_UPPER","FREQ_TREE_MID","FREQ_TREE_GROUND",
                           "XCOV_TREE_UPPER","XCOV_TREE_MID","XCOV_TREE_GROUND",
                           "IMP_TREE_UPPER","IMP_TREE_MID","IMP_TREE_GROUND")

    tcvrOut2 <- subset(merge(tcvrOut1, empty_tcvr, all=TRUE),!is.na(SAMPID))

    tcvrOut3 <- merge(allUIDs,tcvrOut2,by='SAMPID',all.x=T)

    varNames <- names(tcvrOut3)[!names(tcvrOut3) %in% c('SAMPID')]
    tcvrOut4 <- reshape(tcvrOut3, idvar = 'SAMPID', direction = 'long',
                        varying = varNames, times = varNames,
                        timevar = 'METRIC', v.names = 'RESULT')
    
    tcvrOut4$METRIC <- with(tcvrOut4, as.character(METRIC))
    tcvrOut4$RESULT <- with(tcvrOut4, ifelse(is.na(RESULT),0,RESULT))

  }else{
    empty_tcvr <- data.frame(t(rep(NA,27)),stringsAsFactors=F)
    names(empty_tcvr) <- c("FREQ_TALL_TREE","FREQ_HMED_TREE","FREQ_LMED_TREE","FREQ_SMALL_TREE",
                           "FREQ_VSMALL_TREE","FREQ_VTALL_TREE","XCOV_TALL_TREE",
                           "XCOV_HMED_TREE","XCOV_LMED_TREE","XCOV_SMALL_TREE","XCOV_VSMALL_TREE",
                           "XCOV_VTALL_TREE","IMP_TALL_TREE","IMP_HMED_TREE",
                           "IMP_LMED_TREE","IMP_SMALL_TREE","IMP_VSMALL_TREE","IMP_VTALL_TREE",
                           "FREQ_TREE_UPPER","FREQ_TREE_MID","FREQ_TREE_GROUND",
                           "XCOV_TREE_UPPER","XCOV_TREE_MID","XCOV_TREE_GROUND",
                           "IMP_TREE_UPPER","IMP_TREE_MID","IMP_TREE_GROUND")

    tcvrOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)),stringsAsFactors=F), empty_tcvr, all=TRUE)
    tcvrOut1 <- subset(tcvrOut,!is.na(SAMPID))
    
    varNames <- names(tcvrOut1)[!names(tcvrOut1) %in% c('SAMPID')]
    tcvrOut2 <- reshape(tcvrOut1, idvar = 'SAMPID', direction = 'long',
                        varying = varNames, times = varNames,
                        timevar = 'METRIC', v.names = 'RESULT')
    
    tcvrOut2$METRIC <- with(tcvrOut2, as.character(METRIC))
    tcvrOut2$RESULT <- with(tcvrOut2, ifelse(is.na(RESULT),0,RESULT))
    
    tcvrOut4 <- tcvrOut2
  }
  print("Done with tree cover metrics")
  ## Now fill in zeroes in each set of metrics and then combine
  treeOut <- rbind(tsppOut4,tcvrOut4) 
  treeOut$PARAMETER <- treeOut$METRIC
  treeOut$METRIC <- NULL
  
  treeOut.1 <- merge(samples,treeOut,by='SAMPID') 
  treeOut.1 <- subset(treeOut.1, select = -SAMPID)
  
  return(treeOut.1)
}
