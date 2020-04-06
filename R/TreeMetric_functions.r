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
#' nplots <- plyr::ddply(TreesEx,c('UID'),dplyr::summarise
#' ,NPLOTS=length(unique(PLOT)))
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

#  treeIn <- plyr::ddply(treeIn,c('SAMPID'),mutate,NPLOTS=length(unique(PLOT)))
  snags <- subset(treeIn,PARAMETER %in% c('XXTHIN_SNAG','XTHIN_SNAG','THIN_SNAG','JR_SNAG','THICK_SNAG','XTHICK_SNAG','XXTHICK_SNAG'))
  if(nrow(snags)>0){
    snagsOut <- plyr::ddply(snags,c('SAMPID','PARAMETER'),summarise,TOTN=sum(as.numeric(RESULT)),XN=round(TOTN/unique(NPLOTS),2))
    snagsOut1 <- reshape2::melt(snagsOut,id.vars=c('SAMPID','PARAMETER'),variable.name='METRIC',value.name='RESULT')
    snagsOut1 <- plyr::mutate(snagsOut1,METRIC=paste(as.character(METRIC),PARAMETER,sep='_'))

    totsnags <- plyr::ddply(snagsOut,c('SAMPID'),summarise,TOTN_SNAGS=sum(TOTN),XN_SNAGS=round(sum(XN),2))
    totsnags1 <- reshape2::melt(totsnags,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')

    allSnagsOut <- rbind(totsnags1,subset(snagsOut1,select=-PARAMETER))
    allSnagsOut1 <- reshape2::dcast(allSnagsOut,SAMPID~METRIC,value.var='RESULT')

    empty_snags <- data.frame(t(rep(NA,14)),stringsAsFactors=F)

    names(empty_snags) <- c("TOTN_SNAGS","XN_SNAGS","TOTN_JR_SNAG","TOTN_THICK_SNAG","TOTN_THIN_SNAG","TOTN_XTHICK_SNAG","TOTN_XTHIN_SNAG"
                            ,"TOTN_XXTHIN_SNAG","XN_JR_SNAG","XN_THICK_SNAG","XN_THIN_SNAG","XN_XTHICK_SNAG","XN_XTHIN_SNAG","XN_XXTHIN_SNAG")

    allSnagsOut2 <- subset(merge(allSnagsOut1, empty_snags, all=TRUE),!is.na(SAMPID))

    # Merge with the all UIDs and fill in missing values with zeroes
    allSnagsOut3 <- merge(allUIDs,allSnagsOut2,by='SAMPID',all.x=T)

    allSnagsOut4 <- reshape2::melt(allSnagsOut3,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')
    allSnagsOut4 <- plyr::mutate(allSnagsOut4,METRIC=as.character(METRIC),RESULT=ifelse(is.na(RESULT),0,RESULT))
  }else{
    empty_snags <- data.frame(t(rep(NA,14)),stringsAsFactors=F)

    names(empty_snags) <- c("TOTN_SNAGS","XN_SNAGS","TOTN_JR_SNAG","TOTN_THICK_SNAG","TOTN_THIN_SNAG","TOTN_XTHICK_SNAG","TOTN_XTHIN_SNAG"
                            ,"TOTN_XXTHIN_SNAG","XN_JR_SNAG","XN_THICK_SNAG","XN_THIN_SNAG","XN_XTHICK_SNAG","XN_XTHIN_SNAG","XN_XXTHIN_SNAG")

    allSnagsOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)),stringsAsFactors=F), empty_snags, all=TRUE)

    allSnagsOut1 <- subset(allSnagsOut,!is.na(SAMPID))

    allSnagsOut2 <- reshape2::melt(allSnagsOut1,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')

    allSnagsOut3 <- plyr::mutate(allSnagsOut2,METRIC=as.character(METRIC),RESULT=0)

    allSnagsOut4 <- allSnagsOut3
  }
  print("Done with snag metrics")

  treeOut <- merge(samples, allSnagsOut4, by='SAMPID') %>% 
    plyr::rename(c('METRIC'='PARAMETER')) %>%
    dplyr::select(-SAMPID)

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
#' nplots <- plyr::ddply(TreesEx,c('UID'),dplyr::summarise
#' ,NPLOTS=length(unique(PLOT)))
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
  
#  treeIn <- plyr::ddply(treeIn,c('SAMPID'),mutate,NPLOTS=length(unique(PLOT)))
 
  tcnts <- subset(treeIn,PARAMETER %in% c('XXTHIN_TREE','XTHIN_TREE','THIN_TREE','JR_TREE','THICK_TREE','XTHICK_TREE','XXTHICK_TREE'))
  if(nrow(tcnts)>0){
    tcntsOut <- plyr::ddply(tcnts,c('SAMPID','PARAMETER'),summarise,TOTN=sum(as.numeric(RESULT)),XN=round(TOTN/unique(NPLOTS),2))
    tcntsOut1 <- reshape2::melt(tcntsOut,id.vars=c('SAMPID','PARAMETER'),variable.name='METRIC',value.name='RESULT')
    tcntsOut1 <- plyr::mutate(tcntsOut1,METRIC=paste(as.character(METRIC),PARAMETER,sep='_'))

    tottrees <- plyr::ddply(tcntsOut,c('SAMPID'),summarise,TOTN_TREES=sum(TOTN),XN_TREES=round(sum(XN),2))
    tottrees1 <- reshape2::melt(tottrees,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')

    allTreesOut <- rbind(tottrees1,subset(tcntsOut1,select=-PARAMETER))
    allTreesOut1 <- reshape2::dcast(allTreesOut,SAMPID~METRIC,value.var='RESULT')

    empty_trees <- data.frame(t(rep(NA,16)),stringsAsFactors=F)
    names(empty_trees) <- c("TOTN_TREES","XN_TREES","TOTN_JR_TREE","TOTN_THICK_TREE","TOTN_THIN_TREE","TOTN_XTHICK_TREE","TOTN_XTHIN_TREE"
                            ,"TOTN_XXTHICK_TREE","TOTN_XXTHIN_TREE","XN_JR_TREE","XN_THICK_TREE","XN_THIN_TREE","XN_XTHICK_TREE","XN_XTHIN_TREE"
                            ,"XN_XXTHICK_TREE","XN_XXTHIN_TREE")

    allTreesOut2 <- subset(merge(allTreesOut1, empty_trees, all=TRUE),!is.na(SAMPID))

    allTreesOut3 <- merge(allUIDs,allTreesOut2,by='SAMPID',all.x=T)

    allTreesOut4 <- reshape2::melt(allTreesOut3,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')
    allTreesOut4 <- plyr::mutate(allTreesOut4,METRIC=as.character(METRIC),RESULT=ifelse(is.na(RESULT),0,RESULT))
  }else{
    empty_trees <- data.frame(t(rep(NA,16)),stringsAsFactors=F)
    names(empty_trees) <- c("TOTN_TREES","XN_TREES","TOTN_JR_TREE","TOTN_THICK_TREE","TOTN_THIN_TREE","TOTN_XTHICK_TREE","TOTN_XTHIN_TREE"
                            ,"TOTN_XXTHICK_TREE","TOTN_XXTHIN_TREE","XN_JR_TREE","XN_THICK_TREE","XN_THIN_TREE","XN_XTHICK_TREE","XN_XTHIN_TREE"
                            ,"XN_XXTHICK_TREE","XN_XXTHIN_TREE")

    allTreesOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)),stringsAsFactors=F), empty_trees, all=TRUE)

    allTreesOut1 <- subset(allTreesOut,!is.na(SAMPID))

    allTreesOut2 <- reshape2::melt(allTreesOut1,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')

    allTreesOut3 <- plyr::mutate(allTreesOut2,METRIC=as.character(METRIC),RESULT=0)

    allTreesOut4 <- allTreesOut3
  }
  print("Done with tree count metrics")
  treeOut <- merge(samples, allTreesOut4, by='SAMPID') %>% 
    plyr::rename(c('METRIC'='PARAMETER')) %>%
    dplyr::select(-SAMPID)
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
#' nplots <- plyr::ddply(TreesEx,c('UID'),dplyr::summarise
#' ,NPLOTS=length(unique(PLOT)))
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
  
#  treeIn <- plyr::ddply(treeIn,c('SAMPID'),mutate,NPLOTS=length(unique(PLOT)))
  ##### TREE SPECIES METRICS ####################################################################################
  tcvr <- subset(treeIn,PARAMETER %in% c('VSMALL_TREE','SMALL_TREE','LMED_TREE','HMED_TREE','TALL_TREE','VTALL_TREE'))
  if(nrow(tcvr)>0){
    tspp <- reshape2::dcast(subset(treeIn,PARAMETER=='TREE_SPECIES'),SAMPID+PAGE+LINE+PLOT~PARAMETER,value.var='RESULT')

    tcvr1 <- merge(tspp,tcvr,by=c('SAMPID','PLOT','PAGE','LINE'),all=TRUE)
    tcvr1 <- plyr::mutate(tcvr1,PARAM_ALT=car::Recode(PARAMETER,"c('VSMALL_TREE','SMALL_TREE')='TREE_GROUND';c('LMED_TREE','HMED_TREE')='TREE_MID';
                                           c('TALL_TREE','VTALL_TREE')='TREE_UPPER'"))

    totspp <- plyr::ddply(tcvr1,c('SAMPID'),mutate,N_TREESPP=length(unique(TREE_SPECIES)))
    totspp1 <- unique(reshape2::melt(totspp,id.vars='SAMPID',measure.vars='N_TREESPP',variable.name='METRIC',value.name='RESULT'))

    tspp1a <- subset(plyr::ddply(subset(totspp,RESULT!='0'),c('SAMPID','PARAMETER'),summarise,N=length(unique(TREE_SPECIES))),!is.na(PARAMETER))

    tspp1b <- subset(plyr::ddply(subset(totspp,RESULT!='0'),c('SAMPID','PARAM_ALT'),summarise,N=length(unique(TREE_SPECIES))
                           ,PCTN=round((N/unique(N_TREESPP))*100,2)),!is.na(PARAM_ALT))

    tspp2a <- reshape2::melt(tspp1a,id.vars=c('SAMPID','PARAMETER'),variable.name='METRIC',value.name='RESULT')
    tspp2a <- plyr::mutate(tspp2a,METRIC=paste(METRIC,PARAMETER,sep='_'))

    tspp2b <- reshape2::melt(tspp1b,id.vars=c('SAMPID','PARAM_ALT'),variable.name='METRIC',value.name='RESULT')
    tspp2b <- plyr::mutate(tspp2b,METRIC=paste(METRIC,PARAM_ALT,sep='_'))

    tsppOut <- rbind(totspp1,subset(tspp2a,select=-PARAMETER),subset(tspp2b,select=-PARAM_ALT))
    tsppOut1 <- reshape2::dcast(tsppOut,SAMPID~METRIC,value.var='RESULT')

    empty_tspp <- data.frame(t(rep(NA,13)))
    names(empty_tspp) <- c("N_TREESPP","N_TALL_TREE","N_HMED_TREE","N_LMED_TREE","N_SMALL_TREE","N_VSMALL_TREE","N_VTALL_TREE","N_TREE_UPPER"
                           ,"N_TREE_MID","N_TREE_GROUND","PCTN_TREE_UPPER","PCTN_TREE_MID","PCTN_TREE_GROUND")

    tsppOut2 <- subset(merge(tsppOut1, empty_tspp, all=TRUE),!is.na(SAMPID))

    tsppOut3 <- merge(allUIDs,tsppOut2,by='SAMPID',all.x=T)

    tsppOut4 <- reshape2::melt(tsppOut3,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')
    tsppOut4 <- plyr::mutate(tsppOut4,METRIC=as.character(METRIC),RESULT=ifelse(is.na(RESULT),0,RESULT))
  }else{
    empty_tspp <- data.frame(t(rep(NA,13)),stringsAsFactors=F)
    names(empty_tspp) <- c("N_TREESPP","N_TALL_TREE","N_HMED_TREE","N_LMED_TREE","N_SMALL_TREE","N_VSMALL_TREE","N_VTALL_TREE","N_TREE_UPPER"
                           ,"N_TREE_MID","N_TREE_GROUND","PCTN_TREE_UPPER","PCTN_TREE_MID","PCTN_TREE_GROUND")

    tsppOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)),stringsAsFactors=F), empty_tspp, all=TRUE)
    tsppOut1 <- subset(tsppOut,!is.na(SAMPID))
    tsppOut2 <- reshape2::melt(tsppOut1,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')
    tsppOut3 <- plyr::mutate(tsppOut2,METRIC=as.character(METRIC),RESULT=0)
    tsppOut4 <- tsppOut3
  }
  print("Done with tree species metrics")

  ## TREE COVER METRICS ###########################################################################################
  ## Sum by species within plot
  if(nrow(tcvr)>0){
    tcvr2a <- plyr::ddply(tcvr1,c('SAMPID','PLOT','NPLOTS','PARAMETER','TREE_SPECIES'),summarise,COV=sum(as.numeric(RESULT)))
    tcvr2a$COV <- ifelse(tcvr2a$COV>100,100,tcvr2a$COV)
    tcvr2b <- plyr::ddply(tcvr1,c('SAMPID','PLOT','NPLOTS','PARAM_ALT','TREE_SPECIES'),summarise,COV=sum(as.numeric(RESULT)))
    tcvr2b$COV <- ifelse(tcvr2b$COV>100,100,tcvr2b$COV)

    tcvr3a <- plyr::ddply(subset(tcvr2a,COV!=0),c('SAMPID','PARAMETER'),summarise,FREQ=round((length(unique(PLOT))/unique(NPLOTS))*100,2)
                    ,XCOV=round((sum(as.numeric(COV))/unique(NPLOTS)),2),IMP=round((FREQ+XCOV)/2,2))
    tcvr3b <- plyr::ddply(subset(tcvr2b,COV!=0),c('SAMPID','PARAM_ALT'),summarise,FREQ=round((length(unique(PLOT))/unique(NPLOTS))*100,2)
                    ,XCOV=round((sum(as.numeric(COV))/unique(NPLOTS)),2),IMP=round((FREQ+XCOV)/2,2))

    tcvr4a <- reshape2::melt(tcvr3a,id.vars=c('SAMPID','PARAMETER'),variable.name='METRIC',value.name='RESULT')
    tcvr4a <- plyr::mutate(tcvr4a,METRIC=paste(METRIC,PARAMETER,sep='_'))
    tcvr4b <- reshape2::melt(tcvr3b,id.vars=c('SAMPID','PARAM_ALT'),variable.name='METRIC',value.name='RESULT')
    tcvr4b <- plyr::mutate(tcvr4b,METRIC=paste(METRIC,PARAM_ALT,sep='_'))

    tcvrOut <- rbind(subset(tcvr4a,select=-PARAMETER),subset(tcvr4b,select=-PARAM_ALT))
    tcvrOut1 <- reshape2::dcast(tcvrOut,SAMPID~METRIC,value.var='RESULT')

    empty_tcvr <- data.frame(t(rep(NA,27)),stringsAsFactors=F)
    names(empty_tcvr) <- c("FREQ_TALL_TREE","FREQ_HMED_TREE","FREQ_LMED_TREE","FREQ_SMALL_TREE","FREQ_VSMALL_TREE","FREQ_VTALL_TREE","XCOV_TALL_TREE"
                           ,"XCOV_HMED_TREE","XCOV_LMED_TREE","XCOV_SMALL_TREE","XCOV_VSMALL_TREE","XCOV_VTALL_TREE","IMP_TALL_TREE","IMP_HMED_TREE"
                           ,"IMP_LMED_TREE","IMP_SMALL_TREE","IMP_VSMALL_TREE","IMP_VTALL_TREE","FREQ_TREE_UPPER","FREQ_TREE_MID","FREQ_TREE_GROUND"
                           ,"XCOV_TREE_UPPER","XCOV_TREE_MID","XCOV_TREE_GROUND","IMP_TREE_UPPER","IMP_TREE_MID","IMP_TREE_GROUND")

    tcvrOut2 <- subset(merge(tcvrOut1, empty_tcvr, all=TRUE),!is.na(SAMPID))

    tcvrOut3 <- merge(allUIDs,tcvrOut2,by='SAMPID',all.x=T)

    tcvrOut4 <- reshape2::melt(tcvrOut3,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')
    tcvrOut4 <- plyr::mutate( tcvrOut4,METRIC=as.character(METRIC),RESULT=ifelse(is.na(RESULT),0,RESULT))
    
  }else{
    empty_tcvr <- data.frame(t(rep(NA,27)),stringsAsFactors=F)
    names(empty_tcvr) <- c("FREQ_TALL_TREE","FREQ_HMED_TREE","FREQ_LMED_TREE","FREQ_SMALL_TREE","FREQ_VSMALL_TREE","FREQ_VTALL_TREE","XCOV_TALL_TREE"
                           ,"XCOV_HMED_TREE","XCOV_LMED_TREE","XCOV_SMALL_TREE","XCOV_VSMALL_TREE","XCOV_VTALL_TREE","IMP_TALL_TREE","IMP_HMED_TREE"
                           ,"IMP_LMED_TREE","IMP_SMALL_TREE","IMP_VSMALL_TREE","IMP_VTALL_TREE","FREQ_TREE_UPPER","FREQ_TREE_MID","FREQ_TREE_GROUND"
                           ,"XCOV_TREE_UPPER","XCOV_TREE_MID","XCOV_TREE_GROUND","IMP_TREE_UPPER","IMP_TREE_MID","IMP_TREE_GROUND")

    tcvrOut <- merge(data.frame(SAMPID=rep(unique(treeIn$SAMPID)),stringsAsFactors=F), empty_tcvr, all=TRUE)
    tcvrOut1 <- subset(tcvrOut,!is.na(SAMPID))
    tcvrOut2 <- reshape2::melt(tcvrOut1,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')
    tcvrOut3 <- plyr::mutate(tcvrOut2,METRIC=as.character(METRIC),RESULT=0)
    tcvrOut4 <- tcvrOut3
  }
  print("Done with tree cover metrics")
  ## Now fill in zeroes in each set of metrics and then combine
  treeOut <- rbind(tsppOut4,tcvrOut4) %>% 
    plyr::rename(c('METRIC'='PARAMETER'))
  
  treeOut.1 <- merge(samples,treeOut,by='SAMPID') %>%
    dplyr::select(-SAMPID)
  
  return(treeOut.1)
}
