#' @export
#' 
#' @title Calculate metrics based on S&T (Status and Trends) categories
#' 
#' @description This function calculations metrics based on the
#' number of S&T classes found across plots. This function is called by
#' calcVtype_GcovMets().
#' 
#' @param dataIn A data frame containing the following variables:
#' \itemize{
#' \item sampID: Variables identified in \emph{sampID} argument
#'
#' \item PLOT: Sample plot from which data were collected
#'
#' \item PARAMETER: specific measurement type
#'
#' \item RESULT: measured value
#' }
#' The following parameters are used in calculating
#' vegetation type metrics: 'SANDT_CLASS', 'PAL_FARMED'.
#' Additional parameters or variables are ignored.
#' @param nPlot A data frame with the 
#' number of plots sampled associated with
#' each sample, with \emph{sampID} variables and NPLOTS
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID'
#' by default.
#' 
#' @details If any of the parameters are missing, they are assumed
#' to be zeros (if numeric), and metric values associated with any
#' metrics that cannot be calculated due to missing parameters are
#' set to a standardized value.
#' 
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. Data frame contains
#'   \emph{sampID} variables, PARAMETER, RESULT, where values of PARAMETER 
#'   consist of the metric name concatenated with trait value, with valid values
#'   of: N_SANDT, DOM_SANDT, D_SANDT, H_SANDT, J_SANDT. A list of metric
#'   descriptions is provided in the document named 
#'   VegTypes_GrdCover_Metric_Descriptions.pdf included in the help directory
#'   for the package.
#'   
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(Vtype_GrCovEx)
#' # Create data frame with number of plots sampled for each sampID
#' nplots <- data.frame(UID=seq(1:10),NPLOTS=rep(5,10),stringsAsFactors=F)
#' # alternative approach to creating this data frame
#' nplots <- plyr::ddply(Vtype_GrCovEx,c('UID'),dplyr::summarise
#' ,NPLOTS=length(unique(PLOT)))
#'
#' sandtEx <- calcSandTMets(Vtype_GrCovEx,nplots)
#'
#' head(sandtEx)
#' unique(sandtEx$PARAMETER)

calcSandTMets <- function(dataIn,nPlot,sampID='UID'){
  dataIn1 <- merge(dataIn,nPlot,by=sampID)
  
  # Create vector of all samples in dataset
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID,dataIn1[,sampID[i]],sep='.')
  }
  samples <- unique(subset(dataIn1,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID),stringsAsFactors=F)
  
  vhet <- reshape2::dcast(dataIn1,SAMPID+PLOT+NPLOTS~PARAMETER,value.var='RESULT')

  ## First calculate heterogeneity metrics ## Account for datasets without PF parameter so that wide form will not have it
  if('PAL_FARMED' %in% names(vhet)==TRUE){
    vhet <- mutate(vhet,SANDT_CLASS=ifelse(is.na(PAL_FARMED),SANDT_CLASS,'PF'))
  }
  vhet1 <- plyr::ddply(subset(vhet,!is.na(SANDT_CLASS)),c('SAMPID'),mutate,N_SANDT=length(unique(SANDT_CLASS)))
  vhet2 <- plyr::ddply(vhet1,c('SAMPID','NPLOTS','SANDT_CLASS','N_SANDT'),summarise,FREQ=length(SANDT_CLASS)/unique(NPLOTS))
  vhet3 <- plyr::ddply(vhet2,c('SAMPID'),mutate,MAXF=max(FREQ))
  vhet4 <- subset(vhet3,MAXF==FREQ)
  vhet4a <- plyr::ddply(vhet4,c('SAMPID','N_SANDT'),summarise,DOM_SANDT=paste(SANDT_CLASS,collapse='-'))

  vhet5 <- plyr::ddply(vhet2,c('SAMPID','NPLOTS'),summarise,D_SANDT=round(1-sum(FREQ*FREQ),4),H_SANDT=round(-1*(sum(FREQ*log(FREQ))),4)
                 ,J_SANDT=round(ifelse(H_SANDT!=0 & unique(N_SANDT)!=1,H_SANDT/log(unique(N_SANDT)),0),4))

  vhet6 <- merge(vhet4a,vhet5,by='SAMPID')
  # create an empty data frame with all of the metric names as variables to ensure all are included in output
  empty_vhet <- data.frame(t(rep(NA,5)),stringsAsFactors=FALSE)
  names(empty_vhet) <- c('N_SANDT','DOM_SANDT','D_SANDT','H_SANDT','J_SANDT')

  vhet7 <- subset(merge(vhet6, empty_vhet, all=TRUE),!is.na(SAMPID))

  vhet8 <- merge(allUIDs,vhet7,by='SAMPID',all.x=T)

  vhet9 <- reshape2::melt(vhet8,id.vars='SAMPID',measure.vars=c('D_SANDT','H_SANDT','J_SANDT','DOM_SANDT','N_SANDT'),variable.name='METRIC'
                ,value.name='RESULT')
  vhet9 <- mutate(vhet9,METRIC=as.character(METRIC),RESULT=ifelse(METRIC=='DOM_SANDT' & is.na(RESULT),'MISSING'
                                                                  ,ifelse(is.na(RESULT),0,RESULT)))

  print("Done with veg heterogeneity metrics")
  vhetOut <- merge(samples,vhet9,by='SAMPID') %>%
    plyr::rename(c('METRIC'='PARAMETER')) %>%
    dplyr::select(-SAMPID)
  return(vhetOut)
}


#### VASCULAR STRATA
#' @export
#' 
#' @title Calculate vascular strata metrics
#' 
#' @description Calculate vascular strata metrics using
#' vegetation type data. This function is called by
#' calcVtype_GcovMets().
#' 
#' @param dataIn A data frame containing the following
#' variables:
#' \itemize{
#' \item sampID: Variable(s) identified in the \emph{sampID} argument
#' 
#' \item PLOT: Sample plot from which data were collected
#'
#' \item PARAMETER: specific measurement type
#'
#' \item RESULT: measured value
#' }
#' The following parameters are used in
#' calculating tree metrics: 'SUBMERGED_AQ', 'FLOATING_AQ', 
#' 'LIANAS', 'VTALL_VEG', 'TALL_VEG', 'HMED_VEG', 'MED_VEG', 
#' 'SMALL_VEG', 'VSMALL_VEG'.
#' Additional parameters or variables are ignored.
#' @param nPlot A data frame with the 
#' number of plots sampled associated with each
#' sample with \emph{sampID} variables and NPLOTS
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples,
#' 'UID' by default.
#' 
#' @details If any of the parameters are missing, they are
#' assumed to be zeros (if numeric), and metric values associated
#' with any metrics that cannot be calculated due to missing
#' parameters are set to a standardized value.
#' 
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. Data frame contains
#'   \emph{sampID} variables, PARAMETER, RESULT, where values of PARAMETER 
#'   consist of the metric name concatenated with trait value, with valid values
#'   of: FREQ_FLOATING_AQ, FREQ_HMED_VEG, FREQ_LIANAS, FREQ_MED_VEG,
#'   FREQ_SMALL_VEG, FREQ_SUBMERGED_AQ, FREQ_TALL_VEG, FREQ_VSMALL_VEG,
#'   FREQ_VTALL_VEG, IMP_FLOATING_AQ, IMP_HMED_VEG, IMP_LIANAS, IMP_MED_VEG,
#'   IMP_SMALL_VEG, IMP_SUBMERGED_AQ, IMP_TALL_VEG, IMP_VSMALL_VEG,
#'   IMP_VTALL_VEG, XCOV_FLOATING_AQ, XCOV_HMED_VEG, XCOV_LIANAS, XCOV_MED_VEG, 
#'   XCOV_SMALL_VEG, XCOV_SUBMERGED_AQ, XCOV_TALL_VEG, XCOV_VSMALL_VEG,
#'   XCOV_VTALL_VEG, XRCOV_FLOATING_AQ, XRCOV_HMED_VEG, XRCOV_LIANAS,
#'   XRCOV_MED_VEG, XRCOV_SMALL_VEG, XRCOV_SUBMERGED_AQ, XRCOV_TALL_VEG,
#'   XRCOV_VSMALL_VEG, XRCOV_VTALL_VEG, H_VASC_STRATA, J_VASC_STRATA,
#'   D_VASC_STRATA. A list of metric descriptions is provided in the document
#'   VegTypes_GrdCover_Metric_Descriptions.pdf included in the help directory
#'   for the package.
#'   
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(Vtype_GrCovEx)
#' # Create data frame with number of plots sampled for each UID
#' nplots <- data.frame(UID=seq(1:10),NPLOTS=rep(5,10),stringsAsFactors=F)
#' # alternative approach to creating this data frame
#' nplots <- plyr::ddply(Vtype_GrCovEx,c('UID'),dplyr::summarise,
#' NPLOTS=length(unique(PLOT)))
#'
#' stratEx <- calcVascStratMets(Vtype_GrCovEx,nplots)
#'
#' head(stratEx)
#' unique(stratEx$METRIC)
calcVascStratMets <- function(dataIn,nPlot,sampID='UID'){

  dataIn1 <- merge(dataIn,nPlot,by=sampID)

  # Create vector of all samples in dataset
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID,dataIn1[,sampID[i]],sep='.')
  }
  samples <- unique(subset(dataIn1,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID),stringsAsFactors=F)

  vstrat <- subset(dataIn1,PARAMETER %in% c('SUBMERGED_AQ','FLOATING_AQ','LIANAS','VTALL_VEG','TALL_VEG','HMED_VEG','MED_VEG'
                                            ,'SMALL_VEG','VSMALL_VEG') & !is.na(RESULT) & RESULT!='0')
  vstrat <- plyr::ddply(vstrat,c('SAMPID'),mutate,N_VASC_STRATA=length(unique(PARAMETER))
                  ,XTOTCOV_VASC_STRATA=sum(as.numeric(RESULT))/NPLOTS,PLOTSAMP=length(unique(PLOT)))

  vstratPlot <- plyr::ddply(vstrat,c('SAMPID','PLOT','N_VASC_STRATA','NPLOTS','XTOTCOV_VASC_STRATA','PLOTSAMP'),summarise,N_VSTRATA=length(unique(PARAMETER))
                      ,SUM_PLOT=sum(as.numeric(RESULT)))


  vstrat1 <- unique(plyr::ddply(vstratPlot,c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA'),summarise
                                ,XN_VASC_STRATA=sum(N_VSTRATA)/unique(NPLOTS)
                          ,RG_VASC_STRATA=ifelse(unique(PLOTSAMP)==unique(NPLOTS)
                          ,max(N_VSTRATA)-min(N_VSTRATA),max(N_VSTRATA)-0)))

  empty_vstrat <- data.frame(t(rep(NA,4)),stringsAsFactors=FALSE)
  names(empty_vstrat) <- c('N_VASC_STRATA','XTOTCOV_VASC_STRATA','XN_VASC_STRATA','RG_VASC_STRATA')

  vstrat2 <- subset(merge(vstrat1, empty_vstrat, all=TRUE),!is.na(SAMPID))

  vstrat3 <- merge(allUIDs,vstrat2,by='SAMPID',all.x=T)

  vstratMet <- reshape2::melt(vstrat3,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')
  vstratMet <- mutate(vstratMet,METRIC=as.character(METRIC),RESULT=ifelse(is.na(RESULT),0,RESULT))

  ## Calculate frequency, mean cover, relative mean cover, and relative importance by vascular stratum
  indf1 <- plyr::ddply(subset(vstrat,RESULT!=0),c('SAMPID','PARAMETER'),summarise,FREQ=round((length(PLOT)/unique(NPLOTS))*100,2)
                 ,XCOV=round((sum(as.numeric(RESULT))/unique(NPLOTS)),2)
                 ,XRCOV=round((XCOV/unique(XTOTCOV_VASC_STRATA))*100,2),IMP=round((FREQ+XCOV)/2,2))

  outdf <- mutate(reshape2::melt(indf1,id.vars=c('SAMPID','PARAMETER'),value.name='RESULT'),PARAMETER=paste(variable,PARAMETER,sep='_'))

  outdf1 <- reshape2::melt(reshape2::dcast(outdf,SAMPID~PARAMETER,value.var='RESULT'),id.vars='SAMPID'
                           ,variable.name='METRIC',value.name='RESULT')
  outdf1$RESULT[is.na(outdf1$RESULT)] <- 0

  ## Calculate diversity indices
  div1 <- plyr::ddply(indf1,c('SAMPID'),summarise,H_VASC_STRATA=round(-1*sum((XRCOV/100)*log(XRCOV/100)),4)
                ,J_VASC_STRATA=round(H_VASC_STRATA/log(length(PARAMETER)),4)
                ,D_VASC_STRATA=round(1-sum((XRCOV/100)^2),4))
  div1 <- mutate(div1,J_VASC_STRATA=ifelse(!is.na(J_VASC_STRATA),J_VASC_STRATA,0))

  outdf2 <- rbind(outdf1,reshape2::melt(div1,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT'))

  empty_vtype <- data.frame(t(rep(NA,39)),stringsAsFactors=FALSE)
  names(empty_vtype) <- c("FREQ_FLOATING_AQ","FREQ_HMED_VEG","FREQ_LIANAS","FREQ_MED_VEG","FREQ_SMALL_VEG","FREQ_SUBMERGED_AQ"
                          ,"FREQ_TALL_VEG","FREQ_VSMALL_VEG","FREQ_VTALL_VEG","IMP_FLOATING_AQ","IMP_HMED_VEG","IMP_LIANAS"
                          ,"IMP_MED_VEG","IMP_SMALL_VEG","IMP_SUBMERGED_AQ","IMP_TALL_VEG","IMP_VSMALL_VEG","IMP_VTALL_VEG"
                          ,"XCOV_FLOATING_AQ","XCOV_HMED_VEG","XCOV_LIANAS","XCOV_MED_VEG","XCOV_SMALL_VEG","XCOV_SUBMERGED_AQ"
                          ,"XCOV_TALL_VEG","XCOV_VSMALL_VEG","XCOV_VTALL_VEG","XRCOV_FLOATING_AQ","XRCOV_HMED_VEG","XRCOV_LIANAS"
                          ,"XRCOV_MED_VEG","XRCOV_SMALL_VEG","XRCOV_SUBMERGED_AQ","XRCOV_TALL_VEG","XRCOV_VSMALL_VEG","XRCOV_VTALL_VEG"
                          ,"H_VASC_STRATA","J_VASC_STRATA","D_VASC_STRATA")

  outdf3 <- reshape2::dcast(outdf2,SAMPID~METRIC,value.var='RESULT')

  outdf4 <- subset(merge(outdf3, empty_vtype, all=TRUE),!is.na(SAMPID))

  outdf5 <- merge(allUIDs,outdf4,by='SAMPID',all.x=T)

  outdf6 <- reshape2::melt(outdf5,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')
  outdf6 <- mutate(outdf6,METRIC=as.character(METRIC)
                   ,RESULT=ifelse(METRIC %nin% c('D_VASC_STRATA','H_VASC_STRATA') & is.na(RESULT),0,RESULT))

  print("Done with vascular strata metrics")

  # Now combine vtypeMet and vhet8 into a single data frame and widen it
  vtOut <- rbind(outdf6,vstratMet) %>% plyr::rename(c('METRIC'='PARAMETER'))
  vtOut.1 <- merge(samples,vtOut,by='SAMPID') %>%
    dplyr::select(-SAMPID)
  return(vtOut.1)
}



#### NON-VASCULAR STRATA
#' @export
#' 
#' @title Calculate non-vascular metrics
#' 
#' @description This function calculates non-vascular plant
#' metrics based on vegetative and ground cover data. This
#' function is called by calcVtype_GcovMets().
#' 
#' @param dataIn A data frame containing the following variables:
#' \itemize{
#'  \item sampID: variable(s) identified in \emph{sampID} argument
#'
#'  \item PARAMETER: specific measurement type
#'
#'  \item RESULT: measured value
#'  }
#'  The following parameters are used in
#'  calculating tree metrics: 'PEAT_MOSS', 'BRYOPHYTES', 'LICHENS', 
#'  'ARBOREAL', 'ALGAE', 'MACROALGAE'. Additional parameters or
#'  variables are ignored.
#' @param nPlot A data frame with the 
#' number of plots sampled associated with each sample
#'  with sampID variables and NPLOTS.
#' @param sampID  A character vector containing the name(s) of
#'  variable(s) necessary to identify unique samples, 'UID'
#'  by default.
#' 
#' @details If any of the parameters are missing, they are assumed to
#'  be zeros (if numeric), and metric values associated with any metrics
#'  that cannot be calculated due to missing parameters are set to a
#'  standardized value.
#' 
#' @return Data frame containing \emph{sampID} variables, PARAMETER, RESULT,
#'   where values of PARAMETER consist of the metric name concatenated with
#'   trait value, with valid values of: FREQ_ALGAE,  FREQ_ARBOREAL, 
#'   FREQ_BRYOPHYTES, FREQ_LICHENS, FREQ_MACROALGAE, IMP_ALGAE, IMP_ARBOREAL,
#'   IMP_BRYOPHYTES, IMP_LICHENS, IMP_MACROALGAE, XCOV_ALGAE, XCOV_ARBOREAL,
#'   XCOV_BRYOPHYTES, XCOV_LICHENS, XCOV_MACROALGAE, N_PEAT_MOSS_DOM,
#'   FREQ_PEAT_MOSS_DOM. A list of metric descriptions is provided in the
#'   document named VegTypes_GrdCover_Metric_Descriptions.pdf included in the 
#'   help directory for the package.
#' 
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#'  head(Vtype_GrCovEx)
#'  # Create data frame with number of plots sampled for each sampID
#'  nplots <- data.frame(UID=seq(1:10),NPLOTS=rep(5,10),stringsAsFactors=F)
#'  # alternative approach to creating this data frame
#'  nplots <- plyr::ddply(Vtype_GrCovEx,c('UID'),dplyr::summarise,
#'  NPLOTS=length(unique(PLOT)))
#'
#'  nvEx <- calcNonvascMets(Vtype_GrCovEx,nplots)
#'
#'  head(nvEx)
#'  unique(nvEx$PARAMETER)
calcNonvascMets <- function(dataIn,nPlot,sampID='UID'){
  ## Now merge back with input df
  dataIn1 <- merge(dataIn,nPlot,by=sampID)

  # Create vector of all samples in dataset
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID,dataIn1[,sampID[i]],sep='.')
  }
  samples <- unique(subset(dataIn1,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID),stringsAsFactors=F)
  
  nvstrat <- subset(dataIn1,PARAMETER %in% c('PEAT_MOSS','BRYOPHYTES','LICHENS','ARBOREAL','ALGAE','MACROALGAE') &
                      !is.na(RESULT) & RESULT!='0')

  # Now use this data frame to calculate metrics
  ## First calculate metrics for non-vascular veg types, excluding peat moss
  indf1 <- plyr::ddply(subset(nvstrat,PARAMETER!='PEAT_MOSS' & RESULT!=0),c('SAMPID','PARAMETER'),summarise
                       ,FREQ=(length(PLOT)/unique(NPLOTS))*100
                 ,XCOV=(sum(as.numeric(RESULT))/unique(NPLOTS)),IMP=(FREQ+XCOV)/2)

  outdf <- mutate(reshape2::melt(indf1,id.vars=c('SAMPID','PARAMETER'),value.name='RESULT'),PARAMETER=paste(variable,PARAMETER,sep='_'))

  outdf1 <- reshape2::melt(reshape2::dcast(outdf,SAMPID~PARAMETER,value.var='RESULT'),id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')
  outdf1$RESULT[is.na(outdf1$RESULT)] <- 0

  ## now count number and frequency of plots with peat moss dominant (i.e., PEAT_MOSS='Y')
  ## Need to account for situations where PEAT_MOSS not present
  if(nrow(subset(nvstrat,PARAMETER=='PEAT_MOSS' & RESULT!='N'))>0){
    indf2 <- plyr::ddply(subset(nvstrat,PARAMETER=='PEAT_MOSS' & RESULT!='N'),c('SAMPID'),summarise,N_PEAT_MOSS_DOM=length(PLOT)
                   ,FREQ_PEAT_MOSS_DOM=(N_PEAT_MOSS_DOM/unique(NPLOTS))*100)

    outdf2 <- reshape2::melt(indf2,id.vars=c('SAMPID'),variable.name='PARAMETER',value.name='RESULT')
    outdf2a <- reshape2::melt(reshape2::dcast(outdf2,SAMPID~PARAMETER,value.var='RESULT'),id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')
    outdf2a <- mutate(outdf2a,RESULT=ifelse(is.na(RESULT),0,RESULT))

    outdf3 <- rbind(outdf1,outdf2a)
  }else{
    outdf3 <- outdf1
  }

  outdf4 <- reshape2::dcast(outdf3,SAMPID~METRIC,value.var='RESULT')
  empty_nv <- data.frame(t(rep(NA,17)),stringsAsFactors=F)
  names(empty_nv) <- c("FREQ_ALGAE", "FREQ_ARBOREAL","FREQ_BRYOPHYTES","FREQ_LICHENS","FREQ_MACROALGAE","IMP_ALGAE"
                       ,"IMP_ARBOREAL","IMP_BRYOPHYTES","IMP_LICHENS","IMP_MACROALGAE","XCOV_ALGAE","XCOV_ARBOREAL","XCOV_BRYOPHYTES"
                       ,"XCOV_LICHENS","XCOV_MACROALGAE","N_PEAT_MOSS_DOM","FREQ_PEAT_MOSS_DOM")

  outdf5 <- subset(merge(outdf4, empty_nv, all=TRUE),!is.na(SAMPID))

  outdf6 <- merge(allUIDs,outdf5,by='SAMPID',all.x=T)

  outdf7 <- reshape2::melt(outdf6,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')
  outdf7 <- mutate(outdf7,METRIC=as.character(METRIC),RESULT=ifelse(is.na(RESULT),0,RESULT))
  print("Done with non-vascular strata metrics")

  nvOut <- merge(samples,outdf7,by='SAMPID') %>%
    plyr::rename(c('METRIC'='PARAMETER')) %>%
    dplyr::select(-SAMPID)
  return(nvOut)
}


#' @export
#' 
#' @title Calculate water cover metrics
#' 
#' @description This function calculates metrics water
#' cover and depth based on ground cover data. This function is
#' called by calcVtype_GcovMets().
#' 
#' @param dataIn A data frame containing the following variables:
#' \itemize{
#' \item sampID: Variable(s) identified in \emph{sampID} argument
#' 
#' \item PLOT: Sample plot from which data
#' were collected
#'
#' \item PARAMETER: specific measurement type
#'
#' \item RESULT: measured value
#' }
#' The following parameters are used in
#' calculating tree metrics: 'TIME', 'MINIMUM_DEPTH', 
#' 'MAXIMUM_DEPTH', 'PREDOMINANT_DEPTH', 'TOTAL_WATER', 
#' 'WATER_NOVEG', 'WATER_AQVEG', 'WATER_EMERGVEG'.
#' Additional parameters or variables are ignored.
#' @param nPlot A data frame with the 
#' number of plots sampled associated with each
#' sample with \emph{sampID} variables and NPLOTS
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples,
#' 'UID' by default.
#' 
#' @details If any of the parameters are missing, they are assumed
#' to be zeros (if numeric), and metric values associated with any
#' metrics that cannot be calculated due to missing parameters are set
#' to a standardized value.
#' 
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. Data frame contains
#'   \emph{sampID} variables, PARAMETER, RESULT, where values of PARAMETER 
#'   consist of the metric name concatenated with trait value, with valid values
#'   of: MIN_H2O_DEPTH, MAX_H2O_DEPTH, XH2O_DEPTH_AA, MIN_COV_H2O, MAX_COV_H2O,
#'   FREQ_H2O, FREQ_H2O_AQVEG, FREQ_H2O_EMERGVEG, FREQ_H2O_NOVEG, XCOV_H2O, 
#'   XCOV_H2O_AQVEG, XCOV_H2O_EMERGVEG, XCOV_H2O_NOVEG, IMP_H2O, IMP_H2O_AQVEG,
#'   IMP_H2O_EMERGVEG, IMP_H2O_NOVEG, XH2O_DEPTH. A list of metric descriptions
#'   is provided in the document named VegTypes_GrdCover_Metric_Descriptions.pdf
#'   included in the help directory for the package.
#' 
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(Vtype_GrCovEx)
#' # Create data frame with number of plots sampled for each sampID
#' nplots <- data.frame(UID=seq(1:10),NPLOTS=rep(5,10),stringsAsFactors=F)
#' # alternative approach to creating this data frame
#' nplots <- plyr::ddply(Vtype_GrCovEx,c('UID'),dplyr::summarise,
#' NPLOTS=length(unique(PLOT)))
#'
#' wcovEx <- calcWcovMets(Vtype_GrCovEx,nplots)
#'
#' head(wcovEx)
#' unique(wcovEx$METRIC)
calcWcovMets <- function(dataIn,nPlot,sampID='UID'){
  ## Now merge back with input df
  dataIn1 <- merge(dataIn,nPlot,by=sampID)

  # Create vector of all samples in dataset
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID,dataIn1[,sampID[i]],sep='.')
  }
  samples <- unique(subset(dataIn1,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID),stringsAsFactors=F)
  
  ####### WATER COVER AND DEPTH
  wdep <- subset(dataIn1,PARAMETER %in% c('TIME','MINIMUM_DEPTH','MAXIMUM_DEPTH','PREDOMINANT_DEPTH','TOTAL_WATER'
                                          ,'WATER_NOVEG','WATER_AQVEG','WATER_EMERGVEG'))
  wdep1 <- reshape2::dcast(mutate(wdep,RESULT=as.numeric(RESULT)),SAMPID+PLOT+NPLOTS~PARAMETER,value.var='RESULT')
  wdep2 <- subset(wdep,RESULT %nin% c("") & !is.na(RESULT) & PARAMETER %in% c('TOTAL_WATER','WATER_NOVEG','WATER_AQVEG','WATER_EMERGVEG'))
  wdep2a <- reshape2::melt(reshape2::dcast(wdep2,SAMPID+PLOT+NPLOTS~PARAMETER,value.var='RESULT')
                           ,id.vars=c('SAMPID','PLOT','NPLOTS'),variable.name='PARAMETER',value.name='RESULT')
  wdep2a <- mutate(wdep2a,RESULT=ifelse(is.na(RESULT),0,as.numeric(RESULT)))
  
  wat1 <- plyr::ddply(wdep1,c('SAMPID'),summarise
                      ,MIN_H2O_DEPTH=ifelse(any(!is.na(MINIMUM_DEPTH)),min(MINIMUM_DEPTH,na.rm=TRUE),NA)
                      ,MAX_H2O_DEPTH=ifelse(any(!is.na(MAXIMUM_DEPTH)),max(MAXIMUM_DEPTH,na.rm=TRUE),NA)
                      ,XH2O_DEPTH_AA=ifelse(any(!is.na(PREDOMINANT_DEPTH))
                                                ,round(sum(PREDOMINANT_DEPTH,na.rm=TRUE)/unique(NPLOTS),2),NA)
                      ,MIN_COV_H2O=ifelse(any(!is.na(TOTAL_WATER)),min(TOTAL_WATER,na.rm=TRUE),NA)
                      ,MAX_COV_H2O=ifelse(any(!is.na(TOTAL_WATER)),max(TOTAL_WATER,na.rm=TRUE),NA)
  )
  wat1a <- reshape2::melt(wat1,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')

  ## Fix Inf values to 0s
  wat1a <- mutate(wat1a,RESULT=ifelse(RESULT %in% c('Inf','-Inf'),0,ifelse(RESULT=='Inf_-Inf','',RESULT)),METRIC=as.character(METRIC))

  wat2 <- plyr::ddply(subset(wdep2a,RESULT>0),c('SAMPID','PARAMETER'),summarise,FREQ_H2O=round((length(as.numeric(RESULT))/unique(NPLOTS))*100,2)
                ,XCOV_H2O=round(sum(as.numeric(RESULT))/unique(NPLOTS),2)
                ,IMP_H2O=round((FREQ_H2O+XCOV_H2O)/2,2))
  wat2a <- reshape2::melt(wat2,id.vars=c('SAMPID','PARAMETER'),variable.name='METRIC',value.name='RESULT')
  wat2a <- mutate(wat2a,METRIC=as.character(ifelse(PARAMETER=='TOTAL_WATER',as.character(METRIC),paste(METRIC,substring(PARAMETER,7),sep='_')))
                  ,PARAMETER=NULL)

  ## Now pull only PREDOMINANT_DEPTH that is not missing and >0 to calculate XH2O_DEPTH
  wat3 <- plyr::ddply(subset(wdep1,PREDOMINANT_DEPTH>0 & !is.na(PREDOMINANT_DEPTH)),c('SAMPID'),summarise,METRIC='XH2O_DEPTH'
                ,RESULT=mean(PREDOMINANT_DEPTH,na.rm=TRUE))

  watMet <- rbind(wat1a,wat2a,wat3)
  watMet1 <- reshape2::dcast(watMet,SAMPID~METRIC,value.var='RESULT')

  empty_wat <- data.frame(t(rep(NA,18)),stringsAsFactors=F)
  names(empty_wat) <- c("MIN_H2O_DEPTH","MAX_H2O_DEPTH","XH2O_DEPTH_AA","MIN_COV_H2O","MAX_COV_H2O","FREQ_H2O"
                        ,"FREQ_H2O_AQVEG","FREQ_H2O_EMERGVEG","FREQ_H2O_NOVEG","XCOV_H2O","XCOV_H2O_AQVEG","XCOV_H2O_EMERGVEG"
                        ,"XCOV_H2O_NOVEG","IMP_H2O","IMP_H2O_AQVEG","IMP_H2O_EMERGVEG","IMP_H2O_NOVEG","XH2O_DEPTH")

  watMet2 <- subset(merge(watMet1, empty_wat, all=TRUE),!is.na(SAMPID))

  watMet3 <- merge(allUIDs,watMet2,by='SAMPID',all.x=T)

  watMet4 <- reshape2::melt(watMet3,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')
  watMet4 <- mutate(watMet4,METRIC=as.character(METRIC),RESULT=ifelse(is.na(RESULT),'0',RESULT))

  print("Done with water depth metrics")
  watMetOut <- merge(samples,watMet4,by='SAMPID') %>%
    plyr::rename(c('METRIC'='PARAMETER')) %>%
    dplyr::select(-SAMPID)
  
  return(watMetOut)
}



#' @export
#' 
#' @title Calculate bare ground and litter metrics
#' 
#' @description Calculate bare ground and litter metrics based on
#' ground cover data. This function is called by calcVtype_GcovMets().
#' 
#' @param dataIn A data frame containing the following variables:
#' \itemize{
#' \item sampID: Variable(s) identified in \emph{sampID} argument
#' 
#' \item PLOT (Sample plot from which data were collected)
#'
#' \item PARAMETER (specific measurement type)
#'
#' \item RESULT (measured value)
#' }
#' The following parameters are used in
#' calculating tree metrics: 'LITTER_THATCH', 'LITTER_FORB', 
#' 'LITTER_CONIFER', 'LITTER_DECID', 'LITTER_BROADLEAF', 
#' 'LITTER_DEPTH_SW', 'LITTER_DEPTH_NE', 'TOTAL_LITTER', 'WD_FINE', 
#' 'WD_COARSE', 'EXPOSED_SOIL', 'EXPOSED_GRAVEL', 'EXPOSED_ROCK'.
#' Additional parameters or variables are ignored.
#' @param nPlot A data frame with the 
#' number of plots sampled associated with each sample,
#' including \emph{sampID} variables and NPLOTS.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' 
#' @details If any of the parameters are missing, they are assumed to be
#' zeros (if numeric), and metric values associated with any metrics that
#' cannot be calculated due to missing parameters are set to a
#' standardized value.
#' 
#' @return Either a character string containing an error message when metric
#'   calculation is not successful, or a data frame. Data frame contains
#'   \emph{sampID} variables, PARAMETER, RESULT, where values of PARAMETER 
#'   consist of the metric name concatenated with trait value, with valid values
#'   of: FREQ_BAREGD, FREQ_EXPOSED_GRAVEL, FREQ_EXPOSED_ROCK, FREQ_EXPOSED_SOIL,
#'   FREQ_LITTER, FREQ_WD_COARSE, FREQ_WD_FINE, IMP_BAREGD, IMP_EXPOSED_GRAVEL,
#'   IMP_EXPOSED_ROCK, IMP_EXPOSED_SOIL, IMP_LITTER, IMP_WD_COARSE, IMP_WD_FINE,
#'   XCOV_BAREGD, XCOV_EXPOSED_GRAVEL, XCOV_EXPOSED_ROCK, XCOV_EXPOSED_SOIL,
#'   XCOV_LITTER, XCOV_WD_COARSE, XCOV_WD_FINE. A list of metric descriptions is
#'   provided in the document named VegTypes_GrdCover_Metric_Descriptions.pdf 
#'   included in the help directory for the package.
#'   
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(Vtype_GrCovEx)
#' # Create data frame with number of plots sampled for each UID
#' nplots <- data.frame(UID=seq(1:10),NPLOTS=rep(5,10),stringsAsFactors=F)
#' # alternative approach to creating this data frame
#' nplots <- plyr::ddply(Vtype_GrCovEx,c('UID'),dplyr::summarise,
#' NPLOTS=length(unique(PLOT)))
#'
#' bgEx <- calcBareGround_LitterMets(Vtype_GrCovEx,nplots)
#'
#' head(bgEx)
#' unique(bgEx$PARAMETER)

calcBareGround_LitterMets <- function(dataIn,nPlot,sampID='UID'){
  ## Now merge back with input df
  dataIn1 <- merge(dataIn,nPlot,by=sampID)

  # Create vector of all samples in dataset
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID,dataIn1[,sampID[i]],sep='.')
  }
  samples <- unique(subset(dataIn1,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID),stringsAsFactors=F)
  
  ####### LITTER TYPES
  # Need to calculate the number of quadrats sampled using the NE and SW parameters
  litter.sub <- subset(dataIn1,PARAMETER %in% c('LITTER_DEPTH_NE','LITTER_DEPTH_SW'))
  numQuads <- plyr::ddply(litter.sub,c('SAMPID'),summarise,NQUADS=length(RESULT))

  litter1 <- subset(dataIn1,PARAMETER %in% c('LITTER_THATCH','LITTER_FORB','LITTER_CONIFER','LITTER_DECID','LITTER_BROADLEAF'
                                             ,'LITTER_DEPTH_SW','LITTER_DEPTH_NE','TOTAL_LITTER','WD_FINE','WD_COARSE','EXPOSED_SOIL','EXPOSED_GRAVEL'
                                             ,'EXPOSED_ROCK'))

  litter2 <- merge(litter1,numQuads,by='SAMPID')

  littype <- subset(litter2,PARAMETER %in% c('LITTER_THATCH','LITTER_FORB','LITTER_CONIFER','LITTER_DECID','LITTER_BROADLEAF','LITTER_NONE') & RESULT!=0)
  rr1 <- plyr::ddply(littype,c('SAMPID'),mutate,N_LITTER_TYPE=length(unique(PARAMETER)))
  rr2 <- plyr::ddply(rr1,c('SAMPID','PARAMETER','N_LITTER_TYPE'),summarise,NUM=length(PLOT))
  rr3 <- plyr::ddply(rr2,c('SAMPID'),mutate,MAXN=max(NUM))
  rr4 <- subset(rr3,MAXN==NUM)
  rr5 <- plyr::ddply(rr4,c('SAMPID','N_LITTER_TYPE'),summarise,LITTER_TYPE=paste(PARAMETER,collapse='_'))
  rr5 <- mutate(rr5,LITTER_TYPE=gsub('LITTER_','',LITTER_TYPE))
  ## to determine median depth, we must account for any quadrats without depth recorded
  litdep <- subset(litter2,PARAMETER %in% c('LITTER_DEPTH_SW','LITTER_DEPTH_NE'))

  ss1 <- plyr::ddply(litdep,c('SAMPID'),summarise,NSAMP=length(RESULT),XDEPTH_LITTER=round(sum(as.numeric(RESULT))/unique(NQUADS),2))

  litdep1 <- plyr::ddply(litdep,c('SAMPID','NPLOTS','NQUADS'),mutate,NSAMP=length(RESULT),toAdd=NQUADS-NSAMP)

  tt <- plyr::ddply(litdep1,c('SAMPID'),summarise,MEDDEPTH_LITTER=median(as.numeric(RESULT)))

  loutdf <- rbind(reshape2::melt(rr5,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')
                  ,reshape2::melt(ss1,id.vars='SAMPID'
                  ,measure.vars='XDEPTH_LITTER',variable.name='METRIC',value.name='RESULT')
                  ,reshape2::melt(tt,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT'))
  loutdf <- mutate(loutdf,METRIC=as.character(METRIC))
  loutdf1 <- reshape2::dcast(loutdf,SAMPID~METRIC,value.var='RESULT')

  empty_lit <- data.frame(t(rep(NA,4)),stringsAsFactors=FALSE)
  names(empty_lit) <- c('N_LITTER_TYPE','LITTER_TYPE','XDEPTH_LITTER','MEDDEPTH_LITTER')

  loutdf2 <- subset(merge(loutdf1, empty_lit, all=TRUE),!is.na(SAMPID))

  loutdf3 <- merge(allUIDs,loutdf2,by='SAMPID',all.x=T)

  loutdf4 <- reshape2::melt(loutdf3,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')
  loutdf4 <- mutate(loutdf4,METRIC=as.character(METRIC))

  print("Done with litter metrics")
  litterOut <- reshape2::dcast(loutdf4,SAMPID~METRIC,value.var='RESULT')

  ################# BARE GROUND
  bgrd <- subset(dataIn1,PARAMETER %in% c('EXPOSED_SOIL','EXPOSED_GRAVEL','EXPOSED_ROCK','WD_FINE','WD_COARSE','TOTAL_LITTER') & RESULT!=0
                 ,select=c('SAMPID','PLOT','PARAMETER','RESULT','NPLOTS'))
  ## Need to create values for new PARAMETER='BAREGD' based on occurrence of either EXPOSED_SOIL, EXPOSED_GRAVEL, or EXPOSED_ROCK at site
  bgrd1 <- plyr::ddply(subset(bgrd,PARAMETER %in% c('EXPOSED_SOIL','EXPOSED_GRAVEL','EXPOSED_ROCK')
                              ,select=c('SAMPID','PLOT','NPLOTS','RESULT'))
                 ,c('SAMPID','PLOT','NPLOTS'),summarise,PARAMETER='BAREGD',RESULT=as.character(sum(as.numeric(RESULT))))
  bgrdIn <- rbind(bgrd,bgrd1)
  ## Need to fill in zeros if plot sampled and variable is zero
  bgrdIn1 <- reshape2::melt(reshape2::dcast(bgrdIn,SAMPID+PLOT+NPLOTS~PARAMETER,value.var='RESULT'),id.vars=c('SAMPID','PLOT','NPLOTS'),variable.name='PARAMETER'
                  ,value.name='RESULT')
  bgrdIn1 <- mutate(bgrdIn1,RESULT=ifelse(is.na(RESULT),0,as.numeric(RESULT)))

  ## First calculate metrics for non-vascular veg types, excluding peat moss
  bgindf1 <- plyr::ddply(subset(bgrdIn1,RESULT!=0),c('SAMPID','PARAMETER'),summarise,FREQ=round((length(PLOT)/unique(NPLOTS))*100,2)
                   ,XCOV=round((sum(as.numeric(RESULT))/unique(NPLOTS)),2),IMP=round((FREQ+XCOV)/2,2))

  bgoutdf <- mutate(reshape2::melt(bgindf1,id.vars=c('SAMPID','PARAMETER'),value.name='RESULT'),PARAMETER=ifelse(PARAMETER=='TOTAL_LITTER'
                                                                                                    ,paste(variable,'LITTER',sep='_'),paste(variable,PARAMETER,sep='_')))

  bgoutdf1 <- reshape2::melt(reshape2::dcast(bgoutdf,SAMPID~PARAMETER,value.var='RESULT')
                             ,id.vars='SAMPID',variable.name='METRIC',value.name='RESULT')
  bgoutdf1 <- mutate(bgoutdf1,RESULT=ifelse(is.na(RESULT),0,RESULT),METRIC=as.character(METRIC))
  bgoutdf2 <- reshape2::dcast(bgoutdf1,SAMPID~METRIC,value.var='RESULT')

  empty_bg <- data.frame(t(rep(NA,21)),stringsAsFactors=FALSE)
  names(empty_bg) <- c("FREQ_BAREGD","FREQ_EXPOSED_GRAVEL","FREQ_EXPOSED_ROCK","FREQ_EXPOSED_SOIL","FREQ_LITTER","FREQ_WD_COARSE"
                       ,"FREQ_WD_FINE","IMP_BAREGD","IMP_EXPOSED_GRAVEL","IMP_EXPOSED_ROCK","IMP_EXPOSED_SOIL","IMP_LITTER"
                       ,"IMP_WD_COARSE","IMP_WD_FINE","XCOV_BAREGD","XCOV_EXPOSED_GRAVEL","XCOV_EXPOSED_ROCK","XCOV_EXPOSED_SOIL"
                       ,"XCOV_LITTER","XCOV_WD_COARSE","XCOV_WD_FINE")

  bgoutdf3 <- subset(merge(bgoutdf2, empty_bg, all=TRUE),!is.na(SAMPID))

  bgoutdf4 <- merge(allUIDs,bgoutdf3,by='SAMPID',all.x=T)

  bgoutdf5 <- reshape2::melt(bgoutdf4,id.vars=c('SAMPID'),variable.name='METRIC',value.name='RESULT')
  bgoutdf5 <- mutate(bgoutdf5,METRIC=as.character(METRIC),RESULT=ifelse(is.na(RESULT),0,RESULT))

  print("Done with bare ground metrics")
  bgrdOut <- reshape2::dcast(bgoutdf5,SAMPID~METRIC,value.var='RESULT')

  # Now combine bare ground and litter metrics to perform a final check
  combOut <- merge(litterOut,bgrdOut,by='SAMPID',all=T)
  combOut.1 <- reshape2::melt(combOut,id.vars=c('SAMPID'))
    combOut.1 <- mutate(combOut.1,value=ifelse(is.na(value) & variable!='LITTER_TYPE',0,value))
  # Need to determine whether missing LITTER_TYPE value should be ABSENT or NONE
  combOut.wide <- reshape2::dcast(combOut.1,SAMPID~variable,value.var='value')
  combOut.wide <- mutate(combOut.wide,LITTER_TYPE=ifelse(is.na(LITTER_TYPE) & XDEPTH_LITTER==0 & XCOV_LITTER==0 & FREQ_LITTER==0
                                               ,'ABSENT',ifelse(is.na(LITTER_TYPE),'NONE',LITTER_TYPE)))
  combOut.long <- reshape2::melt(combOut.wide,id.vars='SAMPID',variable.name='PARAMETER',value.name='RESULT')
  
  finOut <- merge(samples,combOut.long,by='SAMPID') %>%
    dplyr::select(-SAMPID)
  return(finOut)
}

