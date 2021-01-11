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
#' wetland type metrics: 'SANDT_CLASS' and 'PAL_FARMED', or 
#' 'WETLAND_TYPE' (2016).
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
  
  vhet <- reshape(dataIn1, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'wide',
                  timevar = 'PARAMETER', v.names = 'RESULT')
  names(vhet) <- gsub("RESULT\\.", "", names(vhet))
  

  ## First calculate heterogeneity metrics ## Account for datasets without PF parameter so that wide form will not have it
  if('PAL_FARMED' %in% names(vhet)==TRUE){
    vhet$SANDT_CLASS <- with(vhet, ifelse(is.na(PAL_FARMED),SANDT_CLASS,'PF'))
  }
  vhet <- subset(vhet, !is.na(SANDT_CLASS)) 
  
  vhet1.n <- aggregate(x = list(N_SANDT = vhet$SANDT_CLASS), 
                       by = vhet[c('SAMPID')], FUN = function(x){length(unique(x))})
  vhet1 <- merge(vhet, vhet1.n, by='SAMPID')

  vhet2 <- aggregate(x = list(FREQ = vhet1$PLOT), 
                     by = vhet1[c('SAMPID','NPLOTS','SANDT_CLASS','N_SANDT')],
                     FUN = length)
  
  vhet2$FREQ <- with(vhet2, FREQ/NPLOTS)
  vhet3.maxf <- aggregate(x = list(MAXF = vhet2$FREQ), by = vhet2[c('SAMPID')], FUN = max)
  vhet3 <- merge(vhet2, vhet3.maxf, by='SAMPID')
  vhet4 <- subset(vhet3, MAXF==FREQ)
  # Must aggregate in case there are two dominant classes
  vhet4a <- aggregate(x = list(DOM_SANDT = vhet4$SANDT_CLASS), 
                      by = vhet4[c('SAMPID','N_SANDT')],
                      FUN = function(x){paste(x, collapse = '-')})

  vhet5.d <- aggregate(x = list(D_SANDT = vhet2$FREQ), by = vhet2[c('SAMPID','NPLOTS')],
                       FUN = function(x){round(1 - sum(x*x), 4)})
  
  vhet5.h <- aggregate(x = list(H_SANDT = vhet2$FREQ), by = vhet2[c('SAMPID','NPLOTS')],
                       FUN = function(x){round(-1*sum(x*log(x)),4)})
  
  vhet5.jcalc <- aggregate(x = list(uniqST = vhet2$N_SANDT), by = vhet2[c('SAMPID','NPLOTS')],
                           FUN = function(x){unique(x)}) 
  
  vhet5 <- merge(vhet5.d, vhet5.h, by=c('SAMPID','NPLOTS'))
  vhet5 <- merge(vhet5, vhet5.jcalc, by=c('SAMPID','NPLOTS'))
  vhet5$J_SANDT <- with(vhet5, round(ifelse(H_SANDT!=0 & uniqST!=1,H_SANDT/log(uniqST),0),4))
  
  vhet6 <- merge(vhet4a,vhet5,by='SAMPID')
  # create an empty data frame with all of the metric names as variables to ensure all are included in output
  empty_vhet <- data.frame(t(rep(NA,5)),stringsAsFactors=FALSE)
  names(empty_vhet) <- c('N_SANDT','DOM_SANDT','D_SANDT','H_SANDT','J_SANDT')

  vhet7 <- subset(merge(vhet6, empty_vhet, all=TRUE),!is.na(SAMPID))

  vhet8 <- merge(allUIDs,vhet7,by='SAMPID',all.x=T)

  vhet9 <- reshape(vhet8, idvar = c('SAMPID'), direction = 'long',
                   times = c('D_SANDT','H_SANDT','J_SANDT','DOM_SANDT','N_SANDT'),
                   varying = c('D_SANDT','H_SANDT','J_SANDT','DOM_SANDT','N_SANDT'),
                   timevar = 'METRIC', v.names = 'RESULT')
  vhet9$METRIC <- with(vhet9, as.character(METRIC))
  vhet9$RESULT <- with(vhet9, ifelse(METRIC=='DOM_SANDT' & is.na(RESULT),'MISSING'
                                     ,ifelse(is.na(RESULT),0,RESULT)))

  print("Done with veg heterogeneity metrics")
  vhetOut <- merge(samples, vhet9, by='SAMPID') 
  vhetOut$PARAMETER <- vhetOut$METRIC
  
  vhetOut <- vhetOut[,c(sampID, 'PARAMETER', 'RESULT')]
  
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
#' calculating vascular strata metrics: 'SUBMERGED_AQ', 'FLOATING_AQ', 
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
#' nplots <- aggregate(x = list(NPLOTS = Vtype_GrCovEx$PLOT), 
#' by = Vtype_GrCovEx[c('UID')],
#' FUN = function(x){length(unique(x))})
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

  vstrat <- subset(dataIn1,PARAMETER %in% c('SUBMERGED_AQ','FLOATING_AQ','LIANAS','VTALL_VEG',
                                            'TALL_VEG','HMED_VEG','MED_VEG',
                                            'SMALL_VEG','VSMALL_VEG') & !is.na(RESULT) & RESULT!='0')
  
  vstrat.sum <- aggregate(x = list(XTOTCOV_VASC_STRATA = vstrat$RESULT), 
                          by = vstrat[c('SAMPID','NPLOTS')],
                          FUN = function(x){sum(as.numeric(x))})
  vstrat.sum$XTOTCOV_VASC_STRATA <- with(vstrat.sum, XTOTCOV_VASC_STRATA/NPLOTS)
  
  vstrat.cnt <- aggregate(x = list(N_VASC_STRATA = vstrat$PARAMETER, PLOTSAMP = vstrat$PLOT),
                          by = vstrat[c('SAMPID')], FUN = function(x){length(unique(x))})
  
  vstrat <- merge(vstrat, vstrat.sum, by=c('SAMPID','NPLOTS'))
  vstrat <- merge(vstrat, vstrat.cnt, by='SAMPID')
  
  vstratPlot.n <- aggregate(x = list(N_VSTRATA = vstrat$PARAMETER), 
                            by = vstrat[c('SAMPID','PLOT','N_VASC_STRATA','NPLOTS',
                                          'XTOTCOV_VASC_STRATA','PLOTSAMP')],
                            FUN = function(x){length(unique(x))})
  
  vstratPlot.sum <- aggregate(x = list(SUM_PLOT = vstrat$RESULT), 
                              by = vstrat[c('SAMPID','PLOT','N_VASC_STRATA','NPLOTS',
                                            'XTOTCOV_VASC_STRATA','PLOTSAMP')],
                              FUN = function(x){sum(as.numeric(x))})
  
  vstratPlot <- merge(vstratPlot.n, vstratPlot.sum, by=c('SAMPID','PLOT','N_VASC_STRATA',
                                                         'NPLOTS','XTOTCOV_VASC_STRATA','PLOTSAMP'))

  vstrat1.sum <- aggregate(x = list(XN_VASC_STRATA = vstratPlot$N_VSTRATA),
                           by = vstratPlot[c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP')],
                           FUN = sum)
  
  vstrat1.min <- aggregate(x = list(minN = vstratPlot$N_VSTRATA),
                           by = vstratPlot[c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP')],
                           FUN = min)
  
  vstrat1.max <- aggregate(x = list(maxN = vstratPlot$N_VSTRATA),
                           by = vstratPlot[c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP')],
                           FUN = max)
  
  vstrat1 <- merge(vstrat1.sum, vstrat1.min, 
                   by=c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP'))
  vstrat1 <- merge(vstrat1, vstrat1.max, c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP'))
  
  vstrat1$XN_VASC_STRATA <- with(vstrat1, XN_VASC_STRATA/NPLOTS)
  vstrat1$RG_VASC_STRATA <- with(vstrat1, ifelse(PLOTSAMP==NPLOTS, maxN - minN, maxN - 0))
  
  vstrat1$minN <- NULL
  vstrat1$maxN <- NULL
  vstrat1$PLOTSAMP <- NULL
  vstrat1$NPLOTS <- NULL
  
  empty_vstrat <- data.frame(t(rep(NA,4)),stringsAsFactors=FALSE)
  names(empty_vstrat) <- c('N_VASC_STRATA','XTOTCOV_VASC_STRATA','XN_VASC_STRATA','RG_VASC_STRATA')

  vstrat2 <- subset(merge(vstrat1, empty_vstrat, all=TRUE),!is.na(SAMPID))

  vstrat3 <- merge(allUIDs,vstrat2,by='SAMPID',all.x=T)

  varNames <- names(vstrat3)[!names(vstrat3) %in% c('SAMPID')]
  vstratMet <- reshape(vstrat3, idvar = 'SAMPID', direction = 'long',
                       varying = varNames, times = varNames,
                       timevar = 'METRIC', v.names = 'RESULT')
  
  vstratMet$METRIC <- with(vstratMet, as.character(METRIC))
  vstratMet$RESULT <- with(vstratMet, ifelse(is.na(RESULT),0,RESULT))

  ## Calculate frequency, mean cover, relative mean cover, and relative importance by vascular stratum
  vstrat.pos <- subset(vstrat,RESULT!=0)
  
  indf1.length <- aggregate(x = list(FREQ = vstrat.pos$PLOT), 
                            by = vstrat.pos[c('SAMPID','PARAMETER','NPLOTS','XTOTCOV_VASC_STRATA')],
                            FUN = length)
  
  indf1.sum <- aggregate(x = list(XCOV = vstrat.pos$RESULT), 
                         by = vstrat.pos[c('SAMPID','PARAMETER','NPLOTS','XTOTCOV_VASC_STRATA')],
                         FUN = function(x){sum(as.numeric(x))})
  
  indf1 <- merge(indf1.length, indf1.sum, by = c('SAMPID','PARAMETER','NPLOTS','XTOTCOV_VASC_STRATA'), all=TRUE)
  indf1$FREQ <- with(indf1, round(FREQ/NPLOTS*100, 2))
  indf1$XCOV <- with(indf1, round(XCOV/NPLOTS, 2))
  indf1$XRCOV <- with(indf1, round(XCOV/XTOTCOV_VASC_STRATA*100, 2))
  indf1$IMP <- with(indf1, round((FREQ + XCOV)/2, 2))
  
  # Now drop variables
  indf1$XTOTCOV_VASC_STRATA <- NULL
  indf1$NPLOTS <- NULL
  
  outdf <- reshape(indf1, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                   varying = c('FREQ','XCOV','XRCOV','IMP'), times = c('FREQ','XCOV','XRCOV','IMP'),
                   timevar = 'variable', v.names = 'RESULT')
  
  outdf$PARAMETER <- with(outdf, paste(variable, PARAMETER, sep = '_'))
  outdf$variable <- NULL

  outdf.wide <- reshape(outdf, idvar = c('SAMPID'), direction = 'wide',
                        timevar = 'PARAMETER', v.names = 'RESULT')
  names(outdf.wide) <- gsub("RESULT\\.", "", names(outdf.wide))
  
  varNames <- names(outdf.wide)[!names(outdf.wide) %in% c('SAMPID')]
  outdf1 <- reshape(outdf.wide, idvar = c('SAMPID'), direction = 'long',
                    varying = varNames, times = varNames,
                    timevar = 'METRIC', v.names = 'RESULT')
  outdf1$RESULT[is.na(outdf1$RESULT)] <- 0

  ## Calculate diversity indices
  div1.h <- aggregate(x = list(H_VASC_STRATA = indf1$XRCOV), 
                      by = indf1[c('SAMPID')], FUN = function(x){round(-1*sum((x/100)*log(x/100)), 4)})
  
  div1.j <- aggregate(x = list(J_VASC_STRATA = indf1$PARAMETER), 
                      by = indf1[c('SAMPID')], FUN = length)
  
  div1.d <- aggregate(x = list(D_VASC_STRATA = indf1$XRCOV),
                      by = indf1[c('SAMPID')],
                      FUN = function(x){round(1 - sum((x/100)^2), 4)})
  
  div1 <- merge(div1.h, div1.j, by = 'SAMPID', all=TRUE)
  div1 <- merge(div1, div1.d, by='SAMPID', all=TRUE)
  
  div1$J_VASC_STRATA <- with(div1, round(H_VASC_STRATA/log(J_VASC_STRATA), 4))
  div1$J_VASC_STRATA <- with(div1, ifelse(!is.na(J_VASC_STRATA),J_VASC_STRATA,0))
  
  div1.long <- reshape(div1, idvar = 'SAMPID', direction = 'long',
                       varying = c('H_VASC_STRATA','J_VASC_STRATA','D_VASC_STRATA'), 
                       times = c('H_VASC_STRATA','J_VASC_STRATA','D_VASC_STRATA'),
                       timevar = 'METRIC', v.names = 'RESULT')

  outdf2 <- rbind(outdf1, div1.long)

  empty_vtype <- data.frame(t(rep(NA,39)),stringsAsFactors=FALSE)
  names(empty_vtype) <- c("FREQ_FLOATING_AQ","FREQ_HMED_VEG","FREQ_LIANAS","FREQ_MED_VEG",
                          "FREQ_SMALL_VEG","FREQ_SUBMERGED_AQ",
                          "FREQ_TALL_VEG","FREQ_VSMALL_VEG","FREQ_VTALL_VEG",
                          "IMP_FLOATING_AQ","IMP_HMED_VEG","IMP_LIANAS",
                          "IMP_MED_VEG","IMP_SMALL_VEG","IMP_SUBMERGED_AQ",
                          "IMP_TALL_VEG","IMP_VSMALL_VEG","IMP_VTALL_VEG",
                          "XCOV_FLOATING_AQ","XCOV_HMED_VEG","XCOV_LIANAS",
                          "XCOV_MED_VEG","XCOV_SMALL_VEG","XCOV_SUBMERGED_AQ",
                          "XCOV_TALL_VEG","XCOV_VSMALL_VEG","XCOV_VTALL_VEG",
                          "XRCOV_FLOATING_AQ","XRCOV_HMED_VEG","XRCOV_LIANAS",
                          "XRCOV_MED_VEG","XRCOV_SMALL_VEG","XRCOV_SUBMERGED_AQ",
                          "XRCOV_TALL_VEG","XRCOV_VSMALL_VEG","XRCOV_VTALL_VEG",
                          "H_VASC_STRATA","J_VASC_STRATA","D_VASC_STRATA")

  outdf3 <- reshape(outdf2, idvar = 'SAMPID', direction = 'wide',
                    timevar = 'METRIC', v.names = 'RESULT')
  
  names(outdf3) <- gsub("RESULT\\.", "", names(outdf3))

  outdf4 <- subset(merge(outdf3, empty_vtype, all=TRUE),!is.na(SAMPID))

  outdf5 <- merge(allUIDs,outdf4,by='SAMPID',all.x=T)

  varNames <- names(outdf5)[!names(outdf5) %in% c('SAMPID')]
  outdf6 <- reshape(outdf5, idvar = 'SAMPID', direction = 'long',
                    varying = varNames, times = varNames,
                    timevar = 'METRIC', v.names = 'RESULT')
  outdf6$METRIC <- with(outdf6, as.character(METRIC))
  outdf6$RESULT <- with(outdf6, ifelse(METRIC %nin% c('D_VASC_STRATA','H_VASC_STRATA') & is.na(RESULT),0,RESULT))

  print("Done with vascular strata metrics")

  # Now combine vtypeMet and vhet8 into a single data frame and widen it
  vtOut <- rbind(outdf6,vstratMet) 
  vtOut$PARAMETER <- vtOut$METRIC
  vtOut$METRIC <- NULL
  
  vtOut.1 <- merge(samples,vtOut,by='SAMPID') 
  vtOut.1 <- vtOut.1[, c(sampID, 'PARAMETER', 'RESULT')]
    
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
#'  calculating non-vascular metrics: 'PEAT_MOSS', 'BRYOPHYTES', 'LICHENS', 
#'  'ARBOREAL', 'ALGAE', 'MACROALGAE'. 
#'  Additional parameters or variables are ignored. PEAT_MOSS has values of
#'  'Y/N' for 2011 but only Y or missing for 2016.
#' @param nPlot A data frame with the 
#' number of plots sampled associated with each sample
#'  with sampID variables and NPLOTS.
#' @param sampID  A character vector containing the name(s) of
#'  variable(s) necessary to identify unique samples, 'UID'
#'  by default.
#' @param survyear A string with the survey year. The default is '2011'.
#' Starting in 2016, ARBOREAL began being measured as a categorical 
#' variable, ARBOREAL_ABUNDANCE.
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
#'  nplots <- aggregate(x = list(NPLOTS = Vtype_GrCovEx$PLOT), 
#'  by = Vtype_GrCovEx[c('UID')],
#'  FUN = function(x){length(unique(x))})
#'
#'  nvEx <- calcNonvascMets(Vtype_GrCovEx,nplots)
#'
#'  head(nvEx)
#'  unique(nvEx$PARAMETER)
calcNonvascMets <- function(dataIn,nPlot,sampID='UID', survyear='2011'){
  ## Now merge back with input df
  dataIn1 <- merge(dataIn, nPlot, by=sampID)

  # Create vector of all samples in dataset
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID,dataIn1[,sampID[i]],sep='.')
  }
  samples <- unique(subset(dataIn1,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID),stringsAsFactors=F)
  
  nvstrat <- subset(dataIn1,PARAMETER %in% c('PEAT_MOSS','BRYOPHYTES','LICHENS','ARBOREAL','ALGAE','MACROALGAE') &
                      !is.na(RESULT) & RESULT!='0')
  nvstrat$PARAMETER <- with(nvstrat, gsub('_ABUNDANCE', '', PARAMETER)) # This accounts for change in ARBOREAL to ARBOREAL_ABUNDANCE

  # Now use this data frame to calculate metrics
  ## First calculate metrics for non-vascular veg types, excluding peat moss
  nvstrat.sub <- subset(nvstrat,PARAMETER!='PEAT_MOSS' & RESULT!=0)
  indf1.sum <- aggregate(x = list(XCOV = nvstrat.sub$RESULT),
                         by = nvstrat.sub[c('SAMPID','PARAMETER','NPLOTS')],
                         FUN = function(x){sum(as.numeric(x))})
  
  indf1.length <- aggregate(x = list(FREQ = nvstrat.sub$PLOT), 
                            by = nvstrat.sub[c('SAMPID','PARAMETER','NPLOTS')],
                            FUN = length)
  
  indf1 <- merge(indf1.sum, indf1.length, by = c('SAMPID','PARAMETER','NPLOTS'), all=TRUE)
  
  indf1$FREQ <- with(indf1, FREQ/NPLOTS*100)
  indf1$XCOV <- with(indf1, XCOV/NPLOTS)
  indf1$IMP <- with(indf1, (FREQ + XCOV)/2)
  indf1$NPLOTS <- NULL

  outdf <- reshape(indf1, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                        varying = c('XCOV','FREQ','IMP'), times = c('XCOV','FREQ','IMP'),
                        timevar = 'variable', v.names = 'RESULT') 
  outdf$PARAMETER <- with(outdf, paste(variable, PARAMETER, sep = '_'))
  outdf$variable <- NULL
  
  outdf.wide <- reshape(outdf, idvar = 'SAMPID', direction = 'wide',
                        timevar = 'PARAMETER', v.names = 'RESULT')
  names(outdf.wide) <- gsub("RESULT\\.", "", names(outdf.wide))
  
  varNames <- names(outdf.wide)[!names(outdf.wide) %in% c('SAMPID')]
  outdf1 <- reshape(outdf.wide, idvar = "SAMPID", direction = "long",
                    varying = varNames, times = varNames,
                    timevar = 'METRIC', v.names = 'RESULT')
  
  outdf1$RESULT[is.na(outdf1$RESULT)] <- 0

  ## now count number and frequency of plots with peat moss dominant (i.e., PEAT_MOSS='Y')
  ## Need to account for situations where PEAT_MOSS not present
  if(nrow(subset(nvstrat,PARAMETER=='PEAT_MOSS' & RESULT!='N' & !is.na(RESULT)))>0){
    nvstrat.peat <- subset(nvstrat,PARAMETER=='PEAT_MOSS' & RESULT!='N' & !is.na(RESULT))

    indf2 <- aggregate(x = list(N_PEAT_MOSS_DOM = nvstrat.peat$PLOT),
                       by = nvstrat.peat[c('SAMPID','NPLOTS')],
                       FUN = length)
    
    indf2$FREQ_PEAT_MOSS_DOM <- with(indf2, N_PEAT_MOSS_DOM/NPLOTS*100)
    indf2$NPLOTS <- NULL
    
    outdf2 <- reshape(indf2, idvar = c('SAMPID'), direction = 'long',
                      varying = c('N_PEAT_MOSS_DOM','FREQ_PEAT_MOSS_DOM'),
                      times = c('N_PEAT_MOSS_DOM','FREQ_PEAT_MOSS_DOM'),
                      timevar = 'PARAMETER', v.names = 'RESULT')
    
    outdf2.wide <- reshape(outdf2, idvar = c('SAMPID'), direction = 'wide',
                           timevar = 'PARAMETER', v.names = 'RESULT')
    names(outdf2.wide) <- gsub("RESULT\\.", "", names(outdf2.wide))
    
    
    outdf2a <- reshape(outdf2.wide, idvar = 'SAMPID', direction = 'long',
                       varying = c('N_PEAT_MOSS_DOM','FREQ_PEAT_MOSS_DOM'),
                       times = c('N_PEAT_MOSS_DOM','FREQ_PEAT_MOSS_DOM'),
                       timevar = 'METRIC', v.names = 'RESULT')

    outdf3 <- rbind(outdf1,outdf2a)
  }else{
    outdf3 <- outdf1
  }

  outdf4 <- reshape(outdf3, idvar = c('SAMPID'), direction = 'wide',
                    timevar = 'METRIC', v.names = 'RESULT')
  names(outdf4) <- gsub("RESULT\\.", "", names(outdf4))

  empty_nv <- data.frame(t(rep(NA,17)),stringsAsFactors=F)
  if(survyear=='2011'){
    names(empty_nv) <- c("FREQ_ALGAE", "FREQ_ARBOREAL","FREQ_BRYOPHYTES","FREQ_LICHENS","FREQ_MACROALGAE","IMP_ALGAE"
                         ,"IMP_ARBOREAL","IMP_BRYOPHYTES","IMP_LICHENS","IMP_MACROALGAE","XCOV_ALGAE","XCOV_ARBOREAL","XCOV_BRYOPHYTES"
                         ,"XCOV_LICHENS","XCOV_MACROALGAE","N_PEAT_MOSS_DOM","FREQ_PEAT_MOSS_DOM")
  }else{
    names(empty_nv) <- c("FREQ_ALGAE", "FREQ_BRYOPHYTES","FREQ_LICHENS","FREQ_MACROALGAE","IMP_ALGAE"
                         ,"IMP_BRYOPHYTES","IMP_LICHENS","IMP_MACROALGAE","XCOV_ALGAE","XCOV_BRYOPHYTES"
                         ,"XCOV_LICHENS","XCOV_MACROALGAE","N_PEAT_MOSS_DOM","FREQ_PEAT_MOSS_DOM")
  }
  outdf5 <- subset(merge(outdf4, empty_nv, all=TRUE),!is.na(SAMPID))

  outdf6 <- merge(allUIDs,outdf5,by='SAMPID',all.x=T)

  varNames <- names(outdf6)[!names(outdf6) %in% c('SAMPID')]
  outdf7 <- reshape(outdf6, idvar = 'SAMPID', direction = 'long',
                    varying = varNames, times = varNames,
                    timevar = 'METRIC', v.names = 'RESULT')
  outdf7$METRIC <- with(outdf7, as.character(METRIC))
  outdf7$RESULT <- with(outdf7, ifelse(is.na(RESULT),0,RESULT))
  print("Done with non-vascular strata metrics")

  nvOut <- merge(samples,outdf7,by='SAMPID') 
  nvOut$PARAMETER <- nvOut$METRIC
  nvOut <- nvOut[, c(sampID, 'PARAMETER', 'RESULT')]
  
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
#' calculating water cover metrics: 'TIME', 'MINIMUM_DEPTH', 
#' 'MAXIMUM_DEPTH', 'PREDOMINANT_DEPTH', 'TOTAL_WATER', 
#' 'WATER_NOVEG', 'WATER_AQVEG', 'WATER_EMERGVEG'.
#' Additional parameters or variables are ignored.
#' WATER_NOVEG, WATER_AQVEG, MINIMUM_DEPTH, and 
#' MAXIMUM_DEPTH not measured in 2016. 
#' @param nPlot A data frame with the 
#' number of plots sampled associated with each
#' sample with \emph{sampID} variables and NPLOTS
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples,
#' 'UID' by default.
#' @param survyear A string with the survey year. The default is '2011'.
#' Starting in 2016, only TIME, TOTAL_WATER, and PREDOMINANT_DEPTH is 
#' measured.
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
#' nplots <- aggregate(x = list(NPLOTS = Vtype_GrCovEx$PLOT), 
#' by = Vtype_GrCovEx[c('UID')],
#' FUN = function(x){length(unique(x))})
#'
#' wcovEx <- calcWcovMets(Vtype_GrCovEx,nplots)
#'
#' head(wcovEx)
#' unique(wcovEx$METRIC)
calcWcovMets <- function(dataIn, nPlot, sampID='UID', survyear='2011'){
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
  wdep$RESULT <- with(wdep, as.numeric(RESULT))
  wdep <- subset(wdep, select = c('SAMPID', 'PLOT', 'NPLOTS', 'PARAMETER', 'RESULT'))
  
  wdep1 <- reshape(wdep, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'wide',
                   timevar = 'PARAMETER', v.names = 'RESULT')
  names(wdep1) <- gsub("RESULT\\.", "", names(wdep1))

  wdep2 <- subset(wdep,RESULT %nin% c("") & !is.na(RESULT) & PARAMETER %in% c('TOTAL_WATER','WATER_NOVEG','WATER_AQVEG','WATER_EMERGVEG'))
  
  wdep2a.wide <- reshape(wdep2, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'wide',
                         timevar = 'PARAMETER', v.names = 'RESULT')
  names(wdep2a.wide) <- gsub("RESULT\\.", "", names(wdep2a.wide))
  
  varNames <- names(wdep2a.wide)[!names(wdep2a.wide) %in% c('SAMPID','PLOT','NPLOTS')]
  wdep2a <- reshape(wdep2a.wide, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'long',
                    times = varNames, varying = varNames,
                    timevar = 'PARAMETER', v.names = 'RESULT')
  wdep2a$RESULT <- with(wdep2a, ifelse(is.na(RESULT),0,as.numeric(RESULT)))

  if(survyear=='2011'){
    wat1.min <- aggregate(x = list(MIN_H2O_DEPTH = wdep1$MINIMUM_DEPTH, MIN_COV_H2O = wdep1$TOTAL_WATER),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x)),min(x,na.rm=TRUE),NA)})
    
    wat1.max <- aggregate(x = list(MAX_H2O_DEPTH = wdep1$MAXIMUM_DEPTH, MAX_COV_H2O = wdep1$TOTAL_WATER),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x)),max(x,na.rm=TRUE),NA)})
    
    wat1.sum <- aggregate(x = list(XH2O_DEPTH_AA = wdep1$PREDOMINANT_DEPTH),
                           by = wdep1[c('SAMPID','NPLOTS')],
                           FUN = function(x){ifelse(any(!is.na(x))
                                                    ,sum(x,na.rm=TRUE),NA)})
    wat1.sum$XH2O_DEPTH_AA <- with(wat1.sum, round(XH2O_DEPTH_AA/NPLOTS, 2))
    
    
    wat1 <- merge(wat1.min, wat1.max, by = c('SAMPID','NPLOTS'))
    wat1 <- merge(wat1, wat1.sum, by = c('SAMPID','NPLOTS'))
  }else{
    wat1.min <- aggregate(x = list(MIN_COV_H2O = wdep1$TOTAL_WATER),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x)),min(x,na.rm=TRUE),NA)})
    
    wat1.max <- aggregate(x = list(MAX_COV_H2O = wdep1$TOTAL_WATER),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x)),max(x,na.rm=TRUE),NA)})
    
    wat1.sum <- aggregate(x = list(XH2O_DEPTH_AA = wdep1$PREDOMINANT_DEPTH),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x))
                                                   ,sum(x,na.rm=TRUE),NA)})
    wat1.sum$XH2O_DEPTH_AA <- with(wat1.sum, round(XH2O_DEPTH_AA/NPLOTS, 2))
    
    wat1 <- merge(wat1.min, wat1.max, by = c('SAMPID','NPLOTS'))
    wat1 <- merge(wat1, wat1.sum, by = c('SAMPID','NPLOTS'))
    

  }
  
  varNames <- names(wat1)[!names(wat1) %in% c('SAMPID','NPLOTS')]
  wat1a <- reshape(wat1, idvar = c('SAMPID'), direction = 'long',
                   times = varNames, varying = varNames,
                   timevar = 'METRIC', v.names = 'RESULT')

  ## Fix Inf values to 0s
  wat1a$RESULT <- with(wat1a, ifelse(RESULT %in% c('Inf','-Inf'),0,ifelse(RESULT=='Inf_-Inf','',RESULT)))
  wat1a$METRIC <- with(wat1a, as.character(METRIC))
  wat1a$NPLOTS <- NULL

  wdep2a.sub <- subset(wdep2a, as.numeric(RESULT)>0)
  
  wat2.sum <- aggregate(x = list(XCOV_H2O = wdep2a.sub$RESULT),
                        by = wdep2a.sub[c('SAMPID','PARAMETER','NPLOTS')],
                        FUN = function(x){sum(as.numeric(x))})
  
  wat2.freq <- aggregate(x = list(FREQ_H2O = wdep2a.sub$RESULT),
                         by = wdep2a.sub[c('SAMPID','PARAMETER','NPLOTS')],
                         FUN = function(x){length(as.numeric(x))})
  
  wat2 <- merge(wat2.sum, wat2.freq, by = c('SAMPID','PARAMETER','NPLOTS'), all=TRUE)
  
  wat2$FREQ_H2O <- with(wat2, round(FREQ_H2O/NPLOTS*100, 2))
  wat2$XCOV_H2O <- with(wat2, round(XCOV_H2O/NPLOTS, 2))
  wat2$IMP_H2O <- with(wat2, round((FREQ_H2O + XCOV_H2O)/2, 2))
  wat2$NPLOTS <- NULL

  wat2a <- reshape(wat2, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                   times = c('FREQ_H2O','XCOV_H2O','IMP_H2O'), 
                   varying = c('FREQ_H2O','XCOV_H2O','IMP_H2O'),
                   timevar = 'METRIC', v.names = 'RESULT')
  wat2a$METRIC <- with(wat2a, ifelse(PARAMETER=='TOTAL_WATER', as.character(METRIC), 
                                     as.character(paste(METRIC,substring(PARAMETER,7),sep='_'))))
  
  wat2a$PARAMETER <- NULL

  ## Now pull only PREDOMINANT_DEPTH that is not missing and >0 to calculate XH2O_DEPTH
  wdep1.pos <- subset(wdep1,PREDOMINANT_DEPTH>0 & !is.na(PREDOMINANT_DEPTH))
  
  wat3 <- aggregate(x = list(RESULT = wdep1.pos$PREDOMINANT_DEPTH), 
                    by = wdep1.pos[c('SAMPID')], 
                    FUN = function(x){mean(x, na.rm = TRUE)})
  wat3$METRIC <- 'XH2O_DEPTH'

  watMet <- rbind(wat1a,wat2a,wat3)
  
  watMet1 <- reshape(watMet, idvar = 'SAMPID', direction = 'wide',
                     timevar = 'METRIC', v.names = 'RESULT')
  names(watMet1) <- gsub("RESULT\\.", "", names(watMet1)) 

  if(survyear=='2011'){
    empty_wat <- data.frame(t(rep(NA,18)),stringsAsFactors=F)
    names(empty_wat) <- c("MIN_H2O_DEPTH","MAX_H2O_DEPTH","XH2O_DEPTH_AA","MIN_COV_H2O",
                          "MAX_COV_H2O","FREQ_H2O","FREQ_H2O_AQVEG","FREQ_H2O_EMERGVEG",
                          "FREQ_H2O_NOVEG","XCOV_H2O","XCOV_H2O_AQVEG","XCOV_H2O_EMERGVEG",
                          "XCOV_H2O_NOVEG","IMP_H2O","IMP_H2O_AQVEG","IMP_H2O_EMERGVEG",
                          "IMP_H2O_NOVEG","XH2O_DEPTH")
  }else{
    empty_wat <- data.frame(t(rep(NA,7)),stringsAsFactors=F)
    names(empty_wat) <- c("XH2O_DEPTH_AA","MIN_COV_H2O","MAX_COV_H2O","FREQ_H2O",
                          "XCOV_H2O","IMP_H2O","XH2O_DEPTH")
  }
  watMet2 <- subset(merge(watMet1, empty_wat, all=TRUE),!is.na(SAMPID))

  watMet3 <- merge(allUIDs,watMet2,by='SAMPID',all.x=T)

  varNames <- names(watMet3)[!names(watMet3) %in% c('SAMPID','NPLOTS')]
  watMet4 <- reshape(watMet3, idvar = 'SAMPID', direction = 'long',
                     times = varNames, varying = varNames,
                     timevar = 'METRIC', v.names = 'RESULT')
  
  watMet4$METRIC <- with(watMet4, as.character(METRIC))
  watMet4$RESULT <- with(watMet4, ifelse(is.na(RESULT),'0',RESULT))

  print("Done with water depth metrics")
  watMetOut <- merge(samples,watMet4,by='SAMPID') 
  watMetOut$PARAMETER <- watMetOut$METRIC
  watMetOut <- watMetOut[, c(sampID, 'PARAMETER', 'RESULT')]

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
#' calculating ground cover metrics: 'LITTER_THATCH', 'LITTER_FORB', 
#' 'LITTER_CONIFER', 'LITTER_DECID', 'LITTER_BROADLEAF', 
#' 'LITTER_DEPTH_SW', 'LITTER_DEPTH_NE', 'DEPTH_SW' (2016), 
#' 'DEPTH_NE' (2016), TOTAL_LITTER', 'WD_FINE', 
#' 'WD_COARSE', 'EXPOSED_SOIL', 'EXPOSED_GRAVEL', 'EXPOSED_ROCK',
#'  'PREDOMINANT_LITTER' (2016).
#' 
#' Additional parameters or variables are ignored.
#' @param nPlot A data frame with the 
#' number of plots sampled associated with each sample,
#' including \emph{sampID} variables and NPLOTS.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @param survyear A string with the survey year. The default is '2011'.
#' Starting in 2016, parameter structure changed such that, instead of
#' individual LITTER_XXXX parameters to indicate the preominant litter
#' type, PREDOMINANT_LITTER was used. In addition, DEPTH_NE and DEPTH_SW
#' replaced LITTER_DEPTH_NE and LITTER_DEPTH_SW.
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
#' nplots <- aggregate(x = list(NPLOTS = Vtype_GrCovEx$PLOT), 
#' by = Vtype_GrCovEx[c('UID')],
#' FUN = function(x){length(unique(x))})
#'
#' bgEx <- calcBareGround_LitterMets(Vtype_GrCovEx,nplots)
#'
#' head(bgEx)
#' unique(bgEx$PARAMETER)

calcBareGround_LitterMets <- function(dataIn,nPlot,sampID='UID',survyear='2011'){
  ## Now merge back with input df
  dataIn1 <- merge(dataIn, nPlot, by=sampID)

  # Create vector of all samples in dataset
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID,dataIn1[,sampID[i]],sep='.')
  }
  samples <- unique(subset(dataIn1,select=c(sampID,'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID),stringsAsFactors=F)
  
  ####### LITTER TYPES
  # Need to calculate the number of quadrats sampled using the NE and SW parameters
  litter.sub <- subset(dataIn1,PARAMETER %in% c('LITTER_DEPTH_NE','LITTER_DEPTH_SW','DEPTH_NE','DEPTH_SW')) # Only one set of parameters or the other will be in any dataset
  
  numQuads <- aggregate(x = list(NQUADS = litter.sub$RESULT),
                        by = litter.sub[c('SAMPID')], FUN = length)

  litter1 <- subset(dataIn1,PARAMETER %in% c('LITTER_THATCH','LITTER_FORB','LITTER_CONIFER',
                                             'LITTER_DECID','LITTER_BROADLEAF',
                                             'LITTER_DEPTH_SW','LITTER_DEPTH_NE',
                                             'TOTAL_LITTER','WD_FINE','WD_COARSE',
                                             'EXPOSED_SOIL','EXPOSED_GRAVEL',
                                             'EXPOSED_ROCK','DEPTH_SW','DEPTH_NE',
                                             'PREDOMINANT_LITTER'))

  litter2 <- merge(litter1, numQuads, by='SAMPID')
  litter2 <- subset(litter2, select = c('SAMPID', 'PLOT','PARAMETER','RESULT','NPLOTS','NQUADS'))

  if(survyear=='2011'){
    littype <- subset(litter2, PARAMETER %in% c('LITTER_THATCH','LITTER_FORB','LITTER_CONIFER','LITTER_DECID',
                                               'LITTER_BROADLEAF','LITTER_NONE') & !is.na(RESULT))
    
    rr1 <- aggregate(x = list(N_LITTER_TYPE = littype$PARAMETER),
                     by = littype[c('SAMPID')],
                     FUN = function(x){length(unique(x))})
    
    rr1 <- merge(littype, rr1, by = 'SAMPID')
    rr2 <- aggregate(x = list(NUM = rr1$PLOT), by = rr1[c('SAMPID','PARAMETER','N_LITTER_TYPE')],
                     FUN = length)
    rr3 <- aggregate(x = list(MAXN = rr2$NUM), by = rr2[c('SAMPID')], FUN = max)
    rr3 <- merge(rr2, rr3, by = 'SAMPID')
    rr4 <- subset(rr3, MAXN==NUM)
    rr5 <- aggregate(x = list(LITTER_TYPE = rr4$PARAMETER), 
                     by = rr4[c('SAMPID','N_LITTER_TYPE')],
                     FUN = function(x){paste(x,collapse='_')})
    rr5$LITTER_TYPE <- with(rr5, gsub('LITTER_','',LITTER_TYPE))
  }else{
    littype <- subset(litter2, PARAMETER=='PREDOMINANT_LITTER' & !is.na(RESULT)) 
    
    rr1 <- aggregate(x = list(N_LITTER_TYPE = littype$RESULT),
                     by = littype[c('SAMPID')],
                     FUN = function(x){length(unique(x))})
    
    rr1 <- merge(littype, rr1, by = 'SAMPID')
    rr2 <- aggregate(x = list(NUM = rr1$PLOT), by = rr1[c('SAMPID','RESULT','N_LITTER_TYPE')],
                     FUN = length)
    rr3 <- aggregate(x = list(MAXN = rr2$NUM), by = rr2[c('SAMPID')], FUN = max)
    rr3 <- merge(rr2, rr3, by = 'SAMPID')
    rr4 <- subset(rr3, MAXN==NUM)
    rr5 <- aggregate(x = list(LITTER_TYPE = rr4$RESULT), 
                     by = rr4[c('SAMPID','N_LITTER_TYPE')],
                     FUN = function(x){paste(x,collapse='_')})

  }
  ## to determine median depth, we must account for any quadrats without depth recorded
  litdep <- subset(litter2,PARAMETER %in% c('LITTER_DEPTH_SW','LITTER_DEPTH_NE','DEPTH_NE','DEPTH_SW'))

  ss1.sum <- aggregate(x = list(XDEPTH_LITTER = litdep$RESULT), 
                       by = litdep[c('SAMPID','NQUADS')],
                       FUN = function(x){sum(as.numeric(x))})
  
  ss1.length <- aggregate(x = list(NSAMP = litdep$RESULT), 
                          by = litdep[c('SAMPID','NQUADS')],
                          FUN = length)
  
  ss1 <- merge(ss1.sum, ss1.length, by = c('SAMPID','NQUADS'))
  ss1$XDEPTH_LITTER <- with(ss1, round(XDEPTH_LITTER/NQUADS, 2))

  litdep1 <- aggregate(x = list(NSAMP = litdep$RESULT), 
                       by = litdep[c('SAMPID','NPLOTS','NQUADS')],
                       FUN = length)
  
  litdep1 <- merge(litdep, litdep1, by = c('SAMPID','NPLOTS','NQUADS'))
  
  tt <- aggregate(x = list(MEDDEPTH_LITTER = litdep1$RESULT),
                  by = litdep1[c('SAMPID')], 
                  FUN = function(x){median(as.numeric(x))})

  varNames <- names(rr5)[!names(rr5) %in% 'SAMPID']
  rr5.long <- reshape(rr5, idvar = 'SAMPID', direction = 'long',
                      times = varNames, varying = varNames,
                      timevar = 'METRIC', v.names = 'RESULT')
  
  ss1.long <- reshape(ss1, idvar = 'SAMPID', direction = 'long',
                      times = 'XDEPTH_LITTER', varying = 'XDEPTH_LITTER',
                      timevar = 'METRIC', v.names = 'RESULT') 
  ss1.long <- subset(ss1.long, select = c('SAMPID','METRIC','RESULT'))
  
  varNames <- names(tt)[!names(tt) %in% 'SAMPID']
  tt.long <- reshape(tt, idvar = 'SAMPID', direction = 'long',
                      times = varNames, varying = varNames,
                      timevar = 'METRIC', v.names = 'RESULT')
  
  loutdf <- rbind(rr5.long, ss1.long, tt.long)
  loutdf$METRIC <- with(loutdf, as.character(METRIC))

  loutdf1 <- reshape(loutdf, idvar = 'SAMPID', direction = 'wide',
                     timevar = 'METRIC', v.names = 'RESULT')
  names(loutdf1) <- gsub("RESULT\\.", "", names(loutdf1))

  empty_lit <- data.frame(t(rep(NA,4)),stringsAsFactors=FALSE)
  names(empty_lit) <- c('N_LITTER_TYPE','LITTER_TYPE','XDEPTH_LITTER','MEDDEPTH_LITTER')

  loutdf2 <- subset(merge(loutdf1, empty_lit, all=TRUE),!is.na(SAMPID))

  loutdf3 <- merge(allUIDs,loutdf2,by='SAMPID',all.x=T)

  varNames <- names(loutdf3)[!names(loutdf3) %in% 'SAMPID']
  loutdf4 <- reshape(loutdf3, idvar = 'SAMPID', direction = 'long',
                     times = varNames, varying = varNames,
                     timevar = 'METRIC', v.names = 'RESULT')
  loutdf4$METRIC <- with(loutdf4, as.character(METRIC))

  print("Done with litter metrics")
  litterOut <- reshape(loutdf4, idvar = 'SAMPID', direction = 'wide',
                       timevar = 'METRIC', v.names = 'RESULT') 
  names(litterOut) <- gsub("RESULT\\.", "", names(litterOut))

  ################# BARE GROUND
  bgrd <- subset(dataIn1,PARAMETER %in% c('EXPOSED_SOIL','EXPOSED_GRAVEL','EXPOSED_ROCK',
                                          'WD_FINE','WD_COARSE','TOTAL_LITTER') & 
                   RESULT!=0 & !is.na(RESULT), select=c('SAMPID','PLOT','PARAMETER','RESULT','NPLOTS'))
  ## Need to create values for new PARAMETER='BAREGD' based on occurrence of either EXPOSED_SOIL, EXPOSED_GRAVEL, or EXPOSED_ROCK at site
  bgrd.sub <- subset(bgrd, PARAMETER %in% c('EXPOSED_SOIL','EXPOSED_GRAVEL','EXPOSED_ROCK'))
  
  bgrd1 <- aggregate(x = list(RESULT = bgrd.sub$RESULT), 
                     by = bgrd.sub[c('SAMPID','PLOT','NPLOTS')],
                     FUN = function(x){as.character(sum(as.numeric(x)))})
  bgrd1$PARAMETER <- 'BAREGD'
  bgrdIn <- rbind(bgrd, bgrd1)
  ## Need to fill in zeros if plot sampled and variable is zero
  bgrdIn1.wide <- reshape(bgrdIn, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'wide',
                     timevar = 'PARAMETER', v.names = 'RESULT')
  names(bgrdIn1.wide) <- gsub("RESULT\\.", "", names(bgrdIn1.wide))
  
  varNames <- names(bgrdIn1.wide)[!names(bgrdIn1.wide) %in% c('SAMPID','PLOT','NPLOTS')]
  bgrdIn1 <- reshape(bgrdIn1.wide, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'long',
                     times = varNames, varying = varNames,
                     timevar = 'PARAMETER', v.names = 'RESULT')
  bgrdIn1$RESULT <- with(bgrdIn1, ifelse(is.na(RESULT),0,as.numeric(RESULT)))

  bgrdIn1.sub <- subset(bgrdIn1,RESULT!=0)
  bgindf1.freq <- aggregate(x = list(FREQ = bgrdIn1.sub$RESULT),
                            by = bgrdIn1.sub[c('SAMPID','PARAMETER','NPLOTS')],
                            FUN = length)
    
  bgindf1.sum <- aggregate(x = list(XCOV = bgrdIn1.sub$RESULT),
                           by = bgrdIn1.sub[c('SAMPID','PARAMETER','NPLOTS')],
                           FUN = function(x){sum(as.numeric(x))})
 
  bgindf1 <- merge(bgindf1.sum, bgindf1.freq, by = c('SAMPID','PARAMETER','NPLOTS'), all=TRUE)
  bgindf1$FREQ <- with(bgindf1, round(FREQ/NPLOTS*100, 2))
  bgindf1$XCOV <- with(bgindf1, round(XCOV/NPLOTS, 2))
  bgindf1$IMP <- with(bgindf1, round((FREQ + XCOV)/2, 2))
 
  
  bgoutdf <- reshape(bgindf1, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                     times = c('XCOV','FREQ','IMP'), varying = c('XCOV','FREQ','IMP'),
                     timevar = 'variable', v.names = 'RESULT')
  
  bgoutdf$PARAMETER <- with(bgoutdf, ifelse(PARAMETER=='TOTAL_LITTER',
                                      paste(variable,'LITTER',sep='_'),paste(variable,PARAMETER,sep='_')))
  bgoutdf <- bgoutdf[,c('SAMPID','PARAMETER','RESULT')]
  
  bgoutdf1.wide <- reshape(bgoutdf, idvar = 'SAMPID', direction = 'wide',
                           timevar = 'PARAMETER', v.names = 'RESULT')
  names(bgoutdf1.wide) <- gsub("RESULT\\.", "", names(bgoutdf1.wide))
  
  varNames <- names(bgoutdf1.wide)[!names(bgoutdf1.wide) %in% c('SAMPID')]
  bgoutdf1 <- reshape(bgoutdf1.wide, idvar = c('SAMPID'), direction = 'long',
                      times = varNames, varying = varNames,
                      timevar = 'METRIC', v.names = 'RESULT')
  
  bgoutdf1$RESULT <- with(bgoutdf1, ifelse(is.na(RESULT),0,RESULT))
  bgoutdf1$METRIC <- with(bgoutdf1, as.character(METRIC))
  bgoutdf2 <- reshape(bgoutdf1, idvar = 'SAMPID', direction = 'wide',
                      timevar = 'METRIC', v.names = 'RESULT')
  names(bgoutdf2) <- gsub("RESULT\\.", "", names(bgoutdf2))

  empty_bg <- data.frame(t(rep(NA,21)),stringsAsFactors=FALSE)
  names(empty_bg) <- c("FREQ_BAREGD","FREQ_EXPOSED_GRAVEL","FREQ_EXPOSED_ROCK",
                       "FREQ_EXPOSED_SOIL","FREQ_LITTER","FREQ_WD_COARSE",
                       "FREQ_WD_FINE","IMP_BAREGD","IMP_EXPOSED_GRAVEL",
                       "IMP_EXPOSED_ROCK","IMP_EXPOSED_SOIL","IMP_LITTER",
                       "IMP_WD_COARSE","IMP_WD_FINE","XCOV_BAREGD",
                       "XCOV_EXPOSED_GRAVEL","XCOV_EXPOSED_ROCK","XCOV_EXPOSED_SOIL",
                       "XCOV_LITTER","XCOV_WD_COARSE","XCOV_WD_FINE")

  bgoutdf3 <- subset(merge(bgoutdf2, empty_bg, all=TRUE),!is.na(SAMPID))

  bgoutdf4 <- merge(allUIDs,bgoutdf3,by='SAMPID',all.x=T)

  varNames <- names(bgoutdf4)[!names(bgoutdf4) %in% c('SAMPID')]
  bgoutdf5 <- reshape(bgoutdf4, idvar = 'SAMPID', direction = 'long',
                      times = varNames, varying = varNames,
                      timevar = 'METRIC', v.names = 'RESULT')
  bgoutdf5$METRIC <- with(bgoutdf5, as.character(METRIC))
  bgoutdf5$RESULT <- with(bgoutdf5, ifelse(is.na(RESULT),0,RESULT))

  print("Done with bare ground metrics")
  bgrdOut <- reshape(bgoutdf5, idvar = 'SAMPID', direction = 'wide',
                     timevar = 'METRIC', v.names = 'RESULT')
  names(bgrdOut) <- gsub("RESULT\\.", "", names(bgrdOut))

  # Now combine bare ground and litter metrics to perform a final check
  combOut <- merge(litterOut,bgrdOut,by='SAMPID',all=T)
  
  varNames <- names(combOut)[!names(combOut) %in% c('SAMPID')]
  combOut.1 <- reshape(combOut, idvar = 'SAMPID', direction = 'long',
                       times = varNames, varying = varNames,
                       timevar = 'variable', v.names = 'value')
  combOut.1$value <- with(combOut.1, ifelse(is.na(value) & variable!='LITTER_TYPE',0,value))
  # Need to determine whether missing LITTER_TYPE value should be ABSENT or NONE
  combOut.wide <- reshape(combOut.1, idvar = 'SAMPID', direction = 'wide',
                          timevar = 'variable', v.names = 'value')
  names(combOut.wide) <- gsub("value\\.", "", names(combOut.wide))
  combOut.wide$LITTER_TYPE <- with(combOut.wide, ifelse(is.na(LITTER_TYPE) & XDEPTH_LITTER==0 & XCOV_LITTER==0 & FREQ_LITTER==0, 'ABSENT', ifelse(is.na(LITTER_TYPE),'NONE',LITTER_TYPE)))
  varNames <- names(combOut.wide)[!names(combOut.wide) %in% c('SAMPID')]
  combOut.long <- reshape(combOut.wide, idvar = 'SAMPID', direction = 'long',
                       times = varNames, varying = varNames,
                       timevar = 'PARAMETER', v.names = 'RESULT')

  finOut <- merge(samples, combOut.long, by='SAMPID') 
  finOut <- finOut[c(sampID, 'PARAMETER', 'RESULT')]
    
  return(finOut)
}

