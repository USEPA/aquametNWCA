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
  # Merge dataIn with nPlot
  dataIn1 <- merge(dataIn,nPlot,by=sampID)
  
  # Create SAMPID variable by combining variables in sampID argument
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID,dataIn1[,sampID[i]],sep='.')
  }
  samples <- unique(subset(dataIn1,select=c(sampID,'SAMPID')))
  # Create vector of all samples in dataset
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID),stringsAsFactors=F)
  # Cast data frame and drop prefix added by reshape()
  vhet <- reshape(dataIn1, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'wide',
                  timevar = 'PARAMETER', v.names = 'RESULT')
  names(vhet) <- gsub("RESULT\\.", "", names(vhet))
  
  # First calculate heterogeneity metrics 
  # Account for datasets without PF parameter so that wide form will not have it
  if('PAL_FARMED' %in% names(vhet)==TRUE){
    vhet$SANDT_CLASS <- with(vhet, ifelse(is.na(PAL_FARMED),SANDT_CLASS,'PF'))
  }
  # Drop observations without SANDT_CLASS value
  vhet <- subset(vhet, !is.na(SANDT_CLASS)) 
  # Calculate number of unique S & T classes, merge with data above
  vhet1.n <- aggregate(x = list(N_SANDT = vhet$SANDT_CLASS), 
                       by = vhet[c('SAMPID')], FUN = function(x){length(unique(x))})
  vhet1 <- merge(vhet, vhet1.n, by='SAMPID')
  # Calculate number of plots by class
  vhet2 <- aggregate(x = list(FREQ = vhet1$PLOT), 
                     by = vhet1[c('SAMPID','NPLOTS','SANDT_CLASS','N_SANDT')],
                     FUN = length)
  # Calculate frequency from the above
  vhet2$FREQ <- with(vhet2, FREQ/NPLOTS)
  # Calculate max frequency and then match back to frequencies to find dominant class
  vhet3.maxf <- aggregate(x = list(MAXF = vhet2$FREQ), by = vhet2[c('SAMPID')], FUN = max)
  vhet3 <- merge(vhet2, vhet3.maxf, by='SAMPID')
  vhet4 <- subset(vhet3, MAXF==FREQ)
  # Must aggregate in case there are two dominant classes
  vhet4a <- aggregate(x = list(DOM_SANDT = vhet4$SANDT_CLASS), 
                      by = vhet4[c('SAMPID','N_SANDT')],
                      FUN = function(x){paste(x, collapse = '-')})
  # Calculate D (Simpson diversity) of S&T classes)
  vhet5.d <- aggregate(x = list(D_SANDT = vhet2$FREQ), by = vhet2[c('SAMPID','NPLOTS')],
                       FUN = function(x){round(1 - sum(x*x), 4)})
  # Calculate Shannon Wiener Diversity (H) 
  vhet5.h <- aggregate(x = list(H_SANDT = vhet2$FREQ), by = vhet2[c('SAMPID','NPLOTS')],
                       FUN = function(x){round(-1*sum(x*log(x)),4)})
  # Calculate first step in Pielou's Evenness (J)
  vhet5.jcalc <- aggregate(x = list(uniqST = vhet2$N_SANDT), by = vhet2[c('SAMPID','NPLOTS')],
                           FUN = function(x){unique(x)}) 
  # Merge data frames
  vhet5 <- merge(vhet5.d, vhet5.h, by=c('SAMPID','NPLOTS'))
  vhet5 <- merge(vhet5, vhet5.jcalc, by=c('SAMPID','NPLOTS'))
  # Finish calculation of J
  vhet5$J_SANDT <- with(vhet5, round(ifelse(H_SANDT!=0 & uniqST!=1,H_SANDT/log(uniqST),0),4))
  # Merge with previously calculated metric data frame
  vhet6 <- merge(vhet4a,vhet5,by='SAMPID')
  # create an empty data frame with all of the metric names as variables to ensure all are included in output
  empty_vhet <- data.frame(t(rep(NA,5)),stringsAsFactors=FALSE)
  names(empty_vhet) <- c('N_SANDT','DOM_SANDT','D_SANDT','H_SANDT','J_SANDT')
  # Merge with output data frame and drop row with missing SAMPID
  vhet7 <- subset(merge(vhet6, empty_vhet, all=TRUE),!is.na(SAMPID))
  # Now merge with full list of samples 
  vhet8 <- merge(allUIDs, vhet7, by='SAMPID', all.x=T)
  # Melt data frame and replace missing RESULT values with 0 except for dominant S&T
  vhet9 <- reshape(vhet8, idvar = c('SAMPID'), direction = 'long',
                   times = c('D_SANDT','H_SANDT','J_SANDT','DOM_SANDT','N_SANDT'),
                   varying = c('D_SANDT','H_SANDT','J_SANDT','DOM_SANDT','N_SANDT'),
                   timevar = 'METRIC', v.names = 'RESULT')
  vhet9$METRIC <- with(vhet9, as.character(METRIC))
  vhet9$RESULT <- with(vhet9, ifelse(METRIC=='DOM_SANDT' & is.na(RESULT),'MISSING'
                                     ,ifelse(is.na(RESULT),0,RESULT)))

  print("Done with veg heterogeneity metrics")
  # Merge full samples list with output metrics
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
calcVascStratMets <- function(dataIn, nPlot, sampID='UID'){
  # Merge dataIn with nPlot
  dataIn1 <- merge(dataIn, nPlot, by=sampID)

  # Create SAMPID variable based on combination of variables in sampID
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID, dataIn1[,sampID[i]],sep='.')
  }
  # select unique set of samples along with SAMPID
  samples <- unique(subset(dataIn1, select=c(sampID, 'SAMPID')))
  # Create vector of all samples in dataset
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID), stringsAsFactors=F)
  # Subset to just parameters for vascular strata and drop any with 0 values or missing RESULT
  vstrat <- subset(dataIn1, PARAMETER %in% c('SUBMERGED_AQ','FLOATING_AQ','LIANAS','VTALL_VEG',
                                            'TALL_VEG','HMED_VEG','MED_VEG',
                                            'SMALL_VEG','VSMALL_VEG') & !is.na(RESULT) & RESULT!='0')
  # Sum cover for total cover of vascular strata
  vstrat.sum <- aggregate(x = list(XTOTCOV_VASC_STRATA = vstrat$RESULT), 
                          by = vstrat[c('SAMPID','NPLOTS')],
                          FUN = function(x){sum(as.numeric(x))})
  # Complete calculation of mean cover by dividing by number of plots
  vstrat.sum$XTOTCOV_VASC_STRATA <- with(vstrat.sum, XTOTCOV_VASC_STRATA/NPLOTS)
  # Calculate number of vascular strata and sampled plots
  vstrat.cnt <- aggregate(x = list(N_VASC_STRATA = vstrat$PARAMETER, PLOTSAMP = vstrat$PLOT),
                          by = vstrat[c('SAMPID')], FUN = function(x){length(unique(x))})
  # Merge back with previous data frames
  vstrat <- merge(vstrat, vstrat.sum, by=c('SAMPID','NPLOTS'))
  vstrat <- merge(vstrat, vstrat.cnt, by='SAMPID')
  # Number of vascular strata by plot for later calculation
  vstratPlot.n <- aggregate(x = list(N_VSTRATA = vstrat$PARAMETER), 
                            by = vstrat[c('SAMPID','PLOT','N_VASC_STRATA','NPLOTS',
                                          'XTOTCOV_VASC_STRATA','PLOTSAMP')],
                            FUN = function(x){length(unique(x))})
  # Sum cover values by plot
  vstratPlot.sum <- aggregate(x = list(SUM_PLOT = vstrat$RESULT), 
                              by = vstrat[c('SAMPID','PLOT','N_VASC_STRATA','NPLOTS',
                                            'XTOTCOV_VASC_STRATA','PLOTSAMP')],
                              FUN = function(x){sum(as.numeric(x))})
  # Merge these two calculations together
  vstratPlot <- merge(vstratPlot.n, vstratPlot.sum, by=c('SAMPID','PLOT','N_VASC_STRATA',
                                                         'NPLOTS','XTOTCOV_VASC_STRATA','PLOTSAMP'))
  # Calculate mean number of vascular strata
  vstrat1.sum <- aggregate(x = list(XN_VASC_STRATA = vstratPlot$N_VSTRATA),
                           by = vstratPlot[c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP')],
                           FUN = sum)
  # Calculate minimum number of strata
  vstrat1.min <- aggregate(x = list(minN = vstratPlot$N_VSTRATA),
                           by = vstratPlot[c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP')],
                           FUN = min)
  # Calculate maximum number of strata
  vstrat1.max <- aggregate(x = list(maxN = vstratPlot$N_VSTRATA),
                           by = vstratPlot[c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP')],
                           FUN = max)
  # Merge data sets together
  vstrat1 <- merge(vstrat1.sum, vstrat1.min, 
                   by=c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP'))
  vstrat1 <- merge(vstrat1, vstrat1.max, c('SAMPID','N_VASC_STRATA','XTOTCOV_VASC_STRATA','NPLOTS','PLOTSAMP'))
  # Final calculations for mean and range of number of vascular strata
  vstrat1$XN_VASC_STRATA <- with(vstrat1, XN_VASC_STRATA/NPLOTS)
  vstrat1$RG_VASC_STRATA <- with(vstrat1, ifelse(PLOTSAMP==NPLOTS, maxN - minN, maxN - 0))
  # Drop unnecessary variables
  vstrat1$minN <- NULL
  vstrat1$maxN <- NULL
  vstrat1$PLOTSAMP <- NULL
  vstrat1$NPLOTS <- NULL
  # Create empty data frame
  empty_vstrat <- data.frame(t(rep(NA,4)), stringsAsFactors=FALSE)
  names(empty_vstrat) <- c('N_VASC_STRATA','XTOTCOV_VASC_STRATA','XN_VASC_STRATA','RG_VASC_STRATA')
  # Merge empty data frame output data frame from above, drop row with missing SAMPID
  vstrat2 <- subset(merge(vstrat1, empty_vstrat, all=TRUE), !is.na(SAMPID))
  # Merge full set of samples with output to ensure all samples represented
  vstrat3 <- merge(allUIDs, vstrat2, by='SAMPID', all.x=T)
  # Melt data frame, then replace missing RESULT with 0
  varNames <- names(vstrat3)[!names(vstrat3) %in% c('SAMPID')]
  vstratMet <- reshape(vstrat3, idvar = 'SAMPID', direction = 'long',
                       varying = varNames, times = varNames,
                       timevar = 'METRIC', v.names = 'RESULT')
  
  vstratMet$METRIC <- with(vstratMet, as.character(METRIC))
  vstratMet$RESULT <- with(vstratMet, ifelse(is.na(RESULT),0,RESULT))

  ## Calculate frequency, mean cover, relative mean cover, and relative importance by vascular stratum
  vstrat.pos <- subset(vstrat, RESULT!=0)
  # Calculate precursor for FREQ calculation by parameter
  indf1.length <- aggregate(x = list(FREQ = vstrat.pos$PLOT), 
                            by = vstrat.pos[c('SAMPID','PARAMETER','NPLOTS','XTOTCOV_VASC_STRATA')],
                            FUN = length)
  # Sum cover by parameter
  indf1.sum <- aggregate(x = list(XCOV = vstrat.pos$RESULT), 
                         by = vstrat.pos[c('SAMPID','PARAMETER','NPLOTS','XTOTCOV_VASC_STRATA')],
                         FUN = function(x){sum(as.numeric(x))})
  # Merge calculated values together, then calculate additional metrics 
  indf1 <- merge(indf1.length, indf1.sum, by = c('SAMPID','PARAMETER','NPLOTS','XTOTCOV_VASC_STRATA'), all=TRUE)
  indf1$FREQ <- with(indf1, round(FREQ/NPLOTS*100, 2))
  indf1$XCOV <- with(indf1, round(XCOV/NPLOTS, 2))
  indf1$XRCOV <- with(indf1, round(XCOV/XTOTCOV_VASC_STRATA*100, 2))
  indf1$IMP <- with(indf1, round((FREQ + XCOV)/2, 2))
  
  # Now drop variables
  indf1$XTOTCOV_VASC_STRATA <- NULL
  indf1$NPLOTS <- NULL
  # Melt data frame 
  outdf <- reshape(indf1, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                   varying = c('FREQ','XCOV','XRCOV','IMP'), times = c('FREQ','XCOV','XRCOV','IMP'),
                   timevar = 'variable', v.names = 'RESULT')
  # Update parameter value by combining variable from above with PARAMETER
  outdf$PARAMETER <- with(outdf, paste(variable, PARAMETER, sep = '_'))
  # Drop variable
  outdf$variable <- NULL
  # Cast output data frame to wide format and drop prefix added to variable names by reshape()
  outdf.wide <- reshape(outdf, idvar = c('SAMPID'), direction = 'wide',
                        timevar = 'PARAMETER', v.names = 'RESULT')
  names(outdf.wide) <- gsub("RESULT\\.", "", names(outdf.wide))
  # Melt again, now that we have all metrics represented for all samples
  varNames <- names(outdf.wide)[!names(outdf.wide) %in% c('SAMPID')]
  outdf1 <- reshape(outdf.wide, idvar = c('SAMPID'), direction = 'long',
                    varying = varNames, times = varNames,
                    timevar = 'METRIC', v.names = 'RESULT')
  # Set RESULT to 0 where missing
  outdf1$RESULT[is.na(outdf1$RESULT)] <- 0

  # Calculate diversity indices based on vascular strata
  # Start with Shannon diversity (H)
  div1.h <- aggregate(x = list(H_VASC_STRATA = indf1$XRCOV), 
                      by = indf1[c('SAMPID')], FUN = function(x){round(-1*sum((x/100)*log(x/100)), 4)})
  # Now Pielou's evenness - initial step in calculation
  div1.j <- aggregate(x = list(J_VASC_STRATA = indf1$PARAMETER), 
                      by = indf1[c('SAMPID')], FUN = length)
  # Simpson Diversity index
  div1.d <- aggregate(x = list(D_VASC_STRATA = indf1$XRCOV),
                      by = indf1[c('SAMPID')],
                      FUN = function(x){round(1 - sum((x/100)^2), 4)})
  # Merge data frames
  div1 <- merge(div1.h, div1.j, by = 'SAMPID', all=TRUE)
  div1 <- merge(div1, div1.d, by='SAMPID', all=TRUE)
  # Finish calculation of J, then fill in 0 where missing
  div1$J_VASC_STRATA <- with(div1, round(H_VASC_STRATA/log(J_VASC_STRATA), 4))
  div1$J_VASC_STRATA <- with(div1, ifelse(!is.na(J_VASC_STRATA),J_VASC_STRATA,0))
  # Melt data frame 
  div1.long <- reshape(div1, idvar = 'SAMPID', direction = 'long',
                       varying = c('H_VASC_STRATA','J_VASC_STRATA','D_VASC_STRATA'), 
                       times = c('H_VASC_STRATA','J_VASC_STRATA','D_VASC_STRATA'),
                       timevar = 'METRIC', v.names = 'RESULT')
  # Combine with previously calculated metrics above
  outdf2 <- rbind(outdf1, div1.long)
  # Create empty data frame
  empty_vtype <- data.frame(t(rep(NA,39)), stringsAsFactors=FALSE)
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
  # Cast output data frame from above and drop prefix added by reshape()
  outdf3 <- reshape(outdf2, idvar = 'SAMPID', direction = 'wide',
                    timevar = 'METRIC', v.names = 'RESULT')
  
  names(outdf3) <- gsub("RESULT\\.", "", names(outdf3))
  # Merge with empty data frame and drop row with missing SAMPID
  outdf4 <- subset(merge(outdf3, empty_vtype, all=TRUE), !is.na(SAMPID))
  # Merge with full list of samples to make sure all are represented
  outdf5 <- merge(allUIDs, outdf4, by='SAMPID', all.x=T)
  # Melt output data frame and fill in 0 for appropriate metrics where missing values
  varNames <- names(outdf5)[!names(outdf5) %in% c('SAMPID')]
  outdf6 <- reshape(outdf5, idvar = 'SAMPID', direction = 'long',
                    varying = varNames, times = varNames,
                    timevar = 'METRIC', v.names = 'RESULT')
  outdf6$METRIC <- with(outdf6, as.character(METRIC))
  outdf6$RESULT <- with(outdf6, ifelse(METRIC %nin% c('D_VASC_STRATA','H_VASC_STRATA') & is.na(RESULT), 0, RESULT))

  print("Done with vascular strata metrics")

  # Now combine output data frames, rename METRIC to PARAMETER and drop METRIC
  vtOut <- rbind(outdf6, vstratMet) 
  vtOut$PARAMETER <- vtOut$METRIC
  vtOut$METRIC <- NULL
  # Merge with samples to get back to sampID variables
  vtOut.1 <- merge(samples, vtOut, by='SAMPID') 
  # Select subset of variables
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
calcNonvascMets <- function(dataIn, nPlot, sampID='UID', survyear='2011'){
  # Merge nPlot with dataIn by sampID variables
  dataIn1 <- merge(dataIn, nPlot, by=sampID)

  # Create SAMPID variable based on combination of variables in sampID
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID, dataIn1[,sampID[i]], sep='.')
  }
  # select unique set of samples along with SAMPID
  samples <- unique(subset(dataIn1, select=c(sampID, 'SAMPID')))
  # Create vector of all samples in dataset
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID), stringsAsFactors=F)
  
  # Subset dataset to only non-vascular parameters, keep only non-missing and non-zero cover values
  nvstrat <- subset(dataIn1, PARAMETER %in% c('PEAT_MOSS','BRYOPHYTES','LICHENS','ARBOREAL','ALGAE','MACROALGAE') &
                      !is.na(RESULT) & RESULT!='0')
  # This accounts for change in ARBOREAL to ARBOREAL_ABUNDANCE
  nvstrat$PARAMETER <- with(nvstrat, gsub('_ABUNDANCE', '', PARAMETER)) 

  # Now use this data frame to calculate metrics
  # First calculate metrics for non-vascular veg types, excluding peat moss
  nvstrat.sub <- subset(nvstrat,PARAMETER!='PEAT_MOSS')
  # Sum cover by PARAMETER
  indf1.sum <- aggregate(x = list(XCOV = nvstrat.sub$RESULT),
                         by = nvstrat.sub[c('SAMPID','PARAMETER','NPLOTS')],
                         FUN = function(x){sum(as.numeric(x))})
  # Length of plots by parameter to obtain first step in calculating frequency
  indf1.length <- aggregate(x = list(FREQ = nvstrat.sub$PLOT), 
                            by = nvstrat.sub[c('SAMPID','PARAMETER','NPLOTS')],
                            FUN = length)
  # Merge sum and length calculated values
  indf1 <- merge(indf1.sum, indf1.length, by = c('SAMPID','PARAMETER','NPLOTS'), all=TRUE)
  # Calculate metrics for each parameter
  indf1$FREQ <- with(indf1, FREQ/NPLOTS*100)
  indf1$XCOV <- with(indf1, XCOV/NPLOTS)
  indf1$IMP <- with(indf1, (FREQ + XCOV)/2)
  indf1$NPLOTS <- NULL
  # Melt data frame, then update PARAMETER to combine metric name and PARAMETER value
  outdf <- reshape(indf1, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                        varying = c('XCOV','FREQ','IMP'), times = c('XCOV','FREQ','IMP'),
                        timevar = 'variable', v.names = 'RESULT') 
  outdf$PARAMETER <- with(outdf, paste(variable, PARAMETER, sep = '_'))
  # Drop "variable"
  outdf$variable <- NULL
  # Cast data frame and drop prefix from variable names
  outdf.wide <- reshape(outdf, idvar = 'SAMPID', direction = 'wide',
                        timevar = 'PARAMETER', v.names = 'RESULT')
  names(outdf.wide) <- gsub("RESULT\\.", "", names(outdf.wide))
  # Melt again to ensure all samples are included in output
  varNames <- names(outdf.wide)[!names(outdf.wide) %in% c('SAMPID')]
  outdf1 <- reshape(outdf.wide, idvar = "SAMPID", direction = "long",
                    varying = varNames, times = varNames,
                    timevar = 'METRIC', v.names = 'RESULT')
  # Set all missing RESULT values to 0
  outdf1$RESULT[is.na(outdf1$RESULT)] <- 0

  # Now count number and frequency of plots with peat moss dominant (i.e., PEAT_MOSS='Y')
  # Need to account for situations where PEAT_MOSS not present
  if(nrow(subset(nvstrat,PARAMETER=='PEAT_MOSS' & RESULT!='N' & !is.na(RESULT)))>0){
    # Subset to only the parameter PEAT_MOSS present and non-missing result
    nvstrat.peat <- subset(nvstrat,PARAMETER=='PEAT_MOSS' & RESULT!='N' & !is.na(RESULT))
    # Count number of peat moss-dominated plots
    indf2 <- aggregate(x = list(N_PEAT_MOSS_DOM = nvstrat.peat$PLOT),
                       by = nvstrat.peat[c('SAMPID','NPLOTS')],
                       FUN = length)
    # Use the above to calculate frequency of peat moss
    indf2$FREQ_PEAT_MOSS_DOM <- with(indf2, N_PEAT_MOSS_DOM/NPLOTS*100)
    # Drop NPLOTS
    indf2$NPLOTS <- NULL
    # Melt data frame
    outdf2 <- reshape(indf2, idvar = c('SAMPID'), direction = 'long',
                      varying = c('N_PEAT_MOSS_DOM','FREQ_PEAT_MOSS_DOM'),
                      times = c('N_PEAT_MOSS_DOM','FREQ_PEAT_MOSS_DOM'),
                      timevar = 'PARAMETER', v.names = 'RESULT')
    # Cast back to wide format and drop prefix from reshape()
    outdf2.wide <- reshape(outdf2, idvar = c('SAMPID'), direction = 'wide',
                           timevar = 'PARAMETER', v.names = 'RESULT')
    names(outdf2.wide) <- gsub("RESULT\\.", "", names(outdf2.wide))
    
    # Melt data frame 
    outdf2a <- reshape(outdf2.wide, idvar = 'SAMPID', direction = 'long',
                       varying = c('N_PEAT_MOSS_DOM','FREQ_PEAT_MOSS_DOM'),
                       times = c('N_PEAT_MOSS_DOM','FREQ_PEAT_MOSS_DOM'),
                       timevar = 'METRIC', v.names = 'RESULT')
    # Combine output data frame with previouly calculated metrics
    outdf3 <- rbind(outdf1, outdf2a)
  }else{ # no peat moss - do not calculate
    outdf3 <- outdf1
  }
  # Cast combined data frame and drop prefix added by reshape()
  outdf4 <- reshape(outdf3, idvar = c('SAMPID'), direction = 'wide',
                    timevar = 'METRIC', v.names = 'RESULT')
  names(outdf4) <- gsub("RESULT\\.", "", names(outdf4))
  # Create empty data frame - this will vary by year 
  # 2011 has a few more metrics because of changes in field forms after 2011
  empty_nv <- data.frame(t(rep(NA,17)), stringsAsFactors=F)
  if(survyear=='2011'){
    names(empty_nv) <- c("FREQ_ALGAE", "FREQ_ARBOREAL","FREQ_BRYOPHYTES","FREQ_LICHENS",
                         "FREQ_MACROALGAE","IMP_ALGAE"
                         ,"IMP_ARBOREAL","IMP_BRYOPHYTES","IMP_LICHENS","IMP_MACROALGAE",
                         "XCOV_ALGAE","XCOV_ARBOREAL","XCOV_BRYOPHYTES"
                         ,"XCOV_LICHENS","XCOV_MACROALGAE","N_PEAT_MOSS_DOM","FREQ_PEAT_MOSS_DOM")
  }else{
    names(empty_nv) <- c("FREQ_ALGAE", "FREQ_BRYOPHYTES","FREQ_LICHENS","FREQ_MACROALGAE","IMP_ALGAE"
                         ,"IMP_BRYOPHYTES","IMP_LICHENS","IMP_MACROALGAE","XCOV_ALGAE","XCOV_BRYOPHYTES"
                         ,"XCOV_LICHENS","XCOV_MACROALGAE","N_PEAT_MOSS_DOM","FREQ_PEAT_MOSS_DOM")
  }
  # Merge empty data frame with output data frame, drop row missing SAMPID
  outdf5 <- subset(merge(outdf4, empty_nv, all=TRUE), !is.na(SAMPID))
  # Merge with list of all UIDs to ensure all samples represented
  outdf6 <- merge(allUIDs, outdf5, by='SAMPID', all.x=T)
  # Melt data frame and fill in 0 for missing RESULT values
  varNames <- names(outdf6)[!names(outdf6) %in% c('SAMPID')]
  outdf7 <- reshape(outdf6, idvar = 'SAMPID', direction = 'long',
                    varying = varNames, times = varNames,
                    timevar = 'METRIC', v.names = 'RESULT')
  outdf7$METRIC <- with(outdf7, as.character(METRIC))
  outdf7$RESULT <- with(outdf7, ifelse(is.na(RESULT),0,RESULT))
  print("Done with non-vascular strata metrics")
  # Merge samples with output data frame to bring back sampID variables
  # Rename METRIC to PARAMETER and subset to variables needed
  nvOut <- merge(samples, outdf7, by='SAMPID') 
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
  # Merge nPlot with dataIn by sampID variables
  dataIn1 <- merge(dataIn, nPlot, by=sampID)

  # Create SAMPID variable based on combination of variables in sampID
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID, dataIn1[,sampID[i]], sep='.')
  }
  # select unique set of samples along with SAMPID
  samples <- unique(subset(dataIn1, select=c(sampID, 'SAMPID')))
  # Create vector of all samples in dataset
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID), stringsAsFactors=F)
  
  # WATER COVER AND DEPTH
  # Subset data to just water cover parameters
  wdep <- subset(dataIn1, PARAMETER %in% c('TIME','MINIMUM_DEPTH','MAXIMUM_DEPTH','PREDOMINANT_DEPTH','TOTAL_WATER'
                                          ,'WATER_NOVEG','WATER_AQVEG','WATER_EMERGVEG'))
  # Ensure RESULT is numeric
  wdep$RESULT <- with(wdep, as.numeric(RESULT))
  
  # Subset variables
  wdep <- subset(wdep, select = c('SAMPID', 'PLOT', 'NPLOTS', 'PARAMETER', 'RESULT'))
  
  # Cast data wide and drop prefix from reshape()
  wdep1 <- reshape(wdep, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'wide',
                   timevar = 'PARAMETER', v.names = 'RESULT')
  names(wdep1) <- gsub("RESULT\\.", "", names(wdep1))
  
  # Subset to remove missing and blank RESULT values for total water and aquatic vegetation
  # These are not calculated variables, just direct from the field form
  wdep2 <- subset(wdep,RESULT %nin% c("") & !is.na(RESULT) & PARAMETER %in% c('TOTAL_WATER','WATER_NOVEG','WATER_AQVEG','WATER_EMERGVEG'))
  # Cast wide and drop prefix added by reshape()
  wdep2a.wide <- reshape(wdep2, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'wide',
                         timevar = 'PARAMETER', v.names = 'RESULT')
  names(wdep2a.wide) <- gsub("RESULT\\.", "", names(wdep2a.wide))
  # Melt data frame and fill in missing RESULT values with 0
  varNames <- names(wdep2a.wide)[!names(wdep2a.wide) %in% c('SAMPID','PLOT','NPLOTS')]
  wdep2a <- reshape(wdep2a.wide, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'long',
                    times = varNames, varying = varNames,
                    timevar = 'PARAMETER', v.names = 'RESULT')
  wdep2a$RESULT <- with(wdep2a, ifelse(is.na(RESULT), 0, as.numeric(RESULT)))

  # Variables measured in 2011 and 2016 (and later) changed, so some metrics only apply to 2011 
  if(survyear=='2011'){
    # calculate minimum depth and water cover
    wat1.min <- aggregate(x = list(MIN_H2O_DEPTH = wdep1$MINIMUM_DEPTH, MIN_COV_H2O = wdep1$TOTAL_WATER),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x)), min(x, na.rm=TRUE), NA)})
    # Calculate maximum depth and water cover
    wat1.max <- aggregate(x = list(MAX_H2O_DEPTH = wdep1$MAXIMUM_DEPTH, MAX_COV_H2O = wdep1$TOTAL_WATER),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x)), max(x, na.rm=TRUE), NA)})
    # Calculate first step in mean water depth
    wat1.sum <- aggregate(x = list(XH2O_DEPTH_AA = wdep1$PREDOMINANT_DEPTH),
                           by = wdep1[c('SAMPID','NPLOTS')],
                           FUN = function(x){ifelse(any(!is.na(x))
                                                    ,sum(x, na.rm=TRUE), NA)})
    # Complete mean water depth calculation
    wat1.sum$XH2O_DEPTH_AA <- with(wat1.sum, round(XH2O_DEPTH_AA/NPLOTS, 2))
    
    # Merge metric data frames
    wat1 <- merge(wat1.min, wat1.max, by = c('SAMPID', 'NPLOTS'))
    wat1 <- merge(wat1, wat1.sum, by = c('SAMPID', 'NPLOTS'))
  }else{ # If not 2011, more limited data collected for water
    # Minimum water cover
    wat1.min <- aggregate(x = list(MIN_COV_H2O = wdep1$TOTAL_WATER),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x)), min(x,na.rm=TRUE), NA)})
    # Maximum water cover
    wat1.max <- aggregate(x = list(MAX_COV_H2O = wdep1$TOTAL_WATER),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x)), max(x,na.rm=TRUE), NA)})
    # Sum depths for calculating mean water depth
    wat1.sum <- aggregate(x = list(XH2O_DEPTH_AA = wdep1$PREDOMINANT_DEPTH),
                          by = wdep1[c('SAMPID','NPLOTS')],
                          FUN = function(x){ifelse(any(!is.na(x))
                                                   ,sum(x,na.rm=TRUE),NA)})
    # Finish calculating mean water depth
    wat1.sum$XH2O_DEPTH_AA <- with(wat1.sum, round(XH2O_DEPTH_AA/NPLOTS, 2))
    # Merge metric data frames 
    wat1 <- merge(wat1.min, wat1.max, by = c('SAMPID','NPLOTS'))
    wat1 <- merge(wat1, wat1.sum, by = c('SAMPID','NPLOTS'))
    

  }
  # Melt output metrics
  varNames <- names(wat1)[!names(wat1) %in% c('SAMPID','NPLOTS')]
  wat1a <- reshape(wat1, idvar = c('SAMPID'), direction = 'long',
                   times = varNames, varying = varNames,
                   timevar = 'METRIC', v.names = 'RESULT')

  ## Fix Inf values to 0s
  wat1a$RESULT <- with(wat1a, ifelse(RESULT %in% c('Inf','-Inf'), 0, 
                                     ifelse(RESULT=='Inf_-Inf', '', RESULT)))
  wat1a$METRIC <- with(wat1a, as.character(METRIC))
  wat1a$NPLOTS <- NULL # drop variable

  # Keep only positive values from water parameters
  wdep2a.sub <- subset(wdep2a, as.numeric(RESULT)>0)
  # Sum water cover by type of vegetation
  wat2.sum <- aggregate(x = list(XCOV_H2O = wdep2a.sub$RESULT),
                        by = wdep2a.sub[c('SAMPID','PARAMETER','NPLOTS')],
                        FUN = function(x){sum(as.numeric(x))})
  # Count number of occurrences of each parameter
  wat2.freq <- aggregate(x = list(FREQ_H2O = wdep2a.sub$RESULT),
                         by = wdep2a.sub[c('SAMPID','PARAMETER','NPLOTS')],
                         FUN = function(x){length(as.numeric(x))})
  # Merge data frames
  wat2 <- merge(wat2.sum, wat2.freq, by = c('SAMPID','PARAMETER','NPLOTS'), all=TRUE)
  # Calculate water cover metrics by parameter 
  wat2$FREQ_H2O <- with(wat2, round(FREQ_H2O/NPLOTS*100, 2))
  wat2$XCOV_H2O <- with(wat2, round(XCOV_H2O/NPLOTS, 2))
  wat2$IMP_H2O <- with(wat2, round((FREQ_H2O + XCOV_H2O)/2, 2))
  wat2$NPLOTS <- NULL # drop NPLOTS
  # Melt output data frame
  wat2a <- reshape(wat2, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                   times = c('FREQ_H2O','XCOV_H2O','IMP_H2O'), 
                   varying = c('FREQ_H2O','XCOV_H2O','IMP_H2O'),
                   timevar = 'METRIC', v.names = 'RESULT')
  # Update METRIC by combining PARAMETER and METRIC values, except for TOTAL_WATER parameter
  wat2a$METRIC <- with(wat2a, ifelse(PARAMETER=='TOTAL_WATER', as.character(METRIC), 
                                     as.character(paste(METRIC, substring(PARAMETER, 7), sep='_'))))
  # Drop PARAMETER
  wat2a$PARAMETER <- NULL

  ## Now pull only PREDOMINANT_DEPTH that is not missing and >0 to calculate XH2O_DEPTH
  wdep1.pos <- subset(wdep1,PREDOMINANT_DEPTH>0 & !is.na(PREDOMINANT_DEPTH))
  # Calculate mean water depth based on PREDOMINANT_DEPTH
  # Differs from XH2O_DEPTH_AA because denominator is number of plots sampled in _AA and just 
  # number of values measured >0 for this version
  wat3 <- aggregate(x = list(RESULT = wdep1.pos$PREDOMINANT_DEPTH), 
                    by = wdep1.pos[c('SAMPID')], 
                    FUN = function(x){mean(x, na.rm = TRUE)})
  wat3$METRIC <- 'XH2O_DEPTH'
  # Combine with previously calculated metrics
  watMet <- rbind(wat1a, wat2a, wat3)
  # Cast metric data frame wide and drop prefix added by reshape() to variable names
  watMet1 <- reshape(watMet, idvar = 'SAMPID', direction = 'wide',
                     timevar = 'METRIC', v.names = 'RESULT')
  names(watMet1) <- gsub("RESULT\\.", "", names(watMet1)) 

  # Create empty data frame based on year of survey - more metrics in 2011 than 2016.
  if(survyear=='2011'){
    empty_wat <- data.frame(t(rep(NA,18)), stringsAsFactors=F)
    names(empty_wat) <- c("MIN_H2O_DEPTH","MAX_H2O_DEPTH","XH2O_DEPTH_AA","MIN_COV_H2O",
                          "MAX_COV_H2O","FREQ_H2O","FREQ_H2O_AQVEG","FREQ_H2O_EMERGVEG",
                          "FREQ_H2O_NOVEG","XCOV_H2O","XCOV_H2O_AQVEG","XCOV_H2O_EMERGVEG",
                          "XCOV_H2O_NOVEG","IMP_H2O","IMP_H2O_AQVEG","IMP_H2O_EMERGVEG",
                          "IMP_H2O_NOVEG","XH2O_DEPTH")
  }else{
    empty_wat <- data.frame(t(rep(NA,7)), stringsAsFactors=F)
    names(empty_wat) <- c("XH2O_DEPTH_AA","MIN_COV_H2O","MAX_COV_H2O","FREQ_H2O",
                          "XCOV_H2O","IMP_H2O","XH2O_DEPTH")
  }
  # Merge empty data frame with metric data frame, remove row with missing SAMPID
  watMet2 <- subset(merge(watMet1, empty_wat, all=TRUE), !is.na(SAMPID))
  # Merge with full list of samples to ensure all are included in output file
  watMet3 <- merge(allUIDs,watMet2, by='SAMPID',all.x=T)
  # Melt data frame and substitute missing RESULT values with 0
  varNames <- names(watMet3)[!names(watMet3) %in% c('SAMPID','NPLOTS')]
  watMet4 <- reshape(watMet3, idvar = 'SAMPID', direction = 'long',
                     times = varNames, varying = varNames,
                     timevar = 'METRIC', v.names = 'RESULT')
  
  watMet4$METRIC <- with(watMet4, as.character(METRIC))
  watMet4$RESULT <- with(watMet4, ifelse(is.na(RESULT), '0', RESULT))

  print("Done with water depth metrics")
  # Merge samples with data frame to get sampID variables back, then drop SAMPID
  # Also rename METRIC to PARAMETER
  watMetOut <- merge(samples, watMet4, by='SAMPID') 
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

calcBareGround_LitterMets <- function(dataIn, nPlot, sampID='UID', survyear='2011'){
  ## Now merge back with input df
  dataIn1 <- merge(dataIn, nPlot, by=sampID)

  # Create vector of all samples in dataset
  for(i in 1:length(sampID)){
    if(i==1) dataIn1$SAMPID <- dataIn1[,sampID[i]]
    else dataIn1$SAMPID <- paste(dataIn1$SAMPID, dataIn1[, sampID[i]], sep='.')
  }
  samples <- unique(subset(dataIn1, select=c(sampID, 'SAMPID')))
  
  allUIDs <- data.frame(SAMPID=unique(dataIn1$SAMPID), stringsAsFactors=F)
  
  # LITTER TYPES
  # Need to calculate the number of quadrats sampled using the NE and SW parameters
  # Only one set of parameters or the other will be in any dataset
  litter.sub <- subset(dataIn1,PARAMETER %in% c('LITTER_DEPTH_NE','LITTER_DEPTH_SW','DEPTH_NE','DEPTH_SW')) 
  # Calculate the number of quadrats (SW and NE)
  numQuads <- aggregate(x = list(NQUADS = litter.sub$RESULT),
                        by = litter.sub[c('SAMPID')], FUN = length)
  # Subset data for litter and bare ground parameters
  litter1 <- subset(dataIn1,PARAMETER %in% c('LITTER_THATCH','LITTER_FORB','LITTER_CONIFER',
                                             'LITTER_DECID','LITTER_BROADLEAF',
                                             'LITTER_DEPTH_SW','LITTER_DEPTH_NE',
                                             'TOTAL_LITTER','WD_FINE','WD_COARSE',
                                             'EXPOSED_SOIL','EXPOSED_GRAVEL',
                                             'EXPOSED_ROCK','DEPTH_SW','DEPTH_NE',
                                             'PREDOMINANT_LITTER'))
  # Merge litter subset with number of quadrats
  litter2 <- merge(litter1, numQuads, by='SAMPID')
  # Subset to needed variables
  litter2 <- subset(litter2, select = c('SAMPID', 'PLOT','PARAMETER','RESULT','NPLOTS','NQUADS'))
  # If survey year 2011, differs in parameter names and format of data from later years
  if(survyear=='2011'){
    littype <- subset(litter2, PARAMETER %in% c('LITTER_THATCH','LITTER_FORB','LITTER_CONIFER','LITTER_DECID',
                                               'LITTER_BROADLEAF','LITTER_NONE') & !is.na(RESULT))
    # Count number of litter types
    rr1 <- aggregate(x = list(N_LITTER_TYPE = littype$PARAMETER),
                     by = littype[c('SAMPID')],
                     FUN = function(x){length(unique(x))})
    # Merge this count back to litter data
    rr1 <- merge(littype, rr1, by = 'SAMPID')
    # Count plots by litter type
    rr2 <- aggregate(x = list(NUM = rr1$PLOT), by = rr1[c('SAMPID','PARAMETER','N_LITTER_TYPE')],
                     FUN = length)
    # Find dominant litter type 
    rr3 <- aggregate(x = list(MAXN = rr2$NUM), by = rr2[c('SAMPID')], FUN = max)
    rr3 <- merge(rr2, rr3, by = 'SAMPID')
    rr4 <- subset(rr3, MAXN==NUM)
    # Identify max litter types and combine if more than one with same frequency
    rr5 <- aggregate(x = list(LITTER_TYPE = rr4$PARAMETER), 
                     by = rr4[c('SAMPID','N_LITTER_TYPE')],
                     FUN = function(x){paste(x, collapse='_')})
    rr5$LITTER_TYPE <- with(rr5, gsub('LITTER_', '', LITTER_TYPE))
    
  }else{
    # Subset to predominant litter parameter where not missing
    littype <- subset(litter2, PARAMETER=='PREDOMINANT_LITTER' & !is.na(RESULT)) 
    # Count the total number of litter types
    rr1 <- aggregate(x = list(N_LITTER_TYPE = littype$RESULT),
                     by = littype[c('SAMPID')],
                     FUN = function(x){length(unique(x))})
    # Merge back with data
    rr1 <- merge(littype, rr1, by = 'SAMPID')
    # Count number of each litter type 
    rr2 <- aggregate(x = list(NUM = rr1$PLOT), by = rr1[c('SAMPID','RESULT','N_LITTER_TYPE')],
                     FUN = length)
    # Find max litter type
    rr3 <- aggregate(x = list(MAXN = rr2$NUM), by = rr2[c('SAMPID')], FUN = max)
    rr3 <- merge(rr2, rr3, by = 'SAMPID')
    rr4 <- subset(rr3, MAXN==NUM)
    # Combine litter types if more than one
    rr5 <- aggregate(x = list(LITTER_TYPE = rr4$RESULT), 
                     by = rr4[c('SAMPID','N_LITTER_TYPE')],
                     FUN = function(x){paste(x, collapse='_')})

  }
  # To determine median depth, we must account for any quadrats without depth recorded
  litdep <- subset(litter2,PARAMETER %in% c('LITTER_DEPTH_SW','LITTER_DEPTH_NE','DEPTH_NE','DEPTH_SW'))
  # Sum litter depth
  ss1.sum <- aggregate(x = list(XDEPTH_LITTER = litdep$RESULT), 
                       by = litdep[c('SAMPID','NQUADS')],
                       FUN = function(x){sum(as.numeric(x))})
  # Calculate number of depths
  ss1.length <- aggregate(x = list(NSAMP = litdep$RESULT), 
                          by = litdep[c('SAMPID','NQUADS')],
                          FUN = length)
  # Merge data frames
  ss1 <- merge(ss1.sum, ss1.length, by = c('SAMPID','NQUADS'))
  # Use the merged values to calculate the mean litter depth
  ss1$XDEPTH_LITTER <- with(ss1, round(XDEPTH_LITTER/NQUADS, 2))
  # calculate number of samples with litter depth info
  litdep1 <- aggregate(x = list(NSAMP = litdep$RESULT), 
                       by = litdep[c('SAMPID','NPLOTS','NQUADS')],
                       FUN = length)
  # Merge back with data
  litdep1 <- merge(litdep, litdep1, by = c('SAMPID','NPLOTS','NQUADS'))
  # Calculate median litter depth
  tt <- aggregate(x = list(MEDDEPTH_LITTER = litdep1$RESULT),
                  by = litdep1[c('SAMPID')], 
                  FUN = function(x){median(as.numeric(x))})
  # Melt litter type metrics data frame 
  varNames <- names(rr5)[!names(rr5) %in% 'SAMPID']
  rr5.long <- reshape(rr5, idvar = 'SAMPID', direction = 'long',
                      times = varNames, varying = varNames,
                      timevar = 'METRIC', v.names = 'RESULT')
  # Melt litter depth metrics, select subset of variables 
  ss1.long <- reshape(ss1, idvar = 'SAMPID', direction = 'long',
                      times = 'XDEPTH_LITTER', varying = 'XDEPTH_LITTER',
                      timevar = 'METRIC', v.names = 'RESULT') 
  ss1.long <- subset(ss1.long, select = c('SAMPID','METRIC','RESULT'))
  
  # Melt median litter depth data frame
  varNames <- names(tt)[!names(tt) %in% 'SAMPID']
  tt.long <- reshape(tt, idvar = 'SAMPID', direction = 'long',
                      times = varNames, varying = varNames,
                      timevar = 'METRIC', v.names = 'RESULT')
  # Combine datasets into one data frame
  loutdf <- rbind(rr5.long, ss1.long, tt.long)
  loutdf$METRIC <- with(loutdf, as.character(METRIC))

  # Cast output data frame wide and drop prefix from variable names added by reshape()
  loutdf1 <- reshape(loutdf, idvar = 'SAMPID', direction = 'wide',
                     timevar = 'METRIC', v.names = 'RESULT')
  names(loutdf1) <- gsub("RESULT\\.", "", names(loutdf1))

  # Create empty data frame with all expected metrics
  empty_lit <- data.frame(t(rep(NA,4)), stringsAsFactors=FALSE)
  names(empty_lit) <- c('N_LITTER_TYPE','LITTER_TYPE','XDEPTH_LITTER','MEDDEPTH_LITTER')

  # Merge empty data frame with output data frame, drop row with missing SAMPID
  loutdf2 <- subset(merge(loutdf1, empty_lit, all=TRUE), !is.na(SAMPID))
  # Merge full sample list with output data frame
  loutdf3 <- merge(allUIDs,loutdf2,by='SAMPID',all.x=T)
  # Melt data frame
  varNames <- names(loutdf3)[!names(loutdf3) %in% 'SAMPID']
  loutdf4 <- reshape(loutdf3, idvar = 'SAMPID', direction = 'long',
                     times = varNames, varying = varNames,
                     timevar = 'METRIC', v.names = 'RESULT')
  loutdf4$METRIC <- with(loutdf4, as.character(METRIC))

  print("Done with litter metrics")
  # Cast data frame wide and drop variable name prefix added by reshape()
  litterOut <- reshape(loutdf4, idvar = 'SAMPID', direction = 'wide',
                       timevar = 'METRIC', v.names = 'RESULT') 
  names(litterOut) <- gsub("RESULT\\.", "", names(litterOut))

  # BARE GROUND
  # Subset data to only parameters related to woody debris and bare ground, total litter
  # Keep only positive RESULT values
  bgrd <- subset(dataIn1,PARAMETER %in% c('EXPOSED_SOIL','EXPOSED_GRAVEL','EXPOSED_ROCK',
                                          'WD_FINE','WD_COARSE','TOTAL_LITTER') & 
                   RESULT!=0 & !is.na(RESULT), select=c('SAMPID','PLOT','PARAMETER','RESULT','NPLOTS'))
  # Need to create values for new PARAMETER='BAREGD' based on occurrence of either EXPOSED_SOIL, 
  # EXPOSED_GRAVEL, or EXPOSED_ROCK at site
  bgrd.sub <- subset(bgrd, PARAMETER %in% c('EXPOSED_SOIL','EXPOSED_GRAVEL','EXPOSED_ROCK'))
  # Sum bare ground cover
  bgrd1 <- aggregate(x = list(RESULT = bgrd.sub$RESULT), 
                     by = bgrd.sub[c('SAMPID','PLOT','NPLOTS')],
                     FUN = function(x){as.character(sum(as.numeric(x)))})
  bgrd1$PARAMETER <- 'BAREGD'
  # Add these data back into full dataset
  bgrdIn <- rbind(bgrd, bgrd1)
  
  # Need to fill in zeros if plot sampled and variable is zero
  # Cast data frame, then drop prefix added to variable names by reshape()
  bgrdIn1.wide <- reshape(bgrdIn, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'wide',
                     timevar = 'PARAMETER', v.names = 'RESULT')
  names(bgrdIn1.wide) <- gsub("RESULT\\.", "", names(bgrdIn1.wide))
  # Melt data frame, fill in missing RESULT values with 0
  varNames <- names(bgrdIn1.wide)[!names(bgrdIn1.wide) %in% c('SAMPID','PLOT','NPLOTS')]
  bgrdIn1 <- reshape(bgrdIn1.wide, idvar = c('SAMPID','PLOT','NPLOTS'), direction = 'long',
                     times = varNames, varying = varNames,
                     timevar = 'PARAMETER', v.names = 'RESULT')
  bgrdIn1$RESULT <- with(bgrdIn1, ifelse(is.na(RESULT), 0, as.numeric(RESULT)))
  # Keep only values from above where RESULT is >0
  bgrdIn1.sub <- subset(bgrdIn1,RESULT!=0)
  # Count the number of occurrences of each type of bare ground
  bgindf1.freq <- aggregate(x = list(FREQ = bgrdIn1.sub$RESULT),
                            by = bgrdIn1.sub[c('SAMPID','PARAMETER','NPLOTS')],
                            FUN = length)
  # Sum the cover of each type of bare ground  
  bgindf1.sum <- aggregate(x = list(XCOV = bgrdIn1.sub$RESULT),
                           by = bgrdIn1.sub[c('SAMPID','PARAMETER','NPLOTS')],
                           FUN = function(x){sum(as.numeric(x))})
  # Merge these calculated values, then calculate the metrics based on them
  bgindf1 <- merge(bgindf1.sum, bgindf1.freq, by = c('SAMPID','PARAMETER','NPLOTS'), all=TRUE)
  bgindf1$FREQ <- with(bgindf1, round(FREQ/NPLOTS*100, 2))
  bgindf1$XCOV <- with(bgindf1, round(XCOV/NPLOTS, 2))
  bgindf1$IMP <- with(bgindf1, round((FREQ + XCOV)/2, 2))
 
  # Melt output data frame
  bgoutdf <- reshape(bgindf1, idvar = c('SAMPID','PARAMETER'), direction = 'long',
                     times = c('XCOV','FREQ','IMP'), varying = c('XCOV','FREQ','IMP'),
                     timevar = 'variable', v.names = 'RESULT')
  # Create parameter name based on combinations of variable and PARAMETER
  bgoutdf$PARAMETER <- with(bgoutdf, ifelse(PARAMETER=='TOTAL_LITTER',
                                      paste(variable,'LITTER', sep='_'), 
                                      paste(variable, PARAMETER, sep='_')))
  bgoutdf <- bgoutdf[,c('SAMPID','PARAMETER','RESULT')]
  # Cast data frame and drop prefix to variable name added by reshape()
  bgoutdf1.wide <- reshape(bgoutdf, idvar = 'SAMPID', direction = 'wide',
                           timevar = 'PARAMETER', v.names = 'RESULT')
  names(bgoutdf1.wide) <- gsub("RESULT\\.", "", names(bgoutdf1.wide))
  # Melt data frame
  varNames <- names(bgoutdf1.wide)[!names(bgoutdf1.wide) %in% c('SAMPID')]
  bgoutdf1 <- reshape(bgoutdf1.wide, idvar = c('SAMPID'), direction = 'long',
                      times = varNames, varying = varNames,
                      timevar = 'METRIC', v.names = 'RESULT')
  # Fill in zeroes where missing
  bgoutdf1$RESULT <- with(bgoutdf1, ifelse(is.na(RESULT), 0, RESULT))
  bgoutdf1$METRIC <- with(bgoutdf1, as.character(METRIC))
  # Cast data wide and drop variable name prefix added by reshape()
  bgoutdf2 <- reshape(bgoutdf1, idvar = 'SAMPID', direction = 'wide',
                      timevar = 'METRIC', v.names = 'RESULT')
  names(bgoutdf2) <- gsub("RESULT\\.", "", names(bgoutdf2))
  # Create empty data frame
  empty_bg <- data.frame(t(rep(NA,21)),stringsAsFactors=FALSE)
  names(empty_bg) <- c("FREQ_BAREGD","FREQ_EXPOSED_GRAVEL","FREQ_EXPOSED_ROCK",
                       "FREQ_EXPOSED_SOIL","FREQ_LITTER","FREQ_WD_COARSE",
                       "FREQ_WD_FINE","IMP_BAREGD","IMP_EXPOSED_GRAVEL",
                       "IMP_EXPOSED_ROCK","IMP_EXPOSED_SOIL","IMP_LITTER",
                       "IMP_WD_COARSE","IMP_WD_FINE","XCOV_BAREGD",
                       "XCOV_EXPOSED_GRAVEL","XCOV_EXPOSED_ROCK","XCOV_EXPOSED_SOIL",
                       "XCOV_LITTER","XCOV_WD_COARSE","XCOV_WD_FINE")
  # Merge empty data frame with output data frame and drop row with missing SAMPID
  bgoutdf3 <- subset(merge(bgoutdf2, empty_bg, all=TRUE), !is.na(SAMPID))
  # Merge full sample list with output data frame
  bgoutdf4 <- merge(allUIDs, bgoutdf3, by='SAMPID', all.x=T)
  # Melt resulting data frame and fill in 0 for RESULT where missing
  varNames <- names(bgoutdf4)[!names(bgoutdf4) %in% c('SAMPID')]
  bgoutdf5 <- reshape(bgoutdf4, idvar = 'SAMPID', direction = 'long',
                      times = varNames, varying = varNames,
                      timevar = 'METRIC', v.names = 'RESULT')
  bgoutdf5$METRIC <- with(bgoutdf5, as.character(METRIC))
  bgoutdf5$RESULT <- with(bgoutdf5, ifelse(is.na(RESULT), 0, RESULT))

  print("Done with bare ground metrics")
  # Cast data frame wide and drop variable name prefix added by reshape()
  bgrdOut <- reshape(bgoutdf5, idvar = 'SAMPID', direction = 'wide',
                     timevar = 'METRIC', v.names = 'RESULT')
  names(bgrdOut) <- gsub("RESULT\\.", "", names(bgrdOut))

  # Now combine bare ground and litter metrics to perform a final check
  combOut <- merge(litterOut, bgrdOut, by='SAMPID', all=T)
  
  # Melt data frame and fill in 0 where missing, except if metric is LITTER_TYPE 
  varNames <- names(combOut)[!names(combOut) %in% c('SAMPID')]
  combOut.1 <- reshape(combOut, idvar = 'SAMPID', direction = 'long',
                       times = varNames, varying = varNames,
                       timevar = 'variable', v.names = 'value')
  combOut.1$value <- with(combOut.1, ifelse(is.na(value) & variable!='LITTER_TYPE', 0, value))
  # Need to determine whether missing LITTER_TYPE value should be ABSENT or NONE
  # Cast output data frame and drop variable name prefix added by reshape()
  combOut.wide <- reshape(combOut.1, idvar = 'SAMPID', direction = 'wide',
                          timevar = 'variable', v.names = 'value')
  names(combOut.wide) <- gsub("value\\.", "", names(combOut.wide))
  # If LITTER_TYPE is missing and mean litter depth, mean litter cover, and litter frequency are all 0, 
  # set LITTER_TYPE to 'ABSENT', otherwise set missing LITTER_TYPE to 'NONE'
  combOut.wide$LITTER_TYPE <- with(combOut.wide, ifelse(is.na(LITTER_TYPE) & XDEPTH_LITTER==0 & XCOV_LITTER==0 & FREQ_LITTER==0, 'ABSENT', ifelse(is.na(LITTER_TYPE),'NONE',LITTER_TYPE)))
  # Melt output data frame
  varNames <- names(combOut.wide)[!names(combOut.wide) %in% c('SAMPID')]
  combOut.long <- reshape(combOut.wide, idvar = 'SAMPID', direction = 'long',
                       times = varNames, varying = varNames,
                       timevar = 'PARAMETER', v.names = 'RESULT')
  # Merge data with samples to bring back in sampID variables, then select only needed variables
  finOut <- merge(samples, combOut.long, by='SAMPID') 
  finOut <- finOut[c(sampID, 'PARAMETER', 'RESULT')]
    
  return(finOut)
}

