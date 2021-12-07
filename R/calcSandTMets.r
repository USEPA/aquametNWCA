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


