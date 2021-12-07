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


