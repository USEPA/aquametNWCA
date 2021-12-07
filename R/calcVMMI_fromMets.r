#' @export
#' 
#' @title Calculate 2011 Vegetation MMI and assign condition class
#' 
#' @description This function calculates the NWCA 2011 Vegetation
#' Multimetric Index (VMMI) from metric inputs. If the appropriate
#' variable(s) describing site ecoregion and wetland group type
#' according to NWCA 2011 are included in the input data frame,
#' condition class (Good/Fair/Poor) will also be assigned.
#' 
#' @param metsIn Data frame containing, at a minimum: \itemize{ \item sampID:
#'   Variables identified in \emph{sampID} argument
#'   
#'   \item FQAI_ALL: Floristic Quality Assessment Index based on all species
#'   
#'   \item N_TOL: Number of tolerant species
#'   
#'   \item RIMP_NATSPP: Relative importance of native species
#'   
#'   \item XRCOV_MONOCOTS_NAT: Mean relative cover of native monocots } 
#'   Optional: \itemize{ 
#'   \item NWCA_ECO4: Aggregated ecoregion abbreviations as
#'   used in NWCA 2011. Valid values include: CPL (Coastal Plains), EMU (Eastern
#'   Mountains and Uplands), IPL (Interior Plains), and W (West).
#'   
#'   \item NWCA_WET_GRP: Wetland type groups as used in NWCA 2011. Valid values
#'   include: EH (Estuarine Herbaceous), EW (Estuarine Woody), PRLH (Palustrine,
#'   Riverine, and Lacustrine Herbaceous), PRLW (Palustrine, Riverine, and
#'   Lacustrine Woody).
#'   
#'   \item ECO_X_WETGRP: Combinations of NWCA_ECO4 and NWCA_WET_GRP. Valid 
#'   values include ALL-EH, ALL-EW, CPL-PRLH, CPL-PRLW, EMU-PRLH, EMU-PRLW,
#'   IPL-PRLH, IPL-PRLW, W-PRLH, W-PRLW } To assign condition class, either
#'   NWCA_ECO4 and NWCA_WET_GRP need to be included OR ECO_X_WETGRP needs to be
#'   included in the input data frame.
#' @param sampID A character vector containing the name(s) of variable(s)
#'   necessary to identify unique samples, 'UID' by default
#'   
#' @return  Data frame containing: \itemize{ \item \emph{sampID} Variable(s)
#'   found in the argument sampID
#'   
#'   \item FQAI_ALL_SC: Scored FQAI_ALL, on 0-10 scale
#'   
#'   \item N_TOL_SC: Scored N_TOL, on 0-10 scale
#'   
#'   \item RIMP_NATSPP_SC: Scored RIMP_NATSPP, on 0-10 scale
#'   
#'   \item XRCOV_MONOCOTS_NAT_SC: Scored XRCOV_MONOCOTS_NAT, on 0-10 scale
#'   
#'   \item VMMI: NWCA 2011 Vegetation Multimetric Index on 100-point scale } If
#'   ECO_X_WETGRP is provided or can be determined from NWCA_ECO4 and 
#'   NWCA_WET_GRP, ECO_X_WETGRP and the following is also provided as output: 
#'   \itemize{ \item VEGCOND: Vegetation condition class, based on VMMI score
#'   and site location, as assigned for NWCA 2011. Valid values are Good, Fair,
#'   and Poor. }
#'   
#' @references US Environmental Protection Agency. 2016. National Wetland
#'   Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#'   Environmental Protection Agency, Washington, DC.
#'   
#' @author Karen Blocksom  \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(vmmiMetEx)
#'
#' # Test using ECO_X_WETGRP variable in data frame
#' check.1 <- calcVMMI_fromMets(vmmiMetEx)
#'
#' # Now drop ECO_X_WETGRP and use NWCA_ECO4 and NWCA_WET_GRP
#' check.2 <- calcVMMI_fromMets(subset(vmmiMetEx,select=-ECO_X_WETGRP))
#'
#' # Now run without ECO_X_WETGRP, NWCA_ECO4, and NWCA_WET_GRP
#' check.3 <- calcVMMI_fromMets(subset(vmmiMetEx,select=c('UID','FQAI_ALL'
#'                           ,'N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')))


calcVMMI_fromMets <- function(metsIn, sampID='UID'){
  # Look for UID, metrics in input data frame
  necVars <- c(sampID,'FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')
  # Alert user if any necessary variables are missing
  if(sum(necVars %in% names(metsIn))<5){
    msgVars <- necVars[necVars %nin% names(metsIn)]
    print(paste(paste(msgVars, collapse=', '),"not found in input data frame. Cannot calculate VMMI without them.",sep=' '))
  }

  # Look for region and wetland type variables - cannot assign condition without them
  # Can be individual variables that will be used to combine into ECO_X_WETGRP or ECO_X_WETGRP
  # Currently, we require the variable names to match those from NWCA 2011 to ensure they 
  # contain the correct values
  if(any(c('NWCA_ECO4','NWCA_WET_GRP') %nin% names(metsIn)) & 'ECO_X_WETGRP' %nin% names(metsIn)){
    print("Warning: Must include variables NWCA_WET_GRP and NWCA_ECO4 OR ECO_X_WETGRP to determine condition class. This variable is missing, and condition class will not be assigned, but VMMI will be calculated.")
  }

  # Create ECO_X_WETGRP if not present but NWCA_ECO4 and NWCA_WETGRP are both provided
  if('ECO_X_WETGRP' %nin% names(metsIn) & sum(c('NWCA_ECO4','NWCA_WET_GRP') %in% names(metsIn))==2){
    metsIn$ECO_X_WETGRP <- with(metsIn, ifelse(NWCA_WET_GRP %in% c('EH','EW'), 
                                               paste('ALL', NWCA_WET_GRP, sep='-')
                                               ,paste(NWCA_ECO4, NWCA_WET_GRP, sep='-')))
  }

  # Identify key variables in input dataset
  if('ECO_X_WETGRP' %in% names(metsIn)){
    keyVars <- c(sampID, 'ECO_X_WETGRP')
  }else{
    keyVars <- sampID
  }

  # Make sure all metrics are in numeric format for scoring
  necMets <- c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')
  metsIn[,necMets] <- lapply(metsIn[,necMets], as.numeric)

  # Melt data frame into long format
  metsIn.long <- reshape(metsIn, idvar = keyVars, direction = 'long',
                         varying= c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT'),
                         timevar = 'PARAMETER', v.names = 'RESULT',
                         times = c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT'))
  
  # Set metric scoring thresholds
  metTholds <- data.frame(PARAMETER=c('FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT'), 
                          CEILING=c(38.59,40,100,100),
                          FLOOR=c(6.94,0,44.34,0.06), 
                          DIRECTION=c('POS','NEG','POS','POS'), 
                          stringsAsFactors=FALSE)
  # Merge metric scoring thresholds with metrics in long format
  vMet <- merge(metsIn.long, metTholds, by='PARAMETER')

  # Now apply the function that calculates metric scores by interpolating values between 0 and 10
  scoreMet<-function(dir, x, floor, ceiling){
    if(dir=='POS'){ # if positive metric, interpolate with first approach
      zz<-round(approx(x=c(floor, ceiling), y=c(0,10), xout=x, method='linear', yleft=0, yright=10)$y, 2)
    } else { # if negative metric, interpolate score with this approach
      zz<-round(approx(x=c(floor, ceiling), y=c(10,0), xout=x, method='linear', yleft=10, yright=0)$y, 2)
    }

  }

  # Calculate scores and add scored version of  metric (METRIC_SC) to data frame. SC = rescaled
  # metric score that is used in MMI calculations
  # Keep only variables that are needed
  scored.mets <- vMet[,c(keyVars, 'PARAMETER')]
  # Use scoreMet() function to calculate RESULT value for each metric
  scored.mets$RESULT <- with(scored.mets, with(vMet, mapply(scoreMet, DIRECTION, RESULT, FLOOR, CEILING)))
  # Add _SC to end of each metric name as name of scored metric
  scored.mets$PARAMETER <- with(scored.mets, paste(as.character(PARAMETER), 'SC', sep='_'))

  # Now that we have scored metrics, we can calculate MMI scores and merge with MMI 
  # thresholds to determine condition
  # Sum metrics into MMI
  mmi.1 <- aggregate(x = list(VMMI = scored.mets$RESULT), by = scored.mets[keyVars],
                     FUN = sum)
  # Rescale VMMI to 100-point scale, rounding to 1 digit
  mmi.1$VMMI <- with(mmi.1, round(VMMI*(10/4), 1))
  # Melt data frame into long format
  mmi <- reshape(mmi.1, idvar = keyVars, direction = 'long',
                 varying = 'VMMI', timevar = 'PARAMETER', v.names='RESULT',
                 times = 'VMMI')
  # Combine metric scores and VMMI scores
  mmiOut <- rbind(scored.mets, mmi)

  # Set VMMI condition class thresholds
  mmiTholds <- data.frame(ECO_X_WETGRP=c('CPL-PRLH','CPL-PRLW','ALL-EH','ALL-EW','EMU-PRLH','EMU-PRLW'
                                         ,'IPL-PRLH','IPL-PRLW','W-PRLH','W-PRLW'),
                          p05=c(57.3,52.8,65.0,56.0,41.6,55.8,25.3,40.3,30.0,47.9),
                          p25=c(62.5,58.6,74.1,62.9,63.0,60.5,36.2,49.4,57.4,54.4),
                          stringsAsFactors=FALSE)

  # If ECO_X_WETGRP exists in mmi data frame, assign condition
  if('ECO_X_WETGRP' %in% names(mmi)){
    mmi.1 <- merge(mmi, mmiTholds, by='ECO_X_WETGRP')

    cond <- mmi.1
    cond$VEGCOND <- with(mmi.1, ifelse(RESULT>=p25, 'GOOD', ifelse(RESULT>=p05, 'FAIR', 'POOR')))
    # Melt condition data frame to combine with metric and VMMI scores
    cond.long <- reshape(cond, idvar = keyVars, direction = 'long',
                         varying = 'VEGCOND', times = 'VEGCOND',
                         timevar = 'PARAMETER', v.names = 'RESULT')
    cond.long <- subset(cond.long, select = c(keyVars, 'PARAMETER', 'RESULT'))
    # Combine with metric and VMMI scores
    mmiOut <- rbind(mmiOut, cond.long)
  }
  # Cast output data frame into wide format, then remove prefixes 
  # added by reshape() from variable names
  mmiOut.wide <- reshape(mmiOut, idvar = c(keyVars), direction = 'wide',
                         timevar = 'PARAMETER', v.names = 'RESULT')
  names(mmiOut.wide) <- gsub("RESULT\\.", "", names(mmiOut.wide))
  
return(mmiOut.wide)
}