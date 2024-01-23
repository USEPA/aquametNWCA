# Calculate only duration metrics (both with and without native status, if available)
#' @export
#' 
#' @title Calculate vascular plant duration metrics
#' 
#' @description This function calculates all duration metrics, with
#' additional metrics if NWCA_NATSTAT variable is included
#' in input data frame.
#' 
#' @details Both DURATION and NWCA_NATSTAT variables are recoded to fewer
#' categories. Taxa with 'UND' as native status are excluded.
#' 
#' @param vascIn   Data frame containing cover data summarized by
#' sampID variables and TAXON, with the following fields:
#'  \itemize{
#'     \item sampID: Variable(s) in the argument \emph{sampID}
#'
#'     \item TAXON: Taxon name
#'
#'     \item DURATION: USDA PLANTS DURATION variable, with valid
#'     values "ANNUAL", "BIENNIAL", "PERENNIAL", combinations
#'     of these values, or blank
#'     
#'     \item DUR_ALT (optional): Combinations of DURATION into
#'     ANNUAL, PERENNIAL, ANN_BIEN, and ANN_PEREN, as used in 
#'     NWCA. Created in code if not supplied by user.
#'
#'     \item XABCOV: Mean percent cover of taxon across plots
#'
#'     \item TOTN: Number of taxa in sample
#'
#'     \item sXRCOV: proportion of summed cover across all taxa
#'     (XTOTABCOV) represented by taxon in sample
#'
#'     \item NWCA_NATSTAT (Optional): Native status variable with categories
#'     of 'NAT', 'ADV', 'CRYP', 'INTR', 'UND'
#'     }
#' @param sampID  A character vector containing the name(s) of variable(s)
#'   necessary to identify unique samples, 'UID' by default
#' 
#' @return Data frame containing the \emph{sampID} variables, PARAMETER, RESULT,
#'   where values of PARAMETER consist of the metric name concatenated with
#'   trait value (represented as TRAITNM below): 
#'  \itemize{ \item N_TRAITNM:
#'   Number of taxa with trait
#'   
#'   \item PCTN_TRAITNM: Number of taxa with trait as percentage of \emph{TOTN}
#'   
#'   \item XABCOV_TRAITNM: Sum of \emph{XABCOV} values across taxa with trait
#'   
#'   \item XRCOV_TRAITNM: Sum of \emph{sXRCOV} values across taxa with trait
#' }
#' If NWCA_NATSTAT is in the input data frame, the same set except
#' 'N_' metrics is calculated for native species and
#' alien + cryptogenic species, with metric names suffixes of
#' _NAT and _AC, respectively.
#' 
#' @author Karen Blocksom \email{Blocksom.karen@epa.gov}
#' 
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC.
#' 
#' @examples
#'  head(VascPlantEx)
#'  exPlant <- prepareData(VascPlantEx, taxon_name = 'USDA_NAME', 
#'  inTaxa = taxaNWCA, inNat = ccNatNWCA, inCVal = ccNatNWCA, 
#'  inWIS = wisNWCA, cValReg='STATE')
#'
#'  durEx <- calcDuration(exPlant$byUIDspp)
#'  head(durEx)
#'  unique(durEx$PARAMETER)

calcDuration <- function(vascIn, sampID='UID'){
  vascIn <- as.data.frame(vascIn) # Do this in case read in as a tibble or data.table, which might cause problems
  # First check for necessary variables
  necNames <- c(sampID, 'TAXON', 'DURATION', 'XABCOV', 'TOTN', 'sXRCOV')
  msgNames <- necNames[necNames %nin% names(vascIn)]
  if(length(msgNames)>0){
    print(paste("Missing key variables for metric calculation: ", paste(msgNames,collapse=','),
                ". Try prepareData() function to create necessary input variables.",sep=''))
    return(NULL)
  }
  
  vascIn <- subset(vascIn,!is.na(DURATION) & toupper(DURATION)!='UND')
  # From DURATION, create DUR_ALT variable
  if('DUR_ALT' %in% names(vascIn)){
    vascIn.1 <- vascIn
  }else{
    vascIn.1 <- vascIn
    vascIn.1$DUR_ALT <- vascIn.1$DURATION
    vascIn.1$DUR_ALT[vascIn.1$DURATION %in% c('ANNUAL, BIENNIAL',
                                              'BIENNIAL',
                                              'BIENNIAL, AN')] <- 'ANN_BIEN'
    vascIn.1$DUR_ALT[vascIn.1$DURATION %in% c('ANNUAL, BIENNIAL, PERENNIAL',
                                              'ANNUAL, PERENNIAL',
                                              'PERENNIAL, ANNUAL', 
                                              'BIENNIAL, PERENNIAL',
                                              'ANNUAL, PERENNIAL, BIENNIAL', 
                                              'BIENNIAL, PERENNIAL, AN')] <- 'ANN_PEREN'
    vascIn.1$DUR_ALT[vascIn.1$DURATION %in% c('PERENNIAL, ANNUAL, BIENNIAL',
                                              'PERENNIAL, AN',
                                              'PERENNIAL, BIENNIAL',
                                              'PERENNIAL, BIENNIAL, ANNUAL',
                                              'PERENNIAL, BIENNIAL, AN')] <- 'PERENNIAL'
    
  }
  durOut <- int.calcTraits_MultCat(vascIn.1, 'DUR_ALT', sampID)
  
  empty_base <- data.frame(t(rep(NA,16)), stringsAsFactors=F)
  names(empty_base) <- c("N_ANN_BIEN","N_ANN_PEREN","N_ANNUAL","N_PERENNIAL", "PCTN_ANN_BIEN","PCTN_ANNUAL"
                         ,"PCTN_ANN_PEREN","PCTN_PERENNIAL","XABCOV_ANN_BIEN","XABCOV_ANN_PEREN"
                         ,"XABCOV_ANNUAL","XABCOV_PERENNIAL","XRCOV_ANN_BIEN","XRCOV_ANN_PEREN"
                         ,"XRCOV_ANNUAL","XRCOV_PERENNIAL")
  
  empty_base.nat <- data.frame(t(rep(NA,32)), stringsAsFactors=F)
  names(empty_base.nat) <- c("N_ANN_BIEN_AC","N_ANN_BIEN_NAT","N_ANN_PEREN_AC","N_ANN_PEREN_NAT"
                             ,"N_ANNUAL_AC","N_ANNUAL_NAT","N_PERENNIAL_AC","N_PERENNIAL_NAT"
                             ,"PCTN_ANN_BIEN_AC","PCTN_ANN_BIEN_NAT","PCTN_ANN_PEREN_AC"
                             ,"PCTN_ANN_PEREN_NAT","PCTN_ANNUAL_AC","PCTN_ANNUAL_NAT"
                             ,"PCTN_PERENNIAL_AC","PCTN_PERENNIAL_NAT","XABCOV_ANN_BIEN_AC"
                             ,"XABCOV_ANN_BIEN_NAT","XABCOV_ANN_PEREN_AC","XABCOV_ANN_PEREN_NAT"
                             ,"XABCOV_ANNUAL_AC","XABCOV_ANNUAL_NAT","XABCOV_PERENNIAL_AC"
                             ,"XABCOV_PERENNIAL_NAT","XRCOV_ANN_BIEN_AC","XRCOV_ANN_BIEN_NAT"
                             ,"XRCOV_ANN_PEREN_AC","XRCOV_ANN_PEREN_NAT","XRCOV_ANNUAL_AC"
                             ,"XRCOV_ANNUAL_NAT","XRCOV_PERENNIAL_AC","XRCOV_PERENNIAL_NAT")
  
  if('NWCA_NATSTAT' %in% names(vascIn.1)){
    # Assign native status values to grouped categories
    vascIn.2 <- vascIn.1
    vascIn.2$ALIEN <- with(vascIn.2, ifelse(NWCA_NATSTAT %in% c('INTR','ADV'), 1, 0))
    vascIn.2$NATSTAT_ALT <- with(vascIn.2, ifelse(NWCA_NATSTAT %in% c('INTR','ADV'), 'ALIEN', NWCA_NATSTAT))
    vascIn.2$AC <- with(vascIn.2, ifelse(NWCA_NATSTAT %in% c('INTR','ADV','CRYP'), 1, 0))
    
    vascIn.2$ANNUAL_NAT <- with(vascIn.2, ifelse(DUR_ALT=='ANNUAL' & NATSTAT_ALT=='NAT', 1, 0))
    vascIn.2$ANNUAL_AC <- with(vascIn.2, ifelse(DUR_ALT=='ANNUAL' & AC==1, 1, 0))
    vascIn.2$ANN_BIEN_NAT <- with(vascIn.2, ifelse(DUR_ALT=='ANN_BIEN' & NATSTAT_ALT=='NAT', 1, 0))
    vascIn.2$ANN_BIEN_AC <- with(vascIn.2, ifelse(DUR_ALT=='ANN_BIEN' & AC==1, 1, 0))
    vascIn.2$ANN_PEREN_NAT <- with(vascIn.2, ifelse(DUR_ALT=='ANN_PEREN' & NATSTAT_ALT=='NAT', 1, 0))
    vascIn.2$ANN_PEREN_AC <- with(vascIn.2, ifelse(DUR_ALT=='ANN_PEREN' & AC==1, 1, 0))
    vascIn.2$PERENNIAL_NAT <- with(vascIn.2, ifelse(DUR_ALT=='PERENNIAL' & NATSTAT_ALT=='NAT', 1, 0))
    vascIn.2$PERENNIAL_AC <- with(vascIn.2, ifelse(DUR_ALT=='PERENNIAL' & AC==1, 1, 0))
    
    
    multTraits <- int.combTraits(vascIn.2, c('ANNUAL_NAT','ANNUAL_AC','ANN_BIEN_NAT','ANN_BIEN_AC','ANN_PEREN_NAT'
                                             ,'ANN_PEREN_AC','PERENNIAL_NAT','PERENNIAL_AC'), sampID)
    
    durOut <- rbind(durOut, multTraits)
    
    empty_base <- cbind(empty_base, empty_base.nat)
  }
  
  durOut.1a <- reshape(durOut, idvar = c(sampID), direction = 'wide',
                       timevar = 'PARAMETER', v.names = 'RESULT')
  
  names(durOut.1a) <- gsub("RESULT\\.", "", names(durOut.1a))
  durOut.1a <- merge(durOut.1a, empty_base, all=TRUE)
  
  durOut.1b <- reshape(durOut.1a, idvar = sampID, direction = 'long',
                       varying= names(durOut.1a)[!names(durOut.1a) %in% c(sampID)],
                       timevar = 'PARAMETER', v.names = 'RESULT',
                       times = names(durOut.1a)[!names(durOut.1a) %in% c(sampID)])
  durOut.1b <- subset(durOut.1b, !is.na(eval(as.name(sampID[1]))))
  durOut.1b$RESULT <- with(durOut.1b, ifelse(is.na(RESULT), 0, RESULT))
  durOut.1b$PARAMETER <- with(durOut.1b, as.character(PARAMETER))
  
  durOut.1 <- subset(durOut.1b, PARAMETER %in% names(empty_base))
  
  
  return(durOut.1)
  
}

