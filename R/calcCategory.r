# Category metrics (with and without native status, if available)
#' @export
#' 
#' @title Calculate plant category metrics
#' 
#' @description This function calculates all category-based metrics, including
#' versions with only native species if the variable NWCA_NATSTAT
#' is included in the input data frame.
#' 
#' @details Both CATEGORY and NWCA_NATSTAT variables are recoded to fewer
#' categories. Taxa with 'UND' as native status are excluded.
#' 
#' @param vascIn   Data frame containing cover data summarized by
#' UID and TAXON, with the following fields:
#' \itemize{
#'     \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#'     \item TAXON: Taxon name
#'
#'     \item CATEGORY: USDA PLANTS category variable, with valid values
#'     of DICOT, FERN, GYMNOSPERM, HORSETAIL, LIVERWORT,
#'     LYCOPOD, MONOCOT, QUILLWORT, or blank
#'
#'     \item XABCOV: Mean percent cover of taxon across plots
#'
#'     \item TOTN: Number of taxa in sample
#'
#'     \item sXRCOV: proportion of summed cover across all taxa
#'     (XTOTABCOV) represented by taxon in sample
#'
#'     \item NWCA_NATSTAT (optional): Native status variable with categories
#'     of 'NAT', 'ADV', 'CRYP', 'INTR', 'UND'
#'    }
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#' 
#' @return Data frame containing \emph{sampID} variables, PARAMETER, RESULT,
#'   where values of PARAMETER consist of the metric name concatenated with
#'   trait value (represented as TRAITNM below):
#' \itemize{
#' \item N_TRAITNM: Number of taxa with trait
#'
#' \item PCTN_TRAITNM: Number of taxa with trait as percentage of \emph{TOTN}
#'
#' \item XABCOV_TRAITNM: Sum of \emph{XABCOV} values across taxa with trait
#'
#' \item XRCOV_TRAITNM: Sum of \emph{sXRCOV} values across taxa with trait
#' }
#' For metrics using native status, the metric name has a suffix
#' of NAT, ALIEN, AC, INTR, or CRYP.
#' 
#' @author Karen Blocksom \email{Blocksom.karen@epa.gov}
#' 
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC.
#' 
#' @examples
#' head(VascPlantEx)
#'  exPlant <- prepareData(VascPlantEx, taxon_name = 'USDA_NAME', 
#'  inTaxa = taxaNWCA, inNat = ccNatNWCA, inCVal = ccNatNWCA, 
#'  inWIS = wisNWCA, cValReg='STATE')
#'
#' catEx <- calcCategory(exPlant$byUIDspp)
#'
#' head(catEx)
#' unique(catEx$PARAMETER)

calcCategory <- function(vascIn, sampID='UID'){
  vascIn <- as.data.frame(vascIn) # Do this in case read in as a tibble or data.table, which might cause problems
  # First check for necessary variables
  necNames <- c(sampID, 'TAXON', 'CATEGORY', 'XABCOV', 'TOTN', 'sXRCOV')
  msgNames <- necNames[necNames %nin% names(vascIn)]
  if(length(msgNames)>0){
    print(paste("Missing key variables for metric calculation: ", paste(msgNames, collapse=','),
                ". Try prepareData() function to create necessary input variables.", sep=''))
    return(NULL)
  }
  
  vascIn.1 <- subset(vascIn, !is.na(CATEGORY) & toupper(CATEGORY)!='UND' & CATEGORY!='')
  catOut <- int.calcTraits_MultCat(vascIn.1, 'CATEGORY', sampID)
  
  empty_base <- data.frame(t(rep(NA,32)), stringsAsFactors=F)
  names(empty_base) <- c("N_DICOT","N_FERN","N_GYMNOSPERM","N_MONOCOT","N_LYCOPOD","N_HORSETAIL"
                         ,"N_QUILLWORT","N_LIVERWORT","PCTN_DICOT"
                         ,"PCTN_FERN","PCTN_GYMNOSPERM","PCTN_MONOCOT","PCTN_LYCOPOD","PCTN_HORSETAIL"
                         ,"PCTN_QUILLWORT","PCTN_LIVERWORT","XABCOV_DICOT"
                         ,"XABCOV_FERN","XABCOV_GYMNOSPERM","XABCOV_MONOCOT","XABCOV_LYCOPOD","XABCOV_HORSETAIL"
                         ,"XABCOV_QUILLWORT","XABCOV_LIVERWORT"
                         ,"XRCOV_DICOT","XRCOV_FERN","XRCOV_GYMNOSPERM","XRCOV_MONOCOT"
                         ,"XRCOV_LYCOPOD","XRCOV_HORSETAIL"
                         ,"XRCOV_QUILLWORT","XRCOV_LIVERWORT")
  
  empty_base.nat <- data.frame(t(rep(NA,40)), stringsAsFactors=F)
  names(empty_base.nat) <- c("N_DICOTS_NAT","PCTN_DICOTS_NAT","XABCOV_DICOTS_NAT","XRCOV_DICOTS_NAT"     
                             ,"N_DICOTS_ALIEN","PCTN_DICOTS_ALIEN","XABCOV_DICOTS_ALIEN"
                             ,"XRCOV_DICOTS_ALIEN","N_DICOTS_CRYP","PCTN_DICOTS_CRYP"
                             ,"XABCOV_DICOTS_CRYP","XRCOV_DICOTS_CRYP","N_DICOTS_AC","PCTN_DICOTS_AC"       
                             ,"XABCOV_DICOTS_AC","XRCOV_DICOTS_AC","N_FERNS_NAT","PCTN_FERNS_NAT"
                             ,"XABCOV_FERNS_NAT","XRCOV_FERNS_NAT","N_FERNS_INTR","PCTN_FERNS_INTR"
                             ,"XABCOV_FERNS_INTR","XRCOV_FERNS_INTR","N_MONOCOTS_NAT","PCTN_MONOCOTS_NAT"
                             ,"XABCOV_MONOCOTS_NAT","XRCOV_MONOCOTS_NAT","N_MONOCOTS_ALIEN"
                             ,"PCTN_MONOCOTS_ALIEN","XABCOV_MONOCOTS_ALIEN","XRCOV_MONOCOTS_ALIEN"
                             ,"N_MONOCOTS_CRYP","PCTN_MONOCOTS_CRYP","XABCOV_MONOCOTS_CRYP"
                             ,"XRCOV_MONOCOTS_CRYP","N_MONOCOTS_AC","PCTN_MONOCOTS_AC"
                             ,"XABCOV_MONOCOTS_AC","XRCOV_MONOCOTS_AC")
  
  
  if('NWCA_NATSTAT' %in% names(vascIn)){
    vascIn.2 <- vascIn.1
    vascIn.2$ALIEN <- with(vascIn.2, ifelse(NWCA_NATSTAT %in% c('INTR','ADV'), 1, 0))
    vascIn.2$NATSTAT_ALT <- with(vascIn.2, ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),'ALIEN', NWCA_NATSTAT))
    vascIn.2$AC <- with(vascIn.2, ifelse(NWCA_NATSTAT %in% c('INTR','ADV','CRYP'), 1, 0))
    
    vascIn.2$DICOTS_NAT <- with(vascIn.2, ifelse(CATEGORY=='DICOT' & NATSTAT_ALT=='NAT', 1, 0))
    vascIn.2$DICOTS_ALIEN <- with(vascIn.2, ifelse(CATEGORY=='DICOT' & NATSTAT_ALT=='ALIEN', 1, 0))
    vascIn.2$DICOTS_CRYP <- with(vascIn.2, ifelse(CATEGORY=='DICOT' & NATSTAT_ALT=='CRYP', 1, 0))
    vascIn.2$DICOTS_AC <- with(vascIn.2, ifelse(CATEGORY=='DICOT' & AC==1, 1, 0))
    vascIn.2$FERNS_NAT <- with(vascIn.2, ifelse(CATEGORY=='FERN' & NATSTAT_ALT=='NAT', 1, 0))
    vascIn.2$FERNS_INTR <- with(vascIn.2, ifelse(CATEGORY=='FERN' & NWCA_NATSTAT=='INTR', 1, 0))
    vascIn.2$MONOCOTS_NAT <- with(vascIn.2, ifelse(CATEGORY=='MONOCOT' & NATSTAT_ALT=='NAT', 1, 0))
    vascIn.2$MONOCOTS_ALIEN <- with(vascIn.2, ifelse(CATEGORY=='MONOCOT' & NATSTAT_ALT=='ALIEN', 1, 0))
    vascIn.2$MONOCOTS_CRYP <- with(vascIn.2, ifelse(CATEGORY=='MONOCOT' & NATSTAT_ALT=='CRYP', 1, 0))
    vascIn.2$MONOCOTS_AC <- with(vascIn.2, ifelse(CATEGORY=='MONOCOT' & AC==1, 1, 0))
    
    multTraits <- int.combTraits(vascIn.2,c('DICOTS_NAT','DICOTS_ALIEN','DICOTS_CRYP','DICOTS_AC','FERNS_NAT','FERNS_INTR','MONOCOTS_NAT','MONOCOTS_ALIEN'
                                            ,'MONOCOTS_CRYP','MONOCOTS_AC'), sampID)
    catOut <- rbind(catOut, multTraits)
    
    empty_base <- cbind(empty_base, empty_base.nat)
    
  }
  
  catOut.1a <- reshape(catOut, idvar = c(sampID), direction = 'wide',
                       timevar = 'PARAMETER', v.names = 'RESULT')
  
  names(catOut.1a) <- gsub("RESULT\\.", "", names(catOut.1a))
  catOut.1a <- merge(catOut.1a, empty_base, all=TRUE)
  
  catOut.1b <- reshape(catOut.1a, idvar = sampID, direction = 'long',
                       varying= names(catOut.1a)[!names(catOut.1a) %in% c(sampID)],
                       timevar = 'PARAMETER', v.names = 'RESULT',
                       times = names(catOut.1a)[!names(catOut.1a) %in% c(sampID)])
  catOut.1b <- subset(catOut.1b, !is.na(eval(as.name(sampID[1]))))
  catOut.1b$RESULT <- with(catOut.1b, ifelse(is.na(RESULT), 0, RESULT))
  catOut.1b$PARAMETER <- with(catOut.1b, as.character(PARAMETER))
  
  catOut.1 <- subset(catOut.1b, PARAMETER %in% names(empty_base))
  
  return(catOut.1)
}
