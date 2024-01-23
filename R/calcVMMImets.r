# CalculateMetricSubsets.r


#' @export
#' 
#' @title Calculate only metrics used in NWCA 2011 VMMI
#' 
#' @description This function calculates FQAI_ALL, N_TOL, RIMP_NATSPP,
#' and XRCOV_MONOCOTS_NAT metrics, which are used in the NWCA 2011
#' Vegetation Multimetric Index (VMMI).
#' 
#' @param vascIn Data frame containing cover data summarized at the UID
#' and TAXON level:
#' \itemize{
#'  \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#' \item TAXON: Taxon name
#'
#' \item NWcA_CC: Coefficient of Conservatism as used in NWCA, valid
#' values range between 0 and 10
#'
#' \item NWCA_NATSTAT: Native status as used in NWCA, with native taxa
#' indicated by 'NAT'
#'
#' \item CATEGORY: Category abbreviations from USDA PLANTS database,
#' with monocots indicated by 'MONOCOT'
#'
#' \item sXRCOV: Proportion of summed cover across all taxa
#' (XTOTABCOV) represented by taxon in sample
#'
#' \item sRFREQ: Relative frequency of a taxon, calculated as the
#' percentage of the total frequency of taxon occurrence across
#' all taxa for a UID.
#' }
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#' 
#' @details To calculate the VMMI as used in NWCA 2011, the default
#' taxa lists must be used to create the input data frame.
#' 
#' @return Data frame containing \emph{sampID} variables, FQAI_ALL, N_TOL,
#'   RIMP_NATSPP, and XRCOV_MONOCOTS_NAT, with the following definitions:
#' \itemize{
#' \item FQAI_ALL: Floristic Quality Assessment Index, based on all
#' species
#'
#' \item N_TOL: Number of tolerant species, as defined by C-value
#'
#' \item RIMP_NATSPP: Relative importance of native species
#'
#' \item XRCOV_MONOCOTS_NAT: Mean relative cover of native monocots
#' }
#' 
#' @references US Environmental Protection Agency. 2016.
#' National Wetland Condition Assessment: 2011 Technical Report.
#' EPA-843-R-15-006. US Environmental Protection Agency,
#' Washington, DC.
#' 
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(VascPlantEx)
#'  exPlant <- prepareData(VascPlantEx, taxon_name = 'USDA_NAME', 
#'  inTaxa = taxaNWCA, inNat = ccNatNWCA, inCVal = ccNatNWCA, 
#'  inWIS = wisNWCA, cValReg='STATE')
#'
#' vmmiEx <- calcVMMImets(exPlant$byUIDspp)
#'
#' head(vmmiEx)
#' unique(vmmiEx$PARAMETER)

calcVMMImets <- function(vascIn, sampID='UID'){
  vascIn <- as.data.frame(vascIn) # Do this in case read in as a tibble or data.table, which might cause problems
  # Calculate only the 4 metrics used in the VMMI: FQAI_ALL,N_TOL,RIMP_NATSPP,XRCOV_MONOCOTS_NAT
  # First obtain all of the UIDs in input data frame
  for(i in 1:length(sampID)){
    if(i==1) vascIn$SAMPID <- vascIn[,sampID[i]]
    else vascIn$SAMPID <- paste(vascIn$SAMPID, vascIn[,sampID[i]], sep='.')
  }
  samples <- unique(subset(vascIn, select=c(sampID, 'SAMPID')))
    
  UIDs <- data.frame(SAMPID=unique(subset(vascIn, select='SAMPID')), stringsAsFactors=FALSE)

  # Look at vascIn for necessary variables
  if(sum(c('NWCA_CC','NWCA_NATSTAT','CATEGORY') %in% names(vascIn))<3){
    print("Missing key variable(s) from input data frame! NWCA_CC, NWCA_NATSTAT, and CATEGORY
          required to calculate VMMI metrics")
    return(NULL)
  }
  if(sum(c('sXRCOV','sRFREQ') %in% names(vascIn))<2){
    print("Missing key sums in input data. Must provide sRFREQ and sXRCOV to calculate metrics!
          See help for more info.")
    return(NULL)
  }

  ## Calculate FQAI_ALL
  vascIn.fqai <- subset(vascIn, toupper(NWCA_CC)!='UND' & NWCA_CC!='')
  
  numtaxa <- aggregate(x = list(NUMTAXA = vascIn.fqai$TAXON), by = vascIn.fqai[c('SAMPID')], 
                       FUN = function(x){length(unique(x))})
  
  vascIn.fqai <- merge(vascIn.fqai, numtaxa, by = 'SAMPID')
  vascIn.fqai$NWCA_CC <- as.numeric(vascIn.fqai$NWCA_CC)
  
  fqaiOut <- aggregate(x = list(FQAI_ALL = vascIn.fqai$NWCA_CC), 
            by = vascIn.fqai[c('SAMPID','NUMTAXA')], FUN = sum)
  
  fqaiOut$FQAI_ALL <- with(fqaiOut, round(FQAI_ALL/sqrt(NUMTAXA), 2))

  # Calculate N_TOL
  vascIn.alt <- vascIn
  vascIn.alt$TOL <- with(vascIn.alt, ifelse(NWCA_CC %in% c('4','3','2','1','0'), 1, 0))
  
  vascIn.alt <- subset(vascIn.alt, TOL==1)

  if(nrow(vascIn.alt)>0){
    ntol <- aggregate(x = list(N_TOL = vascIn.alt$TAXON), 
                      by = vascIn.alt[c('SAMPID')], FUN = length)
    
    ntolOut <- merge(UIDs, ntol, by='SAMPID', all.x=TRUE)
    ntolOut$N_TOL <- with(ntolOut, ifelse(is.na(N_TOL), 0, N_TOL))

  }else{
    ntolOut <- data.frame(SAMPID = UIDs$SAMPID, N_TOL = 0, stringsAsFactors = F)
  }

  # Calculate RIMP_NATSPP
  vascIn.nat <- subset(vascIn,NWCA_NATSTAT=='NAT')

  if(nrow(vascIn.nat)>0){
    vascIn.nat.1 <- aggregate(x = list(XRCOV = vascIn.nat$sXRCOV, RFREQ = vascIn.nat$sRFREQ),
                              by = vascIn.nat[c('SAMPID')], FUN = sum)
    
    vascIn.nat.1$XRCOV <- with(vascIn.nat.1, round(XRCOV, 2))
    vascIn.nat.1$RFREQ <- with(vascIn.nat.1, round(RFREQ, 2))
    vascIn.nat.1$RIMP_NATSPP <- with(vascIn.nat.1, round((RFREQ + XRCOV)/2, 2))
    
    natOut <- merge(UIDs, vascIn.nat.1, by = 'SAMPID', all.x=TRUE)
    natOut$RIMP_NATSPP <- with(natOut, ifelse(is.na(RIMP_NATSPP), 0, RIMP_NATSPP))
    
   }else{
    natOut <- data.frame(SAMPID = UIDs$SAMPID, RIMP_NATSPP = 0, stringsAsFactors = F)
  }

  # Calculate XRCOV_MONOCOTS_NAT
  vascIn.mono <- subset(vascIn,CATEGORY=='MONOCOT' & NWCA_NATSTAT=='NAT')

  if(nrow(vascIn.mono)>0){
    vascIn.mono.1 <- aggregate(x = list(XRCOV_MONOCOTS_NAT = vascIn.mono$sXRCOV), 
                               by = vascIn.mono[c('SAMPID')], FUN = sum)
    vascIn.mono.1$XRCOV_MONOCOTS_NAT <- with(vascIn.mono.1, round(XRCOV_MONOCOTS_NAT, 2))

    monoOut <- merge(UIDs, vascIn.mono.1, by = 'SAMPID', all.x=TRUE)
    monoOut$XRCOV_MONOCOTS_NAT <- with(monoOut, ifelse(is.na(XRCOV_MONOCOTS_NAT), 0, XRCOV_MONOCOTS_NAT))

  }else{
    monoOut <- data.frame(SAMPID=UIDs$SAMPID, XRCOV_MONOCOTS_NAT = 0, stringsAsFactors=F)
  }

 # Now combine into one data frame
 allOut <- merge(fqaiOut, natOut, by='SAMPID') 
 allOut <- merge(allOut, monoOut, by='SAMPID')
 allOut <- merge(allOut, ntolOut, by='SAMPID')

 allOut.1 <- merge(allOut, samples, by = 'SAMPID')
 allOut.1$SAMPID <- NULL
 
 allOut.1 <- allOut.1[,c(sampID,'FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')]

 return(allOut.1)

}



