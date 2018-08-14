# CalculateMetricSubsets.r

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
#'     \item XABCOV: Mean percent cover of taxon across plots
#'
#'     \item TOTN: Number of taxa in sample
#'
#'     \item sXRCOV: proportion of summed cover across all taxa
#'     (XTOTABCOV) represented by taxon in sample
#'
#'     \item NWCA_NATSTAT (Optional): Native status variable with categories
#'     of 'NAT','ADV','CRYP','INTR','UND'
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
#'  exPlant <- prepareData(VascPlantEx)
#'
#'  durEx <- calcDuration(exPlant$byUIDspp)
#'  head(durEx)
#'  unique(durEx$PARAMETER)

calcDuration <- function(vascIn,sampID='UID'){
 # First check for necessary variables
  necNames <- c(sampID,'TAXON','DURATION','XABCOV','TOTN','sXRCOV')
  msgNames <- necNames[necNames %nin% names(vascIn)]
  if(length(msgNames)>0){
    print(paste("Missing key variables for metric calculation: ",paste(msgNames,collapse=','),
                ". Try prepareData() function to create necessary input variables.",sep=''))
  return(NULL)
  }

  # From DURATION, create DUR_ALT variable
  vascIn.1 <- plyr::mutate(vascIn,DUR_ALT=car::recode(DURATION,"c('ANNUAL, BIENNIAL','BIENNIAL')='ANN_BIEN';
                  c('ANNUAL, BIENNIAL, PERENNIAL','ANNUAL, PERENNIAL','PERENNIAL, ANNUAL', 'BIENNIAL, PERENNIAL')='ANN_PEREN'"))

  durOut <- int.calcTraits_MultCat(vascIn.1,'DUR_ALT',sampID)
  
  empty_base <- data.frame(t(rep(NA,16)),stringsAsFactors=F)
  names(empty_base) <- c("N_ANN_BIEN","N_ANN_PEREN","N_ANNUAL","N_PERENNIAL", "PCTN_ANN_BIEN","PCTN_ANNUAL"
                         ,"PCTN_ANN_PEREN","PCTN_PERENNIAL","XABCOV_ANN_BIEN","XABCOV_ANN_PEREN"
                         ,"XABCOV_ANNUAL","XABCOV_PERENNIAL","XRCOV_ANN_BIEN","XRCOV_ANN_PEREN"
                         ,"XRCOV_ANNUAL","XRCOV_PERENNIAL")

  if('NWCA_NATSTAT' %in% names(vascIn.1)){
    # Assign native status values to grouped categories
    vascIn.2 <- plyr::mutate(vascIn.1,ALIEN=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),1,0)
                              ,NATSTAT_ALT=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),'ALIEN',NWCA_NATSTAT)
                              ,AC=ifelse(NWCA_NATSTAT %in% c('INTR','ADV','CRYP'),1,0))


    vascIn.2 <- plyr::mutate(vascIn.2,ANNUAL_NAT=ifelse(DUR_ALT=='ANNUAL' & NATSTAT_ALT=='NAT',1,0)
                                 ,ANNUAL_AC=ifelse(DUR_ALT=='ANNUAL' & AC==1,1,0)
                                 ,ANN_BIEN_NAT=ifelse(DUR_ALT=='ANN_BIEN' & NATSTAT_ALT=='NAT',1,0)
                                 ,ANN_BIEN_AC=ifelse(DUR_ALT=='ANN_BIEN' & AC==1,1,0)
                                 ,ANN_PEREN_NAT=ifelse(DUR_ALT=='ANN_PEREN' & NATSTAT_ALT=='NAT',1,0)
                                 ,ANN_PEREN_AC=ifelse(DUR_ALT=='ANN_PEREN' & AC==1,1,0)
                                 ,PERENNIAL_NAT=ifelse(DUR_ALT=='PERENNIAL' & NATSTAT_ALT=='NAT',1,0)
                                 ,PERENNIAL_AC=ifelse(DUR_ALT=='PERENNIAL' & AC==1,1,0)
                           
    
      )

    multTraits <- int.combTraits(vascIn.2,c('ANNUAL_NAT','ANNUAL_AC','ANN_BIEN_NAT','ANN_BIEN_AC','ANN_PEREN_NAT'
                                            ,'ANN_PEREN_AC','PERENNIAL_NAT','PERENNIAL_AC'),sampID)
    
    durOut <- rbind(durOut,multTraits)
    
    empty_base.nat <- data.frame(t(rep(NA,32)),stringsAsFactors=F)
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


  }

  return(durOut)

}

# Growth habit metrics (with and without native status, if available)
#' @export
#' @title Calculate growth habit metrics
#' @description This function calculates all growth habit metrics, with
#' additional metrics if NWCA_NATSTAT variable is included
#' in input data frame.
#' @details Both GROWTH_HABIT and NWCA_NATSTAT variables are recoded to fewer
#' categories. Taxa with 'UND' as native status are excluded.
#' @param vascIn   Data frame containing cover data summarized by
#' UID and TAXON, with the following fields:
#' \itemize{
#'     \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#'     \item TAXON: Taxon name
#'
#'     \item GROWTH_HABIT: USDA PLANTS GROWTH_HABIT variable with valid
#'     values FORB, GRAMINOID, HERB, SHRUB, SUBSHRUB, TREE, VINE,
#'     NONVASCULAR, combinations of these, or blank
#'
#'     \item XABCOV: Mean percent cover of taxon across plots
#'
#'     \item TOTN: Number of taxa in sample
#'
#'     \item sXRCOV: proportion of summed cover across all taxa
#'     (XTOTABCOV) represented by taxon in sample
#'
#'     \item NWCA_NATSTAT (optional): Native status variable with categories
#'     of 'NAT','ADV','CRYP','INTR','UND'
#'    }
#' @param sampID  A character vector containing the name(s) of variable(s)
#'   necessary to identify unique samples, 'UID' by default
#' @return Data frame containing \emph{sampID} variables, PARAMETER, RESULT,
#'   where values of PARAMETER consist of the metric name concatenated with
#'   trait value (represented as TRAITNM below): 
#'\itemize{
#' \item N_TRAITNM: Number of taxa with trait
#'
#' \item PCTN_TRAITNM: Number of taxa with trait as percentage of \emph{TOTN}
#'
#' \item XABCOV_TRAITNM: Sum of \emph{XABCOV} values across taxa with trait
#'
#' \item XRCOV_TRAITNM: Sum of \emph{sXRCOV} values across taxa with trait
#' }
#' If NWCA_NATSTAT is in the input data frame, the same set except
#' 'N_' metrics is calculated for native species and
#' alien + cryptogenic species, with metric names suffixes of
#' _NAT and _AC, respectively.
#' @author Karen Blocksom \email{Blocksom.karen@epa.gov}
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC.
#' @examples
#' head(VascPlantEx)
#' exPlant <- prepareData(VascPlantEx)
#'
#' ghEx <- calcGrowthHabit(exPlant$byUIDspp)
#'
#' head(ghEx)
#' unique(ghEx$PARAMETER)

calcGrowthHabit <- function(vascIn,sampID='UID'){
  # First check for necessary variables
  necNames <- c(sampID,'TAXON','GROWTH_HABIT','XABCOV','TOTN','sXRCOV')
  msgNames <- necNames[necNames %nin% names(vascIn)]
  if(length(msgNames)>0){
    print(paste("Missing key variables for metric calculation: ",paste(msgNames,collapse=','),
                ". Try prepareData() function to create necessary input variables.",sep=''))
    return(NULL)
  }

  vascIn.1 <- plyr::mutate(vascIn,GRH_ALT=car::recode(GROWTH_HABIT,"c('FORB/HERB','FORB/HERB, SHRUB','FORB/HERB, SHRUB, SUBSHRUB'
                          ,'FORB/HERB, SUBSHRUB')='FORB';c('SUBSHRUB, FORB/HERB','SUBSHRUB, SHRUB, FORB/HERB')='SSHRUB_FORB'
                          ;c('SUBSHRUB','SUBSHRUB, SHRUB','SHRUB, SUBSHRUB','SUBSHRUB, SHRUB, TREE')='SSHRUB_SHRUB'
                          ;c('SHRUB','SHRUB, TREE','TREE, SUBSHRUB, SHRUB','SHRUB, SUBSHRUB, TREE')='SHRUB'
                          ;c('TREE, SHRUB','TREE, SHRUB, VINE')='TREE_SHRUB'
                          ;c('VINE, FORB/HERB','SUBSHRUB, FORB/HERB, VINE','FORB/HERB, VINE')='VINE'
                          ;c('VINE, SHRUB','VINE, SUBSHRUB','SUBSHRUB, VINE','SHRUB, VINE','SHRUB, FORB/HERB, SUBSHRUB, VINE'
                          ,'SHRUB, SUBSHRUB, VINE')='VINE_SHRUB'
                          ;c('GRAMINOID','SUBSHRUB, SHRUB, GRAMINOID')='GRAMINOID'")
                           ,TREE_COMB=ifelse(GRH_ALT %in% c('TREE','TREE_SHRUB'),1,0)
                           ,SHRUB_COMB=ifelse(GRH_ALT %in% c('SSHRUB_SHRUB','SSHURB_FORB','SHRUB'),1,0)
                           ,HERB=ifelse(GRH_ALT %in% c('GRAMINOID','FORB'),1,0))

  sppGRH <- int.calcTraits_MultCat(vascIn.1,'GRH_ALT',sampID)
  sppGRH.1 <- int.combTraits(vascIn.1,c('TREE_COMB','SHRUB_COMB','HERB'),sampID)
  grhOut <- rbind(sppGRH,sppGRH.1)

  if('NWCA_NATSTAT' %in% names(vascIn.1)){
    vascIn.2 <- plyr::mutate(vascIn.1,ALIEN=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),1,0)
                              ,NATSTAT_ALT=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),'ALIEN',NWCA_NATSTAT)
                              ,AC=ifelse(NWCA_NATSTAT %in% c('INTR','ADV','CRYP'),1,0))

    vascIn.2 <- plyr::mutate(vascIn.2,GRAMINOID_AC=ifelse(GRH_ALT=='GRAMINOID' & AC==1,1,0)
                         ,GRAMINOID_NAT=ifelse(GRH_ALT=='GRAMINOID' & NATSTAT_ALT=='NAT',1,0)
                         ,FORB_AC=ifelse(GRH_ALT=='FORB' & AC==1,1,0),FORB_NAT=ifelse(GRH_ALT=='FORB' & NATSTAT_ALT=='NAT',1,0)
                         ,HERB_AC=ifelse(HERB==1 & AC==1,1,0),HERB_NAT=ifelse(HERB==1 & NATSTAT_ALT=='NAT',1,0)
                         ,SHRUB_COMB_AC=ifelse(SHRUB_COMB==1 & AC==1,1,0),SHRUB_COMB_NAT=ifelse(SHRUB_COMB==1 & NATSTAT_ALT=='NAT',1,0)
                         ,TREE_COMB_AC=ifelse(TREE_COMB==1 & AC==1,1,0),TREE_COMB_NAT=ifelse(TREE_COMB==1 & NATSTAT_ALT=='NAT',1,0)
                         ,VINE_AC=ifelse(GRH_ALT=='VINE' & AC==1,1,0),VINE_NAT=ifelse(GRH_ALT=='VINE' & NATSTAT_ALT=='NAT',1,0)
                         ,VINE_SHRUB_AC=ifelse(GRH_ALT=='VINE_SHRUB' & AC==1,1,0)
                         ,VINE_SHRUB_NAT=ifelse(GRH_ALT=='VINE_SHRUB' & NATSTAT_ALT=='NAT',1,0)
      )

    multTraits <- int.combTraits(vascIn.2,c('GRAMINOID_AC','GRAMINOID_NAT','FORB_AC','FORB_NAT','HERB_AC','HERB_NAT','SHRUB_COMB_AC','SHRUB_COMB_NAT'
                                        ,'TREE_COMB_AC','TREE_COMB_NAT','VINE_AC','VINE_NAT','VINE_SHRUB_AC','VINE_SHRUB_NAT')
                              ,sampID)
    grhOut <- rbind(grhOut,multTraits)

  }

  return(grhOut)
}

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
#'     of DICOT, FERN, GYMNOSPERM, HORSETAIL, LICHEN, LIVERWORT,
#'     LYCOPOD, MONOCOT, MOSS, or blank
#'
#'     \item XABCOV: Mean percent cover of taxon across plots
#'
#'     \item TOTN: Number of taxa in sample
#'
#'     \item sXRCOV: proportion of summed cover across all taxa
#'     (XTOTABCOV) represented by taxon in sample
#'
#'     \item NWCA_NATSTAT (optional): Native status variable with categories
#'     of 'NAT','ADV','CRYP','INTR','UND'
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
#' exPlant <- prepareData(VascPlantEx)
#'
#' catEx <- calcCategory(exPlant$byUIDspp)
#'
#' head(catEx)
#' unique(catEx$PARAMETER)

calcCategory <- function(vascIn,sampID='UID'){
  # First check for necessary variables
  necNames <- c(sampID,'TAXON','CATEGORY','XABCOV','TOTN','sXRCOV')
  msgNames <- necNames[necNames %nin% names(vascIn)]
  if(length(msgNames)>0){
    print(paste("Missing key variables for metric calculation: ",paste(msgNames,collapse=','),
                ". Try prepareData() function to create necessary input variables.",sep=''))
    return(NULL)
  }

  vascIn.1 <- subset(vascIn,!is.na(CATEGORY))
  catOut <- int.calcTraits_MultCat(vascIn.1,'CATEGORY',sampID)

  if('NWCA_NATSTAT' %in% names(vascIn)){
    vascIn.2 <- plyr::mutate(vascIn.1,ALIEN=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),1,0)
                           ,NATSTAT_ALT=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),'ALIEN',NWCA_NATSTAT)
                           ,AC=ifelse(NWCA_NATSTAT %in% c('INTR','ADV','CRYP'),1,0))
    vascIn.2 <- plyr::mutate(vascIn.2,DICOTS_NAT=ifelse(CATEGORY=='DICOT' & NATSTAT_ALT=='NAT',1,0)
                               ,DICOTS_ALIEN=ifelse(CATEGORY=='DICOT' & NATSTAT_ALT=='ALIEN',1,0)
                               ,DICOTS_CRYP=ifelse(CATEGORY=='DICOT' & NATSTAT_ALT=='CRYP',1,0)
                               ,DICOTS_AC=ifelse(CATEGORY=='DICOT' & AC==1,1,0)
                               ,FERNS_NAT=ifelse(CATEGORY=='FERN' & NATSTAT_ALT=='NAT',1,0)
                               ,FERNS_INTR=ifelse(CATEGORY=='FERN' & NWCA_NATSTAT=='INTR',1,0)
                               ,MONOCOTS_NAT=ifelse(CATEGORY=='MONOCOT' & NATSTAT_ALT=='NAT',1,0)
                               ,MONOCOTS_ALIEN=ifelse(CATEGORY=='MONOCOT' & NATSTAT_ALT=='ALIEN',1,0)
                               ,MONOCOTS_CRYP=ifelse(CATEGORY=='MONOCOT' & NATSTAT_ALT=='CRYP',1,0)
                               ,MONOCOTS_AC=ifelse(CATEGORY=='MONOCOT' & AC==1,1,0))


    multTraits <- int.combTraits(vascIn.2,c('DICOTS_NAT','DICOTS_ALIEN','DICOTS_CRYP','DICOTS_AC','FERNS_NAT','FERNS_INTR','MONOCOTS_NAT','MONOCOTS_ALIEN'
                                    ,'MONOCOTS_CRYP','MONOCOTS_AC'),sampID)
    catOut <- rbind(catOut,multTraits)

    }

  return(catOut)
}

# Wetland indicator status metrics (with and without native status, if available)
#' @export
#' 
#' @title Calculate Wetland Indicator Status metrics
#' 
#' @description This function calculates Wetland Indicator Status (WIS)
#' metrics, including variations based on native status,
#' if NWCA_NATSTAT is present in the input data frame.
#' 
#' @param vascIn Data frame containing cover data summarized by
#' UID and TAXON, with the following fields:
#' \itemize{
#'     \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#'     \item TAXON: Taxon name
#'
#'     \item WIS: Wetland Indicator Status, with possible values
#'     of FAC, FACU, FACW, OBL, UPL, NL, TBD, or missing.
#'     NL is recoded to UPL in function, TBD is recoded
#'     to missing in function.
#'
#'     \item FREQ: Frequency of taxon among plots at site
#'
#'     \item XABCOV: Mean percent cover of taxon across plots
#'
#'     \item TOTN: Number of taxa in sample
#'
#'     \item sXRCOV: proportion of summed cover across all taxa
#'     (XTOTABCOV) represented by taxon in sample
#'
#'     \item NWCA_NATSTAT (optional): Native status variable with
#'       categories of 'NAT','ADV','CRYP','INTR','UND'.
#'       UND taxa are ignored.
#'    }
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default

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
#' In addition WIS indices based on all and native species only (if
#' NWCA_NATSTAT is provided in the input data frame), with the
#' suffixes ALL and NAT, respectively. WIS values recoded as numeric
#' with OBL=1, FACW=2, FAC=3, FACU=4, UPL=5:
#' \itemize{ \item WETIND_COV_ALL, WETIND_COV_NAT: Wetland Index, weighted by
#' taxon cover
#' 
#' \item WETIND_FREQ_ALL, WETIND_FREQ_NAT: Wetland Index, weight by frequency } 
#' 
#' If NWCA_NATSTAT is provided in the input data frame, the following metrics
#' are calculated based on Alien + Cryptogenic species: 
#'\itemize{ \item
#' N_OBLFACW_AC: Number of alien and cryptogenic taxa with WIS of either OBL or
#' FACW.
#' 
#' \item XABCOV_OBLFACW_AC: Mean absolute cover of alien and cryptogenic taxa 
#' with WIS of either OBL or FACW.
#' 
#' \item XRCOV_OBLFACW_AC: Mean relative cover of alien and cryptogenic taxa 
#' with WIS of either OBL or FACW. }

#' @author Karen Blocksom \email{Blocksom.karen@epa.gov}
#' 
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC.
#' 
#' @examples
#' head(VascPlantEx)
#' exPlant <- prepareData(VascPlantEx)
#'
#' wisEx <- calcWIS(exPlant$byUIDspp)
#'
#' head(wisEx)
#' unique(wisEx$PARAMETER)

calcWIS <- function(vascIn,sampID='UID'){
  # First check for necessary variables
  necNames <- c(sampID,'TAXON','WIS','FREQ','XABCOV','TOTN','sXRCOV')
  msgNames <- necNames[necNames %nin% names(vascIn)]
  if(length(msgNames)>0){
    print(paste("Missing key variables for metric calculation: ",paste(msgNames,collapse=','),
                ". Try prepareData() function to create necessary input variables.",sep=''))
    return(NULL)
  }

  vascIn <- plyr::mutate(vascIn,ECOIND=car::recode(WIS,"c('FACU')=4;c('FAC')=3;c('FACW')=2;c('OBL')=1;
                                  c('UPL','NL')=5;c('TBD',NA)=NA")
                        ,WIS=car::recode(WIS,"c('UPL','NL')='UPL';c(NA,'TBD')=NA"))

  # Overall metric calculations
  vascIn.1 <- subset(vascIn,!is.na(WIS))
  sppWIS <- int.calcTraits_MultCat(vascIn.1,'WIS',sampID)

  ## Calculate Wetland indicator status metrics and melt into long format
  vascIn.2 <- subset(vascIn.1,!is.na(ECOIND)) %>%
    plyr::mutate(ECOIND=as.numeric(ECOIND)) %>%
    plyr::ddply(c(sampID),summarise
                        ,WETIND_COV_ALL=round(sum(XABCOV*ECOIND)/sum(XABCOV),2)
                        ,WETIND_FREQ_ALL=round(sum(FREQ*ECOIND)/sum(FREQ),2)) %>%
    reshape2::melt(id.vars=c(sampID),variable.name='PARAMETER',value.name='RESULT')

  wisOut <- rbind(sppWIS,vascIn.2)

  # Metrics using only subsets of data based on NATSTAT_ALT
  if('NWCA_NATSTAT' %in% names(vascIn)){
    vascIn.nat <- plyr::mutate(vascIn,ALIEN=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),1,0)
                           ,NATSTAT_ALT=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),'ALIEN',NWCA_NATSTAT)
                           ,AC=ifelse(NWCA_NATSTAT %in% c('INTR','ADV','CRYP'),1,0))

    vascIn.nat.1 <- subset(vascIn.nat,NWCA_NATSTAT=='NAT')

    ## Calculate Wetland indicator status metrics and melt into long format
    wisOut.nat <- subset(vascIn.nat.1,!is.na(ECOIND)) %>%
      plyr::mutate(ECOIND=as.numeric(ECOIND)) %>%
      plyr::ddply(c(sampID),summarise,WETIND_COV_NAT=round(sum(XABCOV*ECOIND)/sum(XABCOV),2)
                          ,WETIND_FREQ_NAT=round(sum(FREQ*ECOIND)/sum(FREQ),2)) %>%
            reshape2::melt(id.vars=c(sampID),variable.name='PARAMETER',value.name='RESULT') %>%
              plyr::mutate(PARAMETER=as.character(PARAMETER))

    # Obligate and facultative wet alien and cryptogenic species
    vascIn.obl <- plyr::mutate(vascIn.nat,OBLFACW_AC=ifelse(WIS %in% c('OBL','FACW') & NATSTAT_ALT %in% c('ALIEN','CRYP'),1,0))
    ofOut <- int.calcTraits_Indicator(vascIn.obl,'OBLFACW_AC',sampID) %>%
      plyr::mutate(PARAMETER=as.character(PARAMETER)) %>%
      subset(PARAMETER %nin% c('PCTN_OBLFACW_AC'))

    wisOut <- rbind(wisOut,wisOut.nat,ofOut)

  }

  outdf <- plyr::mutate(wisOut,RESULT=ifelse(is.na(RESULT)|is.infinite(RESULT),0,RESULT))

  return(outdf)
}


# Metrics using CC values
#' @export
#' 
#' @title Calculate Wetland Indicator Status metrics
#' 
#' @description This function calculates Wetland Indicator Status (WIS)
#' metrics, including variations based on native status,
#' if NWCA_NATSTAT is present in the input data frame.
#' 
#' @param vascIn Data frame containing cover data summarized by
#' UID and TAXON, with the following fields:
#' \itemize{
#'     \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#'     \item TAXON: Taxon name
#'
#'     \item XABCOV: Mean percent cover of taxon across plots
#'
#'     \item FREQ: Number plots in which taxon occurs
#'
#'     \item NWCA_CC: Coefficient of conservatism values by taxon
#'     from NWCA
#'
#'     \item NWCA_NATSTAT (optional): Native status variable with
#'       categories of 'NAT','ADV','CRYP','INTR','UND'.
#'       UND taxa are ignored.
#'    }
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default

#' @return   Data frame containing \emph{sampID} variables, PARAMETER, RESULT,
#'   where values of PARAMETER are:
#' \itemize{
#'     \item XC: Mean coefficient of conservatism (unweighted)
#'
#'     \item FQAI: Floral Quality Assessment Index (unweighted)
#'
#'     \item XC_FREQ: Mean coefficient of conservatism weighted by
#'     relative frequency
#'
#'     \item FQAI_FREQ: Floral Quality Assessment Index weighted by
#'     relative frequency
#'
#'     \item XC_COV: Mean coefficient of conservatism weighted by
#'     relative absolute cover
#'
#'     \item FQAI_COV: Floral Quality Assessment Index weighted by
#'     relative absolute cover
#'
#'     \item For metrics based native species, the metric name has a
#'     suffix of NAT.
#'    }
#' 
#' @author Karen Blocksom \email{Blocksom.karen@epa.gov}
#' 
#' @references US Environmental Protection Agency. 2016. National Wetland 
#'   Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US 
#'   Environmental Protection Agency, Washington, DC.
#'   
#' @examples
#' head(VascPlantEx)
#' exPlant <- prepareData(VascPlantEx)
#'
#' ccEx <- calcCC(exPlant$byUIDspp)
#'
#' head(ccEx)
#' unique(ccEx$PARAMETER)
calcCC <- function(vascIn,sampID='UID'){
  if('NWCA_CC' %nin% names(vascIn)){
    print("Missing NWCA_CC from input data frame - cannot calculate metrics!")
    return(NULL)
  }

  ## Calculate mean CC and FQAI indices
  totals <- plyr::ddply(vascIn,c(sampID),summarise,SUBTOTFREQ=sum(FREQ),SUBXTOTABCOV=sum(XABCOV))
  vascIn.1 <- merge(subset(vascIn,select=names(vascIn) %nin% c('TOTFREQ','XTOTABCOV')),totals,by=sampID)
  vascIn.2 <- plyr::ddply(subset(vascIn.1,NWCA_CC!='Und'),c(sampID),summarise,XC=round(sum(as.numeric(NWCA_CC))/length(unique(TAXON)),2)
                        ,FQAI=round(sum(as.numeric(NWCA_CC))/sqrt(length(unique(TAXON))),2)
                        ,XC_FREQ=round(sum((FREQ/SUBTOTFREQ)*100*as.numeric(NWCA_CC))/length(unique(TAXON)),2)
                        ,FQAI_FREQ=round(sum((FREQ/SUBTOTFREQ)*100*as.numeric(NWCA_CC))/sqrt(length(unique(TAXON))),2)
                        ,XC_COV=round(sum((XABCOV/SUBXTOTABCOV)*100*as.numeric(NWCA_CC))/length(unique(TAXON)),2)
                        ,FQAI_COV=round(sum((XABCOV/SUBXTOTABCOV)*100*as.numeric(NWCA_CC))/sqrt(length(unique(TAXON))),2)
  )

  ccOut <- reshape2::melt(vascIn.2,id.vars=c(sampID),variable.name='PARAMETER',value.name='RESULT') %>%
    plyr::mutate(PARAMETER=paste(as.character(PARAMETER),'ALL',sep='_'))

  # Now create recoded indicator variables based on NWCA_CC values
  vascIn.alt <- plyr::mutate(vascIn,SEN=ifelse(NWCA_CC %in% c('7','8','9','10'),1,0),ISEN=ifelse(NWCA_CC %in% c('5','6'),1,0)
                          ,TOL=ifelse(NWCA_CC %in% c('4','3','2','1','0'),1,0),HTOL=ifelse(NWCA_CC %in% c('0','1','2'),1,0)
                          ,HSEN=ifelse(NWCA_CC %in% c('9','10'),1,0))

  multTraits <- int.combTraits(vascIn.alt,c('SEN','TOL','ISEN','HTOL','HSEN'),sampID)

  ccOut <- rbind(ccOut,multTraits)

  # Now, if NATSTAT_ALT available, calculate additional metrics
  if('NWCA_NATSTAT' %in% names(vascIn)){
    vascIn.nat <- subset(vascIn,NWCA_NATSTAT=='NAT')

    totals.nat <- plyr::ddply(vascIn.nat,c(sampID),summarise,SUBTOTFREQ=sum(FREQ),SUBXTOTABCOV=sum(XABCOV))
    vascIn1.nat <- merge(subset(vascIn.nat,select=names(vascIn.nat) %nin% c('TOTFREQ','XTOTABCOV')),totals.nat,by=sampID)
    vascIn1a.nat <- plyr::ddply(subset(vascIn1.nat,NWCA_CC!='Und'),c(sampID),summarise
                              ,XC=round(sum(as.numeric(NWCA_CC))/length(unique(TAXON)),2)
                          ,FQAI=round(sum(as.numeric(NWCA_CC))/sqrt(length(unique(TAXON))),2)
                          ,XC_FREQ=round(sum((FREQ/SUBTOTFREQ)*100*as.numeric(NWCA_CC))/length(unique(TAXON)),2)
                          ,FQAI_FREQ=round(sum((FREQ/SUBTOTFREQ)*100*as.numeric(NWCA_CC))/sqrt(length(unique(TAXON))),2)
                          ,XC_COV=round(sum((XABCOV/SUBXTOTABCOV)*100*as.numeric(NWCA_CC))/length(unique(TAXON)),2)
                          ,FQAI_COV=round(sum((XABCOV/SUBXTOTABCOV)*100*as.numeric(NWCA_CC))/sqrt(length(unique(TAXON))),2))

    vascIn1b.nat <- reshape2::melt(vascIn1a.nat,id.vars=c(sampID),variable.name='PARAMETER',value.name='RESULT') %>%
      plyr::mutate(PARAMETER=paste(as.character(PARAMETER),'NAT',sep='_'))

    ccOut <- rbind(ccOut,vascIn1b.nat)

  }

  outdf <- plyr::mutate(ccOut,RESULT=ifelse(is.na(RESULT)|is.infinite(RESULT),0,RESULT))

  return(outdf)
}


# Metrics using only native status
#' @export
#' 
#' @title Calculate metrics based only on native status
#' 
#' @description This function calculates all metrics based
#' only on native status.
#' 
#' @param vascIn Data frame containing cover data summarized by
#' UID and TAXON, with the following fields:
#' \itemize{
#'     \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#'     \item TAXON: Taxon name
#'
#'     \item XABCOV: Mean percent cover of taxon across plots
#'
#'     \item TOTN: Number of taxa in sample
#'
#'     \item sXRCOV: proportion of summed cover across all taxa
#'     (XTOTABCOV) represented by taxon in sample
#'
#'     \item sRFREQ: Relative frequency of a taxon, calculated
#'     as the percentage of the total frequency of taxon
#'     occurrence across all taxa for a UID
#'
#'     \item NWCA_NATSTAT: Native status variable with
#'     categories of 'NAT','ADV','CRYP','INTR','UND'
#'  }
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default

#' @return     Data frame containing \emph{sampID} variables, PARAMETER, RESULT,
#'   where values of PARAMETER consist of the metric name concatenated with
#'   trait value (represented as TRAITNM below):
#' \itemize{
#'   \item PCTN_TRAITNM: Number of taxa with trait as percentage of \emph{TOTN}
#'
#'   \item XABCOV_TRAITNM: Sum of \emph{XABCOV} values across taxa with trait
#'
#'   \item XRCOV_TRAITNM: Sum of \emph{sXRCOV} values across taxa with trait
#'
#'   \item RFREQ_TRAITNM: Sum of \emph{sRFREQ} values across taxa with trait value
#'
#'   \item RIMP_TRAITNM: Relative importance ((RFREQ_TRAITVAL + XRCOV_TRAITVAL)/2)
#'   of taxa with trait value
#'   }
#' 
#' @author Karen Blocksom \email{Blocksom.karen@epa.gov}
#' 
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC.
#' 
#' @examples
#' head(VascPlantEx)
#' exPlant <- prepareData(VascPlantEx)
#'
#' natEx <- calcCC(exPlant$byUIDspp)
#'
#' head(natEx)
#' unique(natEx$PARAMETER)

calcNative <- function(vascIn,sampID='UID'){
  if('NWCA_NATSTAT' %nin% names(vascIn)){
    print("Missing NWCA_NATSTAT from input data frame - cannot calculate metrics! If NWCA_NATSTAT exists,
          run prepareData() function to create input data frame.")
    return(NULL)
  }

  vascIn <- plyr::mutate(vascIn,ALIEN=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),1,0)
                          ,AC=ifelse(NWCA_NATSTAT %in% c('INTR','ADV','CRYP'),1,0))

  sppNATSTAT <- int.calcTraits_MultCat.alt(vascIn,'NWCA_NATSTAT',sampID)

  alienTrait <- int.calcTraits_Indicator.alt(vascIn,'ALIEN',sampID) %>% plyr::mutate(PARAMETER=paste(PARAMETER,'SPP',sep=''))
  acTrait <- int.calcTraits_Indicator.alt(vascIn,'AC',sampID)

  natstatOut <- rbind(sppNATSTAT,alienTrait,acTrait)

  return(natstatOut)
}


#' @export
#' 
#' @title Calculate richness metrics
#' 
#' @description This function calculates richness metrics using plot- and
#' UID-based datasets at the species, genus, and family levels.
#' 
#' @param byUIDspp Data frame containing species data summarized to the UID
#' and TAXON level, with the following variables:
#' \itemize{
#'  \item sampID: Variables identified by \emph{sampID} argument
#'
#' \item STATE: State two-letter abbreviation for site
#'
#' \item USAC_REGION: U.S. Army Corps of Engineers region name for site
#'
#' \item TAXON: Taxon name
#'
#' \item COVER: Total percent cover for taxon
#'
#' \item DISTINCT: Distinctness (0/1) of taxon within sample, where
#' taxa above species level are not distinct (0) if taxa
#' at a lower level are also included in sample. 
#' }
#' @param byPlotspp Data frame containing species data summarized
#' to the UID, PLOT, and TAXON level, with the following variables:
#' \itemize{
#'  \item sampID: Variables identified in sampID argument
#'
#' \item PLOT: Plot from which data were collected
#'
#' \item STATE: State two-letter abbreviation for site
#'
#' \item USAC_REGION: U.S. Army Corps of Engineers region name for site
#'
#' \item TAXON: Taxon name
#'
#' \item COVER: Total percent cover for taxon
#'
#' \item DISTINCT: Distinctness (0/1) of taxon within sample, where
#' taxa above species level are not distinct (0) if taxa
#' at a lower level are also included in sample.
#' }
#' @param byUIDgen Data frame containing genus-level data summarized to the
#'   sampID variables and TAXON level and containing the same variables as 
#'   byUIDspp (described above).
#' @param byPlotgen Data frame containing genus-level data summarized to the
#'   sampID variables, PLOT, and TAXON level and containing the same variables
#'   as byPlotspp (described above).
#' @param byUIDfam  Data frame containing family-level data summarized to the
#'   sampID variables and TAXON level and containing the same variables as
#'   byUIDspp (described above).
#' @param byPlotfam Data frame containing family-level data summarized to the
#'   sampID variables, PLOT, and TAXON level and containing the same variables
#'   as byPlotspp (described above).
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#' 
#' @details The prepareData() function creates a list object with
#' all of the necessary input data frames. For each taxonomic level,
#' the function createDFs() creates a list with a data frame summarized
#' by UID and one by sampID variables and PLOT.
#' 
#' @return   Data frame containing \emph{sampID} variables, PARAMETER, and
#'   RESULT, with one row of results per parameter and \emph{sampID}. The values
#'   for PARAMETER consist of the metric name concatenated with taxonomic level
#'   (represented as SPP, GEN, and FAM below):
#' \itemize{
#' \item TOTN_SPP, TOTN_GEN, TOTN_FAM: Number of unique taxa in sample (UID)
#'
#' \item XN_SPP, XN_GEN, XN_FAM: Mean number of taxa per plot
#'
#' \item MEDN_SPP, MEDN_GEN, MEDN_FAM: Median number of taxa per plot
#'
#' \item SDN_SPP, SDN_GEN, SDN_FAM: Standard deviation of number of taxa per plot
#'
#' \item N_PLOTS: Number of plots sampled for UID
#' }
#' 
#' If NWCA_NATSTAT is present in the input data frame, the following
#' metrics are calculated. GRP in the metric names below represents
#' native status values of NAT, ADV, CRYP, INTR, and subsets ALIEN
#' (INTR + ADV) and AC (ALIEN and CRYP):
#' \itemize{
#' \item TOTN_GRP: Number of unique taxa in sample (UID)
#'
#' \item XN_GRP: Mean number of taxa per plot
#'
#' \item MEDN_GRP: Median number of taxa per plot
#'
#' \item SDN_GRP: Standard deviation of number of taxa per plot
#' }
#' @references US Environmental Protection Agency. 2016. National
#' Wetland Condition Assessment: 2011 Technical Report. EPA-843-R-15-006.
#' US Environmental Protection Agency, Washington, DC.
#' 
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#'   head(VascPlantEx)
#'   exPlant <- prepareData(VascPlantEx)
#'
#'   richEx <- calcRichness(exPlant$byUIDspp,exPlant$byPlotspp,
#'          exPlant$byUIDgen,exPlant$byPlotgen,exPlant$byUIDfam,exPlant$byPlotfam)
#'
#'   head(richEx)
#'   unique(richEx$PARAMETER)

calcRichness <- function(byUIDspp,byPlotspp,byUIDgen,byPlotgen,byUIDfam,byPlotfam,sampID='UID'){

  sppRich <- int.calcRich(byUIDspp,byPlotspp,'SPP',sampID)
  genRich <- int.calcRich(byUIDgen,byPlotgen,'GEN',sampID)
  famRich <- int.calcRich(byUIDfam,byPlotfam,'FAM',sampID)

  richOut <- rbind(sppRich,genRich,famRich)

  if('NWCA_NATSTAT' %in% names(byUIDspp)){
    natRich <- int.calcRichNS(byUIDspp,byPlotspp,c('NAT'),'NATSPP',sampID)
    advRich <- int.calcRichNS(byUIDspp,byPlotspp,c('ADV'),'ADVSPP',sampID)
    crypRich <- int.calcRichNS(byUIDspp,byPlotspp,c('CRYP'),'CRYPSPP',sampID)
    intrRich <- int.calcRichNS(byUIDspp,byPlotspp,c('INTR'),'INTRSPP',sampID)
    alienRich <- int.calcRichNS(byUIDspp,byPlotspp,c('ADV','INTR'),'ALIENSPP',sampID)
    acRich <- int.calcRichNS(byUIDspp,byPlotspp,c('INTR','ADV','CRYP'),'AC',sampID)

    # Combine all into a single df
    allNSrich <- rbind(natRich,advRich,crypRich,intrRich,alienRich,acRich)

    richOut <- rbind(richOut,allNSrich)
  }
    # Must fill in missing categories for all sites with zeros
  formula <- paste(paste(sampID,collapse='+'),'~PARAMETER',sep='')
  allRichOut <- reshape2::dcast(richOut,eval(formula),value.var='RESULT') %>%
    reshape2::melt(id.vars=sampID,variable.name='PARAMETER',value.name='RESULT') %>%
    plyr::mutate(RESULT=ifelse(is.na(RESULT),0,RESULT),PARAMETER=as.character(PARAMETER))

  return(allRichOut)
}


#' @export
#' 
#' @title Calculate diversity indices
#' 
#' @description Calculate Simpson Diversity, Shannon-Wiener
#' Diversity, and Pielou's Evenness indices, with variations based on
#' native status if NWCA_NATSTAT is included in the input data frame.
#' 
#' @param vascIn Data frame containing cover data summarized by
#' UID and TAXON, with the following fields:
#' \itemize{
#'  \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#' \item TAXON: Taxon name
#'
#' \item CATEGORY: USDA PLANTS category variable
#'
#' \item XABCOV: Mean percent cover of taxon across plots
#'
#' \item Optional: NWCA_NATSTAT: Native status variable with categories of
#' 'NAT','ADV','CRYP','INTR','UND'
#' }
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#' 
#' @return Data frame containing \emph{sampID} variables, PARAMETER, RESULT,
#'   where values of PARAMETER are:
#' \itemize{
#' \item D_ALL: Simpson diversity index based on all taxa
#'
#' \item H_ALL: Shannon-Wiener diversity index based on all taxa
#'
#' \item J_ALL: Pielou's evenness index based on all taxa
#' }
#' If NWCA_NATSTAT is included in input data frame, the
#' following variations based on subsets of taxa are
#' calculated:
#' \itemize{
#' \item D_NAT: Simpson diversity based only on native species
#'
#' \item D_ALIEN: Simpson diversity based only on alien species
#'
#' \item D_AC: Simpson diversity based only on alien and cryptogenic
#' species
#'
#' \item H_NAT: Shannon-Wiener diversity index based on native species
#'
#' \item H_ALIEN: Shannon-Wiener diversity index based on alien species
#'
#' \item H_AC: Shannon-Wiener diversity index based on alien and
#' cryptogenic species
#'
#' \item J_NAT: Pielou's evenness index based on native taxa
#'
#' \item J_ALIEN: Pielou's evenness index based on alien taxa
#'
#' \item J_AC: Pielou's evenness index based on alien and cryptogenic taxa
#' }
#' 
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC.
#' 
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(VascPlantEx)
#' exPlant <- prepareData(VascPlantEx)
#'
#' divEx <- calcDiversity(exPlant$byUIDspp)
#'
#' head(divEx)
#' unique(divEx$PARAMETER)

calcDiversity <- function(vascIn,sampID='UID'){
  divOut <- int.calcIndices(vascIn,'ALL',sampID)

 if('NWCA_NATSTAT' %in% names(vascIn)){

   vascIn.1 <- plyr::mutate(vascIn,NATSTAT_ALT=ifelse(NWCA_NATSTAT %in% c('INTR','ADV'),'ALIEN',NWCA_NATSTAT)
                          ,AC=ifelse(NWCA_NATSTAT %in% c('INTR','ADV','CRYP'),1,0))

   nsvalues <- c('NAT','ALIEN')
   for(i in 1:length(nsvalues)){
     vascIn.2 <- subset(vascIn.1,NATSTAT_ALT==nsvalues[i])

     nsOut <- int.calcIndices(vascIn.2,nsvalues[i],sampID)
     divOut <- rbind(divOut,nsOut)
   }

   # AC metrics
   vascIn.ac <- subset(vascIn.1,AC=='1')
   acOut <- int.calcIndices(vascIn.ac,'AC',sampID)

   divOut <- rbind(divOut,acOut)

 }
 formula <- paste(paste(sampID,collapse='+'),'~PARAMETER',sep='')
 dfOut <- reshape2::dcast(divOut,eval(formula),value.var='RESULT') %>%
   reshape2::melt(id.vars=sampID,variable.name='PARAMETER',value.name='RESULT') %>%
   plyr::mutate(dfOut,RESULT=ifelse(is.na(RESULT),0,RESULT),PARAMETER=as.character(PARAMETER))

 return(dfOut)
}


#' @export
#' 
#' @title Calculate Bray-Curtis metrics
#' 
#' @description This function calculates the mean Bray-Curtis
#' distances among plots for all species, and includes a version using only
#' native species if the variable NWCA_NATSTAT is included in the input
#' data frame. This variable is found in the ccNatNWCA dataset.
#' 
#' @param vascIn Data frame containing cover data summarized by 
#' \emph{sampID} variables, PLOT, and
#' DISTINCT at the species level. Must also contain at least one of the
#' following: USDA_NAME (taxon name) or SPECIES_NAME_ID (numeric code
#' for taxon). If NWCA_NATSTAT is included, a native species version is
#' also calculated.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#' 
#' @details This function calculates metrics based on all species, and on
#' native species if the appropriate variable is included in the
#' input data frame.
#' 
#' @return Data frame containing \emph{sampID} variables, PARAMETER, and RESULT,
#'   with one row of results per parameter and UID. The values for PARAMETER
#'   consist of the metric name concatenated with taxonomic level (represented
#'   as TAXLEVEL below):
#' \itemize{
#' \item XBCDIST_SPP: Mean Bray-Curtis distance among plots in sample
#'
#' \item XBCDIST_NATSPP: Mean Bray-Curtis distance among plots in sample
#' based only on native species.
#' }

calcBCmets <- function(vascIn,sampID='UID'){
  # Need to account for cases where there is no SPECIES_NAME_ID by creating one just for this
  # calculation
  if('SPECIES_NAME_ID' %nin% names(vascIn)){
    uniqNames <- unique(vascIn$USDA_NAME)
    uniqNames <- plyr::mutate(uniqNames,SPECIES_NAME_ID=seq(1,length(uniqNames)))
    vascIn <- merge(vascIn,uniqNames,by='UsDA_NAME')
  }

  forDist <- plyr::ddply(vascIn,c(sampID,'PLOT','SPECIES_NAME_ID')
                         ,summarise,COVER=sum(as.numeric(COVER)),.progress='tk')
  # This df needs to be in wide format
  forDist <- plyr::mutate(forDist,SPECIES=paste('s',SPECIES_NAME_ID,sep=''))

  meanBC <- int.calcXBC(forDist,sampID)

  xbcOut <- meanBC

  if('NWCA_NATSTAT' %in% names(vascIn)){
    forDist.nat <- plyr::ddply(vascIn,c(sampID,'PLOT','SPECIES_NAME_ID','NWCA_NATSTAT')
                               ,summarise,COVER=sum(as.numeric(COVER)),.progress='tk') %>%
      plyr::mutate(SPECIES=paste('s',SPECIES_NAME_ID,sep='')) %>%
      dplyr::filter(NWCA_NATSTAT=='NAT')

    meanBC_nat <- int.calcXBC(forDist.nat,sampID) %>% plyr::rename(c('XBCDIST_SPP'='XBCDIST_NATSPP'))

    xbcOut <- merge(xbcOut,meanBC_nat,by='UID',all.x=T)
  }

  xbcOut.1 <- reshape2::melt(xbcOut,id.vars=sampID,variable.name='PARAMETER',value.name='RESULT') %>%
    plyr::mutate(RESULT=ifelse(is.na(RESULT),0,RESULT),PARAMETER=as.character(PARAMETER))

  return(xbcOut.1)
}

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
#' exPlant <- prepareData(VascPlantEx)
#'
#' vmmiEx <- calcVMMImets(exPlant$byUIDspp)
#'
#' head(vmmiEx)
#' unique(vmmiEx$PARAMETER)

calcVMMImets <- function(vascIn,sampID='UID'){
  # Calculate only the 4 metrics used in the VMMI: FQAI_ALL,N_TOL,RIMP_NATSPP,XRCOV_MONOCOTS_NAT
  # First obtain all of the UIDs in input data frame
  for(i in 1:length(sampID)){
    if(i==1) vascIn$SAMPID <- vascIn[,sampID[i]]
    else vascIn$SAMPID <- paste(vascIn$SAMPID,vascIn[,sampID[i]],sep='.')
  }
  samples <- unique(subset(vascIn,select=c(sampID,'SAMPID')))
    
  UIDs <- data.frame(SAMPID=unique(subset(vascIn,select='SAMPID')),stringsAsFactors=FALSE)

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
  fqaiOut <- plyr::ddply(subset(vascIn,NWCA_CC!='Und'),'SAMPID',summarise,PARAMETER='FQAI_ALL'
                      ,RESULT=round(sum(as.numeric(NWCA_CC))/sqrt(length(unique(TAXON))),2))

  # Calculate N_TOL
  vascIn.alt <- plyr::mutate(vascIn,TOL=ifelse(NWCA_CC %in% c('4','3','2','1','0'),1,0)) %>%
    subset(TOL==1)

  if(nrow(vascIn.alt)>0){
    ntol <- plyr::ddply(vascIn.alt,'SAMPID',summarise,N_TOL=length(TAXON))
    ntolOut <- merge(UIDs,ntol,by='SAMPID',all.x=T) %>%
      reshape2::melt(id.vars='SAMPID',variable.name='PARAMETER',value.name='RESULT') %>%
      plyr::mutate(RESULT=ifelse(is.na(RESULT),0,RESULT))

  }else{
    ntolOut <- data.frame(SAMPID=UIDs$SAMPID,PARAMETER='N_TOL',RESULT=0,stringsAsFactors=F)
  }

  # Calculate RIMP_NATSPP
  vascIn.nat <- subset(vascIn,NWCA_NATSTAT=='NAT')

  if(nrow(vascIn.nat)>0){
    vascIn.nat.1 <- plyr::ddply(vascIn.nat,'SAMPID',plyr::summarise,XRCOV=round(sum(sXRCOV),2)
                              ,RFREQ=round(sum(sRFREQ),2),RIMP_NATSPP=round((RFREQ+XRCOV)/2,2))

    natOut <- merge(UIDs,vascIn.nat.1,by='SAMPID',all.x=T) %>%
      reshape2::melt(id.vars='SAMPID',measure.vars='RIMP_NATSPP',variable.name='PARAMETER'
                     ,value.name='RESULT') %>%
      plyr::mutate(RESULT=ifelse(is.na(RESULT),0,RESULT))

   }else{
    natOut <- data.frame(SAMPID=UIDs$SAMPID,PARAMETER='RIMP_NATSPP',RESULT=0,stringsAsFactors=F)
  }

  # Calculate XRCOV_MONOCOTS_NAT
  vascIn.mono <- subset(vascIn,CATEGORY=='MONOCOT' & NWCA_NATSTAT=='NAT')

  if(nrow(vascIn.mono)>0){
    vascIn.mono.1 <- plyr::ddply(vascIn.mono,c('SAMPID'),plyr::summarise,XRCOV_MONOCOTS_NAT=round(sum(sXRCOV),2))

    monoOut <- merge(UIDs,vascIn.mono.1,by='SAMPID',all.x=T) %>%
      reshape2::melt(id.vars=c('SAMPID'),measure.vars='XRCOV_MONOCOTS_NAT',variable.name='PARAMETER'
                    ,value.name='RESULT') %>%
      plyr::mutate(RESULT=ifelse(is.na(RESULT),0,RESULT))

  }else{
    monoOut <- data.frame(SAMPID=UIDs$SAMPID,PARAMETER='XRCOV_MONOCOTS_NAT',RESULT=0,stringsAsFactors=F)
  }

 # Now combine into one data frame
 allOut <- rbind(fqaiOut,natOut,monoOut,ntolOut)
 allOut.1 <- reshape2::dcast(allOut,SAMPID~PARAMETER,value.var='RESULT') %>%
   merge(samples,by='SAMPID') %>% 
   dplyr::select(-SAMPID)
 
 allOut.1 <- allOut.1[,c(sampID,'FQAI_ALL','N_TOL','RIMP_NATSPP','XRCOV_MONOCOTS_NAT')]

 return(allOut.1)

}



