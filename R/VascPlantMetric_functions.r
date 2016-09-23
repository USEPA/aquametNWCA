#' @export
#' @title Create data frames for function input
#' @description This internal function merges the inputs of
#' plant cover data and taxalist and summarizes by the taxonomic
#' level specified. Not intended for use on its own.
#' @param tvar String with the level of taxonomy
#' ('USDA_NAME','GENUS','FAMILY')
#' @param indf Data frame with vegetation cover data, having
#' the following variables:
#' \itemize{
#'
#' \item sampID - A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#'
#' \item PLOT - Plot number of data (1 to 5 possible)
#'
#' \item STATE - State postal code of site
#'
#' \item USAC_REGION - U.S. Army Corps of Engineers Region Name
#'
#' \item USDA_NAME - Taxon name, must match with taxa data frame
#'
#' \item COVER - Percentage estimated vegetation cover by TAXON in PLOT
#' }
#' @param taxa Data frame containing USDA_NAME, GENUS, CATEGORY,
#' GROWTH_HABIT, and DURATION variables. Dataset taxaNWCA is the
#' default if not specified. These variables are assumed to be
#' populated with abbreviations as found in UsDA PLANTS database.
#' @return A list containing two data frames. The first data frame
#' summarizes data by UID and TAXON and contains UID, STATE,
#' USAC_REGION, TAXON, NUM, XABCOV, DISTINCT, and NPLOTS. NUM is the
#' number of plots in which taxon occurs, and XABCOV is the mean
#' absolute COVER across plots. DISTINCT is the value 1 assigned to
#' each row. NPLOTS is the number of plots sampled.
#'
#' The second summarizes by UID, PLOT, and TAXON. Each data frame
#' contains UID, PLOT, STATE, USAC_REGION, TAXON, COVER, and
#' DISTINCT. DISTINCT assigns the value for each taxon as 1 if the
#' taxon has COVER>0 or 0 if not. COVER is the sum of the COVER
#' variable by UID, PLOT, and TAXON.
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' @examples
#' head(VascPlantEx)
#' data(taxaNWCA)
#'
#' outEx <- createDFs('GENUS',VascPlantEx,taxaNWCA)
#' head(outEx$byUID)
#' head(outEx$byPlot)
createDFs <- function(sampID='UID',tvar,indf,taxa){
  
  # First merge the taxa list with the cover data by USDA_NAME
  indf.1 <- merge(indf,taxa,by='USDA_NAME')
  indf.1$tobj <- indf.1[,tvar] # Set tobj as the value of tvar
  indf.1 <- plyr::mutate(indf.1,COVER=as.numeric(COVER)) # Make sure COVER is numeric
  # Set value for TAXON as either specified taxon level above species, or as USDA_NAME
  byPlot <- plyr::mutate(indf.1,TAXON=ifelse(!is.na(tobj) & tobj!='',tobj,USDA_NAME))
  # Sum COVER by UID, PLOT, and TAXON to ensure there is only one row per species in input data, set DISTINCT as 1 to be
  # taxon counter
  byPlot1 <- plyr::ddply(byPlot,c(sampID,'STATE','USAC_REGION','PLOT','TAXON'),summarise,COVER=sum(COVER)
                         ,DISTINCT=ifelse(COVER>0,1,0),.progress='tk')
  byPlot1$COVER[byPlot1$COVER>100] <- 100 # Cap sums at 100 percent
  print("Done with plot summing")

  ## Calculate frequency and relative frequency variables by taxon
  byUID1 <- plyr::ddply(byPlot1,c(sampID,'STATE','USAC_REGION','TAXON'),summarise,NUM=sum(DISTINCT),XABCOV=mean(COVER)
                        ,NPLOTS=length(unique(PLOT)))
  byUID2 <- plyr::mutate(byUID1,DISTINCT=1)

  print("Done with sAMPID summing")

  byDF <- list(byUID=byUID2,byPlot=byPlot1)
  return(byDF)
}

#' @export
#' @title Prepare data for calculating metrics
#' @description Assemble datasets from data frame containing
#' vegetation cover by UID and PLOT, and taxa lists containing
#' variables used in metric calculation. Species-, genus-, and
#' family-level datasets are summarized by plot and UID. For the
#' species-level output data summarized by UID, various traits
#' are added to the output data set.
#' @param indf Data frame containing cover data summarized at
#' the UID, PLOT, and taxon value, with the following variable
#' names:
#'  \itemize{
#' \item sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#'
#' \item PLOT: Plot number of data
#'
#' \item USDA_NAME: Taxon name, based on USDA PLANTS names,
#' supplemented with NWCA names, if not in USDA PLANTS
#'
#' \item COVER: Percent cover of taxon in plot, including zeros
#' for plots in which a taxon does not occur.
#' }
#' @param inTaxa Data frame with all taxa in indf, with variables:
#' \itemize{
#' \item USDA_NAME: Taxon name, based on USDA PLANTS names,
#' supplemented with NWCA names, if necessary.
#'
#' \item FAMILY: Family name of taxon
#'
#' \item GENUS: Genus name of taxon
#'
#' \item CATEGORY (optional): USDA PLANTS category variable,
#' necessary to calculate category metrics.
#'
#' \item DURATION (optional): USDA PLANTS duration variable,
#' necessary to calculate duration metrics.
#'
#' \item GROWTH_HABIT (optional): USDA PLANTS growth habit
#' variable,necessary to calculate growth habit metrics
#'
#' \item SPECIES_NAME_ID (optional): Taxonomic ID number,
#' which will be used in Bray-Curtis distance metrics
#' if available.
#'
#' \item STATE: Postal abbreviation for state of sample
#'
#' \item USAC_REGION: U.S. Army Corps of Engineers region
#' abbreviation for sample, to correspond to those
#' in inWIS data frame.
#' }
#' @param inNatCC Data frame with C-values and native status:
#'
#' \itemize{
#' \item USDA_NAME: Taxon name
#'
#' \item GEOG_ID: Postal abbreviation for STATE of taxon
#'
#' \item NWCA_CC (optional): Coefficient of conservatism,
#' as used in NWCA, necessary to calculate CC metrics
#'
#' \item NWCA_NATSTAT (optional): Native status variable,
#' as used in NWCA, necessary to calculate native status metrics.
#' }
#' @param inWIS Data frame with Wetland Indicator Status, from
#' U.S. Army Corps of Engineers (USAC):
#' \itemize{
#' \item USDA_NAME: Taxon name
#'
#' \item GEOG_ID: USAC region, abbreviated to match those used in input
#' data frame
#'
#' \item WIS: Wetland Indicator Status as provided by USAC
#' or added for NWCA
#' }
#' @details This function calls the createDFs() function, which
#' sums cover by UID, PLOT, TAXON, with sums > 100 truncated to
#' 100 percent.
#' @return A list containing six data frames:
#' \itemize{
#'  \item
#'    byUIDspp: Data frame with data summarized by sampID variables and TAXON
#'    at the species level and contains:
#'    \itemize{
#'      \item sampID: A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#'
#'      \item STATE: State abbreviation of site location, necessary for merging
#' with inNatCC
#'
#'      \item USAC_REGION: U.S. Army Corps of Engineers region abbreviation
#' necessary for merging with inWIS.
#'
#'      \item TAXON: Taxon name
#'
#'      \item NUM: Number of occurrences of taxon across plots
#'
#'      \item XABCOV: Mean percent absolute cover of taxon across plots
#'
#'      \item DISTINCT: Distinctness value for each taxon, 1 if
#' the taxon has COVER>0 and 0 if not.
#'
#'      \item NPLOTS: Number of plots in sample (1-5)
#'
#'      \item TOTN: Total number of taxa in sample
#'
#'      \item XTOTABCOV: Sum of XABCOV across all taxa in sample
#'
#'      \item sXRCOV: taxon mean relative cover (XABCOV/XTOTABCOV)*100
#'
#'      \item FREQ: Relative number of plots in which taxon occurs
#' (NUM/NPLOTS)*100
#'
#'      \item TOTFREQ: Sum of FREQ across all taxa in sample
#'
#'      \item SRFREQ: Relative frequency of taxon relative to total
#' frequency (FREQ/TOTFREQ)*100
#'    }
#' This data frame is also merged with the input taxa data
#' frames and contains, in addition, GENUS and FAMILY, but
#' also contains (depending on the input taxa lists): CATEGORY,
#' DURATION, GROWTH_HABIT, NWCA_CC, NWCA_NATSTAT, WIS.
#'
#' \item byPlotspp: Data frame with data summarized by UID, PLOT, and
#' TAXON at the species level. Each data frame contains:
#'    \itemize{
#'      \item sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#'
#'      \item PLOT: Plot number
#'
#'      \item STATE: State of sample location
#'
#'      \item USAC_REGION: USAC region code
#'
#'      \item TAXON: Taxon name
#'
#'      \item COVER: Sum of cover by TAXON within plot
#'
#'      \item DISTINCT: Distinctness of taxon, value of 1 assigned to each row
#' }
#' 
#' \item byUIDgen: Data frame with data summarized by UID and TAXON
#' at the genus level and contains \emph{sampID}, STATE, USAC_REGION,
#' TAXON, NUM, XABCOV, and DISTINCT. NUM is the number of
#' plots in which taxon occurs, and XABCOV is the mean
#' absolute COVER across plots. DISTINCT is the value 1
#' assigned to each row.
#'
#' \item byPlotgen: Data frame with data summarized by UID, PLOT, and
#' TAXON at the genus level. Each data frame contains \emph{sampID},
#' PLOT, STATE, USAC_REGION, TAXON, COVER, and DISTINCT.
#' DISTINCT assigns the value for each taxon as 1 if the
#' taxon has COVER>0 or 0 if not. COVER is the sum of the
#' COVER variable.
#'
#' \item byUIDfam: Data frame with data summarized by UID and TAXON
#' at the family level and contains \emph{sampID}, STATE, USAC_REGION,
#' TAXON, NUM, XABCOV, and DISTINCT. NUM is the number of plots
#' in which taxon occurs, and XABCOV is the mean absolute COVER
#' across plots. DISTINCT is the value 1 assigned to each row.
#'
#' \item byPlotfam: Data frame with data summarized by \emph{sampID}
#' , PLOT, and TAXON at the family level. Each data frame 
#' contains \emph{sampID}, PLOT, STATE, USAC_REGION, TAXON, 
#' COVER, and DISTINCT. DISTINCT assigns the value for each 
#' taxon as 1 if the taxon has COVER>0 or 0 if not. COVER is the 
#' sum of the COVER variable.
#' }
#' @references US Environmental Protection Agency. 2016.
#' National Wetland Condition Assessment: 2011 Technical Report.
#' EPA-843-R-15-006. US Environmental Protection Agency,
#' Washington, DC.
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' @examples
#' head(VascPlantEx)
#' prepEx <- prepareData(VascPlantEx)
#'
#' str(prepEx)
prepareData <- function(indf,sampID='UID',inTaxa=taxaNWCA,inNatCC=ccNatNWCA,inWIS=wisNWCA){
  # Read in various input datasets, and create output dataset based on available
  # types of data - must have cover data and taxonomic data at the very least
  # If
  datNames <- c(sampID,'PLOT','USDA_NAME','COVER','STATE','USAC_REGION')
  if(any(datNames %nin% names(indf))){
    print(paste("Missing key variables! Should be ",sampID," PLOT, USDA_NAME, and COVER.",sep='')
    return(NULL)
  }
  # Input taxa list with taxonomy and life history traits
  taxNames <- c('USDA_NAME','FAMILY','GENUS')
  altNames <- c('CATEGORY','GROWTH_HABIT','DURATION')
  if(any(taxNames %nin% names(inTaxa))){
    print("Missing key variables! Need at least USDA_NAME, FAMILY, GENUS to calculate metrics.")
    return(NULL)
  }
  if(any(altNames %nin% names(inTaxa))){
    msgNames <- altNames[altNames %nin% names(inTaxa)]
    print(paste("Will not be able to calculate metrics using ",paste(msgNames,collapse=','),
                " without these parameters in inTaxa",sep=''))
  }

  # Coefficients of conservatism and native status values
    ccNames <- c('USDA_NAME','GEOG_ID')
    ccVars <- c('NWCA_CC','NWCA_NATSTAT')
  # This only applies if someone specifies a taxalist not included in the package
    if(any(ccNames %nin% names(inNatCC))){
      print("Missing key variables! Need variables named USDA_NAME, GEOG_ID (State) to match up
            with cover data. This taxa list cannot be used in calculations. Either revise file
            or use default taxa list by not specifying the inNatCC argument.")
      return(NULL)
    }
    if(any(ccVars %nin% names(inNatCC))){
      msgNames <- ccVars[ccVars %nin% names(inNatCC)]
      print(paste("Warning: Will not be able to calculate metrics using ",paste(msgNames,collapse=','),
                  " without these parameter in inNatCC."))
    }


  # Wetland Indicator Status values
  if(!is.null(inWIS)){
    wisNames <- c('USDA_NAME','GEOG_ID','WIS')
    if(any(wisNames %nin% names(inWIS))){
      print("Missing key variables! Need USDA_NAME, GEOG_ID, WIS to match up with cover data.
            This taxa list cannot be used in calculations. Either revise file or use default
            taxa list by not specifying inWIS argument.")
      return(NULL)
    }
  }

  ## Create dfs for species level, genus, family, and order
  # First construct list object with by plot and by sampID summarizations
  dfSPP <- createDFs('USDA_NAME',indf,inTaxa)
  # For species-level data, run additional checks and add additional information
  # Merge cover data with taxalist
  dfSPP.byUID.1a <- merge(dfSPP[[1]],inTaxa,by.x='TAXON',by.y='USDA_NAME')

  # If any taxa in the cover data do not match up with the taxalist, return
  # missing names and end function
  if(nrow(dfSPP.byUID.1a)!=nrow(dfSPP[[1]])){
    print("Not all taxa in dfIn match up with names in taxaIn!")
    check1 <- merge(dfSPP.byUID.1a,dfSPP[[1]],by=c(sampID,'TAXON'),all.y=T)
    checkout <- unique(subset(check1,is.na(SPECIES_NAME_ID),select=c('TAXON')))
    print(checkout)
    return(NULL)
  }

  # If all taxa match taxalist, merge now with CC/native status by STATE
  if(!is.null(inNatCC)){
    dfSPP.byUID.1b <- merge(dfSPP.byUID.1a,inNatCC,by.x=c('TAXON','STATE'),by.y=c('USDA_NAME','GEOG_ID'))
  }else{
    dfSPP.byUID.1b <- dfSPP.byUID.1a
  }
  # Merge with WIS values by USAC_REGION
  if(!is.null(inWIS)){
    dfSPP.byUID.1c <- merge(dfSPP.byUID.1b,inWIS,by.x=c('TAXON','USAC_REGION'),by.y=c('USDA_NAME','GEOG_ID'),all.x=T)
  }else{
    dfSPP.byUID.1c <- dfSPP.byUID.1b
  }

  # Calculate totals and add them to output data frame
  dfSPP.byUID.fin <- plyr::ddply(dfSPP.byUID.1c,sampID,mutate,TOTN=length(TAXON)
                                 ,XTOTABCOV=sum(XABCOV),TOTFREQ=sum(NUM/NPLOTS)*100) %>%
    plyr::mutate(sXRCOV=XABCOV/XTOTABCOV*100, FREQ=(NUM/NPLOTS)*100
                 ,sRFREQ=(FREQ/TOTFREQ)*100)

  dfSPP[[1]] <- dfSPP.byUID.fin

  ## Also want to add NWCA_NATSTAT to dfSPP[[2]], byPlot
  dfSPP.byPlot <- merge(dfSPP[[2]],inNatCC,by.x=c('STATE','TAXON'),by.y=c('GEOG_ID','USDA_NAME'))
  dfSPP[[2]] <- dfSPP.byPlot

  # Create datasets for genus and family levels which will only be used for richness metrics
  dfGEN <- createDFs('GENUS',indf,inTaxa)
  dfFAM <- createDFs('FAMILY',indf,inTaxa)

  outDF <- list(byUIDspp=dfSPP[[1]],byPlotspp=dfSPP[[2]],byUIDgen=dfGEN[[1]],byPlotgen=dfGEN[[2]],byUIDfam=dfFAM[[1]],byPlotfam=dfFAM[[2]])
  print("Done preparing datasets")
  return(outDF)
}

#' @export
#' @title Calculate vascular plant richness metrics
#' @description This internal function calculates taxa richness
#' of sample at a specified level of taxonomy, based on output
#' from CreateDFs(). Not intended for use on its own.
#' @param x Data frame containing cover data summarized by sampID,
#' TAXON, DISTINCT at the specified taxonomic level (tlevel)
#' @param y Data frame containing cover data summarized by sampID,
#' PLOT, TAXON, and DISTINCT at the specified taxonomic level.
#' @param tlevel Taxonomic level of input dataset, abbreviated
#' to serve as suffix to metric names ('SPP','GEN','FAM')
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @return Data frame containing UID, PARAMETER, and RESULT,
#' with one row of results per parameter and UID. The values for
#' PARAMETER consist of the metric name concatenated with taxonomic
#' level (represented as TAXLEVEL below):
#'
#' TOTN_TAXLEVEL: Number of unique taxa in sample 
#'
#' XN_TAXLEVEL: Mean number of taxa per plot
#'
#' MEDN_TAXLEVEL: Median number of taxa per plot
#'
#' SDN_TAXLEVEL: Standard deviation of number of taxa per plot
#'
#' N_PLOTS: Number of plots sampled for sampID
#'
#' @author Karen Blocksom
.calcRich <- function(x,y,tlevel,sampID='UID') {
  xx1 <- plyr::ddply(x,sampID,summarise,TOTN_TAXA=sum(DISTINCT))
  ## Now calculate richness by plot to obtain average richness per plot
  yy1 <- merge(subset(y,select=c(sampID,'PLOT','DISTINCT')),xx1,by=sampID)
  yy2 <- plyr::ddply(yy1,c(sampID,'PLOT','TOTN_TAXA'),summarise,N_TAXA=sum(DISTINCT))
  yy3 <- plyr::ddply(yy2,c(sampID,'TOTN_TAXA'),summarise,XN_TAXA=round(mean(N_TAXA),2),MEDN_TAXA=median(N_TAXA)
               ,SDN_TAXA=round(sd(N_TAXA),2))

  outdf <- reshape2::melt(yy3,id.vars=sampID,variable.name='PARAMETER',value.name='RESULT')
  outdf$PARAMETER <- gsub('TAXA',tlevel,outdf$PARAMETER)
  if(tlevel!='SPP'){
    outdf <- subset(outdf,PARAMETER!='NPLOTS')
  }
  return(outdf)
}

#' @export 
#' @title Calculate metrics for traits with >2 categories
#' @description This internal function calculates metrics using
#' traits with more than two categories as values. Used in
#' calculating metrics for growth habit and duration.
#' @param indf Input data frame with UID, TAXON, TOTN, XABCOV,
#' sXRCOV, and variable with name in argument trait. TOTN is the
#' total number of taxa at the lowest level for a UID. XABCOV is
#' the mean absolute cover of taxon. sXRCOV is the percentage of
#' the total sum of absolute cover across all taxa for a UID.
#' @param trait Character string with name of variable containing
#' traits of interest
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @return Data frame containing UID, PARAMETER, RESULT, where values
#' of PARAMETER consist of the metric name concatenated with each
#' value of the trait (represented as TRAITVAL below):
#'
#' N_TRAITVAL: Number of taxa with trait value
#'
#' PCTN_TRAITVAL: Number of taxa with trait value as percentage of TOTN
#'
#' XABCOV_TRAITVAL: Sum of XABCOV values across taxa with trait value
#'
#' XRCOV_TRAITVAL: Sum of sXRCOV values across taxa with trait value
#' @author Karen Blocksom
.calcTraits_MultCat <- function(indf,trait,sampID='UID'){
  indf1 <- subset(indf,!is.na(eval(as.name(trait))) & eval(as.name(trait))!='')

  indf2 <- plyr::ddply(indf1,c(sampID,trait),summarise,N=length(TAXON),PCTN=round(N/unique(TOTN)*100,2)
                 ,XABCOV=round(sum(XABCOV),2),XRCOV=round(sum(sXRCOV),2))
  outdf <- reshape2::melt(indf2,id.vars=c(sampID,trait))
  
  formula <- paste(paste(paste(sampID,collapse='+'),'~variable',sep=''),trait,sep='+')
  outdf1 <- reshape2::melt(reshape2::dcast(outdf,eval(formula),value.var='value'),id.vars=sampID,variable.name='PARAMETER'
                           ,value.name='RESULT')
  outdf1 <- plyr::mutate(outdf1,RESULT=ifelse(is.na(RESULT),0,RESULT),PARAMETER=as.character(PARAMETER))
  return(outdf1)
}

#' @export 
#' @title Alternate metric calculations for traits with >2 categories
#' @description This internal function calculates metrics using traits
#' having >2 categories as values. Used for native status variables.
#' @param indf Input data frame with UID, TAXON, TOTN, XABCOV, sXRCOV,
#' sRFREQ, and variable with name in argument trait. TOTN is the total
#' number of taxa at the lowest level for a UID. XABCOV is the mean
#' absolute cover of taxon. sXRCOV is the percentage of the total
#' sum of absolute cover across all taxa for a UID. sRFREQ is the relative
#' frequency of a taxon, calculated as the percentage of the total frequency
#' of taxon occurrence across all taxa for a UID.
#' @param trait Character string with name of variable containing
#' traits of interest
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @return Data frame containing UID, PARAMETER, RESULT, where values
#' of PARAMETER consist of the metric name concatenated with each
#' value of the trait (represented as TRAITVAL below):
#'
#' N_TRAITVAL: Number of taxa with trait value
#'
#' PCTN_TRAITVAL: Number of taxa with trait value as percentage of TOTN
#'
#' XABCOV_TRAITVAL: Sum of XABCOV values across taxa with trait value
#'
#' XRCOV_TRAITVAL: Sum of sXRCOV values across taxa with trait value
#'
#' RFREQ_TRAITVAL: Sum of sRFREQ values across taxa with trait value
#'
#' RIMP_TRAITVAL: Relative importance ((RFREQ_TRAITVAL + XRCOV_TRAITVAL)/2)
#' of taxa with trait value
#' @author Karen Blocksom

.calcTraits_MultCat.alt <- function(indf,trait,sampID='UID'){
  indf1 <- subset(indf,!is.na(eval(as.name(trait))) & eval(as.name(trait))!='')

  indf2 <- plyr::ddply(indf1,c(sampID,trait),summarise,PCTN=round(length(TAXON)/unique(TOTN)*100,2)
                 ,XABCOV=round(sum(XABCOV),2),XRCOV=round(sum(sXRCOV),2),RFREQ=round(sum(sRFREQ),2)
                 ,RIMP=round((RFREQ+XRCOV)/2,2))
  outdf <- reshape2::melt(indf2,id.vars=c(sampID,trait))
  
  formula <- paste(paste(paste(sampID,collapse='+'),'~variable',sep=''),trait,sep='+')
  outdf1 <- reshape2::melt(reshape2::dcast(outdf,eval(formula),value.var='value'),id.vars=sampID,variable.name='PARAMETER'
                           ,value.name='RESULT') %>%
    plyr::mutate(PARAMETER=as.character(PARAMETER)) %>%
    dplyr::filter(substring(PARAMETER,nchar(PARAMETER)-3,nchar(PARAMETER))!='_UND') %>%
    plyr::mutate(RESULT=ifelse(is.na(RESULT),0,RESULT)
                         ,PARAMETER=paste(as.character(PARAMETER),'SPP',sep=''))

  return(outdf1)
}

#' @export 
#' @title Calculate metrics using traits with only two values (0/1)
#' @description This internal function calculates metrics for traits
#' that are indicator values (0 or 1). Not intended to be used on its
#' own. Output feeds into combTraits function.
#' @param indf Input data frame with variables UID, TAXON, TOTN, XABCOV,
#' sXRCOV, and a binary variable with name in trait argument
#' and possible values of 0 and 1, with 1 indicating trait
#' present for taxon.
#' @param trait Character string with name of variable containing
#' traits of interest
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @return Data frame containing UID, PARAMETER, RESULT, where values
#' of PARAMETER consist of the metric name concatenated with trait name
#' (represented as TRAITNM below):
#'
#' N_TRAITNM: Number of taxa with trait
#'
#' PCTN_TRAITNM: Number of taxa with trait as percentage of TOTN
#'
#' XABCOV_TRAITNM: Sum of XABCOV values across taxa with trait
#'
#' XRCOV_TRAITNM: Sum of sXRCOV values across taxa with trait
#' @author Karen Blocksom
.calcTraits_Indicator <- function(indf,trait,sampID='UID'){
  
  for(i in 1:length(sampID)){
    if(i==1) inCts$SAMPID <- inCts[,sampID[i]]
    else inCts$SAMPID <- paste(inCts$SAMPID,inCts[,sampID[i]],sep='.')
  }
  samples <- unique(subset(inCts,select=c(sampID,'SAMPID')))
    
  UIDs <- data.frame(SAMPID=unique(subset(indf,select='SAMPID')),stringsAsFactors=FALSE)
  indf1 <- subset(indf,eval(as.name(trait))==1)
  if(nrow(indf1)>0){
    indf2 <- plyr::ddply(indf1,sampID,summarise,N=length(TAXON),PCTN=round(N/unique(TOTN)*100,2)
                   ,XABCOV=round(sum(XABCOV),2),XRCOV=round(sum(sXRCOV),2))

    outdf <- reshape2::melt(indf2,id.vars='SAMPID')
    outdf$variable <- paste(outdf$variable,trait,sep='_')

    
    outdf1 <- reshape2::melt(merge(UIDs,dcast(outdf,SAMPID~variable),by='SAMPID',all.x=TRUE),id.vars='UID',variable.name='PARAMETER'
                   ,value.name='RESULT')
    outdf1 <- plyr::mutate(outdf1,RESULT=ifelse(is.na(RESULT),0,RESULT)) 
  }else{
    numUIDs <- length(unique(indf$SAMPID))
    outdf1 <- data.frame(SAMPID=rep(unique(indf$SAMPID),4),PARAMETER=c(rep('N',numUIDs),rep('PCTN',numUIDs),rep('XABCOV',numUIDs)
                                                                ,rep('XRCOV',numUIDs)),RESULT=0,stringsAsFactors=F)
      outdf1$PARAMETER <- paste(outdf1$PARAMETER,trait,sep='_')
  }
  outdf2 <- merge(samples,outdf1,by='SAMPID') %>% dplyr::select(-SAMPID)
  return(outdf2)
}

#' @export 
#' @title Combine trait metric calculations
#' @description This internal function calls calcTraits_Indicator()
#' repeatedly a set of traits provided in a character vector. Not
#' intended for use on its own.
#' @param indf Input data frame containing UID, TAXON, TOTN, XABCOV,
#' sXRCOV, and variables with all names matching those in traits argument.
#' @param traits Character vector containing one or more traits variable
#' names.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @return Data frame containing UID, PARAMETER, RESULT, where values of
#' PARAMETER consist of the metric name concatenated with each trait name
#' (represented as TRAITNM below):
#'
#' N_TRAITNM: Number of taxa with trait
#'
#' PCTN_TRAITNM: Number of taxa with trait as percentage of TOTN
#'
#' XABCOV_TRAITNM: Sum of XABCOV values across taxa with trait
#'
#' XRCOV_TRAITNM: Sum of sXRCOV values across taxa with trait
#' @author Karen Blocksom
.combTraits <- function(indf,traits,sampID='UID'){
  for(i in 1:length(traits)){
    tmpOut <- .calcTraits_Indicator(indf,traits[i],sampID)
    if(i==1){
      outdf <- tmpOut
    }else{
      outdf <- rbind(outdf,tmpOut)
    }
  }
  return(outdf)
}

#' @export 
#' @title Calculate metrics using traits with only two values (0/1)
#' @description This internal function calculates metrics for traits
#' that are indicator values (0 or 1), specifically for alien and
#' cryptogenic species groups. Not intended to be used on its own.
#' @param indf Input data frame with UID, TAXON, TOTN, XABCOV, sXRCOV,
#' sRFREQ, and variable with name in argument trait. TOTN is the total
#' number of taxa at the lowest level for a UID. XABCOV is the mean
#' absolute cover of taxon. sXRCOV is the percentage of the total
#' sum of absolute cover across all taxa for a UID. sRFREQ is the relative
#' frequency of a taxon, calculated as the percentage of the total
#' frequency of taxon occurrence across all taxa for a UID.
#' @param trait Character string containing name of binary trait variable.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @return Data frame containing UID, PARAMETER, RESULT, where values of
#' PARAMETER consist of the metric name concatenated with trait name
#' (represented as TRAITNM below):
#'
#' N_TRAITNM: Number of taxa with trait value
#'
#' PCTN_TRAITNM: Number of taxa with trait value as percentage of TOTN
#'
#' XABCOV_TRAITNM: Sum of XABCOV values across taxa with trait value
#'
#' XRCOV_TRAITNM: Sum of sXRCOV values across taxa with trait value
#'
#' RFREQ_TRAITNM: Sum of sRFREQ values across taxa with trait value
#'
#' RIMP_TRAITNM: Relative importance ((RFREQ_TRAITVAL + XRCOV_TRAITVAL)/2)
#' of taxa with trait value
#' @author Karen Blocksom
.calcTraits_Indicator.alt <- function(indf,trait,sampID='UID'){
  
  for(i in 1:length(sampID)){
    if(i==1) inCts$SAMPID <- inCts[,sampID[i]]
    else inCts$SAMPID <- paste(inCts$SAMPID,inCts[,sampID[i]],sep='.')
  }
  samples <- unique(subset(inCts,select=c(sampID,'SAMPID')))
  
  UIDs <- data.frame(SAMPID=unique(subset(indf,select='SAMPID')),stringsAsFactors=FALSE)
  indf1 <- subset(indf,eval(as.name(trait))==1)
  indf2 <- plyr::ddply(indf1,c('SAMPID'),summarise,PCTN=round(length(TAXON)/unique(TOTN)*100,2)
                 ,XABCOV=round(sum(XABCOV),2),XRCOV=round(sum(sXRCOV),2)
                 ,RFREQ=round(sum(sRFREQ),2),RIMP=round((RFREQ+XRCOV)/2,2))

  outdf <- reshape2::melt(indf2,id.vars='SAMPID')
  outdf$variable <- paste(outdf$variable,trait,sep='_')
  outdf1 <- reshape2::melt(merge(UIDs,dcast(outdf,SAMPID~variable),by='UID',all.x=TRUE),id.vars='SAMPID',variable.name='PARAMETER'
                 ,value.name='RESULT')
  outdf1 <- plyr::mutate(outdf1,RESULT=ifelse(is.na(RESULT),0,RESULT))
  
  outdf2 <- merge(samples,outdf1,by='SAMPID') %>% dplyr::select(-SAMPID)
  return(outdf2)
}

#' @export 
#' @title Calculate Mean C and FQAI indices
#' @description This internal function calculates several variations of
#' Mean C and FQAI, based on cover, frequency, and number of taxa. Also
#' calculated are Not intended for use on its own.
#' @param indf Input data frame containing:
#'
#' sampID: A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#'
#' TAXON: Taxon name
#'
#' XABCOV: Mean absolute cover of taxon
#' @param subgrp  Character string of subgroup abbreviation to add
#' as suffix to metric name. Default value is NULL.
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @return Data frame containing UID, PARAMETER, and RESULT, where values
#' for PARAMETER consist of name below concatenated with subgrp value as
#' suffix (represented as SUBGRP below):
#'
#' H_SUBGRP: Shannon-Wiener Diversity Index H' = -1*sum(pi*ln(pi)), where pi
#' is proportion of species i
#'
#' J_SUBGRP: Eveness (Pielou). J = H'/ln(S), where S is number of species
#' observed
#'
#' D_SUBGRP: Simpson Diversity Index. D = 1 - sum(pi^2), where pi is
#' proportion of species i
#' @author Karen Blocksom

.calcIndices <- function(indf,subgrp=NULL,sampID='UID'){
  
  ## Calculate mean CC and FQAI indices
  indf.1 <- plyr::ddply(indf,sampID,mutate,SUBXTOTABCOV=sum(XABCOV))
  ## Calculate diversity indices
  indf.2 <- plyr::ddply(indf.1,sampID,summarise,H=round(-1*sum((XABCOV/SUBXTOTABCOV)*log(XABCOV/SUBXTOTABCOV)),4)
                 ,J=round(H/log(length(TAXON)),4),D=round(1-sum((XABCOV/SUBXTOTABCOV)^2),4))

  ## Combine all three data frames into one
  outdf <- reshape2::melt(indf.2,id.vars=sampID,variable.name='PARAMETER',value.name='RESULT') %>%
        plyr::mutate(PARAMETER=paste(as.character(PARAMETER),subgrp,sep='_')
                     ,RESULT=ifelse(is.na(RESULT)|is.infinite(RESULT),0,RESULT))

  return(outdf)
}

#' @export 
#' @title Calculate richness metrics for native status subsets
#' @description This internal function calculates species richness
#' of sample for subset based on native status values specified.
#' Not intended for use on its own.
#' @param x Data frame containing cover data summarized by UID, TAXON,
#' NWCA_NATSTAT, and DISTINCT
#' @param y Data frame containing cover data summarized by UID, PLOT,
#' TAXON, NWCA_NATSTAT, and DISTINCT
#' @param natstat Character vector containing Values of NWCA_NATSTAT
#' variable to include in rich metrics
#' @param grpname String containing suffix to add to metric name to
#' represent this group
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @return Data frame containing UID, PARAMETER, and RESULT, with one
#' row of results per parameter and UID. The values for PARAMETER
#' consist of the metric name concatenated with grpname value
#' (represented as GRP below):
#'
#' TOTN_GRP: Number of unique taxa in sample 
#'
#' XN_GRP: Mean number of taxa per plot
#'
#' MEDN_GRP: Median number of taxa per plot
#'
#' SDN_GRP: Standard deviation of number of taxa per plot
#' @author Karen Blocksom
.calcRichNS <- function(x,y,natstat,grpname,sampID='UID') {
  xx1 <- plyr::ddply(subset(x,NWCA_NATSTAT %in% natstat),sampID,summarise,TOTN_TAXA=sum(DISTINCT))
  ## Now calculate richness by plot to obtain average richness per plot
  yy1 <- merge(subset(y,NWCA_NATSTAT %in% natstat,select=c(sampID,'PLOT','DISTINCT')),xx1,by='UID')
  yy2 <- plyr::ddply(yy1,c(sampID,'PLOT','TOTN_TAXA'),summarise,N_TAXA=sum(DISTINCT))
  yy3 <- plyr::ddply(yy2,c(sampID,'TOTN_TAXA'),summarise,XN_TAXA=round(mean(N_TAXA),2),MEDN_TAXA=median(N_TAXA)
                     ,SDN_TAXA=round(sd(N_TAXA),2))

  outdf <- reshape2::melt(yy3,id.vars=sampID,variable.name='PARAMETER',value.name='RESULT')
  outdf$PARAMETER <- gsub('TAXA',grpname,outdf$PARAMETER)
  return(outdf)
}

#' @export 
#' @title Calculate mean Bray-Curtis distance among plots
#' @description This internal function calculates the mean Bray-Curtis
#' distance among   plots based on species composition.
#' @param x Data frame containing:
#'
#' sampID: A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#'
#' PLOT: Plot number of sample
#'
#' SPECIES: Character code for taxon with no spaces
#'
#' COVER: Percentage vegetative cover for a given taxon
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' @return Data frame consisting of UID and XCDIST_SPP, defined as:
#'
#' Within AA dissimilarity based on species composition = Mean of
#' between plot Bray-Curtis Distance (Dissimilarity)
#' @author Karen Blocksom
.calcXBC <- function(x,sampID='UID'){
  for(i in 1:length(sampID)){
    if(i==1) inCts$SAMPID <- inCts[,sampID[i]]
    else inCts$SAMPID <- paste(inCts$SAMPID,inCts[,sampID[i]],sep='.')
  }
  samples <- unique(subset(inCts,select=c(sampID,'SAMPID')))
  
  
  uidlist <- data.frame(SAMPID=unique(x$SAMPID),stringsAsFactors=FALSE)
  outdf <- data.frame(SAMPID=numeric(0),XBCDIST=numeric(0),stringsAsFactors=FALSE)

  for(i in 1:nrow(uidlist)){
    x1 <- subset(x,SAMPID==uidlist[i,])
    x2 <- dcast(x1,PLOT~SPECIES,value.var='COVER')
    x3 <- ecodist::distance(x2[,2:length(x2)],'bray-curtis')
    outx <- data.frame(SAMPID=uidlist[i,],XBCDIST_SPP=round(mean(x3),4),stringsAsFactors=FALSE)
    outdf <- rbind(outdf,outx)
  }
  outdf.1 <- merge(samples,outdf,by='SAMPID') %>% dplyr::select(-SAMPID)

  return(outdf.1)
}

