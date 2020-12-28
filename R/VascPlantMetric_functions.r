#' @export
#' 
#' @title Create data frames for function input
#' 
#' @description This internal function merges the inputs of
#' plant cover data and taxalist and summarizes by the taxonomic
#' level specified. Not intended for use on its own.
#' 
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples, 'UID' by default
#' @param tvar String with the level of taxonomy
#' ('USDA_NAME','GENUS','FAMILY')
#' @param vascIn Data frame with vegetation cover data, having
#' the following variables:
#' \itemize{
#'
#' \item sampID: Variable(s) found in the argument \emph{sampID}
#'
#' \item PLOT:  Plot number of data (1 to 5 possible)
#'
#' \item STATE: Two-letter state postal code of site, used to link
#' native status to taxa in native status taxalist
#'
#' \item USAC_REGION: U.S. Army Corps of Engineers Region abbreviation for
#'   sample, to correspond to those in Wetland Indicator Status taxalist
#' 
#' \item Variable named in \emph{cValReg}: NWCA C-value region abbreviation for sample, to correspond
#' those in C-value taxalist
#'
#' \item USDA_NAME: Taxon name, must match with taxa data frame
#'
#' \item COVER: Percentage estimated vegetation cover by TAXON in PLOT
#' }
#' @param taxa Data frame containing USDA_NAME, GENUS, CATEGORY,
#' GROWTH_HABIT, and DURATION variables. Dataset taxaNWCA is the
#' default if not specified. These variables are assumed to be
#' populated with abbreviations as found in USDA PLANTS database.
#' @param cValReg String containing the name of the variable in \emph{vascIn}
#'  which specifies the C-value region. 
#' 
#' @return A list containing two data frames. The first data frame summarizes
#'   data by \emph{sampID} variables and TAXON and contains sampID variables,
#'   STATE, USAC_REGION, TAXON, NUM, XABCOV, DISTINCT, and NPLOTS. NUM is the
#'   number of plots in which taxon occurs, and XABCOV is the mean absolute
#'   COVER across plots. DISTINCT is the value 1 assigned to each row. NPLOTS is
#'   the number of plots sampled.
#'   
#'   The second summarizes by \emph{sampID} variables, PLOT, and TAXON. Each
#'   data frame contains sampID variables, PLOT, STATE, USAC_REGION, TAXON,
#'   COVER, and DISTINCT. DISTINCT assigns the value for each taxon as 1 if the 
#'   taxon has COVER>0 or 0 if not. COVER is the sum of the COVER variable by
#'   sampID variables, PLOT, and TAXON.
#'   
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#' 
#' @examples
#' head(VascPlantEx)
#' data(taxaNWCA)
#'
#' outEx <- createDFs(sampID='UID','GENUS',VascPlantEx,taxaNWCA,cValReg='STATE')
#' head(outEx$byUID)
#' head(outEx$byPlot)
createDFs <- function(sampID='UID',tvar,vascIn,taxa, cValReg='NWC_CREG'){
  
  # Drop SPECIES_NAME_ID if present here
  if('SPECIES_NAME_ID' %in% names(taxa)){
    taxa <- subset(taxa, select=-SPECIES_NAME_ID)
  }
  
  # First merge the taxa list with the cover data by USDA_NAME
  vascIn.1 <- merge(vascIn,taxa,by='USDA_NAME')
  vascIn.1$tobj <- vascIn.1[,tvar] # Set tobj as the value of tvar
  vascIn.1$COVER <- with(vascIn.1, as.numeric(COVER)) # Make sure COVER is numeric
  
  # Set value for TAXON as either specified taxon level above species, or as USDA_NAME
  byPlot <- vascIn.1
  byPlot$TAXON <- with(byPlot, ifelse(!is.na(tobj) & tobj!='',tobj,USDA_NAME))
  
  # Sum COVER by SAMPID, PLOT, and TAXON to ensure there is only one row per species in input data, set DISTINCT as 1 to be
  # taxon counter
  if(cValReg=='STATE'){
    byPlot1 <- aggregate(x = list(COVER = byPlot$COVER), 
                         by = byPlot[,c(sampID, 'STATE', 'USAC_REGION', 'PLOT', 'TAXON')],
                         FUN = sum)
    byPlot1$DISTINCT <- with(byPlot1, ifelse(COVER>0,1,0))
    
    # byPlot1 <- plyr::ddply(byPlot,c(sampID,'STATE','USAC_REGION','PLOT','TAXON'),
    #                        summarise,COVER=sum(COVER)
    #                        ,DISTINCT=ifelse(COVER>0,1,0),.progress='tk')
  }else{
    byPlot1 <- aggregate(x = list(COVER = byPlot$COVER),
                         by = byPlot[,c(sampID, 'STATE', 'USAC_REGION', cValReg, 'PLOT', 'TAXON')],
                         FUN = sum)
    byPlot1$DISTINCT <- with(byPlot1, ifelse(COVER>0,1,0))
    
    # byPlot1 <- plyr::ddply(byPlot,c(sampID,'STATE','USAC_REGION',cValReg,'PLOT','TAXON'),
    #                        summarise,COVER=sum(COVER)
    #                        ,DISTINCT=ifelse(COVER>0,1,0),.progress='tk')
  }
  
  byPlot1$COVER[byPlot1$COVER>100] <- 100 # Cap sums at 100 percent
  print("Done with plot summing")

  ## Calculate frequency and relative frequency variables by taxon
  if(cValReg=='STATE'){
    byUID.sum <- aggregate(x = list(NUM = byPlot1$DISTINCT),
                           by = byPlot1[c(sampID,'STATE','USAC_REGION','TAXON')],
                           FUN = sum)
    
    byUID.mean <- aggregate(x = list(XABCOV = byPlot1$COVER),
                            byPlot1[,c(sampID,'STATE','USAC_REGION','TAXON')],
                            FUN = mean)
    
    byUID.length <- aggregate(x = list(NPLOTS = byPlot1$PLOT),
                              byPlot1[,c(sampID,'STATE','USAC_REGION','TAXON')],
                              FUN = function(x){length(unique(x))})
    
    byUID1a <- merge(byUID.sum, byUID.mean, by = c(sampID,'STATE','USAC_REGION','TAXON'))
    byUID1 <- merge(byUID1a, byUID.length, by = c(sampID,'STATE','USAC_REGION','TAXON'))
    
    # byUID1 <- plyr::ddply(byPlot1,c(sampID,'STATE','USAC_REGION','TAXON'),
    #                       summarise,NUM=sum(DISTINCT),XABCOV=mean(COVER)
    #                       ,NPLOTS=length(unique(PLOT)))
  }else{
    byUID.sum <- aggregate(x = list(NUM = byPlot1$DISTINCT),
                           by = byPlot1[,c(sampID,'STATE','USAC_REGION',cValReg,'TAXON')],
                           FUN = sum)
    
    byUID.mean <- aggregate(x = list(XABCOV = byPlot1$COVER),
                            by = byPlot1[,c(sampID,'STATE','USAC_REGION',cValReg,'TAXON')],
                            FUN = mean)
    
    byUID.length <- aggregate(x = list(NPLOTS = byPlot1$PLOT),
                              by = byPlot1[,c(sampID,'STATE','USAC_REGION',cValReg,'TAXON')],
                              FUN = function(x){length(unique(x))})
    
    byUID1a <- merge(byUID.sum, byUID.mean, by = c(sampID,'STATE','USAC_REGION',cValReg,'TAXON'))
    byUID1 <- merge(byUID1a, byUID.length, by = c(sampID,'STATE','USAC_REGION',cValReg,'TAXON'))
    
    # byUID1 <- plyr::ddply(byPlot1,c(sampID,'STATE','USAC_REGION',cValReg,'TAXON'),
    #                     summarise,NUM=sum(DISTINCT),XABCOV=mean(COVER)
    #                     ,NPLOTS=length(unique(PLOT)))
  }
  byUID1$DISTINCT <- 1
  
  print("Done with sampID summing")

  byDF <- list(byUID=byUID1,byPlot=byPlot1)
  return(byDF)
}

#' @export
#' 
#' @title Prepare data for calculating metrics
#'   
#' @description Assemble datasets from data frame containing vegetation cover by
#'   sampID variables and PLOT, and taxa lists containing variables used in
#'   metric calculation. Species-, genus-, and family-level datasets are
#'   summarized by plot and sample. For the species-level output data summarized
#'   by sampID variables, various traits are added to the output data set.
#'   
#' @param vascIn Data frame containing cover data summarized by
#'   \emph{sampID} variables, PLOT, and taxon value, with the following variable
#'   names: \itemize{ 
#'   \item sampID A character vector containing the name(s) of 
#'   variable(s) necessary to identify unique samples, 'UID' by default
#'   
#'   \item PLOT: Plot number of data
#'   
#'   \item USDA_NAME: Taxon name, based on USDA PLANTS names, supplemented with
#'   NWCA names, if not in USDA PLANTS
#'   
#'   \item COVER: Percent cover of taxon in plot, including zeros for plots in
#'   which a taxon does not occur. 
#'   
#'   \item STATE: Two-letter state postal code of site, used to link
#' native status to taxa in native status taxalist (inNat)
#'
#'   \item USAC_REGION: U.S. Army Corps of Engineers region abbreviation for
#'   sample, to correspond to GEOG_ID in Wetland Indicator Status taxalist (inWIS)
#'   
#'   \item Variable named in \emph{cValReg}: NWCA C-value regions: values must match GEOG_ID
#'   in C-value taxalist (inCVal)
#'   }
#'
#' @param sampID A character vector containing the name(s) of variable(s)
#'   necessary to identify unique samples, 'UID' by default
#' @param inTaxa Data frame with all taxa in vascIn, with variables: \itemize{ 
#'   \item USDA_NAME: Taxon name, based on USDA PLANTS names, supplemented with
#'   NWCA names, if necessary.
#'   
#'   \item FAMILY: Family name of taxon
#'   
#'   \item GENUS: Genus name of taxon
#'   
#'   \item CATEGORY (optional): USDA PLANTS category variable, necessary to
#'   calculate category metrics.
#'   
#'   \item DURATION (optional): USDA PLANTS duration variable, necessary to
#'   calculate duration metrics.
#'   
#'   \item GROWTH_HABIT (optional): USDA PLANTS growth habit variable,necessary
#'   to calculate growth habit metrics
#'   
#'   \item SPECIES_NAME_ID (optional): Taxonomic ID number, which will be used
#'   in Bray-Curtis distance metrics if available. }
#' @param inNat Data frame with native status:
#'   
#'   \itemize{ \item USDA_NAME: Taxon name
#'   
#'   \item GEOG_ID: Postal abbreviation for STATE of taxon
#'   
#'   \item NWCA_NATSTAT: Native status variable, as used in NWCA,
#'   necessary to calculate native status metrics. }
#' @param inCVal Data frame with coefficient of conservatism values:
#'  \itemize{ \item USDA_NAME: Taxon name
#'  
#'  \item GEOG_ID: Code indicating C region for site, as supplied with cover 
#'  data.
#'  
#'  \item NWCA_CC: Coefficient of conservatism (C-value), as used in
#'  NWCA, necessary to calculate metrics based on C-values.
#'  }
#' @param inWIS Data frame with Wetland Indicator Status, from U.S. Army Corps
#'   of Engineers (USAC): \itemize{ \item USDA_NAME: Taxon name
#'   
#'   \item GEOG_ID: USAC region, abbreviated to match USAC_REGION used in 
#'   input data frame
#'   
#'   \item WIS: Wetland Indicator Status as provided by USAC or added for NWCA }
#' @param cValReg String containing the name of the variable in \emph{vascIn}
#'  which specifies the C-value region. 
#'   
#' @details This function calls the createDFs() function, which sums cover by
#'   \emph{sampID} variables, PLOT, TAXON, with sums > 100 truncated to 100
#'   percent.
#'   
#' @return A list containing six data frames: \itemize{ \item byUIDspp: Data
#'   frame with data summarized by \emph{sampID} variables and TAXON at the
#'   species level and contains: \itemize{ \item sampID: Variable(s) identified
#'   in \emph{sampID} argument
#'   
#'   \item STATE: State abbreviation of site location, necessary for merging 
#'   with inNat
#'   
#'   \item USAC_REGION: U.S. Army Corps of Engineers region abbreviation 
#'   necessary for merging with inWIS.
#'   
#'   \item TAXON: Taxon name
#'   
#'   \item NUM: Number of occurrences of taxon across plots
#'   
#'   \item XABCOV: Mean percent absolute cover of taxon across plots
#'   
#'   \item DISTINCT: Distinctness value for each taxon, 1 if the taxon has
#'   COVER>0 and 0 if not.
#'   
#'   \item NPLOTS: Number of plots in sample (1-5)
#'   
#'   \item TOTN: Total number of taxa in sample
#'   
#'   \item XTOTABCOV: Sum of \emph{XABCOV} across all taxa in sample
#'   
#'   \item sXRCOV: taxon mean relative cover (XABCOV/XTOTABCOV)*100
#'   
#'   \item FREQ: Relative number of plots in which taxon occurs (NUM/NPLOTS)*100
#'   
#'   \item TOTFREQ: Sum of \emph{FREQ} across all taxa in sample
#'   
#'   \item SRFREQ: Relative frequency of taxon relative to total frequency
#'   (FREQ/TOTFREQ)*100 } This data frame is also merged with the input taxa
#'   data frames and contains, in addition, GENUS and FAMILY, but also contains
#'   (depending on the input taxa lists): CATEGORY, DURATION, GROWTH_HABIT,
#'   NWCA_CC, NWCA_NATSTAT, WIS.
#'   
#'   \item byPlotspp: Data frame with data summarized by \emph{sampID}
#'   variables, PLOT, and TAXON at the species level. Each data frame contains: 
#'   \itemize{ \item sampID Variables identified in \emph{sampID} argument
#'   
#'   \item PLOT: Plot number
#'   
#'   \item STATE: State of sample location
#'   
#'   \item USAC_REGION: USAC region code
#'   
#'   \item \emph{cValReg}: C-value region code
#'   
#'   \item TAXON: Taxon name
#'   
#'   \item COVER: Sum of cover by TAXON within plot
#'   
#'   \item DISTINCT: Distinctness of taxon, value of 1 assigned to each row }
#'   
#'   \item byUIDgen: Data frame with data summarized by sampID variables and
#'   TAXON at the genus level and contains \emph{sampID}, STATE, USAC_REGION,
#'   NWC_CREG,
#'   TAXON, NUM, XABCOV, and DISTINCT. NUM is the number of plots in which taxon
#'   occurs, and XABCOV is the mean absolute COVER across plots. DISTINCT is the
#'   value 1 assigned to each row.
#'   
#'   \item byPlotgen: Data frame with data summarized by sampID variables, PLOT,
#'   and TAXON at the genus level. Each data frame contains \emph{sampID}, PLOT,
#'   STATE, USAC_REGION, NWC_CREG, TAXON, COVER, and DISTINCT. DISTINCT assigns 
#'   the value
#'   for each taxon as 1 if the taxon has COVER>0 or 0 if not. COVER is the sum
#'   of the COVER variable.
#'   
#'   \item byUIDfam: Data frame with data summarized by sampID variables and
#'   TAXON at the family level and contains \emph{sampID}, STATE, USAC_REGION,
#'   NWC_CREG, 
#'   TAXON, NUM, XABCOV, and DISTINCT. NUM is the number of plots in which taxon
#'   occurs, and XABCOV is the mean absolute COVER across plots. DISTINCT is the
#'   value 1 assigned to each row.
#'   
#'   \item byPlotfam: Data frame with data summarized by \emph{sampID} , PLOT,
#'   and TAXON at the family level. Each data frame contains \emph{sampID},
#'   PLOT, STATE, USAC_REGION, NWC_CREG, TAXON, COVER, and DISTINCT. DISTINCT 
#'   assigns the
#'   value for each taxon as 1 if the taxon has COVER>0 or 0 if not. COVER is
#'   the sum of the COVER variable. }
#'   
#' @references US Environmental Protection Agency. 2016. National Wetland
#'   Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#'   Environmental Protection Agency, Washington, DC.
#'   
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#'   
#' @examples
#' head(VascPlantEx)
#' prepEx <- prepareData(VascPlantEx, cValReg='STATE')
#' 
#' str(prepEx)
prepareData <- function(vascIn,sampID='UID', inTaxa=taxaNWCA, inNat=ccNatNWCA, 
                        inCVal=ccNatNWCA, inWIS=wisNWCA, cValReg='NWC_CREG'){
  # Read in various input datasets, and create output dataset based on available
  # types of data - must have cover data and taxonomic data at the very least
  # If
  datNames <- c(sampID,'PLOT','USDA_NAME','COVER','STATE','USAC_REGION',cValReg)
  if(any(datNames %nin% names(vascIn))){
    print(paste("Missing key variables! Should be ",sampID," PLOT, USDA_NAME, COVER, STATE, USAC_REGION, and ", cValReg, ".",sep=''))
    return(NULL)
  }
  vascIn <- subset(vascIn, select=names(vascIn) %in% c(sampID, 'PLOT','USDA_NAME','COVER','STATE','USAC_REGION',cValReg))
  
  # Input taxa list with taxonomy and life history traits
  if(!is.null(inTaxa)){
    taxNames <- c('USDA_NAME','FAMILY','GENUS')
    altNames <- c('CATEGORY','GROWTH_HABIT','DURATION','DUR_ALT','GRH_ALT','HERB','TREE_COMB','SHRUB_COMB','VINE_ALL')
    if(any(taxNames %nin% names(inTaxa))){
      print("Missing key variables! Need at least USDA_NAME, FAMILY, GENUS to calculate metrics.")
      return(NULL)
    }
    if(any(altNames %nin% names(inTaxa))){
      msgNames <- altNames[altNames %nin% names(inTaxa)]
      print(paste("Will not be able to calculate metrics using ",paste(msgNames,collapse=','),
                  " without these parameters in inTaxa",sep=''))
    }
    inTaxa <- subset(inTaxa, select=names(inTaxa) %in% c('USDA_NAME', taxNames, altNames,'SPECIES_NAME_ID'))
  }
  
  # Coefficients of conservatism and native status values
  if(!is.null(inNat)){
    natNames <- c('USDA_NAME','GEOG_ID')
    natVars <- c('NWCA_NATSTAT')
  # This only applies if someone specifies a taxalist not included in the package
    if(any(natNames %nin% names(inNat))){
      print("Missing key variables! Need variables named USDA_NAME, GEOG_ID to match up
            with cover data. This taxa list cannot be used in calculations. Either revise file
            or use default taxa list by not specifying the inNat argument.")
      return(NULL)
    }
    if(any(natVars %nin% names(inNat))){
      msgNames <- natVars[natVars %nin% names(inNat)]
      print(paste("Warning: Will not be able to calculate metrics using ",paste(msgNames,collapse=','),
                  " without these parameter in inNat."))
    }
    inNat <- subset(inNat, select=names(inNat) %in% c('USDA_NAME','GEOG_ID','NWCA_NATSTAT'))
  }
 
  if(!is.null(inCVal)){
    ccNames <- c('USDA_NAME','GEOG_ID')
    ccVars <- c('NWCA_CC')
    # This only applies if someone specifies a taxalist not included in the package
    if(any(ccNames %nin% names(inCVal))){
      print("Missing key variables! Need variables named USDA_NAME, GEOG_ID to match up
            with cover data. This taxa list cannot be used in calculations. Either revise file
            or use default taxa list by not specifying the inCVal argument.")
      return(NULL)
    }
    if(any(ccVars %nin% names(inCVal))){
      msgNames <- ccVars[ccVars %nin% names(inCVal)]
      print(paste("Warning: Will not be able to calculate metrics using ",paste(msgNames,collapse=','),
                  " without these parameters in inCVal."))
    }
    inCVal <- subset(inCVal, select=names(inCVal) %in% c('USDA_NAME','GEOG_ID','NWCA_CC'))
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
    inWIS <- subset(inWIS, select=names(inWIS) %in% c('USDA_NAME','GEOG_ID','WIS','ECOIND'))
  }

  ## Create dfs for species level, genus, family, and order
  # First construct list object with by plot and by sampID summarizations
  dfSPP <- createDFs(sampID,'USDA_NAME',vascIn,inTaxa,cValReg)
  print(names(dfSPP))
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

  # If all taxa match taxalist, merge now with CC/native status by cValReg
  if(!is.null(inCVal)){
    dfSPP.byUID.1b <- merge(dfSPP.byUID.1a,inCVal,by.x=c('TAXON',cValReg),by.y=c('USDA_NAME','GEOG_ID'), all.x=TRUE)
    print(names(dfSPP.byUID.1b))
  }else{
    dfSPP.byUID.1b <- dfSPP.byUID.1a
  }
  # merge with native status by STATE
  if(!is.null(inNat)){
    dfSPP.byUID.1c <- merge(dfSPP.byUID.1b,inNat,by.x=c('TAXON','STATE'),by.y=c('USDA_NAME','GEOG_ID'), all.x=TRUE)
  }else{
    dfSPP.byUID.1c <- dfSPP.byUID.1b
  }
  # Merge with WIS values by USAC_REGION
  if(!is.null(inWIS)){
    dfSPP.byUID.1d <- merge(dfSPP.byUID.1c,inWIS,by.x=c('TAXON','USAC_REGION'),by.y=c('USDA_NAME','GEOG_ID'),all.x=T)
  }else{
    dfSPP.byUID.1d <- dfSPP.byUID.1c
  }
 print(names(dfSPP.byUID.1d))
 
  # Calculate totals and add them to output data frame
 dfSPP.byUID.length <- aggregate(x = list(TOTN = dfSPP.byUID.1d$TAXON),
                                 by = dfSPP.byUID.1d[c(sampID)],
                                 FUN = length)
 
 dfSPP.byUID.1d$TOTFREQ <- with(dfSPP.byUID.1d, NUM/NPLOTS)
 dfSPP.byUID.sum <- aggregate(x = list(TOTFREQ = dfSPP.byUID.1d$TOTFREQ, XTOTABCOV = dfSPP.byUID.1d$XABCOV),
                              by = dfSPP.byUID.1d[c(sampID)], 
                              FUN = sum)
 dfSPP.byUID.sum$TOTFREQ <- dfSPP.byUID.sum$TOTFREQ*100
 
 dfSPP.byUID.1d$TOTFREQ <- NULL
 
 dfSPP.byUID.fin <- merge(dfSPP.byUID.1d, dfSPP.byUID.length, by = sampID)
 dfSPP.byUID.fin <- merge(dfSPP.byUID.fin, dfSPP.byUID.sum, by = sampID)
 
 dfSPP.byUID.fin$sXRCOV <- with(dfSPP.byUID.fin, XABCOV/XTOTABCOV*100)
 dfSPP.byUID.fin$FREQ <- with(dfSPP.byUID.fin, NUM/NPLOTS*100)
 dfSPP.byUID.fin$sRFREQ <- with(dfSPP.byUID.fin, (FREQ/TOTFREQ)*100)
 
  # dfSPP.byUID.fin <- plyr::ddply(dfSPP.byUID.1d,sampID,mutate,TOTN=length(TAXON)
  #                                ,XTOTABCOV=sum(XABCOV),TOTFREQ=sum(NUM/NPLOTS)*100) %>%
  #   plyr::mutate(sXRCOV=XABCOV/XTOTABCOV*100, FREQ=(NUM/NPLOTS)*100
  #                ,sRFREQ=(FREQ/TOTFREQ)*100)

  dfSPP[[1]] <- dfSPP.byUID.fin

  ## Also want to add NWCA_NATSTAT to dfSPP[[2]], byPlot
  dfSPP.byPlot <- merge(dfSPP[[2]],inNat,by.x=c('STATE','TAXON'),by.y=c('GEOG_ID','USDA_NAME'))
  dfSPP[[2]] <- dfSPP.byPlot

  # Create datasets for genus and family levels which will only be used for richness metrics
  dfGEN <- createDFs(sampID,'GENUS',vascIn,inTaxa,cValReg)
  dfFAM <- createDFs(sampID,'FAMILY',vascIn,inTaxa,cValReg)

  outDF <- list(byUIDspp=dfSPP[[1]],byPlotspp=dfSPP[[2]],byUIDgen=dfGEN[[1]],byPlotgen=dfGEN[[2]],byUIDfam=dfFAM[[1]],byPlotfam=dfFAM[[2]])
  print("Done preparing datasets")
  return(outDF)
}

#' @export
#' 
#' @title Calculate vascular plant richness metrics
#'   
#' @description This internal function calculates taxa richness of sample at a
#'   specified level of taxonomy, based on output from CreateDFs(). This 
#'   function used in \code{calcRichness()} function.
#'   
#' @section Warning: This function not intended for use on its own
#'   
#' @param x Data frame containing cover data summarized by \emph{sampID} 
#'   variables, TAXON, DISTINCT at the specified taxonomic level (tlevel)
#' @param y Data frame containing cover data summarized by \emph{sampID} 
#'   variables, PLOT, TAXON, and DISTINCT at the specified taxonomic level.
#' @param tlevel Taxonomic level of input dataset, abbreviated to serve as
#'   suffix to metric names ('SPP','GEN','FAM')
#' @param sampID A character vector containing the name(s) of variable(s)
#'   necessary to identify unique samples
#'   
#' @return Data frame containing \emph{sampID} variables, PARAMETER, and RESULT,
#'   with one row of results per parameter and sample. The values for PARAMETER 
#'   consist of the metric name concatenated with taxonomic level (represented 
#'   as TAXLEVEL below):
#'   
#'   TOTN_TAXLEVEL: Number of unique taxa in sample
#'   
#'   XN_TAXLEVEL: Mean number of taxa per plot
#'   
#'   MEDN_TAXLEVEL: Median number of taxa per plot
#'   
#'   SDN_TAXLEVEL: Standard deviation of number of taxa per plot
#'   
#'   N_PLOTS: Number of plots sampled for sampID
#'   
#' @author Karen Blocksom
#'   
#'   
int.calcRich <- function(x,y,tlevel,sampID) {
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
#' 
#' @description This internal function calculates metrics using traits with more
#'   than two categories as values. Used in \code{calcDuration()}, \code{calcGrowthHabit()}, 
#'   \code{calcCategory()}, \code{calcWIS()} functions.
#'   
#' @section Warning: This function not intended for use on its own
#'   
#' @param vascIn Input data frame with \emph{sampID} variables, TAXON, TOTN, 
#'   XABCOV, sXRCOV, and variable with name in \emph{trait} argument. TOTN is 
#'   the total number of taxa at the lowest level for a \emph{sampID} variables.
#'   XABCOV is the mean absolute cover of taxon. sXRCOV is the percentage of the
#'   total sum of absolute cover across all taxa for \emph{sampID} variables.
#' @param trait Character string with name of variable containing traits of 
#'   interest
#' @param sampID A character vector containing the name(s) of variable(s) 
#'   necessary to identify unique samples
#' @return Data frame containing \emph{sampID} variables, PARAMETER, RESULT, 
#'   where values of PARAMETER consist of the metric name concatenated with each
#'   value of the trait (represented as TRAITVAL below):
#'   
#'   N_TRAITVAL: Number of taxa with trait value
#'   
#'   PCTN_TRAITVAL: Number of taxa with trait value as percentage of TOTN
#'   
#'   XABCOV_TRAITVAL: Sum of XABCOV values across taxa with trait value
#'   
#'   XRCOV_TRAITVAL: Sum of sXRCOV values across taxa with trait value
#' @author Karen Blocksom
int.calcTraits_MultCat <- function(vascIn,trait,sampID){
  vascIn1 <- subset(vascIn,!is.na(eval(as.name(trait))) & eval(as.name(trait))!='')

  vascIn2 <- plyr::ddply(vascIn1,c(sampID,trait),summarise,N=length(TAXON),PCTN=round(N/unique(TOTN)*100,2)
                 ,XABCOV=round(sum(XABCOV),2),XRCOV=round(sum(sXRCOV),2))
  outdf <- reshape2::melt(vascIn2,id.vars=c(sampID,trait))
  
  formula <- paste(paste(paste(sampID,collapse='+'),'~variable',sep=''),trait,sep='+')
  outdf1 <- reshape2::melt(reshape2::dcast(outdf,eval(formula),value.var='value'),id.vars=sampID,variable.name='PARAMETER'
                           ,value.name='RESULT')
  outdf1 <- plyr::mutate(outdf1,RESULT=ifelse(is.na(RESULT),0,RESULT),PARAMETER=as.character(PARAMETER))
  return(outdf1)
}

#' @export 
#' 
#' @title Alternate metric calculations for traits with >2 categories
#' 
#' @description This internal function calculates metrics using traits having >2
#'   categories as values. Used for native status variables. Used in 
#'   \code{calcNative()} function.
#'   
#' @section Warning: This function not intended for use on its own
#'   
#' @param vascIn Input data frame with \emph{sampID} variables, TAXON, TOTN, 
#'   XABCOV, sXRCOV, sRFREQ, and variable with name in argument trait. 
#'   \emph{TOTN} is the total number of taxa at the lowest level for 
#'   \emph{sampID} variables. \emph{XABCOV} is the mean absolute cover of taxon.
#'   \emph{sXRCOV} is the percentage of the total sum of absolute cover across 
#'   all taxa for a sample. \emph{sRFREQ} is the relative frequency of a taxon, 
#'   calculated as the percentage of the total frequency of taxon occurrence 
#'   across all taxa for a sample.
#' @param trait Character string with name of variable containing traits of 
#'   interest
#' @param sampID A character vector containing the name(s) of variable(s) 
#'   necessary to identify unique samples
#'   
#' @return Data frame containing \emph{sampID} variables, PARAMETER, RESULT, 
#'   where values of PARAMETER consist of the metric name concatenated with each
#'   value of the trait (represented as TRAITVAL below):
#'   
#'   N_TRAITVAL: Number of taxa with trait value
#'   
#'   PCTN_TRAITVAL: Number of taxa with trait value as percentage of TOTN
#'   
#'   XABCOV_TRAITVAL: Sum of XABCOV values across taxa with trait value
#'   
#'   XRCOV_TRAITVAL: Sum of sXRCOV values across taxa with trait value
#'   
#'   RFREQ_TRAITVAL: Sum of sRFREQ values across taxa with trait value
#'   
#'   RIMP_TRAITVAL: Relative importance ((RFREQ_TRAITVAL + XRCOV_TRAITVAL)/2)
#' of taxa with trait value
#' 
#' @author Karen Blocksom

int.calcTraits_MultCat.alt <- function(vascIn,trait,sampID){
  vascIn1 <- subset(vascIn,!is.na(eval(as.name(trait))) & eval(as.name(trait))!='')

  vascIn2 <- plyr::ddply(vascIn1,c(sampID,trait),summarise,PCTN=round(length(TAXON)/unique(TOTN)*100,2)
                 ,XABCOV=round(sum(XABCOV),2),XRCOV=round(sum(sXRCOV),2),RFREQ=round(sum(sRFREQ),2)
                 ,RIMP=round((RFREQ+XRCOV)/2,2))
  outdf <- reshape2::melt(vascIn2,id.vars=c(sampID,trait))
  
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
#' 
#' @title Calculate metrics using traits with only two values (0/1)
#' 
#' @description This internal function calculates metrics for traits
#' that are indicator values (0 or 1). Output feeds into 
#' \code{int.combTraits()} function. Used in 
#' \code{calcWIS()} function.
#' 
#' @section Warning: This function not intended for use on its own
#' 
#' @param vascIn Input data frame with \emph{sampID} variables, TAXON, TOTN,
#'   XABCOV, sXRCOV, and a binary variable with name in \emph{trait} argument 
#'   and possible values of 0 and 1, with 1 indicating trait present for taxon.
#' @param trait Character string with name of variable containing traits of
#'   interest
#' @param sampID A character vector containing the name(s) of variable(s)
#'   necessary to identify unique samples
#'   
#' @return Data frame containing \emph{sampID} variables, PARAMETER, RESULT,
#'   where values of PARAMETER consist of the metric name concatenated with
#'   trait name (represented as TRAITNM below):
#'
#' N_TRAITNM: Number of taxa with trait
#'
#' PCTN_TRAITNM: Number of taxa with trait as percentage of \emph{TOTN}
#'
#' XABCOV_TRAITNM: Sum of \emph{XABCOV} values across taxa with trait
#'
#' XRCOV_TRAITNM: Sum of \emph{sXRCOV} values across taxa with trait
#' @author Karen Blocksom
int.calcTraits_Indicator <- function(vascIn,trait,sampID){
  
  for(i in 1:length(sampID)){
    if(i==1) vascIn$SAMPID <- vascIn[,sampID[i]]
    else vascIn$SAMPID <- paste(vascIn$SAMPID,vascIn[,sampID[i]],sep='.')
  }
  samples <- unique(subset(vascIn,select=c(sampID,'SAMPID')))
    
  UIDs <- data.frame(SAMPID=unique(subset(vascIn,select='SAMPID')),stringsAsFactors=FALSE)
  vascIn1 <- subset(vascIn,eval(as.name(trait))==1)
  if(nrow(vascIn1)>0){
    vascIn2 <- plyr::ddply(vascIn1,c('SAMPID'),summarise,N=length(TAXON),PCTN=round(N/unique(TOTN)*100,2)
                   ,XABCOV=round(sum(XABCOV),2),XRCOV=round(sum(sXRCOV),2))

    outdf <- reshape2::melt(vascIn2,id.vars='SAMPID')
    outdf$variable <- paste(outdf$variable,trait,sep='_')

    
    outdf1 <- reshape2::melt(merge(UIDs,dcast(outdf,SAMPID~variable),by='SAMPID',all.x=TRUE),id.vars='SAMPID',variable.name='PARAMETER'
                   ,value.name='RESULT')
    outdf1 <- plyr::mutate(outdf1,RESULT=ifelse(is.na(RESULT),0,RESULT)) 
  }else{
    numUIDs <- length(unique(vascIn$SAMPID))
    outdf1 <- data.frame(SAMPID=rep(unique(vascIn$SAMPID),4),PARAMETER=c(rep('N',numUIDs),rep('PCTN',numUIDs),rep('XABCOV',numUIDs)
                                                                ,rep('XRCOV',numUIDs)),RESULT=0,stringsAsFactors=F)
      outdf1$PARAMETER <- paste(outdf1$PARAMETER,trait,sep='_')
  }
  outdf2 <- merge(samples,outdf1,by='SAMPID') %>% dplyr::select(-SAMPID)
  return(outdf2)
}

#' @export 
#' 
#' @title Combine trait metric calculations
#' 
#' @description This internal function calls calcTraits_Indicator()
#' repeatedly a set of traits provided in a character vector. Used in 
#' \code{calcDuration()}, \code{calcGrowthHabit()}, \code{calcCategory()},
#' \code{calcCC()} functions.
#' 
#' @section Warning: This function not intended for use on its own
#' 
#' @param vascIn Input data frame containing \emph{sampID} variables, TAXON, 
#'   TOTN, XABCOV, sXRCOV, and variables with all names matching those in traits
#'   argument.
#' @param traits Character vector containing one or more traits variable names.
#' @param sampID A character vector containing the name(s) of variable(s) 
#'   necessary to identify unique samples
#'   
#' @return Data frame containing \emph{sampID} variables, PARAMETER, RESULT, 
#'   where values of PARAMETER consist of the metric name concatenated with each
#'   trait name (represented as TRAITNM below):
#'   
#'   N_TRAITNM: Number of taxa with trait
#'   
#'   PCTN_TRAITNM: Number of taxa with trait as percentage of \emph{TOTN}
#'   
#'   XABCOV_TRAITNM: Sum of \emph{XABCOV} values across taxa with trait
#'   
#'   XRCOV_TRAITNM: Sum of \emph{sXRCOV} values across taxa with trait
#'   
#' @author Karen Blocksom
int.combTraits <- function(vascIn,traits,sampID){
  for(i in 1:length(traits)){
    tmpOut <- int.calcTraits_Indicator(vascIn,traits[i],sampID)
    if(i==1){
      outdf <- tmpOut
    }else{
      outdf <- rbind(outdf,tmpOut)
    }
  }
  return(outdf)
}

#' @export 
#' 
#' @title Calculate metrics using traits with only two values (0/1)
#' 
#' @description This internal function calculates metrics for traits that are 
#'   indicator values (0 or 1), specifically for alien and cryptogenic species 
#'   groups. Used by \code{calcNative()} function.
#'   
#' @section Warning: This function not intended for use on its own
#'   
#' @param vascIn Input data frame with \emph{sampID} variables, TAXON, TOTN, 
#'   XABCOV, sXRCOV, sRFREQ, and variable with name in argument trait. 
#'   \emph{TOTN} is the total number of taxa at the lowest level for each sample
#'   (combination of sampID variables). \emph{XABCOV} is the mean absolute cover
#'   of taxon. sXRCOV is the percentage of the total sum of absolute cover 
#'   across all taxa for a sample. sRFREQ is the relative frequency of a taxon, 
#'   calculated as the percentage of the total frequency of taxon occurrence 
#'   across all taxa for a sample.
#' @param trait Character string containing name of binary trait variable.
#' @param sampID A character vector containing the name(s) of variable(s) 
#'   necessary to identify unique samples
#'   
#' @return Data frame containing \emph{sampID} variables, PARAMETER, RESULT, 
#'   where values of PARAMETER consist of the metric name concatenated with 
#'   trait name (represented as TRAITNM below):
#'   
#'   N_TRAITNM: Number of taxa with trait value
#'
#' PCTN_TRAITNM: Number of taxa with trait value as percentage of \emph{TOTN}
#'
#' XABCOV_TRAITNM: Sum of \emph{XABCOV} values across taxa with trait value
#'
#' XRCOV_TRAITNM: Sum of \emph{sXRCOV} values across taxa with trait value
#'
#' RFREQ_TRAITNM: Sum of \emph{sRFREQ} values across taxa with trait value
#'
#' RIMP_TRAITNM: Relative importance ((RFREQ_TRAITVAL + XRCOV_TRAITVAL)/2)
#' of taxa with trait value
#' 
#' @author Karen Blocksom
int.calcTraits_Indicator.alt <- function(vascIn,trait,sampID){
  
  for(i in 1:length(sampID)){
    if(i==1) vascIn$SAMPID <- vascIn[,sampID[i]]
    else vascIn$SAMPID <- paste(vascIn$SAMPID,vascIn[,sampID[i]],sep='.')
  }
  samples <- unique(subset(vascIn,select=c(sampID,'SAMPID')))
  
  UIDs <- data.frame(SAMPID=unique(subset(vascIn,select='SAMPID')),stringsAsFactors=FALSE)
  vascIn1 <- subset(vascIn,eval(as.name(trait))==1)
  vascIn2 <- plyr::ddply(vascIn1,c('SAMPID'),summarise,PCTN=round(length(TAXON)/unique(TOTN)*100,2)
                 ,XABCOV=round(sum(XABCOV),2),XRCOV=round(sum(sXRCOV),2)
                 ,RFREQ=round(sum(sRFREQ),2),RIMP=round((RFREQ+XRCOV)/2,2))

  outdf <- reshape2::melt(vascIn2,id.vars='SAMPID')
  outdf$variable <- paste(outdf$variable,trait,sep='_')
  outdf1 <- reshape2::melt(merge(UIDs,dcast(outdf,SAMPID~variable),by='SAMPID',all.x=TRUE),id.vars='SAMPID',variable.name='PARAMETER'
                 ,value.name='RESULT')
  outdf1 <- plyr::mutate(outdf1,RESULT=ifelse(is.na(RESULT),0,RESULT))
  
  outdf2 <- merge(samples,outdf1,by='SAMPID') %>% dplyr::select(-SAMPID)
  return(outdf2)
}

#' @export 
#' 
#' @title Calculate Mean C and FQAI indices
#' 
#' @description This internal function calculates several variations of Mean C 
#'   and FQAI, based on cover, frequency, and number of taxa. Used in 
#'   \code{calcDiversity()} function.
#'   
#' @section Warning: This function not intended for use on its own
#'   
#' @param vascIn Input data frame containing:
#'   
#'   sampID: Variables identified in \emph{sampID} argument
#'   
#'   TAXON: Taxon name
#'   
#'   XABCOV: Mean absolute cover of taxon
#' @param subgrp  Character string of subgroup abbreviation to add as suffix to 
#'   metric name. Default value is NULL.
#' @param sampID A character vector containing the name(s) of variable(s) 
#'   necessary to identify unique samples
#'   
#' @return Data frame containing sampID variables, PARAMETER, and RESULT, where 
#'   values for PARAMETER consist of name below concatenated with subgrp value 
#'   as suffix (represented as SUBGRP below):
#'   
#'   H_SUBGRP: Shannon-Wiener Diversity Index H' = -1*sum(pi*ln(pi)), where pi 
#'   is proportion of species i
#'   
#'   J_SUBGRP: Eveness (Pielou). J = H'/ln(S), where S is number of species 
#'   observed
#'   
#' D_SUBGRP: Simpson Diversity Index. D = 1 - sum(pi^2), where pi is
#' proportion of species i
#' 
#' @author Karen Blocksom

int.calcIndices <- function(vascIn,subgrp=NULL,sampID){
  
  ## Calculate mean CC and FQAI indices
  vascIn.1 <- plyr::ddply(vascIn,sampID,mutate,SUBXTOTABCOV=sum(XABCOV))
  ## Calculate diversity indices
  vascIn.2 <- plyr::ddply(vascIn.1,sampID,summarise,H=round(-1*sum((XABCOV/SUBXTOTABCOV)*log(XABCOV/SUBXTOTABCOV)),4)
                 ,J=round(H/log(length(TAXON)),4),D=round(1-sum((XABCOV/SUBXTOTABCOV)^2),4))

  ## Combine all three data frames into one
  outdf <- reshape2::melt(vascIn.2,id.vars=sampID,variable.name='PARAMETER',value.name='RESULT') %>%
        plyr::mutate(PARAMETER=paste(as.character(PARAMETER),subgrp,sep='_')
                     ,RESULT=ifelse(is.na(RESULT)|is.infinite(RESULT),0,RESULT))

  return(outdf)
}

#' @export 
#' 
#' @title Calculate richness metrics for native status subsets
#' 
#' @description This internal function calculates species richness
#' of sample for subset based on native status values specified. 
#' Used by \code{calcRichness()}
#' 
#' @section Warning: This function not intended for use on its own
#' 
#' @param x Data frame containing cover data summarized by sampID 
#' variables, TAXON, NWCA_NATSTAT, and DISTINCT
#' @param y Data frame containing cover data summarized by sampID 
#' variables, PLOT, TAXON, NWCA_NATSTAT, and DISTINCT
#' @param natstat Character vector containing Values of NWCA_NATSTAT
#' variable to include in rich metrics
#' @param grpname String containing suffix to add to metric name to
#' represent this group
#' @param sampID A character vector containing the name(s) of
#' variable(s) necessary to identify unique samples
#' 
#' @return Data frame containing \emph{sampID} variables, PARAMETER, and RESULT,
#'   with one row of results per parameter and sample. The values for PARAMETER 
#'   consist of the metric name concatenated with grpname value (represented as 
#'   GRP below):
#'   
#' TOTN_GRP: Number of unique taxa in sample 
#'
#' XN_GRP: Mean number of taxa per plot
#'
#' MEDN_GRP: Median number of taxa per plot
#'
#' SDN_GRP: Standard deviation of number of taxa per plot
#' 
#' @author Karen Blocksom
int.calcRichNS <- function(x,y,natstat,grpname,sampID) {
  xx1 <- plyr::ddply(subset(x,NWCA_NATSTAT %in% natstat),sampID,summarise,TOTN_TAXA=sum(DISTINCT))
  ## Now calculate richness by plot to obtain average richness per plot
  yy1 <- merge(subset(y,NWCA_NATSTAT %in% natstat,select=c(sampID,'PLOT','DISTINCT')),xx1,by=sampID)
  yy2 <- plyr::ddply(yy1,c(sampID,'PLOT','TOTN_TAXA'),summarise,N_TAXA=sum(DISTINCT))
  yy3 <- plyr::ddply(yy2,c(sampID,'TOTN_TAXA'),summarise,XN_TAXA=round(mean(N_TAXA),2),MEDN_TAXA=median(N_TAXA)
                     ,SDN_TAXA=round(sd(N_TAXA),2))

  outdf <- reshape2::melt(yy3,id.vars=sampID,variable.name='PARAMETER',value.name='RESULT')
  outdf$PARAMETER <- gsub('TAXA',grpname,outdf$PARAMETER)
  return(outdf)
}

#' @export
#' 
#' @title Calculate mean Bray-Curtis distance among plots
#' 
#' @description This internal function calculates the mean Bray-Curtis distance 
#'   among plots based on species composition. Used in \code{calcBCmets()} 
#'   function.
#'   
#' @section Warning: This function not intended for use on its own
#'   
#' @param x Data frame containing:
#'   
#'   sampID: A character vector containing the name(s) of variable(s) necessary 
#'   to identify unique samples
#'   
#'   PLOT: Plot number of sample
#'   
#'   SPECIES: Character code for taxon with no spaces
#'   
#'   COVER: Percentage vegetative cover for a given taxon
#' @param sampID A character vector containing the name(s) of variable(s) 
#'   necessary to identify unique samples
#'   
#' @return Data frame consisting of \emph{sampID} variables and XCDIST_SPP, 
#'   defined as:
#'   
#'   Within AA dissimilarity based on species composition = Mean of between plot
#'   Bray-Curtis Distance (Dissimilarity)
#' @author Karen Blocksom
int.calcXBC <- function(x,sampID){
  for(i in 1:length(sampID)){
    if(i==1) x$SAMPID <- x[,sampID[i]]
    else x$SAMPID <- paste(x$SAMPID,x[,sampID[i]],sep='.')
  }
  samples <- unique(subset(x,select=c(sampID,'SAMPID')))
  
  
  uidlist <- data.frame(SAMPID=unique(x$SAMPID),stringsAsFactors=FALSE)
  outdf <- data.frame(SAMPID=numeric(0),XBCDIST=numeric(0),stringsAsFactors=FALSE)

  for(i in 1:nrow(uidlist)){
    x1 <- subset(x,SAMPID==uidlist[i,])
    x2 <- dcast(x1,PLOT~SPECIES,value.var='COVER')
      x2[is.na(x2)] <- 0
    x3 <- ecodist::distance(x2[,2:length(x2)],'bray-curtis')
    outx <- data.frame(SAMPID=uidlist[i,],XBCDIST_SPP=round(mean(x3),4),stringsAsFactors=FALSE)
    outdf <- rbind(outdf,outx)
  }
  outdf.1 <- merge(samples,outdf,by='SAMPID') %>% dplyr::select(-SAMPID)

  return(outdf.1)
}

