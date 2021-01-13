#' @export
#' 
#' @title Prepare data for calculating metrics - alternative version
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
#'   \item Variable named in \emph{state}: Two-letter state postal code of site, used to link
#' native status to taxa in native status taxalist (inNat)
#'
#'   \item Variable named in \emph{coeReg}: U.S. Army Corps of Engineers region abbreviation for
#'   sample, to correspond to GEOG_ID in Wetland Indicator Status taxalist (inWIS)
#'   
#'   \item Variable named in \emph{cValReg}: NWCA C-value regions: values must match GEOG_ID
#'   in C-value taxalist (inCVal)
#'   }   
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
#'  \item NWCA_CVAL: Coefficient of conservatism (C-value), as used in
#'  NWCA, necessary to calculate metrics based on C-values.
#'  }
#' @param inWIS Data frame with Wetland Indicator Status, from U.S. Army Corps
#'   of Engineers (USAC): \itemize{ \item USDA_NAME: Taxon name
#'   
#'   \item GEOG_ID: USAC region, abbreviated to match those used in input data
#'   frame
#'   
#'   \item WIS: Wetland Indicator Status as provided by USAC or added for NWCA }
#' 
#' @param state String containing the name of the state in \emph{vascIn}, 
#' with default value of 'STATE'  
#' 
#' @param coeReg String containing the name of the U.S. Army Corps of Engineers
#' region in \emph{vascIn} associated with Wetland Indicator Status, 
#' with default value of 'USAC_REGION'  
#' 
#' @param cValReg String containing the name of the Coefficient of 
#' Conservatism region in \emph{vascIn}, with default value of 'STATE'
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
#'   \item \emph{cValReg} (as specified in function arguments): Coefficient of 
#'   Conservatism region associated with site, necessary for merging with 
#'   inCVal, with default of 'STATE'
#'   
#'   \item \emph{state}: Two-letter state code value
#'   (will match name specified in function arguments)
#' 
#'   \item \emph{coeReg}: U.S. Army Corps of Engineers region code value 
#'   (will match name specified in function arguments)
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
#'   \item \emph{state}: Two-letter state code value associated with native status
#'   (will match name specified in function arguments)
#' 
#'   \item \emph{coeReg}: U.S. Army Corps of Engineers region code value 
#'   (will match name specified in function arguments)
#'   
#'   \item \emph{cValReg}: Site region associated with Coefficient of Conservatism 
#'   value (will match name specified in function arguments)
#'   
#'   \item TAXON: Taxon name
#'   
#'   \item COVER: Sum of cover by TAXON within plot
#'   
#'   \item DISTINCT: Distinctness of taxon, value of 1 assigned to each row }
#'   
#'   \item byUIDgen: Data frame with data summarized by sampID variables and
#'   TAXON at the genus level and contains \emph{sampID}, \emph{cValReg} variable, 
#'  \emph{coeReg},  
#'   TAXON, NUM, XABCOV, and DISTINCT. NUM is the number of plots in which taxon
#'   occurs, and XABCOV is the mean absolute COVER across plots. DISTINCT is the
#'   value 1 assigned to each row.
#'   
#'   \item byPlotgen: Data frame with data summarized by sampID variables, PLOT,
#'   and TAXON at the genus level. Each data frame contains \emph{sampID}, PLOT,
#'   \emph{cValReg} variable, \emph{coeReg}, TAXON, COVER, and DISTINCT. DISTINCT 
#'   assigns the value for each taxon as 1 if the taxon has COVER>0 or 0 if not. 
#'   COVER is the sum of the COVER variable.
#'   
#'   \item byUIDfam: Data frame with data summarized by sampID variables and
#'   TAXON at the family level and contains \emph{sampID}, \emph{cValReg} 
#'   variable, \emph{coeReg}, TAXON, NUM, XABCOV, and DISTINCT. NUM is the number of 
#'   plots in which taxon occurs, and XABCOV is the mean absolute COVER across plots. 
#'   DISTINCT is the value 1 assigned to each row.
#'   
#'   \item byPlotfam: Data frame with data summarized by \emph{sampID} , PLOT,
#'   and TAXON at the family level. Each data frame contains \emph{sampID},
#'   PLOT, \emph{cValReg} variable, \emph{coeReg}, TAXON, COVER, and DISTINCT. 
#'   DISTINCT assigns the value for each taxon as 1 if the taxon has COVER>0 or 
#'   0 if not. COVER is the sum of the COVER variable. }
#'   
#' @references US Environmental Protection Agency. 2016. National Wetland
#'   Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#'   Environmental Protection Agency, Washington, DC.
#'   
#' @author Karen Blocksom \email{blocksom.karen@epa.gov}
#'   
#' @examples
#' head(VascPlantEx)
#' prepEx <- nwcaVegData(vascIn=VascPlantEx, cValReg='STATE')
#' 
#' str(prepEx)
nwcaVegData <- function(vascIn, sampID='UID', inTaxa=taxaNWCA, inNat=ccNatNWCA, 
                        inCVal=ccNatNWCA, inWIS=wisNWCA, state = 'STATE',
                        coeReg = 'USAC_REGION', cValReg = "STATE"){
  # Read in various input datasets, and create output dataset based on available
  # types of data - must have cover data and taxonomic data at the very least
  # If
  datNames <- c(sampID, 'PLOT', 'USDA_NAME', 'COVER', state, cValReg, coeReg) 
  if(any(unique(datNames) %nin% names(vascIn))){
    print(paste("Missing key variables! Should be ", unique(datNames), ".", sep=''))
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
    print(paste("Will not be able to calculate metrics that use ",paste(msgNames, collapse=','),
                " without these parameters in inTaxa",sep=''))
  }
  
  # Coefficients of conservatism and native status values
  ccNatNames <- c('USDA_NAME','GEOG_ID')
  ccVars <- c('NWCA_CC')
  natVars <- c('NWCA_NATSTAT')
  # This only applies if someone specifies a taxalist not included in the package
  if(any(ccNatNames %nin% names(inNat))|any(ccNatNames %nin% names(inCVal))){
    print("Missing key variables! Need variables named USDA_NAME, GEOG_ID to match up
          with cover data. This taxa list cannot be used in calculations. Either revise file
          or use default taxa list by not specifying the inNat argument.")
    return(NULL)
  }
  if(any(ccVars %nin% names(inCVal))){
    msgNames <- ccVars[ccVars %nin% names(inCVal)]
    print(paste("Warning: Will not be able to calculate metrics using ", 
                paste(msgNames, collapse=','),
                " without these parameter in inCVal."))
  }
  if(any(natVars %nin% names(inNat))){
    msgNames <- natVars[natVars %nin% names(inNat)]
    print(paste("Warning: Will not be able to calculate metrics using ", 
                paste(msgNames, collapse=','),
                " without these parameter in inCVal."))
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
  dfSPP <- nwcaVegInput(sampID, 'USDA_NAME', vascIn, inTaxa, state, coeReg, cValReg)
  # For species-level data, run additional checks and add additional information
  # Merge cover data with taxalist
  dfSPP.byUID.1a <- merge(dfSPP[[1]], inTaxa, by.x='TAXON', by.y='USDA_NAME')
  
  # If any taxa in the cover data do not match up with the taxalist, return
  # missing names and end function
  if(nrow(dfSPP.byUID.1a)!=nrow(dfSPP[[1]])){
    print("Not all taxa in dfIn match up with names in taxaIn!")
    check1 <- merge(dfSPP.byUID.1a, dfSPP[[1]], by=c(sampID,'TAXON'), all.y=T)
    checkout <- unique(subset(check1, is.na(SPECIES_NAME_ID), select=c('TAXON')))
    print(checkout)
    return(NULL)
  }
  
  # If all taxa match taxalist, merge now with CC/native status by C of C region
  if(!is.null(inNat)){
    dfSPP.byUID.1b <- merge(dfSPP.byUID.1a, inNat, by.x=c('TAXON', state),
                            by.y=c('USDA_NAME','GEOG_ID'))
  }else{
    dfSPP.byUID.1b <- dfSPP.byUID.1a
  }
  if(!is.null(inCVal)){
    dfSPP.byUID.1b <- merge(dfSPP.byUID.1a, inCVal, by.x=c('TAXON', cValReg),
                            by.y=c('USDA_NAME','GEOG_ID'))
  }else{
    dfSPP.byUID.1b <- dfSPP.byUID.1a
  }
  # Merge with WIS values by USAC_REGION
  if(!is.null(inWIS)){
    dfSPP.byUID.1c <- merge(dfSPP.byUID.1b, inWIS, by.x=c('TAXON', coeReg), 
                            by.y=c('USDA_NAME','GEOG_ID'),all.x=T)
  }else{
    dfSPP.byUID.1c <- dfSPP.byUID.1b
  }
  
  # Calculate totals and add them to output data frame
  dfSPP.byUID.1c$TOTFREQ <- with(dfSPP.byUID.1c, NUM/NPLOTS)
  
  dfSPP.byUID.sum <- aggregate(x = list(XTOTABCOV = dfSPP.byUID.1c$XABCOV, 
                                        TOTFREQ = dfSPP.byUID.1c$TOTFREQ),
                               by = dfSPP.byUID.1c[c(sampID)], FUN = sum)
  
  dfSPP.byUID.length <- aggregate(x = list(TOTN = dfSPP.byUID.1c$TAXON),
                                  by = dfSPP.byUID.1c[c(sampID)], FUN = length)
  
  dfSPP.byUID.1c$TOTFREQ <- NULL
  
  dfSPP.byUID.fin <- merge(dfSPP.byUID.1c, dfSPP.byUID.sum, by = sampID)
  dfSPP.byUID.fin <- merge(dfSPP.byUID.fin, dfSPP.byUID.length, by = sampID)
  
  dfSPP.byUID.fin$TOTFREQ <- with(dfSPP.byUID.fin, TOTFREQ*100)
  dfSPP.byUID.fin$sXRCOV <- with(dfSPP.byUID.fin, XABCOV/XTOTABCOV*100)
  dfSPP.byUID.fin$FREQ <- with(dfSPP.byUID.fin, NUM/NPLOTS*100)
  dfSPP.byUID.fin$sRFREQ <- with(dfSPP.byUID.fin, (FREQ/TOTFREQ)*100)
  
  dfSPP[[1]] <- dfSPP.byUID.fin
  
  ## Also want to add NWCA_NATSTAT to dfSPP[[2]], byPlot
  dfSPP.byPlot <- merge(dfSPP[[2]], inCVal, by.x=c(cValReg, 'TAXON'), by.y=c('GEOG_ID', 'USDA_NAME'))
  dfSPP[[2]] <- dfSPP.byPlot
  
  # Create datasets for genus and family levels which will only be used for richness metrics
  dfGEN <- nwcaVegInput(sampID, 'GENUS', vascIn, inTaxa, cValReg)
  dfFAM <- nwcaVegInput(sampID, 'FAMILY', vascIn, inTaxa, cValReg)
  
  outDF <- list(byUIDspp=dfSPP[[1]], byPlotspp=dfSPP[[2]], byUIDgen=dfGEN[[1]],
                byPlotgen=dfGEN[[2]], byUIDfam=dfFAM[[1]], byPlotfam=dfFAM[[2]])
  print("Done preparing datasets")
  
  # class(outDF) <- append(class(outDF),"nwcaVegData")
  return(outDF)
  }





#' @export
#' 
#' @title Create data frames for function input - alternative version
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
#' \item \emph{state}: Two-letter state code value associated with native status
#'   (will match name specified in function arguments)
#' 
#' \item \emph{coeReg}: U.S. Army Corps of Engineers region code value 
#'   (will match name specified in function arguments)
#' 
#' \item \emph{cValReg} (as specified in function arguments): Site region associated with Coefficient of 
#' Conservatism value (will match name specified in function arguments)
#'
#' \item USDA_NAME: Taxon name, must match with taxa data frame
#'
#' \item COVER: Percentage estimated vegetation cover by TAXON in PLOT
#' }
#' @param taxa Data frame containing USDA_NAME, GENUS, CATEGORY,
#' GROWTH_HABIT, and DURATION variables. Dataset taxaNWCA is the
#' default if not specified. These variables are assumed to be
#' populated with abbreviations as found in USDA PLANTS database.
#' 
#' @param state String containing the name of the state in \emph{vascIn}, 
#' with default value of 'STATE'  
#' 
#' @param coeReg String containing the name of the U.S. Army Corps of Engineers
#' region in \emph{vascIn} associated with Wetland Indicator Status, 
#' with default value of 'USAC_REGION'  
#' 
#' @param cValReg Character string with name of variable containing 
#' region associated with Coefficients of Conservatism.
#' 
#' @return A list containing two data frames. The first data frame summarizes
#'   data by \emph{sampID} variables and TAXON and contains sampID variables,
#'   \emph{cValReg} variable, \emph{coeReg}, TAXON, NUM, XABCOV, DISTINCT, and NPLOTS. NUM is the
#'   number of plots in which taxon occurs, and XABCOV is the mean absolute
#'   COVER across plots. DISTINCT is the value 1 assigned to each row. NPLOTS is
#'   the number of plots sampled.
#'   
#'   The second summarizes by \emph{sampID} variables, PLOT, and TAXON. Each
#'   data frame contains sampID variables, PLOT, \emph{cValReg} variable, \emph{coeReg}, TAXON,
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
#' outEx <- nwcaVegInput(sampID='UID','GENUS',VascPlantEx,taxaNWCA,cValReg='STATE')
#' head(outEx$byUID)
#' head(outEx$byPlot)
nwcaVegInput <- function(sampID='UID', tvar, vascIn, taxa, state = 'STATE', 
                         coeReg = 'USAC_REGION', cValReg='STATE'){
  
    # First merge the taxa list with the cover data by USDA_NAME
  vascIn.1 <- merge(vascIn, taxa, by='USDA_NAME')
  vascIn.1$tobj <- vascIn.1[,tvar] # Set tobj as the value of tvar
  vascIn.1$COVER <- with(vascIn.1, as.numeric(COVER)) # Make sure COVER is numeric

  # Set value for TAXON as either specified taxon level above species, or as USDA_NAME
  byPlot <- vascIn.1
  byPlot$TAXON <- with(vascIn.1, ifelse(!is.na(tobj) & tobj!='',tobj, USDA_NAME))

  # Sum COVER by SAMPID, PLOT, and TAXON to ensure there is only one row per species in input data, set DISTINCT as 1 to be
  # taxon counter
  byVars <- unique(c(sampID, cValReg, coeReg,'TAXON', state))
  
  byPlot1.sum <- aggregate(x = list(COVER = byPlot$COVER), by = byPlot[c(byVars, 'PLOT')],
                       FUN = sum)
  
  byPlot1.dist <- aggregate(x = list(DISTINCT = byPlot1.sum$COVER), by = byPlot1.sum[c(byVars, 'PLOT')],
                            FUN = function(x){ifelse(x >0, 1, 0)})
  
  byPlot1 <- merge(byPlot1.sum, byPlot1.dist, by = c(byVars, 'PLOT'))
  byPlot1 <- subset(byPlot1, select = c(byVars, 'PLOT', 'COVER', 'DISTINCT'))
  
  byPlot1$COVER[byPlot1$COVER>100] <- 100 # Cap sums at 100 percent
  print("Done with plot summing")
  
  ## Calculate frequency and relative frequency variables by taxon
  byUID1.sum <- aggregate(x = list(NUM = byPlot1$DISTINCT), by = byPlot1[c(byVars)],
                          FUN = sum) 
  
  byUID1.mean <- aggregate(x = list(XABCOV = byPlot1$COVER), by = byPlot1[c(byVars)],
                           FUN = mean)
  
  byUID1.length <- aggregate(x = list(NPLOTS = byPlot1$PLOT), by = byPlot1[c(byVars)],
                             FUN = function(x){length(unique(x))})
  
  byUID2 <- merge(byUID1.sum, byUID1.mean, by = byVars)
  byUID2 <- merge(byUID2, byUID1.length, by = byVars)
  byUID2$DISTINCT <- 1
  byUID2 <- subset(byUID2, select = c(byVars, 'NUM', 'XABCOV', 'NPLOTS', 'DISTINCT'))
  
  print("Done with sampID summing")
  
  byDF <- list(byUID=byUID2, byPlot=byPlot1)
  
  return(byDF)
}