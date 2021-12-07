


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
int.calcRich <- function(x, y, tlevel, sampID) {
  # Calculate number of taxa by sample 
  xx1 <- aggregate(x = list(TOTN_TAXA = x$DISTINCT), by = x[c(sampID)], FUN = sum)
  # Now calculate richness by plot to obtain average richness per plot
  # Merge plot data with number of taxa
  yy1 <- merge(subset(y, select=c(sampID, 'PLOT', 'DISTINCT')), xx1, by=sampID)
  # Count number of taxa by plot
  yy2 <- aggregate(x = list(N_TAXA = yy1$DISTINCT), by = yy1[c(sampID, 'PLOT', 'TOTN_TAXA')], FUN = sum)
  # From this, calculate mean number of taxa in plot
  yy3 <- aggregate(x = list(XN_TAXA = yy2$N_TAXA), by = yy2[c(sampID, 'TOTN_TAXA')], 
                   FUN = function(z){round(mean(z),2)})
  # Median number of taxa by sample
  yy4 <- aggregate(x = list(MEDN_TAXA = yy2$N_TAXA), by = yy2[c(sampID, 'TOTN_TAXA')], FUN = median)
  # Standard deviation of taxa by sample
  yy5 <- aggregate(x = list(SDN_TAXA = yy2$N_TAXA), by = yy2[c(sampID, 'TOTN_TAXA')], 
                   FUN = function(z){round(sd(z),2)})
  # Merge data frames
  zz1 <- merge(yy3, yy4, by = c(sampID, 'TOTN_TAXA')) 
  zz2 <- merge(zz1, yy5, by = c(sampID, 'TOTN_TAXA'))
  # Melt data frame 
  outdf <- reshape(zz2, idvar = sampID, direction = 'long',
                   varying= names(zz2)[!names(zz2) %in% c(sampID)],
                   timevar = 'PARAMETER', v.names = 'RESULT',
                   times = names(zz2)[!names(zz2) %in% c(sampID)])
  # Substitute tlevel values for TAXA in PARAMETER values 
  # (i.e., TOTN_GEN instead of TOTN_TAXA for tlevel='GENUS')
  outdf$PARAMETER <- gsub('TAXA', tlevel, outdf$PARAMETER)
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
int.calcTraits_MultCat <- function(vascIn, trait, sampID){
  # Subset input data to only non-missing and non-blank values of trait
  vascIn1 <- subset(vascIn,!is.na(eval(as.name(trait))) & eval(as.name(trait))!='')
  # Calculate number of taxa with trait
  vascIn.length <- aggregate(x = list(N = vascIn1$TAXON), by = vascIn1[c(sampID, trait)],
                             FUN = length)
  # Merge back to input data frame
  vascIn1a <- merge(vascIn1, vascIn.length, by = c(sampID, trait))
  # Way to obtain total number of taxa in sample 
  vascIn1.pctn <- aggregate(x = list(uniqN = vascIn1a$TOTN), by = vascIn1a[c(sampID, trait, 'N')],
                             FUN = unique)
  # Calculate PCTN from this information
  vascIn1.pctn$PCTN <- with(vascIn1.pctn, round(N/uniqN*100, 2))
  # Drop uniqN variable
  vascIn1.pctn$uniqN <- NULL
  # Calculate mean absolute and relative cover by trait
  vascIn.sum <- aggregate(x = list(XABCOV = vascIn1$XABCOV, XRCOV = vascIn1$sXRCOV),
                          by = vascIn1[c(sampID, trait)], 
                          FUN = function(x){round(sum(x), 2)})
  # Merge metric output data frames 
  vascIn2 <- merge(vascIn1.pctn, vascIn.sum, by = c(sampID, trait))
  # Melt data frame
  outdf <- reshape(vascIn2, idvar = c(sampID, trait), direction = 'long',
                   varying= names(vascIn2)[!names(vascIn2) %in% c(sampID, trait)],
                   timevar = 'variable', v.names = 'value',
                   times = names(vascIn2)[!names(vascIn2) %in% c(sampID, trait)])
  # Update metric name by combining variable and trait value
  outdf$variable <- paste(outdf$variable, outdf[,trait], sep='_')
  # Drop trait variable
  outdf[,trait] <- NULL
  # Cast data frame wide again, dropping prefix added by reshape() from variable names
  outdf.wide <- reshape(outdf, idvar = c(sampID), direction = 'wide',
                        timevar = 'variable', v.names='value')

  names(outdf.wide) <- gsub("value\\.", "", names(outdf.wide))
  # Melt data frame and replace missing values with 0
  outdf1 <- reshape(outdf.wide, idvar = sampID, direction = 'long',
                    varying = names(outdf.wide)[!names(outdf.wide) %in% c(sampID)],
                    timevar = 'PARAMETER', v.names = 'RESULT',
                    times = names(outdf.wide)[!names(outdf.wide) %in% c(sampID)])
  
  outdf1$RESULT <- with(outdf1, ifelse(is.na(RESULT), 0, RESULT))
  outdf1$PARAMETER <- with(outdf1, as.character(PARAMETER))
  
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

int.calcTraits_MultCat.alt <- function(vascIn, trait, sampID){
  # Subset input data to keep only non-missing and non-blank traits
  vascIn1 <- subset(vascIn,!is.na(eval(as.name(trait))) & eval(as.name(trait))!='')
  # calculate number of taxa with trait
  vascIn1.ntax <- aggregate(x = list(NTAX = vascIn1$TAXON), by = vascIn1[c(sampID, trait)],
                            FUN = length)
  # Merge back to input data
  vascIn1a <- merge(vascIn1, vascIn1.ntax, by = c(sampID, trait))
  # Obtain TOTN (total number of taxa) from input
  vascIn1.pctn <- aggregate(x = list(uniqN = vascIn1a$TOTN), by = vascIn1a[c(sampID, trait,'NTAX')],
                            FUN = unique)
  # calculate percent taxa for trait from above output
  vascIn1.pctn$PCTN <- with(vascIn1.pctn, round(NTAX/uniqN*100, 2))
  # Drop unnecessary variables
  vascIn1.pctn$uniqN <- NULL
  vascIn1.pctn$NTAX <- NULL
  # Calculate mean relative and absolute cover for trait, then relative importance
  vascIn1.sum <- aggregate(x = list(XABCOV = vascIn1$XABCOV, XRCOV = vascIn1$sXRCOV,
                           RFREQ = vascIn1$sRFREQ), by = vascIn1[c(sampID, trait)],
                           FUN = function(z){round(sum(z),2)})
  vascIn1.sum$RIMP <- with(vascIn1.sum, round((RFREQ+XRCOV)/2, 2))
  # Merge metric outputs 
  vascIn2 <- merge(vascIn1.pctn, vascIn1.sum, by = c(sampID, trait))
  # Melt data frame 
  outdf <- reshape(vascIn2, idvar = c(sampID, trait), direction = 'long',
                   varying= names(vascIn2)[!names(vascIn2) %in% c(sampID, trait)],
                   timevar = 'variable', v.names = 'value',
                   times = names(vascIn2)[!names(vascIn2) %in% c(sampID, trait)])
  # Update metric name by combining variable and trait name, drop trait variable
  outdf$variable <- paste(outdf$variable, outdf[,trait], sep='_')
  outdf[, trait] <- NULL
  # Cast data frame wide and drop prefix added by reshape() from variable names
  outdf.wide <- reshape(outdf, idvar = c(sampID), direction = 'wide',
                        timevar = 'variable', v.names='value')
  
  names(outdf.wide) <- gsub("value\\.", "", names(outdf.wide))
  # Melt data frame
  outdf1 <- reshape(outdf.wide, idvar = sampID, direction = 'long',
                    varying = names(outdf.wide)[!names(outdf.wide) %in% c(sampID)],
                    timevar = 'PARAMETER', v.names = 'RESULT',
                    times = names(outdf.wide)[!names(outdf.wide) %in% c(sampID)])
  # Drop metrics for undetermined trait values
  outdf1 <- subset(outdf1, substring(PARAMETER,nchar(PARAMETER)-3,nchar(PARAMETER))!='_UND')
  # Set missing result values to 0
  outdf1$RESULT <- with(outdf1, ifelse(is.na(RESULT), 0, RESULT))
  # Add _SPP to metric names
  outdf1$PARAMETER <- with(outdf1, paste(as.character(PARAMETER), 'SPP', sep=''))
  
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
int.calcTraits_Indicator <- function(vascIn, trait, sampID){
  # Create SAMPID based on variables in sampID
  for(i in 1:length(sampID)){
    if(i==1) vascIn$SAMPID <- vascIn[,sampID[i]]
    else vascIn$SAMPID <- paste(vascIn$SAMPID,vascIn[,sampID[i]],sep='.')
  }
  # Create list of unique samples 
  samples <- unique(subset(vascIn,select=c(sampID,'SAMPID')))
  # Create list of unique SAMPID values  
  UIDs <- data.frame(SAMPID=unique(subset(vascIn, select='SAMPID')), stringsAsFactors=FALSE)
  # Subset input data frame to only keep rows where trait = 1
  vascIn1 <- subset(vascIn, eval(as.name(trait))==1)
  # If there are any rows in this data frame, perform first part
  if(nrow(vascIn1)>0){
    # calculate number of taxa with trait
    vascIn.length <- aggregate(x = list(N = vascIn1$TAXON), by = vascIn1[c('SAMPID')],
                               FUN = length)
    # Merge back to input data
    vascIn1a <- merge(vascIn1, vascIn.length, by = c('SAMPID'))
    # Obtain total number of taxa from input
    vascIn1.pctn <- aggregate(x = list(uniqN = vascIn1a$TOTN), by = vascIn1a[c('SAMPID','N')],
                              FUN = unique)
    # Calculate Percent taxa with trait
    vascIn1.pctn$PCTN <- with(vascIn1.pctn, round(N/uniqN*100, 2))
    # Drop uniqN variable
    vascIn1.pctn$uniqN <- NULL
    # calculate mean absolute and relative cover by trait
    vascIn.sum <- aggregate(x = list(XABCOV = vascIn1$XABCOV, XRCOV = vascIn1$sXRCOV),
                            by = vascIn1[c('SAMPID')], 
                            FUN = function(x){round(sum(x), 2)})
    # Merge metric data frames
    vascIn2 <- merge(vascIn1.pctn, vascIn.sum, by = c('SAMPID'))
    # Melt data
    outdf <- reshape(vascIn2, idvar = c('SAMPID'), direction = 'long',
            varying= names(vascIn2)[!names(vascIn2) %in% c('SAMPID')],
            timevar = 'variable', v.names = 'value',
            times = names(vascIn2)[!names(vascIn2) %in% c('SAMPID')])
    # Rename metric by adding trait name to variable
    outdf$variable <- paste(outdf$variable, trait, sep='_')
    # Cast wide again, remove prefix added by reshape() to variable names
    outdf.wide <- reshape(outdf, idvar = c('SAMPID'), direction = 'wide',
                          timevar = 'variable', v.names='value')
    
    names(outdf.wide) <- gsub("value\\.", "", names(outdf.wide))
    # Melt data frame and replace missing RESULT with 0
    outdf1 <- reshape(outdf.wide, idvar = 'SAMPID', direction = 'long',
                      varying = names(outdf.wide)[!names(outdf.wide) %in% c('SAMPID')],
                      timevar = 'PARAMETER', v.names = 'RESULT',
                      times = names(outdf.wide)[!names(outdf.wide) %in% c('SAMPID')])
    
    outdf1$RESULT <- with(outdf1, ifelse(is.na(RESULT), 0, RESULT))
    outdf1$PARAMETER <- with(outdf1, as.character(PARAMETER))
    
  }else{ # If not rows with trait
    # Calculate number of unique samples
    numUIDs <- length(unique(vascIn$SAMPID))
    # Create data frame with appropriate metrics and 0 values for RESULT
    outdf1 <- data.frame(SAMPID=rep(unique(vascIn$SAMPID), 4), 
                         PARAMETER=c(rep('N',numUIDs),rep('PCTN', numUIDs), rep('XABCOV',numUIDs)
                                , rep('XRCOV',numUIDs)), RESULT=0, stringsAsFactors=F)
      outdf1$PARAMETER <- paste(outdf1$PARAMETER,trait,sep='_')
  }
  # Merge samples with output data frame to get back sampID variables, drop SAMPID
  outdf2 <- merge(samples, outdf1, by='SAMPID') 
  outdf2$SAMPID <- NULL
  
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
int.combTraits <- function(vascIn, traits, sampID){
  # For each value in traits variable, run int.calcTraits_Indicator function
  # then combine all trait value outputs together
  for(i in 1:length(traits)){
    tmpOut <- int.calcTraits_Indicator(vascIn, traits[i], sampID)
    if(i==1){
      outdf <- tmpOut
    }else{
      outdf <- rbind(outdf, tmpOut)
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
int.calcTraits_Indicator.alt <- function(vascIn, trait, sampID){
  # Create SAMPID variable by combining values of variables in sampID argument
  for(i in 1:length(sampID)){
    if(i==1) vascIn$SAMPID <- vascIn[, sampID[i]]
    else vascIn$SAMPID <- paste(vascIn$SAMPID, vascIn[, sampID[i]], sep='.')
  }
  # Create list of samples by sampID variables and SAMPID
  samples <- unique(subset(vascIn, select=c(sampID,'SAMPID')))
  # Create list of unique SAMPID values
  UIDs <- data.frame(SAMPID=unique(subset(vascIn, select='SAMPID')), stringsAsFactors=FALSE)
  # Subset input data to only keep trait values of 1
  vascIn1 <- subset(vascIn, eval(as.name(trait))==1)
  # Calculate number of taxa with trait
  vascIn1.ntax <- aggregate(x = list(NTAX = vascIn1$TAXON), by = vascIn1[c('SAMPID')],
                            FUN = length)
  # Merge back to subsetted input data
  vascIn1a <- merge(vascIn1, vascIn1.ntax, by = c('SAMPID'))
  # Obtain total number of taxa from input data
  vascIn1.pctn <- aggregate(x = list(uniqN = vascIn1a$TOTN), by = vascIn1a[c('SAMPID','NTAX')],
                            FUN = unique)
  # Use that value to calculate % taxa with trait, then drop unnecessary variables
  vascIn1.pctn$PCTN <- with(vascIn1.pctn, round(NTAX/uniqN*100, 2))
  vascIn1.pctn$uniqN <- NULL
  vascIn1.pctn$NTAX <- NULL
  # Calculate mean absolute and relative cover by trait
  vascIn1.sum <- aggregate(x = list(XABCOV = vascIn1$XABCOV, XRCOV = vascIn1$sXRCOV,
                                    RFREQ = vascIn1$sRFREQ), by = vascIn1[c('SAMPID')],
                           FUN = function(z){round(sum(z),2)})
  # Calculate relative importance from above metrics
  vascIn1.sum$RIMP <- with(vascIn1.sum, round((RFREQ+XRCOV)/2, 2))
  # Merge metric results
  vascIn2 <- merge(vascIn1.pctn, vascIn1.sum, by = c('SAMPID'))
  # Melt data frame and create metric name from trait name and variable 
  outdf <- reshape(vascIn2, idvar = c('SAMPID'), direction = 'long',
                   varying= names(vascIn2)[!names(vascIn2) %in% c('SAMPID')],
                   timevar = 'variable', v.names = 'value',
                   times = names(vascIn2)[!names(vascIn2) %in% c('SAMPID')])
  
  outdf$variable <- paste(outdf$variable, trait, sep='_')
  # Cast data frame wide and drop prefix added by reshape() to variable names
  outdf.wide <- reshape(outdf, idvar = c('SAMPID'), direction = 'wide',
                        timevar = 'variable', v.names='value')
  
  names(outdf.wide) <- gsub("value\\.", "", names(outdf.wide))
  # Melt data frame and remove rows with metrics based on UND (undetermined) trait values 
  outdf1 <- reshape(outdf.wide, idvar = 'SAMPID', direction = 'long',
                    varying = names(outdf.wide)[!names(outdf.wide) %in% c('SAMPID')],
                    timevar = 'PARAMETER', v.names = 'RESULT',
                    times = names(outdf.wide)[!names(outdf.wide) %in% c('SAMPID')])
  
  outdf1 <- subset(outdf1, substring(PARAMETER, nchar(PARAMETER)-3, nchar(PARAMETER))!='_UND')
  # set RESULT to 0 where missing
  outdf1$RESULT <- with(outdf1, ifelse(is.na(RESULT), 0, RESULT))
  # Merge samples to get sampID variables back, then drop SAMPID
  outdf2 <- merge(samples, outdf1, by='SAMPID')
  outdf2$SAMPID <- NULL
  
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

int.calcIndices <- function(vascIn, subgrp=NULL, sampID){
  
  ## Calculate mean C and FQAI indices
  # First calculate total absolute cover for input data frame (by sampID variables)
  vascIn.sum <- aggregate(x = list(SUBXTOTABCOV = vascIn$XABCOV),
                        by = vascIn[c(sampID)], FUN = sum)
  # Merge input dataset with totals calculated above
  vascIn.1 <- merge(vascIn, vascIn.sum, by = c(sampID))
  # Calculate diversity indices
  # Use above calculations to do further calculations needed for indices 
  vascIn.1$xrcov <- with(vascIn.1, (XABCOV/SUBXTOTABCOV))
  vascIn.1$hcalc <- with(vascIn.1, xrcov*log(xrcov))
  vascIn.1$dcalc <- with(vascIn.1, xrcov^2)
  # Use the above calculations to further calculate subparts of H and D indices
  vascIn.1.sum <- aggregate(x = list(Hsub = vascIn.1$hcalc, Dsub = vascIn.1$dcalc),
                            by = vascIn.1[c(sampID)], 
                            FUN = sum)
  # Count number of taxa
  vascIn.1.jcalc <- aggregate(x = list(jcalc = vascIn.1$TAXON), 
                              by = vascIn.1[c(sampID)],
                              FUN = length)
  # Merge metrics together
  vascIn.2 <- merge(vascIn.1.sum, vascIn.1.jcalc, by = sampID)
  # Finalize calculations of H, J, and D indices
  vascIn.2$H <- with(vascIn.2, round(-1*Hsub, 4))
  vascIn.2$J <- with(vascIn.2, round(H/log(jcalc), 4))
  vascIn.2$D <- with(vascIn.2, round(1 - Dsub, 4))
  # Subset data to select appropriate variables
  vascIn.3 <- subset(vascIn.2, select=c(sampID, 'H', 'J', 'D'))
  # Melt data frame
  outdf <- reshape(vascIn.3, idvar = sampID, direction = 'long',
          varying = c('H','J','D'),
          timevar = 'PARAMETER', v.names = 'RESULT',
          times = c('H','J','D'))
  # Update value of PARAMETER with subgrp added to end of name
  outdf$PARAMETER <- with(outdf, paste(as.character(PARAMETER), subgrp,sep='_'))
  # Fill in missing or infinite RESULT values with 0
  outdf$RESULT <- with(outdf, ifelse(is.na(RESULT)|is.infinite(RESULT), 0, RESULT))

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
int.calcRichNS <- function(x, y, natstat, grpname, sampID) {
  # Subset first input data frame, keeping only values of NWCA_NATSTAT included in
  # natstat argument - sample-level data
  xin <- subset(x, NWCA_NATSTAT %in% natstat)
  # Sum DISTINCT to get total taxa 
  xx1 <- aggregate(x = list(TOTN_TAXA = xin$DISTINCT), by = xin[c(sampID)],
                   FUN = sum)
  
  ## Now calculate richness by plot to obtain average richness per plot
  # Again, subset 2nd input data frame to keep only NWCA_NATSTAT values in natstat
  # Plot-level data
  # Then merge with xx1 by sampID
  yy1 <- merge(subset(y,NWCA_NATSTAT %in% natstat,
                      select=c(sampID,'PLOT','DISTINCT')), xx1, by=sampID)
  # Calculate number of taxa by plot
  yy2 <- aggregate(x = list(N_TAXA = yy1$DISTINCT), 
                   by = yy1[c(sampID, 'PLOT','TOTN_TAXA')],
                   FUN = sum)
  # Calculate mean number of taxa per plot
  yy3 <- aggregate(x = list(XN_TAXA = yy2$N_TAXA), by = yy2[c(sampID, 'TOTN_TAXA')], 
                   FUN = function(z){round(mean(z),2)})
  # Calculate median number of taxa per plot
  yy4 <- aggregate(x = list(MEDN_TAXA = yy2$N_TAXA), by = yy2[c(sampID, 'TOTN_TAXA')], FUN = median)
  # Calculate standard deviation in number of taxa per plot
  yy5 <- aggregate(x = list(SDN_TAXA = yy2$N_TAXA), by = yy2[c(sampID, 'TOTN_TAXA')], 
                   FUN = function(z){round(sd(z),2)})
  # Merge datasets together
  zz1 <- merge(yy3, yy4, by = c(sampID, 'TOTN_TAXA')) 
  zz2 <- merge(zz1, yy5, by = c(sampID, 'TOTN_TAXA'))
  # Melt data and update PARAMETER name with grpname (taxonomic level- SPP, GEN, FAM)
  outdf <- reshape(zz2, idvar = sampID, direction = 'long',
                   varying= names(zz2)[!names(zz2) %in% c(sampID)],
                   timevar = 'PARAMETER', v.names = 'RESULT',
                   times = names(zz2)[!names(zz2) %in% c(sampID)])
  
  outdf$PARAMETER <- gsub('TAXA', grpname, outdf$PARAMETER)
  
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
int.calcXBC <- function(x, sampID){
  # Create SAMPID based on combining variables in sampID argumetn
  for(i in 1:length(sampID)){
    if(i==1) x$SAMPID <- x[, sampID[i]]
    else x$SAMPID <- paste(x$SAMPID, x[,sampID[i]], sep='.')
  }
  # List of unique samples in dataset
  samples <- unique(subset(x, select=c(sampID,'SAMPID')))
  # Create list of unique SAMPID values
  uidlist <- data.frame(SAMPID=unique(x$SAMPID), stringsAsFactors=FALSE)
  # Create empty output data frame
  outdf <- data.frame(SAMPID=numeric(0), XBCDIST=numeric(0), stringsAsFactors=FALSE)
  # For each sample in data set, subset to just that sample, then cast wide
  # Then, remove prefix added by reshape() to variable names
  for(i in 1:nrow(uidlist)){
    x1 <- subset(x,SAMPID==uidlist[i,], select = c('PLOT','SPECIES','COVER'))
    x2 <- reshape(x1, idvar = c('PLOT'), direction = 'wide',
                  timevar = 'SPECIES', v.names='COVER')
    
    names(x2) <- gsub("COVER\\.", "", names(x2))
    # Fill in missing values in matrix with 0
    x2[is.na(x2)] <- 0
    # Calculate Bray-Curtis distance among plots sampled 
    x3 <- ecodist::distance(x2[,2:length(x2)],'bray-curtis')
    # From this, calculate mean B-C distance
    outx <- data.frame(SAMPID=uidlist[i,], XBCDIST_SPP=round(mean(x3),4), stringsAsFactors=FALSE)
    # Combine this with output data frame, adding a row at a time
    outdf <- rbind(outdf, outx)
  }
  # Merge back with samples to get sampID variables back, then drop SAMPID
  outdf.1 <- merge(samples, outdf, by='SAMPID') 
  outdf.1$SAMPID <- NULL

  return(outdf.1)
}

