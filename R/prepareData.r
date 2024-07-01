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
#'   \item variable name provided in \emph{taxon_name}: Taxon name, based on
#'   USDA PLANTS names, supplemented with
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
#'
#' @param sampID A character vector containing the name(s) of variable(s)
#'   necessary to identify unique samples, 'UID' by default
#' @param taxon_name String containing the name of variable for taxon name in
#' \emph{vascIn} and in taxalists.
#' @param inTaxa Data frame with all taxa in vascIn, with variables: \itemize{
#'   \item variable name provided in \emph{taxon_name}: Taxon name, consistent
#'   with name in vascIn data frame.
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
#'
#' @param inNat Data frame with native status:\itemize{
#'   \item variable name provided in \emph{taxon_name}: Taxon name, consistent
#'   with name in vascIn data frame.
#'
#'   \item GEOG_ID: Postal abbreviation for STATE of taxon
#'
#'   \item NWCA_NATSTAT: Native status variable, as used in NWCA,
#'   necessary to calculate native status metrics. }
#' @param inCVal Data frame with coefficient of conservatism values:
#'  \itemize{
#'  \item variable name provided in \emph{taxon_name}: Taxon name, consistent
#'   with name in vascIn data frame.
#'
#'  \item GEOG_ID: Code indicating C region for site, as supplied with cover
#'  data.
#'
#'  \item NWCA_CC: Coefficient of conservatism (C-value), as used in
#'  NWCA, necessary to calculate metrics based on C-values.
#'  }
#' @param inWIS Data frame with Wetland Indicator Status, from U.S. Army Corps
#'   of Engineers (USAC): \itemize{
#'   \item variable name provided in \emph{taxon_name}: Taxon name, consistent
#'   with name in vascIn data frame.
#'
#'   \item GEOG_ID: USAC region, abbreviated to match USAC_REGION used in
#'   input data frame
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
#' \item Variable named in \emph{state}: Two-letter state postal code of site, used to link
#' native status to taxa in native status taxalist (inNat)
#'
#' \item Variable named in \emph{coeReg}: U.S. Army Corps of Engineers region abbreviation for
#'   sample, to correspond to GEOG_ID in Wetland Indicator Status taxalist (inWIS)
#'
#' \item Variable named in \emph{cValReg}: NWCA C-value regions: values must match GEOG_ID
#'   in C-value taxalist (inCVal)
#' \item TAXON: Taxon name
#'
#' \item NUM: Number of occurrences of taxon across plots
#'
#' \item XABCOV: Mean percent absolute cover of taxon across plots
#'
#' \item DISTINCT: Distinctness value for each taxon, 1 if the taxon has
#'   COVER>0 and 0 if not.
#'
#' \item NPLOTS: Number of plots in sample (1-5)
#'
#' \item TOTN: Total number of taxa in sample
#'
#' \item XTOTABCOV: Sum of \emph{XABCOV} across all taxa in sample
#'
#' \item sXRCOV: taxon mean relative cover (XABCOV/XTOTABCOV)*100
#'
#' \item FREQ: Relative number of plots in which taxon occurs (NUM/NPLOTS)*100
#'
#' \item TOTFREQ: Sum of \emph{FREQ} across all taxa in sample
#'
#' \item SRFREQ: Relative frequency of taxon relative to total frequency
#'   (FREQ/TOTFREQ)*100 } This data frame is also merged with the input taxa
#'   data frames and contains, in addition, GENUS and FAMILY, but also contains
#'   (depending on the input taxa lists): CATEGORY, DURATION, GROWTH_HABIT,
#'   NWCA_CC, NWCA_NATSTAT, WIS.
#'
#' \item byPlotspp: Data frame with data summarized by \emph{sampID}
#'   variables, PLOT, and TAXON at the species level. Each data frame contains:
#'   \itemize{ \item sampID Variables identified in \emph{sampID} argument
#'
#'   \item PLOT: Plot number
#'
#'   \item \emph{state}: State of sample location
#'
#'   \item \emph{coeReg}: USAC region code
#'
#'   \item \emph{cValReg}: C-value region code
#'
#'   \item variable name provided in \emph{taxon_name}: Taxon name, consistent
#'   with name in vascIn data frame.
#'
#'   \item COVER: Sum of cover by TAXON within plot
#'
#'   \item DISTINCT: Distinctness of taxon, value of 1 assigned to each row }
#'
#'   \item byUIDgen: Data frame with data summarized by sampID variables and
#'   TAXON at the genus level and contains \emph{sampID}, \emph{state},
#'   \emph{coeReg}, \emph{cValReg},
#'   TAXON, NUM, XABCOV, and DISTINCT. NUM is the number of plots in which taxon
#'   occurs, and XABCOV is the mean absolute COVER across plots. DISTINCT is the
#'   value 1 assigned to each row.
#'
#'   \item byPlotgen: Data frame with data summarized by sampID variables, PLOT,
#'   and TAXON at the genus level. Each data frame contains \emph{sampID}, PLOT,
#'   \emph{state}, \emph{coeReg}, \emph{cValReg}, TAXON, COVER, and DISTINCT.
#'   DISTINCT assigns the value
#'   for each taxon as 1 if the taxon has COVER>0 or 0 if not. COVER is the sum
#'   of the COVER variable.
#'
#'   \item byUIDfam: Data frame with data summarized by sampID variables and
#'   TAXON at the family level and contains \emph{sampID}, \emph{state},
#'   \emph{coeReg}, \emph{cValReg},
#'   TAXON, NUM, XABCOV, and DISTINCT. NUM is the number of plots in which taxon
#'   occurs, and XABCOV is the mean absolute COVER across plots. DISTINCT is the
#'   value 1 assigned to each row.
#'
#'   \item byPlotfam: Data frame with data summarized by \emph{sampID} , PLOT,
#'   and TAXON at the family level. Each data frame contains \emph{sampID},
#'   PLOT, \emph{state}, \emph{coeReg}, \emph{cValReg}, TAXON, COVER,
#'   and DISTINCT. DISTINCT assigns the
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
#' prepEx <- prepareData(VascPlantEx,
#'   taxon_name = "USDA_NAME", state = "STATE",
#'   coeReg = "USAC_REGION", cValReg = "STATE"
#' )
#'
#' str(prepEx)
prepareData <- function(vascIn, sampID = "UID", taxon_name, inTaxa = taxaNWCA_2016,
                        inNat = nativeNWCA_2016,
                        inCVal = cvalNWCA_2016, inWIS = wisNWCA_2016, state = "STATE",
                        coeReg = "USAC_REGION", cValReg = "NWC_CREG") {
  # Read in various input datasets, and create output dataset based on available
  # types of data - must have cover data and taxonomic data at the very least
  datNames <- c(sampID, "PLOT", taxon_name, "COVER", state, coeReg, cValReg)
  # Alert user if variables are missing and stop function
  if (any(datNames %nin% names(vascIn))) {
    print(paste("Missing key variables! Should be ", sampID, " PLOT, ", taxon_name, ", COVER,",
      state, coeReg, " and ", cValReg, ".",
      sep = ""
    ))
    return(NULL)
  }
  # Subset data to only keep relevant variables
  vascIn <- subset(vascIn, select = names(vascIn) %in% c(
    sampID, "PLOT", taxon_name, "COVER",
    state, coeReg, cValReg
  ))

  # Input taxa list with taxonomy and life history traits
  if (!is.null(inTaxa)) {
    # Need certain taxonomy levels, other traits are optional
    taxNames <- c(taxon_name, "FAMILY", "GENUS")
    altNames <- c("CATEGORY", "GROWTH_HABIT", "DURATION", "DUR_ALT", "GRH_ALT", "HERB", "TREE_COMB", "SHRUB_COMB", "VINE_ALL")
    # Alert user if any necessary names are missing from input taxalist and stop function
    if (any(taxNames %nin% names(inTaxa))) {
      print("Missing key variables! Need at least USDA_NAME, FAMILY, GENUS to calculate metrics.")
      return(NULL)
    }
    # Alert user if optional traits are missing that certain metrics may not be calculated
    if (any(altNames %nin% names(inTaxa))) {
      msgNames <- altNames[altNames %nin% names(inTaxa)]
      print(paste("Will not be able to calculate metrics using ", paste(msgNames, collapse = ","),
        " without these parameters in inTaxa",
        sep = ""
      ))
    }
    # Subset taxalist to keep only relevant variables
    inTaxa <- subset(inTaxa, select = names(inTaxa) %in% c(taxon_name, taxNames, altNames, "SPECIES_NAME_ID"))
  }

  # Check for names in native status taxalist, both necessary and optional
  if (!is.null(inNat)) {
    natNames <- c(taxon_name, "GEOG_ID", "NWCA_NATSTAT")
    # This only applies if someone specifies a taxalist not included in the package
    if (any(natNames %nin% names(inNat))) {
      print(paste("Missing key variables! Need variables named ", taxon_name, " GEOG_ID, and
            NWCA_NATSTAT to match up
            with cover data. This taxa list cannot be used in calculations. Either revise file
            or use default taxa list by not specifying the inNat argument."))
      return(NULL)
    }
    inNat <- subset(inNat, select = names(inNat) %in% c(taxon_name, "GEOG_ID", "NWCA_NATSTAT"))
  }

  # Check for names in C-values taxalist, both necessary
  if (!is.null(inCVal)) {
    if ("NWCA_CC" %nin% names(inCVal) & "NWCA_CVAL" %in% names(inCVal)) {
      inCVal$NWCA_CC <- inCVal$NWCA_CVAL
    }
    ccNames <- c(taxon_name, "GEOG_ID", "NWCA_CC")
    # This only applies if someone specifies a taxalist not included in the package
    if (any(ccNames %nin% names(inCVal))) {
      print(paste("Missing key variables! Need variables named ", taxon_name, " GEOG_ID, and
            NWCA_CC to match up
            with cover data. This taxa list cannot be used in calculations. Either revise file
            or use default taxa list by not specifying the inCVal argument."))
      return(NULL)
    }
    inCVal <- subset(inCVal, select = names(inCVal) %in% c(taxon_name, "GEOG_ID", "NWCA_CC"))
  }

  # Wetland Indicator Status values
  if (!is.null(inWIS)) {
    wisNames <- c(taxon_name, "GEOG_ID", "WIS")
    if (any(wisNames %nin% names(inWIS))) {
      print(paste("Missing key variables! Need variables named ", taxon_name, " GEOG_ID, and
            WIS to match up
            with cover data. This taxa list cannot be used in calculations. Either revise file
            or use default taxa list by not specifying the inWIS argument."))
      return(NULL)
    }
    inWIS <- subset(inWIS, select = names(inWIS) %in% c(taxon_name, "GEOG_ID", "WIS", "ECOIND"))
  }

  ## Create dfs for species level, genus, family, and order
  # First construct list object with by plot and by sampID summarizations
  dfSPP <- createDFs(sampID, taxon_name, vascIn, taxon_name, inTaxa, state, coeReg, cValReg)
  print(names(dfSPP))

  # For species-level data, run additional checks and add additional information
  # Merge cover data with taxalist
  dfSPP.byUID.1a <- merge(dfSPP[[1]], inTaxa, by.x = "TAXON", by.y = taxon_name)

  # If any taxa in the cover data do not match up with the taxalist, return
  # missing names and end function
  if (nrow(dfSPP.byUID.1a) != nrow(dfSPP[[1]])) {
    print("Not all taxa in dfIn match up with names in taxaIn!")
    check1 <- merge(dfSPP.byUID.1a, dfSPP[[1]], by = c(sampID, "TAXON"), all.y = T)
    checkout <- unique(subset(check1, is.na(SPECIES_NAME_ID), select = c("TAXON")))
    print(checkout)
    return(NULL)
  }

  # If all taxa match taxalist, merge now with C-value status by cValReg variable
  if (!is.null(inCVal)) {
    dfSPP.byUID.1b <- merge(dfSPP.byUID.1a, inCVal,
      by.x = c("TAXON", cValReg),
      by.y = c(taxon_name, "GEOG_ID"), all.x = TRUE
    )
    print(names(dfSPP.byUID.1b))
  } else {
    dfSPP.byUID.1b <- dfSPP.byUID.1a
  }
  # merge with native status by state variable
  if (!is.null(inNat)) {
    dfSPP.byUID.1c <- merge(dfSPP.byUID.1b, inNat,
      by.x = c("TAXON", state),
      by.y = c(taxon_name, "GEOG_ID"), all.x = TRUE
    )
  } else {
    dfSPP.byUID.1c <- dfSPP.byUID.1b
  }
  # Merge with WIS values by coeReg variable
  if (!is.null(inWIS)) {
    dfSPP.byUID.1d <- merge(dfSPP.byUID.1c, inWIS,
      by.x = c("TAXON", coeReg),
      by.y = c(taxon_name, "GEOG_ID"), all.x = T
    )
  } else {
    dfSPP.byUID.1d <- dfSPP.byUID.1c
  }
  print(names(dfSPP.byUID.1d))

  # Calculate totals and add them to output data frame
  # Number of taxa
  dfSPP.byUID.length <- aggregate(
    x = list(TOTN = dfSPP.byUID.1d$TAXON),
    by = dfSPP.byUID.1d[c(sampID)],
    FUN = length
  )
  # Number of plots for taxon/number of plots sampled (frequency)
  dfSPP.byUID.1d$TOTFREQ <- with(dfSPP.byUID.1d, NUM / NPLOTS)
  # Now sum freqency and absolute cover across all taxa
  dfSPP.byUID.sum <- aggregate(
    x = list(
      TOTFREQ = dfSPP.byUID.1d$TOTFREQ,
      XTOTABCOV = dfSPP.byUID.1d$XABCOV
    ),
    by = dfSPP.byUID.1d[c(sampID)],
    FUN = sum
  )
  # Multiply frequency by 100 for percent over all taxa
  dfSPP.byUID.sum$TOTFREQ <- dfSPP.byUID.sum$TOTFREQ * 100
  # Drop TOTFREQ from interim data frame
  dfSPP.byUID.1d$TOTFREQ <- NULL
  # Merge calculated data frames
  dfSPP.byUID.fin <- merge(dfSPP.byUID.1d, dfSPP.byUID.length, by = sampID)
  dfSPP.byUID.fin <- merge(dfSPP.byUID.fin, dfSPP.byUID.sum, by = sampID)
  # Now calculate relative cover and frequency variables for each taxon
  dfSPP.byUID.fin$sXRCOV <- with(dfSPP.byUID.fin, XABCOV / XTOTABCOV * 100)
  dfSPP.byUID.fin$FREQ <- with(dfSPP.byUID.fin, NUM / NPLOTS * 100)
  dfSPP.byUID.fin$sRFREQ <- with(dfSPP.byUID.fin, (FREQ / TOTFREQ) * 100)
  # This inserts data frame into first position in list of species-level data
  dfSPP[[1]] <- dfSPP.byUID.fin

  ## Also want to add NWCA_NATSTAT to dfSPP[[2]], byPlot for species-level data
  dfSPP.byPlot <- merge(dfSPP[[2]], inNat, by.x = c(state, "TAXON"), by.y = c("GEOG_ID", taxon_name))
  dfSPP[[2]] <- dfSPP.byPlot

  # Create datasets for genus and family levels which will only be used for richness metrics
  dfGEN <- createDFs(sampID, "GENUS", vascIn, taxon_name, inTaxa, state, coeReg, cValReg)
  dfFAM <- createDFs(sampID, "FAMILY", vascIn, taxon_name, inTaxa, state, coeReg, cValReg)
  # Create full output list using species-, genus-, and family-level data frames by UID and plot
  outDF <- list(
    byUIDspp = dfSPP[[1]], byPlotspp = dfSPP[[2]], byUIDgen = dfGEN[[1]],
    byPlotgen = dfGEN[[2]], byUIDfam = dfFAM[[1]], byPlotfam = dfFAM[[2]]
  )

  print("Done preparing datasets")
  return(outDF)
}
