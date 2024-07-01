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
#' (taxon_name,'GENUS','FAMILY')
#' @param vascIn Data frame with vegetation cover data, having
#' the following variables:
#' \itemize{
#'
#' \item sampID: Variable(s) found in the argument \emph{sampID}
#'
#' \item PLOT:  Plot number of data (1 to 5 possible)
#'
#' \item Variable named in \emph{state}: Two-letter state postal code of site, used to link
#' native status to taxa in native status taxalist (inNat)
#'
#' \item Variable named in \emph{coeReg}: U.S. Army Corps of Engineers region abbreviation for
#'   sample, to correspond to GEOG_ID in Wetland Indicator Status taxalist (inWIS)
#'
#' \item Variable named in \emph{cValReg}: NWCA C-value region abbreviation for sample, to correspond
#' those in C-value taxalist
#'
#' \item Variable named in \emph{taxon_name}: Taxon name, must match with taxa data frame
#'
#' \item COVER: Percentage estimated vegetation cover by TAXON in PLOT
#' }
#' @param taxon_name String containing the name of variable for taxon name in
#' \emph{vascIn} and in taxalists.
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
#' outEx <- createDFs(
#'   sampID = "UID", "GENUS", VascPlantEx, taxon_name = "USDA_NAME",
#'   taxaNWCA, cValReg = "STATE"
#' )
#' head(outEx$byUID)
#' head(outEx$byPlot)
createDFs <- function(sampID = "UID", tvar, vascIn, taxon_name, taxa, state = "STATE",
                      coeReg = "USAC_REGION", cValReg = "NWC_CREG") {
  # Drop SPECIES_NAME_ID if present here
  if ("SPECIES_NAME_ID" %in% names(taxa)) {
    taxa <- subset(taxa, select = -SPECIES_NAME_ID)
  }

  # First merge the taxa list with the cover data by USDA_NAME
  vascIn.1 <- merge(vascIn, taxa, by = taxon_name)
  vascIn.1$tobj <- vascIn.1[, tvar] # Set tobj as the value of tvar
  vascIn.1$COVER <- with(vascIn.1, as.numeric(COVER)) # Make sure COVER is numeric

  # Set value for TAXON as either specified taxon level above species, or as USDA_NAME
  byPlot <- vascIn.1
  byPlot$TAXON <- with(byPlot, ifelse(!is.na(tobj) & tobj != "", tobj, eval(as.name(taxon_name))))

  # Sum COVER by SAMPID, PLOT, and TAXON to ensure there is only one row per species in
  # input data, set DISTINCT as 1 to be taxon counter
  if (cValReg == state) { # If cValReg variable is same as state variable, need to only specify once
    # Sum cover by taxon
    byPlot1 <- aggregate(
      x = list(COVER = byPlot$COVER),
      by = byPlot[, c(sampID, state, coeReg, "PLOT", "TAXON")],
      FUN = sum
    )
    # Set DISTINCT to 1 if cover>0
    byPlot1$DISTINCT <- with(byPlot1, ifelse(COVER > 0, 1, 0))
  } else { # If cValReg and state are not same variable, keep both
    # Sum cover by taxon
    byPlot1 <- aggregate(
      x = list(COVER = byPlot$COVER),
      by = byPlot[, c(sampID, state, coeReg, cValReg, "PLOT", "TAXON")],
      FUN = sum
    )
    # Sum cover by taxon
    byPlot1$DISTINCT <- with(byPlot1, ifelse(COVER > 0, 1, 0))
  }

  byPlot1$COVER[byPlot1$COVER > 100] <- 100 # Cap sums at 100 percent
  print("Done with plot summing")

  ## Calculate frequency and relative frequency variables by taxon
  if (cValReg == state) {
    # First calculate number of plots where taxon occurs
    byUID.sum <- aggregate(
      x = list(NUM = byPlot1$DISTINCT),
      by = byPlot1[c(sampID, state, coeReg, "TAXON")],
      FUN = sum
    )
    # Calculate mean absolute cover by taxon
    byUID.mean <- aggregate(
      x = list(XABCOV = byPlot1$COVER),
      byPlot1[, c(sampID, state, coeReg, "TAXON")],
      FUN = mean
    )
    # Count number of plots sampled (assumes zeros are filled in where taxon missing)
    byUID.length <- aggregate(
      x = list(NPLOTS = byPlot1$PLOT),
      byPlot1[, c(sampID, state, coeReg, "TAXON")],
      FUN = function(x) {
        length(unique(x))
      }
    )
    # Merge data frames
    byUID1a <- merge(byUID.sum, byUID.mean, by = c(sampID, state, coeReg, "TAXON"))
    byUID1 <- merge(byUID1a, byUID.length, by = c(sampID, state, coeReg, "TAXON"))
  } else {
    # Same as above but keeping cValReg and state variables
    byUID.sum <- aggregate(
      x = list(NUM = byPlot1$DISTINCT),
      by = byPlot1[, c(sampID, state, coeReg, cValReg, "TAXON")],
      FUN = sum
    )

    byUID.mean <- aggregate(
      x = list(XABCOV = byPlot1$COVER),
      by = byPlot1[, c(sampID, state, coeReg, cValReg, "TAXON")],
      FUN = mean
    )

    byUID.length <- aggregate(
      x = list(NPLOTS = byPlot1$PLOT),
      by = byPlot1[, c(sampID, state, coeReg, cValReg, "TAXON")],
      FUN = function(x) {
        length(unique(x))
      }
    )

    byUID1a <- merge(byUID.sum, byUID.mean, by = c(sampID, state, coeReg, cValReg, "TAXON"))
    byUID1 <- merge(byUID1a, byUID.length, by = c(sampID, state, coeReg, cValReg, "TAXON"))
  }
  # Set DISTINCT as 1 for all taxa
  byUID1$DISTINCT <- 1

  print("Done with sampID summing")

  # Create list of byUID and byPlot data frames
  byDF <- list(byUID = byUID1, byPlot = byPlot1)
  return(byDF)
}
