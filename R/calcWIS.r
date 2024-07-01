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
#'       categories of 'NAT', 'ADV', 'CRYP', 'INTR', 'UND'.
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
#' with OBL=1, FACW=2, FAC=3, FACU=4, UPL=5 for WETIND metrics and as
#' OBL=5, FACW=4, FAC=3 FACU=2, and UPD=1 for WETIND2 metrics:
#' \itemize{ \item WETIND_COV_ALL, WETIND_COV_NAT: Wetland Index, weighted by
#' taxon cover, based on numeric conversion of WIS, with lower numbers
#' indicating wetter conditions
#'
#' \item WETIND_FREQ_ALL, WETIND_FREQ_NAT: Wetland Index, weight by frequency,
#' based on numeric conversion of WIS, with lower numbers
#' indicating wetter conditions }
#'
#' \itemize{ \item WETIND2_COV_ALL, WETIND2_COV_NAT: Wetland Index, weighted by
#' taxon cover, based on numeric conversion of WIS, with higher numbers
#' indicating wetter conditions
#'
#' \item WETIND2_FREQ_ALL, WETIND2_FREQ_NAT: Wetland Index, weight by frequency,
#' based on numeric conversion of WIS, with higher numbers
#' indicating wetter conditions }
#'
#' If NWCA_NATSTAT is provided in the input data frame, the following metrics
#' are calculated based on Alien + Cryptogenic species:
#' \itemize{ \item
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
#' exPlant <- prepareData(VascPlantEx,
#'   taxon_name = "USDA_NAME",
#'   inTaxa = taxaNWCA, inNat = ccNatNWCA, inCVal = ccNatNWCA,
#'   inWIS = wisNWCA, cValReg = "STATE"
#' )
#'
#' wisEx <- calcWIS(exPlant$byUIDspp)
#'
#' head(wisEx)
#' unique(wisEx$PARAMETER)
calcWIS <- function(vascIn, sampID = "UID") {
  vascIn <- as.data.frame(vascIn) # Do this in case read in as a tibble or data.table, which might cause problems
  # First check for necessary variables
  necNames <- c(sampID, "TAXON", "WIS", "FREQ", "XABCOV", "TOTN", "sXRCOV")
  msgNames <- necNames[necNames %nin% names(vascIn)]
  if (length(msgNames) > 0) {
    print(paste("Missing key variables for metric calculation: ", paste(msgNames, collapse = ","),
      ". Try prepareData() function to create necessary input variables.",
      sep = ""
    ))
    return(NULL)
  }

  if ("ECOIND" %in% names(vascIn)) { # NEED TO DEAL WITH NOL species

    vascIn$WIS[vascIn$WIS %in% c("UPL", "NL")] <- "UPL"
    vascIn$WIS[is.na(vascIn$WIS) | vascIn$WIS %in% c("TBD", "UND", "NOL")] <- NA # This places taxa with UND with missing
    vascIn$ECOIND1 <- vascIn$ECOIND
    vascIn$ECOIND2 <- with(vascIn, ifelse(is.na(ECOIND), NA, 6 - as.integer(ECOIND))) # Alternate version of ECOIND numbering so higher numbers are for wetter conditions
  } else {
    vascIn$WIS[vascIn$WIS %in% c("UPL", "NL")] <- "UPL"
    vascIn$WIS[is.na(vascIn$WIS) | vascIn$WIS %in% c("TBD", "UND", "NOL")] <- NA
    vascIn$ECOIND1 <- NA
    vascIn$ECOIND1[vascIn$WIS == "OBL"] <- 1
    vascIn$ECOIND1[vascIn$WIS == "FACW"] <- 2
    vascIn$ECOIND1[vascIn$WIS == "FAC"] <- 3
    vascIn$ECOIND1[vascIn$WIS == "FACU"] <- 4
    vascIn$ECOIND1[vascIn$WIS %in% c("UPL", "NL")] <- 5
    vascIn$ECOIND1[is.na(vascIn$WIS) | vascIn$WIS %in% c("TBD", "UND", "NOL")] <- NA
    vascIn$ECOIND2 <- with(vascIn, ifelse(is.na(ECOIND1), NA, 6 - as.integer(ECOIND1)))
  }


  # Overall metric calculations
  vascIn.1 <- subset(vascIn, !is.na(WIS))

  vascIn.1$OBL_FACW <- ifelse(vascIn.1$WIS %in% c("OBL", "FACW"), 1, 0)
  vascIn.1$OBL_FACW_FAC <- ifelse(vascIn.1$WIS %in% c("OBL", "FACW", "FAC"), 1, 0)
  vascIn.1$FAC_FACU <- ifelse(vascIn.1$WIS %in% c("FACU", "FAC"), 1, 0)

  sppWIS <- int.calcTraits_MultCat(vascIn.1, "WIS", sampID)

  multTraits <- int.combTraits(vascIn.1, c("OBL_FACW", "OBL_FACW_FAC", "FAC_FACU"), sampID)

  ## Calculate Wetland indicator status metrics and melt into long format
  vascIn.2a <- subset(vascIn.1, !is.na(ECOIND1))

  vascIn.2a[, c("ECOIND1", "ECOIND2")] <- lapply(vascIn.2a[, c("ECOIND1", "ECOIND2")], as.numeric)

  vascIn.2b <- aggregate(
    x = list(
      coveco1 = vascIn.2a$XABCOV * vascIn.2a$ECOIND1,
      freqeco1 = vascIn.2a$FREQ * vascIn.2a$ECOIND1,
      coveco2 = vascIn.2a$XABCOV * vascIn.2a$ECOIND2,
      freqeco2 = vascIn.2a$FREQ * vascIn.2a$ECOIND2,
      sumabcov = vascIn.2a$XABCOV,
      sumfreq = vascIn.2a$FREQ
    ),
    by = vascIn.2a[c(sampID)],
    FUN = sum
  )

  vascIn.2b$WETIND_COV_ALL <- with(vascIn.2b, round(coveco1 / sumabcov, 2))
  vascIn.2b$WETIND_FREQ_ALL <- with(vascIn.2b, round(freqeco1 / sumfreq, 2))
  vascIn.2b$WETIND2_COV_ALL <- with(vascIn.2b, round(coveco2 / sumabcov, 2))
  vascIn.2b$WETIND2_FREQ_ALL <- with(vascIn.2b, round(freqeco2 / sumfreq, 2))
  vascIn.2b <- subset(vascIn.2b, select = c(sampID, "WETIND_COV_ALL", "WETIND_FREQ_ALL", "WETIND2_COV_ALL", "WETIND2_FREQ_ALL"))

  vascIn.2 <- reshape(vascIn.2b,
    idvar = sampID, direction = "long",
    varying = c("WETIND_COV_ALL", "WETIND_FREQ_ALL", "WETIND2_COV_ALL", "WETIND2_FREQ_ALL"),
    timevar = "PARAMETER", v.names = "RESULT",
    times = c("WETIND_COV_ALL", "WETIND_FREQ_ALL", "WETIND2_COV_ALL", "WETIND2_FREQ_ALL")
  )

  wisOut <- rbind(sppWIS, multTraits, vascIn.2)

  empty_base <- data.frame(t(rep(NA, 36)), stringsAsFactors = F)
  names(empty_base) <- c(
    "N_FAC", "N_FACU", "N_FACW", "N_OBL", "N_UPL",
    "N_OBL_FACW", "N_OBL_FACW_FAC", "N_FAC_FACU", "PCTN_FAC",
    "PCTN_FACU", "PCTN_FACW", "PCTN_OBL", "PCTN_UPL", "PCTN_OBL_FACW", "PCTN_OBL_FACW_FAC",
    "PCTN_FAC_FACU", "XABCOV_FAC", "XABCOV_FACU",
    "XABCOV_FACW", "XABCOV_OBL", "XABCOV_UPL", "XABCOV_OBL_FACW", "XABCOV_OBL_FACW_FAC",
    "XABCOV_FAC_FACU", "XRCOV_FAC", "XRCOV_FACU", "XRCOV_FACW",
    "XRCOV_OBL", "XRCOV_UPL", "XRCOV_OBL_FACW", "XRCOV_OBL_FACW_FAC", "XRCOV_FAC_FACU",
    "WETIND_COV_ALL", "WETIND_FREQ_ALL", "WETIND2_COV_ALL",
    "WETIND2_FREQ_ALL"
  )

  empty_base.nat <- data.frame(t(rep(NA, 7)), stringsAsFactors = F)
  names(empty_base.nat) <- c(
    "WETIND_COV_NAT", "WETIND_FREQ_NAT",
    "WETIND2_COV_NAT", "WETIND2_FREQ_NAT",
    "N_OBLFACW_AC", "XABCOV_OBLFACW_AC", "XRCOV_OBLFACW_AC"
  )

  # Metrics using only subsets of data based on NATSTAT_ALT
  if ("NWCA_NATSTAT" %in% names(vascIn)) {
    vascIn.nat <- vascIn
    vascIn.nat$ALIEN <- with(vascIn.nat, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV"), 1, 0))
    vascIn.nat$NATSTAT_ALT <- with(vascIn.nat, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV"), "ALIEN", NWCA_NATSTAT))
    vascIn.nat$AC <- with(vascIn.nat, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV", "CRYP"), 1, 0))

    vascIn.nat.1 <- subset(vascIn.nat, NWCA_NATSTAT == "NAT")

    ## Calculate Wetland indicator status metrics and melt into long format
    wisOut.nat.a <- subset(vascIn.nat.1, !is.na(ECOIND1))
    wisOut.nat.a[, c("ECOIND1", "ECOIND2")] <- lapply(wisOut.nat.a[, c("ECOIND1", "ECOIND2")], as.numeric)

    wisOut.nat.b <- aggregate(
      x = list(
        coveco1 = wisOut.nat.a$XABCOV * wisOut.nat.a$ECOIND1,
        freqeco1 = wisOut.nat.a$FREQ * wisOut.nat.a$ECOIND1,
        coveco2 = wisOut.nat.a$XABCOV * wisOut.nat.a$ECOIND2,
        freqeco2 = wisOut.nat.a$FREQ * wisOut.nat.a$ECOIND2,
        sumabcov = wisOut.nat.a$XABCOV,
        sumfreq = wisOut.nat.a$FREQ
      ),
      by = wisOut.nat.a[c(sampID)],
      FUN = sum
    )

    wisOut.nat.b$WETIND_COV_NAT <- with(wisOut.nat.b, round(coveco1 / sumabcov, 2))
    wisOut.nat.b$WETIND_FREQ_NAT <- with(wisOut.nat.b, round(freqeco1 / sumfreq, 2))
    wisOut.nat.b$WETIND2_COV_NAT <- with(wisOut.nat.b, round(coveco2 / sumabcov, 2))
    wisOut.nat.b$WETIND2_FREQ_NAT <- with(wisOut.nat.b, round(freqeco2 / sumfreq, 2))
    wisOut.nat.b <- subset(wisOut.nat.b, select = c(
      sampID, "WETIND_COV_NAT", "WETIND_FREQ_NAT",
      "WETIND2_COV_NAT", "WETIND2_FREQ_NAT"
    ))

    wisOut.nat.2 <- reshape(wisOut.nat.b,
      idvar = sampID, direction = "long",
      varying = c("WETIND_COV_NAT", "WETIND_FREQ_NAT", "WETIND2_COV_NAT", "WETIND2_FREQ_NAT"),
      timevar = "PARAMETER", v.names = "RESULT",
      times = c("WETIND_COV_NAT", "WETIND_FREQ_NAT", "WETIND2_COV_NAT", "WETIND2_FREQ_NAT")
    )



    vascIn.obl <- vascIn.nat
    vascIn.obl$OBLFACW_AC <- with(vascIn.obl, ifelse(WIS %in% c("OBL", "FACW") & NATSTAT_ALT %in% c("ALIEN", "CRYP"), 1, 0))

    ofOut <- int.calcTraits_Indicator(vascIn.obl, "OBLFACW_AC", sampID)
    ofOut$PARAMETER <- as.character(ofOut$PARAMETER)

    ofOut <- subset(ofOut, PARAMETER %nin% c("PCTN_OBLFACW_AC"))

    wisOut <- rbind(wisOut, wisOut.nat.2, ofOut)

    empty_base <- cbind(empty_base, empty_base.nat)
  }

  wisOut.1a <- reshape(wisOut,
    idvar = c(sampID), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(wisOut.1a) <- gsub("RESULT\\.", "", names(wisOut.1a))

  wisOut.1a <- merge(wisOut.1a, empty_base, all = TRUE)

  wisOut.1b <- reshape(wisOut.1a,
    idvar = sampID, direction = "long",
    varying = names(wisOut.1a)[!names(wisOut.1a) %in% c(sampID)],
    timevar = "PARAMETER", v.names = "RESULT",
    times = names(wisOut.1a)[!names(wisOut.1a) %in% c(sampID)]
  )
  wisOut.1b <- subset(wisOut.1b, !is.na(eval(as.name(sampID[1]))))
  wisOut.1b$RESULT <- with(wisOut.1b, ifelse(is.na(RESULT), 0, RESULT))
  wisOut.1b$PARAMETER <- with(wisOut.1b, as.character(PARAMETER))

  wisOut.1 <- subset(wisOut.1b, PARAMETER %in% names(empty_base))

  return(wisOut.1)
}
