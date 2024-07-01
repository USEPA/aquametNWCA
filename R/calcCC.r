# Metrics using C values
#' @export
#'
#' @title Calculate Coefficient of Conservatism metrics
#'
#' @description This function calculates C of C
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
#'       categories of 'NAT', 'ADV', 'CRYP', 'INTR', 'UND'.
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
#' exPlant <- prepareData(VascPlantEx,
#'   taxon_name = "USDA_NAME",
#'   inTaxa = taxaNWCA, inNat = ccNatNWCA, inCVal = ccNatNWCA,
#'   inWIS = wisNWCA, cValReg = "STATE"
#' )
#'
#' ccEx <- calcCC(exPlant$byUIDspp)
#'
#' head(ccEx)
#' unique(ccEx$PARAMETER)
calcCC <- function(vascIn, sampID = "UID") {
  vascIn <- as.data.frame(vascIn) # Do this in case read in as a tibble or data.table, which might cause problems
  if ("NWCA_CC" %nin% names(vascIn)) {
    print("Missing NWCA_CC from input data frame - cannot calculate metrics!")
    return(NULL)
  }

  ## Calculate mean CC and FQAI indices
  vascIn.sub <- subset(vascIn, toupper(NWCA_CC) != "UND")

  totals <- aggregate(
    x = list(SUBTOTFREQ = vascIn$FREQ, SUBXTOTABCOV = vascIn$XABCOV), by = vascIn[c(sampID)],
    FUN = sum
  )
  numtaxa <- aggregate(
    x = list(NUMTAXA = vascIn.sub$TAXON), by = vascIn.sub[c(sampID)],
    FUN = function(x) {
      length(unique(x))
    }
  )

  vascIn.1a <- merge(subset(vascIn.sub, select = names(vascIn.sub) %nin% c("TOTFREQ", "XTOTABCOV")),
    totals,
    by = sampID
  )
  vascIn.1 <- merge(vascIn.1a, numtaxa, by = sampID)

  vascIn.2 <- vascIn.1
  vascIn.2$NWCA_CC <- as.numeric(vascIn.2$NWCA_CC)

  vascIn.2$XC_FREQ <- with(vascIn.2, ((FREQ / SUBTOTFREQ) * NWCA_CC * 100))
  vascIn.2$FQAI_FREQ <- with(vascIn.2, ((FREQ / SUBTOTFREQ) * NWCA_CC * 100))
  vascIn.2$XC_COV <- with(vascIn.2, (XABCOV / SUBXTOTABCOV) * NWCA_CC * 100)
  vascIn.2$FQAI_COV <- with(vascIn.2, (XABCOV / SUBXTOTABCOV) * NWCA_CC * 100)

  vascIn.2a <- aggregate(
    x = list(
      XC = vascIn.2$NWCA_CC, XC_FREQ = vascIn.2$XC_FREQ,
      XC_COV = vascIn.2$XC_COV
    ),
    by = vascIn.2[c(sampID, "NUMTAXA")],
    FUN = sum
  )

  vascIn.2a$XC <- with(vascIn.2a, round(XC / NUMTAXA, 2))
  vascIn.2a$XC_FREQ <- with(vascIn.2a, round(XC_FREQ / NUMTAXA, 2))
  vascIn.2a$XC_COV <- with(vascIn.2a, round(XC_COV / NUMTAXA, 2))
  vascIn.2a$NUMTAXA <- NULL

  vascIn.2b <- aggregate(
    x = list(
      FQAI = vascIn.2$NWCA_CC,
      FQAI_FREQ = vascIn.2$FQAI_FREQ,
      FQAI_COV = vascIn.2$FQAI_COV
    ),
    by = vascIn.2[c(sampID, "NUMTAXA")], FUN = sum
  )

  vascIn.2b$FQAI <- with(vascIn.2b, round(FQAI / sqrt(NUMTAXA), 2))
  vascIn.2b$FQAI_FREQ <- with(vascIn.2b, round(FQAI_FREQ / sqrt(NUMTAXA), 2))
  vascIn.2b$FQAI_COV <- with(vascIn.2b, round(FQAI_COV / sqrt(NUMTAXA), 2))
  vascIn.2b$NUMTAXA <- NULL

  vascIn.2c <- merge(vascIn.2a, vascIn.2b, by = sampID)

  ccOut <- reshape(vascIn.2c,
    idvar = sampID, direction = "long",
    varying = names(vascIn.2c)[!names(vascIn.2c) %in% sampID],
    timevar = "variable", v.names = "RESULT",
    times = names(vascIn.2c)[!names(vascIn.2c) %in% sampID]
  )

  ccOut$PARAMETER <- paste(ccOut$variable, "ALL", sep = "_")
  ccOut$variable <- NULL

  # Create empty data frame
  empty_base <- data.frame(t(rep(NA, 26)), stringsAsFactors = F)
  names(empty_base) <- c(
    "XC_ALL", "FQAI_ALL", "XC_FREQ_ALL", "FQAI_FREQ_ALL", "XC_COV_ALL", "FQAI_COV_ALL",
    "N_HTOL", "PCTN_HTOL", "XABCOV_HTOL", "XRCOV_HTOL",
    "N_HSEN", "PCTN_HSEN", "XABCOV_HSEN", "XRCOV_HSEN",
    "N_TOL", "PCTN_TOL", "XABCOV_TOL", "XRCOV_TOL",
    "N_ISEN", "PCTN_ISEN", "XABCOV_ISEN", "XRCOV_ISEN",
    "N_SEN", "PCTN_SEN", "XABCOV_SEN", "XRCOV_SEN"
  )

  empty_base.nat <- data.frame(t(rep(NA, 6)), stringsAsFactors = F)
  names(empty_base.nat) <- c("XC_NAT", "FQAI_NAT", "XC_FREQ_NAT", "FQAI_FREQ_NAT", "XC_COV_NAT", "FQAI_COV_NAT")


  # Now create recoded indicator variables based on NWCA_CC values
  vascIn.alt <- vascIn
  vascIn.alt$SEN <- with(vascIn.alt, ifelse(NWCA_CC %in% c("7", "8", "9", "10"), 1, 0))
  vascIn.alt$ISEN <- with(vascIn.alt, ifelse(NWCA_CC %in% c("5", "6"), 1, 0))
  vascIn.alt$TOL <- with(vascIn.alt, ifelse(NWCA_CC %in% c("4", "3", "2", "1", "0"), 1, 0))
  vascIn.alt$HTOL <- with(vascIn.alt, ifelse(NWCA_CC %in% c("0", "1", "2"), 1, 0))
  vascIn.alt$HSEN <- with(vascIn.alt, ifelse(NWCA_CC %in% c("9", "10"), 1, 0))

  multTraits <- int.combTraits(vascIn.alt, c("SEN", "TOL", "ISEN", "HTOL", "HSEN"), sampID)

  ccOut <- rbind(ccOut, multTraits)

  # Now, if NATSTAT_ALT available, calculate additional metrics
  if ("NWCA_NATSTAT" %in% names(vascIn)) {
    vascIn.nat <- subset(vascIn, NWCA_NATSTAT == "NAT")

    vascIn.nat.sub <- subset(vascIn.nat, toupper(NWCA_CC) != "UND")

    totals.nat <- aggregate(
      x = list(SUBTOTFREQ = vascIn.nat$FREQ, SUBXTOTABCOV = vascIn.nat$XABCOV),
      by = vascIn.nat[c(sampID)],
      FUN = sum
    )

    numtaxa.nat <- aggregate(
      x = list(NUMTAXA = vascIn.nat.sub$TAXON), by = vascIn.nat.sub[c(sampID)],
      FUN = function(x) {
        length(unique(x))
      }
    )

    vascIn.1a.nat <- merge(subset(vascIn.nat.sub, select = names(vascIn.nat.sub) %nin% c("TOTFREQ", "XTOTABCOV")),
      totals.nat,
      by = sampID
    )
    vascIn.1.nat <- merge(vascIn.1a.nat, numtaxa.nat, by = sampID)
    vascIn.2.nat <- vascIn.1.nat
    vascIn.2.nat$NWCA_CC <- as.numeric(vascIn.2.nat$NWCA_CC)

    vascIn.2.nat$XC_FREQ <- with(vascIn.2.nat, (FREQ / SUBTOTFREQ) * NWCA_CC * 100)
    vascIn.2.nat$FQAI_FREQ <- with(vascIn.2.nat, (FREQ / SUBTOTFREQ) * NWCA_CC * 100)
    vascIn.2.nat$XC_COV <- with(vascIn.2.nat, (XABCOV / SUBXTOTABCOV) * NWCA_CC * 100)
    vascIn.2.nat$FQAI_COV <- with(vascIn.2.nat, (XABCOV / SUBXTOTABCOV) * NWCA_CC * 100)

    vascIn.2a.nat <- aggregate(
      x = list(
        XC = vascIn.2.nat$NWCA_CC, XC_FREQ = vascIn.2.nat$XC_FREQ,
        XC_COV = vascIn.2.nat$XC_COV
      ),
      by = vascIn.2.nat[c(sampID, "NUMTAXA")],
      FUN = sum
    )

    vascIn.2a.nat$XC <- with(vascIn.2a.nat, round(XC / NUMTAXA, 2))
    vascIn.2a.nat$XC_FREQ <- with(vascIn.2a.nat, round(XC_FREQ / NUMTAXA, 2))
    vascIn.2a.nat$XC_COV <- with(vascIn.2a.nat, round(XC_COV / NUMTAXA, 2))
    vascIn.2a.nat$NUMTAXA <- NULL

    vascIn.2b.nat <- aggregate(
      x = list(
        FQAI = vascIn.2.nat$NWCA_CC,
        FQAI_FREQ = vascIn.2.nat$FQAI_FREQ,
        FQAI_COV = vascIn.2.nat$FQAI_COV
      ),
      by = vascIn.2.nat[c(sampID, "NUMTAXA")],
      FUN = sum
    )

    vascIn.2b.nat$FQAI <- with(vascIn.2b.nat, round(FQAI / sqrt(NUMTAXA), 2))
    vascIn.2b.nat$FQAI_FREQ <- with(vascIn.2b.nat, round(FQAI_FREQ / sqrt(NUMTAXA), 2))
    vascIn.2b.nat$FQAI_COV <- with(vascIn.2b.nat, round(FQAI_COV / sqrt(NUMTAXA), 2))
    vascIn.2b.nat$NUMTAXA <- NULL

    vascIn.2c.nat <- merge(vascIn.2a.nat, vascIn.2b.nat, by = sampID)

    ccOut.nat <- reshape(vascIn.2c.nat,
      idvar = sampID, direction = "long",
      varying = names(vascIn.2c.nat)[!names(vascIn.2c.nat) %in% sampID],
      timevar = "variable", v.names = "RESULT",
      times = names(vascIn.2c.nat)[!names(vascIn.2c.nat) %in% sampID]
    )

    ccOut.nat$PARAMETER <- paste(ccOut.nat$variable, "NAT", sep = "_")
    ccOut.nat$variable <- NULL


    ccOut <- rbind(ccOut, ccOut.nat)

    empty_base <- cbind(empty_base, empty_base.nat)
  }


  outdf <- reshape(ccOut,
    idvar = c(sampID), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(outdf) <- gsub("RESULT\\.", "", names(outdf))


  outdf <- merge(outdf, empty_base, all = TRUE)

  outdf.1 <- reshape(outdf,
    idvar = sampID, direction = "long",
    varying = names(outdf)[!names(outdf) %in% c(sampID)],
    timevar = "PARAMETER", v.names = "RESULT",
    times = names(outdf)[!names(outdf) %in% c(sampID)]
  )

  outdf.1 <- subset(outdf.1, !is.na(eval(as.name(sampID[1]))))
  outdf.1$RESULT <- with(outdf.1, ifelse(is.na(RESULT), 0, RESULT))
  outdf.1$PARAMETER <- with(outdf.1, as.character(PARAMETER))

  outdf.1 <- subset(outdf.1, PARAMETER %in% names(empty_base))

  return(outdf.1)
}
