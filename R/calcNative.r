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
#'     categories of 'NAT', 'ADV', 'CRYP', 'INTR', 'UND'
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
#' A list of metric descriptions is provided in the document named 
#' \href{https://github.com/USEPA/aquametNWCA/blob/main/inst/VascPlant_Metric_Descriptions.pdf}{VascPlant_Metric_Descriptions.pdf}
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
#' natEx <- calcNative(exPlant$byUIDspp)
#'
#' head(natEx)
#' unique(natEx$PARAMETER)
calcNative <- function(vascIn, sampID = "UID") {
  vascIn <- as.data.frame(vascIn) # Do this in case read in as a tibble or data.table, which might cause problems
  if ("NWCA_NATSTAT" %nin% names(vascIn)) {
    print("Missing NWCA_NATSTAT from input data frame - cannot calculate metrics! If NWCA_NATSTAT exists,
          run prepareData() function to create input data frame.")
    return(NULL)
  }

  vascIn$ALIEN <- with(vascIn, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV"), 1, 0))
  vascIn$AC <- with(vascIn, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV", "CRYP"), 1, 0))

  sppNATSTAT <- int.calcTraits_MultCat.alt(vascIn, "NWCA_NATSTAT", sampID)

  alienTrait <- int.calcTraits_Indicator.alt(vascIn, "ALIEN", sampID)
  alienTrait$PARAMETER <- with(alienTrait, paste(PARAMETER, "SPP", sep = ""))

  acTrait <- int.calcTraits_Indicator.alt(vascIn, "AC", sampID)

  natstatOut <- rbind(sppNATSTAT, alienTrait, acTrait)

  empty_base <- data.frame(t(rep(NA, 30)), stringsAsFactors = F)
  names(empty_base) <- c(
    "PCTN_ADVSPP", "PCTN_CRYPSPP", "PCTN_INTRSPP", "PCTN_NATSPP", "XABCOV_ADVSPP",
    "XABCOV_CRYPSPP", "XABCOV_INTRSPP", "XABCOV_NATSPP", "XRCOV_ADVSPP",
    "XRCOV_CRYPSPP", "XRCOV_INTRSPP", "XRCOV_NATSPP", "RFREQ_ADVSPP", "RFREQ_CRYPSPP",
    "RFREQ_INTRSPP", "RFREQ_NATSPP", "RIMP_ADVSPP", "RIMP_CRYPSPP", "RIMP_INTRSPP",
    "RIMP_NATSPP", "PCTN_ALIENSPP", "XRCOV_ALIENSPP", "RFREQ_ALIENSPP",
    "RIMP_ALIENSPP", "XABCOV_ALIENSPP", "PCTN_AC", "XRCOV_AC", "RFREQ_AC",
    "RIMP_AC", "XABCOV_AC"
  )


  outdf <- reshape(natstatOut,
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
  outdf.1$RESULT <- with(outdf.1, ifelse(is.na(RESULT) | is.infinite(RESULT), 0, RESULT))
  outdf.1$PARAMETER <- with(outdf.1, as.character(PARAMETER))

  outdf.1 <- subset(outdf.1, PARAMETER %in% names(empty_base))

  return(outdf.1)
}
