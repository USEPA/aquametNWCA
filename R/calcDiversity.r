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
#' 'NAT', 'ADV', 'CRYP', 'INTR', 'UND'
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
#' exPlant <- prepareData(VascPlantEx,
#'   taxon_name = "USDA_NAME",
#'   inTaxa = taxaNWCA, inNat = ccNatNWCA, inCVal = ccNatNWCA,
#'   inWIS = wisNWCA, cValReg = "STATE"
#' )
#'
#' divEx <- calcDiversity(exPlant$byUIDspp)
#'
#' head(divEx)
#' unique(divEx$PARAMETER)
calcDiversity <- function(vascIn, sampID = "UID") {
  vascIn <- as.data.frame(vascIn) # Do this in case read in as a tibble or data.table, which might cause problems
  divOut <- int.calcIndices(vascIn, "ALL", sampID)

  if ("NWCA_NATSTAT" %in% names(vascIn)) {
    vascIn.1 <- vascIn
    vascIn.1$NATSTAT_ALT <- with(vascIn.1, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV"), "ALIEN", NWCA_NATSTAT))
    vascIn.1$AC <- with(vascIn.1, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV", "CRYP"), 1, 0))

    nsvalues <- c("NAT", "ALIEN")
    for (i in 1:length(nsvalues)) {
      vascIn.2 <- subset(vascIn.1, NATSTAT_ALT == nsvalues[i])

      nsOut <- int.calcIndices(vascIn.2, nsvalues[i], sampID)
      divOut <- rbind(divOut, nsOut)
    }

    # AC metrics
    vascIn.ac <- subset(vascIn.1, AC == "1")
    acOut <- int.calcIndices(vascIn.ac, "AC", sampID)

    divOut <- rbind(divOut, acOut)
  }

  divOut.wide <- reshape(divOut,
    idvar = c(sampID), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(divOut.wide) <- gsub("RESULT\\.", "", names(divOut.wide))

  outdf <- reshape(divOut.wide,
    idvar = sampID, direction = "long",
    varying = names(divOut.wide)[!names(divOut.wide) %in% c(sampID)],
    timevar = "PARAMETER", v.names = "RESULT",
    times = names(divOut.wide)[!names(divOut.wide) %in% c(sampID)]
  )

  outdf$RESULT <- with(outdf, ifelse(is.na(RESULT), 0, RESULT))
  outdf$PARAMETER <- as.character(outdf$PARAMETER)

  return(outdf)
}
