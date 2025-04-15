# Growth habit metrics (with and without native status, if available)
#' @export
#' @title Calculate growth habit metrics
#' @description This function calculates all growth habit metrics, with
#' additional metrics if NWCA_NATSTAT variable is included
#' in input data frame.
#' @details Both GROWTH_HABIT and NWCA_NATSTAT variables are recoded to fewer
#' categories. Taxa with 'UND' as native status are excluded.
#' @param vascIn   Data frame containing cover data summarized by
#' UID and TAXON, with the following fields:
#' \itemize{
#'     \item sampID: Variable(s) identified in \emph{sampID} argument
#'
#'     \item TAXON: Taxon name
#'
#'     \item GROWTH_HABIT: USDA PLANTS GROWTH_HABIT variable with valid
#'     values FORB, GRAMINOID, HERB, SHRUB, SUBSHRUB, TREE, VINE,
#'     NONVASCULAR, combinations of these, or blank
#'
#'     \item GRH_ALT (optional): Combinations of GROWTH_HABIT variable
#'     as used in NWCA. Created in code if not supplied by user. Valid
#'     values include FORB, GRAMINOID, SHRUB, SSHRUB_FORB, SSHRUB_SHRUB,
#'     TREE, TREE_SHRUB, VINE, VINE_SHRUB.
#'
#'     \item TREE_COMB (optional): Indicator (1/0) for values of
#'     TREE or TREE_SHRUB in GRH_ALT.
#'     If not supplied by user, will be created in code
#'
#'     \item SHRUB_COMB (optional): Indicator value (0/1) for GRH_ALT
#'     values of SHRUB, SSHRUB_FORB, or SSHRUB_SHRUB. Created in code
#'     if not supplied by user.
#'
#'     \item HERB (optional): Indicator value (0/1) for GRH_ALT values
#'     of GRAMINOID or FORB. Created in code if not supplied by user.
#'
#'     \item XABCOV: Mean percent cover of taxon across plots
#'
#'     \item TOTN: Number of taxa in sample
#'
#'     \item sXRCOV: proportion of summed cover across all taxa
#'     (XTOTABCOV) represented by taxon in sample
#'
#'     \item NWCA_NATSTAT (optional): Native status variable with categories
#'     of 'NAT', 'ADV', 'CRYP', 'INTR', 'UND'
#'    }
#' @param sampID  A character vector containing the name(s) of variable(s)
#'   necessary to identify unique samples, 'UID' by default
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
#' If NWCA_NATSTAT is in the input data frame, the same set except
#' 'N_' metrics is calculated for native species and
#' alien + cryptogenic species, with metric names suffixes of
#' _NAT and _AC, respectively.
#' 
#' #' A list of metric descriptions is provided in the document named 
#' \href{https://github.com/USEPA/aquametNWCA/blob/main/inst/VascPlant_Metric_Descriptions.pdf}{VascPlant_Metric_Descriptions.pdf}
#' @author Karen Blocksom \email{Blocksom.karen@epa.gov}
#' @references US Environmental Protection Agency. 2016. National Wetland
#' Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#' Environmental Protection Agency, Washington, DC.
#' @examples
#' head(VascPlantEx)
#' exPlant <- prepareData(VascPlantEx,
#'   taxon_name = "USDA_NAME",
#'   inTaxa = taxaNWCA, inNat = ccNatNWCA, inCVal = ccNatNWCA,
#'   inWIS = wisNWCA, cValReg = "STATE"
#' )
#'
#' ghEx <- calcGrowthHabit(exPlant$byUIDspp)
#'
#' head(ghEx)
#' unique(ghEx$PARAMETER)
calcGrowthHabit <- function(vascIn, sampID = "UID") {
  vascIn <- as.data.frame(vascIn) # Do this in case read in as a tibble or data.table, which might cause problems
  # First check for necessary variables
  necNames <- c(sampID, "TAXON", "GROWTH_HABIT", "XABCOV", "TOTN", "sXRCOV")
  msgNames <- necNames[necNames %nin% names(vascIn)]
  if (length(msgNames) > 0) {
    print(paste("Missing key variables for metric calculation: ", paste(msgNames, collapse = ","),
      ". Try prepareData() function to create necessary input variables.",
      sep = ""
    ))
    return(NULL)
  }

  vascIn <- subset(vascIn, !is.na(GROWTH_HABIT) & toupper(GROWTH_HABIT) != "UND" & GROWTH_HABIT != "")

  # Check for GRH_ALT variable
  if ("GRH_ALT" %in% names(vascIn)) {
    vascIn.1 <- vascIn

    vascIn.1$GRH_ALT <- with(vascIn.1, gsub("SUBSHRUB", "SSHRUB", GRH_ALT))
  } else {
    vascIn.1 <- vascIn
    vascIn.1$GRH_ALT <- vascIn.1$GROWTH_HABIT
    vascIn.1$GRH_ALT[vascIn.1$GROWTH_HABIT %in% c(
      "FORB/HERB",
      "FORB/HERB, SHRUB",
      "FORB/HERB, SHRUB, SUBSHRUB",
      "FORB/HERB, SUBSHRUB",
      "FORB/HERB, SUBSHRUB, SHRUB",
      "FORB/HERB, SUBSHRUB, TREE"
    )] <- "FORB"
    vascIn.1$GRH_ALT[vascIn.1$GROWTH_HABIT %in% c(
      "SUBSHRUB, FORB/HERB",
      "SUBSHRUB, SHRUB, FORB/HERB",
      "SUBSHRUB, FORB/HERB, SHRUB",
      "SHRUB, SUBSHRUB, FORB/HERB"
    )] <- "SSHRUB_FORB"
    vascIn.1$GRH_ALT[vascIn.1$GROWTH_HABIT %in% c(
      "SUBSHRUB",
      "SUBSHRUB, SHRUB",
      "SHRUB, SUBSHRUB",
      "SUBSHRUB, SHRUB, TREE",
      "SHRUB, FORB/HERB, SUBSHRUB"
    )] <- "SSHRUB_SHRUB"
    vascIn.1$GRH_ALT[vascIn.1$GROWTH_HABIT %in% c(
      "SHRUB",
      "SHRUB, TREE",
      "TREE, SUBSHRUB, SHRUB",
      "SHRUB, SUBSHRUB, TREE",
      "SUBSHRUB, FORB/HERB, SHRUB, TREE",
      "FORB/HERB, SHRUB, SUBSHRUB, TREE",
      "SHRUB,"
    )] <- "SHRUB"
    vascIn.1$GRH_ALT[vascIn.1$GROWTH_HABIT %in% c(
      "TREE, SHRUB",
      "TREE, SHRUB, VINE",
      "TREE, SHRUB, SUBSHRUB"
    )] <- "TREE_SHRUB"
    vascIn.1$GRH_ALT[vascIn.1$GROWTH_HABIT %in% c(
      "VINE, FORB/HERB",
      "SUBSHRUB, FORB/HERB, VINE",
      "FORB/HERB, VINE",
      "FORB/HERB, VINE, SUBSHRUB",
      "VINE, FORB/HERB, SUBSHRUB",
      "VINE, HERBACEOUS"
    )] <- "VINE"
    vascIn.1$GRH_ALT[vascIn.1$GROWTH_HABIT %in% c(
      "VINE, SHRUB",
      "VINE, SUBSHRUB",
      "SUBSHRUB, VINE",
      "SHRUB, VINE",
      "SHRUB, FORB/HERB, SUBSHRUB, VINE",
      "SHRUB, SUBSHRUB, VINE",
      "VINE, TREE, SHRUB",
      "SHRUB, VINE, FORB/HERB",
      "SUBSHRUB, VINE, SHRUB",
      "SUBSHRUB, VINE, FORB/HERB",
      "SHRUB, VINE, SUBSHRUB",
      "SHRUB, FORB/HERB, VINE",
      "VINE, WOODY",
      "FORB/HERB, SHRUB, SUBSHRUB, TREE, VINE",
      "SHRUB, SUBSHRUB, TREE, VINE",
      "FORB/HERB, SHRUB, SUBSHRUB, VINE"
    )] <- "VINE_SHRUB"
    vascIn.1$GRH_ALT[vascIn.1$GROWTH_HABIT %in% c(
      "GRAMINOID",
      "SUBSHRUB, SHRUB, GRAMINOID",
      "GRAMINOID, SHRUB, SUBSHRUB",
      "GRAMINOID, SHRUB, VINE",
      "SUBSHRUB, GRAMINOID, SHRUB"
    )] <- "GRAMINOID"
  }
  # Check for TREE_COMB variable
  if (("TREE_COMB" %in% names(vascIn)) == FALSE) {
    vascIn.1$TREE_COMB <- with(vascIn.1, ifelse(GRH_ALT %in% c("TREE", "TREE_SHRUB"), 1, 0))
  }
  # Check for SHRUB_COMB variable
  if (("SHRUB_COMB" %in% names(vascIn)) == FALSE) {
    vascIn.1$SHRUB_COMB <- with(vascIn.1, ifelse(GRH_ALT %in% c("SSHRUB_SHRUB", "SHRUB"), 1, 0))
  }
  # Check for HERB variable
  if (("HERB" %in% names(vascIn)) == FALSE) {
    vascIn.1$HERB <- with(vascIn.1, ifelse(GRH_ALT %in% c("GRAMINOID", "FORB", "SSHRUB_FORB"), 1, 0))
  }
  # Combine VINE and VINE_SHRUB
  if (("VINE_ALL" %in% names(vascIn)) == FALSE) {
    vascIn.1$VINE_ALL <- with(vascIn.1, ifelse(GRH_ALT %in% c("VINE", "VINE_SHRUB"), 1, 0))
  }

  sppGRH <- int.calcTraits_MultCat(vascIn.1, "GRH_ALT", sampID)
  sppGRH.1 <- int.combTraits(vascIn.1, c("TREE_COMB", "SHRUB_COMB", "HERB", "VINE_ALL"), sampID)
  grhOut <- rbind(sppGRH, sppGRH.1)

  empty_base <- data.frame(t(rep(NA, 52)), stringsAsFactors = F)
  names(empty_base) <- c(
    "N_FORB", "N_GRAMINOID", "N_SHRUB", "N_SSHRUB_FORB", "N_SSHRUB_SHRUB",
    "N_TREE", "N_TREE_SHRUB", "N_VINE", "N_VINE_SHRUB", "PCTN_FORB",
    "PCTN_GRAMINOID", "PCTN_SHRUB", "PCTN_SSHRUB_FORB", "PCTN_SSHRUB_SHRUB", "PCTN_TREE",
    "PCTN_TREE_SHRUB", "PCTN_VINE", "PCTN_VINE_SHRUB", "XABCOV_FORB", "XABCOV_GRAMINOID",
    "XABCOV_SHRUB", "XABCOV_SSHRUB_FORB", "XABCOV_SSHRUB_SHRUB", "XABCOV_TREE", "XABCOV_TREE_SHRUB",
    "XABCOV_VINE", "XABCOV_VINE_SHRUB", "XRCOV_FORB", "XRCOV_GRAMINOID", "XRCOV_SHRUB",
    "XRCOV_SSHRUB_FORB", "XRCOV_SSHRUB_SHRUB", "XRCOV_TREE", "XRCOV_TREE_SHRUB", "XRCOV_VINE",
    "XRCOV_VINE_SHRUB", "N_TREE_COMB", "PCTN_TREE_COMB", "XABCOV_TREE_COMB", "XRCOV_TREE_COMB",
    "N_SHRUB_COMB", "PCTN_SHRUB_COMB", "XABCOV_SHRUB_COMB", "XRCOV_SHRUB_COMB", "N_HERB",
    "PCTN_HERB", "XABCOV_HERB", "XRCOV_HERB", "N_VINE_ALL",
    "PCTN_VINE_ALL", "XABCOV_VINE_ALL", "XRCOV_VINE_ALL"
  )

  empty_base.nat <- data.frame(t(rep(NA, 64)), stringsAsFactors = F)
  names(empty_base.nat) <- c(
    "N_GRAMINOID_AC", "PCTN_GRAMINOID_AC",
    "XABCOV_GRAMINOID_AC", "XRCOV_GRAMINOID_AC", "N_GRAMINOID_NAT", "PCTN_GRAMINOID_NAT", "XABCOV_GRAMINOID_NAT",
    "XRCOV_GRAMINOID_NAT", "N_FORB_AC", "PCTN_FORB_AC", "XABCOV_FORB_AC", "XRCOV_FORB_AC",
    "N_FORB_NAT", "PCTN_FORB_NAT", "XABCOV_FORB_NAT", "XRCOV_FORB_NAT", "N_HERB_AC",
    "PCTN_HERB_AC", "XABCOV_HERB_AC", "XRCOV_HERB_AC", "N_HERB_NAT", "PCTN_HERB_NAT",
    "XABCOV_HERB_NAT", "XRCOV_HERB_NAT", "N_SHRUB_COMB_AC", "PCTN_SHRUB_COMB_AC", "XABCOV_SHRUB_COMB_AC",
    "XRCOV_SHRUB_COMB_AC", "N_SHRUB_COMB_NAT", "PCTN_SHRUB_COMB_NAT", "XABCOV_SHRUB_COMB_NAT", "XRCOV_SHRUB_COMB_NAT",
    "N_TREE_COMB_AC", "PCTN_TREE_COMB_AC", "XABCOV_TREE_COMB_AC", "XRCOV_TREE_COMB_AC", "N_TREE_COMB_NAT",
    "PCTN_TREE_COMB_NAT", "XABCOV_TREE_COMB_NAT", "XRCOV_TREE_COMB_NAT", "N_VINE_AC", "PCTN_VINE_AC",
    "XABCOV_VINE_AC", "XRCOV_VINE_AC", "N_VINE_NAT", "PCTN_VINE_NAT", "XABCOV_VINE_NAT",
    "XRCOV_VINE_NAT", "N_VINE_SHRUB_AC", "PCTN_VINE_SHRUB_AC", "XABCOV_VINE_SHRUB_AC", "XRCOV_VINE_SHRUB_AC",
    "N_VINE_SHRUB_NAT", "PCTN_VINE_SHRUB_NAT", "XABCOV_VINE_SHRUB_NAT", "XRCOV_VINE_SHRUB_NAT",
    "N_VINE_ALL_AC", "PCTN_VINE_ALL_AC", "XABCOV_VINE_ALL_AC", "XRCOV_VINE_ALL_AC",
    "N_VINE_ALL_NAT", "PCTN_VINE_ALL_NAT", "XABCOV_VINE_ALL_NAT", "XRCOV_VINE_ALL_NAT"
  )


  if ("NWCA_NATSTAT" %in% names(vascIn.1)) {
    vascIn.2 <- vascIn.1
    vascIn.2$ALIEN <- with(vascIn.2, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV"), 1, 0))
    vascIn.2$NATSTAT_ALT <- with(vascIn.2, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV"), "ALIEN", NWCA_NATSTAT))
    vascIn.2$AC <- with(vascIn.2, ifelse(NWCA_NATSTAT %in% c("INTR", "ADV", "CRYP"), 1, 0))

    vascIn.2$GRAMINOID_AC <- with(vascIn.2, ifelse(GRH_ALT == "GRAMINOID" & AC == 1, 1, 0))
    vascIn.2$GRAMINOID_NAT <- with(vascIn.2, ifelse(GRH_ALT == "GRAMINOID" & NATSTAT_ALT == "NAT", 1, 0))
    vascIn.2$FORB_AC <- with(vascIn.2, ifelse(GRH_ALT == "FORB" & AC == 1, 1, 0))
    vascIn.2$FORB_NAT <- with(vascIn.2, ifelse(GRH_ALT == "FORB" & NATSTAT_ALT == "NAT", 1, 0))
    vascIn.2$HERB_AC <- with(vascIn.2, ifelse(HERB == 1 & AC == 1, 1, 0))
    vascIn.2$HERB_NAT <- with(vascIn.2, ifelse(HERB == 1 & NATSTAT_ALT == "NAT", 1, 0))
    vascIn.2$SHRUB_COMB_AC <- with(vascIn.2, ifelse(SHRUB_COMB == 1 & AC == 1, 1, 0))
    vascIn.2$SHRUB_COMB_NAT <- with(vascIn.2, ifelse(SHRUB_COMB == 1 & NATSTAT_ALT == "NAT", 1, 0))
    vascIn.2$TREE_COMB_AC <- with(vascIn.2, ifelse(TREE_COMB == 1 & AC == 1, 1, 0))
    vascIn.2$TREE_COMB_NAT <- with(vascIn.2, ifelse(TREE_COMB == 1 & NATSTAT_ALT == "NAT", 1, 0))
    vascIn.2$VINE_AC <- with(vascIn.2, ifelse(GRH_ALT == "VINE" & AC == 1, 1, 0))
    vascIn.2$VINE_NAT <- with(vascIn.2, ifelse(GRH_ALT == "VINE" & NATSTAT_ALT == "NAT", 1, 0))
    vascIn.2$VINE_SHRUB_AC <- with(vascIn.2, ifelse(GRH_ALT == "VINE_SHRUB" & AC == 1, 1, 0))
    vascIn.2$VINE_SHRUB_NAT <- with(vascIn.2, ifelse(GRH_ALT == "VINE_SHRUB" & NATSTAT_ALT == "NAT", 1, 0))
    vascIn.2$VINE_ALL_AC <- with(vascIn.2, ifelse(VINE_ALL == 1 & AC == 1, 1, 0))
    vascIn.2$VINE_ALL_NAT <- with(vascIn.2, ifelse(VINE_ALL == 1 & NATSTAT_ALT == "NAT", 1, 0))


    multTraits <- int.combTraits(
      vascIn.2, c(
        "GRAMINOID_AC", "GRAMINOID_NAT", "FORB_AC", "FORB_NAT",
        "HERB_AC", "HERB_NAT", "SHRUB_COMB_AC", "SHRUB_COMB_NAT",
        "TREE_COMB_AC", "TREE_COMB_NAT", "VINE_AC", "VINE_NAT",
        "VINE_SHRUB_AC", "VINE_SHRUB_NAT",
        "VINE_ALL_AC", "VINE_ALL_NAT"
      ),
      sampID
    )

    grhOut <- rbind(grhOut, multTraits)

    empty_base <- cbind(empty_base, empty_base.nat)
  }

  grhOut.1a <- reshape(grhOut,
    idvar = c(sampID), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )

  names(grhOut.1a) <- gsub("RESULT\\.", "", names(grhOut.1a))
  grhOut.1a <- merge(grhOut.1a, empty_base, all = TRUE)

  grhOut.1b <- reshape(grhOut.1a,
    idvar = sampID, direction = "long",
    varying = names(grhOut.1a)[!names(grhOut.1a) %in% c(sampID)],
    timevar = "PARAMETER", v.names = "RESULT",
    times = names(grhOut.1a)[!names(grhOut.1a) %in% c(sampID)]
  )
  grhOut.1b <- subset(grhOut.1b, !is.na(eval(as.name(sampID[1]))))
  grhOut.1b$RESULT <- with(grhOut.1b, ifelse(is.na(RESULT), 0, RESULT))
  grhOut.1b$PARAMETER <- with(grhOut.1b, as.character(PARAMETER))

  grhOut.1 <- subset(grhOut.1b, PARAMETER %in% names(empty_base))

  return(grhOut.1)
}
