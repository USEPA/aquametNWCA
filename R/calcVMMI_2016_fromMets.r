#' @export
#'
#' @title Calculate Vegetation MMI and assign condition class
#'
#' @description This function calculates the NWCA 2016 Vegetation
#' Multimetric Index (VMMI) from metric and wetland class group
#' inputs. If the appropriate
#' variable describing NWCA 2016 reporting units (RPT_UNIT)
#' is included in the input data frame,
#' condition class (Good/Fair/Poor) will also be assigned.
#'
#' @param metsIn Data frame containing, at a minimum:
#' \itemize{
#'   \item sampID: Variables identified in \emph{sampID} argument
#'
#'   \item wetcls_grp: Variable identified in \emph{wetcls_grp} argument,
#'   valid values of this variable are: EH, EW, PRLH, PRLW
#'
#'   \item Metrics required, depending on wetland class
#'   groups in data:
#'    \itemize{
#'     \item EH: N_ANNUAL, PCTN_NATSPP, XRCOV_FORB, XRCOV_HTOL,
#'     XRCOV_MONOCOTS_NAT, XRCOV_SEN
#'
#'     \item EW: PCTN_ISEN, PCTN_MONOCOT, PCTN_NATSPP, RIMP_NATSPP,
#'     XCOV_WD_FINE, XRCOV_GRAMINOID
#'
#'     \item PRLH: FQAI_ALL, N_TOL, PCTN_OBL_FACW, XRCOV_NATSPP
#'
#'     \item RFREQ_NATSPP, XC_ALL, XRCOV_MONOCOTS_NAT, XRCOV_NATSPP
#'
#'    }
#'   \item RPT_UNIT (optional): NWCA 2016 reporting units with the following
#'   valid values: ARW, ATL, GFC, GPL, ICP, NCE, PAC, SAP, TPL, WVM.
#'   Required to assign condition class based on VMMI_2016.
#'   }
#'
#' @param wetcls_grp String containing name of variable containing wetland
#' class group, with valid values of EH (estuarine herbaceous),
#' EW (estuarine woody), PRLH (inland herbaceous), and
#' PRLW (inland woody).
#'
#' @param sampID A character vector containing the name(s) of variable(s)
#'   necessary to identify unique samples, 'UID' by default
#'
#' @return  Data frame containing: \itemize{
#' \item \emph{sampID} Variable(s)
#'   found in the argument sampID
#'
#'   \item All metrics used in wetland class-specific MMIs, with
#'   _SC16 on the end of the metric
#'
#'   \item VMMI_2016: NWCA 2016 Vegetation Multimetric Index on
#'   100-point scale } If
#'   RPT_UNIT is provided, the following is also provided as output:
#'   \itemize{
#'   \item VEGCOND_2016: Vegetation condition class, based on VMMI score
#'   and site location, as assigned for NWCA 2016. Valid values are Good, Fair,
#'   and Poor. }
#'
#' @references US Environmental Protection Agency. 2016. National Wetland
#'   Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
#'   Environmental Protection Agency, Washington, DC.
#'
#' @author Karen Blocksom  \email{blocksom.karen@epa.gov}
#'
#' @examples
#' head(vmmiMetEx)
#'
#' # Test using ECO_X_WETGRP variable in data frame
#' check.1 <- calcVMMI_fromMets(vmmiMetEx)
#'
#' # Now drop ECO_X_WETGRP and use NWCA_ECO4 and NWCA_WET_GRP
#' check.2 <- calcVMMI_fromMets(subset(vmmiMetEx, select = -ECO_X_WETGRP))
#'
#' # Now run without ECO_X_WETGRP, NWCA_ECO4, and NWCA_WET_GRP
#' check.3 <- calcVMMI_fromMets(subset(vmmiMetEx, select = c(
#'   "UID", "FQAI_ALL",
#'   "N_TOL", "RIMP_NATSPP", "XRCOV_MONOCOTS_NAT"
#' )))
calcVMMI_2016_fromMets <- function(metsIn, sampID = "UID", wetcls_grp = NULL) {
  # Look for UID, metrics in input data frame

  necVars <- c(sampID, wetcls_grp)
  # Alert user if any necessary variables are missing
  if (sum(necVars %in% names(metsIn)) < 2) {
    msgVars <- necVars[necVars %nin% names(metsIn)]
    print(paste(paste(msgVars, collapse = ", "), "not found in input data frame. Cannot calculate VMMI without them.", sep = " "))
  }

  necMets <- data.frame(
    wet_grp = c(rep("EH", 6), rep("EW", 6), rep("PRLH", 4), rep("PRLW", 4)),
    METRIC = c(
      "N_ANNUAL", "PCTN_NATSPP", "XRCOV_FORB", "XRCOV_HTOL",
      "XRCOV_MONOCOTS_NAT", "XRCOV_SEN",
      "PCTN_ISEN", "PCTN_MONOCOT", "PCTN_NATSPP", "RIMP_NATSPP",
      "XCOV_WD_FINE", "XRCOV_GRAMINOID",
      "FQAI_ALL", "N_TOL", "PCTN_OBL_FACW", "XRCOV_NATSPP",
      "RFREQ_NATSPP", "XC_ALL", "XRCOV_MONOCOTS_NAT", "XRCOV_NATSPP"
    ),
    stringsAsFactors = F
  )
  names(necMets)[names(necMets) == "wet_grp"] <- wetcls_grp

  # How many different wetcls_grp values are there in the input dataset and then are all
  # of the metrics present?
  wetcls_pres <- unique(metsIn[, wetcls_grp])

  necMets_sub <- unique(subset(necMets, eval(as.name(wetcls_grp)) %in% wetcls_pres, select = "METRIC"))

  if (sum(necMets_sub$METRIC %in% names(metsIn)) < nrow(necMets_sub)) {
    msgMets <- necMets_sub$METRIC[necMets_sub$METRIC %nin% names(metsIn)]
    print(paste(
      paste(msgMets, collapse = ", "),
      "not found in input data frame. Cannot calculate VMMI without them."
    ))
  }

  # Look for region and wetland type variables - cannot assign condition without them
  # Can be individual variables that will be used to combine into ECO_X_WETGRP or ECO_X_WETGRP
  # Currently, we require the variable names to match those from NWCA 2011 to ensure they
  # contain the correct values
  if ("RPT_UNIT" %nin% names(metsIn) & "PRLW" %in% wetcls_pres) {
    print("Warning: Must include variable RPT_UNIT for NWCA 2016 to determine condition class for PRLW sites. This variable is missing, and condition class will not be assigned, but VMMI will be calculated.")
  }

  # Use RPT_UNIT
  if ("RPT_UNIT" %in% names(metsIn) & "PRLW" %in% wetcls_pres) {
    metsIn$RPT_FINAL <- with(metsIn, ifelse(eval(as.name(wetcls_grp)) %in% c("PRLW") &
      RPT_UNIT %in% c("ARW", "GPL", "TPL"),
    "PRLW-PLN_ARIDW",
    ifelse(eval(as.name(wetcls_grp)) == "PRLW" &
      RPT_UNIT %nin% c("ARW", "GPL", "TPL"), "PRLW-OTHER",
    eval(as.name(wetcls_grp))
    )
    ))
  } else {
    metsIn$RPT_FINAL <- metsIn[, wetcls_grp]
  }

  # Identify key variables in input dataset
  if ("RPT_FINAL" %in% names(metsIn)) {
    keyVars <- c(sampID, wetcls_grp, "RPT_FINAL")
  } else {
    keyVars <- c(sampID, wetcls_grp)
  }

  # Make sure all metrics are in numeric format for scoring
  metsIn[, necMets_sub$METRIC] <- lapply(metsIn[, necMets_sub$METRIC], as.numeric)

  metsIn <- metsIn[, c(keyVars, necMets_sub$METRIC)]

  # Melt data frame into long format
  metsIn.long <- reshape(metsIn,
    idvar = keyVars, direction = "long",
    varying = c(necMets_sub$METRIC),
    timevar = "PARAMETER", v.names = "RESULT",
    times = c(necMets_sub$METRIC)
  )

  # Set metric scoring thresholds
  metTholds <- data.frame(
    wetcls_grp = c(
      rep("EH", 6), rep("EW", 6),
      rep("PRLH", 4), rep("PRLW", 4)
    ),
    PARAMETER = c(
      "N_ANNUAL", "PCTN_NATSPP", "XRCOV_FORB",
      "XRCOV_HTOL", "XRCOV_MONOCOTS_NAT", "XRCOV_SEN",
      "PCTN_ISEN", "PCTN_MONOCOT", "PCTN_NATSPP",
      "RIMP_NATSPP", "XCOV_WD_FINE", "XRCOV_GRAMINOID",
      "FQAI_ALL", "N_TOL", "PCTN_OBL_FACW", "XRCOV_NATSPP",
      "RFREQ_NATSPP", "XC_ALL", "XRCOV_MONOCOTS_NAT",
      "XRCOV_NATSPP"
    ),
    CEILING = c(
      2, 100, 69.367, 84.5735, 100, 100,
      45.45, 55.6805, 100, 100, 13.85, 90.3625,
      35.7735, 41, 100, 100,
      100, 6.19, 48.414, 100
    ),
    FLOOR = c(
      0, 62.959, 0, 0, 0.288, 0,
      7.569, 0, 66.982, 68.562, 0, 0,
      4.9045, 3, 17.2055, 12.4155,
      62.832, 2.524, 0.17, 53.042
    ),
    DIRECTION = c(
      "NEG", "POS", "NEG", "NEG", "POS", "POS",
      "POS", "POS", "POS", "POS", "NEG", "POS",
      "POS", "NEG", "POS", "POS",
      rep("POS", 4)
    ),
    stringsAsFactors = FALSE
  )

  # Merge metric scoring thresholds with metrics in long format
  vMet <- merge(metsIn.long, metTholds,
    by.x = c(wetcls_grp, "PARAMETER"),
    by.y = c("wetcls_grp", "PARAMETER")
  )

  # Now apply the function that calculates metric scores by interpolating values between 0 and 10
  scoreMet <- function(dir, x, floor, ceiling) {
    if (dir == "POS") { # if positive metric, interpolate with first approach
      zz <- round(approx(x = c(floor, ceiling), y = c(0, 10), xout = x, method = "linear", yleft = 0, yright = 10)$y, 2)
    } else { # if negative metric, interpolate score with this approach
      zz <- round(approx(x = c(floor, ceiling), y = c(10, 0), xout = x, method = "linear", yleft = 10, yright = 0)$y, 2)
    }
  }

  # Calculate scores and add scored version of  metric (METRIC_SC) to data frame. SC = rescaled
  # metric score that is used in MMI calculations
  # Keep only variables that are needed
  scored.mets <- vMet[, c(keyVars, "PARAMETER")]
  # Use scoreMet() function to calculate RESULT value for each metric
  scored.mets$RESULT <- with(scored.mets, with(vMet, mapply(scoreMet, DIRECTION, RESULT, FLOOR, CEILING)))
  # Add _SC to end of each metric name as name of scored metric
  scored.mets$PARAMETER <- with(scored.mets, paste(as.character(PARAMETER), "SC16", sep = "_"))

  # Now that we have scored metrics, we can calculate MMI scores and merge with MMI
  # thresholds to determine condition
  # Sum metrics into MMI
  mmi.1 <- aggregate(
    x = list(VMMI_2016 = scored.mets$RESULT), by = scored.mets[keyVars],
    FUN = sum
  )

  numMets <- data.frame(
    wetcls_grp = c("EH", "EW", "PRLH", "PRLW"),
    nummets = c(6, 6, 4, 4), stringsAsFactors = F
  )

  mmi.1 <- merge(mmi.1, numMets, by.x = wetcls_grp, by.y = "wetcls_grp")
  # Rescale VMMI to 100-point scale, rounding to 1 digit
  mmi.1$VMMI_2016 <- with(mmi.1, round(VMMI_2016 * (10 / nummets), 1))
  # Melt data frame into long format
  mmi <- reshape(mmi.1,
    idvar = keyVars, direction = "long",
    varying = "VMMI_2016", timevar = "PARAMETER", v.names = "RESULT",
    times = "VMMI_2016"
  )

  mmi$nummets <- NULL
  # Combine metric scores and VMMI scores
  mmiOut <- rbind(scored.mets, mmi)

  # Set VMMI condition class thresholds
  mmiTholds <- data.frame(
    RPT_FINAL = c("EH", "EW", "PRLH", "PRLW-OTHER", "PRLW-PLN_ARIDW"),
    p05 = c(73.6, 64.6, 63.8, 53.7, 43.7),
    p25 = c(86.4, 69.8, 74.2, 65.5, 49.9),
    stringsAsFactors = FALSE
  )

  # If ECO_X_WETGRP exists in mmi data frame, assign condition
  if ("RPT_FINAL" %in% names(mmi)) {
    mmi.1 <- merge(mmi, mmiTholds, by = "RPT_FINAL")

    cond <- mmi.1
    cond$VEGCOND_2016 <- with(mmi.1, ifelse(RESULT >= p25, "Good", ifelse(RESULT >= p05, "Fair", "Poor")))
    # Melt condition data frame to combine with metric and VMMI scores
    cond.long <- reshape(cond,
      idvar = keyVars, direction = "long",
      varying = "VEGCOND_2016", times = "VEGCOND_2016",
      timevar = "PARAMETER", v.names = "RESULT"
    )
    cond.long <- subset(cond.long, select = c(keyVars, "PARAMETER", "RESULT"))
    # Combine with metric and VMMI scores
    mmiOut <- rbind(mmiOut, cond.long)
  }
  # Cast output data frame into wide format, then remove prefixes
  # added by reshape() from variable names
  mmiOut.wide <- reshape(mmiOut,
    idvar = c(keyVars), direction = "wide",
    timevar = "PARAMETER", v.names = "RESULT"
  )
  names(mmiOut.wide) <- gsub("RESULT\\.", "", names(mmiOut.wide))

  return(mmiOut.wide)
}
