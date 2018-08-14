#' aquametNWCA: A package to calculate vegetation metrics and the 
#' vegetation multimetric index (VMMI) used in NWCA 2011
#' 
#' @docType package
#' @name aquametNWCA
#' 
#' @importFrom Hmisc "%nin%"
#' @importFrom gtools smartbind
#' @importFrom dplyr filter select "%>%"
#' @importFrom plyr ddply mutate summarise summarize rename
#' @importFrom reshape2 dcast melt
#' @importFrom car recode Recode
#' @importFrom ecodist distance
#' @importFrom stats approx median sd
#' 
#' @keywords package
#' @title aquametNWCA
#' 
#' @section Vascular plant functions:
#' These functions calculate either subsets of metrics or the
#' full set of metrics, based on vascular plant data:
#' \itemize{
#' \item calcBareGround_LitterMets()
#' \item calcBCmets()
#' \item calcCategory()
#' \item calcCC()
#' \item calcDiversity()
#' \item calcDuration()
#' \item calcGrowthHabit()
#' \item calcNative()
#' \item calcRichness()
#' \item calcVascPlantMets()
#' \item calcWIS()
#' }
#' 
#' @section Tree metric functions:
#' These functions calculate either subsets of metrics or the 
#' full set of metrics, based on tree data:
#' \itemize{ 
#' \item calcSnagMets()
#' \item calcTreeCntMets()
#' \item calcTreeCoverMets()
#' \item calcTreeMets()
#' }
#' 
#' @section Vascular strata and ground cover functions:
#' These functions calculate either subsets of metrics or the
#' full set of metrics, based on ground cover or vascular strata 
#' data:
#' \itemize{
#' \item calcBareGround_LitterMets()
#' \item calcNonvascMets()
#' \item calcSandTMets()
#' \item calcVascStratMets()
#' \item calcVtype_GcovMets()
#' \item calcWcovMets()
#' }
#' 
#' @section MMI calculation:
#' These functions calculate the MMI from metrics, and 
#' also only MMI metrics:
#' \itemize{
#' \item calcVMMI_fromMets()
#' \item calcVMMImets()
#' }
#' 
#' @section Data preparation:
#' These functions prepare input datasets for the functions
#' listed above, in the expected formats, from raw data:
#' \itemize{
#' \item createDFs()
#' \item prepareData()
#' \item nwcaVegData() (alternative to prepareData() that allows
#' use of a C-value region other than 'STATE')
#' \item nwcaVegInput() (alternative to createDFs() that allows
#' use of a C-value region other than 'STATE')
#' }
#' 
#' @section Included datasets used in examples:
#' \itemize{
#' \item taxaNWCA - NWCA 2011 plant taxa list, including category,
#' duration, and growth habit, and taxonomy
#'
#' \item ccNatNWCA - NWCA 2011 plant list of Coefficient of Conservatism
#' values by state
#' 
#' \item wisNWCA - NWCA 2011 plant list of Wetland Indicator Status by
#' U.S. Army Corps of Engineers region
#' 
#' \item VascPlantEx - Example vascular plant dataset
#' 
#' \item TreesEx - Example tree dataset
#' 
#' \item Vtyp_GrCovEx - Example vegetation type (strata) and ground cover dataset
#' 
#' \item vmmiMetsEx - Example dataset containing metrics necessary to calculate 
#' the VMMI 
#' }
#' 
NULL
## NULL