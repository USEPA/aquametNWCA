#' aquametNWCA: A package to calculate vegetation metrics and the 
#' vegetation multimetric index (VMMI) used in NWCA
#' 
#' @docType package
#' @name aquametNWCA
#' 
#' @importFrom Hmisc "%nin%"
#' @importFrom gtools smartbind
#' @importFrom ecodist distance
#' @importFrom stats approx median sd aggregate reshape
#' @importFrom rlang .data
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
if(getRversion() >= "3.0") utils::globalVariables(c('COVER'
 ,'SPECIES_NAME_ID','NWCA_NATSTAT','PARAMETER','RESULT'
 ,'AC','CATEGORY','COV','DISTINCT','DURATION','DUR_ALT'
 ,'ECOIND','ECOIND1','ECOIND2','FREQ','FREQ_H2O'
 ,'FREQ_LITTER','GRH_ALT'
 ,'GROWTH_HABIT','H','HERB','H_SANDT','H_VASC_STRATA'
 ,'J_VASC_STRATA','LITTER_TYPE','MAXF','MAXIMUM_DEPTH'
 ,'MAXN','METRIC','MINIMUM_DEPTH','N','NATSTAT_ALT','NPLOTS'
 ,'NQUADS','NSAMP','NUM','NWCA_CC','NWCA_ECO4'
 ,'NWCA_WET_GRP','N_PEAT_MOSS_DOM','N_SANDT','N_TAXA'
 ,'N_TREESPP','N_VSTRATA','PAL_FARMED','PARAM_ALT','PLOT'
 ,'PLOTSAMP','PREDOMINANT_DEPTH','RFREQ','SAMPID'
 ,'SANDT_CLASS','SHRUB_COMB','SUBTOTFREQ','SUBXTOTABCOV'
 ,'TAXON','TOL','TOTAL_WATER','TOTFREQ','TOTN','TREE_COMB'
 ,'TREE_SPECIES','USDA_NAME','WIS','XABCOV','XCOV','XCOV_H2O'
 ,'XCOV_LITTER','XDEPTH_LITTER','XN','XRCOV','XTOTABCOV'
 ,'XTOTCOV_VASC_STRATA'
 ,'ccNatNWCA','p05','p25','sRFREQ','sXRCOV','taxaNWCA'
 ,'tobj','value','variable','wisNWCA'))
 
