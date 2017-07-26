#' NWCA 2011 Plant CC and Native Status
#' 
#' A dataset containing NWCA 2011 state-specific Coefficient 
#' of Conservatism and native status values for vascular plant taxa
#' 
#' @name ccNatNWCA
#' @format A data frame with 13203 observations on the following 5 variables:
#' \describe{
#'  \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA 
#'  taxonomic ID number}
#'  \item{USDA_NAME}{Accepted USDA PLANTS name, or the official NWCA name 
#'  where not available}
#'  \item{GEOG_ID}{Geographic unit to which value applies, which is STATE 
#'  in this dataset}
#'  \item{NWCA_CC}{State-specific Coefficient of Conservatism value as 
#'  determined for NWCA 2011}
#'  \item{NWCA_NATSTAT}{State-specific Native status as determined for NWCA 2011}
#' }
#' @note This dataset is the taxa list used for state-specific CC and native 
#' status values used in NWCA 2011.
#' @keywords datasets
#' @examples 
#' head(ccNatNWCA)
#' str(ccNatNWCCA)
#' 
"ccNatNWCA"

#' NWCA 2011 Vascular Plant Taxa List
#' 
#' A dataset containing taxonomy and basic characteristics for vascular plant 
#' taxa as used in NWCA 2011.
#' 
#' @name taxaNWCA
#' 
#' @format A data frame with 3814 observations on the following 12 variables.
#' \describe{
#'   \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA taxonomic 
#'   ID number}
#'   \item{USDA_NAME}{Accepted USDA PLANTS name, or the official NWCA name 
#'     where not available}
#'   \item{DIVISION}{Division level of taxonomy}  
#'   \item{ORDER}{Order level of taxonomy}  
#'   \item{FAMILY}{Family level of taxonomy}
#'   \item{GENUS}{Genus level of taxonomy}  
#'   \item{SPECIES}{Species level of taxonomy}
#'   \item{SUBSPECIES}{Subspecies level of taxonomy}
#'   \item{VARIETY}{Variety level of taxonomy}
#'   \item{CATEGORY}{Plant category of taxon, including MONOCOT, GYMNOSPERM, 
#'   DICOT, FERN, LICHEN, HORSETAIL, MOSS, LYCOPOD, and LIVERWORT}
#'   \item{GROWTH_HABIT}{Growth habit of taxon}
#'   \item{DURATION}{Lifecycle duration of taxon}
#'   }
#'   @note This dataset is the taxa list used for state-specific CC and 
#'   native status values used in NWCA 2011.
#'   
#'   @examples 
#'   head(taxaNWCA)
#'   str(taxaNWCA)
#'   
#'   @keywords datasets
"taxaNWCA"

#' NWCA 2011 Wetland Indicator Status 
#' 
#' A dataset containing NWCA 2011 USAC-specific Wetland Indicator Status 
#' values for vascular plant taxa
#' 
#' @name wisNWCA
#' @format A data frame with 19818 observations on the following 5 variables.
#' \describe{
#'     \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA taxonomic 
#'     ID number}
#'     \item{USDA_NAME}{Accepted USDA PLANTS name, or the official NWCA name 
#'     where not available}
#'     \item{GEOG_ID}{Geographic unit to which value applies, which is US Army 
#'     Corps region in this dataset}
#'     \item{WIS}{USAC-specific Wetland Indicator Status as used for NWCA 2011}
#' }
#' @note This dataset is the taxa list used for USAC-specific Wetland Indicator 
#' Status values used in NWCA 2011.
#' 
#' @examples
#' head(wisNWCA)
#' str(wisNWCA)
#' 
#' @keywords datasets
"wisNWCA"

#' Example Tree data
#' 
#' A dataset containing tree data for use in treeMets() example.
#' 
#' @name TreesEx
#' @format A data frame with 1528 observations on the following 6 variables.
#'   \describe{
#'   \item{UID}{A unique site visit ID}
#'   \item{PLOT}{Plot number of measurement}
#'   \item{PAGE}{Page of field form for measurement}
#'   \item{LINE}{Line number of measurement}
#'   \item{PARAMETER}{Name of tree parameter measured}
#'   \item{RESULT}{Measured value of parameter}
#'   }
#'   
#' @note This dataset is a small subset of actual NWCA 2011 data, with UIDs 
#' recoded.
#' @examples
#'   data(TreesEx)
#'   str(TreesEx)
#' @keywords datasets
"TreesEx"

#' Example Vascular Plant Cover data
#' 
#' A dataset containing vascular plant data for use in vascPlantMets() example.
#' 
#' @name VascPlantEx
#' @format A data frame with 1200 observations on the following 6 variables.
#'   \describe{
#'   \item{UID}{A unique site visit ID}
#'   \item{STATE}{Two letter state code of sample site}
#'   \item{USAC_REGION}{US Army Corps region of sample site}
#'   \item{PLOT}{Plot number of measurement}
#'   \item{USDA_NAME}{Taxon name corresponding to input taxa tables}
#'   \item{COVER}{A numeric vector representing the percent cover of a 
#'   taxon within a plot}
#'   }
#' @note This dataset is a small subset of actual NWCA 2011 data, with UIDs 
#' recoded.
#' @examples
#' head(VascPlantEx)
#' str(VascPlantEx) 
#' @keywords datasets
"VascPlantEx"

#' Example VMMI metric data
#' 
#' A dataset containing metrics used in VMMI calculation for use in example.
#' 
#' @name vmmiMetEx
#' @format A data frame with 10 observations on the following 8 variables.
#'  \describe{
#'  \item{UID}{A unique site visit ID}
#'  \item{NWCA_ECO4}{Aggregated ecoregions as used in NWCA 2011}
#'  \item{NWCA_WET_GRP}{Wetland type group as used in NWCA 2011}
#'  \item{ECO_X_WETGRP}{Combination of NWCA_ECO4 and NWCA_WET_GRP as used in 
#'  NWCA 2011}
#'  \item{FQAI_ALL}{Floristic Quality Assessment Index based on all species}
#'  \item{N_TOL}{Number of tolerant species, based on C-Value (NWCA_CC)}
#'  \item{RIMP_NATSPP}{Relative importance of native species}
#'  \item{XRCOV_MONOCOTS_NAT}{Relative mean cover of native monocots}
#'  }
#' @note This dataset is a small subset of actual NWCA 2011 data, with UIDs 
#' recoded.
#' @examples
#' head(vmmiMetEx)
#' str(vmmiMetEx)
#' @keywords datasets
"vmmiMetEx"

#' Example Vegetation Type and Ground Surface Attributes data
#' 
#' A dataset containing vegetation type and ground surface attributes data 
#' for use in vtype_gcovMets() example.
#' 
#' @name Vtype_GrCovEx
#' @format A data frame with 1775 observations on the following  variables.
#'  \describe{
#'  \item{UID}{A unique site visit ID}
#'  \item{PLOT}{Plot number of measurement}
#'  \item{PARAMETER}{Name of parameter measured}
#'  \item{RESULT}{Measured value of parameter}
#'  }
#'  
#' @note This dataset is a small subset of actual NWCA 2011 data, with UIDs 
#'  recoded.
#'  
#' @examples
#' head(Vtype_GrCovEx)
#' str(Vtype_GrCovEx)
#' @keywords datasets
"Vtype_GrCovEx"




