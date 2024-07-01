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
#' str(ccNatNWCA)
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
#' @examples
#' head(taxaNWCA)
#' str(taxaNWCA)
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
#' data(TreesEx)
#' str(TreesEx)
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


#' NWCA 2016 Coefficient of Conservation Values (C-values)
#'
#' A data frame containing NWCA 2016 state- or region-specific
#' Coefficient of Conservatism values for vascular plant taxa
#'
#' @name cvalNWCA_2016
#' @format A data frame containing 24206 observations on the
#' following variables
#' \describe{
#'  \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA
#'  taxonomic ID number}
#'  \item{NWCA_NAME}{Accepted USDA PLANTS name, or the official NWCA name
#'  where not available}
#'  \item{GEOG_TYPE}{Type of Geographic Unit}
#'  \item{GEOG_ID}{Geographic unit to which value applies, which is NWC_CREG16
#'  in this dataset}
#'  \item{NWCA_CVAL}{C-value region-specific Coefficient of
#'  Conservatism value as used in NWCA 2016}
#' }
#' @note This dataset is the list of C-values used for site-specific and
#' taxon-specific assignments for NWCA 2016.
#' @keywords datasets
#' @examples
#' head(cvalNWCA_2016)
#' str(cvalNWCA_2016)
#' @keywords datasets
"cvalNWCA_2016"

#' NWCA 2016 Native Status Values
#'
#' A data frame containing NWCA 2016 state-specific
#' native status values for vascular plant taxa
#'
#' @name nativeNWCA_2016
#' @format A data frame containing 21356 observations on the
#' following variables
#' \describe{
#'  \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA
#'  taxonomic ID number}
#'  \item{NWCA_NAME}{Accepted USDA PLANTS name, or the official NWCA name
#'  where not available}
#'  \item{GEOG_ID}{Geographic unit to which value applies, which is
#'  the two-letter state abbreviation in this dataset}
#'  \item{NWCA_NATSTAT}{Native status of taxon. Valid values:
#'  ADV (Adventive), CRYP (Cryptogenic), INTR (introduced),
#'  NAT (native), UND (Undetermined).}
#'  \item{NATSTAT_ALT}{Combined native status values as used in calculations.
#'  Valid values: ALIEN (INTR + ADV), CRYP, NAT, UND}
#'  \item{ALIEN}{Indicator value (1/0) for whether taxon is considered alien.}
#'  \item{AC}{Indicator value (1/0) for whether taxon is considered alien or
#'  cryptogenic}
#' }
#' @note This dataset is the list of native status values used
#' for site-specific and #' taxon-specific assignments for
#' NWCA 2016.
#' @keywords datasets
#' @examples
#' head(nativeNWCA_2016)
#' str(nativeNWCA_2016)
#' @keywords datasets
"nativeNWCA_2016"

#' NWCA 2016 Vascular Plant Taxa List
#'
#' A dataset containing taxonomy and basic characteristics for vascular plant
#' taxa as used in NWCA 2016.
#'
#' @name taxaNWCA_2016
#'
#' @format A data frame with 5044 observations on the following 28 variables.
#' \describe{
#'   \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA taxonomic
#'   ID number}
#'   \item{NWCA_NAME}{Accepted USDA PLANTS name, or the official NWCA name
#'     where not available}
#'   \item{ACCEPTED_SYMBOL}{Accepted symbol for taxon from USDA
#'   PLANTS or NWCA 2016 if not in USDA PLANTS}
#'   \item{TAXON_LEVEL}{Taxonomic level of name: CL (class), F (Family),
#'   G (Genus), GH (Growth habit), H (), NV (Nonvascular),
#'   SC (), SF (),
#'   SSP (Subspecies), SP (Species), T (), U (),
#'   VAR (Variety)}
#'   \item{DIVISION}{Division level of taxonomy}
#'   \item{CLASS}{Class level of taxonomy}
#'   \item{SUBCLASS}{Subclass level of taxonomy}
#'   \item{ORDER}{Order level of taxonomy}
#'   \item{FAMILY}{Family level of taxonomy}
#'   \item{GENUS}{Genus level of taxonomy}
#'   \item{CATEGORY}{Plant category of taxon, including MONOCOT, GYMNOSPERM,
#'   DICOT, FERN, LICHEN, HORSETAIL, MOSS, LYCOPOD, and LIVERWORT}
#'   \item{GROWTH_HABIT}{Growth habit of taxon as provided by USDA PLANTS}
#'   \item{DURATION}{Lifecycle duration of taxon as provided by USDA PLANTS}
#'   \item{VINE_ALL}{Indicator value for all vine types (1/0)}
#'   \item{HERB}{Indicator value for herbaceous vascular plant (1/0)}
#'   \item{SHRUB_COMB}{Indicator value for all types of shrubs (1/0)}
#'   \item{TREE_COMB}{Indicator value for all types of trees (1/0)}
#'   }
#'   @note This dataset is the taxa list used for state-specific CC and
#'   native status values used in NWCA 2016.
#'
#' @examples
#' head(taxaNWCA_2016)
#' str(taxaNWCA_2016)
#'
#'   @keywords datasets
"taxaNWCA_2016"

#' NWCA 2016 Wetland Indicator Status
#'
#' A dataset containing NWCA 2016 USAC-specific Wetland Indicator Status
#' values for vascular plant taxa
#'
#' @name wisNWCA_2016
#' @format A data frame with 9584 observations on the following 5 variables.
#' \describe{
#'     \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA taxonomic
#'     ID number}
#'     \item{NWCA_NAME}{Accepted USDA PLANTS name, or the official NWCA name
#'     where not available}
#'     \item{GEOG_ID}{Geographic unit to which value applies, which is US Army
#'     Corps region in this dataset}
#'     \item{WIS}{USAC-specific Wetland Indicator Status as used for NWCA 2016}
#'     \item{ECOIND1}{Wetland indicator status converted to numeric value, with
#'     lower numbers more characteristic of wetlands}
#'     \item{ECOIND2}{Wetland indicator status converted to numeric value, with
#'     higher numbers more characteristic of wetlands}
#' }
#' @note This dataset is the taxa list used for USAC-specific Wetland Indicator
#' Status values used in NWCA 2016.
#'
#' @examples
#' head(wisNWCA_2016)
#' str(wisNWCA_2016)
#'
#' @keywords datasets
"wisNWCA_2016"
