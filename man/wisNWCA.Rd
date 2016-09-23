\name{wisNWCA}
\alias{wisNWCA}
\concept{WIS Values}
\concept{Wetland Indicator Status Values}
\docType{data}
\title{
NWCA 2011 Wetland Indicator Status 
}
\description{
A dataset containing NWCA 2011 USAC-specific Wetland Indicator Status values for vascular plant taxa
}
\usage{data(wisNWCA)}
\format{
  A data frame with 19818 observations on the following 5 variables.
  \describe{
    \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA taxonomic ID number}
    \item{USDA_NAME}{Accepted USDA PLANTS name, or the official NWCA name where not available}
    \item{GEOG_ID}{Geographic unit to which value applies, which is US Army Corps region in this dataset}
    \item{WIS}{USAC-specific Wetland Indicator Status as used for NWCA 2011}
  }
}
\note{
This dataset is the taxa list used for USAC-specific Wetland Indicator Status values used in NWCA 2011.
}

\examples{
data(wisNWCA)
str(wisNWCA)
}
\keyword{datasets}
