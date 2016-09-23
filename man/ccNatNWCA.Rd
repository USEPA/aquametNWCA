\name{ccNatNWCA}
\alias{ccNatNWCA}
\alias{CC Values}
\concept{Native Status Values}
\docType{data}
\title{
NWCA 2011 Plant CC and Native Status 
}
\description{
A dataset containing NWCA 2011 state-specific Coefficient of Conservatism and native status values for vascular plant taxa
}
\usage{data("ccNatNWCA")}
\format{
  A data frame with 13203 observations on the following 5 variables.
  \describe{
    \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA taxonomic ID number}
    \item{USDA_NAME}{Accepted USDA PLANTS name, or the official NWCA name where not available}
    \item{GEOG_ID}{Geographic unit to which value applies, which is STATE in this dataset}
    \item{NWCA_CC}{State-specific Coefficient of Conservatism value as determined for NWCA 2011}
    \item{NWCA_NATSTAT}{State-specific Native status as determined for NWCA 2011}
  }
}
\note{
This dataset is the taxa list used for state-specific CC and native status values used in NWCA 2011.
}
\examples{
data(ccNatNWCA)
str(ccNatNWCA)
}
\keyword{datasets}
