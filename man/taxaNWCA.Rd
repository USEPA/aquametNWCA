\name{taxaNWCA}
\alias{taxaNWCA}
\concept{NWCA Taxa List}
\concept{Category Values}
\concept{Duration Values}
\concept{Growth Habit Values}
\docType{data}
\title{
NWCA 2011 Vascular Plant Taxa List
}
\description{
A dataset containing taxonomy and basic characteristics for vascular plant taxa as used in NWCA 2011.
}
\usage{data(taxaNWCA)}
\format{
  A data frame with 3814 observations on the following 12 variables.
  \describe{
    \item{SPECIES_NAME_ID}{A numeric vector containing the NWCA taxonomic ID number}
    \item{USDA_NAME}{Accepted USDA PLANTS name, or the official NWCA name where not available}
    \item{DIVISION}{Division level of taxonomy}
    \item{ORDER}{Order level of taxonomy}
    \item{FAMILY}{Family level of taxonomy}
    \item{GENUS}{Genus level of taxonomy}
    \item{SPECIES}{Species level of taxonomy}
    \item{SUBSPECIES}{Subspecies level of taxonomy}
    \item{VARIETY}{Variety level of taxonomy}
    \item{CATEGORY}{Plant category of taxon, including MONOCOT, GYMNOSPERM, DICOT, 
    FERN, LICHEN, HORSETAIL, MOSS, LYCOPOD, and LIVERWORT}
    \item{GROWTH_HABIT}{Growth habit of taxon}
    \item{DURATION}{Lifecycle duration of taxon}
  }
}
\note{
This dataset is the taxa list used for taxonomy and basic characteristics used in NWCA 2011.
}
\examples{
data(taxaNWCA)
str(taxaNWCA)
}
\keyword{datasets}
