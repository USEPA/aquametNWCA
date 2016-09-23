\name{VascPlantEx}
\alias{VascPlantEx}
\concept{Example Vascular Plant Data}
\docType{data}
\title{Example Vascular Plant Cover data
}
\description{
A dataset containing vascular plant data for use in vascPlantMets() example.
}
\usage{data(VascPlantEx)}
\format{
  A data frame with 1200 observations on the following 6 variables.
  \describe{
    \item{UID}{A unique site visit ID}
    \item{STATE}{Two letter state code of sample site}
    \item{USAC_REGION}{US Army Corps region of sample site}
    \item{PLOT}{Plot number of measurement}
    \item{USDA_NAME}{Taxon name corresponding to input taxa tables}
    \item{COVER}{A numeric vector representing the percent cover of a taxon within a plot}
  }
}
\note{
This dataset is a small subset of actual NWCA 2011 data, with UIDs recoded.
}
\examples{
data(VascPlantEx)
str(VascPlantEx) 
head(VascPlantEx)
}
\keyword{datasets}
