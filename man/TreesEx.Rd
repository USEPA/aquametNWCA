\name{TreesEx}
\alias{TreesEx}
\concept{Example Tree Data}
\docType{data}
\title{
Example Tree data
}
\description{
A dataset containing tree data for use in treeMets() example.
}
\usage{data(TreesEx)}
\format{
  A data frame with 1528 observations on the following 6 variables.
  \describe{
    \item{UID}{A unique site visit ID}
    \item{PLOT}{Plot number of measurement}
    \item{PAGE}{Page of field form for measurement}
    \item{LINE}{Line number of measurement}
    \item{PARAMETER}{Name of tree parameter measured}
    \item{RESULT}{Measured value of parameter}
  }
}
\note{
This dataset is a small subset of actual NWCA 2011 data, with UIDs recoded.
}
\examples{
data(TreesEx)
str(TreesEx)
}
\keyword{datasets}
