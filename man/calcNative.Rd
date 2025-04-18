% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcNative.r
\name{calcNative}
\alias{calcNative}
\title{Calculate metrics based only on native status}
\usage{
calcNative(vascIn, sampID = "UID")
}
\arguments{
\item{vascIn}{Data frame containing cover data summarized by
UID and TAXON, with the following fields:
\itemize{
    \item sampID: Variable(s) identified in \emph{sampID} argument

    \item TAXON: Taxon name

    \item XABCOV: Mean percent cover of taxon across plots

    \item TOTN: Number of taxa in sample

    \item sXRCOV: proportion of summed cover across all taxa
    (XTOTABCOV) represented by taxon in sample

    \item sRFREQ: Relative frequency of a taxon, calculated
    as the percentage of the total frequency of taxon
    occurrence across all taxa for a UID

    \item NWCA_NATSTAT: Native status variable with
    categories of 'NAT', 'ADV', 'CRYP', 'INTR', 'UND'
 }}

\item{sampID}{A character vector containing the name(s) of
variable(s) necessary to identify unique samples, 'UID' by default}
}
\value{
Data frame containing \emph{sampID} variables, PARAMETER, RESULT,
  where values of PARAMETER consist of the metric name concatenated with
  trait value (represented as TRAITNM below):
\itemize{
  \item PCTN_TRAITNM: Number of taxa with trait as percentage of \emph{TOTN}

  \item XABCOV_TRAITNM: Sum of \emph{XABCOV} values across taxa with trait

  \item XRCOV_TRAITNM: Sum of \emph{sXRCOV} values across taxa with trait

  \item RFREQ_TRAITNM: Sum of \emph{sRFREQ} values across taxa with trait value

  \item RIMP_TRAITNM: Relative importance ((RFREQ_TRAITVAL + XRCOV_TRAITVAL)/2)
  of taxa with trait value
  }
A list of metric descriptions is provided in the document named 
\href{https://github.com/USEPA/aquametNWCA/blob/main/inst/VascPlant_Metric_Descriptions.pdf}{VascPlant_Metric_Descriptions.pdf}
}
\description{
This function calculates all metrics based
only on native status.
}
\examples{
head(VascPlantEx)
exPlant <- prepareData(VascPlantEx,
  taxon_name = "USDA_NAME",
  inTaxa = taxaNWCA, inNat = ccNatNWCA, inCVal = ccNatNWCA,
  inWIS = wisNWCA, cValReg = "STATE"
)

natEx <- calcNative(exPlant$byUIDspp)

head(natEx)
unique(natEx$PARAMETER)
}
\references{
US Environmental Protection Agency. 2016. National Wetland
Condition Assessment: 2011 Technical Report. EPA-843-R-15-006. US
Environmental Protection Agency, Washington, DC.
}
\author{
Karen Blocksom \email{Blocksom.karen@epa.gov}
}
