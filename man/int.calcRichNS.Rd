% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VascPlantMetric_functions.r
\name{int.calcRichNS}
\alias{int.calcRichNS}
\title{Calculate richness metrics for native status subsets}
\usage{
int.calcRichNS(x, y, natstat, grpname, sampID)
}
\arguments{
\item{x}{Data frame containing cover data summarized by sampID
variables, TAXON, NWCA_NATSTAT, and DISTINCT}

\item{y}{Data frame containing cover data summarized by sampID
variables, PLOT, TAXON, NWCA_NATSTAT, and DISTINCT}

\item{natstat}{Character vector containing Values of NWCA_NATSTAT
variable to include in rich metrics}

\item{grpname}{String containing suffix to add to metric name to
represent this group}

\item{sampID}{A character vector containing the name(s) of
variable(s) necessary to identify unique samples}
}
\value{
Data frame containing \emph{sampID} variables, PARAMETER, and RESULT,
  with one row of results per parameter and sample. The values for PARAMETER
  consist of the metric name concatenated with grpname value (represented as
  GRP below):

TOTN_GRP: Number of unique taxa in sample

XN_GRP: Mean number of taxa per plot

MEDN_GRP: Median number of taxa per plot

SDN_GRP: Standard deviation of number of taxa per plot
}
\description{
This internal function calculates species richness
of sample for subset based on native status values specified.
Used by \code{calcRichness()}
}
\section{Warning}{
 This function not intended for use on its own
}

\author{
Karen Blocksom
}
