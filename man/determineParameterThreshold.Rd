% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7_Clustering.R
\name{determineParameterThreshold}
\alias{determineParameterThreshold}
\title{Plot density of a given parameter on an interactive plot}
\usage{
determineParameterThreshold(
  data = NULL,
  parameter = NULL,
  displayedCells = 1e+05
)
}
\arguments{
\item{data}{A string defining which dataset to use. It should be the output of the \code{openUMAPData()} function. Defaults to \code{NULL}.}

\item{parameter}{A string defining which parameter to plot Defaults to \code{NULL}.}

\item{displayedCells}{A string defining the maximum number of cells to be plotted. Defaults to \code{100000}.}
}
\value{
A \code{ggplotly} plot. Please note that after the threshold selection is performed, the graph should be close with \code{dev.off()} call.
}
\description{
This function allows to display the density of a desired parameter on an interactive \code{ggplotly} plot to choose the best threshold to use for negative and positive cells separation.
}
