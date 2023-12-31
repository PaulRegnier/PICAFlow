% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_Downsampling.R
\name{poolData}
\alias{poolData}
\title{Reduce and downsample dataset}
\usage{
poolData(
  flowSet = NULL,
  groupVector = NULL,
  parametersToKeep = NULL,
  downsampleMinEvents = 10000,
  rescale = TRUE,
  rescale_min = 1,
  rescale_max = 10,
  maxCellsNb = 5e+05,
  estimateThreshold = FALSE,
  coresNumber = 1
)
}
\arguments{
\item{flowSet}{A string defining the \code{flowSet} to use. Defaults to \code{NULL}.}

\item{groupVector}{A character vector defining the groups for each sample. Defaults to \code{NULL}.}

\item{parametersToKeep}{A character vector defining the parameters to keep in the subsequent dataset. Typically, one wants to remove parameters that were already used to gate cells. Defaults to \code{NULL}.}

\item{downsampleMinEvents}{An integer defining the minimum number of cells per sample to be included in the downsampled dataset. Defaults to \code{10000}.}

\item{rescale}{A boolean defining if the data should be rescaled parameter-by-parameter. Defaults to \code{TRUE}.}

\item{rescale_min}{An integer defining the minimum value for the rescaled dataset. Used only if \code{rescale} is set to \code{TRUE}. Defaults to \code{1}.}

\item{rescale_max}{An integer defining the maximum value for the rescaled dataset. Used only if \code{rescale} is set to \code{TRUE}. Defaults to \code{10}.}

\item{maxCellsNb}{An integer defining the maximum number of cells for the downsampled dataset. If the actual number of cells is greater than \code{maxCellsNb}, a global downsampling factor will be applied to keep the number of cells of the downsampled dataset below \code{maxCellsNb}. Defaults to \code{500000}.}

\item{estimateThreshold}{A boolean defining if a given \code{maxCellsNb} threshold should rather be replaced by the test of all thresholds possible according to the number of cells in each sample (if set to \code{TRUE}). Be careful, this function is expensive in terms of computational resources, so do not use it if your dataset is rather large. Defaults to \code{FALSE}.}

\item{coresNumber}{An integer defining the number of cores to use to analyze the peaks. Be careful, this function is expensive in terms of computational resources, so do not increase this value too much (typically 1 to 4 maximum). Defaults to \code{1}.}
}
\value{
If \code{estimateThreshold = FALSE}, the function will return a list of 2 elements named \code{data} (containing the actual downsampled data) and \code{message} (containing the log messages generated during the downsampling process). If \code{estimateThreshold = TRUE}, the function will rather return a plot in PDF format showing the relation between a given threshold and the final number of cells after downsampling.
}
\description{
This function allows to reduce the dataset by removing useless paramaters as well as tag cells for the construction of a downsampled dataset in which all groups contribute equally to the dataset and in which each sample contributes equally within a given group. Please note that from this step, \code{rds} files will contain table-like data and not \code{flowFrame} nor \code{flowSet} objects anymore.
}
