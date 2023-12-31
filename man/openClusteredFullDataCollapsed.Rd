% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7_Clustering.R
\name{openClusteredFullDataCollapsed}
\alias{openClusteredFullDataCollapsed}
\title{Open the desired clustered dataset}
\usage{
openClusteredFullDataCollapsed(datasetToUse)
}
\arguments{
\item{datasetToUse}{A string defining which dataset to open. It can be either \code{full} for the full downsampled dataset, \code{training} for the training downsampled subdataset or \code{validation} for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to \code{full}.}
}
\value{
A data frame containing the selected dataset to use.
}
\description{
This function allows to open the desired clustered dataset for subsequent analysis.
}
