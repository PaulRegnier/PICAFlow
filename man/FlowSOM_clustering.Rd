% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7_Clustering.R
\name{FlowSOM_clustering}
\alias{FlowSOM_clustering}
\title{Perform FlowSOM clustering}
\usage{
FlowSOM_clustering(
  data = NULL,
  parametersToUse = NULL,
  seed = NULL,
  maxMeta = 90,
  datasetFolder = "full"
)
}
\arguments{
\item{data}{A data frame containing the dataset to use. Defaults to \code{NULL}.}

\item{parametersToUse}{A vector defining the parameters to use for the clustering. Defaults to \code{NULL}.}

\item{seed}{An integer defining the seed to use to obtain reproducible results. Defaults to \code{NULL}.}

\item{maxMeta}{An integer defining the maximum number of clusters to generate. Please note that FlowSOM algorithm cannot exceed 90 clusters maximum. However, FlowSOM will surely output way less than this number, depending on the dataset of interest. Defaults to \code{90}.}

\item{datasetFolder}{A string defining in which folder to save the results. It can be either \code{full} for the full downsampled dataset, \code{training} for the training downsampled subdataset or \code{validation} for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to \code{full}.}
}
\value{
Generated rds file is saved to \code{rds} directory.
}
\description{
This function allows to perform FlowSOM clustering algorithm on the dataset and parameters of interest.
}
