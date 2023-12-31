% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7_Clustering.R
\name{clustersPercentagesHeatmap}
\alias{clustersPercentagesHeatmap}
\title{Export information about clusters abundance (group- or sample-wise)}
\usage{
clustersPercentagesHeatmap(
  prefix = NULL,
  metricUsed = "median",
  datasetFolder = "full",
  mode = "group"
)
}
\arguments{
\item{prefix}{A string defining the prefix of the clusters abundance file to use. Defaults to \code{NULL}.}

\item{metricUsed}{A string defining the metric to use for statistics computation. It can be either \code{mean} or \code{median}. Defaults to \code{median}.}

\item{datasetFolder}{A string defining in which folder to save the results. It can be either \code{full} for the full downsampled dataset, \code{training} for the training downsampled subdataset or \code{validation} for the validation downsampled subdataset. The value must match the origin of the data used. Defaults to \code{full}.}

\item{mode}{A string defining the mode to use for heatmap generation. It can be either \code{group} or \code{sample}. Defaults to \code{group}.}
}
\value{
Generated PDF files are saved to \verb{output > 7_Clustering > datasetFolder}.
}
\description{
This function allows to export heatmaps about clusters abundance per group or per sample.
}
