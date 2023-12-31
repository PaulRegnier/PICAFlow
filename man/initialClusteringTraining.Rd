% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7_Clustering.R
\name{initialClusteringTraining}
\alias{initialClusteringTraining}
\title{Perform initial cell clustering training}
\usage{
initialClusteringTraining(
  data = NULL,
  parametersToAnalyze = NULL,
  subsetDownsampled = 25000,
  clusterMinPercentage = 0.5
)
}
\arguments{
\item{data}{A data frame containing the dataset to use. Defaults to \code{NULL}.}

\item{parametersToAnalyze}{A character vector defining the parameters to analyze in the subsequent dataset. Defaults to \code{NULL}.}

\item{subsetDownsampled}{An integer defining the maximum number of cells to use for the analysis. Defaults to \code{25000}.}

\item{clusterMinPercentage}{An integer defining the minimum frequency of a cluster to be kept as is. Values represent percentages and range from \code{0} to \code{100}. Defaults to \code{0.5}.}
}
\value{
A list of 4 elements: \code{dataTraining} which contains the actual data used for training, \code{clustering} which contains the clusters associated to each cell, \code{treeCuts} which contains the number of clusters for different hierarchical tree cutting and \code{hSequence} which contains the different thresholds used to cut the hierarchical tree.
}
\description{
This function allows perform the initial step of the cell clustering training.
}
