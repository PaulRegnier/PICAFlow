% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7_Clustering.R
\name{finalClusteringTraining}
\alias{finalClusteringTraining}
\title{Perform final cell clustering training}
\usage{
finalClusteringTraining(
  initialClusteringData = NULL,
  clusterMinPercentage = NULL,
  cutoff = NULL
)
}
\arguments{
\item{initialClusteringData}{An object generated using the \code{initialClusteringTraining()} function. Defaults to \code{NULL}.}

\item{clusterMinPercentage}{An integer defining the minimum frequency of a cluster to be kept as is. Values represent percentages and range from \code{0} to \code{100}. Defaults to \code{0.5}.}

\item{cutoff}{An integer defining the actual cutoff to use. Defaults to \code{NULL}.}
}
\value{
A list of 2 elements: \code{dataTraining} which contains the actual data used for training and \code{cutoff} which contains the cutoff value used for the tree cutting.
}
\description{
This function allows to apply the selected cutoff (for hierarchical clustering tree cutting) to the training dataset.
}
