% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7_Clustering.R
\name{plotInitialClusteringTraining}
\alias{plotInitialClusteringTraining}
\title{Display an interactive plot for choosing an optimal cutoff value}
\usage{
plotInitialClusteringTraining(initialClustering = NULL)
}
\arguments{
\item{initialClustering}{An object generated using the \code{initialClusteringTraining()} function. Defaults to \code{NULL}.}
}
\value{
A \code{ggplotly} plot. Please note that after the cutoff selection is performed, the graph should be close with \code{dev.off()} call.
}
\description{
This function allows to display an interactive \code{ggplotly} decision plot to help the user to choose the optimal cutoff for the hierarchical tree cutting.
}
