% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0_Preprocessing.R
\name{exportRDSFilesFromPool}
\alias{exportRDSFilesFromPool}
\title{Export new \code{rds} files from a previously generated pooled version}
\usage{
exportRDSFilesFromPool(RDSFileToUse = NULL, coresNumber = 2)
}
\arguments{
\item{RDSFileToUse}{A vector of all parameters to plot. Defaults to \code{NULL}.}

\item{coresNumber}{An integer defining the number of cores to use to export the files. It corresponds to the number of samples to be treated concomitantly. Defaults to \code{2}.}
}
\value{
Generated \code{rds} files are saved to \code{rds} directory. Please note that they will overwrite any preexisting \code{rds} file with the same name.
}
\description{
This function exports new \code{rds} files (one per sample) from a previously generated pooled version
}
