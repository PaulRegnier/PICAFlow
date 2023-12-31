% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_Analysis.R
\name{openFiles}
\alias{openFiles}
\title{Open clusters abundance files}
\usage{
openFiles(pattern = NULL, sampleNamesColumn = NULL)
}
\arguments{
\item{pattern}{A string defining the pattern used to determine which files should be opened. Defaults to \code{NULL}.}

\item{sampleNamesColumn}{A string defining the column name the be used as common comparator for each file to merge. This column should be present in each file to be opened and should present identical values (for instance, individualized sample names). Defaults to \code{NULL}.}
}
\value{
Return value is directly passed to \code{mergeData()} function.
}
\description{
This function allows to open one or several clusters abundance files. As this function is called by \code{mergeData()} function, it should not be called as is.
}
