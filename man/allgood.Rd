% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgood.R
\name{allgood}
\alias{allgood}
\title{A small function that removes stations with missing values from a group}
\usage{
allgood(x, miss = 0.1, verbose = FALSE)
}
\arguments{
\item{x}{a \code{station} object}

\item{miss}{fraction of data that may be missing, e.g., if miss=.1 then stations with than 10\% missing data are removed}

\item{verbose}{a boolean; if TRUE print information about progress}
}
\description{
Useful before performing PCA.
}
