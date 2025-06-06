% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.dsensemble.pca}
\alias{plot.dsensemble.pca}
\title{plot dsensemble pca results}
\usage{
\method{plot}{dsensemble.pca}(
  x,
  ...,
  pts = FALSE,
  target.show = TRUE,
  map.show = TRUE,
  it = 0,
  ip = 1,
  envcol = rgb(1, 0, 0, 0.2),
  legend.show = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{x}{input object to be plotted}

\item{\dots}{additional arguments}

\item{pts}{a boolean; if TRUE plot points?}

\item{target.show}{a boolean; if TRUE show diagnostics as a target (see \code{\link{diagnose}})}

\item{map.show}{a boolean; if TRUE show map of stations}

\item{legend.show}{a boolean; if TRUE show legend}

\item{verbose}{a boolean; if TRUE print information about progress}
}
\description{
plot dsensemble pca results
}
