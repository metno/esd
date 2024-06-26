% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.mvr}
\alias{plot.mvr}
\title{Plot esd objects}
\usage{
\method{plot}{mvr}(x, verbose = FALSE, ...)
}
\arguments{
\item{x}{the object to be plotted}

\item{verbose}{a boolean; if TRUE print information about progress}

\item{\dots}{additional arguments}
}
\value{
None
}
\description{
The plot functions are S3 methods for esd objects, based on \code{plot}. 
The function \code{plot.mvr} produces a plot for an \code{mvr} object which
is the output of the multi-variate regression function \code{MVR},
}
\seealso{
plot.station plot.eof plot.ds
}
\keyword{graphics}
\keyword{plot}
