\name{gridmap}
\alias{gridmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates a griddded map
%%  ~~function to do ... ~~
}
\description{A function that uses \code{LatticeKrieg} and elevation data to grid station based data and present a map.}
\usage{
gridmap(Y, FUN = "mean", colbar = list(col='t2m'), project = "lonlat", xlim = NULL, ylim = NULL, zlim = NULL, verbose = FALSE, plot = FALSE, new=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A station object}
  \item{FUN}{A function or name of a function, e.g, "mean" or "trend"}
  \item{colbar}{A list specifying the color bar, e.g., list(col="precip", breaks=seq(1,10), rev=FALSE)}
  \item{project}{projection: "lonlat" or "sphere"}
  \item{xlim}{range of x-axis}
  \item{ylim}{range of y-axis}
  \item{zlim}{range of color axis}
  \item{verbose}{if TRUE print information about progress}
  \item{plot}{if TRUE display results as plots}
}

\examples{
data("precip.NORDKLIM")
precip.gp <- gridmap(precip.NORDKLIM, plot=TRUE)
}
