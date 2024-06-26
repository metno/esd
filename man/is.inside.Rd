% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.inside.R
\name{is.inside}
\alias{is.inside}
\alias{test.is.inside}
\title{A test to see if a point is inside a polygon}
\usage{
is.inside(x, y, verbose = FALSE, plot = FALSE)
}
\arguments{
\item{x}{an esd-object or a list/data.frame with the elements \code{x$x} and
\code{x$y} containing the coordinates.}

\item{y}{A polygon in the shape of a list/data.frame with the elements
\code{x$x} and \code{x$y} containing the coordinates.}

\item{verbose}{\code{TRUE} prints out diagnostics for the code.}

\item{plot}{\code{TRUE} provides a graphical disgnostic.}

\item{N}{Number of tests with random coordinates}
}
\value{
a boolean; \code{TRUE} if the point(s) \code{x} is/are inside the polygon \code{y}.
}
\description{
\code{is.inside} checks whether a point or a set of points is inside a polygon, e.g., borders of a country.
}
\examples{

\dontrun{ 
library(readINAMdata)
data("Moz")

osmoz <- !is.inside(data.frame(x=lons,y=lats),Moz)
plot(lons,lats)
points(lons[osmoz],lats[osmoz],pch=19,col='red')
points(lons[!osmoz],lats[!osmoz],pch=19,col='green')
points(lons[is.na(osmoz)],lats[is.na(osmoz)],pch=19,col='black')
lines(Moz,type='b')
}

}
