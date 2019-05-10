\name{is.inside}
\alias{is.inside}
\alias{test.is.inside}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{is.inside checks whether a point or a set of points is inside a polygon}
\description{
\code{is.inside} returns \code{TRUE} if the point(s) \code{x} is/are inside a polygon \code{y}. 
}
\usage{
is.inside(x, y, verbose = FALSE, plot = FALSE)
test.is.inside(N=5000,verbose=TRUE,plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an esd-object or a list/data.frame with the elements \code{x$x} and \code{x$y} containing the coordinates.}
  \item{y}{A polygon in the shape of a list/data.frame with the elements \code{x$x} and \code{x$y} containing the coordinates.}
  \item{verbose}{
\code{TRUE} prints out diagnostics for the code.
}
  \item{plot}{
\code{TRUE} provides a graphical disgnostic.
}
\item{N}{Number of tests with random coordinates}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
