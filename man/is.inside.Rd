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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{ 
## Example 
library(esd)
library(readINAMdata)
## Get a polygon for Mozambique
data("Moz")

osmoz <- !is.inside(data.frame(x=lons,y=lats),Moz)
plot(lons,lats)
points(lons[osmoz],lats[osmoz],pch=19,col='red')
points(lons[!osmoz],lats[!osmoz],pch=19,col='green')
points(lons[is.na(osmoz)],lats[is.na(osmoz)],pch=19,col='black')
lines(Moz,type='b')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
