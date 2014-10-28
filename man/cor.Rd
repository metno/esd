\name{corfield}
\alias{corfield}
\alias{corfield.default}
\alias{corfield.zoo}
\alias{corfield.field}
\alias{corfield.field.station}
\alias{corfield.station}
\alias{corfield.eof}
\title{Correlation}
\description{
Compute the correlation between field objects and station/field.
}
\usage{
corfield(x,y, use='pairwise.complete.obs',...)
corfield.default(x,y,...)
corfield.zoo(x,y,plot=TRUE,use='pairwise.complete.obs',...) 
corfield.field(x,y,plot=TRUE,use='pairwise.complete.obs',...)
corfield.field.station(x,y,plot=TRUE, use='pairwise.complete.obs',...)
corfield.station(x,y,plot=TRUE, use='pairwise.complete.obs',...)
corfield.eof(x,y,pattern=1,plot=TRUE,
             use='pairwise.complete.obs',na.action='na.omit',...)
}
\arguments{
    \item{x}{data object}
    \item{y}{data object}
    \item{plot}{TRUE: plot the results}
    \item{use}{see \code{\link{cor}}}.
 }
\value{Map of correlation}
\author{R.E. Benestad and A. Mezghani}
\examples{
x <- t2m.ERAINT(lon=c(-40,30),lat=c(0,50))
y <- t2m.NCEP(lon=c(-40,30),lat=c(0,50))
r <- corfield(annual(x),annual(y))

data(Oslo)
t2m <- t2m.ERAINT()
x <- subset(Oslo,it=1)
y <- subset(t2m,it=1)
r <- corfield(x,y)
}
\keyword{manip}
