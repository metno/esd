\name{aggregate}
\alias{aggregate.comb}
\alias{aggregate.field}
\alias{aggregate.station}
\alias{aggregate.area}
\title{aggregate}
\description{
  The aggregation is based on the S3 method for \code{zoo} objects, but
  takes care of extra house keeping, such as attributes with meta data.
}
\usage{
aggregate.station(x,by,FUN = 'mean', ..., regular = NULL, frequency = NULL) 
aggregate.field(x,by,FUN = 'mean', ..., regular = NULL, frequency = NULL)
aggregate.comb(x,by,FUN = 'mean', ...,
               regular = NULL, frequency = NULL)
aggregate.area(x,is=NULL,it=NULL,FUN='sum',na.rm=TRUE,smallx=FALSE)
}
\arguments{
\item{x}{A \code{\link{station}}, \code{\link{spell}} or a \code{\link{field}} object}
\item{by}{see \code{\link{aggregate.zoo}}}
\item{FUN}{see \code{\link{aggregate.zoo}}}
\item{regular}{see \code{\link{aggregate.zoo}}}
\item{frequency}{see \code{\link{aggregate.zoo}}}
\item{is}{spatial selection - see \code{\link{subset.field}}}
\item{na.rm}{TRUE: ignore NA - see see \code{\link{mean}}}
}
\value{The call returns a station object}

\seealso{
 \code{\link{spatial.avg.field}} \code{\link{as.4seasons}}, \code{\link{annual}}
}

\author{R.E. Benestad}
\examples{
# Example: use aggregate to compute annual mean temperature for Svalbard:
data(Svalbard)
y <- aggregate(Svalbard, year, FUN='mean', na.rm = FALSE) 
plot(y)

# Example for getting seasonal aggregate of standard deviation of
data(Oslo)
ym <- as.4seasons(Oslo,FUN='mean')
ys <- as.4seasons(Oslo,FUN='sd')
y <- combine(ym,ys)
plot(y)

x <- t2m.NCEP()
z <- aggregate.area(x,FUN='mean')
plot(z)
}
\keyword{utilities}
