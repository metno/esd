\name{aggregate}
\alias{aggregate.comb}
\alias{aggregate.field}
\alias{aggregate.station}
\alias{aggregate.area}
\alias{aggregate.size}
\alias{aggregate.size.matrix}
\alias{aggregate.size.field}
\title{aggregate}
\description{
  The aggregation is based on the S3 method for \code{zoo} objects, but
  takes care of extra house keeping, such as attributes with meta data.

  \code{aggregate.area} is used for aggregating spatial statistics, such
  as the global mean or the global area of some
  phenomenon. \code{aggregate.size} is similar to
  \code{aggregate.area}, but returns the size statistics (square meters) for individual
  events (defined as gridboxes touching each other).
}
\usage{
aggregate.station(x,by,FUN = 'mean', ..., regular = NULL, frequency = NULL) 
aggregate.field(x,by,FUN = 'mean', ..., regular = NULL, frequency = NULL)
aggregate.comb(x,by,FUN = 'mean', ...,
               regular = NULL, frequency = NULL)
aggregate.area(x,is=NULL,it=NULL,FUN='sum',na.rm=TRUE,smallx=FALSE,
                           a== 6378, x0=NULL)
aggregate.size.matrix(x,x0,plot=FALSE,verbose=FALSE,a=6.378e06,...)
aggregate.size.field(x,x0,plot=FALSE,verbose=FALSE,...)
}
\arguments{
\item{x}{A \code{\link{station}}, \code{\link{spell}} or a \code{\link{field}} object}
\item{by}{see \code{\link{aggregate.zoo}}}
\item{FUN}{see \code{\link{aggregate.zoo}}.  In addition, 'area','exceedance','lessthan'.}
\item{regular}{see \code{\link{aggregate.zoo}}}
\item{frequency}{see \code{\link{aggregate.zoo}}}
\item{is}{spatial selection - see \code{\link{subset.field}}}
\item{na.rm}{TRUE: ignore NA - see see \code{\link{mean}}}
\item{a}{Earth's radius (km)}
\item{x0}{threshold for FUN=c('exceedance','lessthan')}
}
\value{The call returns a station object}

\seealso{
 \code{\link{spatial.avg.field}} \code{\link{as.4seasons}}, \code{\link{annual}}
}

\author{R.E. Benestad}
\examples{
## Example: use aggregate to compute annual mean temperature for Svalbard:
data(Svalbard)
y <- aggregate(Svalbard, year, FUN='mean', na.rm = FALSE) 
plot(y)

## Example for getting seasonal aggregate of standard deviation of
data(Oslo)
ym <- as.4seasons(Oslo,FUN='mean')
ys <- as.4seasons(Oslo,FUN='sd')
y <- combine(ym,ys)
plot(y)

## Aggregate the global mean temperature
x <- t2m.NCEP()
z <- aggregate.area(x,FUN='mean')
plot(z)

## Aggregate statistics for the size of high-pressure systems
slp <- slp.DNMI()
highs <- aggregate.size(slp,x0=1015)
plot(highs)

lows <-  aggregate.size(-slp,x0=-1015)
plot(lows)
}
\keyword{utilities}
