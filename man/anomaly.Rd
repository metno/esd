\name{anomaly}
\alias{anomaly}
\alias{anomaly.default}
\alias{anomaly.comb}
\alias{anomaly.field}
\alias{anomaly.station}
\alias{anomaly.annual}
\alias{anomaly.month}
\alias{anomaly.season}
\alias{anomaly.day}
\alias{climatology}
\alias{climatology.default}
\alias{climatology.field}
\alias{climatology.station}
\alias{clim2pca}
\alias{clim2pca.default}
\alias{clim2pca.month}
\alias{clim2pca.day}
\title{Anomaly and Climatology}
\description{
  S3-method that computes anomalies and/or climatology for time series and fields.
  \code{clim2pca} is unfinished
}
\usage{
anomaly(x,...)
anomaly.default(x,...)
anomaly.comb(x,...)
anomaly.field(x,...)
anomaly.station(x,...)
anomaly.dsensemble(x,ref=NULL,...)
anomaly.annual(x,ref=1961:1990)
anomaly.month(x,ref=NULL)
anomaly.season(x,ref=NULL) 
anomaly.day(x,ref=NULL)
climatology(x,...)
climatology.default(x)
climatology.field(x)
climatology.station(x)
clim2pca(x,...)
clim2pca.default(x)
clim2pca.month(x)
clim2pca.day(x)
}
\arguments{
  \item{x}{A station or field object}
  \item{ref}{vector defining the reference interval}
}
\value{The call returns a similar object as x containing anomalies and
  data climatology.

In `anomaly.dsensemble', the default value of the reference period is
taken as the available time period from observations, i.e. same time
period as in attribute `station' is used as baseperiod to compute anomalies of GCM downscaled results.}

\seealso{
 \code{\link{as.anomaly}}, \code{\link{as.climatology}}
}

\author{R.E. Benestad}
\examples{
data(ferder)
plot(anomaly(ferder))

}
\keyword{utilities}
