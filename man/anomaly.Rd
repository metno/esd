% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anomaly.R
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
\alias{anomaly.dsensemble}
\alias{as.anomaly}
\alias{as.anomaly.default}
\alias{as.anomaly.zoo}
\alias{as.anomaly.list}
\alias{as.anomaly.station}
\alias{as.anomaly.field}
\alias{climatology}
\alias{as.climatology}
\title{Anomaly and Climatology}
\usage{
anomaly(x, ...)
}
\arguments{
\item{x}{A station or field object}

\item{\dots}{additional arguments}

\item{ref}{vector defining the reference interval}

\item{na.rm}{a boolean; if TRUE remove NA values}

\item{verbose}{a boolean; if TRUE print information about progress}
}
\value{
a similar object as x containing anomalies and climatology
}
\description{
S3-method that computes anomalies and/or climatology for time series and
fields.
}
\details{
In 'anomaly.dsensemble', the default value of the reference period is taken
as the available time period from observations, i.e., same time period as in
attribute `station' is used as base period to compute anomalies of GCM
downscaled results.
}
\examples{
data(ferder)
plot(anomaly(ferder))


}
\seealso{
as.stand
}
\keyword{utilities}
