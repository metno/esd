\name{transforms}
\alias{C.C.eq}
\alias{precip.vul}
\alias{t2m.vul}
\alias{precip.rv}
\alias{nv}
\alias{precip.Pr}
\alias{t2m.Pr}
\alias{NE}

\title{Various formulas, equations and transforms.}
\description{Computes different formulas.

\code{C.C.eq}: Clapeyron-Clausius equation (saturation evaporation
pressure) where  \code{x} is a data object holding the temperature.

\code{precip.vul}: and index for the vulerability to precipitation
defined as wetmean(x)/wetfreq(x). High when the mean intensity is high
and/or the frequency is low (it rains seldom, but when it rains, it
really pours down).

\code{t2m.vul}: and index for the vulerability to temperature defined as
the mean spell length for heat waves with temperatures exceeding 30C
(default).

\code{precip.rv}: a rough estimate of the return value for precipitation
under the assumption that it is exponentially distributed. Gives
apprximate answers for low return levels (less than 20
years). Advantage, can be predicted given wet-day mean and frequency.

\code{nv}: number of valid data points.

\code{precip.Pr}: rough estimate of the probability of more than x0 of
rain based on an exponential distribution.

\code{t2m.Pr}: rough estimate of the probability of more than x0 of rain
based on a normal distribution.

\code{NE}: predicts the number of events given the probability Pr. 
}
\usage{
C.C.eq(x)
precip.vul(x)
t2m.vul(x,x0=30,is=1)
precip.rv(x,tau=10)
nv(x)
precip.Pr(x,x0=10)
t2m.Pr(x,x0=10,na.rm=TRUE)
NE(p)
}

\arguments{
  \item{x}{a data object}
  \item{p}{a probability}
  \item{x0}{a threshold value}
  \item{tau}{time scale (years)}
  \item{is}{which of the spell results [1,2]}
  \item{na.rm}{See \code{\link{mean}}.}
}

\value{The right hand side of the equation}


\author{R. Benestad, MET Norway}
\examples{
t2m <- t2m.DNMI(lon=c(-70,-10),lat=c(20,60))
es <- C.C.eq(t2m)
map(es)
}
\keyword{parameter,element}