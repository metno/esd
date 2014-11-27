\name{WG}
\alias{Weather generators}
\alias{WG}
\alias{WG.station}
\alias{}
\alias{}
\title{Weather generators for conditional simulation.}
\description{ 
Weather generators for conditional simulation of daily temperature
and/or precipitation, given mean and/or standard deviation. The family
of WG functions procude stochastic time series with similar
characteristics as the station series provided (if none if provided, it
will use either ferder or bjornholt provided by the esd-package). Here
characteristics means similar mean value, standard deviation, and
spectral properties. \code{FTscramble} takes the Fourier components
(doing a Fourier Transform - FT) of a series and reassigns random phase
to each frequency and then returns a new series through an inverse
FT. The FT scrambling is used for temperature, but not for precipitation
that is non-Gaussian and involves sporadic events with rain. For
precipitation, a different approach is used, taking the wet-day
frequency of each year and using the wet-day mean and ranomly generated
exponentially distributed numbers to provide similar aggregated annual
statistics as the station or predicted though downscaling. The
precipitation WG can also take into account the number of consequtive
number-of-dry-days statistics, using either a Poisson or a gemoetric
distribution.

The weather generater produces a series with similar length as the
provided sample data, but with shifted dates according to specified
scenarios for annual mean mean/standard deviation/wet-day mean/wet-day
frequency.

\code{WG.FT.day.t2m} generates daily temperature from seasonal means and
standard deviations. It is given a sample station series, and uses
\code{FTscramble} to generate a series with random phase but similar (or
predicted - in the future) spectral characteristics. It then uses
a quantile transform to prescribe predicted mean and standard deviation,
assuming the distributions are normal. The temperal structure (power
spectrum) is therefore similar as the sample provided.

\code{WG.fw.day.precip} uses the annual wet-day mean and the wet-day frequency
as input, and takes a sample station of daily values to stochastically simulate
number consequtive wet days based on its annual mean number. If not
specified, it is taken from the sample data after being phase scrambeled
(\code{FTscramble}) The number of wet-days per year is estimated from
the wed-day frequency, it too taken to be phase scrambled estimates from
the sample data unless specifically specified. The
daily amount is taken from stochastic values generated with
\code{\link{rexp}}. The number of consequtive wet days can be
approximated by a geometric distribution (\code{\link{rgeom}}), and the
annual mean number was estimated from the sample series. 
}
\usage{
FTscramble(x,interval=NULL,spell.stats=FALSE,wetfreq.pred=FALSE)
WG(x,...)
WG.station(x,option='default',...)
WG.FT.day.t2m(x=NULL,amean=NULL,asd=NULL,t=NULL,eofs=1:4,
          select=NULL,lon=c(-20,20),lat=c(-20,20),
          plot=TRUE,biascorrect=TRUE,verbose=FALSE)
WG.fw.day.precip(x=NULL,mu=NULL,fw=NULL,ndd=NULL,t=NULL,
             threshold=1,select=NULL,
             lon=c(-10,10),lat=c(-10,10),
             plot=TRUE,biascorrect=TRUE,
             method='rgeom')
WG.pca.day.t2m.precip(t2m=NULL,precip=NULL,threshold=1,select=NULL)
}
\arguments{
    \item{x}{station object}
    \item{option}{Define the type of WG}
    \item{amean}{annual mean values. If NULL, use those estimated from
      x; if NA, estimate using \code{\link{DSensemble.t2m}}, or if
      provided, assume a 'dsensemble' object.}
    \item{asd}{annual standard deviation. If NULL, use those estimated from
      x; if NA, estimate using \code{\link{DSensemble.t2m}}, or if
      provided, assume a 'dsensemble' object.}
    \item{t}{Time axis. If null, use the same as x or the last interval
      of same length as x from downscaled results.}
    \item{eofs}{passed on to \code{\link{DSensemble.t2m}}}
    \item{select}{passed on to \code{\link{DSensemble.t2m}}}
    \item{lon}{passed on to \code{\link{DSensemble.t2m}}}
    \item{lat}{passed on to \code{\link{DSensemble.t2m}}}
    \item{plot}{}
    \item{biascorrect}{passed on to \code{\link{DSensemble.t2m}}}
    \item{verbose}{passed on to \code{\link{DSensemble.t2m}}}
    \item{mu}{annual wet-mean values. If NULL, use those estimated from
      x; if NA, estimate using \code{\link{DSensemble.t2m}}, or if
      provided, assume a 'dsensemble' object.}
    \item{fw}{annual wet-day frequency. If NULL, use those estimated from
      x; if NA, estimate using \code{\link{DSensemble.t2m}}, or if
      provided, assume a 'dsensemble' object.}
    \item{ndd}{annual mean dry spell length. If NULL, use those estimated from
      x; if NA, estimate using \code{\link{DSensemble.t2m}}, or if
      provided, assume a 'dsensemble' object.}
    \item{threshold}{Definition of a rainy day.}
    \item{method}{Assume a gemoetric or a poisson distribution. Can also
    define ownth methods.}
    \item{t2m}{station object with temperature}
    \item{precip}{station object with precipitation.}
 }
\value{}
\author{R.E. Benestad}
\examples{
data(ferder)
t2m <- WG(ferder)
data(bjornholt)
pr <- WG(bjornholt)
}
\keyword{manip}
