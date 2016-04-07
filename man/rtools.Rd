\name{rtools}
\alias{as.decimal}
\alias{nv}
\alias{cv}
\alias{q5}
\alias{q95}
\alias{q995}
\alias{lag.station}
\alias{lag.field}
\alias{filt}
\alias{filt.default}
\alias{exit}
\alias{figlab}
\alias{ndig}
\alias{attrcp}
\alias{ensemblemean}
\alias{propchange}
\alias{TGW}
\alias{geostrophicwind}
\alias{geostrophicwind.station}
\alias{geostrophicwind.field}
\alias{stand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simple and handy functions.
}
\description{
\code{lag.station} and \code{lag.station} are wrap-around functions for
\code{\link[zoo]{lag.zoo}} that maintains all the attributes.
\code{attrcp(x,y)} passes on attributes from \code{x} to \code{y} and
returns the \code{y} with new attributes.
\code{ensemblemean} returns the ensemble mean for dsensemble
objects. The argument \code{FUN} can also be used to estiamte other
statistics, e.g \code{FUN='q95'} where \code{q95=function(x) apply(x,1,quantile,probs=0.95)}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{TGW} uses triangulation of pressure measurements to estimate
geostrophic wind based on Alexandersson et al. (1998), Glob. Atm. and
Oce. Sys.
\code{stand} gives a standardised time series.
}
\usage{
as.decimal(x = NULL)
cv(x,na.rm=TRUE)
nv(x)
q5(x)
q95(x)
q995(x)
trend.coef(x)
trend.pval(x)
lag.station(x,...)
lag.field(x,...)
exit()
ndig(x)
filt(x,...)
filt.default(x,n,type='ma',lowpass=TRUE)
figlab(x,xpos=0.001,ypos=0.001)
attrcp(x, y, ignore = c("name", "model", "n.apps",
"appendix","dimnames"))
ensemblemean(x,FUN='rowMeans')
propchange(x,it0=c(1979,2013))
TGW(triangle,f=1.25e-4,rho=1.25,verbose=FALSE)
geostrophicwind.station(x,f=1.25e-4,rho=1.25,verbose=FALSE)
geostrophicwind.field(x,f=1.25e-4,rho=1.25,verbose=FALSE)
stand(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data.frame or a coredata zoo object.}
   \item{na.rm}{If TRUE, remove NA's from data}
   \item{type}{'ma' for moving average (box-car), 'gauss' for Gaussian,
     'binom' for binomial filter,' parzen' for Parzen filter, 'hanning'
     for Hanning filter, or 'welch' for Welch filter.}
   \item{lowpass}{True for smoothing, otherwise the highpass results is
     returned}
   \item{triangle}{a group of three stations with sea-level
     pressure', e.g. from ECA\&D.}
}
\details{
  'as.decimal' converts between degree-minute-second into decimal value.
  
  'cv' computes the coefficient of variation.
  
  'nv' count the number of valid data points.
  
  'q5','q95' and 'q995' are shortcuts to the 5\%, 95\%, and 99.5\%
  percentiles.
  
  'trend.coef' and 'trend.pval' return the coefficient and the p-value
  of the linear trend.

  'exit' is a handy function for exiting the R session without saving.

  'figlab' is a handy function for labelling figures (e.g. 'Figure 1')

  'ndig' estimates the number of digits for round(x,ndig), e.g. in
  scales for plotting.
}
\value{
 \item{as.decimal }{Decimal value}
 \item{trend.coef }{Linear trend per decade}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{A. Mezghani}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Monthly mean temperature at 'Oslo - Blindern' station 
data(Oslo)
## Compute the linear trend and the p-value on annual aggregated values 
tr <- trend.coef(coredata(annual(Oslo)))
pval <- trend.pval(coredata(annual(Oslo)))
\dontrun{
pp <- station(param='slp',cntr='Denmark',src='ecad')
wind <- TGW(subset(pp,is=c(1,3,10))
plot(wind)
ws <- sqrt(wind[,1]^2 + wind[,2]^2)
plot(ws)
hist(ws)

## Estimate wind for a larger group of stations
wind <- geostrophicwind(pp,nmax=10)
u <- subset(wind,param='u')
v <- subset(wind,param='u')
ws <- sqrt(u^2+v^2)
ws <- attrcp(v,ws)
class(ws) <- class(v)
attr(ws,'variable')='windspeed'
attr(ws,'longname')='geostrophic wind speed'
map(ws,FUN='quantile',probs=0.98)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{rtools}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
