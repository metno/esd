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
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simple and handy functions.

\code{lag.station} and \code{lag.station} are wrap-around functions for
\code{\link[zoo]{lag.zoo}} that maintains all the attributes.
\code{attrcp(x,y)} passes on attributes from \code{x} to \code{y} and
returns the \code{y} with new attributes.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
attrcp(x, y, ignore = c("name", "model", "n.apps", "appendix","dimnames"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data.frame or a coredata zoo object.}
   \item{na.rm}{If TRUE, remove NA's from data}
   \item{type}{'ma' for moving average (box-car), 'gauss' for Gaussian,
     'binom' for binomial filter,' parzen' for Parzen filter, 'hanning'
     for Hanning filter, or 'welch' for Welch filter.}
   \item{lowpass}{True for smoothing, otherwise the highpass results is returned}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{rtools}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
