\name{rtools}
\alias{as.decimal}
\alias{nv}
\alias{cv}
\alias{q5}
\alias{q95}
\alias{q995}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simple and handy functions.
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
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data.frame or a coredata zoo object.}
   \item{na.rm}{If TRUE, remove NA's from data}
}
\details{
  'as.decimal' converts between degree-minute-second into decimal value.
  
  'cv' computes the coefficient of variation.
  
  'nv' count the number of valid data points.
  
  'q5','q95' and 'q995' are shortcuts to the 5\%, 95\%, and 99.5\%
  percentiles.
  
  'trend.coef' and 'trend.pval' return the coefficient and the p-value
  of the linear trend.
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
