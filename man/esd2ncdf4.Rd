\name{esd2ncdf4}
\alias{esd2ncdf4}
\alias{esd2ncdf4.default}
\alias{esd2ncdf4.field}
\alias{esd2ncdf4.station}
\alias{esd2ncdf4.eof}
\alias{esd2ncdf4.pca}
\alias{esd2ncdf4.dsensemble}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Export esd-objects to netCDF files}
\description{
\code{esd2ncdf4} saves the data and its main structure (metadata) as
netCDF files according to the CF-convention \url{http://cfconventions.org/}.  
}
\usage{
esd2ncdf4(x, ...)
esd2ncdf4.default(x,...)
esd2ncdf4.field(x,fname='field.nc',scale=10,offset=NULL,torg="1970-01-01")
esd2ncdf4.station(x,fname='station.nc',scale=10,offset=NULL,torg="1970-01-01")
esd2ncdf4.eof(x,fname='eof.nc',scale=10,offset=NULL,torg="1970-01-01")
esd2ncdf4.pca(x,fname='pca.nc',scale=10,offset=NULL,torg="1970-01-01")
esd2ncdf4.dsensemble(x,fname='dsensemble.nc',scale=10,offset=NULL,torg="1970-01-01")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{esd object}
  \item{fname}{file name}
  \item{scale}{scaling factor}
  \item{offset}{offset}
  \item{torg}{time origin}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, ...) 
UseMethod("esd2ncdf4")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
