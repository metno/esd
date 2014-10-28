\name{write2ncdf}
\alias{write2ncdf}
\alias{write2ncdf.station}
\title{Saves climate data as netCDF.}
\description{Method to save station data as netCDF, making sure to include
  the data structure and meta-data (attributes). The code tries to
  follow the netCDf 'CF' convention. The method is built on the
  \code{\link{ncdf4}} package.
}
\usage{
write2ncdf(x,file,...)
write2ncdf.station(x,file,prec='short',missval=-99.9,offs=0,scalf=0.1,
                   torg='1899-12-31',verbose=FALSE)
}

\arguments{
  \item{x}{data object}
  \item{file}{file name}
  \item{prec}{Precision: see \code{\link[ncdf4]{ncvar_def}}}
  \item{missval}{Missing value: see \code{\link[ncdf4]{ncvar_def}}}
  \item{offs}{Sets the attribute 'add_offset' which is added to the
  values stored (to save space may be represented as 'short').}
  \item{scalf}{Sets the atttribute 'scale_factor' which is used to scale
  (multiply) the values stored (to save space may be represented as 'short').}
  \item{torg}{Time origin}
  \item{verbose}{TRUE - clutter the screen.}
}
\value{A "zoo" "field" object with additional attributes used for further processing.}

\seealso{
       \code{\link{test.retrieve.ncdf4}}.
}

\author{R.E. Benestad}
\examples{
nacd <- station(src='nacd')
X <- annual(nacd)
write2ncdf(X,file='test.nc')
}
\keyword{?}
