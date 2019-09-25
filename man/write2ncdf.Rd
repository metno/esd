\name{write2ncdf4}
\alias{write2ncdf4}
\alias{write2ncdf4.station}
\alias{write2ncdf4.default}
\alias{write2ncdf4.field}
\alias{write2ncdf4.list}
\alias{write2ncdf4.station}
\alias{write2ncdf4.eof}
\alias{write2ncdf4.pca}
\alias{write2ncdf4.dsensemble}
\title{Saves climate data as netCDF.}
\description{Method to save station data as netCDF, making sure to include
  the data structure and meta-data (attributes). The code tries to
  follow the netCDf 'CF' convention. The method is built on the
  \code{\link{ncdf4}} package.
}
\usage{
write2ncdf4(x,...) 
write2ncdf4.default(x,...)
write2ncdf4.list(x,file='field.nc',prec='short',scale=0.1,offset=NULL,
                 torg="1970-01-01",missval=-999,verbose=FALSE)
		 
write2ncdf4.field(x,file='field.nc',prec='short',scale=10,offset=NULL,
                 torg="1970-01-01",missval=-999,verbose=FALSE)
write2ncdf4.eof(x,file='eof.nc',prec='short',scale=10,offset=NULL,
                torg="1970-01-01",missval=-999)
write2ncdf4.pca(x,file='pca.nc',prec='short',scale=10,offset=NULL,
                torg="1970-01-01",missval=-999) 
write2ncdf4.dsensemble(x,file='dsensemble.nc',prec='short',scale=10,offset=NULL,
                       torg="1970-01-01",missval=-999) 
write2ncdf4.station(x,file='station.nc',prec='short',offset=0, missval=-999,
                    scale=0.1,torg='1899-12-31',nmin=30,verbose=FALSE) 
}

\arguments{
  \item{x}{data object}
  \item{file}{file name}
  \item{prec}{Precision: see \code{\link[ncdf4]{ncvar_def}}}
  \item{missval}{Missing value: see \code{\link[ncdf4]{ncvar_def}}}
  \item{offset}{Sets the attribute 'add_offset' which is added to the
  values stored (to save space may be represented as 'short').}
  \item{scale}{Sets the atttribute 'scale_factor' which is used to scale
  (multiply) the values stored (to save space may be represented as 'short').}
  \item{torg}{Time origin}
  \item{nmin}{Minimum number of years for saving valid summary statistics}
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
write2ncdf4(X,file='test.nc')
}
\keyword{netcdf}
\keyword{ncdf4}
\keyword{save}