\name{Data}
\alias{station.meta}
\alias{NACD}
\alias{NARP}
\alias{Oslo}
\alias{Svalbard}
\alias{t2m.NORDKLIM}
\alias{precip.NORDKLIM}
\alias{geoborders}
\alias{t2m.ERAINT}
\alias{t2m.ERA40}
\alias{t2m.NCEP}
\alias{t2m.MERRA}
\alias{sst.NCEP}
\alias{precip.ERAINT}
\alias{slp.ERAINT}
\alias{slp.MERRA}
\alias{slp.NCEP}
\alias{t2m.NorESM.M}
\alias{t2m.DNMI}
\alias{sst.DNMI}
\alias{slp.DNMI}
\alias{eof.t2m.ERAINT}
\alias{eof.t2m.ERA40}
\alias{eof.t2m.NCEP}
\alias{eof.t2m.MERRA}
\alias{eof.sst.NCEP}
\alias{eof.precip.ERAINT}
\alias{eof.slp.ERAINT}
\alias{eof.slp.MERRA}
\alias{eof.slp.NCEP}
\alias{eof.t2m.NorESM.M}
\alias{eof.t2m.DNMI}
\alias{eof.sst.DNMI}
\alias{eof.slp.DNMI}
\alias{NAOI}
\alias{sunspots}
\alias{NINO3.4}
\alias{bjornholt}
\alias{vardo}
\alias{ferder}
\alias{HadCRUT4}
\title{Sample data.}
\description{Different data sets: station data from northern Europe
  (NACD, NARP) and historic reconstructions (Oslo, Svalbard) from
  Dr. Nordli, Met Norway.

  The object \code{station.meta} contains station information, used in
  the methods \code{\link{station}}. 
  
  Also reduced representation of re-analyses, where the data have been
  sampled by skipping grid points to reduce the spatial dimensions and
  stored as 20 EOFS (30 for precipitation). The data compression
  facilitated by the EOFs can provide 80-90\% of the variance in the
  data. ESD uses the large-scale features from these reanalyses, and
  hence this information loss may be acceptable for downscaling work.

  A reduced copy of the NorESM (M RCP 4.5) is also provided for the examples
  and demonstrations on how the downscaling can be implemented. Note:
  downscaling for end-users should never be based on one GCM simulation
  alone. 

  The object \code{geoborders} contains data on coastlines and borders, used in
  the methods \code{\link{map}}.

}
\usage{
data(bjornholt)
data(ferder)
data(vardo)
data(eca.meta)
data(station.meta)
data(NACD)
data(NARP)
data(Oslo)
data(Svalbard)
data(eof.t2m.ERAINT)
data(eof.t2m.ERA40)
data(eof.t2m.NCEP)
data(eof.precip.ERAINT)
data(eof.slp.ERAINT)
data(eof.slp.MERRA)
data(eof.slp.NCEP)
data(eof.t2m.NorESM.M)
data(eof.t2m.DNMI)
data(eof.sst.DNMI)
data(eof.slp.DNMI)
data(geoborders)
slp.MERRA(lon=NULL,lat=NULL,anomaly=FALSE)
t2m.MERRA(lon=NULL,lat=NULL,anomaly=FALSE)
t2m.NCEP(lon=NULL,lat=NULL,anomaly=FALSE)
sst.NCEP(lon=NULL,lat=NULL,anomaly=FALSE)
slp.NCEP(lon=NULL,lat=NULL,anomaly=FALSE)
t2m.ERAINT(lon=NULL,lat=NULL,anomaly=FALSE)
precip.ERAINT(lon=NULL,lat=NULL,anomaly=FALSE)
slp.ERAINT(lon=NULL,lat=NULL,anomaly=FALSE)
t2m.ERA40(lon=NULL,lat=NULL,anomaly=FALSE)
t2m.DNMI(lon=NULL,lat=NULL,anomaly=FALSE)
slp.DNMI(lon=NULL,lat=NULL,anomaly=FALSE)
sst.DNMI(lon=NULL,lat=NULL,anomaly=FALSE)
t2m.NorESM.M(lon=NULL,lat=NULL,anomaly=FALSE)
data(NAOI)
data(sunspots)
data(NINO3.4)
HadCRUT4(url="http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.2.0.0.monthly_ns_avg.txt",plot=FALSE)
}

\arguments{
\item{lon}{longitude range c(lin.min,lon.max)}
\item{lat}{latitude range}
\item{anomaly}{TRUE: return anomaly}
\item{url}{source of data}
\item{plot}{TRUE:plot}
}

\seealso{
 \code{\link{spatial.avg.field}} \code{\link{as.4seasons}}, \code{\link{annual}}
}


\value{Numeric vectors/matrices with a set of attributes describing the data.}
\author{R.E. Benestad}
\examples{
data(Oslo)
year <- as.numeric( format(index(Oslo), '\%Y') ) 
plot(aggregate(Oslo, by=year, mean, na.rm = FALSE))
}
\keyword{datasets}
