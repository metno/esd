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
\alias{SOI}
\alias{GSL}
\alias{bjornholt}
\alias{vardo}
\alias{ferder}
\alias{dse.ferder}
\alias{HadCRUT4}
\alias{NASAgiss}
\alias{dse.Oslo}
\alias{glossstations}
\alias{IPCC.AR5.Table.9.A.1}
\alias{global.t2m.gcm}
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

  \code{glossstations} contains META-data for GLOSS stations taken from
  the table in \url{http://www.gloss-sealevel.org/station_handbook/stations/#.Vqtc6kL4phg}
  
  Some data sets, such as NINO3.4 and NAOI come with a 'frozen' version
  in the package, but there are also functions that read the most recent
  version of these indeces from the Internet. GSL reads the global mean sea level.

}
\usage{
data(bjornholt)
data(ferder)
data(dse.ferder)
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
data(IPCC.AR5.Table.9.A.1)
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
data(sunspots)
data(NINO3.4)
data(NAOI)
NINO3.4(url='ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/ersst3b.nino.mth.ascii',header=TRUE)
NAO(url='http://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table',header=FALSE)
SOI(url='ftp://ftp.bom.gov.au/anon/home/ncc/www/sco/soi/soiplaintext.html',header=FALSE)
data(dse.Oslo)
HadCRUT4(url="http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.2.0.0.monthly_ns_avg.txt",plot=FALSE)
NASAgiss(url='http://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.txt')
GSL <-
function(url='http://www3.epa.gov/climatechange/images/indicator_downloads/sea-level_fig-1.csv')
data(global.t2m.gcm)
}

\arguments{
\item{lon}{longitude range c(lin.min,lon.max)}
\item{lat}{latitude range}
\item{anomaly}{TRUE: return anomaly}
\item{url}{source of data}
\item{plot}{TRUE:plot}
}

\seealso{
 \code{\link{aggregate.area}} \code{\link{as.4seasons}}, \code{\link{annual}}
}


\value{Numeric vectors/matrices with a set of attributes describing the data.}
\author{R.E. Benestad}
\examples{
data(Oslo)
year <- as.numeric( format(index(Oslo), '\%Y') ) 
plot(aggregate(Oslo, by=year,FUN='mean', na.rm = FALSE))

data(etopo5)
z <- subset(etopo5,is=list(lon=c(-10,30),lat=c(40,60)))
map(z)

}
\keyword{datasets}
