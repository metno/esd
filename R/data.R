#' Sample data
#'
#' The object \code{geoborders} contains data on coastlines and borders, used
#' in the methods \code{\link{map}}.
#'
#' \code{etopo5} is a 5-minute gridded elevation data set provided by NOAA as described in 
#' "Data Announcement 88-MGG-02, Digital relief of the Surface of the Earth. NOAA, National Geophysical Data Center, Boulder, Colorado, 1988."
#'
#' The object \code{station.meta} contains meta data for various sources of station data
#' (NACD, NARP, NORDKLIM, ECAD, GHCN, and METNO) used in the methods \code{\link{station}}.
#'
#' \code{NACD}, \code{NARP}, and \code{nordklim.data} contain staion data from Northern Europe from
#' the North Atlantic Climatological Dataset (NACD),
#' the Nordic Arctic Research Programme (NARP),
#' and the NORDKLIM project, respectively,
#' which are used in the methods \code{\link{station}}.
#'
#' The temperature and precipitation data from NORDKLIM are also avaialble
#' as \code{station} objects in \code{t2m.NORDKLIM} and \code{precip.NORDKLIM}.
#'
#' \code{Oslo} and \code{Svalbard} are historic reconstructions of temperature from Oslo (1837-2018)
#' and Svalbard (1898-2018) provided by Dr. Nordli, Met Norway.
#'
#' \code{ferder} and \code{vardo} are time series of temperature and \code{bjornholt} of precipitation from stations in Norway,
#' downloaded from the Met Norway data archive. Ferder and Bjornholt are located near Olso
#' while Vardo is located in Northern Norway.
#' 
#' Samples of re-analyses are provided but to reduce the data size they have been stored as 20 EOFS (30 for precipitation)
#' To reconstruct the fields, use the functions \code{precip.ERAINT}, \code{slp.NCEP}, \code{t2m.NCEP}, \code{sst.NCEP},
#' \code{slp.DNMI}, \code{sst.DNMI}, and \code{t2m.DNMI}. 
#' The data compression facilitated by the EOFs can provide 80-90\% of the variance in the data.
#' ESD uses the large-scale features from these reanalyses, and hence this information loss may be acceptable for downscaling work.
#'
#' A reduced copy of the NorESM (M RCP 4.5) is also provided for the examples
#' and demonstrations on how the downscaling can be implemented. Note:
#' downscaling for end-users should never be based on one GCM simulation alone.
#'
#' \code{slp.ERA5} provides a small sample of 6-hourly ERA5 sea level pressure data from the North Atlantic
#' from September 30 to October 10 of 2016, used to test the \code{\link{CCI}} method.
#'
#' Some data sets (\code{NINO3.4}, \code{NAOI} ) come with a 'frozen' version in the package,
#' but there are also functions that read the most recent version of these indeces from the Internet
#' with functions \code{\link{NINO3.4}} and \code{\link{NAO}}.
#'
#' @aliases geoborders etopo5 station.meta NACD NARP Oslo Svalbard bjornholt ferder Svalbard Oslo vardo
#' nordklim.data station.meta t2m.NORDKLIM precip.NORDKLIM
#' eof.precip.ERAINT eof.slp.NCEP eof.sst.NCEP eof.t2m.NCEP
#' eof.slp.DNMI eof.sst.DNMI eof.t2m.DNMI
#' eof.t2m.NorESM.M
#' precip.ERAINT slp.NCEP sst.NCEP t2m.NCEP
#' slp.DNMI sst.DNMI t2m.DNMI
#' t2m.NorESM.M
#' slp.ERA5
#' NINO3.4 NAOI sunspots
#' arctic.t2m.cmip5   arctic.t2m.cmip3 global.t2m.cmip3 
#' global.t2m.gcm global.t2m.cmip5 scandinavia.t2m.cmip3 scandinavia.t2m.cmip5
#' IPCC.AR5.Table.9.A.1 gcmresolution
#' dse.ferder dse.Svalbard dse.Oslo
#' imilast.M03 storms 
#' mu.eq.f.tx
#'
#' @param lon longitude range c(lin.min,lon.max)
#' @param lat latitude range
#' @param anomaly TRUE: return anomaly
#' @param url source of data
#' @param plot TRUE:plot
#' @return Numeric vectors/matrices with a set of attributes describing the
#' data.
#' @author R.E. Benestad
#' @seealso \code{\link{aggregate.area}} \code{\link{as.4seasons}},
#' \code{\link{annual}}
#' @keywords datasets
#' @examples
#' 
#' data(Oslo)
#' year <- as.numeric( format(index(Oslo), '%Y') ) 
#' plot(aggregate(Oslo, by=year,FUN='mean', na.rm = FALSE), new=FALSE)
#' 
#' data(etopo5)
#' z <- subset(etopo5,is=list(lon=c(-10,30),lat=c(40,60)))
#' map(z, new=FALSE)
#' 
#' @export t2m.NCEP
t2m.NCEP <- function(lon=NULL,lat=NULL,anomaly=FALSE,latest=FALSE,
                     url='ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc',
                     verbose=FALSE) {
  if(verbose) print("t2m.NCEP")
  if (latest) {
    ## If the URL exists, check if a newer version should be downloaded
    if (file.exists('air.mon.mean.nc')) {
      ## If it exists locally, check if the two are the same version
      finfo.url <- file.info(url)
      finfo.loc <- file.info('air.mon.mean.nc')
      if (as.Date(finfo.loc$ctime) < as.Date(Sys.time()) - 7) {
          if (verbose) print(paste('download from',url))
          file.remove('air.mon.mean.nc')
          download.file(url,'air.mon.mean.nc')
      }
    } else download.file(url,'air.mon.mean.nc')
    t2m.NCEP <- retrieve('air.mon.mean.nc',lon=lon,lat=lat)
    if (anomaly) t2m.NCEP <- anomaly(t2m.NCEP)
  } else {
    data("eof.t2m.NCEP",envir=environment())
    t2m.NCEP<- eof2field(eof.t2m.NCEP,is=list(lon=lon,lat=lat),anomaly=anomaly)
  }
  return(t2m.NCEP)
}

#' @export sst.NCEP
sst.NCEP <- function(lon=NULL,lat=NULL,anomaly=FALSE,verbose=FALSE) {
  if(verbose) print("sst.NCEP")
  data("eof.sst.NCEP",envir=environment())
  sst.NCEP<- eof2field(eof.sst.NCEP,is=list(lon=lon,lat=lat),anomaly=anomaly)
  return(sst.NCEP)
}

#' @export slp.NCEP
slp.NCEP <- function(lon=NULL,lat=NULL,anomaly=FALSE,verbose=FALSE) {
  if(verbose) print("slp.NCEP")
  data("eof.slp.NCEP",envir=environment())
  slp.NCEP <- eof2field(eof.slp.NCEP,is=list(lon=lon,lat=lat),anomaly=anomaly)
  return(slp.NCEP)
}

#' @export precip.ERAINT
precip.ERAINT <- function(lon=NULL,lat=NULL,anomaly=FALSE,verbose=FALSE) {
  if(verbose) print("precip.ERAINT")
  data("eof.precip.ERAINT",envir=environment())
  precip.ERAINT <- eof2field(eof.precip.ERAINT,is=list(lon=lon,lat=lat),anomaly=anomaly)
  return(precip.ERAINT)
}

#' @export t2m.DNMI
t2m.DNMI <- function(lon=NULL,lat=NULL,anomaly=FALSE,verbose=FALSE) {
  if(verbose) print("t2m.DNMI")
  data("eof.t2m.DNMI",envir=environment())
  t2m.DNMI <- eof2field(eof.t2m.DNMI,is=list(lon=lon,lat=lat),anomaly=anomaly)
  return(t2m.DNMI)
}

#' @export slp.DNMI
slp.DNMI <- function(lon=NULL,lat=NULL,anomaly=FALSE,verbose=FALSE) {
  if(verbose) print("slp.DNMI")
  data("eof.slp.DNMI",envir=environment())
  slp.DNMI <- eof2field(eof.slp.DNMI,is=list(lon=lon,lat=lat),anomaly=anomaly)
  return(slp.DNMI)
}

#' @export sst.DNMI
sst.DNMI <- function(lon=NULL,lat=NULL,anomaly=FALSE,verbose=FALSE) {
  if(verbose) print("sst.DNMI")
  data("eof.sst.DNMI",envir=environment())
  sst.DNMI <- eof2field(eof.sst.DNMI,is=list(lon=lon,lat=lat),anomaly=anomaly)
  return(sst.DNMI)
}

#' @export t2m.NorESM.M
t2m.NorESM.M <- function(lon=NULL,lat=NULL,anomaly=FALSE,verbose=FALSE) {
  if(verbose) print("t2m.NorESM.M")
  data("eof.t2m.NorESM.M",envir=environment())
  t2m.NorESM.M <- eof2field(eof.t2m.NorESM.M,is=list(lon=lon,lat=lat),anomaly=anomaly)
  return(t2m.NorESM.M)
}
