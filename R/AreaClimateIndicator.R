#' Read daily station data of the Norwegian Meteorological Institute from thredds netCDF
#' using OpenDAP
#' 
#' \code{AreaClimateIndicator} returns area-based climate indicators provided through the Met Norway thredds server,
#' some of which are described in Benestad et al (2022)  
#' "Global hydro-climatological indicators and changes in the global hydrological cycle and rainfall patterns", 
#' PLOS Climate, PCLM-D-21-00079R1, DOI: 10.1371/journal.pclm.0000029, and an extension of those to temperature and other thresholds.
#' 
#' @aliases AreaClimateIndicator
#' @seealso NINO3.4, NAO, GSL, QBO, CET, CO2, AMO, IOD,  aggregate
#'
#' @param param The element to read c('t2m','tmax','tmin','precip')
#' @param is the area or index space to select domain c('Global','50S-50N','Arctic')
#' @param it Index time to select times scale c('as.annual','as.4seasons','as.monthly','day')
#' @param threshold select area of precipitation exceeding 1 mm/day or temperature above 0C (the freezing point)
#' @param FUN is the function used in aggregation 
#' @examples
#' ## Get the daily minimum temperature for Oslo-Blindern (station ID 18700)
#' x <- AreaClimateIndicator()
#' y <- AreaClimateIndicator(param='t2m',is='Arctic',threshold=0)
#' z <- AreaClimateIndicator(param='precip',is='50S-50N',it='day',threshold=1)
#' w <- AreaClimateIndicator(param='t2m',it='max',threshold=40)
#'
#' @export AreaClimateIndicator
AreaClimateIndicator <- function(param='precip',is='Global',it='annual',FUN='mean',threshold=1,plot=TRUE,verbose=FALSE,
                                 url='https://thredds.met.no/thredds/dodsC/metusers/rasmusb/') { 
  if (verbose) print(match.call())
  
  varid <- switch(tolower(param), 't2m'='TAS', 'precip'='PT', 'tmax'='TAX','tmin'='TAN')
  unit <- switch(tolower(param), 't2m'='C', 'precip'='mm', 'tmax'='C','tmin'='C')
  url <- paste0('https://thredds.met.no/thredds/dodsC/metusers/rasmusb/era5-area-',is,'-',
                varid,'.gt.',threshold,unit,'.nc')
  if (verbose) print(url) 
  
  ncid <- nc_open(url)
  vars <- names(ncid$var)
  varid <- vars[-grep('time_bnds',vars)]
  if (verbose) print(paste('Reading',varid))
  y <- ncvar_get(ncid,varid)
  t <- as.Date(ncvar_get(ncid,'time')/24,origin='1900-01-01')
  lname <- ncatt_get(ncid,varid,'long_name')$value
  nc_close(ncid)
  X <- zoo(y,order.by=t)
  #X <- subset(X,it=c(1981,2019))
  if (verbose) print(range(index(X)))
  if (is.character(it)) {
    if (tolower(it)=='annual') X <- annual(X,FUN=FUN) else 
      if (tolower(it)=='season') X <- as.4seasons(X,FUN=FUN) else
        if (tolower(it)=='month') X <- as.monthly(X,FUN=FUN)
  }
  
  if (plot) { 
    par(bty='n')
    plot(X,main=paste('ERA5:',is,'surface area fraction with',param,'exceeding',threshold,unit,'on a daily basis'),lwd=2,
         ylab='fraction',xlab=paste(range(index(X)),collapse=' - '),new=FALSE); grid()
  }
  attr(X,'url') <-url
  attr(X,'unit') <-unit
  attr(X,'long_name') <-lname
  attr(X,'threshold') <- threshold
  attr(X,'domain') <- is
  attr(X,'history') <- history.stamp()
  invisible(X)
}
