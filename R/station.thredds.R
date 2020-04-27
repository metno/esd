#' Read daily station data of the Norwegian Meteorological Institute from thredds netCDF
#' using OpenDAP
#' 
#' \code{station.thredds} is a wrapper that uses retrieve.station combined with information about the
#' files stored on Thredds. The analysis can also be applied to either EOFs or fields.
#'
#' \code{meta.thredds} retrieves meta data from Thredds.  
#' 
#' @aliases station.thredds meta.thredds
#' @seealso retrieve.station, station, radar
#'
#' @param param The element to read c('t2m','tmax','tmin','precip','slp','sd','fx','fg','dd')
#' @param is Index space to select station (list)
#' @param stid Station ID to select
#' @param loc Name of location to select
#' @param lon Range of longitudes to select
#' @param lat Range of latitudes to select
#' @param alt Range of altitudes to select of stations above (positive) or below (negative) a threshold
#' @param it Range of times to select
#' @param cntr Countries to select
#' @param start.year.before Select stations with record starting before a given year
#' @param end.year.after Select stations with record ending after a given year
#' @param nmin Select stations with minimum number of valid data points
#' @param onebyone If many stations, select them first individually and then combine
#' @param verbose If TRUE print information about progress.
#' @param ... Other arguments.
#' @return A station object.
#'
#' @importFrom stats cov cor
#'
#' @examples
#' ## Get the daily minimum temperature for Oslo-Blindern (station ID 18700)
#' tmin <- station.thredds(param='tmax',stid=18700)
#'
#' meta <- meta.thredds(param='precip')
#' precip <- station.thredds(meta[1,])
#'
#' @export station.thredds
station.thredds <- function(param='t2m',is = NULL, stid = NULL, 
                            loc = NULL, lon = NULL, lat = NULL, it = NULL, alt = NULL, 
                            cntr = NULL, start.year.before = NULL, end.year.after = NULL, 
                            nmin = NULL, verbose = FALSE, onebyone = FALSE, ...) {
  if (verbose) t0 <- Sys.time()
  if (is.character(param)) { 
    if (param=='slp') param <- 'pp'
    url <- paste0('https://thredds.met.no/thredds/dodsC/metusers/rasmusb/',
                  param,'.metnod.nc')
    if (param=='pp') param <- 'slp'
  } else if (inherits(param,'stationsummary')) {
    stid <- param$station.id
    url <- attr(param,'url')
    param <- attr(param,'variable')
  }
  if (verbose) print(url)
  y <- retrieve.station(url,param=param,is=is,stid=stid,loc=loc,
                        lon=lon,lat=lat,it=it,alt=alt,cntr=cntr,
                        start.year.before = start.year.before, 
                        end.year.after = end.year.after,
                        nmin = nmin, verbose = verbose, onebyone = onebyone, ...)
  if (verbose) print(paste('Time taken for reading the data was', 
                           round(Sys.time() - t0),'s'))
  attr(y,'url') <- url
  invisible(y)
}

#' @export meta.thredds
meta.thredds <- function(param='t2m',verbose=FALSE) {
  Y <- NULL
  url <- NULL
  if (param=='slp') param <- 'pp'
  url <- paste0('https://thredds.met.no/thredds/dodsC/metusers/rasmusb/',
                param,'.metnod.nc')
  if (param=='pp') param <- 'slp'
  if (verbose) print(url)
  Y <- retrieve.stationsummary(url)
  attr(Y,'url') <- url
  attr(Y,'variable') <- param
  invisible(Y)
}

## wraparound
select.station.thredds <- function(x=NULL, ..., loc=NULL, param=NULL,  ele=NULL, stid=NULL, 
                                   lon=NULL, lat=NULL, alt=NULL, cntr=NULL, src=NULL, it=NULL, 
                                   nmin=NULL, user='external', verbose=FALSE) {
  if (is.null(param)) param <- c('t2m','tmin','tmax','precip','slp')
  for (i in 1:length(param)) {
    meta <- meta.thredds(param=param[i])
    if (i ==1) ss <- meta else ss <- cbind(ss,meta)
  }
  invisible(ss)
}
