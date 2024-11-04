#' Retrieve station record from a given data source.
#'
#' @seealso [func()] station station.default station.ecad station.nacd
#' station.narp station.nordklim station.metnod station.metnom station.ghcnd
#' station.ghcnm station.ghcnm station.sonel station.gloss station.newlyn
#' station.giss clean.station allgood station.thredds map.station select.station
#'
#' @param x a data.frame with GHCND metadata
#' @param param Parameter or element type or variable identifier. There are
#' several core parameters or elements as well as a number of additional
#' parameters. The parameters or elements are: precip = Precipitation (mm) tas,
#' tavg = 2m-surface temperature (in degrees Celcius) tmax, tasmax = Maximum
#' temperature (in degrees Celcius) tmin, tasmin = Minimum temperature (in
#' degrees Celcius)
#' @param url URL of the GHCND database
#' @param stid A string of characters as an identifier of the weather/climate
#' station.
#' @param lon Numeric value of longitude (in decimal degrees East) for the
#' reference point (e.g. weather station) as a vector
#' containing the range of longitude values in the form of c(lon.min,lon.max)
#' @param lat Numeric value of latitudes for the range of latitude values in the
#' form of c(lat.min,lat.max)
#' @param alt Numeric value of altitude (in meters a.s.l.) used for selection.
#' Positive value, select all stations above this altitude; for negative
#' values, select all stations below this latitude.
#' @param cntr A string or a vector of strings of the full name of the country:
#' Select the stations from a specified country or a set of countries.
#' @param plot Logical value. If, TRUE provides a plot.
#' @param verbose Logical value defaulting to FALSE. If FALSE, do not display
#' comments (silent mode). If TRUE, displays extra information on progress.
#' @return A time series of "zoo" "station" class with additional attributes
#' used for further processing.
#'
#' @author Rasmus Benestad
#' @keywords station.GHCND
#' @examples
#' ## Read rain gauge data from Mozambique
#' ## It may take a little time reading all the metadata over the Internet
#' meta.GHCND(cntr='Mozambique',verbose=TRUE) -> mz
#' Y <- station.GHCND(mz,param='precip')
#' plot(Y)
#'
#' meta.GHCND(cntr='Ghana',verbose=TRUE) -> gz
#' ## Get all available variables: returned as a list object where 
#' X <- station.GHCND(gz)
#' print(names(gz))
#' plot(X$precip)
#' 
#' @exportS3Method
#' @export station.GHCND
station.GHCND <- function(x=NULL,cntr=NULL,param=NULL,lon=NULL,lat=NULL,
                          url='https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access',
                          sep=',',verbose=FALSE) {

  if (verbose) print('station.GHCND')
  if (is.null(x)) {
    if (is.null(cntr)) cntr <- 'Mozambique'
    if (is.null(param)) param <- 'precip' 
    x <- meta.GHCND(cntr=cntr,param,verbose=verbose)
  }
  filenames <- paste0(url,'/',gsub(' ','',x$station_id),'.csv')
  Precip <- NULL; Tmax <- NULL; Tmin <- NULL; T2m <- NULL
  iso2 <- ISO2cntrcode()
  ii <- 1
  for (file2get in filenames) {
    if (verbose) print(sub(paste0(url,'/'),'',file2get))
    ghcnd <- read.table(file2get,sep=sep,header=TRUE)
    content <- names(ghcnd)
    if (length(grep('PRCP',content))>0) { 
      precip <- zoo(x=ghcnd$PRCP,order.by=as.Date(ghcnd$DATE))
      precip <- as.station(precip,stid=ghcnd$STATION[1],loc=ghcnd$NAME[1],
                           lon=ghcnd$LONGITUDE[1],lat=ghcnd$LATITUDE[1],
                           alt=ghcnd$ELEVATION,cntr=x$country[ii],
                           param='precip',unit='mm',longname='24-hr precipitation',
                           src='GHCN',url=url)
      if (is.null(Precip)) Precip <- precip else Precip <- combine.stations(Precip,precip)
    }
    if (length(grep('TMAX',content))>0) { 
      tmax <- zoo(x=ghcnd$TMAX-273.15,order.by=as.Date(ghcnd$DATE))
      tmax <- as.station(tmax,stid=ghcnd$STATION[1],loc=ghcnd$NAME[1],
                         lon=ghcnd$LONGITUDE[1],lat=ghcnd$LATITUDE[1],
                         alt=ghcnd$ELEVATION,cntr=x$country[ii],
                         param='tmax',unit='degC',longname='daily maximum temperature',
                         src='GHCN',url=url)
      if (is.null(Tmax)) Tmax <- tmax else Tmax <- combine.stations(Tmax,tmax)
    }
    if (length(grep('TMIN',content))>0) { 
      tmin <- zoo(x=ghcnd$TMIN-273.15,order.by=as.Date(ghcnd$DATE))
      tmin <- as.station(tmax,stid=ghcnd$STATION[1],loc=ghcnd$NAME[1],
                         lon=ghcnd$LONGITUDE[1],lat=ghcnd$LATITUDE[1],
                         alt=ghcnd$ELEVATION,cntr=x$country[ii],
                         param='tmin',unit='degC',longname='daily minimum temperature',
                         src='GHCN',url=url)
      if (is.null(Tmin)) Tmin <- tmin else Tmin <- combine.stations(Tmin,tmin)
    }
    if (length(grep('TAVG',content))>0) { 
      t2m <- zoo(x=ghcnd$TAVG-273.15,order.by=as.Date(ghcnd$DATE))
      t2m <- as.station(t2m,stid=ghcnd$STATION[1],loc=ghcnd$NAME[1],
                        lon=ghcnd$LONGITUDE[1],lat=ghcnd$LATITUDE[1],
                        alt=ghcnd$ELEVATION,cntr=x$country[ii],
                        param='t2m',unit='degC',longname='daily average temperature',
                        src='GHCN',url=url)
      if (is.null(T2m)) T2m <- t2m else T2m <- combine.stations(T2m,t2m)
    }
    ii <- ii + 1
  }
  if (is.null(param)) result <- list(precip=Precip,tmax=Tmax,tmin=Tmin,t2m=T2m) else
    result <- switch(tolower(param),'precip'=Precip,'tmax'=Tmax,'tmin'=Tmin,'t2m'=T2m)
  invisible(result)
}

#' @exportS3Method
#' @export meta.GHCND
meta.GHCND <- function(url='https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt',
                       urlinv='https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt',
                       urlcntr='https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-countries.txt',
                       cntr=NULL,param='precip',lon=NULL,lat=NULL,alt=NULL,
                       widths=c(12,9,10,7,34,4,10),
                       metaID=c('station_id','latitude','longitude','altitude',
                               'location','GSN','ID'),verbose=FALSE,plot=FALSE) {
  require(dplyr)
  require(tidyverse)
  if (verbose) print('meta.GHCND')
  ## Read the metadata
  meta <- read.fwf(url,widths=widths,comment.char = "")
  names(meta) <- metaID
  meta$station_id <- gsub(' ','',meta$station_id)
  inventory <- read.table(urlinv)
  ## Read station inventory
  names(inventory) <- c('station_id','latitude','longitude','variable','start','end')
  ncs <- max(nchar(readLines(urlcntr))) - 2
  cntrcode <- read.fwf(urlcntr,widths=c(2,ncs))
  ## Add country information to the meta data list
  names(cntrcode) <- c('cntrcode','country')
  cntrs <- substr(meta$station_id,1,2); country=cntrs
  Cntrs <- table(cntrs)
  if (verbose) print(Cntrs)
  for (ic in rownames(Cntrs)) {
    if (verbose) print(paste(ic,cntrcode$country[cntrcode$cntrcode==ic]))
    country[cntrs==ic] <- cntrcode$country[cntrcode$cntrcode==ic]
  }
  if (verbose) print(table(country))
  meta <- cbind(meta,country)
  ## Add parameter and start/end to the meta data list
  ns <- length(meta$station_id)
  inv <- data.frame(params=rep('NA',ns),start=rep(NA,ns),end=rep(NA,ns))
  ii <- 0
  if (verbose) print(paste(ns,'GHCND stations'))
  # for (id in meta$station_id) {
  #   idmatch <- is.element(inventory$station_id,id)
  #   iis <- match(id,meta$station_id)
  #   inv$params[iis] <- paste(inventory$variable,collapse='|')
  #   inv$start[iis] <- min(inventory$start,na.rm=TRUE)
  #   inv$end[iis] <- min(inventory$end,na.rm=TRUE)
  #   if (ii %% 100 == 0) cat('.')
  #   ii <- ii + 1
  # }
  #meta <- cbind(meta,inv)
  
  # Grouping the data by station_id and summarise the data 
  inventory_summary <- inventory %>%
    group_by(station_id) %>%
    summarise(params = paste(variable, collapse = '|'),
              start = min(start, na.rm = TRUE),
              end = max(end, na.rm = TRUE))
  # Joining the meta data with the inventory_summary
  meta <- left_join(meta, inventory_summary, by = "station_id")
  if (verbose) print(names(meta))
  if (!is.null(cntr)) {
    if (verbose) print(paste('Select',cntr))
    sel <- grep(tolower(cntr),tolower(meta$country))
    meta <- meta[sel,]
  }
  if (!is.null(lon)) {
    if (verbose) print(paste('Select longitudes',paste(lon,collapse='-')))
    sel <- (meta$longitude >= min(lon)) & (meta$longitude <= max(lon))
    meta <- meta[sel,]
  }
  if (!is.null(lat)) {
    if (verbose) print(paste('Select longitudes',paste(lat,collapse='-')))
    sel <- (meta$latitude >= min(lat)) & (meta$latitude <= max(lat))
    meta <- meta[sel,]
  }
  if (!is.null(alt)) {
    if (verbose) print(paste('Select altitudes',paste(alt,collapse='-')))
    if (length(alt)>1) sel <- (meta$altitude >= min(alt)) & (meta$altitude <= max(alt)) else
      if (alt <0) sel <- (meta$altitude <= abs(alt)) else
        sel <- (meta$altitude >= alt)
    meta <- meta[sel,]
  }
  
  if (plot) {
    plot(meta$longitude,meta$latitude,pch=19,col=rgb(0.5,0,0,0.2),
                 main='GHCND',cex=0.5)
    data("geoborders")
    lines(geoborders,col='grey')
  }
  if (verbose) print(dim(meta))
  class(meta) <- c("stationmeta","data.frame")
  invisible(meta)
}