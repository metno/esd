#' Fetch station data from GHND and convert into \code{station} objects.
#' 
#' Transform an input object into the esd class \code{station}. 
#' \code{GHCND.meta} reads the inventory and metadata of the stations (slow due to big data volume - requires Internet access).
#'
#' \code{GHCND.station} reads the raw data (slow due to big data volume - requires Internet access).

#' @aliases GHCND.meta GHCND.station
#' @seealso station as.station ghcnd.meta ghcnd.station
#' 
#' @param x a 'stationmeta' object
#' @param url lURL of the database
#' 
#' @examples
#' meta <- ghcnd.meta()
#' m <- subset(meta,lon=c(24,35),lat=c(-30,-25),verbose=TRUE)
#' y <- ghcnd.station(subset(m,is=1:10),verbose=TRUE)

#' @export
GHCND.meta <- function(url='https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/doc/ghcnd-stations.txt') {
  # meta <- read.fwf(url,widths = c(11,9,10,7,3,31,4,4,6),
  #                  col.names=c('ID','latitude','longitude','elevation','state','location','GSN-flag','HCN','WMO-ID'))
  ## GHCND
  # IV. FORMAT OF "ghcnd-stations.txt"
  # 
  # ------------------------------
  #   Variable   Columns   Type
  # ------------------------------
  #   ID            1-11   Character 
  # LATITUDE     13-20   Real
  # LONGITUDE    22-30   Real
  # ELEVATION    32-37   Real
  # STATE        39-40   Character
  # NAME         42-71   Character
  # GSN FLAG     73-75   Character
  # HCN/CRN FLAG 77-79   Character
  # WMO ID       81-85   Character
  # -----------------------------
  # https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
  # FIPS country codes https://en.wikipedia.org/wiki/List_of_FIPS_country_codes
  fips <- readLines('http://efele.net/maps/fips-10/data/fips-414.txt')
  fips <- fips[grep('country',fips)]
  fips <- gsub('00_414_414_country____',';',fips)
  fips <- gsub('__','',fips)
  fips <- unlist(strsplit(fips,';')); dim(fips) <- c(2,length(fips)/2)
  ## Get the GHCND metadata
  lines <- readLines(url)
  lines <- lines[nchar(lines)==85]
  id <- substr(lines,1,11)
  lat <- as.numeric(substr(lines,13,20))
  lon <- as.numeric(substr(lines,22,30))
  alt <- as.numeric(substr(lines,32,37))
  ## Extract country information from the station ID and the embedded FIPS codes
  cntr <- fips[2, match(substr(id,1,2),fips[1,])]
  qlty <- substr(lines,73,79)
  ele <- rep(100.901,length(id))
  param <- rep('precip+t2m',length(id))
  loc <- substr(lines,42,71)
  wmo.id <- substr(lines,81,85)
  meta <- data.frame(station_id=id,longitude=lon,latitude=lat,altitude=alt,country=cntr,location=loc,
                     wmo=wmo.id,source=rep('GHCND',length(id)),quality=qlty,element=ele,variable=param,
                     start=rep(NA,length(id)),end=rep(NA,length(id)))
  class(meta) <- c('stationmeta','data.frame')
  return(meta)
}

#' @export
GHCND.station <- function(x,url='https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/',verbose=FALSE) {
  n <- length(x$station_id)
  if (verbose) print(paste('Read',n,'stations'))
  Precip <- NULL; Tmax <- NULL; Tmin <- NULL; T2m <- NULL
  for (is in 1:n) {
    csvfile <- paste0(url,'/',x$station_id[is],'.csv')
    if (verbose) {print(x$location[is]); print(csvfile)}
    y <- read.csv(csvfile)
    t <- as.Date(y$DATE)
    if (!is.null(y$PRCP)) {
      if (sum(is.finite(y$PRCP))>0) { 
        if (verbose) print(paste('Precipitation:',paste(range(y$PRCP,na.rm=TRUE),collapse=' - ')))
        
        precip <- as.station(zoo(0.1*y$PRCP,order.by=t),loc=x$location[is],lon=x$longitude[is],lat=x$latitude[is],alt=x$altitude[is],
                             stid=x$station_id[is],wmo.is=x$wmo.id,variable='precip',unit='mm',cntr=x$country[is],
                             src='GHCND',url=csvfile)
        if (is.null(Precip)) Precip <- precip else Precip <- combine.stations(Precip,precip)
      }
    }
    if (!is.null(y$TMAX)) {
      if (sum(is.finite(y$TMAX))>0) { 
        if (verbose) print(paste('Tmax:',paste(range(y$TMAX,na.rm=TRUE),collapse=' - ')))
        tmax <- as.station(zoo(0.1*y$TMAX,order.by=t),loc=x$location[is],lon=x$longitude[is],lat=x$latitude[is],alt=x$altitude[is],
                           stid=x$station_id[is],wmo.is=x$wmo.id,variable='tmax',unit='degC',cntr=x$country[is],
                           src='GHCND',url=csvfile)
        if (is.null(Tmax)) Tmax <- tmax else Tmax <- combine.stations(Tmax,tmax)
      }
    }
    if (!is.null(y$TMIN)) {
      if (sum(is.finite(y$TMIN))>0) { 
        if (verbose) print(paste('Tmin:',paste(range(y$TMIN,na.rm=TRUE),collapse=' - ')))
        tmin <- as.station(zoo(0.1*y$TMIN,order.by=t),loc=x$location[is],lon=x$longitude[is],lat=x$latitude[is],alt=x$altitude[is],
                           stid=x$station_id[is],wmo.is=x$wmo.id,variable='tmin',unit='degC',cntr=x$country[is],
                           src='GHCND',url=csvfile)
        if (is.null(Tmin)) Tmin <- tmin else Tmin <- combine.stations(Tmin,tmin)
      }
    }
    if (!is.null(0.1*y$TAVG)) { 
      if (sum(is.finite(y$TAVG))>0) { 
        if (verbose) print(paste('T2m:',paste(range(y$TAVG,na.rm=TRUE),collapse=' - ')))
        t2m <- as.station(zoo(y$TAVG,order.by=t),loc=x$location[is],lon=x$longitude[is],lat=x$latitude[is],alt=x$altitude[is],
                          stid=x$station_id[is],wmo.is=x$wmo.id,variable='t2m',unit='degC',cntr=x$country[is],
                          src='GHCND',url=csvfile)
        if (is.null(T2m)) T2m <- t2m else T2m <- combine.stations(T2m,t2m)
      }
    }
  }
  return(list(precip=Precip,tmax=Tmax,tmin=Tmin,t2m=T2m))
}