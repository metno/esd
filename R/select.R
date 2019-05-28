#' Select from meta data base
#'
#' Function that searches the meta data base for the requested station data
#' Search priority: ID, name, coordinates, altitude, country,...
#' Can return several matches 
'
#' @export
select <- function(stid=NULL,param=NULL,lon=NULL,lat=NULL,alt=NULL,cntr=NULL,
                    ...) UseMethod("select")

select.default <- function() {
}


# Selects a station time series based on station number or location (stid) or on the basis
# of location or altitude. Returns a pointer to selected objects in the meta-data list. Need to use
# this list to determine source.

# author : Abdelkader Mezghani
# date   : 30-05-2013
#' @export
select.station <- function (stid = NULL, param = NULL, lon = NULL, lat = NULL, 
    alt = NULL, cntr = NULL, src = NULL,silent=FALSE) 
{
    load("esd/data/station.meta.rda")
    n <- length(station.meta$stid)
    #STID <- station.meta$stid
    #browser()
    if (!is.null(param)) {
        ele.c <- switch(tolower(param), t2m = "101", tg = "101", 
            rr = "601", slp = "401", cloud = "801", t2 = "101", 
            precip = "601", `101` = "101", `401` = "401", `601` = "601",`801` = "801")
        station.meta <- subset(station.meta,element==ele.c)
    }
    if (!is.null(stid)) {
        if (is.numeric(stid)) {
            id <- stid
            station.meta <- subset(station.meta,stid==id)
        }
        else if (is.character(stid)) {
            id <- grep(tolower(stid),tolower(station.meta$location))
            station.meta <- station.meta[id,]
        }
    }
    if (!is.null(lon) & !is.null(lat)) {
        if ((length(lon) == 1) & (length(lat) == 1)) {
            d <- distAB(lon, lat, station.meta$lon, station.meta$lat)
            id <- d==min(d,na.rm=TRUE) ; id[is.na(id)] <- FALSE # Ak some of the lon values are NA's
            station.meta <- subset(station.meta,id)
        }
        else if ((length(lon) == 2) & (length(lat) == 2)) {
            lon.rng <- lon ; lat.rng <- lat 
            station.meta <- subset(station.meta,((lon >= lon.rng[1]) & (lon <= lon.rng[2])) & (lat >= lat.rng[1]) & (lat <= lat.rng[2]))           
        }
        else return(NULL)
    }
    if (!is.null(alt)) {
        if (length(alt) == 1) {
           alt.rng <- c(alt-0.1*alt,alt+.1*alt) # set the altitude range to +/-10% of alt
        }
        else if (length(alt) == 2) {
           alt.rng <- alt      
        }
        station.meta <- subset(station.meta,(alt >= alt.rng[1]) & (alt <= alt.rng[2]))   
    }
    if (!is.null(cntr)) {
        id <- is.element(tolower(station.meta$country),tolower(cntr))
        station.meta <- subset(station.meta,id)   
    }
    if (!is.null(param)) {
        id <- is.element(station.meta$element,ele.c)
        station.meta <- subset(station.meta,id)   
    }
    if (!is.null(src)) {
        id <- is.element(station.meta$src, src) 
        station.meta <- subset(station.meta,id)
    }
    if (dim(station.meta)[1]!=0)
       return(station.meta)
    else {
       if (silent) print("No available stations for your selection")
       return(NULL)
    }	
}

test.select.station <- function() {

# RUN !
# Available ECA&D stations for the range of longitude between 0 and 10 deg. East, the range of latitude between 20 and 30 deg. North and for altitudes around 400m +/- 40 m
available.station <- select.station2(lon = c(0,10),lat=c(20,30),src="ECA&D",alt=560,silent=FALSE)
print(available.station)

# Available stations for ECA&D data
available.station <- select.station(src="ECA&D")
print(available.station)

# Available stations for NACD data
available.station <- select.station(src="NACD")
print(available.station)

# Available stations for GHCN data
# Not Run !
available.station <- select.station(src="GHCN")
print(available.station)

# Available stations for Norway
available.station <- select.station(cntr="NORWAY")
print(available.station)

# Available stations for FRANCE
available.station <- select.station(cntr="FRANCE")
print(available.station)

# Available stations recording 2m-surface temperature
available.station <- select.station(param="t2m")
print(available.station)

# Available data sources for Oslo station
available.station <- select.station(stid=18700)
print(available.station)
# ....
}
