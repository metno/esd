# There is another version in select.R. Which one is correct?
select.station <- function (x=NULL,loc=NULL , param = NULL,  ele = NULL,
                            stid = NULL ,lon = NULL, lat = NULL, 
                            alt = NULL, cntr = NULL, ..., src = NULL ,
                            it = NULL , nmin = NULL , verbose=FALSE) {
  if (verbose) print('select.station')
  if (is.null(x)) {
    data("station.meta",envir=environment())
    station.meta$end[is.na(station.meta$end)] <- format(Sys.time(),'%Y')
    station.meta <- as.data.frame(station.meta,stringsAsFactors=FALSE)
    if (!is.null(param)) ele <- apply(as.matrix(param),1,esd2ele) 
  } else {
    if (inherits(x,"station")) {      
      station_id <- attr(x, "station_id")
      location <- attr(x,"location")
      country <- attr(x,"country")
      longitude <- attr(x,"longitude")
      latitude <- attr(x,"latitude")
      altitude <- attr(x,"altitude")
      element <- apply(as.matrix(attr(x,"variable")),1,esd2ele)
      start <- rep(year(x)[1],length(station_id))
      end <- rep(year(x)[length(index(x))],length(station_id))
      source <- attr(x,"source")
      quality <- attr(x,"quality")

      station.meta <- data.frame(station_id = station_id,location = location, country = country,
                                 longitude = longitude, latitude = latitude,altitude = altitude,
                                 element = element, start = start,end = end,source = source,quality= quality)

      # update ele using element
      ## ele <- element
      ## param <- esd2ele(ele)
    } else stop("x must be an object of calss 'station'") 
  }
  ##
  if (!is.null(param) & is.null(ele)) {
    print("No variable found for your selection or the param identifier has not been set correctly.")
    print("Please refresh your selection based on the list below")
    print(as.matrix(ele2param(src=src))[,c(2,5,6)])
  }  
  ## get the lenght of the data base
  n <- length(station.meta$station_id)
  ## Search by station identifier
  if (!is.null(stid)) {
    if (is.numeric(stid)) {
      id <- is.element(station.meta$station_id,stid)
      station.meta <- station.meta[id,]
    } else if (is.character(stid)) {
      id <- is.element(station.meta$station_id,stid) # grep(stid,station.meta$station_id)
      station.meta <- station.meta[id,]
    }
  }
  ## Search by the closest station to longitude and latitude values
  if (length(lon)==1 & length(lat)==1) {
    ## AM 25.09.2013 STILL NEED TO TEST
    d <- distAB(lon, lat, station.meta$longitude, station.meta$latitude)
    id <- d==min(d,na.rm=TRUE)
    ##id[is.na(id)] <- FALSE # Ak some of the lon values are NA's
    station.meta <- station.meta[id,]
  } else {## Search by longitude values or range of values
    if (!is.null(lon)) {##search by longitude values or a range of values
      lon.rng <- range(lon,na.rm=TRUE) 
      id <- (station.meta$longitude >= lon.rng[1]) & (station.meta$longitude <= lon.rng[2])
      station.meta <- station.meta[id,]  
    }
    ## Search by latitude values or range of values
    if (!is.null(lat)) {
      lat.rng <- range(lat) 
      id <- (station.meta$latitude >= lat.rng[1]) & (station.meta$latitude <= lat.rng[2])
      station.meta <- station.meta[id,]
    }
  }
  ## Search by altitude values or range of values
  if (!is.null(alt)) {
    if (length(alt) == 1) {
      if (alt > 0) alt.rng <- c(alt,10000)
      else alt.rng <- c(0,abs(alt))
    }  else if (length(alt) == 2) {
      alt.rng <- alt      
    }
    id <- (station.meta$altitude >= alt.rng[1]) & (station.meta$altitude <= alt.rng[2])
    station.meta <- station.meta[id,]   
  }
  ## Search by country name
  if (!is.null(cntr)) {
    id <- is.element(tolower(station.meta$country),tolower(cntr))
    station.meta <- station.meta[id,]   
  }
  ##
  ## Search by data source
  if (!is.null(src)) {
    id <- is.element(tolower(station.meta$source),tolower(src))
    station.meta <- station.meta[id,]
  }
    ## 
    ## Search by location
  if (!is.null(loc)) {
    ## id <- is.element(tolower(station.meta$location),tolower(loc))
      pattern <- paste(loc,collapse='|')
      id <- grep(pattern=pattern,station.meta$location,ignore.case=TRUE,...)
      station.meta <- station.meta[id,]
  }
    ##
    ## Search by starting and ending years
    if (!is.null(it)) { 
        it.rng <- range(as.numeric(it),na.rm=TRUE)
        id <- (as.numeric(station.meta$start) <= it.rng[1]) & (as.numeric(station.meta$end) >= it.rng[2])
        ##
        if (sum(id,na.rm=TRUE)==0) {
            print(paste('No records that cover the period ',it.rng[1],'-',it.rng[2],'. Earliest observation from ',
                        min(as.numeric(station.meta$start)),' and latest observation from ',
                        max(as.numeric(station.meta$end)),sep=''))
        }
        station.meta <- station.meta[id,]
        station.meta$start <- rep(it.rng[1],length(station.meta$loc))
        ## paste('01-01-',rep(it.rng[1],length(station.meta$loc)),sep='')
        station.meta$end <- rep(it.rng[2],length(station.meta$loc))
        ## paste('31-12-',rep(it.rng[2],length(station.meta$loc)),sep='')
    }

  ## Search by minimum number of years
  if (!is.null(nmin)) { 
    ny <- as.numeric(station.meta$end) - as.numeric(station.meta$start) + 1
    id <- (ny >= nmin)
    station.meta <- station.meta[id,]
  } 
  ## Search by esd element
  if (!is.null(ele)) {
    id <- is.element(station.meta$element,ele)
    station.meta <- station.meta[id,]
  }
  ## Outputs
  if (dim(station.meta)[1]!=0) {
    station.meta$station_id <- as.character(station.meta$station_id)
    station.meta$location <- as.character(station.meta$location)
    station.meta$country <- as.character(station.meta$country)
    station.meta$source <- as.character(station.meta$source)
    class(station.meta) <- c("stationmeta","data.frame")
    return(station.meta)
  } else {
    print("No available stations found for your selection")
    return(NULL)
  }
}

## test.select.station performs a series of tests and print results for visual checks
# NOT EXPORTED
test.select.station <- function() {
  available.station <- select.station(param="t2m",lon = c(0,10),lat=c(20,30),src="ECAD",alt=c(500,1000),verbose=FALSE)
  summary(available.station)

  available.station <- select.station(param="t2m",src="ECAD")
  str(available.station)

  available.station <- select.station(src="NACD")
  str(available.station)
  map(available.station)
  
  available.station <- select.station(param="t2m",src="GHCND")
  str(available.station)
  map(available.station)

  available.station <- select.station(cntr="NORWAY")
  str(available.station)

  available.station <- select.station(param="t2m")
  map(available.station)

  available.station <- select.station(param="precip",lon=c(0,30),lat=c(50,70))
  str(available.station)
  map(available.station)

  available.station <- select.station(loc="oslo")
  str(available.station)
  map(available.station)

  available.station <- select.station(loc="oslo")
  src <- rownames(table(available.station$source))
  print(src)
}
