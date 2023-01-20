
#' Retrieve records from frost.met.no, for one element and one or more stations
#' 
#' Fetch data from frost.met.no and return station object
#'
#' @aliases metno.frost.data metno.frost.station
#' 
#' @param stid A string of characters as an identifier of the weather/climate
#' station.
#' @param param Parameter or element type or variable identifier. There are
#' several core parameters or elements as well as a number of additional
#' parameters. The parameters or elements are: precip = Precipitation (mm) tas,
#' tavg = 2m-surface temperature (in degrees Celcius) tmax, tasmax = Maximum
#' temperature (in degrees Celcius) tmin, tasmin = Minimum temperature (in
#' degrees Celcius)
#' @param it A vector of two dates (from and to) in the form "2014-01-01".
#' @param \dots additional arguments  
#' @param keyfile The path where the frost.met.no credentials keyfile is stored locally.
#' @param url The URL of the data portal or webpage for requesting new client credentials. 
#' @param lon Numeric value of longitude (in decimal degrees East) for the
#' reference point (e.g. weather station) as a single value or a vector
#' containing the range of longitude values in the form of c(lon.min,lon.max)
#' @param lat Numeric value of latitude for the reference point (in decimal
#' degrees North) or a vector containing the range of latitude values in the
#' form of c(lat.min,lat.max)
#' @param loc A string of characters as the name of the location
#' (weather/climate station) or an object of class "stationmeta".
#' @param alt Numeric value of altitude (in meters a.s.l.) used for selection.
#' Positive value, select all stations above this altitude; for negative
#' values, select all stations below this latitude.
#' @param cntr A string or a vector of strings of the full name of the country:
#' Select the stations from a specified country or a set of countries.
#' @param timeresolutions A string of the desired time resolution. Accepted strings
#' are the ISO 8601 durations P1M, P1D and PT1M, as well as the strings MONTHLY,
#' MONTH, DAILY, DAY, MINUTE and MIN. If no value is set, this defaults to P1M.
#' @param levels Numeric value of sensor height used for selection. If no value
#' is set, then standard heights will be selected.
#' @param timeoffsets A string of the ISO 8601 duration desired for the offset
#' value of timeseries. If no value is set, then the best available option will 
#' be selected.
#' @param performancecategories A string with a comma-separated list of categories
#' to accept for sensor performance. If no value is set, then the default A,B,C
#' will be used.
#' @param exposurecategories A string with a comma-separated list of categories to
#' accept for siting quality of the sensor. If no value is set, then the default
#' 1,2 will be used.
#' @param qualities A string with a comma-separated list of quality codes to
#' accept for the individual measurements. If no value is set, then the default
#' 0,1,2,3,4,5 will be used.
#' @param fetch.meta Logical value defaulting to TRUE. If TRUE, download timeseries
#' metadata from frost.met.no. If FALSE, try to load data file station.meta instead.
#' @param path The path where the data are stored. Can be a symbolic link.
#' @param browser A string signifying which internet browser to open to ask
#' user to create new credentials for frost.met.no
#' @param save2file Logical value defaulting to FALSE. If TRUE, save a file
#' with the output data.
#' @param verbose Logical value defaulting to FALSE. If FALSE, do not display
#' comments (silent mode). If TRUE, displays extra information on progress.
#' @return A time series of "zoo" "station" class with additional attributes
#' used for further processing.
#'
#' @author K. Tunheim
#'
#' @keywords parameter,data,metno,norway,frost
#' 
#' @examples
#'
#' \dontrun{ 
#' metno.frost.data(param='t2m', stid=18700, it=c('2020-01-01','2020-02-01'), timeresolutions='P1D')
#' }
#'
#' @export metno.frost.data
metno.frost.data <- function(keyfile='~/.FrostAPI.key', url='https://frost.met.no/auth/requestCredentials.html',
                             stid=NULL, param=NULL, it=NULL,
                             lon=NULL, lat=NULL, loc=NULL, alt=NULL, cntr=NULL,
                             timeresolutions='P1M', levels="default", timeoffsets="default", 
                             performancecategories="A,B,C", exposurecategories="1,2", 
                             qualities='0,1,2,3,4,5', fetch.meta=TRUE, path=NULL, 
                             browser="firefox", save2file=FALSE, verbose=FALSE, ...) {
  
  if (verbose) print(match.call())

  ## Check requirements
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("metno.frost.data: Package 'jsonlite' needed to use 'metno.frost.station'. Please install it.")
  }
  if (is.null(param)) {
    stop("metno.frost.data: param must be defined")
  }

  ## Ensure there is a key file and get frost keys
  frostID <- metno.frost.keyfile(keyfile, verbose)

  ## Translate intuitive time resolution names into ISO 8601 duration
  timeresolutions <- switch(toupper(timeresolutions), 
                            "MONTHLY"="P1M", "MONTH"="P1M",
                            "DAILY"="P1D", "DAY"="P1D",
                            "MINUTE"="PT1M", "MIN"="PT1M",
                            timeresolutions)
  
  ## Get parameter metadata for desired param and time resolution
  elementid <- esd2ele(param)
  frostparaminfo <- ele2param(elementid, src="metno.frost")
  frostcfname <- gsub('*', timeresolutions, frostparaminfo$param, fixed=TRUE)
  if (verbose) print( paste('metno.frost.data: Frost parameters:', frostcfname) )

  ## Get all station metadata for desired param and time resolution
  meta <- metno.frost.getmetadata(fetch.meta, param, timeresolutions)

  ## Get list of relevant stations / source IDs
  frostsources <- metno.frost.stations(stid, lat, lon, meta)
  if (verbose) print( paste("metno.frost.data: Frost stations:", paste(frostsources, collapse=",")) )

  ## Update or set start and end dates using metadata
  it <- metno.frost.dates(it, stid, elementid, meta)
  if (verbose) print( paste0("metno.frost.data: Frost dates: ", it[1], '/', it[2]) )

  ## Fetch the data in as many API queries as frost.met.no requires
  data <- metno.frost.querysplitter(frostID, frostsources, c(frostcfname), it,
                                    timeresolutions, levels, timeoffsets, performancecategories, 
                                    exposurecategories, qualities, verbose)

  ## If empty result then return here
  if(is.null(data)) {
    return(invisible(NULL))
  }
  
  ## Rearrange data and transform into zoo and station object

  ## Get list of stations
  sourceId <- unique(data$sourceId)
  ## Extract data column
  var <- data[[3]]
  ## Parse time column into appropriate type
  if(timeresolutions=="PT1M") {
    time <- as.POSIXct(data$referenceTime)
  } else {
    time <- as.Date(data$referenceTime)
  }
  ## Get unique time steps
  tvec <- seq(
    min(time),
    max(time), 
    by = switch(timeresolutions, "P1D"="day", "P1M"="month", "PT1M"="min")
  )
  ## Set up matrix with a row per time step and a column per station
  X <- matrix(NA, nrow=length(tvec), ncol=length(sourceId))
  for(i in 1:ncol(X)) {
    j <- sapply(time[data$sourceId==sourceId[i]], function(x) which(tvec==x))
    X[j,i] <- var[data$sourceId==sourceId[i]]
  }
  ## Turn the matrix into a zoo object
  varzoo <- zoo(X, order.by=tvec)

  # Transform to station object and attach attributes
  stid <- gsub("[A-Z]|:.*","",toupper(sourceId))
  inds <- sapply(stid, function(x) which(meta$station_id==x)[1])
  METNO.FROST <- as.station(varzoo, stid=stid, loc=meta$location[inds],
                            param=param, quality=qualities, 
                            cntr=meta$country[inds],
                            lon=meta$longitude[inds], 
                            lat=meta$latitude[inds], 
                            alt=meta$altitude[inds],
                            src=switch(timeresolutions, 
                                       'PT1M'='METNO.FROST.MINUTE',
                                       'P1D'='METNOD.FROST', 
                                       'P1M'='METNOM.FROST'),
                            url="http://frost.met.no",
                            longname=frostparaminfo$longname,
                            unit=frostparaminfo$unit,
                            aspect="original",
                            reference="Frost API (http://frost.met.no)",
                            info="Frost API (http://frost.met.no)"
  )
  attr(METNO.FROST,'history') <- history.stamp(METNO.FROST)

  ## Save to file if requested
  if(save2file) {
    ## add rownames and colnames for titles
    rownames(X) <- as.character(tvec)
    colnames(X) <- substring(sourceId, 3, nchar(sourceId)-2)
    if (is.null(path)) path <- 'data.METNO'
    dir.create(path, showWarnings=FALSE, recursive=TRUE)
    ext <- switch(timeresolutions, 'PT1M'='obs', 'P1D'='dly', 'P1M'='mon')
    filename <- paste(attr(METNO.FROST,"source"),param,ext,sep=".")
    write.table(X, file=file.path(path,filename), col.names = NA, sep=",", quote = FALSE)
  }

  invisible(METNO.FROST)
}

## Split a data request into multiple frost API requests as needed
metno.frost.querysplitter <- function(frostID, frostsources=NULL, frostelements=NULL, it=NULL,
                                   timeresolutions='P1M', levels="default", timeoffsets="default", 
                                   performancecategories="A,B,C", exposurecategories="1,2", 
                                   qualities='0,1,2,3,4,5', verbose=FALSE) {
  
  ## Set up frost.met.no API query parameters
  sourcestring <- paste(paste0('SN', frostsources), collapse=',')
  elementstring <- paste(frostelements, collapse=',')
  parameters <- list(
    sources = sourcestring,
    elements = elementstring,
    referencetime = paste0(it[1], '/', it[2]),
    timeresolutions = timeresolutions,
    levels = levels,
    timeoffsets = timeoffsets,
    performancecategories = performancecategories,
    exposurecategories = exposurecategories,
    qualities = qualities,
    fields = "referenceTime,sourceId,elementId,value"
  )

  ## Attempt an initial request
  result <- metno.frost.dataquery(frostID, parameters, verbose)

  ## If status code is 200 (OK), then all is well; if not, split the query and do recursive calls
  data <- NULL
  if (result$statuscode == 200) {
    data <- result$data
  } else if (result$statuscode == 414) {
    ## Split longest list (frostsources or frostelements) and try again
    if (nchar(sourcestring) > nchar(elementstring)) {
      if (verbose) print("metno.frost.querysplitter: Splitting station list for frost.met.no query")
      halflength = floor( length(frostsources)/2 )
      sourcespart1 <- frostsources[1:halflength]
      sourcespart2 <- tail(frostsources, halflength)

      resultpart1 <- metno.frost.querysplitter(frostID, sourcespart1, frostelements, it, timeresolutions, 
                                            levels, timeoffsets, performancecategories, exposurecategories, 
                                            qualities, verbose)
      resultpart2 <- metno.frost.querysplitter(frostID, sourcespart2, frostelements, it, timeresolutions, 
                                            levels, timeoffsets, performancecategories,
                                            exposurecategories, qualities, verbose)
    } else {
      if (verbose) print("metno.frost.querysplitter: Splitting param list for frost.met.no query")
      halflength = floor( length(frostelements)/2 )
      elementpart1 <- frostelements[1:halflength]
      elementpart2 <- tail(frostelements, halflength)

      resultpart1 <- metno.frost.querysplitter(frostID, frostsources, elementpart1, it, timeresolutions, 
                                            levels, timeoffsets, performancecategories, exposurecategories, 
                                            qualities, verbose)
      resultpart2 <- metno.frost.querysplitter(frostID, frostsources, elementpart2, it, timeresolutions, 
                                            levels, timeoffsets, performancecategories, exposurecategories, 
                                            qualities, verbose)
    }
    ## Merge and return the two results
    data <- merge(resultpart1, resultpart2, all=TRUE)
  } else if (result$statuscode == 403) {
    if (verbose) print("metno.frost.querysplitter: Splitting time interval for frost.met.no query")
    t1 <- as.numeric(as.POSIXct(it[1], tz='UTC'))
    t2 <- as.numeric(as.POSIXct(it[2], tz='UTC'))
    tmid <- as.POSIXct((t1 + t2) / 2, origin = '1970-01-01', tz='UTC')
    itmid <- strftime(tmid, "%Y-%m-%dT%H:%M:%S%z", tz='UTC')
    itpart1 <- c(it[1], itmid)
    itpart2 <- c(itmid, it[2])

    resultpart1 <- metno.frost.querysplitter(frostID, frostsources, frostelements, itpart1, timeresolutions, 
                                          levels, timeoffsets, performancecategories, exposurecategories, 
                                          qualities, verbose)
    resultpart2 <- metno.frost.querysplitter(frostID, frostsources, frostelements, itpart2, timeresolutions, 
                                          levels, timeoffsets, performancecategories, exposurecategories, 
                                          qualities, verbose)
    ## Merge and return the two results
    data <- merge(resultpart1, resultpart2, all=TRUE)
  } else if (result$statuscode == 412) {
    ## TODO: this one should be accepted, with a legitimate NULL return that is handled upstream
    stop('metno.frost.querysplitter: No data in frost.met.no for these criteria')
  } else if (result$statuscode == 400) {
    stop('metno.frost.querysplitter: Something went wrong in esd (invalid frost.met.no API request)')
  } else if (result$statuscode == 500) {
    stop('metno.frost.querysplitter: Something went wrong in frost.met.no')
  } else if (result$statuscode == 429 | result$statuscode == 503) {
    stop('metno.frost.querysplitter: frost.met.no is under high load, try again later')
  }
  data
}

## Perform a single frost query and return parsed JSON results
metno.frost.dataquery <- function(frostID, parameters, verbose=FALSE) {
  ## Perform request
  endpoint <- "https://frost.met.no/observations/v0.csv"
  resp <- httr::GET(
    URLencode(endpoint),
    query = parameters,
    httr::authenticate(frostID[1], '')
  )
  statuscode <- httr::status_code(resp)

  ## Print URL if desired
  if (verbose) print(paste('metno.frost.dataquery:', resp$url))
  
  ## Parse the CSV (hopefully good enough)
  if (statuscode == 200) {
    xs <- httr::content(resp, encoding='UTF-8')
    table <- read.table(text=xs, sep=',', header=TRUE, check.names=FALSE)
    colnames(table) = sub("\\(-\\)", "", colnames(table))
  } else {
    table <- NULL
  }

  ## Return both the table and the status code
  list(data=table, statuscode=statuscode)
}

## Ensure a keyfile exists, ask user to generate one if not
metno.frost.keyfile <- function(keyfile, verbose) {
  if (file.exists(keyfile)) {
    if (verbose) print(paste('metno.frost.keyfile: Read client ID from',keyfile))
    frostID <- readLines(keyfile)
  } else { 
    if (verbose) print(paste('metno.frost.keyfile: Generate new client ID',url))  
    system(paste(browser,url))
    frostID <- rep("",2)
    frostID[1] <- readline('Please give me the first key:')
    frostID[2] <- readline('Please give me the second key:')
    writeLines(frostID,con=keyfile)
  }
  frostID
}

## Return metadata data frame for station metadata
metno.frost.getmetadata <- function(fetch.meta, param, timeresolutions, verbose=FALSE) {
  if(fetch.meta | timeresolutions=="PT1M") {
    meta.function <- switch(toupper(timeresolutions),
                            "PT1M"=metno.frost.meta.minute,
                            "P1D"=metno.frost.meta.day, 
                            "P1M"=metno.frost.meta.month)
    station.meta <- meta.function(param=c(param), save2file=FALSE, verbose=verbose)
  } else {
    data("station.meta", envir=environment())
  }
  id <- station.meta$source==switch(toupper(timeresolutions), 
                                    "PT1M"="METNO.FROST.MINUTE",
                                    "P1D"="METNOD.FROST",
                                    "P1M"="METNOM.FROST")
  meta <- station.meta[id,]
}

## Decide initial station list to use for frost query
metno.frost.stations <- function(stid, lat, lon, meta) {
  ## If stid is NULL, try to use (lat,lon) to find station(s) in metadata table
  if(is.null(stid) && (!is.null(lon) & !is.null(lat))) {
    if(length(lon)==1 & length(lat)==1) {
      d <- distAB(lon, lat, unlist(meta$lon), unlist(meta$lat))
      stid <- meta.stid[which.min(d)]
    } else {
      ok <- mapply(function(x,y) !is.na(x) & x<=max(lon) & x>=min(lon) & 
                     !is.na(y) & y<=max(lat) & y>=min(lat), meta$lon, meta$lat)
      stid <- unique(meta$station_id[ok])
    }
  }
  ## If stid is still NULL, fetch all data for given period
  if(is.null(stid)) {
    stid <- unique(meta$station_id)
  }
  stid
}

## Decide from and to dates to use for frost query
metno.frost.dates <- function(it, stid, elementid, meta) {
  ## Set initial start and end from it variable
  start <- NULL; end <- NULL
  if (!is.null(it)) {
    start <- it[1]; end <- it[2]
  }

  ## Select index for subset of stations that have the desired parameter
  i <- meta$station_id %in% stid & meta$element %in% elementid
  if(sum(i)==0) {
    stop('metno.frost.dates: Found no stations with given criteria')
  }
  
  ## Get earliest start and latest end of the relevant timeseries
  meta.start <- meta$start[i]
  meta.end <- meta$end[i]
  if(is.dates(meta.end)) {
    meta.end[is.na(meta.end)] <- strftime(Sys.time(), "%Y-%m-%d")
  } else {
    meta.end[!is.na(meta.end)] <- paste0(meta.end[!is.na(meta.end)],"-12-31")
    meta.end[is.na(meta.end)] <- strftime(Sys.time(), "%Y-%m-%d")
  }
  if(!is.dates(meta.start)) {
    meta.start <- paste0(meta.start,"-01-01")
  }
  if(is.null(start)) {
    start <- min(meta.start)
  } else if(!is.dates(start)) {
    start <- paste0(start,"-01-01")
  }
  if(is.null(end)) {
    end <- max(meta.end)
  } else if(!is.dates(end)) {
    end <- paste0(end,"-12-31")
  }

  ## Shrink date interval in variable 'it' if real start and end is smaller
  if(start<min(meta.start)) {
    start <- min(meta.start)
  }
  if(end>max(meta.end)) {
    end <- max(meta.end)
  }

  ## Are there any stations left?
  hasdata <- meta.start<=end & meta.end>=start
  ## Check if there are any stations left
  if(sum(hasdata)==0) {
    stop('metno.frost.dates: Found no stations with given criteria')
  }

  ## Return an updated it variable
  it <- c(start, end)
}


## Add alias metno.frost.station
metno.frost.station <- metno.frost.data
