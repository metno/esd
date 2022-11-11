#' Developing a new frost query splitter/data fetcher here

# 1. Receive lists of sources, elements, from/totime
# 2. Form URL, check length. Too long? Halve the longest of the sources term and the elements term. Call recursively (once for each part).
# 3. URL not too long? Try call and check if error says asking for too much data. If so, halve the from/to and call recursively (once for each part).

## main function, this is the only one that should be exported
#' @export metno.frost.data
metno.frost.data <- function(keyfile='~/.FrostAPI.key', url='https://frost.met.no/auth/requestCredentials.html',
                             stid=NULL, param=NULL, it=NULL,
                             lon=NULL, lat=NULL, loc=NULL, alt=NULL, cntr=NULL,
                             timeresolutions='P1M', levels="default", timeoffsets="default", 
                             performancecategories="A,B,C", exposurecategories="1,2", 
                             qualities='0,1,2,3,4,5', fetch.meta=TRUE, path=NULL, 
                             browser="firefox", save2file=FALSE, verbose=FALSE) {
  
  if(verbose) print("Fetch data from the Frost API (http://frost.met.no)")

  ## Check requirements
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' needed to use 'metno.frost.station'. Please install it.")
  }
  if (is.null(param)) {
    stop("param must be defined")
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
  frostelements <- metno.frost.elements(elementid, timeresolutions)
  if (verbose) print( paste('Frost parameters = ', paste(frostelements, collapse=",")) )

  ## Get all station metadata for desired param and time resolution
  meta <- metno.frost.getmetadata(fetch.meta, param, timeresolutions)

  ## Get list of relevant stations / source IDs
  frostsources <- metno.frost.stations(stid, lat, lon, meta)
  if (verbose) print( paste("Frost stations:", paste(frostsources, collapse=",")) )

  ## Update or set start and end dates using metadata
  it <- metno.frost.dates(it, stid, elementid, meta)
  if (verbose) print( paste("Frost dates =", it[1], it[2]) )

  ## Fetch the data in as many API queries as frost.met.no requires
  data <- metno.frost.querysplitter(frostID, frostsources, frostelements, it,
                                    timeresolutions, levels, timeoffsets, performancecategories, 
                                    exposurecategories, qualities, verbose)

  ## TODO: create zoo object and station object
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
    qualities = qualities
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
      if (verbose) print("Splitting station list for frost.met.no query")
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
      if (verbose) print("Splitting param list for frost.met.no query")
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
    if (verbose) print("Splitting time interval for frost.met.no query")
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
    stop('No data in frost.met.no for these criteria')
  } else if (result$statuscode == 400) {
    stop('Something went wrong in esd (invalid frost.met.no API request)')
  } else if (result$statuscode == 500) {
    stop('Something went wrong in frost.met.no')
  } else if (result$statuscode == 429 | result$statuscode == 503) {
    stop('frost.met.no is under high load, try again later')
  }
  data
}

## Perform a single frost query and return parsed JSON results
metno.frost.dataquery <- function(frostID, parameters, verbose=FALSE) {
  ## Perform request
  endpoint <- "https://frost.met.no/observations/v0.jsonld"
  resp <- httr::GET(
    URLencode(endpoint),
    query = parameters,
    httr::authenticate(frostID[1], '')
  )

  ## Print URL if desired
  if (verbose) print(resp$url)

  # Handle any returned error
  statuscode <- httr::status_code(resp)
  table <- NULL
  if (statuscode == 200) {
    ## Parse the JSON
    xs <- httr::content(resp)
    marshalled <- jsonlite::toJSON(xs$data)
    table <- jsonlite::fromJSON(marshalled, flatten=TRUE)
    ## TODO: flatten in a better way without using tidyr?
  }
  ## Return both the table and the status code
  list(data=table, statuscode=statuscode)
}

# Ensure a keyfile exists, ask user to generate one if not
metno.frost.keyfile <- function(keyfile, verbose) {
  if (file.exists(keyfile)) {
    if (verbose) print(paste('Read client ID from',keyfile))
    frostID <- readLines(keyfile)
  } else { 
    if (verbose) print(paste('Generate new client ID',url))  
    system(paste(browser,url))
    frostID <- rep("",2)
    frostID[1] <- readline('Please give me the first key:')
    frostID[2] <- readline('Please give me the second key:')
    writeLines(frostID,con=keyfile)
  }
  frostID
}

## Return metadata data frame for station metadata
metno.frost.getmetadata <- function(fetch.meta, param, timeresolutions) {
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

## Return the CF name corresponding to the given esd element id
metno.frost.elements <- function(elementid, timeresolutions) {
  frostparaminfo <- ele2param(elementid, src="metno.frost")
  frostcfname <- gsub('*', timeresolutions, frostparaminfo$param, fixed=TRUE)
  frostelements <- c(frostcfname)
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
    stop('Found no stations with given criteria')
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
    stop('Found no stations with given criteria')
  }

  ## Return an updated it variable
  it <- c(start, end)
}
