#' Download METNO station metadata using frost.met.no
#' 
#' Where there are multiple measuring periods registered for the parameter,
#' only the earliest start time and the latest end time are used.
#' 
#' @aliases metno.frost.meta.day metno.frost.meta.month metno.frost.meta.minute
#' 
#' @param param Vector of parameters
#' @param save2file if TRUE, save metadata in a local file
#' @param verbose if TRUE, print diagnostics
#' @param \dots additional arguments  
#' @param url The URL to the webpage to request new client credentials at Frost API
#'
#' @return A meta data matrix object for all stations in METNO's collection
#' that have measured any of the given parameters. Start and end time are included. 
#'
#' @author K. Tunheim
#'
#' @keywords parameter,metadata,metno,norway,frost
#'
#' @examples
#'
#' \dontrun{
#' # Fetch metadata for all stations and measuring periods of the t2m parameter
#' metno.frost.meta.day(param=c('t2m'))
#' # Fetch metadata for all stations and measuring periods of all available parameters
#' metno.frost.meta.month()
#' }
#'
#' @export metno.frost.meta.day
metno.frost.meta.day <- function(param=c("t2m","precip","tmin","tmax","slp","pon","pox","fg","fx","pp","sd"), 
                                 save2file=FALSE, path=NULL, verbose=FALSE, ...) {
  if (verbose) print(match.call())
  X <- metno.frost.meta.default(param=param, timeresolutions="P1D", verbose=verbose, ...)
  filename <- "meta.metno.frost.day.rda"
  attr(X, "source") <- "METNOD.FROST"
  attr(X, "version") <- NA
  attr(X, "URL") <- "http://frost.met.no"
  attr(X, "file") <- filename
  attr(X, "cite") <- ""
  attr(X, "date") <- date()
  attr(X,"call") <- match.call()
  attr(X, "history") <- history.stamp(X)
  if (save2file) {
    meta.metno.frost.day <- X
    if(!is.null(path)) filename <- file.path(path,filename)
    ## KMP 2020-01-24: Added version because the rda files produced with R 3.5.1 (version 3) 
    ## is of a vectorized format that can't be read by earlier R version.
    save(meta.metno.frost.day, file=filename, version=2)
    rm("meta.metno.frost.day")
  }
  invisible(X)
}

# Get metadata for monthly timeseries - removed DD, DD06, DD12, DD18, SD
#' @export metno.frost.meta.month
metno.frost.meta.month <- function(param=c("t2m","precip","tmin","tmax","slp","pon","pox","fg","fx"), 
                                   save2file=FALSE, path=NULL, verbose=FALSE,...) {
  if (verbose) print(match.call())
  X <- metno.frost.meta.default(param=param, timeresolutions="P1M", verbose=verbose, ...)
  filename <- "meta.metno.frost.month.rda"
  attr(X, "source") <- "METNOM.FROST"  
  attr(X, "version") <- NA
  attr(X, "URL") <- "http://frost.met.no"
  attr(X, "file") <- "metno.frost.meta.month.rda"
  attr(X, "cite") <- ""
  attr(X, "date") <- date()
  attr(X, "call") <- match.call()
  attr(X, "history") <- history.stamp(X)
  if (save2file) {
    meta.metno.frost.month <- X
    if(!is.null(path)) filename <- file.path(path,filename)
    ## KMP 2020-01-24: Added version because the rda files produced with R 3.5.1 (version 3) 
    ## is of a vectorized format that can't be read by earlier R version.
    save(meta.metno.frost.month, file=filename, version=2)
    rm("meta.metno.frost.month")
  }
  invisible(X)
}

# Get metadata for minute timeseries
#' @export metno.frost.meta.minute
metno.frost.meta.minute <- function(param=c("t2m","precip","tmin","tmax","slp","pon","pox","fg","fx"), 
                                    save2file=FALSE, path=NULL, verbose=FALSE, ...) {
  if (verbose) print(match.call())
  X <- metno.frost.meta.default(param=param, timeresolutions="PT1M", verbose=verbose, ...)
  filename <- "meta.metno.frost.minute.rda"
  attr(X, "source") <- "METNO.FROST.MINUTE"
  attr(X, "version") <- NA
  attr(X, "URL") <- "http://frost.met.no"
  attr(X, "file") <- filename
  attr(X, "cite") <- ""
  attr(X, "date") <- date()
  attr(X,"call") <- match.call()
  attr(X, "history") <- history.stamp(X)
  if (save2file) {
    meta.metno.frost.min <- X
    if(!is.null(path)) filename <- file.path(path,filename)
    save(meta.metno.frost.min, file=filename, version=2)
    rm("meta.metno.frost.min")
  }
  invisible(X)
}


metno.frost.meta.default <- function(keyfile='~/.FrostAPI.key', param=c("t2m"), 
                                     timeresolutions="P1M", levels="default", timeoffsets="default", 
                                     performancecategories="A,B,C", exposurecategories="1,2", 
                                     url='https://frost.met.no/auth/requestCredentials.html',
                                     browser="firefox", verbose = FALSE) {
  if (verbose) print(match.call())
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("metno.frost.meta.default: Package 'jsonlite' needed, please install it.")
  } else {
    
    # KMP 2020-01-22: enable timeresolutions notation monthly and daily
    timeresolutions <- switch(toupper(timeresolutions), 
                              "MONTHLY"="P1M", "MONTH"="P1M",
                              "DAILY"="P1D", "DAY"="P1D",
                              "MINUTE"="PT1M", "MIN"="PT1M",
                              timeresolutions)
    
    # convert all param to local param names
    getparam1 <- function(x) {
      withstar <- ele2param(x, src="METNO.FROST")$param
      gsub('*', timeresolutions, withstar, fixed=TRUE)
    }
    ele <- sapply(param, esd2ele)
    param1s <- sapply(ele, getparam1)
    names(param1s) <- ele
    strparam <- paste0(param1s, collapse=",")
    if (verbose) print(paste('metno.frost.meta.default: params:', strparam))
    
    # Get a client_id using function from metno.frost.data.R
    frostID <- metno.frost.keyfile(keyfile, verbose)

    url1 <- paste0(
      "https://", 
      frostID[1],
      #client_id, 
      "@frost.met.no/",
      "sources/v0.jsonld",
      "?types=SensorSystem",
      "&country=Norge",
      "&validtime=0000-01-01/9999-01-01",
      "&fields=id,name,masl,country,county,countyId,municipality,municipalityId,geometry"
    )
    # KT 2020-05-25 - fetch all stations in and around Svalbard too
    url_sj <- paste0(
      "https://",
      frostID[1],
      "@frost.met.no/",
      "sources/v0.jsonld",
      "?types=SensorSystem",
      "&country=Svalbard og Jan Mayen",
      "&validtime=0000-01-01/9999-01-01",
      "&fields=id,name,masl,country,county,countyId,municipality,municipalityId,geometry"
    )
    url2 <- paste0(
      "https://",
      frostID[1],
      "@frost.met.no/",
      "observations/availableTimeSeries/v0.jsonld",
      "?elements=", strparam,
      "&timeresolutions=", timeresolutions,
      "&levels=", levels,
      "&timeoffsets=", timeoffsets,
      "&performancecategories=", performancecategories,
      "&exposurecategories=", exposurecategories,
      "&fields=sourceId,elementId,validFrom,validTo"
    )
    if (verbose) {
      print(paste('metno.frost.meta.default:', url1))
      print(paste('metno.frost.meta.default:', url_sj))
      print(paste('metno.frost.meta.default:', url2))
    }
    # KT 2020-05-26: getting data from both Norge and Svalbard and Jan Mayen
    xs_no <- jsonlite::fromJSON(URLencode(url1), flatten=TRUE)
    xs_sj <- jsonlite::fromJSON(URLencode(url_sj), flatten=TRUE)
    xs1 <- rbind(xs_no$data, xs_sj$data)
    xs1$lon = sapply(xs1$geometry.coordinates, function(x) x[1])
    xs1$lon[sapply(xs1$lon, is.null)] <- NA
    xs1$lon <- unlist(xs1$lon)
    xs1$lat = sapply(xs1$geometry.coordinates, function(x) x[2])
    xs1$lat[sapply(xs1$lat, is.null)] <- NA
    xs1$lat <- unlist(xs1$lat)
    df1 <- xs1[c("id","name","country","lon","lat","masl","municipality",
                 "municipalityId","county","countyId")]

    xs2 <- jsonlite::fromJSON(URLencode(url2), flatten=TRUE)
    df2 <- xs2$data
    df2$sourceId = substring(df2$sourceId, 1, nchar(df2$sourceId)-2)
    df <- data.frame(NULL)
    for (i in which(param1s %in% df2$elementId)) {
      # KT 2020-05-26: preserve NA as latest validTo
      dfparam = df2[df2$elementId == param1s[i], ]
      dfparam$validTo[is.na(dfparam$validTo)] = "9999-12-31T00:00:00.000Z"
      validFrom = try(aggregate(validFrom ~ sourceId, data=dfparam, min), silent=TRUE)
      validTo = try(aggregate(validTo ~ sourceId, data=dfparam, max), silent=TRUE)
      validTo[validTo=="9999-12-31T00:00:00.000Z"] <- NA
      if (class(validFrom) != "try-error" & length(validFrom) > 0) {
        validFrom$validFrom <- as.Date(validFrom$validFrom)
        validTo$validTo <- as.Date(validTo$validTo)

        period = merge(validFrom, validTo, by='sourceId', all.x=TRUE)
        stperiod = merge(df1, period, by.x="id", by.y="sourceId")

        colnames(stperiod) = c("station_id","location","country","lon","lat","altitude",
                               "municipality","municipalityid","county","countyid","start","end")
        stperiod$element <- rep(names(param1s[i]),
                                      length(stperiod$station_id))

        # convert to UTM
        utmZone <- 33
        XY <- LatLon2UTM(lat=stperiod$lat, lon=stperiod$lon, zone=utmZone)
        stperiod$utm_east  <- XY[[1]]
        stperiod$utm_north <- XY[[2]]
        stperiod$utm_zone  <- rep(utmZone, length(stperiod$station_id))

        df <- rbind(stperiod, df, stringsAsFactors=FALSE)
      }
    }
    
    #invisible(df)
    ## Same format as station.meta
    var <- df$element
    for(element in unique(df$element)) {
      var[df$element==element] <- esd2ele(element)
    }
    cntr <- sapply(df$country, function(x) switch(x, "Norge"="NORWAY", x))
    X <- data.frame("station_id"=gsub("[^0-9]","",df$station_id),
                    "location"=df$location,
                    "country"=cntr,
                    "longitude"=df$lon,
                    "latitude"=df$lat,
                    "altitude"=df$altitude,
                    "element"=df$element,
                    "start"=strftime(df$start, format="%Y"),
                    "end"=strftime(df$end, format="%Y"),
                    "source"=switch(timeresolutions,
                                    "P1D"="METNOD.FROST",
                                    "P1M"="METNOM.FROST",
                                    "PT1M"="METNO.FROST.MINUTE"),
                    "wmo"=rep(NA,length(df$station_id)),
                    "quality"=rep(NA,length(df$station_id)),
                    "variable"=var, stringsAsFactors=FALSE)
    attr(X,"metnoURLs") <- "http://frost.met.no"
    attr(X,"author") <- "K. Tunheim & K. Parding"
    attr(X,"date") <- Sys.time()
    attr(X,"history") <- history.stamp(X)
    class(X) <- c("stationmeta", class(X))
    invisible(X)
  }
}