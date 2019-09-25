## This downloads a full list of METNO station metadata
## relevant to the specified param, using the data API frost.met.no.
##
## Where there are multiple measuring periods registered for the parameter,
## only the earliest start time and the latest end time are used.
##
## Author: K. Tunheim

# source("~/esd/R/dictionary.R")
library(jsonlite) # remove this line in final version?

# get diurnal timeseries - removed DD, DD06, DD12, DD18, SD
metno.frost.meta.diurnal <- function(param=c("t2m","precip","tmin","tmax","slp","pon","pox","fg","fx"), save=TRUE,...) {
  X <- metno.frost.meta.default(param=param, timeresolutions="P1D", ...)

  attr(X, "source") <- "METNO.FROST.DIURNAL"  
  attr(X, "version") <- NA
  attr(X, "URL") <- "http://frost.met.no"
  attr(X, "file") <- "metno.frost.meta.diurnal.rda"
  attr(X, "cite") <- ""
  attr(X, "date") <- date()
  attr(X,"call") <- match.call()

  if (save) {
    metno.frost.meta.diurnal <- X
    save(metno.frost.meta.diurnal, file="metno.frost.meta.diurnal.rda")
    rm("metno.frost.meta.diurnal")
  }

  invisible(X)
}

# get monthly timeseries - removed DD, DD06, DD12, DD18, SD
metno.frost.meta.month <- function(param=c("t2m","precip","tmin","tmax","slp","pon","pox","fg","fx"), save=TRUE,...) {
  X <- metno.frost.meta.default(param=param, timeresolutions="P1M", ...)

  attr(X, "source") <- "METNO.FROST.MONTH"  
  attr(X, "version") <- NA
  attr(X, "URL") <- "http://frost.met.no"
  attr(X, "file") <- "metno.frost.meta.month.rda"
  attr(X, "cite") <- ""
  attr(X, "date") <- date()
  attr(X,"call") <- match.call()

  if (save) {
    metno.frost.meta.month <- X
    save(metno.frost.meta.month, file="metno.frost.meta.month.rda")
    rm("metno.frost.meta.month")
  }

  invisible(X)
}

metno.frost.meta.default <- function(param=c("t2m"), timeresolutions="P1M", levels="default", timeoffsets="default", 
                        performancecategories="A,B,C", exposurecategories="1,2", verbose = FALSE) {
  # TODO: get a client_id
  client_id <- '0763dab1-d398-4a56-ba5d-601d7d352999'

  # convert all param to local param names
  getparam1 <- function(x) {
    withstar <- ele2param(x, src="metno.frost")$param
    gsub('*', timeresolutions, withstar, fixed=TRUE)
  }
  ele <- sapply(param, esd2ele)
  param1s <- sapply(ele, getparam1)
  names(param1s) <- ele
  strparam <- paste0(param1s, collapse=",")

  if (verbose) print(strparam)

  url1 <- paste0(
    "https://", client_id, "@frost.met.no/",
    "sources/v0.jsonld",
    "?types=SensorSystem",
    "&country=Norge",
    "&fields=id,name,masl,country,county,countyId,municipality,municipalityId,geometry"
  )
  url2 <- paste0(
    "https://", client_id, "@frost.met.no/",
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
    print(url1)
    print(url2)
  }

  xs1 <- jsonlite::fromJSON(URLencode(url1), flatten=T)
  xs1$data$lon = sapply(xs1$data$geometry.coordinates, function(x) x[1])
  xs1$data$lat = sapply(xs1$data$geometry.coordinates, function(x) x[2])
  df1 <- xs1$data[c("id","name","country","lon","lat","masl","municipality","municipalityId","county","countyId")]

  xs2 <- jsonlite::fromJSON(URLencode(url2), flatten=T)
  df2 <- xs2$data
  df2$sourceId = substring(df2$sourceId, 1, nchar(df2$sourceId)-2)

  df <- data.frame(NULL)
  for (i in 1:length(param1s)) {
    dfparam = df2[df2$elementId == param1s[i], ]
    validFrom = try(aggregate(validFrom ~ sourceId, data=dfparam, min), silent=TRUE)
    validTo = try(aggregate(validTo ~ sourceId, data=dfparam, max), silent=TRUE)

    if (class(validFrom) != "try-error" & length(validFrom) > 0) {
      validFrom$validFrom <- as.Date(validFrom$validFrom)
      validTo$validTo <- as.Date(validTo$validTo)

      period = merge(validFrom, validTo, by='sourceId', all.x=TRUE)
      stperiod = merge(df1, period, by.x="id", by.y="sourceId")

      colnames(stperiod) = c("station_id","location","country","lon","lat","altitude",
                             "municipality","municipalityid","county","countyid","start","end")

      stperiod$element <- rep(names(param1s[i]),length(stperiod$station_id))

      # convert to UTM
      utmZone <- 33
      XY <- LatLon2UTM(lat=stperiod$lat, lon=stperiod$lon, zone=utmZone)
      stperiod$utm_east  <- XY[[1]]
      stperiod$utm_north <- XY[[2]]
      stperiod$utm_zone  <- rep(utmZone, length(stperiod$station_id))

      df <- rbind(stperiod, df, stringsAsFactors=FALSE)
    }
  }
  
  invisible(df)
}
