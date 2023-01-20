
## This function is used to check wether there are errors in the programming
# NOT EXPORTED
test.station <- function(ss=NULL,stid=NULL,alt=NULL,lat=c(50,70),lon=c(0,30),param="precip",
                         src=c("GHCND","GHCNM","NORDKLIM","NACD","METNOM","METNOD","ECAD"),verbose=FALSE) {
  for (i in 1:length(src)) {
    y <- station(stid=stid,alt=alt,lat=lat,lon=lon,param=param,src=src[i],nmin=100,verbose=verbose)
    print(summary(y))
  }
}

#' Retrieve station record from a given data source.
#' 
#' \code{allgood} and \code{clean.station} provide two filters for extracting
#' stations with good data (discarding missing values). \code{allgood} will not
#' leave any NA's whereas \code{clean.station} provides a more 'gentle'
#' filtering.
#' 
#' \code{station.sonel}, \code{station.gloss}, and \code{station.newlyn} read
#' sea level from tidal gauges in France (SONEL), on a global scale (GLOSS) and
#' for a single station (sub-daily data)in the UK (Newlyn).
#' 
#' @aliases station station.default station.ecad station.nacd
#' station.narp station.nordklim station.metnod station.metnom station.ghcnd
#' station.ghcnm station.ghcnm station.sonel station.gloss station.newlyn
#' station.giss
#'
#' @seealso clean.station allgood station.thredds map.station select.station
#'
#' @param loc A string of characters as the name of the location
#' (weather/climate station) or an object of class "stationmeta".
#' @param param Parameter or element type or variable identifier. There are
#' several core parameters or elements as well as a number of additional
#' parameters. The parameters or elements are: precip = Precipitation (mm) tas,
#' tavg = 2m-surface temperature (in degrees Celcius) tmax, tasmax = Maximum
#' temperature (in degrees Celcius) tmin, tasmin = Minimum temperature (in
#' degrees Celcius)
#' @param src Source: limit the downscaling to a specific data set ("NARP",
#' "NACD", "NORDKLIMA", "GHCNM", "METNOM", "ECAD", "GHCND" and "METNOD")
#' @param stid A string of characters as an identifier of the weather/climate
#' station.
#' @param lon Numeric value of longitude (in decimal degrees East) for the
#' reference point (e.g. weather station) as a single value or a vector
#' containing the range of longitude values in the form of c(lon.min,lon.max)
#' @param lat Numeric value of latitude for the reference point (in decimal
#' degrees North) or a vector containing the range of latitude values in the
#' form of c(lat.min,lat.max)
#' @param alt Numeric value of altitude (in meters a.s.l.) used for selection.
#' Positive value, select all stations above this altitude; for negative
#' values, select all stations below this latitude.
#' @param cntr A string or a vector of strings of the full name of the country:
#' Select the stations from a specified country or a set of countries.
#' @param it A single integer or a vector of integers or Dates. An integer in
#' the range of [1:12] for months, an integer of 4 digits for years (e.g.
#' 2014), or a vector of Dates in the form "2014-01-01").
#' @param is Index space: integer (station ID) or character (location name) or
#' list with lon/lat ranges.
#' @param nmin Select only stations with at least nmin number of years, months
#' or days depending on the class of object x (e.g. 30 years).
#' @param plot Logical value. If, TRUE provides a plot.
#' @param verbose Logical value defaulting to FALSE. If FALSE, do not display
#' comments (silent mode). If TRUE, displays extra information on progress.
#' @param path The path where the data are stored. Can be a symbolic link.
#' @return A time series of "zoo" "station" class with additional attributes
#' used for further processing.
#' @param url The URL of the data portal or webpage for requesting new client credentials. 
#'
#' @author A. Mezghani
#'
#' @keywords select.station
#' @examples
#' 
#'  \dontrun{
#' # Get daily and monthly mean temperature for "Oslo" station ("18700") from METNO data source
#' t2m.dly <- station.metnod(stid="18700",param="t2m")
#' t2m.mon <- station.metnom(stid="18700",param="t2m")
#' 
#' # Get daily data from the ECA&D data source:
#' # If called for the first time, the script will download a huge chunk of
#' # data and store it locally.
#' # select meta for "De Bilt" station into ss,
#' ss <- select.station(loc = "de bilt",param="t2m",src="ECAD")
#' # Retrieve the data from the local directory specified in path based on
#' # previous selected station 
#' t2m.dly <- station.ecad(loc=ss,path="~/data.ECAD")
#' # or directly retrieve the data without a prior selection 
#' t2m.dly <- station.ecad(loc = "oslo - blindern",param="t2m",path="~/data.ECAD")
#' plot(t2m.dly)
#' # Aggregate to monthly and annual mean temperature values and plot the results
#' t2m.mon <- as.monthly(t2m.dly, FUN="mean") ; plot(t2m.mon)
#' t2m.ann <- as.annual(t2m.mon, FUN = "mean") ; plot(t2m.ann)
#' # specify one station from ECAD, and this time get daily mean precipitation
#' precip.dly <- station.ecad(loc="Oxford",param="precip") ; plot(precip.dly)
#' # Aggregate to annual accumulated precipitation values and plot the result
#' precip.ann <- as.annual(precip.dly,FUN="sum") ; plot(precip.ann)
#' 
#' # Get daily data from the GHCND data source
#' # Select a subset of stations across Norway with a minimum number of
#' # 130 years using "GHCND" as a data source, retrieve the data and show its
#' # structure.
#' ss <- select.station(cntr="NORWAY",param="precip",src="GHCND",nmin=130)
#' y <- station.ghcnd(loc=ss , path="~/data.GHCND",plot=TRUE)
#' str(y)
#' # Subselect one station and display the geographical location of both selected 
#' # stations and highlight the subselected station (is=2).
#' y1 <- subset(y,is=2)
#' map(y, xlim = c(-10,30), ylim = c(50,70), cex=1, select=y1, cex.select=2, showall=TRUE)
#' }
#'
#' @importFrom utils download.file head read.csv read.csv2 read.fwf read.table write.table untar unzip
#'
#' @export station
station <- function(...) UseMethod("station")

#' @exportS3Method
#' @export station.ecad
station.ecad <- function(...) {
  y <- station(src="ecad",...)
  invisible(y)
}

#' @exportS3Method
#' @export station.ghcnd
station.ghcnd <- function(...) {
  y <- station(src="ghcnd",...)
  invisible(y)
}

#' @exportS3Method
#' @export station.nacd
station.nacd <- function(...) {
  y <- station(src="nacd",...)
  invisible(y)
}

#' @exportS3Method
#' @export station.narp
station.narp <- function(...) {
  y <- station(src="narp",...)
  invisible(y)
}

#' @exportS3Method
#' @export station.nordklim
station.nordklim <- function(...) {
  ## 
  y <- station(src="nordklim",...)
  invisible(y)
}

#' @exportS3Method
#' @export station.metnom
station.metnom <- function(...) {
  y <- station(src="metnom",...)
  invisible(y)
}

#' @exportS3Method
#' @export station.metnod
station.metnod <- function(...) {
  y <- station(src="metnod",...)
  invisible(y)
}

#' @exportS3Method
#' @export station.ghcnm
station.ghcnm <- function(...) {
  y <- station(src="ghcnm",...)
  invisible(y)
}

# NOT EXPORTED
station.daily <- function(...,src=NULL) {
  
  SRC <- c("METNOD","ECAD","GHCND")  
  if (is.null(src)) src <- SRC else src <- intersect(tolower(src),tolower(SRC))
  y <- station(src=src,...)
  attr(y,"call") <- match.call()
  invisible(y)
  
}

# NOT EXPORTED
station.monthly <- function(...,src=NULL) {
  
  SRC <- c("NACD","NARP","NORDKLIM","METNOM","GHCNM")
  if (is.null(src)) src <- SRC else src <- intersect(src,SRC)
  y <- station(src=src,...)
  attr(y,"call") <- match.call()
  invisible(y)
  
}

#' @exportS3Method
#' @export station.default
station.default <- function(..., loc=NULL, param='t2m', src=NULL, path=NULL,
                            qual=NULL, url = NULL, stid=NULL, lon=NULL, lat=NULL, alt=NULL,
                            cntr=NULL, it= NULL,nmin=NULL, plot=FALSE, verbose=FALSE,
                            path.ecad=NULL,url.ecad=NULL,
                            path.ghcnm=NULL,url.ghcnm=NULL,
                            path.ghcnd=NULL,url.ghcnd=NULL,
                            path.metnom=NULL,url.metnom=NULL,
                            path.metnod=NULL,url.metnod=NULL,
                            user='external',save2file=FALSE) {
  
  if (verbose) print(match.call())
  ## REB 2020-05-26: a hack - the old code took forever...
  if (!is.null(src)) if (src=='metnod.thredds') {
    y <- station.thredds(locs=loc,param=param,stid=stid,lon=lon,lat=lat,
                         alt=alt,cntr=cntr,it=it,nmin=nmin,verbse=verbose)
    return(y)
  }
  ## end of hack meant to streamline the use REB...
  ## REB 2021-05-11
  l <- list(...)
  dots <- names(l)
  if (length(l) > 0) { 
    if (verbose) print("station.default: Taking '...' to select station (station.meta/data.frame)")
    ic <- match('ss',names(l))
    if (!is.na(ic)) ss <- l$ss else ss <- l[[1]] 
  } else {
    if (verbose) print("station.default: ss <- NULL")
    ss <- NULL
  }
  if (inherits(loc,"stationmeta")) {
    ss <- loc
  } else if (is.character(loc)) {
    loc <- loc
    ss <- NULL
  } else if(length(list(...))>0) {
    if(inherits(list(...)[[1]],"stationmeta")) {
      ss <- list(...)[[1]]
    }
  }
  
  if (is.null(ss)) {
    ss <- select.station(stid=stid,loc=loc,lon=lon,lat=lat,alt=alt,cntr=cntr,
                         param=param,src=src,it=it,nmin=nmin,user=user,
                         verbose=verbose)
  }
  
  if ((param=="t2m") & is.null(ss)) {
    param0 <- param
    ssn <- select.station(param="tmin",stid=stid,loc=loc,lon=lon,lat=lat,alt=alt,
                          cntr=cntr,src=src,it=it,nmin=nmin,user=user,verbose=verbose)
    ssx <- select.station(param="tmax",stid=stid,loc=loc,lon=lon,lat=lat,alt=alt,
                          cntr=cntr,src=src,it=it,nmin=nmin,user=user,verbose=verbose)
    if (!is.null(ssn) & !is.null(ssx)) {
      class(ssn) <- class(ssx) <- "data.frame"
    } else {
      print('station.default: Found no stations with given criteria')
      return(NULL)
    }
    ss <- subset(ssx,ssx$station_id==ssn$station_id) # keep only stations recording both min and max
    if (is.null(ss)) {
      return(NULL)
    } else {
      rl <- readline(paste0("station.default: T2m is not available for your selection but TMIN and TMAX have been found",
                            " - Would you like to continue using the averaged values? (y or n): "))
      if ((rl=="y") | rl==("ye") | (rl=="yes")) {
        ss$element <- rep(esd2ele(param),length(ss$station_id)) ## update element with param
      } else {
        stop("station.default: Process stopped")
      }
    }
  } 

  if (verbose) {
    print("station.default: Station ID:")
    str(ss$station_id)
  }

  X <- NULL
  src <- as.character(ss$source)
  sources <- unique(src)

  ## Add .FROST to METNOM and METNOD source names for backwards compatibility – unless internal user!
  if(any(grepl("METNO",sources))) {
    if(user!="metno") {
      sources <- unique(sapply(sources, function(x) {
        switch(toupper(x), "METNOM"="METNOM.FROST", "METNOD"="METNOD.FROST", x)}))
    }
  }

  ## Loop through requested data sources
  for(s in sources) {
    if(verbose) print(paste("station.default: Retrieving data from source",s))

    ## Set stid, param and a default retrieval path
    stid <- ss$station_id[src==s]
    param <- apply(as.matrix(ss$element),1,esd2ele)[src==s]
    if(!is.null(it)) {
      start <- min(it)
      end <- max(it)
    } else {
      args <- list(...)
      if("start" %in% names(args)) start <- args$start else start <- NULL
      if("end" %in% names(args)) end <- args$end else end <- NULL
    }
    path <- paste0("data.",toupper(s))

    ## Set specific path and other defaults depending on source
    if(grepl("METNOM.FROST",toupper(s))) {
      if(!is.null(path.metnom)) path <- path.metnom
      timeres <- "P1M"
    } else if(grepl("METNO.FROST.MINUTE",toupper(s))) {
      timeres <- "PT1M"
    } else if(grepl("METNOD.FROST",toupper(s))) {
      if(!is.null(path.metnod)) path <- path.metnod
      timeres <- "P1D"
    } else if(grepl("METNOD.THREDDS",toupper(s))) {
      if(!is.null(path.metnod)) path <- path.metnod
    } else if(s=="ECAD") {
      if (!is.null(path.ecad)) path <- path.ecad
      if (is.null(url.ecad)) {
        url="https://knmi-ecad-assets-prd.s3.amazonaws.com/download/ECA_nonblend"
        #url="http://www.ecad.eu/utils/downloadfile.php?file=download/ECA_nonblend"
      } else url <- url.ecad
    } else if(s=="GHCNM") {
      if(!is.null(path.ghcnm)) path <- path.ghcnm
      if(is.null(url.ghcnm)) url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn" else url <- url.ghcnm
    } else if(s=="GHCND") {
      ## https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/
      if(!is.null(path.ghcnd)) path <- path.ghcnd
      if(is.null(url.ghcnd)) url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all" else url <- url.ghcnd
    }

    ## Loop through requested parameters
    for(param0 in unique(param)) {
      if(grepl("FROST",toupper(s))) {
        ## Special treatment for FROST, because multiple stations should be fetched in one query for efficiency
        ## (TODO: If X can contain multiple parameters, metno.frost.station can be modified to allow this as well)
        x <- metno.frost.station(timeresolutions=timeres, stid=stid, 
                                 param=param0, verbose=verbose, path=path, 
                                 it=c(start, end), save2file=save2file)
        if(!is.null(x)) if(is.null(X)) X <- x else X <- combine.station(X,x)
      } else {
        ## Common variables for all source types
        j <- ss$source==s & ss$variable==param0
        stid <- ss$station_id[j]
        loc <- ss$location[j]
        lon <- ss$longitude[j]
        lat <- ss$latitude[j]
        alt <- ss$altitude[j]
        cntr <- ss$country[j]
        ele <- ss$element[j]
        qual <- ss$quality[j]
        start <- ss$start[j]
        end <- ss$end[j]
        if(verbose) print(paste("station.default: Retrieving data from",
                                length(stid), "records ..."))

        ## Fetch data from one station at a time
        for (i in 1:length(stid)) {
          if(verbose) print(paste('station.default:',i,toupper(param0),stid[i],loc[i],cntr[i],s))

          ## Custom treatment for the different source types
          if (grepl("METNOD",toupper(s))) {
            ## NOTE: this block will only trigger for internal users
            if (param0!='dd') {
              param1 <- esd2ele(param0)
              x <- station.thredds(param=param0,stid=stid[i],loc=loc[i],lon=lon[i],lat=lat[i],
                                   it=as.numeric(c(start[i],end[i])),alt=alt[i],cntr=cntr[i],verbose=verbose)
            } else {
              dd06 <- station.thredds(param='dd06',stid=stid[i],loc=loc[i],lon=lon[i],lat=lat[i],
                                      it=c(start[i],end[i]),alt=alt[i],cntr=cntr[i],verbose=verbose)
              dd12 <- station.thredds(param='dd12',stid=stid[i],loc=loc[i],lon=lon[i],lat=lat[i],
                                      it=c(start[i],end[i]),alt=alt[i],cntr=cntr[i],verbose=verbose)
              dd18 <- station.thredds(param='dd18',stid=stid[i],loc=loc[i],lon=lon[i],lat=lat[i],
                                      it=c(start[i],end[i]),alt=alt[i],cntr=cntr[i],verbose=verbose)
              x <- (dd06 + dd12 + dd18) / 3
              x <- attrcp(dd06,x)
              class(x) <- class(dd06)
              rm(dd06,dd12,dd18)
              if(verbose) print("station.default: WARNING : Averaged wind direction values computed from 06, 12,and 18 UTC")
              param1 <- "DD"
              ele <-  "502"
              attr(x,'variable') <- param1
              attr(x,'element') <- ele
              attr(x,'longname') <- "Average of wind directions at 06 12 and 18 utc" 
            }
          } else if(s=="ECAD") {
            x <- ecad.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
                              cntr=cntr[i],qual=qual[i],param=param0,verbose=verbose,
                              path=path, url=url)
          } else if (s=="NACD") {
            x <- nacd.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
                              cntr=cntr[i],qual=qual[i],param=param0,verbose=verbose)
          } else if (s=="NARP") {
            x <- narp.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
                              cntr=cntr[i],qual=qual[i],param=param0,verbose=verbose)
          } else if (s=="NORDKLIM") {
            x <- nordklim.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],
                                  loc=loc[i],cntr=cntr[i],qual=qual[i],param=param0,
                                  verbose=verbose)
          } else if (s=="GHCNM") {
            x <- ghcnm.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
                               cntr=cntr[i],qual=qual[i],param=param0,verbose=verbose,
                               path = path,url=url)
          } else if (s=="GHCND") {
            if(param0=="t2m") { ## compute the avg 
              ghcnd.tmin <- ghcnd.station.int(param="tmin",stid=stid[i],lon=lon[i],lat=lat[i],
                                          alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],
                                          verbose=verbose,path = path,url=url)
              ghcnd.tmax <- ghcnd.station.int(param="tmax",stid=stid[i],lon=lon[i],lat=lat[i],
                                          alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],
                                          verbose=verbose, path=path,url=url)
              if (is.null(ghcnd.tmax) |  is.null(ghcnd.tmin)) {
                x <- NULL
              } else {
                x <- (ghcnd.tmin + ghcnd.tmax) / 2
                x <- attrcp(ghcnd.tmin,x)
                class(x) <- class(ghcnd.tmin)
                rm(ghcnd.tmin,ghcnd.tmax)
                if(verbose) print("station.default: WARNING : Average temperature values have been computed from TMIN and TMAX values")
                param1 <- "TAVG"
                ele <-  "101"
                attr(x,'variable') <- param
                attr(x,'element') <- ele
                attr(x,'longname') <- "Mean temperature"
              }
            } else {
              x <- ghcnd.station.int(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
                                 cntr=cntr[i],qual=qual[i],param=param0,verbose=verbose,
                                 path = path,url=url)
            }
          }

          ## Let user know if no data was found, otherwise set data to variable X
          if (is.null(x) | (sum(is.na(coredata(x)))==length(coredata(x))) ) {
            if(verbose) print("station.default: Warning : No values found in the time series for-> This station will be ignored")
            if(verbose) print(paste('station.default: stid=',stid[i],'lon=',lon[i],'lat=',lat[i],'alt=',alt[i],
                                    'loc=',loc[i],'cntr=',cntr[i],'param=',param0,'path=',path,
                                    'url=',url)) # REB 2016-07-26
            x <- NULL
          } else {
            if (verbose) {print("station.default: obs"); str(x)}
            if (verbose) print(paste('station.default: Combine the station records for i=',i))
            if (is.null(X)) X <- x else X <- combine.station(X,x)
          }
        }
      }
    }
  }
  if (is.null(X)) return(NULL)
  if (plot & !is.null(X)) map.station(X,col="darkgreen",bg="green", cex=0.7)
  invisible(X)
}

# NOT EXPORTED
t2m.ghcnd.avg <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param="t2m",
                          path = path , ver = "v3.2.0.20130120",adj = TRUE,force=FALSE,flag = FALSE, off = FALSE, verbose = FALSE) {
  ghcnd.tmin <- station(param="tmin",stid = stid, src = "GHCND", verbose = verbose , path = path)
  if (is.null(ghcnd.tmin)) return(NULL)
  ghcnd.tmax <- station(param="tmax",stid = stid, src = "GHCND", verbose = verbose , path = path) 
  if (is.null(ghcnd.tmax)) return(NULL)
  ghcnd <- (ghcnd.tmin + ghcnd.tmax) / 2
  ghcnd <- attrcp(ghcnd.tmin,ghcnd)
  if(verbose) print("station.default: WARNING : Average temperature values have been computed from TMIN and TMAX values")
  param1 <- "TAVG"
  ele <-  "101"
  attr(ghcnd,'variable') <- switch(toupper(param1),'TAVG'=expression(T[2*m]),
                                   'TMAX'=expression(paste("max ",T[2*m])),
                                   'TMIN'=expression(paste("min ",T[2*m])))
  attr(ghcnd,'element') <- ele
  invisible(ghcnd)
}

# NOT EXPORTED - internal function
## KMP 2021-09-22: I updated the url as the old path seems not to work anymore. 
ecad.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,
                         param=NULL,qual=NULL,path="data.ECAD",remove.suspect=FALSE,
                         url="https://knmi-ecad-assets-prd.s3.amazonaws.com/download/ECA_nonblend",
                         #url="http://www.ecad.eu/utils/downloadfile.php?file=download/ECA_nonblend",
                         verbose=FALSE) {
  
  if (verbose) print('station.default: ecad.station...')
  ele <- esd2ele(param=param)
  if (is.null(ele)) {
    param1 <- as.character(ele2param(ele=param,src="ECAD")$param)
  } else {
    param1 <- as.character(ele2param(ele=ele,src="ECAD")$param)
  }
  if (!is.null(param) & (!is.null(dim(param1)[1]))) {
    if (dim(param1)==0) { 
      stop('Please refresh your selection, element not found in meta data')
    }
  }
  scale <-as.numeric(ele2param(ele=ele,src="ECAD")$scale)
  if (verbose) print(paste('station.default: scale=',scale))
  fdata <- paste(url,"_",tolower(param1),".zip",sep="") 
  #text  <- unlist(strsplit(fdata,split="/"))
  #text2 <- text[length(text)]
  #destfile <- file.path(path,text2,fsep= .Platform$file.sep)
  destfile <- file.path(path, basename(fdata), fsep= .Platform$file.sep)
  destfile2 <- file.path(path, paste0('ECA_nonblend_',tolower(param1)),
                         fsep= .Platform$file.sep)
  if (file.exists(destfile) & !file.exists(destfile2)) {
    unzip(destfile,exdir=substr(destfile,1,nchar(destfile)-4))
  } 
  ## If folder does not exist, then download and unzip
  if (!file.exists(destfile2)) {
    print(paste('station.default: Could not find',destfile2))
    print('station.default: Download the ECA&D data... please be patient.')
    if(!file.exists(path)) dir.create(path,showWarnings = FALSE,recursive=TRUE)
    download.file(fdata,destfile,method = "wget", quiet = !verbose, mode = "w", cacheOK = TRUE,
                  extra = getOption("download.file.extra"))
    unzip(destfile,exdir=substr(destfile,1,nchar(destfile)-4))
  }
  if (verbose) print(paste("station.default: station.ecad: the folder has been found",destfile2))
  newpath <- substr(destfile,1,nchar(destfile)-4) 
  stid <- gsub(' ','',stid)
  for (i in 1:length(stid)) {
    while(nchar(stid[i]) < 6) stid[i] <- paste('0',stid[i],sep="")
  }
  fnames <- paste(toupper(param1),'_SOUID',stid,'.txt',sep="")
  fnames <- file.path(newpath,fnames,fsep = .Platform$file.sep)
  
  ipick <- file.exists(fnames)
  if (verbose) print(paste('station.default: Looking for',fnames[1]))
  if (sum(ipick)==0)  return(NULL)
  if (sum(ipick)!=1) {
    warning('station.default: More than one matches - I choose the first!')
    ipick <- (1:length(ipick))[ipick][1]
  }
  fname <- fnames[ipick]
  if (verbose) print(paste('station.default: FILENAME:',fname))
  
  x <- read.table(fname,header=TRUE,skip=18,sep=",")
  
  eval(parse(text=paste0("ecad <- scale * x$",param1)))
  year <- substr(as.character(x$DATE),1,4)
  L <- length(year)
  month <- substr(as.character(x$DATE),5,6)
  day <- substr(as.character(x$DATE),7,8)
  if (verbose) {
    print('station.default: Make use of the quality flags')
    str(x)
  }        
  if (dim(x)[2]==5) dataQ <- x[[5]] else dataQ <- ecad*0               
  ecad[ecad <=-999] <- NA                       
  if (remove.suspect) ecad[dataQ == 1] <- NA    
  ecad[dataQ == 9] <- NA                        
  
  if (verbose) {
    print(paste('station.default:',c(year[1],month[1],day[1])))
    print(paste('station.default:',c(year[L],month[L],day[L])))
    print(paste('station.default:',summary(ecad)))
    print(paste('station.default:',length(ecad)))
    print(paste('station.default:',stid))
  }
  
  if (!check.bad.dates(year,month,day)) {
    ECAD <- zoo(ecad, order.by = as.Date(paste(year, month, day, sep = "-"), 
                                         by='day', length.out = L))
  } else {
    print("station.default: Warning : Bad dates were found for this station -> Ignored")
    return(NULL)
  }
  
  if (sum(ECAD,na.rm=TRUE)==0) {
    print("station.default: Warning : No recorded values are found for this station -> Ignored") ; return(NULL)
  } 
  
  ECAD <- as.station(ECAD, stid=stid, lon=lon, lat=lat, alt=alt,
                     ## ele=esd2ele(param), freq=1,calendar='gregorian',
                     quality=qual, cntr=cntr, loc=loc, src='ECAD', url=fdata,
                     #url=paste0("http://eca.knmi.nl/utils/downloadfile.php?file=download/ECA_nonblend_",as.character(param),".zip"),
                     param=param, aspect="original",
                     unit=switch(param1,'TG'='degree Celsius','TX'='deg C',
                                 'TN'='deg C','CC'='oktas', 'DD'='degrees',
                                 'FG'='m/s','FX'='m/s','HU'='%','PP'='hPa',
                                 'SS'='hours','RR'='mm/day'),
                     longname=as.character(ele2param(ele=ele,src="ECAD")$longname[2]),
                     reference=paste0("Klein Tank, A.M.G. and Coauthors, 2002.",
                                      " Daily dataset of 20th-century surface air temperature and precipitation series for the European Climate Assessment.",
                                      " Int. J. of Climatol., 22, 1441-1453."),
                     info= "Data and metadata available at http://eca.knmi.nl")
  attr(ECAD,'history') <- c(match.call(),date())
  attr(ECAD,'history') <- history.stamp(ECAD)
  if (verbose) print('station.default: --- exit ecad.station ---')
  invisible(ECAD)
}

# NOT EXPORTED - internal function
nacd.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,verbose=FALSE) {
  if(verbose) print("nacd.station")
  if (is.null(stid)) return(NULL)
  
  ele <- esd2ele(param=param)
  
  data("NACD", envir = environment())
  loc <- gsub("-",".",loc)
  loc <- gsub("/",".",loc)
  #id.dot <- grep('.',loc)
  #if (substr(loc,nchar(loc),id.dot[length(id.dot)])) loc <- substr(loc,1,nchar(loc)-1) # AM 29.07.2013 added 
  #elem<-switch(ele.c,'101'='t2m','601'='precip','401'='slp','801'='cloud')
  
  if (tolower(param)=='cc') param <- 'cloud'
  txt <- paste("x <- NACD$",loc,".",param,sep="") # AM 27.08.2013 line updated : "elem" is replaced by "ele"
  eval(parse(text=txt))
  
  if (is.null(x)) {
    print("station.default: No recorded values are found for this station")
    return(NULL)
  }
  
  x.name <- as.character(ele2param(ele=ele,src="NACD")$param[2])
  unit <- as.character(ele2param(ele=ele,src="NACD")$unit[4])
  unit <- switch(ele,
                 '101'='degree Celsius',
                 '111'='degree Celsius',
                 '112'='degree Celsius',
                 '122'='degree Celsius','401'='hPa',
                 '601'='mm/month','701'='days','801'='%')
  ny <- length(attr(x,'year'))
  year <- rep(attr(x,'year'),each=12); L <- length(year)
  month <- rep(1:12,ny)
  day <- rep(1,length(month))
  x[x < -99] <- NA
  
  NACD <- zoo(c(t(x)),order.by = as.Date(paste(year, month, day, sep = "-"), by='month', length.out = L))
  names(NACD) <- loc # AM 30.07.2013 added           
  
  unit <- switch(ele,
                 '101'='degree Celsius','111'='degree Celsius',
                 '112'='degree Celsius','122'='degree Celsius',
                 '401'='hPa', '601'='mm','701'='days','801'='%')
  
  NACD <- as.station(NACD, stid=stid, loc=loc, cntr=cntr,
                     lon=round(lon,digits = 2),
                     lat=round(lat,digits = 2), alt=alt,
                     #calendar='gregorian',
                     quality=qual, src='NACD', url=NA, param=esd2ele(ele),
                     aspect="original", unit=unit, longname=x.name,
                     reference="Frich et al. (1996), DMI scientific report 96-1",
                     info="Data and metadata adopted from clim.pact")
  
  attr(NACD,'history') <- c(match.call(),date())
  attr(NACD,'history') <- history.stamp(NACD)
  invisible(NACD)
}

# NOT EXPORTED
narp.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,verbose=FALSE) {
  if(verbose) print("station.default: narp.station")
  data(NARP,envir=environment())
  ele <- esd2ele(param=param)
  
  x.name <- as.character(ele2param(ele=ele,src="NARP")$param[2])
  unit <-  as.character(ele2param(ele=ele,src="NARP")$unit[4])
  scale <- as.numeric(ele2param(ele=ele,src="NARP")$scale_factor[3])
  iii <- is.element(NARP[,1],stid) & is.element(NARP[,2],as.numeric(ele))
  
  if (sum(iii) ==0) {
    print("station.default: Warning : No recorded values are found for this station -> Ignored")
    return(NULL)
  }
  x <- NARP[iii,4:15]*scale
  x[x <= -99] <- NA
  ny <- length(NARP[iii,3])
  year <- sort(rep(NARP[iii,3],12))
  L <- length(year)
  month <- rep(1:12,ny)
  day <- rep(1,length(year))
  
  NARP <- zoo(c(t(x)), order.by = as.Date(paste(year, month, day, sep = "-"), by='month', length.out = L))
  
  NARP <- as.station(NARP, stid=stid, loc=loc, lon=lon, lat=lat, alt=alt,
                     ## frequency=1, calendar='gregorian',
                     cntr=cntr, quality=NA, src='NARP', url=NA, param=param,
                     unit=unit, longname=x.name,
                     aspect="original",
                     reference="Nordic Arctic Research Programme",
                     info="narp2esd.R")
  attr(NARP,'history') <- c(match.call(),date())
  attr(NARP,'history') <- history.stamp(NARP)
  invisible(NARP)
}

# NOT EXPORTED - internal function
nordklim.station <- function(stid=NULL,loc=NULL,lon=NULL,lat=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,verbose=FALSE,path=NULL) {
  if(verbose) print("station.default: nordklim.station") 
  if (is.null(stid)) return(NULL)
  data(nordklim.data,envir=environment())
  x <- nordklim.data ; rm(nordklim.data)
  station_id <- stid
  x <- subset(x, station_id==stid)
  ele <- esd2ele(param) 
  x <- subset(x,element==ele)
  
  if (dim(x)[1] < 1) {
    print("station.default: Warning : No recorded values are found for this station -> Ignored")
    return(NULL)
  }
  x[x < -999] <- NA
  scale_factor <- as.numeric(as.matrix(ele2param(ele=ele,src="NORDKLIM")$scale_factor[3]))
  xx <- as.matrix(x[,4:15]*scale_factor)
  unit <- as.character(as.matrix(ele2param(ele=ele,src="NORDKLIM")$unit[4]))
  ## Vector of dates
  year <- sort(rep(x$start,12))
  ny <- length(x$start)
  month <- rep(1:12,ny)
  day <- rep(1,length(year))
  ## Longname
  lname <- as.character(ele2param(ele=ele,src="NORDKLIM")$longname[2])
  
  ## format data into a zoo object
  NORDKLIM <- zoo(c(t(xx)),order.by = as.Date(paste(year, month, day, sep = "-"), by='month', length.out = ny))
  names(NORDKLIM) <- loc # AM 30.07.2013 added   
  NORDKLIM <- as.station(NORDKLIM,stid=stid, quality=qual, lon=lon,lat=lat,alt=alt, ##frequency=1,calendar='gregorian',
                         cntr=cntr,loc=loc,src='NORDKLIM', url="http://www.smhi.se/hfa_coord/nordklim/index.php?page=dataset",
                         longname=lname,unit=unit, param=param, aspect="original",
                         reference="Tuomenvirta et al. (2001), NORDKLIM data set 1.0 - Description and illustrations. DNMI report no. 08/01.",
                         info="nordklim2esd.R")
  
  attr(NORDKLIM,'call') <- match.call()
  attr(NORDKLIM,'history') <- c(match.call(),date())
  attr(NORDKLIM,'history') <- history.stamp(NORDKLIM)
  invisible(NORDKLIM)
}

# NOT EXPORTED - internal function
ghcnm.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,ver="v3",
                          path="data.GHCNM",url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn",adj = "qca",force=FALSE,
                          flag=FALSE,verbose=FALSE) {
  
  ele <-esd2ele(param=param) 
  
  ## extract param1 and scale variables
  param1 <-as.character(ele2param(ele=ele,src="GHCNM")$param[5])
  scale <-as.numeric(ele2param(ele=ele,src="GHCNM")$scale_factor[3])
  
  if (verbose) print("station.default: station.GHCNM")
  
  ghcnm <- ghcnm.data(ele=ele,stid=stid,src="ghcnm",ver=ver,adj=adj,path=path,url=url,force=force,flag=flag,verbose=verbose)
  x <- c(t(ghcnm[,5:16]))*scale
  
  year <- sort(rep(ghcnm$year,12)) ; L <- length(year)
  month <-rep(1:12,length(ghcnm$year)) 
  day <- rep("01",L)
  
  GHCNM <- zoo(x,order.by = as.Date(paste(year, month, day, sep = "-"), by='month', length.out = L))
  
  if (sum(GHCNM,na.rm=TRUE)==0) {
    print("station.default: Warning : No recorded values are found for this station -> Ignored")
    return(NULL)
  } 
  ##print("station.default: attributes")
  GHCNM <- as.station(GHCNM,stid=stid, quality=qual, lon=lon,lat=lat,alt=alt,##frequency=1,calendar='gregorian',
                      cntr=cntr, loc=loc, src='GHCNM', url=paste(url,ver,sep="/"),longname=as.character(ele2param(ele=ele,src="GHCNM")$longname[2]),
                      unit=switch(param1,'TAVG'='degree Celsius','TMAX'='degree Celsius','TMIN'='degree Celsius'), param=param, aspect="original",
                      reference=paste0("J. H. Lawrimore, M. J. Menne, B. E. Gleason, C. N. Williams, D. B. Wuertz, R. S. Vose, and J. Rennie ",
                                       "(2011), An overview of the Global Historical Climatology Network monthly mean temperature data set",
                                       ", version 3, J. Geophys. Res., 116, D19121, doi:10.1029/2011JD016187."),
                      info="Data and metadata available at the ftp://ftp.ncdc.noaa.gov/pub/data/ghcn")
  
  
  ## attr(GHCNM,'call') <- match.call()
  attr(GHCNM,'history') <- c(match.call(),date())
  attr(GHCNM,'history') <- history.stamp(GHCNM)
  invisible(GHCNM)
}

# NOT EXPORTED - internal function
ghcnd.station.int <- function(stid=NULL, lon=NULL, lat=NULL, loc=NULL, alt=NULL, cntr=NULL, qual=NULL, param=NULL,
                          path="data.GHCND", url=NULL, adj=TRUE, force=FALSE, flag=FALSE, off=FALSE, verbose=FALSE) {
  
  if (verbose) print("station.default: station.GHCND")
  ele <- esd2ele(param=param) 
  param1 <- as.character(ele2param(ele=ele,src="GHCND")$param[5])
  ## REB 2021-05-11 fix
  #scale <- as.numeric(ele2param(ele=ele,src="GHCND")[3])
  scale <- as.numeric(ele2param(ele=ele,src="GHCND")$scale_factor[3])
  ghcnd <- ghcnd.data(param=param1, stid=stid, src="ghcnd", path=path, url=url,
                      force=force, flag=flag, verbose=verbose, rm.file=FALSE)
  
  if (is.null(ghcnd)) return(NULL)
  
  x <- c(t(ghcnd[,6:36]))*scale
  
  day_id <-substr(names(ghcnd[,6:dim(ghcnd)[2]]),4,nchar(names(ghcnd[,6:dim(ghcnd)[2]]))) 
  year <- rep(as.character(ghcnd$YEAR),each=length(day_id)) ; L <- length(year)
  month <- rep(ghcnd$MONTH,each=length(day_id))
  day <- rep(day_id,dim(ghcnd)[1])
  
  vdate <- as.Date(paste(year, month, day, sep = "-"), by='day', length.out = L)
  id <- !is.na(vdate)
  vdate <- vdate[id]
  x <- x[id]
  
  if (is.null(x)) return(NULL)
  
  GHCND <- zoo(x,order.by = vdate)
  
  if (sum(GHCND,na.rm=TRUE)==0) {
    print("station.default: Warning : No recorded values are found for this station -> Ignored")
    return(NULL)
  } 
  
  GHCND <- as.station(GHCND,stid=stid, quality=qual, lon=lon,lat=lat,alt=alt,##frequency=1,calendar='gregorian',
                      cntr=cntr, loc=loc,src='GHCND', url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn",
                      longname=as.character(ele2param(ele=ele,src="GHCND")$longname[2]),
                      unit=as.character(ele2param(ele=ele,src="GHCND")$unit[4]), param=param, aspect="original",
                      reference=paste0("J. H. Lawrimore, M. J. Menne, B. E. Gleason, C. N. Williams, D. B. Wuertz, R. S. Vose, and J. Rennie ",
                                       "(2011), An overview of the Global Historical Climatology Network monthly mean temperature data set",
                                       ", version 3, J. Geophys. Res., 116, D19121, doi:10.1029/2011JD016187."),
                      info="Data and metadata available at the ftp://ftp.ncdc.noaa.gov/pub/data/ghcn")
  
  attr(GHCND,'call') <- match.call()
  attr(GHCND,'history') <- c(match.call(),date())
  attr(GHCND,'history') <- history.stamp(GHCND)
  ## class(GHCND) <- c("station","day","zoo")
  invisible(GHCND)
}

#' @export station.giss
station.giss <- function(...,url=NULL) {
  t2m <- read.table(url,skip=2,header=TRUE)
  meta <- read.table(url,nrows=1,as.is=TRUE)
  lat <- as.numeric(substr(meta$V3,2,nchar(meta$V3)-2))
  lon <- as.numeric(sub(')','',meta$V4))
  yr <- t2m$YEAR
  t2m[t2m >= 999] <- NA
  t2m <- zoo(c(t(as.matrix(t2m[2:13]))),
             order.by=as.Date(paste(sort(rep(yr,12)),rep(1:12,length(yr)),'01',sep='-')))
  t2m <- as.station(t2m,loc=meta$V1,param='t2m',unit='degC',
                    lon=lon,lat=lat,alt=NA,
                    cntr=meta$V2,longname='Surface temperature',
                    stid=meta$V5,quality=meta$V10,src=paste(meta$V6,meta$V7),url=url,
                    reference='Parker, et. al. (1992), Int. J. Clim.',info=NA, method= NA)
  return(t2m)
}
