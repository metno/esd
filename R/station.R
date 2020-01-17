
library(data.table) # TODO: this should not be here in production

## This function is used to check wether there are errors in the programming
# NOT EXPORTED
test.station <- function(ss=NULL,stid=NULL,alt=NULL,lat=c(50,70),lon=c(0,30),param="precip",
    src=c("GHCND","GHCNM","NORDKLIM","NACD","METNOM","METNOD","METNO.FROST","ECAD"),verbose=FALSE) {
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
#' @aliases station select.station station.default station.ecad station.nacd
#' station.narp station.nordklim station.metnod station.metnom station.ghcnd
#' station.ghcnm station.ghcnm station.sonel station.gloss station.newlyn
#' station.giss
#'
#' @seealso clean.station allgood station.thredds
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
#' @author A. Mezghani
#' @seealso \code{\link{map.station}}.
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

#' @export station.ecad
station.ecad <- function(...) {
  y <- station(src="ecad",...)
  invisible(y)
}

#' @export station.ghcnd
station.ghcnd <- function(...) {
  y <- station(src="ghcnd",...)
  invisible(y)
}

#' @export station.nacd
station.nacd <- function(...) {
  y <- station(src="nacd",...)
  invisible(y)
}

#' @export station.narp
station.narp <- function(...) {
  y <- station(src="narp",...)
  invisible(y)
}

#' @export station.nordklim
station.nordklim <- function(...) {
  ## 
  y <- station(src="nordklim",...)
  invisible(y)
}

#' @export station.metnom
station.metnom <- function(...) {
  y <- station(src="metnom",...)
  invisible(y)
}

#' @export station.metnod
station.metnod <- function(...) {
  y <- station(src="metnod",...)
  invisible(y)
}

#' @export station.metno.frost
station.metno.frost <- function(...) {
  y <- station(src="metno.frost",...)
  invisible(y)
}

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

#' @export station.default
station.default <- function(..., loc=NULL, param='t2m', src=NULL, path=NULL,
                            qual=NULL, url = NULL, stid=NULL, lon=NULL, lat=NULL, alt=NULL,
			    cntr=NULL, it= NULL,nmin=NULL, plot=FALSE, verbose=FALSE,
                            path.ecad=NULL,url.ecad=NULL,
                            path.ghcnm=NULL,url.ghcnm=NULL,
                            path.ghcnd=NULL,url.ghcnd=NULL,
                            path.metnom=NULL,url.metnom=NULL,
                            path.metnod=NULL,url.metnod=NULL,
                            path.metno.frost=NULL,url.metno.frost=NULL,
                            user='external',save2file=TRUE) {
  
  if (verbose) {
    print('station.default')
    print(match.call())
  }
  if (inherits(loc,"stationmeta")) {
    ss <- loc
  } else if (is.character(loc)) {
    loc <- loc
    ss <- NULL
  } else ss <- NULL
  
  X <- NULL
  SRC <- src
  PATH <- path
  URL <- url
  
  if (is.null(ss)) {
    if (verbose) print('select.station')
    ss <- select.station(stid=stid,loc=loc,lon=lon,lat=lat,alt=alt,cntr=cntr,param=param,src=src,it=it,nmin=nmin)
  }  
  if ((param=="t2m") & is.null(ss)) {
    param0 <- param
    ssn <- select.station(param="tmin",stid=stid,loc=loc,lon=lon,lat=lat,alt=alt,cntr=cntr,src=src,it=it,nmin=nmin)
    ssx <- select.station(param="tmax",stid=stid,loc=loc,lon=lon,lat=lat,alt=alt,cntr=cntr,src=src,it=it,nmin=nmin)
    if (!is.null(ssn) & !is.null(ssx)) {
      class(ssn) <- class(ssx) <- "data.frame"
    } else {
      print('Found no stations with given criteria')
      return(NULL)
    }
    ss <- subset(ssx,ssx$station_id==ssn$station_id) # keep only stations recording both min and max
    if (is.null(ss)) {
      return(NULL)
    } else {
      rl <- readline("T2m is not available for your selection but TMIN and TMAX have been found - Would you like to continue using the averaged values? (y or n): ")
      if ((rl=="y") | rl==("ye") | (rl=="yes")) {
        ss$element <- rep(esd2ele(param),length(ss$station_id)) ## update element with param
      } else {
        stop("Process stopped")
      }
    }
  } 

  id <- ss$station_id
  if (verbose) {
    print("Station ID:")
    str(id)
  }

  stid <- ss$station_id
  src <- as.character(ss$source)
  loc <- as.character(ss$location)
  lon <- ss$longitude
  lat <- ss$latitude
  alt <- ss$altitude
  cntr <- as.character(ss$country)
  ele <- ss$element
  qual <- ss$quality
  start <- ss$start
  end <- ss$end
  param <-  apply(as.matrix(ss$element),1,esd2ele)
  rm(ss)
  
  print(paste("Retrieving data from",length(id),"records ..."))
  
  X <- NULL
  
  for (i in 1:length(id)) {
    if (!is.null(param[i]) & !is.null(src[i])) {
      if ((param[i] == "t2m") & (src[i]=="GHCND")) {
        param0 <- param[i] 
      } else param0 <-  NULL
    }
    
    if(verbose) {
      if (!is.null(param0)) {
        print(paste(i,toupper(param0),stid[i],loc[i],cntr[i],src[i]))
      } else {
        print(paste(i,toupper(param[i]),stid[i],loc[i],cntr[i],src[i]))
      }
    }
    
    if (src[i]=="METNOD") { #AM-29.08.2013 added for metno data
      if (is.null(path.metnod)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.metnod ## default path
      if (is.null(url.metnod)) {
          if (user == 'metno') url="http://klapp/metnopub/production/"
          else url= 'ftp://ftp.met.no/projects/chasepl/test'
      } else {
        url <- url.metnod ## default url
      }
      
      if (param[i] == 'dd') {
        param1 <- esd2ele(param[i])
      } else {
        param1 <- NULL
      }
      if (is.null(param1)) {
        x <- metnod.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],
	  start=start[i],end=end[i],qual=qual[i],param=param[i],verbose=verbose, path=path,url=url,user=user,save2file = save2file)
      } else { ## compute the avg 
        dd06 <- metnod.station(param='dd06',stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
	     cntr=cntr[i],start=start[i],end=end[i],qual=qual[i],verbose=verbose, path=path,url=url,user=user,save2file = save2file)
        dd12 <- metnod.station(param='dd12',stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
	     cntr=cntr[i],start=start[i],end=end[i],qual=qual[i],verbose=verbose, path=path,url=url,user=user,save2file = save2file)
        dd18 <- metnod.station(param='dd18',stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
	     cntr=cntr[i],start=start[i],end=end[i],qual=qual[i],verbose=verbose, path=path,url=url,user=user,save2file = save2file)
        x <- (dd06 + dd12 + dd18) / 3
        x <- attrcp(dd06,x)
        class(x) <- class(dd06)
        rm(dd06,dd12,dd18)
        print("WARNING : Averaged wind direction values computed from 06, 12,and 18 UTC")
        param1 <- "DD"
        ele <-  "502"
        attr(x,'variable') <- param
        attr(x,'element') <- ele
        attr(x,'longname') <- "Average of wind directions at 06 12 and 18 utc"
      }
      if (verbose) {print("obs"); str(x)}
      if ( is.null(x)| (sum(is.na(coredata(x)))==length(coredata(x))) ) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
    } else if (src[i]=="METNOM") { #AM-29.08.2013 added for metno data
      if (is.null(path.metnom)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.metnom ## default path
      if (is.null(url.metnom)) url="http://klapp/metnopub/production/" else url <- url.metnom ## default url
      x <- metnom.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],
      	start=start[i],end=end[i],qual=qual[i],param=param[i],verbose=verbose, path=path,url=url,user=user) ## ,path=path, url=url
      if (verbose) {print("obs"); str(x)}
      if ( is.null(x) | (sum(is.na(coredata(x)))==length(coredata(x))) ) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
      ##if (!is.null(x)) X <- combine.stations(X,x)
    } else if (src[i]=="METNO.FROST.DIURNAL") {
      x <- metno.frost.station(timeresolutions='P1D',...)
    } else if (src[i]=="METNO.FROST.MONTH") {
      x <- metno.frost.station(timeresolutions='P1M',...)
    } else if (src[i]=="ECAD") {
      if (is.null(path.ecad)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.ecad ## default path
      if (is.null(url.ecad)) url="http://www.ecad.eu/utils/downloadfile.php?file=download/ECA_nonblend" else url <- url.ecad ## default url
      x <- ecad.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],
                        qual=qual[i],param=param[i],verbose=verbose,path=path, url=url)
      if (verbose) {print("obs"); str(x)}
      
      if ( is.null(x) | (sum(is.na(coredata(x)))==length(coredata(x))) ) {
        print("Warning : No values found in the time series for-> This station will be ignored")
        print(paste('stid=',stid[i],'lon=',lon[i],'lat=',lat[i],'alt=',alt[i],'loc=',loc[i],'cntr=',
                    cntr[i],'param=',param[i],'path=',path,'url=',url)) # REB 2016-07-26
        x <- NULL
      }
    } else if (src[i]=="NACD") {
      x <- nacd.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose)
      if (verbose) {print("obs"); str(x)}
      if ( is.null(x) | (sum(is.na(coredata(x)))==length(coredata(x))) ) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
    } else if (src[i]=="NARP") {
      x <- narp.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose)
      if (verbose) {print("obs"); str(x)}
      if ( is.null(x) | (sum(is.na(coredata(x)))==length(coredata(x))) ) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
    } else if (src[i]=="NORDKLIM") {
      x <- nordklim.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose)
      if (verbose) {print("obs"); str(x)}
      if ( is.null(x) | (sum(is.na(coredata(x)))==length(coredata(x))) ) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
    } else if (src[i]=="GHCNM") {
      if (is.null(path.ghcnm)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.ghcnm ## default path
      if (is.null(url.ghcnm)) url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn" else url <- url.ghcnm ## default url
      x <- ghcnm.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
      	cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose,path = path,url=url)
      if (verbose) {print("obs"); str(x)}
      if ( is.null(x) | (sum(is.na(coredata(x)))==length(coredata(x))) ) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
    } else if (src[i]=="GHCND") {
      if (is.null(path.ghcnd)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.ghcnd## default path
      if (is.null(url.ghcnd)) url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all" else url <- url.ghcnd ## default url
      if (!is.null(param0)) { ## compute the avg 
        ghcnd.tmin <- ghcnd.station(param="tmin",stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
		   cntr=cntr[i],qual=qual[i],verbose=verbose,path = path,url=url)
        ghcnd.tmax <- ghcnd.station(param="tmax",stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],
		   cntr=cntr[i],qual=qual[i],verbose=verbose,path = path,url=url)
        if (is.null(ghcnd.tmax) |  is.null(ghcnd.tmin)) {
          x <- NULL
        } else {
          x <- (ghcnd.tmin + ghcnd.tmax) / 2
          x <- attrcp(ghcnd.tmin,x)
          class(x) <- class(ghcnd.tmin)
          rm(ghcnd.tmin,ghcnd.tmax)
          print("WARNING : Average temperature values have been computed from TMIN and TMAX values")
          param1 <- "TAVG"
          ele <-  "101"
          attr(x,'variable') <- param
          attr(x,'element') <- ele
          attr(x,'longname') <- "Mean temperature"
        }
      }
      else {
        x <- ghcnd.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],
	  qual=qual[i],param=param[i],verbose=verbose,path = path,url=url)
        if (verbose) {print("obs"); str(x)}
        if ( is.null(x) | (sum(is.na(coredata(x)))==length(coredata(x))) ) {
          print("Warning : No values found in the time series -> This station will be ignored")
          x <- NULL
        }
      }
    }
    
    if (verbose) print(paste('Combine the station records for i=',i))
    if (!is.null(x)) {
      if (is.null(X)) X <- x else X <- combine.stations(X,x)
    }
  }
  if (is.null(X)) return(NULL)
  
  if (plot & !is.null(X)) map.station(X,col="darkgreen",bg="green", cex=0.7)
  
  if (is.null(stid)) {
    x <- NULL
    print(match.call())
  }
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
  print("WARNING : Average temperature values have been computed from TMIN and TMAX values")
  param1 <- "TAVG"
  ele <-  "101"
  attr(ghcnd,'variable') <- switch(param1,'TAVG'=expression(T[2*m]),'TMAX'=expression(paste("max ",T[2*m])),'TMIN'=expression(paste("min ",T[2*m])))
  attr(ghcnd,'element') <- ele
  invisible(ghcnd)
}

# NOT EXPORTED - internal function
ecad.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,
                         param=NULL,qual=NULL,path="data.ECAD",remove.suspect=FALSE,
                         url="http://www.ecad.eu/utils/downloadfile.php?file=download/ECA_nonblend",verbose=FALSE) {  ## it=it,nmin=nmin

  if (verbose) print('ecad.station')
  ele <- esd2ele(param=param)
  if (is.null(ele)) {
    param1 <-as.character(ele2param(ele=param,src="ECAD")[5])
  } else {
    param1 <- as.character(ele2param(ele=ele,src="ECAD")[5])
  }
  if (!is.null(param) & (!is.null(dim(param1)[1]))) {
    if (dim(param1)==0) { 
      stop('Please refresh your selection, element not found in meta data')
    }
  }
  scale <-as.numeric(ele2param(ele=ele,src="ECAD")[3])
  if (verbose) print(paste('scale=',scale))                      
 
  fdata <- paste(url,"_",tolower(param1),".zip",sep="") 
  text  <- unlist(strsplit(fdata,split="/"))
  text2 <- text[length(text)]
  destfile <- file.path(path,text2,fsep= .Platform$file.sep)
  text3 <- paste('ECA','nonblend',tolower(param1),sep='_')
  destfile2 <- file.path(path,text3,fsep= .Platform$file.sep)

  if (file.exists(destfile) & !file.exists(destfile2)) {
    unzip(destfile,exdir=substr(destfile,1,nchar(destfile)-4))
  } 
  ## If folder does not exist, then download and unzip
  if (!file.exists(destfile2)) {
    if(!file.exists(path)) dir.create(path,showWarnings = FALSE,recursive=TRUE)
    download.file(fdata,destfile,method = "wget", quiet = !verbose, mode = "w", cacheOK = TRUE,
                  extra = getOption("download.file.extra"))
    unzip(destfile,exdir=substr(destfile,1,nchar(destfile)-4))
  }
  if (verbose) print("station.ecad")
  newpath <- substr(destfile,1,nchar(destfile)-4) 
  stid <- gsub(' ','',stid)
  for (i in 1:length(stid)) {
    while(nchar(stid[i]) < 6) stid[i] <- paste('0',stid[i],sep="")
  }
  fnames <- paste(toupper(param1),'_SOUID',stid,'.txt',sep="")
  fnames <- file.path(newpath,fnames,fsep = .Platform$file.sep)

  ipick <- file.exists(fnames)
  if (sum(ipick)==0)  return(NULL)
  if (sum(ipick)!=1) {
    warning('More than one matches - I choose the first!')
    ipick <- (1:length(ipick))[ipick][1]
  }
  fname <- fnames[ipick]
  if (verbose) print(fname)
  
  x <- read.table(fname,header=TRUE,skip=18,sep=",")

  eval(parse(text=paste("ecad <- scale * x$",param1,sep="")))
  year <- substr(as.character(x$DATE),1,4);L <- length(year)
  month <- substr(as.character(x$DATE),5,6)
  day <- substr(as.character(x$DATE),7,8)
  if (verbose) {
    print('Make use of the quality flags')
    str(x)
  }        
  if (dim(x)[2]==5) dataQ <- x[[5]] else dataQ <- ecad*0               
  ecad[ecad <=-999] <- NA                       
  if (remove.suspect) ecad[dataQ == 1] <- NA    
  ecad[dataQ == 9] <- NA                        
  
  if (verbose) {
    print(c(year[1],month[1],day[1]))
    print(c(year[L],month[L],day[L]))
    print(summary(ecad))
    print(length(ecad))
    print(stid)
  }
 
  if (!check.bad.dates(year,month,day))  
             ECAD <- zoo(ecad, order.by = as.Date(paste(year, month, day, sep = "-"), 
                                                  by='day', length.out = L)) else {
             print("Warning : Bad dates were found for this station -> Ignored")
             return(NULL)                                          
  }
  
  if (sum(ECAD,na.rm=TRUE)==0) {
    print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)
  } 
  ECAD <- as.station(ECAD, stid=stid, lon=lon, lat=lat, alt=alt,
                     ## ele=esd2ele(param), freq=1,calendar='gregorian',
                     quality=qual, cntr=cntr, loc=loc, src='ECAD',
                     url=paste("http://eca.knmi.nl/utils/downloadfile.php?file=download/ECA_nonblend_",as.character(param),".zip",sep=''),
			param=param, aspect="original",
                     unit=switch(param1,'TG'='degree Celsius','TX'='deg C','TN'='deg C', 'CC'='oktas',
		     				     'DD'='degrees','FG'='m/s', 'FX'='m/s','HU'='%','PP'='hPa', 'SS'='hours','RR'='mm/day'),
                     longname=as.character(ele2param(ele=ele,src="ECAD")[2]),
                     reference="Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of 20th-century surface air temperature and precipitation series for the European Climate Assessment. Int. J. of Climatol., 22, 1441-1453.",
                     info= "Data and metadata available at http://eca.knmi.nl")
  attr(ECAD,'history') <- c(match.call(),date())
  attr(ECAD,'history') <- history.stamp(ECAD)
  if (verbose) print('--- exit ecad.station ---')
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
    print("No recorded values are found for this station")
    return(NULL)
  }
  
  x.name <- as.character(ele2param(ele=ele,src="NACD")[2])
  unit <- as.character(ele2param(ele=ele,src="NACD")[4])
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
  if(verbose) print("narp.station")
  data(NARP,envir=environment())
  ele <- esd2ele(param=param)
  
  x.name <- as.character(ele2param(ele=ele,src="NARP")[2])
  unit <-  as.character(ele2param(ele=ele,src="NARP")[4])
  scale <- as.numeric(ele2param(ele=ele,src="NARP")[3])
  iii <- is.element(NARP[,1],stid) & is.element(NARP[,2],as.numeric(ele))
 
  if (sum(iii) ==0) {
    print("Warning : No recorded values are found for this station -> Ignored")
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
  if(verbose) print("nordklim.station") 
  if (is.null(stid)) return(NULL)
  data(nordklim.data,envir=environment())
  x <- nordklim.data ; rm(nordklim.data)
  station_id <- stid
  x <- subset(x, station_id==stid)
  ele <- esd2ele(param) 
  x <- subset(x,element==ele)

  if (dim(x)[1] < 1) {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)}
  x[x < -999] <- NA
  scale_factor <- as.numeric(as.matrix(ele2param(ele=ele,src="NORDKLIM")[3]))
  xx <- as.matrix(x[,4:15]*scale_factor)
  unit <- as.character(as.matrix(ele2param(ele=ele,src="NORDKLIM")[4]))
  ## Vector of dates
  year <- sort(rep(x$start,12))
  ny <- length(x$start)
  month <- rep(1:12,ny)
  day <- rep(1,length(year))
  ## Longname
  lname <- as.character(ele2param(ele=ele,src="NORDKLIM")[2])
  
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
ghcnm.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,ver="v3",path="data.GHCNM",url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn",adj = "qca",force=FALSE,flag = FALSE, verbose = FALSE) {

  ele <-esd2ele(param=param) 
  
  ## extract param1 and scale variables
  param1 <-as.character(ele2param(ele=ele,src="GHCNM")[5])
  scale <-as.numeric(ele2param(ele=ele,src="GHCNM")[3])
  
  if (verbose) print("station.GHCNM")

  ghcnm <- ghcnm.data(ele=ele,stid = stid, src = "ghcnm", ver = ver , adj = adj, path = path, url=url,force = force, flag = flag, verbose = verbose)
  x <- c(t(ghcnm[,5:16]))*scale

  year <- sort(rep(ghcnm$year,12)) ; L <- length(year)
  month <-rep(1:12,length(ghcnm$year)) 
  day <- rep("01",L)
  
  GHCNM <- zoo(x,order.by = as.Date(paste(year, month, day, sep = "-"), by='month', length.out = L))
  
  if (sum(GHCNM,na.rm=TRUE)==0) {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)} 
  ##print("attributes")
  GHCNM <- as.station(GHCNM,stid=stid, quality=qual, lon=lon,lat=lat,alt=alt,##frequency=1,calendar='gregorian',
                      cntr=cntr, loc=loc, src='GHCNM', url=paste(url,ver,sep="/"),longname=as.character(ele2param(ele=ele,src="GHCNM")[2]),
                      unit=switch(param1,'TAVG'='degree Celsius','TMAX'='degree Celsius','TMIN'='degree Celsius'), param=param, aspect="original",
                      reference="J. H. Lawrimore, M. J. Menne, B. E. Gleason, C. N. Williams, D. B. Wuertz, R. S. Vose, and J. Rennie (2011), An overview of the Global Historical Climatology Network monthly mean temperature data set, version 3, J. Geophys. Res., 116, D19121, doi:10.1029/2011JD016187.",
                      info="Data and metadata available at the ftp://ftp.ncdc.noaa.gov/pub/data/ghcn")
  
  
  ## attr(GHCNM,'call') <- match.call()
  attr(GHCNM,'history') <- c(match.call(),date())
  attr(GHCNM,'history') <- history.stamp(GHCNM)
  invisible(GHCNM)
}

# NOT EXPORTED - internal function
ghcnd.station <- function(stid=NULL, lon=NULL, lat=NULL, loc=NULL, alt=NULL, cntr=NULL, qual=NULL, param=NULL,
	      path="data.GHCND", url=NULL, adj=TRUE, force=FALSE, flag=FALSE, off=FALSE, verbose=FALSE) {
  
  if (verbose) print("station.GHCND")
  
  ele <- esd2ele(param=param) 
  param1 <- as.character(ele2param(ele=ele,src="GHCND")[5])
  scale <- as.numeric(ele2param(ele=ele,src="GHCND")[3])
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
  
  if (sum(GHCND,na.rm=TRUE)==0) {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)} 
  
  GHCND <- as.station(GHCND,stid=stid, quality=qual, lon=lon,lat=lat,alt=alt,##frequency=1,calendar='gregorian',
                      cntr=cntr, loc=loc,src='GHCND', url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn",longname=as.character(ele2param(ele=ele,src="GHCND")[2]),
                      unit=as.character(ele2param(ele=ele,src="GHCND")[4]), param=param, aspect="original",
                      reference="J. H. Lawrimore, M. J. Menne, B. E. Gleason, C. N. Williams, D. B. Wuertz, R. S. Vose, and J. Rennie (2011), An overview of the Global Historical Climatology Network monthly mean temperature data set, version 3, J. Geophys. Res., 116, D19121, doi:10.1029/2011JD016187.",
                      info="Data and metadata available at the ftp://ftp.ncdc.noaa.gov/pub/data/ghcn")
  
  attr(GHCND,'call') <- match.call()
  attr(GHCND,'history') <- c(match.call(),date())
  attr(GHCND,'history') <- history.stamp(GHCND)
  ## class(GHCND) <- c("station","day","zoo")
  invisible(GHCND)
}

# NOT EXPORTED - internal function
metnom.station <-  function(re=15,stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,
                            start=NULL,end=NULL,param=NULL,verbose=FALSE, h = NULL, nmt = 0,
                            path = NULL, dup = "A", user='metno',url = "http://klapp/metnopub/production/",save2file = TRUE) {
  
  if (user=='metno') {
    y <- metno.station.internal(re=re,stid=stid,lon=lon,lat=lat,loc=loc,alt=alt,cntr=cntr,qual=qual,
                       start=start,end=end,param=param,verbose=verbose,h = h, nmt = nmt,
                       path = path, dup = dup, url = url,save2file = save2file)
  } else {
    y <- metno.station(re=re,stid=stid,lon=lon,lat=lat,loc=loc,alt=alt,cntr=cntr,qual=qual,
                       start=start,end=end,param=param,verbose=verbose,h = h, nmt = nmt,
                       path = path, dup = dup, url = url)
  }
  if (!is.null(y))
    attr(y,"source") <- "METNOM"
  invisible(y)
}

# NOT EXPORTED
metnod.station <-  function(re=14, url = 'ftp://ftp.met.no/projects/chasepl/test', user='else',save2file = TRUE,...) {
  ## 
  ## url <- "ftp://ftp.met.no/projects/chasepl/test"
  if (user=='metno') {
    y <- metno.station.internal(re=re,url=url,save2file = save2file,...)
  } else { 
    y <- metno.station(re=re,url=url,save2file = save2file,...)
  }
  if (!is.null(y)) 
    attr(y,"source") <- "METNOD"
  invisible(y)
}

# NOT EXPORTED - internal function
metno.station.internal <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,
                          qual=NULL,start=NULL,end=NULL,param=NULL,verbose=FALSE,
                          re = 14,h = NULL, nmt = 0,  path = NULL, dup = "A",qa = 4,
                          url = "http://klapp/metnopub/production/",save2file=FALSE) {
  if (verbose) print("http://eklima.met.no")
  if (!is.na(end)) end1 <-format(as.Date(paste("31.12.",as.character(end),sep=""),format='%d.%m.%Y'),'%d.%m.%Y')
  if (!is.na(start)) start1 <- format(as.Date(paste("01.01.",as.character(start),sep=""),format='%d.%m.%Y'),'%d.%m.%Y')
  if (!is.null(url)) {
    Filnavn <- paste(url, "metno?re=", re, "&ct=text/plain&del=space&ddel=dot&nod=NA&split=1", sep = "")
    if (!is.null(h)) {
      Filnavn <- paste(Filnavn, "&h=", h, sep = "")
    }
    param1 <- ele2param(ele=esd2ele(param),src="metno")$param ## switch(param,"t2m"="TAM","precip"="RR")
    for (i in 1:length(param1)) Filnavn <- paste(Filnavn,"&p=", param1[i], sep = "")
    if (!is.null(h)) {
      Filnavn <- paste(Filnavn, "&nmt=", nmt, sep = "")
    }
    Filnavn <- paste(Filnavn, "&fd=", start1, "&td=", end1, sep = "")
    Filnavn <- paste(Filnavn, "&s=", stid, sep = "")
    Filnavn <- paste(Filnavn, "&qa=", qa, sep = "")
    if (!is.null(h)) 
      Filnavn <- paste(Filnavn, "&dup=", dup, sep = "")
  } else stop("The url must be specified")
  
  if (verbose) print(Filnavn)
  firstline <- readLines(Filnavn, n = 1, encoding = "latin1")
  if (substr(firstline, 1, 3) == "***") {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)}

  Datasett <- as.list(read.table(Filnavn,dec = ".", header = TRUE, as.is = TRUE, fileEncoding = "latin1"))
  
  if (param1=='RR') {
    Datasett$RR[Datasett$RR == "."] <- "0"
    Datasett$RR[Datasett$RR == "x"] <- NA
  } else if (param1 == 'SA') {
    Datasett$SA[Datasett$SA == "."] <- "0"
  }
  ext <- switch(as.character(re), '14' = 'dly', '17' = 'obs', '15' = 'mon')
  
  if (save2file) {
    if (is.null(path)) 
      path <- 'data.METNO'
    dir.create(path,showWarnings = FALSE,recursive = TRUE)
    #if (nchar(stid)<=5) 
    stid <- sprintf("%05d", as.numeric(stid))
    write.table(Datasett,file=file.path(path,paste(param1,'_',stid,'.',ext,sep='')),row.names = FALSE,col.names = names(Datasett))
  }
  eval(parse(text = paste("y <- as.numeric(Datasett$", param1,")", sep = "")))
  if (sum(y,na.rm=TRUE)==0) {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)}
  
  type <- switch(re, '14' = "daily values", '17' = "observations", '15' = "Monthly means")
  if (is.na(end)) end <- format(Sys.time(),'%Y')
  year <- Datasett$Year ## sort(rep(c(start:end),12))
  ny <- length(year)
  month <- Datasett$Month ## rep(1:12,length(c(start:end)))
  if (re==14) day <- Datasett$Day else day <- "01" ## rep(1,length(year))
  
  METNO <- zoo(y,order.by = as.Date(paste(year, month, day, sep = "-")))                                  
  
  if (sum(METNO,na.rm=TRUE)==0) {
    print("Warning : No recorded values are found for this station -> Ignored")
    return(NULL)
  } 

  METNO <- as.station(METNO,stid=stid, quality=qual, lon=lon,lat=lat,alt=alt,
                      ##frequency=1,calendar='gregorian',
                      cntr=cntr,loc=loc,src='METNO', url=Filnavn,
                      longname=as.character(ele2param(ele=esd2ele(param),src="METNO")[2]),
                      unit=as.character(ele2param(ele=esd2ele(param),src="METNO")[4]),
                      param=param, aspect="original",
                      reference="Klimadata Vare Huset archive (http://eklima.met.no)",
                      info="Klima Data Vare Huset archive (http://eklima.met.no)")
  
  attr(METNO,'history') <- c(match.call(),date())
  attr(METNO,'history') <- history.stamp(METNO)
  if (re==14) class(METNO) <- c("station","day","zoo") else if (re==15) class(METNO) <- c("station","month","zoo")
  invisible(METNO)
}

#' MetNo meta data function
#'
#' @export
stnr <- function (navn = NULL, lon = NULL, lat = NULL, max.dist = 10, 
                  alt = NULL, Fylke = NULL, Kommune = NULL, fy = NULL, ty = NULL, 
                  ny = NULL, param = "TAM", plot = FALSE, print = FALSE) {
  met.no.meta <- MET.no.meta(param = param, print = print)
  iue <- nchar(met.no.meta$TODATE) == 2
  met.no.meta$TODATE[iue] <- Sys.time()
  i9c <- (nchar(met.no.meta$TODATE) == 9)
  met.no.meta$TODATE[i9c] <- paste("0", met.no.meta$TODATE[i9c], 
                                   sep = "")
  nyrs <- as.numeric(substr(met.no.meta$TODATE, 7, 10)) -
    as.numeric(substr(met.no.meta$FROMDATE, 7, 10)) + 1
  if (!is.null(ny)) {
    keep <- (nyrs >= ny) & (is.finite(nyrs))
    print(summary(nyrs))
    print(paste("Only stations with", ny, "years of data:", 
                sum(keep), "in total"))
    met.no.meta <- met.no.meta[keep, ]
  }
  met.no.meta$Lon[!is.finite(met.no.meta$Lon)] <- -90
  met.no.meta$Lat[!is.finite(met.no.meta$Lat)] <- -90
  ii <- 1:length(met.no.meta$STNR)
  if (!is.null(navn)) {
    ii <- grep(toupper(navn), met.no.meta$Navn)
    print(navn)
    print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii]))
  }
  II <- ii # Is II supposed to be the same as ii?
  if (plot) {
    data("geoborders",envir=environment())
    plot(c(0, 32), c(57, 73), type = "n", xlab = "lon", ylab = "lat")
    lines(geoborders$x,geoborders$y,col="darkblue")
    lines(attr(geoborders,'borders')$x,attr(geoborders,'borders')$y,col="pink")
    lines(geoborders$x+360,geoborders$y,col="darkblue")
    points(met.no.meta$Lon, met.no.meta$Lat, col = "grey", cex = 0.8)
  }
  if (xor(is.null(lon), is.null(lat))) 
    stop("both or none of lon/lat must be specified")
  if (!is.null(lon)) {
    if (length(lon) == 1) {
      if (plot) points(lon, lat, pch = "+", col = "blue", cex = 0.7)
      d <- round(distAB(lon,lat,met.no.meta$Lon, met.no.meta$Lat)/1000, 3)
      print(length(d))
      ii <- II[(d <= max.dist)]
      print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii], 
                  met.no.meta$Lon[ii], met.no.meta$Lat[ii], d[ii]))
    } else if (length(lon) == 2) {
      if (plot) polygon(c(lon[1], lon[2], lon[2], lon[1], lon[1]), 
                        c(lat[1], lat[1], lat[2], lat[2], lat[1]), 
                        border = "blue", lwd = 2)
      if (length(lat)==1) stop("both or none of lon/lat must have two entries")
      ii <- II[(met.no.meta$Lon >= min(lon)) & 
               (met.no.meta$Lon <= max(lon)) & 
               (met.no.meta$Lat >= min(lat)) & 
               (met.no.meta$Lat <= max(lat))]
      print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii], 
                  met.no.meta$Lon[ii], met.no.meta$Lat[ii]))
      met.no.meta <- met.no.meta[ii, ]
    }
  }
  if (!is.null(alt)) {
    if (length(alt) == 1) {
      if (alt > 0) {
        ii <- (met.no.meta$Hoh >= alt) & is.finite(met.no.meta$Hoh)
      } else {
        ii <- (met.no.meta$Hoh <= abs(alt)) &
               is.finite(met.no.meta$Hoh)
      }
    } else {
      ii <- (met.no.meta$Hoh >= min(alt)) & 
            (met.no.meta$Hoh <= max(alt))
    }
    ii[is.na(met.no.meta$Stnr[ii])] <- FALSE
    print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii],
                met.no.meta$Hoh[ii]))
    met.no.meta <- met.no.meta[ii, ]
  }
  if (!is.null(Fylke)) {
    ii <- is.element(toupper(met.no.meta$Fylke), toupper(Fylke))
    print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii],
                met.no.meta$Fylke[ii], 
                met.no.meta$Kommune[ii]))
    met.no.meta <- met.no.meta[ii, ]
  }
  if (!is.null(Kommune)) {
    ii <- is.element(toupper(met.no.meta$Kommune), toupper(Kommune))
    print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii],
                met.no.meta$Fylke[ii], 
                met.no.meta$Kommune[ii]))
    met.no.meta <- met.no.meta[ii, ]
  }
  if (plot) {
    points(met.no.meta$Lon, met.no.meta$Lat, pch = 19, col = "red", 
           cex = 0.6)
    text(met.no.meta$Lon, met.no.meta$Lat, met.no.meta$Stnr, cex = 0.5)
  }
  invisible(met.no.meta)
}

# internal function - no need to export
MET.no.meta <- function (param = "TAM", print = FALSE) {
  url <- paste("http://klapp/metnopub/production/metno?re=27&ct=text/plain&del=semicolon&tab=T_ELEM_MONTH&p=", 
               param, "&geo=lat&geo=utm&geo=amsl&geo=name&geo=cnr&geo=muni&nod=NA", 
               sep = "")
  dnmi.meta <- read.table(url, header = TRUE, sep = ";", as.is = TRUE, 
                          fileEncoding = "latin1")
  if (print) {
    print(url)
    print(summary(dnmi.meta))
  }
  dnmi.meta$Stnr <- as.numeric(dnmi.meta$STNR)
  dnmi.meta$Lon <- as.numeric(dnmi.meta$LON)
  dnmi.meta$Lat <- as.numeric(dnmi.meta$LAT)
  dnmi.meta$Hoh <- as.numeric(dnmi.meta$AMSL)
  dnmi.meta$Navn <- dnmi.meta$ST_NAME
  dnmi.meta$Fylke <- dnmi.meta$COUNTY
  dnmi.meta$Kommune <- dnmi.meta$MUNICIPALITY
  invisible(dnmi.meta)
}

# NOT EXPORTED - internal function
metno.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,
                          qual=NULL,start=NA,end=NA,param=NULL,verbose=FALSE,
                          re = 14,h = NULL, nmt = 0,  path = NULL, dup = "A",
                          url = "ftp://ftp.met.no/projects/chasepl/test",save2file=FALSE) {
  if (verbose) print("metno.station - access data from ftp.met.no")
  ## 
  ## if (!is.na(end)) end1 <- format(Sys.time(),'%d.%m.%Y')
  #if (!is.na(end)) end1 <-format(as.Date(paste("31.12.",as.character(end),sep=""),format='%d.%m.%Y'),'%d.%m.%Y')
  #if (!is.na(start)) start1 <- format(as.Date(paste("01.01.",as.character(start),sep=""),format='%d.%m.%Y'),'%d.%m.%Y')
  param1 <- ele2param(ele=esd2ele(param),src="metno")$param ## switch(param,"t2m"="TAM","precip"="RR")
  ext <- switch(as.character(re), '14' = 'dly', '17' = 'obs', '15' = 'mon')
  if (ext=='dly') {
    path <- 'data.METNOD'
  } else if (ext=='mon') {
    path <- 'data.METNOM'
  }
  if (!is.null(url)) {
    Filnavn <- file.path(url, path, paste(param1,'_',sprintf('%05d',as.numeric(stid)),'.',ext, sep = ""))
  } else stop("The url must be specified")
  if (verbose) print(Filnavn)
  Datasett <- as.list(read.table(Filnavn,dec = ".", header = TRUE, as.is = TRUE, fileEncoding = "latin1"))
  if (param1=='RR') 
    Datasett$RR[Datasett$RR == "."] <- "0"
  
  if (save2file) {
    if (is.null(path)) 
      path <- 'data.METNO'
    dir.create(path,showWarnings = FALSE, recursive = TRUE)
    #if (nchar(stid)<=5) 
    stid <- sprintf("%05d", as.numeric(stid))
    write.table(Datasett,file=file.path(path,paste(param1,'_',stid,'.',ext,sep='')),row.names = FALSE,col.names = names(Datasett))
  }
  
  eval(parse(text = paste("y <- as.numeric(Datasett$", param1, ")",sep = "")))
  if (sum(y,na.rm=TRUE)==0) {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)}
  
  type <- switch(re, '14' = "daily values", '17' = "observations", '15' = "Monthly means")
  ## 
  if (is.na(end)) end <- format(Sys.time(),'%Y')
  year <- Datasett$Year ## sort(rep(c(start:end),12))
  ny <- length(year)
  month <- Datasett$Month ## rep(1:12,length(c(start:end)))
  if (re==14) day <- Datasett$Day else day <- "01" ## rep(1,length(year))
  
  METNO <- zoo(y,order.by = as.Date(paste(year, month, day, sep = "-")))                                  
  
  if (sum(METNO,na.rm=TRUE)==0) {
    print("Warning : No recorded values are found for this station -> Ignored")
    return(NULL)
  } 
  ##print("attributes")
  
  ## Add meta data as attributes:
  METNO <- as.station(METNO,stid=stid, quality=qual, lon=lon,lat=lat,alt=alt,
                      ##frequency=1,calendar='gregorian',
                      cntr=cntr,loc=loc,src='METNO', url=Filnavn,
                      longname=as.character(ele2param(ele=esd2ele(param),src="METNO")[2]),
                      unit=as.character(ele2param(ele=esd2ele(param),src="METNO")[4]),
                      param=param, aspect="original",
                      reference="Klimadata Vare Huset archive (http://eklima.met.no)",
                      info="Klima Data Vare Huset archive (http://eklima.met.no)")
  
  attr(METNO,'history') <- c(match.call(),date())
  attr(METNO,'history') <- history.stamp(METNO)
  if (re==14) class(METNO) <- c("station","day","zoo") else if (re==15) class(METNO) <- c("station","month","zoo")
  #close(Filenavn)
  invisible(METNO)
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

metno.frost.station <- function(stid=NULL, param=NULL, start=NULL, end=NULL,
                          lon=NULL, lat=NULL, loc=NULL, alt=NULL, cntr=NULL,
                          timeresolutions='P1M', levels="default", timeoffsets="default", 
                          performancecategories="A,B,C", exposurecategories="1,2", qualities='0,1,2,3,4,5',
                          path=NULL, save2file=FALSE, verbose=FALSE) {
  ## This downloads METNO observational data,
  ## using the data API frost.met.no.
  ##
  ## Author: K. Tunheim
  
  ## tut <- metno.frost.station(stid='18700', param='t2m', start='2017-01-01', end='2018-01-01', cntr='Norway', timeresolutions='P1M', verbose=TRUE)

  if (verbose) print("http://frost.met.no")
  
  ## TODO 1: if stnr not defined but (lat,lon) is defined, then try to find nearest by looking at the metadata table

  if (is.null(stid) || is.null(param) || is.null(start) || is.null(end))
    stop("stid, param, start and end must be defined")
  
  ## TODO 2: lat, lon, alt and cntr are not used in the call, that will also require a lookup in the metadata
  if (timeresolutions=="P1M") {
    ## FETCH FROM:
    ## metno.frost.meta.month.rda
    
  } else if (timeresolutions=="P1D") {
    ## FETCH FROM:
    ## metno.frost.meta.diurnal.rda
  }

  ## TODO? if stid == NULL, actually use station meta to fetch all data for given period, assuming period is not too long?

  ## TODO: get a client_id
  client_id <- '0763dab1-d398-4a56-ba5d-601d7d352999'
  
  param1info <- ele2param(esd2ele(param), src="metno.frost")
  param1 <- gsub('*', timeresolutions, param1info$param, fixed=TRUE)

  url <- paste0(
      "https://", client_id, "@frost.met.no/observations/v0.jsonld",
      "?sources=", paste0('SN',stid),
      "&referencetime=", paste0(start,"/",end),
      "&timeresolutions=", timeresolutions,
      "&elements=", param1,
      "&levels=", levels,
      "&timeoffsets=", timeoffsets,
      "&performancecategories=", performancecategories,
      "&exposurecategories=", exposurecategories,
      "&qualities=", qualities
  )
  
  if (verbose) print(url)

  xs <- try(jsonlite::fromJSON(URLencode(url),flatten=TRUE))
  if (class(xs) == 'try-error') stop("Data retrieval from frost.met.no was not successful")
  data <- xs$data

  df <- data.table()
  for (i in 1:length(data$observations)) {
      row <- data$observations[[i]]
      row$sourceId <- data$sourceId[[i]]
      row$referenceTime <- data$referenceTime[[i]]
      df <- data.table::rbindlist(list(df, row), fill=TRUE)
  }

  df <- as.data.frame(df)[c("referenceTime","value")]
  df$referenceTime <- as.Date(df$referenceTime)
  
  # TODO: how to (get and) use zoo?
  var <- zoo(df)

  ## TODO: will this work?
  METNO.FROST <- as.station(df,
    stid=stid,
    loc=loc,
    param=param,
    quality=qualities,
    cntr=cntr,
    lon=lon, lat=lat, alt=alt,
    src='METNO.FROST',
    url=url,
    longname=param1info$longname,
    unit=param1info$unit,
    aspect="original",
    reference="Frost API (http://frost.met.no)",
    info="Frost API (http://frost.met.no)"
  )

  attr(METNO.FROST,'history') <- c(match.call(),date())
  attr(METNO.FROST,'history') <- history.stamp(METNO.FROST)

  invisible(METNO.FROST)
}

