#' Fetch station data from GHND and convert into \code{station} objects.
#' 
#' Transform an input object into the esd class \code{station}. 
#' \code{ghcnd.meta} reads the inventory and metadata of the stations (slow due to big data volume - requires Internet access).
#'
#' \code{ghcnd.station} reads the raw data (slow due to big data volume - requires Internet access).

#' @aliases ghcnd.meta ghcnd.station
#' @seealso station as.station ghcnd.meta ghcnd.station
#' 
#' @param x a 'stationmeta' object
#' @param url lURL of the database
#' 
#' @examples
#' \dontrun{
#' meta <- ghcnd.meta(is=list(cntr=c('Mozambique','Zimbabwe','Malawi')))
#' m <- subset(meta,lon=c(24,35),lat=c(-30,-25),verbose=TRUE)
#' y <- ghcnd.station(subset(m,is=1:10),verbose=TRUE)
#' }
#'
#' @export
ghcnd.meta <- function(is=NULL,url='https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/doc/ghcnd-stations.txt') {
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
  if (!is.null(is)) meta <- subset(meta,is=is)
  return(meta)
}

#' @export
ghcnd.station <- function(x,url='https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/',verbose=FALSE) {
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
        t2m <- as.station(zoo(0.1*y$TAVG,order.by=t),loc=x$location[is],lon=x$longitude[is],lat=x$latitude[is],alt=x$altitude[is],
                          stid=x$station_id[is],wmo.is=x$wmo.id,variable='t2m',unit='degC',cntr=x$country[is],
                          src='GHCND',url=csvfile)
        if (is.null(T2m)) T2m <- t2m else T2m <- combine.stations(T2m,t2m)
      }
    }
  }
  return(list(precip=Precip,tmax=Tmax,tmin=Tmin,t2m=T2m))
}

#' @export
ghcnd.data <- function(param="PRCP", stid="ACW00011604", src="ghcnd" , path="data.GHCND",
                       url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all", flag=FALSE, miss2na=TRUE, 
                       force=TRUE, verbose=TRUE, save.file=FALSE, rm.file =TRUE) {
  if(verbose) print("ghcnd.data")
  if(is.null(url)) url <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all"
  src <- tolower(src)

  if (!file.exists(path) & save.file) {
    test <- readline(paste("Directory,",path," does not exist ! Would you like to create it (y or n)",sep=""))
    if ((tolower(test) == "yes") | (tolower(test) == "ye") | (tolower(test) == "y")) dir.create(path)
  } 
  
  if (verbose) print(paste("PARAM", param ,sep=" <<-->> "))
  if (src == "ghcnd"){
    if (verbose) print(paste("Source",src))     
    ## set destination full path filename to destfile
    file <- paste(stid,"dly",sep=".")
    if (!save.file) path <- url
    destfile <- file.path(path, file, fsep=.Platform$file.sep)
    ##
    if (!file.exists(destfile) | (file.info(destfile)$size==0) | !save.file) {  
      if (verbose) print("Reading data from ftp.ncdc.noaa.gov")
      url = paste(url,file,sep="/")
      if (!save.file) {
        destfile <- url(url)
      } else {
        test <- try(eval(download.file(url, destfile, method = "wget", quiet=FALSE, 
                                       mode="w", cacheOK=TRUE, extra=getOption("download.file.extra"))))
        if (test>0) {return(NULL)} ##;  setwd(oldpath)}
      }
    }
    ##
    if (save.file) {
      if (file.info(destfile)$size==0) {
        print(paste("Warning : File",destfile,"has null size",sep=" "))
        return(NULL)
      }
    }

    ## Reading data as text ...
    ## setwd(newpath)	
    ##	
    ## fdata <- paste(stid,"dly",sep=".")  	
    ##fdata <- "ghcnd.tavg.v3.2.0.20130120.qca.dat"   
    ## if (!is.null(stid) & save.file) datatext = readLines(destfile) ##readLines(fdata)
    ghcnd.data <- try(read.fwf(destfile, widths=c(3,8,4,2,4,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,
                                              5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3),
                           col.names=c("COUNTRY.CODE","ID","YEAR","MONTH","ELEMENT","DAY1","MQSFLAG1","DAY2","MQSFLAG2",
			                      "DAY3","MQSFLAG3","DAY4","MQSFLAG4","DAY5","MQSFLAG5","DAY6","MQSFLAG6","DAY7","MQSFLAG7",
				                    "DAY8","MQSFLAG8","DAY9","MQSFLAG9","DAY10","MQSFLAG10","DAY11","MQSFLAG11","DAY12","MQSFLAG12",
				                    "DAY13","MQSFLAG13","DAY14","MQSFLAG14","DAY15","MQSFLAG15","DAY16","MQSFLAG16","DAY17",
				                    "MQSFLAG17","DAY18","MQSFLAG18","DAY19","MQSFLAG19","DAY20","MQSFLAG20","DAY21","MQSFLAG21",
				                    "DAY22","MQSFLAG22","DAY23","MQSFLAG23","DAY24","MQSFLAG24","DAY25","MQSFLAG25","DAY26",
				                    "MQSFLAG26","DAY27","MQSFLAG27","DAY28","MQSFLAG28","DAY29","MQSFLAG29","DAY30","MQSFLAG30",
				                    "DAY31","MQSFLAG31"), sep="\t", as.is=TRUE))
    if (inherits(ghcnd.data,"try-error")) {#} | dim(ghcnd.data)[1]==0) {
      print("Warning : No data found for that station")
      return(NULL)
    }
    ## attach(ghcnd.data,warn.conflicts = FALSE)
    ## 
    ## Remove flags from values
    ghcnd.data <- subset(ghcnd.data,select = seq(-7,-dim(ghcnd.data)[2],-2))
    
    ## Replacement of missing value -9999 by "NA"
    if (miss2na) ghcnd.data[ghcnd.data == -9999] = NA	
    
    ## Scale the values 
    ## ghcnd.data[,5:16] = ghcnd.data[,5:16] / 100 
    # 
    ## Select variable param or element
    id <- ghcnd.data$ELEMENT == param
    ghcnd.data <- ghcnd.data[id,]
    
    ## Save into file
    ## save(ghcnd.data,file="ghcnd.data.rda")
    
    ## Remove downloaded files if necessary to save disk space
    if (rm.file) file.remove(destfile)
    
    #} ## else {if (verbose) print("Reading data from R-data file ...")
    ##      load("ghcnd.data.rda")
    ##      if (verbose) print("Done !")
    ##    }
    ## reset to home directory
    ## setwd(oldpath)
    return(ghcnd.data)
  }
}


## OLD VERSIONS OF GHCND FUNCTIONS: 
# Functions to fetch data from the Global Historical Climatology Network (GHCN) data base
#
# \code{ghcnd.meta} and \code{ghncm.meta} read and organize metadata of daily (ghcnd) and monthly (ghcnm) GHCN data.
# \code{ghncd.data} and \code{ghcnd.data} read daily and monthly mean GHCN data from NOAA (\url{ftp.ncdc.noaa.gov}).
#
# @aliases ghcnd.meta ghcnd.data ghcnm.meta ghcnm.data
#
# @param param climate variable
# @param src source of data
# @param path path to directory where to save data
# @param url url to data source
# @param save.file a boolean; if TRUE save data or metadata
# @param verbose a boolean; if TRUE print information about progress
# @param force a boolean; if TRUE overwrite old file
#
# @export
ghcnd.meta.old <- function(param=NULL, src="ghcnd", path="data.GHCND",
                       url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily",
                       save.file=FALSE, verbose=TRUE, force=TRUE) {
  ## Get old path
  ## oldpath <- getwd()
  if (!verbose) print(paste("Param", param, sep = " <-- "))
  ## Setting Path to Work directory
  ##newpath <- paste(path,src,"/",sep = "")
  if (!verbose) print(paste("Setting Work Directory to ",path, sep=" -> "))
  if (!file.exists(path) & save.file) {
    test <- readline(paste("Directory",path," does not exist ! Would you like to create it (yes or no)",sep=""))
    if ((tolower(test) == "yes") |(tolower(test) == "ye") | (tolower(test) == "y")) dir.create(path)
  }
  if(file.exists(path) & save.file) setwd(path)
  
  if (src == "ghcnd"){
    if (verbose) print(paste("Source",src))
    
    if (!file.exists("ghcnd.meta.rda") | force) {  
      ## Reading data source - version
      destfile <- paste(src,"-version.txt",sep="")
      if (!file.exists(destfile)) {
        url <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/",destfile,sep="")
        download.file(url, destfile, method = "wget", quiet = TRUE, mode = "w", 
                      cacheOK = TRUE, extra = getOption("download.file.extra")) 
      } 
      text <- readLines(destfile)
      tunlist <- unlist(strsplit(text,' '))
      version.por <- tunlist[grep("-por-",tunlist)]
      if (verbose) print(paste("Original Version",version.por, sep = " <-- "))  
      version.upd <- tunlist[grep("-upd-",tunlist)]
      if (verbose) print(paste("Updated Version",version.upd, sep = " <-- "))
      
      ## Reading data source - inventory  
      finventory <- paste(src,"-inventory.txt",sep="")
      ##destfile = paste(source,"-inventory.txt",sep="")
      url = paste(url,finventory,sep="/")
      if (!file.exists(destfile)) {
        download.file(url, finventory, method = "wget", quiet = FALSE, mode = "w",
                      cacheOK = TRUE, extra = getOption("download.file.extra"))
      }
      
      if (verbose) print("Reading metadata...") 	
      #
      meta1 <- read.fwf(finventory,widths=c(11,9,11,4,5,5),
                        col.names=c("stid","lat","lon","param","start","end"), sep = "\t",
                        as.is=TRUE,header=FALSE)
      ## meta1 <- read.fwf(finventory,widths=c(3,8,9,11,4,5,5),col.names=c("country.code","wmo.code","lat","lon","param","start","end"), sep = "\t",as.is=TRUE) 
      meta <- read.fwf("ghcnd-stations.txt",widths=c(11,9,10,7,4,30,4,4,6),
                       col.names=c("stid","lat","lon","alt","state","stnm","gsnflag","hcnflag","wmo_id"), 
                       sep = "\t",as.is=TRUE,header=FALSE)
     
      ## meta <- read.fwf("ghcnd-stations.txt",widths=c(2,1,8,9,10,7,4,30,4,4,6),
      ##                  col.names=c("cntr.abb","cntr.netw.c","stid","lat","lon","alt","state","stnm","gsnflag","hcnflag","wmo_id"), sep = "\t",as.is=TRUE
      
      ## Replace element name by element id 
      meta1$param <- sub("TMIN",meta1$param,replacement="121")
      meta1$param <- sub("TMAX",meta1$param,replacement="111")      
      meta1$param <- sub("PRCP",meta1$param,replacement="601")
      meta1$param <- sub("SNOW",meta1$param,replacement="999")
      meta1$param <- sub("SNWD",meta1$param,replacement="901")                                    
      
      ## remove auxilary elements and keep only core element mentioned above
      id <- !is.na(as.numeric(meta1$param))
      meta1 <- subset(meta1,id)
     
      ## compute the number of stations
      cc <- as.integer(table(factor(meta1$stid)))

      ## extract the country abbreviated name
      meta$cntr.abb <- substr(meta$stid,1,2)

      ## Generate full country name
      cname <- ghcnd.country.abb()$name
      meta$cntr <- cname[as.integer(factor(meta$cntr.abb))]
                                        # remove extra empty spaces
      meta$stnm <- as.character(gsub("  ",x=meta$stnm,replacement=""))
      n <- nchar(meta$stnm)
      id <- (substr(meta$stnm,n,n)==" ")
      meta$stnm[id] <- substr(meta$stnm[id],1,n-1)  
      
      ## Remove meta data of missing stations in the inventory
      meta <- meta[-which(is.element(meta$stid,meta1$stid)==FALSE),]
      ##                                 # esd output format
      ghcnd.meta <- list(station_id = as.character(meta1$stid), 
                         location = rep(as.character(meta$stnm),cc), country = as.character(rep(meta$cntr,cc)) , 
                         longitude = as.numeric(meta1$lon), latitude = as.numeric(meta1$lat),  
                         altitude = as.numeric(rep(meta$alt,cc)), element = as.numeric(meta1$param), 
                         start = as.integer(meta1$start), end = as.integer(meta1$end))
                                        # Add attributes
      attr(ghcnd.meta, "source") <- src  
      attr(ghcnd.meta, "version") <- c(version.por,version.upd)
      attr(ghcnd.meta, "history") <- c("ghcnd-inventory.txt","ghcnd-stations.txt")
      attr(ghcnd.meta, "URL") <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/"
      attr(ghcnd.meta, "file") <- "ghcnd.meta.rda"
      attr(ghcnd.meta, "cite") <- "Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of 20th-century surface air temperature and precipitation series for the European Climate Assessment. Int. J. of Climatol., 22, 1441-1453."
      attr(ghcnd.meta, "date") <- date()
      attr(ghcnd.meta, "call") <- match.call()
                                        # save ouptuts in file
      save(ghcnd.meta,file="ghcnd.meta.rda")
    } else {
      if (verbose) print("Reading Meta data from rda file ...")
      load("ghcnd.meta.rda")
      if (verbose) print("Done !")
    }
    ## reset the path to user path
    ## setwd(oldpath)
    return(ghcnd.meta)
  }
}

