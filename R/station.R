## Author A. Mezghani git-test
## Can return one station or a group of stations. A cluster of regional stations can be
## reppresented in terms of PCA, maintaining the spatial consistencies between them.
## It is also possible to use a set of station objects to provide multi-information from
## one location: e.g. T(2m), precip; mean, variance, skewness, & kurtosis...

## Major updates : 29.07.2013 ; 29.08.2013 ; 03.09.2013 ; 25.09.2013 ; 18.10.2013
## station.metno(ok) ; station.nordklim(ok) ; station.nacd(ok) ; station.ecad(ok) ; station.narp(in progress) ; station.ghcnm(ok) ; station.ghcnd(almost done - checking for t2m)  
## require(zoo)

## ecad (updated) , 

## This function is used to check wether there are errors in the programming !
test.station <- function(ss=NULL,stid=NULL,alt=NULL,lat=c(50,70),lon=c(0,30),param="precip",src=c("GHCND","GHCNM","NORDKLIM","NACD","METNOM","METNOD","ECAD"),verbose=FALSE) {
  for (i in 1:length(src)) {
    y <- station(stid=stid,alt=alt,lat=lat,lon=lon,param=param,src=src[i],nmin=100,verbose=verbose)
    print(summary(y))
  }
}

## Define method
station <- function(stid=NULL,...) UseMethod("station")

station.ecad <- function(...) {
  ## 
  y <- station(src="ecad",...)
  ## attr(y,"call") <- match.call()
  invisible(y)
}

station.ghcnd <- function(...) {
  ## 
  y <- station(src="ghcnd",...)
  ## attr(y,"call") <- match.call()
  invisible(y)
}

station.nacd <- function(...) {
  ## 
  y <- station(src="nacd",...)
  ## attr(y,"call") <- match.call()
  invisible(y)
}

station.narp <- function(...) {
  ## 
  y <- station(src="narp",...)
  ## attr(y,"call") <- match.call()
  invisible(y)
}

station.nordklim <- function(...) {
  ## 
  y <- station(src="nordklim",...)
  invisible(y)
}

station.metnom <- function(...) {
  ## 
  y <- station(src="metnom",...)
  ## attr(y,"call") <- match.call()
  invisible(y)
}

station.metnod <- function(...) {
  ## 
  y <- station(src="metnod",...)
  ## attr(y,"call") <- match.call()
  invisible(y)
}

station.ghcnm <- function(...) {
  ## 
  y <- station(src="ghcnm",...)
  ## attr(y,"call") <- match.call()
  invisible(y)
}

station.daily <- function(src=NULL,...) {
  
  SRC <- c("METNOD","ECAD","GHCND")  
  if (is.null(src)) src <- SRC else src <- intersect(tolower(src),tolower(SRC))
  y <- station(src=src,...)
  attr(y,"call") <- match.call()
  invisible(y)
  
}

station.monthly <- function(src=NULL,...) {
  
  SRC <- c("NACD","NARP","NORDKLIM","METNOM","GHCNM")
  if (is.null(src)) src <- SRC else src <- intersect(src,SRC)
  y <- station(src=src,...)
  attr(y,"call") <- match.call()
  invisible(y)
  
}

## default / retrieve several stations into a zoo object from one or several data sources
station.default <- function(loc=NULL, param='t2m',src = NULL, path=NULL, qual=NULL,url = NULL,
                            stid=NULL, lon=NULL, lat=NULL, alt=NULL, cntr=NULL,
                            it= NULL,nmin=NULL, plot=FALSE,verbose=FALSE,
                            path.ecad=NULL,url.ecad=NULL,
                            path.ghcnm=NULL,url.ghcnm=NULL,
                            path.ghcnd=NULL,url.ghcnd=NULL,
                            path.metnom=NULL,url.metnom=NULL,
                            path.metnod=NULL,url.metnod=NULL) {
  ##
  ## check wether x is a 'location' or a 'stationmeta' object
  
  if (inherits(loc,"stationmeta")) {
    ss <- loc
  } else if (is.character(loc)) {
    loc <- loc
    ss <- NULL
  } else ss <- NULL
  ## else stop("x must be either a stationmeta or location object")
  
  ## Initialize X
  X <- NULL
  SRC <- src
  PATH <- path
  URL <- url
  
  ## Select one or a set of stations based on the metadata
  if (is.null(ss)) { 
    ss <- select.station(stid=stid,loc=loc,lon=lon,lat=lat,alt=alt,cntr=cntr,param=param,src=src,it=it,nmin=nmin) # AM-29.07.2013 "loc" added into the arguments 
  } 
  ## 
  ##if (!is.null(ss)) {
  ##  SRC <- src
  ##}
  if ((param=="t2m") & is.null(ss)) {
    param0 <- param
    ssn <- select.station(param="tmin",stid=stid,loc=loc,lon=lon,lat=lat,alt=alt,cntr=cntr,src=src,it=it,nmin=nmin)
    ssx <- select.station(param="tmax",stid=stid,loc=loc,lon=lon,lat=lat,alt=alt,cntr=cntr,src=src,it=it,nmin=nmin)
    if (!is.null(ssn) & !is.null(ssx))
      class(ssn) <- class(ssx) <- "data.frame" else
      {print('Found no stations with given criteria'); return(NULL)}
    ss <- subset(ssx,ssx$station_id==ssn$station_id) # keep only stations recording both min and max
    if (is.null(ss))
      return(NULL)
    else {
      rl <- readline("T2m is not available for your selection but TMIN and TMAX have been found - Would you like to continue using the averaged values? (y or n): ")
      if ((rl=="y") | rl==("ye") | (rl=="yes"))
        ss$element <- rep(esd2ele(param),length(ss$station_id)) ## update element with param
      else stop("Process stopped")
    }
  } 
  ## 
  
  ## Update attributes based on selected stations (ss)
  id <- ss$station_id
  if (verbose) {print("Station ID:"); str(id) }
  ## 
  ## Extract attributes of all selected stations
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
  param <- apply(as.matrix(ss$element),1,esd2ele)
  rm(ss)
  
  print(paste("Retrieving data from",length(id),"records ..."))
  
  
  ## start loap on available stations
  for (i in 1:length(id)) {
    ##
    ## For now param0 is only used for GHCND dataset
    if (!is.null(param[i]) & !is.null(src[i])) {
      if ((param[i] == "t2m") & (src[i]=="GHCND")) {
        param0 <- param[i] 
        ## param <- apply(as.matrix(ss$element),1,esd2ele)
      } else param0 <-  NULL
      ## used only for ghcnd#update param value from ss
    }
    
    if (!is.null(param0))
      print(paste(i,toupper(param0),stid[i],loc[i],cntr[i],src[i]))
    else
      print(paste(i,toupper(param[i]),stid[i],loc[i],cntr[i],src[i]))
    
    if (src[i]=="METNOD") { #AM-29.08.2013 added for metno data
      ## 
      if (is.null(path.metnod)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.metnod ## default path
      if (is.null(url.metnod)) url="http://klapp/metnopub/production/" else url <- url.metnod ## default url
      x <- metnod.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],start=start[i],end=end[i],qual=qual[i],param=param[i],verbose=verbose, path=path,url=url) ## ,path=path, url=url
      if (verbose) {print("obs"); str(x)}
      if (sum(is.na(coredata(x)))==length(coredata(x))) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
      ## else if (!is.null(x)) X <- combine.stations(X,x)
    } else if (src[i]=="METNOM") { #AM-29.08.2013 added for metno data
      ## 
      if (is.null(path.metnom)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.metnom ## default path
      if (is.null(url.metnom)) url="http://klapp/metnopub/production/" else url <- url.metnom ## default url
      x <- metnom.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],start=start[i],end=end[i],qual=qual[i],param=param[i],verbose=verbose, path=path,url=url) ## ,path=path, url=url
      if (verbose) {print("obs"); str(x)}
      if (sum(is.na(coredata(x)))==length(coredata(x))) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
      ##if (!is.null(x)) X <- combine.stations(X,x)
    } else if (src[i]=="ECAD") { #AM-29.07.2013 added "|(src[i]=="ECAD")"
      ## 
      if (is.null(path.ecad)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.ecad ## default path
      if (is.null(url.ecad)) url="http://www.ecad.eu/utils/downloadfile.php?file=download/ECA_blend" else url <- url.ecad ## default url
      x <- ecad.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose,path=path, url=url)
      if (verbose) {print("obs"); str(x)}
      if (sum(is.na(coredata(x)))==length(coredata(x))) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
      ## else if (!is.null(x)) X <- combine.stations(X,x)
    } else if (src[i]=="NACD") {
      ##if (is.null(path.nacd)) path <- paste("data.",tolower(src[i]),sep="") ## default path
      ##if (is.null(url.nacd)) url="" ## default url
      x <- nacd.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose)
      ##
      if (verbose) {print("obs"); str(x)}
      if (sum(is.na(coredata(x)))==length(coredata(x))) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
      ## else if (!is.null(x)) X <- combine.stations(X,x) 
    } else if (src[i]=="NARP") {
      ##
      x <- narp.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose)
      if (verbose) {print("obs"); str(x)}
      if (sum(is.na(coredata(x)))==length(coredata(x))) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
      ## else if (!is.null(x)) X <- combine.stations(X,x)  
    } else if (src[i]=="NORDKLIM") {
      ##
      x <- nordklim.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose)
      if (verbose) {print("obs"); str(x)}
      if (sum(is.na(coredata(x)))==length(coredata(x))) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
      ## if (!is.null(x)) X <- combine.stations(X,x)  
      #AM-29.07.2013 added begin
    } else if (src[i]=="GHCNM") {
      ## 
      if (is.null(path.ghcnm)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.ghcnm ## default path
      if (is.null(url.ghcnm)) url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn" else url <- url.ghcnm ## default url
      x <- ghcnm.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose,path = path,url=url)
      if (verbose) {print("obs"); str(x)}
      if (sum(is.na(coredata(x)))==length(coredata(x))) {
        print("Warning : No values found in the time series -> This station will be ignored")
        x <- NULL
      }
      ## if (!is.null(x)) X <- combine.stations(X,x)      
    } else if (src[i]=="GHCND") {
      if (is.null(path.ghcnd)) path <- paste("data.",toupper(src[i]),sep="") else path <- path.ghcnd## default path
      if (is.null(url.ghcnd)) url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all" else url <- url.ghcnd ## default url
      ##
      if (!is.null(param0)) { ## compute the avg 
        ghcnd.tmin <- ghcnd.station(param="tmin",stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],verbose=verbose,path = path,url=url)
        ghcnd.tmax <- ghcnd.station(param="tmax",stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],verbose=verbose,path = path,url=url)
        ## scale.tmin <-as.numeric(ele2param(ele="121",src="GHCND")[3])
        if (is.null(ghcnd.tmax) |  is.null(ghcnd.tmin)) 
          x <- NULL
        else {
          ## scale.tmax <-as.numeric(ele2param(ele="111",src="GHCND")[3])
          ## compute the average values from tmin and tmax
          x <- (ghcnd.tmin + ghcnd.tmax) / 2
          ## copy all attributes
          x <- attrcp(ghcnd.tmin,x)
          class(x) <- class(ghcnd.tmin)
          rm(ghcnd.tmin,ghcnd.tmax)
          print("WARNING : Average temperature values have been computed from TMIN and TMAX values")
          param1 <- "TAVG"
          ele <-  "101"
          ## update attributes variable and ele
          attr(x,'variable') <- param
          ## attr(x,'variable') <- switch(param1,'TAVG'=expression(T[2*m]),'TMAX'=expression(paste("max ",T[2*m])),'TMIN'=expression(paste("min ",T[2*m])))
          attr(x,'element') <- ele
          attr(x,'longname') <- "Mean temperature"
        }
        ## if (!is.null(x)) X <- combine.stations(X,x)
      }
      else {
        x <- ghcnd.station(stid=stid[i],lon=lon[i],lat=lat[i],alt=alt[i],loc=loc[i],cntr=cntr[i],qual=qual[i],param=param[i],verbose=verbose,path = path,url=url)
        if (verbose) {print("obs"); str(x)}
        ## 
        if (sum(is.na(coredata(x)))==length(coredata(x))) {
          print("Warning : No values found in the time series -> This station will be ignored")
          x <- NULL
        }
        ## if (!is.null(x)) X <- combine.stations(X,x)
      }
    }
    
    if (!is.null(x)) {
      if (i==1)
        X <- x 
      else
        X <- combine.stations(X,x)
    }
  }
  ## 
  if (is.null(X)) {return(NULL)} ## ;setwd(oldpath)}
  
  ##if (length(id)>1) colnames(X) <- attr(X,"location") #AM-30.07.2013 added line to update colnames in the zoo object
  ## Plot the locations:
  if (plot & !is.null(X)) map.station(X,col="darkgreen",bg="green", cex=0.7)
  
  if (is.null(stid)) {
    x <- NULL
    print(match.call())
  }
  invisible(X)
}

t2m.ghcnd.avg <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param="t2m", path = path , ver = "v3.2.0.20130120",adj = TRUE,force=FALSE,flag = FALSE, off = FALSE, verbose = FALSE) {
  ## compute the average if param is "t2m" and dispplay a warning 
  ## 
  ghcnd.tmin <- station(param="tmin",stid = stid, src = "GHCND", verbose = verbose , path = path)
  if (is.null(ghcnd.tmin)) return(NULL)
  ## scale.tmin <-as.numeric(ele2param(ele="121",src="GHCND")[3])
  ghcnd.tmax <- station(param="tmax",stid = stid, src = "GHCND", verbose = verbose , path = path) 
  if (is.null(ghcnd.tmax)) return(NULL)
  ## scale.tmax <-as.numeric(ele2param(ele="111",src="GHCND")[3])
  ## compute the average values from tmin and tmax
  ghcnd <- (ghcnd.tmin + ghcnd.tmax) / 2
  ## copy all attributes
  ghcnd <- attrcp(ghcnd.tmin,ghcnd)
  print("WARNING : Average temperature values have been computed from TMIN and TMAX values")
  param1 <- "TAVG"
  ele <-  "101"
  ## update attributes variable and ele
  attr(ghcnd,'variable') <- switch(param1,'TAVG'=expression(T[2*m]),'TMAX'=expression(paste("max ",T[2*m])),'TMIN'=expression(paste("min ",T[2*m])))
  attr(ghcnd,'element') <- ele
  invisible(ghcnd)
}

ecad.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,param=NULL,qual=NULL,path="data.ECAD",url="http://www.ecad.eu/utils/downloadfile.php?file=download/ECA_blend",verbose=FALSE) {  ## it=it,nmin=nmin
  ## ECAD basic function to retrieve data for one station
  ## http://eca.knmi.nl/
  ## ECA&D was initiated by the ECSN in 1998 
  ## Retrieve data from the ECSN (ECA&D) data:
  ## data("station.meta",envir=environment())
  
  ##ss <- select.station(src="ecad",stid=stid,lon=lon,lat=lat,loc=loc,alt=alt,cntr=cntr,qual=qual,param=NULL,it=it,nmin=nmin,verbose=verbose)
  ##if (is.null(ss)) return(NULL)
  ##stid <- ss$station_id
  ##loc <- as.character(ss$location)
  ##lon <- ss$longitude
  ##lat <- ss$latitude
  ##alt <- ss$altitude
  ##cntr <- as.character(ss$country)
  ##ele <- ss$element
  ##qual <- ss$quality
  ##start <- ss$start
  ##end <- ss$end
  ##param <- apply(as.matrix(ss$element),1,esd2ele)
  ##rm(ss)
  ##browser()
  ele <- esd2ele(param=param)
  if (is.null(ele)) 
    param1 <-as.character(ele2param(ele=param,src="ECAD")[5])
  else
    param1 <- as.character(ele2param(ele=ele,src="ECAD")[5])
  
  if (!is.null(param) & (!is.null(dim(param1)[1])))
    if (dim(param1)==0) 
      stop('Please refrech your selection, element not found in meta data')
  
  scale <-as.numeric(ele2param(ele=ele,src="ECAD")[3])
  ##param.gp <- substr(param1,1,nchar(param1)-1)
  fdata <- paste(url,"_",tolower(param1),".zip",sep="") 
  text  <- unlist(strsplit(fdata,split="/"))
  text2 <- text[length(text)]
  destfile <- file.path(path,text2,fsep= .Platform$file.sep)
  text3 <- paste('ECA','blend',tolower(param1),sep='_')
  destfile2 <- file.path(path,text3,fsep= .Platform$file.sep)
  ## 
  ## If zip file exist and not the data folder, then unzip
  if (file.exists(destfile) & !file.exists(destfile2)) {
    unzip(destfile,exdir=substr(destfile,1,nchar(destfile)-4))
  } 
  ## if folder does not exist, then download and unzip
  if (!file.exists(destfile2)) { 
    download.file(fdata,destfile,method = "wget", quiet = FALSE, mode = "w", cacheOK = TRUE, extra = getOption("download.file.extra"))
    unzip(destfile,exdir=substr(destfile,1,nchar(destfile)-4))
  }
  
  if (verbose) print("station.ecad")
  newpath <- substr(destfile,1,nchar(destfile)-4) 
  stid <- gsub(' ','',stid)
  for (i in 1:length(stid)) 
    while(nchar(stid[i]) < 6) stid[i] <- paste('0',stid[i],sep="")
  
  
  fnames <- paste(toupper(param1),'_STAID',stid,'.txt',sep="")
  fnames <- file.path(newpath,fnames,fsep = .Platform$file.sep)
  ##print(fnames)
  ipick <- file.exists(fnames)
  if (sum(ipick)==0)  return(NULL)
  if (sum(ipick)!=1) {
    warning('More than one matches - I choose the first!')
    ipick <- (1:length(ipick))[ipick][1]
  }
  fname <- fnames[ipick]
  if (verbose) print(fname)
  
  x <- read.table(fname,header=TRUE,skip=20,sep=",")
  ##
  eval(parse(text=paste("ecad <- scale * x$",param1,sep="")))
  year <- substr(as.character(x$DATE),1,4);L <- length(year)
  month <- substr(as.character(x$DATE),5,6)
  day <- substr(as.character(x$DATE),7,8)
  
  ecad[ecad < -99] <- NA
  
  if (verbose) {
    print(c(year[1],month[1],day[1]))
    print(c(year[L],month[L],day[L]))
    print(summary(ecad))
    print(length(ecad))
    print(stid)
  }
  
  ECAD <- zoo(ecad, order.by = as.Date(paste(year, month, day, sep = "-"), by='day', length.out = L))
  
  if (sum(ECAD,na.rm=TRUE)==0) {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)} 
  
  ## create a station object
  ECAD <- as.station(ECAD, stid=stid, lon=lon, lat=lat, alt=alt,
                     ## ele=esd2ele(param), freq=1,calendar='gregorian',
                     quality=qual, cntr=cntr, loc=loc, src='ECAD',
                     url="http://eca.knmi.nl/utils/downloadfile.php?file=download/ECA_blend_all.zip", param=param, aspect="original",
                     unit=switch(param1,'TG'='degree Celsius','TX'='deg C','TN'='deg C', 'CC'='oktas','DD'='degrees','FG'='m/s', 'FX'='m/s','HU'='%','PP'='hPa', 'SS'='hours','RR'='mm/day'),
                     longname=as.character(ele2param(ele=ele,src="ECAD")[2]),
                     reference="Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of 20th-century surface air temperature and precipitation series for the European Climate Assessment. Int. J. of Climatol., 22, 1441-1453.",
                     info= "Data and metadata available at http://eca.knmi.nl")
  ## additional attributes
  attr(ECAD,'history') <- c(match.call(),date())
  attr(ECAD,'history') <- history.stamp(ECAD)
  ## class(ECAD) <- c("station","day","zoo")
  invisible(ECAD)
}

nacd.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,verbose=FALSE) {
  ## This R routine reads the NACD data. The code
  ## will not work with the original NACD files: a space
  ## must be inserted between the December value and the
  ## country tag, and the missing values must be changed
  ## from '-9999' to ' -999'.
  ##
  ## Arguments:
  ## 'location' determines the time series.
  ## 'ele.c' determines the element (default=T2m).
  ##
  ## R.E. Benestad
  ## AM 28.08.2013 Has been reviewed and seems to work properly and fast !
  ## AM 26.08.2013 station.nacd() cannot deal with several stations ! So stid should be finite number !
  ## print("station.nacd")
  if (is.null(stid)) return(NULL) ## AM 26.08.2013 added
  
  ele <- esd2ele(param=param)
  #  ele.c<-switch(tolower(param),'t2m'='101','tg'='101','rr'='601','slp'='401','cloud'='801','t2'='101','precip'='601','101'='101','401'='401','601'='601','801'='801')
  
  ## 
  ## load("esd/data/NACD.rda")
  data("NACD")
  loc <- gsub("-",".",loc) # AM replace.char() replaced by gsub()
  loc <- gsub("/",".",loc) # AM replace.char() replaced by gsub()
  if (substr(loc,nchar(loc),nchar(loc)) ==".") loc <- substr(loc,1,nchar(loc)-1) # AM 29.07.2013 added 
  #elem<-switch(ele.c,'101'='t2m','601'='precip','401'='slp','801'='cloud')
  string <- paste("x <- NACD$",loc,".",param,sep="") # AM 27.08.2013 line updated : "elem" is replaced by "ele"
  ##print(string)
  ##
  eval(parse(text=string))
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
  ## print(as.character(meta$V16))
  #  quality<-switch(as.character(quality),
  #                ' H'='Homogenous, rigorously tested & adjusted',
  #                'h'='Homogenous, rigorously tested & adjusted',
  #                ' T'='Tested, maybe adjusted but not perfectly H.',
  #                't'='Tested, maybe adjusted but not perfectly H.',
  #                ' N'='Not tested for inhomogenouity',
  #                'n'='Not tested for inhomogenouity',
  #                ' E'='Environm. changes prevents clim.change studies',
  #                'e'='Environm. changes prevents clim.change studies',
  #                ' I'='Inhomogenous series which presently are unadjustable',
  #                'i'='Inhomogenous series which presently are unadjustable')
  
  ny <- length(attr(x,'year'))
  year <- sort(rep(attr(x,'year'),12)); L <- length(year)
  month <- rep(1:12,length(ny))
  day <- rep(1,length(month))
  x[x < -99] <- NA
  ## 
  NACD <- zoo(c(t(x)),order.by = as.Date(paste(year, month, day, sep = "-"), by='month', length.out = L))
  names(NACD) <- loc # AM 30.07.2013 added           
  #print("attributes")
  # Add meta data as attributes:
  
  unit <- switch(ele,
                 '101'='degree Celsius','111'='degree Celsius',
                 '112'='degree Celsius','122'='degree Celsius',
                 '401'='hPa', '601'='mm','701'='days','801'='%')
  
  NACD <- as.station(NACD, stid=stid, loc=loc, cntr=cntr, lon=lon, lat=lat, alt=alt,
                     ##freq=1, calendar='gregorian',
                     quality=qual, src='NACD', url=NA, param=param,
                     aspect="original", unit=unit, longname=x.name,
                     reference="Frich et al. (1996), DMI scientific report 96-1",
                     info="Data and metadata adopted from clim.pact")
  
  ## Additional attributes
  attr(NACD,'history') <- c(match.call(),date())
  attr(NACD,'history') <- history.stamp(NACD)
  ## class(NACD) <- c("station","month","zoo")
  invisible(NACD)
}


narp.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,verbose=FALSE) {
  ## This R function reads the NARP data.
  ## R.E. Benestad
  ## Major changes by A. Mezghani
  ##print("station.narp")
  ##data("station.meta",envir=environment())
  ##load("esd/data/station.meta.rda")
  ## if (is.null(stid)) return(NULL) 
  ## 
  ## load the original data
  data(NARP,envir=environment())
  ## Get esd ele from param argument.
  ele <- esd2ele(param=param)
  ## 
  x.name <- as.character(ele2param(ele=ele,src="NARP")[2])
  unit <-  as.character(ele2param(ele=ele,src="NARP")[4])
  scale <- as.numeric(ele2param(ele=ele,src="NARP")[3])
  iii <- is.element(NARP[,1],stid) & is.element(NARP[,2],as.numeric(ele))
  ## 
  if (sum(iii) ==0) {
    print("Warning : No recorded values are found for this station -> Ignored")
    return(NULL)}
  x <- NARP[iii,4:15]*scale
  x[x <= -99] <- NA
  ny <- length(NARP[iii,3])
  year <- sort(rep(NARP[iii,3],12))
  L <- length(year)
  month <- rep(1:12,ny)
  day <- rep(1,length(year))
  
  NARP <- zoo(c(t(x)), order.by = as.Date(paste(year, month, day, sep = "-"), by='month', length.out = L))
  
  ## Format as station object
  NARP <- as.station(NARP, stid=stid, loc=loc, lon=lon, lat=lat, alt=alt,
                     ## frequency=1, calendar='gregorian',
                     cntr=cntr, quality=NA, src='NARP', url=NA, param=param,
                     unit=unit, longname=x.name,
                     aspect="original",
                     reference="Nordic Arctic Research Programme",
                     info="narp2esd.R")
  ## Additional attributes
  attr(NARP,'history') <- c(match.call(),date())
  attr(NARP,'history') <- history.stamp(NARP)
  ## class(NARP) <- c("station","month","zoo")
  invisible(NARP)
}

nordklim.station <- function(stid=NULL,loc=NULL,lon=NULL,lat=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,verbose=FALSE,path=NULL) {
  ## This R function reads the NARP data.
  ## R.E. Benestad
  ## Major changes by A. Mezghani
  ##print("station.narp")
  ##data("station.meta",envir=environment())
  ##load("esd/data/station.meta.rda")
  ## 
  if (is.null(stid)) return(NULL)
  ## load the data - should be replaced by the comment below
  ## load("nordklim.data.rda") ## data(nordklim.data,envir=environment())
  data(nordklim.data,envir=environment())
  x <- nordklim.data ; rm(nordklim.data)
  station_id <- stid
  x <- subset(x, station_id==stid)
  ele <- esd2ele(param) 
  x <- subset(x,element==ele)
  # 
  if (dim(x)[1] < 1) {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)}
  ## set missing values to NA
  x[x < -999] <- NA
  ## convert values to deg. C and get data to x1
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
  ## class(NORDKLIM) <- c("station","month","zoo")
  invisible(NORDKLIM)
}

ghcnm.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,ver="v3",path="data.GHCNM",url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn",adj = "qca",force=FALSE,flag = FALSE, verbose = FALSE) {
  ## AM-29.07.2013 "path" argument set to new folder  
  ## /function
  ## GHCNM was initiated by the  
  ## Retrieve data from the GHCNM data:
  ## data("station.meta",envir=environment())
  ## 
  
  ## Convert to esd element
  ele <-esd2ele(param=param) 
  
  ## extract param1 and scale variables
  param1 <-as.character(ele2param(ele=ele,src="GHCNM")[5])
  scale <-as.numeric(ele2param(ele=ele,src="GHCNM")[3])
  
  if (verbose) print("station.GHCNM")
  
  ghcnm <- ghcnm.data(ele=ele,stid = stid, src = "ghcnm", ver = ver , adj = adj, path = path, url=url,force = force, flag = flag, verbose = verbose)
  x <- c(t(ghcnm[,5:16]))*scale
  ## 
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
  ## class(GHCNM) <- c("station","month","zoo")
  invisible(GHCNM)
}

ghcnd.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,param=NULL,path="data.GHCND", url=NULL,adj = TRUE,force=FALSE,flag = FALSE, off = FALSE, verbose = FALSE) {
  ## AM-29.07.2013 "path" argument set to new folder  
  ## /function
  ## GHCND was initiated by the  
  ## Retrieve data from the GHCND data:
  ##  data("station.meta",envir=environment())
  ## depends on ghcnd.data(), 
  
  ## Convert param to esd element
  ele <-esd2ele(param=param) 
  param1 <-as.character(ele2param(ele=ele,src="GHCND")[5])
  
  if (verbose) print("station.GHCND")
  
  scale <-as.numeric(ele2param(ele=ele,src="GHCND")[3])
  ## Get the data
  ghcnd <- ghcnd.data(param = param1,stid = stid, src = "ghcnd", adj = adj, path = path, url=url, force = force, flag = flag, verbose = verbose,rm.file=FALSE)
  
  if (is.null(ghcnd)) return(NULL)
  ## 
  
  ## Convert unit by scale factor 
  x <- c(t(ghcnd[,6:36]))*scale
  
  ## Time vector
  day_id <-substr(names(ghcnd[,6:dim(ghcnd)[2]]),4,nchar(names(ghcnd[,6:dim(ghcnd)[2]]))) 
  year <- rep(as.character(ghcnd$YEAR),each=length(day_id)) ; L <- length(year)
  month <- rep(ghcnd$MONTH,each=length(day_id))
  day <- rep(day_id,dim(ghcnd)[1])
  
  ## Remove erroneous dates
  vdate <- as.Date(paste(year, month, day, sep = "-"), by='day', length.out = L)
  id <- !is.na(vdate)
  vdate <- vdate[id]
  x <- x[id]
  
  if (is.null(x)) return(NULL)
  
  ## format results as a zoo object 
  GHCND <- zoo(x,order.by = vdate)
  
  ## GHCND <-GHCND[!is.na(GHCND)] 
  if (sum(GHCND,na.rm=TRUE)==0) {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)} 
  
  ## Format the results as station object
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

## The following three functions are proper to METNO data source. The default function is station.metno() 
## and two derivative functions that are station.metnod() and station.metnom() for daily and monthly velues, respectively.
## Author : A. Mezghani
## adapted from stnr() function 
metnom.station <-  function(re=15,stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,qual=NULL,
                            start=NULL,end=NULL,param=NULL,verbose=FALSE, h = NULL, nmt = 0,
                            path = NULL, dup = "A", url = "http://klapp/metnopub/production/") {
  
  y <- metno.station(re=re,stid=stid,lon=lon,lat=lat,loc=loc,alt=alt,cntr=cntr,qual=qual,
                     start=start,end=end,param=param,verbose=verbose,h = h, nmt = nmt,
                     path = path, dup = dup, url = url)
  if (!is.null(y))
    attr(y,"source") <- "METNOM"
  invisible(y)
}
metnod.station <-  function(re=14, ...) {
  ## 
  y <- metno.station(re=re,...)
  if (!is.null(y))
    attr(y,"source") <- "METNOD"
  invisible(y)
}
metno.station <- function(stid=NULL,lon=NULL,lat=NULL,loc=NULL,alt=NULL,cntr=NULL,
                          qual=NULL,start=NULL,end=NULL,param=NULL,verbose=FALSE,
                          re = 14,h = NULL, nmt = 0,  path = NULL, dup = "A",
                          url = "http://klapp/metnopub/production/") {
  
  if (verbose) print("http://eklima.met.no")
  ## 
  ## if (!is.na(end)) end1 <- format(Sys.time(),'%d.%m.%Y')
  if (!is.na(end)) end1 <-format(as.Date(paste("31.12.",as.character(end),sep=""),format='%d.%m.%Y'),'%d.%m.%Y')
  if (!is.na(start)) start1 <- format(as.Date(paste("01.01.",as.character(start),sep=""),format='%d.%m.%Y'),'%d.%m.%Y')
  if (!is.null(url)) {
    Filnavn <- paste(url, "metno?re=", re, "&ct=text/plain&del=space&ddel=dot&nod=NA&split=1", sep = "")
    if (!is.null(h)) 
      Filnavn <- paste(Filnavn, "&h=", h, sep = "")
    param1 <- ele2param(ele=esd2ele(param),src="metno")$param ## switch(param,"t2m"="TAM","precip"="RR")
    for (i in 1:length(param1)) Filnavn <- paste(Filnavn,"&p=", param1[i], sep = "")
    if (!is.null(h)) 
      Filnavn <- paste(Filnavn, "&nmt=", nmt, sep = "")
    Filnavn <- paste(Filnavn, "&fd=", start1, "&td=", end1, sep = "")
    Filnavn <- paste(Filnavn, "&s=", stid, sep = "")
    ##for (i in 1:length(StNr)) Filnavn <- paste(Filnavn, "&s=", StNr[i], sep = "")
    if (!is.null(h)) 
      Filnavn <- paste(Filnavn, "&dup=", dup, sep = "")
  } else stop("The url must be specified")
  
  if (verbose) print(Filnavn)
  
  firstline <- readLines(Filnavn, n = 1, encoding = "latin1")
  if (substr(firstline, 1, 3) == "***") {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)}
  ## 
  Datasett <- as.list(read.table(Filnavn,dec = ".", header = TRUE, as.is = TRUE, fileEncoding = "latin1"))
  Datasett$RR[Datasett$RR == "."] <- "0"
  eval(parse(text = paste("y <- as.numeric(Datasett$", param1, ")",sep = "")))
  if (sum(y,na.rm=TRUE)==0) {print("Warning : No recorded values are found for this station -> Ignored") ; return(NULL)}
  
  type <- switch(re, `14` = "daily values", `17` = "observations", '15' = "Monthly means")
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
  invisible(METNO)
}



stnr <- function (navn = NULL, lon = NULL, lat = NULL, max.dist = 10, 
                  alt = NULL, Fylke = NULL, Kommune = NULL, fy = NULL, ty = NULL, 
                  ny = NULL, param = "TAM", plot = FALSE, print = FALSE) 
{
  met.no.meta <- MET.no.meta(param = param, print = print)
  iue <- nchar(met.no.meta$TODATE) == 2
  met.no.meta$TODATE[iue] <- now()
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
    ii <- grep(upper.case(navn), met.no.meta$Navn)
    print(navn)
    print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii]))
  }
  if (plot) {
    plot(c(0, 32), c(57, 73), type = "n", xlab = "lon", ylab = "lat")
    addland()
    points(met.no.meta$Lon, met.no.meta$Lat, col = "grey", cex = 0.8)
  }
  if (xor(is.null(lon), is.null(lat))) 
    stop("both or none of lon/lat must be specified")
  if (!is.null(lon)) {
    if (length(lon) == 1) {
      if (plot) 
        points(lon, lat, pch = "+", col = "blue", cex = 0.7)
      d <- round(distAB(lon,lat,met.no.meta$Lon, met.no.meta$Lat)/1000, 
                 3)
      print(length(d))
      ii <- II[(d <= max.dist)]
      print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii], 
                  met.no.meta$Lon[ii], met.no.meta$Lat[ii], d[ii]))
    }
    else if (length(lon) == 2) {
      if (plot) 
        polygon(c(lon[1], lon[2], lon[2], lon[1], lon[1]), 
                c(lat[1], lat[1], lat[2], lat[2], lat[1]), 
                border = "blue", lwd = 2)
      if (length(lat) == 1) 
        stop("both or none of lon/lat must have two entries")
      ii <- II[(met.no.meta$Lon >= min(lon)) & (met.no.meta$Lon <= 
                                                  max(lon)) & (met.no.meta$Lat >= min(lat)) & (met.no.meta$Lat <= 
                                                                                                 max(lat))]
      print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii], 
                  met.no.meta$Lon[ii], met.no.meta$Lat[ii]))
      met.no.meta <- met.no.meta[ii, ]
    }
  }
  if (!is.null(alt)) {
    if (length(alt) == 1) {
      if (alt > 0) 
        ii <- (met.no.meta$Hoh >= alt) & is.finite(met.no.meta$Hoh)
      else ii <- (met.no.meta$Hoh <= abs(alt)) &
          is.finite(met.no.meta$Hoh)
    }
    else ii <- (met.no.meta$Hoh >= min(alt)) & (met.no.meta$Hoh <= 
                                                  max(alt))
    ii[is.na(met.no.meta$Stnr[ii])] <- FALSE
    print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii],
                met.no.meta$Hoh[ii]))
    met.no.meta <- met.no.meta[ii, ]
  }
  if (!is.null(Fylke)) {
    ii <- is.element(upper.case(met.no.meta$Fylke), upper.case(Fylke))
    print(rbind(met.no.meta$Navn[ii], met.no.meta$Stnr[ii],
                met.no.meta$Fylke[ii], 
                met.no.meta$Kommune[ii]))
    met.no.meta <- met.no.meta[ii, ]
  }
  if (!is.null(Kommune)) {
    ii <- is.element(upper.case(met.no.meta$Kommune), upper.case(Kommune))
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


replace.char <- function (c, s, ny.c)  {
  if (c == ny.c) return(s)
  nc <- nchar(c); ns <- nchar(s)
  is <- 1
  tries <- instring(c, s)
  if (length(tries)==0) return(s)
  #print(nc)
  while ( (instring(c, s)[1] > 0) & (is <= length(tries)) ) {        
    ii <- instring(c, s)[1]
    #print(ii); print(c);print(s)
    if (ii > 1) {
      s <- paste(substr(s, 1, ii - 1), ny.c,
                 substr(s, ii + nc, nchar(s)), sep = "")
    } else if (ii==1) s <- paste(ny.c,
                                 substr(s, ii + nc + 1, nchar(s)), sep = "")
    is <- is + 1
  }
  s
}

