## Description  : Retrieving meta and data from the GHCNM data source 
## Author 	: Abdelkader Mezghani (MET)
## Created 	: 22-03-2013 
## Last update	: 14-06-2013
## Functions    : test.ghcnm(completed) ; t2m.GHCNM(completed) ; metaghcnm(completed) ; dataghcnm(completed)

## TEST FUCNTION
                                        #  source("~/SHARED/R_scripts/esd/R/ghcnm.R")
test.ghcnm <- function(silent=FALSE) {
  stid <- 89001000
  obs <- t2m.ghcnm(stid = stid , verbose = verbose)
  str(obs)
  plot(obs)
}

## MAIN FUNCTION
                                        # Get and format output data for esd further processing
t2m.GHCNM <- function(stid = "10160403000", param = "TAVG", src = "ghcnm", ver = "v3.2.0", adj = TRUE, path = "/klimadata/work/abdelkaderm/data/", force = TRUE, flag = FALSE, verbose = FALSE) {
                                        # Get the meta data
  m <- ghcnm.meta(stid = stid , verbose = verbose)
  x <- ghcnm.data(stid = stid , ver = ver , verbose = verbose)
                                        # convert units
  x[,5:16] = x[,5:16]/100
                                        # extract data
                                        # x = data.frame(subset(x,select=c(-1:-3))) # remove duplicate meta info
                                        # create zoo object
  z <- zoo(cbind(x[5:16]),order.by = as.numeric(rownames(table(x$YEAR))))
  
  attr(z, "station_id") <- m$station_id
  attr(z, "variable")	<- "t2m"
  attr(z, "longname") 	<- "Temperature"
  attr(z, "unit") 	<- "degC"
  attr(z, "quality") 	<- "adjusted"
  attr(z, "longitude") 	<- m$longitude
  attr(z, "latitude") 	<- m$latitude
  attr(z, "altitude") 	<- m$altitude
  attr(z, "frequency") 	<- 1
  attr(z, "calendar") 	<- "gregorian"
  attr(z, "country") 	<- m$country
  attr(z, "location") 	<- m$location
  attr(z, "source") 	<- "GHCNM"
  attr(z, "version")  	<- attr(m,"version") 
  attr(z, "url") 	<- paste(url,ver,sep="/")
  attr(z, "history") 	<- c("ghcnm.data.rda","ghcnm.meta.rda") 
  attr(z, "date-stamp") <- date()
  attr(z, "inv_filename") <- attr(m,"inv_filename")
  attr(z, "call") 	<- match.call()

  class(z) <- c("station","monthly")
  
  invisible(z)

## obsfinalformat <- list(PARAM = rownames(table(obs.data$ELEMENT)),ID = paste(obs.meta$COUNTRY.CODE,obs.meta$STN,sep=""),NAME = obs.meta$NAME, LAT =obs.meta$LATITUDE,LON =obs.meta$LONGITUDE,ELEV = obs.meta$STNELEV,EXTRA = obs.meta$EXTRA,VAL = obs,src = upper.case(src), OR_VERSION = version.por , UPD_VERSION = version.upd, ADJUSTED = upper.case(adj), INV_FILE = invfile, DATA_FILE = fdata, ftp_src_link = "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn")
}

## SUB-FUNCTION "metaghcnm"
ghcnm.meta <- function(stid = NULL, ele = 101,src = "ghcnm",ver = "v3",adj = "qca", path = "data.GHCNM",url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/",force = TRUE,verbose=FALSE) {
  ## SETTING AND/OR CREATING WORK DIRECTORY
  ## AM 02.09.2013 This function could not handle several variables in the same time, each variable should be loaded seperately. i.e. the argument ele should be one single value and not a vector of values !!! 
  
  path <- file.path(path,ver,fsep= .Platform$file.sep) 
  if (verbose) print(paste("Working Directory is :",path, sep=" -> "))
  if (!file.exists(path)) {
    test <- readline("Directory does not exist ! Would you like to create it (yes or no)")
    if (tolower(test) == "yes") dir.create(path,recursive=TRUE)
  }
  ##setwd(newpath)
                                        # SET PARAMETER
  param = tolower(as.character(ele2param(ele=ele,src=toupper(src))[1,5]))
  if (verbose) print("Reading meta data ...")
  if (verbose) print(paste("Param", param,sep=" <<-->> "))
                                        #browser()     
  if (tolower(src) == "ghcnm"){
    if (verbose) print(paste("src",src))
    if (!file.exists("ghcnm.meta.rda") | force) {
     # newpath2 <- paste(path,src,sep = "")   
     # setwd(newpath2)
                                #destpath = paste(newpath,src,"/",sep="")
      destfile <- zipfile <- paste(src,param,"latest",adj,"tar","gz",sep=".")
                                        #destfile = paste(zipfile,sep="") 
      url <- paste(url,ver,"/",sep="")
      fullurl = paste(url,zipfile,sep="")
      if (!file.exists(destfile)) {  
        if (verbose) print("Downloading metadata from 'ftp.ncdc.noaa.gov'")    	
        download.file(fullurl, destfile, method = "wget", quiet = FALSE, mode = "w",cacheOK = TRUE, extra = getOption("download.file.extra"))
      }
      files <- untar(destfile,list="TRUE")
      inv.file <- files[grep("inv",files)]
      dat.file <- files[grep("dat",files)]
      upd.ver <- substr(inv.file, nchar(inv.file)-22,nchar(inv.file)-8)
      untar(destfile, files = c(inv.file,dat.file), list = FALSE,exdir =path)	
                                        # Cleaning the file from "#" character
      text <- readLines(inv.file,encoding="US-ASCII")
      text1 <- gsub("#",replacement="X",x=text)
      writeLines(text1,inv.file)							
      	
      if (verbose) print(paste("Version",upd.ver))
      if (verbose) print(paste("Adjusted", adj))
                                        # Updating the filename	
      #invfile <- paste(src,param,version,adj,"inv",sep=".")  	
      #setwd(newpath)  	
                                        #finventory <- paste(newpath,invfile,sep="")
                                        # Reading meta data ...
      meta <- read.fwf(inv.file,widths=c(11,9,10,8,30,38),col.names=c("stid","lat","lon","alt","loc","extra"),sep = "\t",as.is=TRUE)
      #attach(meta,warn.conflicts=FALSE)
                                        #if (!is.null(stid)) { 
                                        #ID <- paste(ghcnm.meta$COUNTRY.CODE,ghcnm.meta$STN,sep="")	  
                                        #ghcnm.meta <- subset(ghcnm.meta,ID==stid)
                                        #}
                                        # extract country code
      cntr <- substr(meta$stid,1,3)
             				# Generate country full name
      cc <- as.integer(table(factor(cntr)))
      meta$country <- rep(ghcn.country.code()$name[is.element(ghcn.country.code()$code,levels(factor(cntr)))],cc)
                                        # remove consecutive blanks from location name
      meta$loc <- gsub(pattern="  ",x=meta$loc,replacement="")
                                        # remove "extra" coloumn (number 6) from data frame
      meta <- meta[,-6]
                                        
      txt <- readLines(con=dat.file,encoding="US-ASCII")
      stid <- substr(txt,1,11)
      year <- substr(txt,12,15)
      cc <- as.integer(table(stid))
      # get start (first) and end (last) years      
      start <- year[c(1,cumsum(cc)[1:length(cc)-1]+1)] 
      end <- year[cumsum(cc)]
      ele <- rep(ele,length(meta$stid))
      
      ghcnm.meta <- list(station_id =meta$stid , location = as.character(meta$loc) , country = as.character(meta$country) , longitude = as.numeric(meta$lon) , latitude = as.numeric(meta$lat) ,  altitude = as.numeric(meta$alt) , element = as.numeric(ele) , start = as.integer(start) , end = as.integer(end))
                                        # Add attributes to dataframe  	
      attr(ghcnm.meta,"source") <- src
      attr(ghcnm.meta,"version") <- upd.ver	
      attr(ghcnm.meta,"history") <- c(inv.file,dat.file)
      attr(ghcnm.meta,"URL") <- paste(url,ver,sep="/")
      attr(ghcnm.meta,"cite") <- "J. H. Lawrimore, M. J. Menne, B. E. Gleason, C. N. Williams, D. B. Wuertz, R. S. Vose, and J. Rennie (2011), An overview of the Global Historical Climatology Network monthly mean temperature data set, version 3, J. Geophys. Res., 116, D19121, doi:10.1029/2011JD016187." 
      attr(ghcnm.meta, "file") <- "ghcnm.meta.rda"
      attr(ghcnm.meta, "call") <- match.call()
                                        # save in rda file
      save(ghcnm.meta,file="ghcnm.meta.rda")
    } else {
      if (verbose) print("Reading Meta data from rda file ...")
      load("ghcnm.meta.rda")
      if (verbose) print("Done !")
    }
    ## reset to home directory
    setwd("~")
    return(ghcnm.meta)
  }
  
} # END OF metaghcnm function

                                        # SUB-FUNCTION "dataghcnm"
ghcnm.data <- function(ele = 101, stid = "10160403000", ver = "v3", adj = "qca", src = "ghcnm", path = "data.GHCNM",url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn",flag = FALSE, force = FALSE, verbose=FALSE) {

  ## oldpath <- getwd()
  ## path <- file.path(path,ver,fsep= .Platform$file.sep) 
  ## convert ele to param
  param = tolower(as.character(ele2param(ele=ele,src="GHCNM")[5]))
  ## Print work directory
  if (verbose) print(paste("Working Directory ",path, sep=" -> "))
  ## browser()                                           # check path exist or not, otherwise, create it
  if (!file.exists(path)) {
      test <- readline(paste("No such directory",path,"! Would you like to create it (y or n)",sep =" "))
      if ((tolower(test) == "y") | (tolower(test) == "yes")|(tolower(test) == "ye")) {dir.create(path,recursive=TRUE)}
      else return(NULL)
  }
  ## set work directory to new path
  ## setwd(newpath)
  ## print comments
  if (verbose) print(paste("PARAM", param ,sep=" <<-->> "))
  if (verbose) print(paste('source',src))
  
  ## if (!file.exists("ghcnm.data.rda") | force == TRUE) {
  destfile <- file.path(path,paste(src,param,"latest",adj,"tar","gz",sep="."))
  
  gzip.name <- paste(src,param,"latest",adj,"tar","gz",sep=".")
  gzip.path <- file.path(path,gzip.name,fsep = .Platform$file.sep)
  ## browser()
  ## Check if file exists 
  if (!file.exists(gzip.path) | force) {  
    url = paste(url,ver,gzip.name,sep="/")
    if (verbose) print(paste("Reading data from",url))    	
    test <- try(eval(download.file(url, gzip.path, method = "wget", quiet = FALSE, mode = "w", cacheOK = TRUE, extra = getOption("download.file.extra"))))
    if (test>0) {
      if (verbose) print(paste(destfile,"does not exist",sep=" "))
      return(NULL)
    }
  }
  ## KMP 2016-09-16: This part has to be done regardless of if gzip.path exists or not.
  #} else {
      untar(gzip.path, list = FALSE, exdir=path) ## files = c(inv.file,gzip.path),
      ## do we really have to untar for every station? we could search for file first
      files <- untar(gzip.path,list="TRUE")
      inv.file <- files[grep("inv",files)]
      data.file <- files[grep("dat",files)]
      upd.ver <- substr(inv.file, nchar(inv.file)-22,nchar(inv.file)-8)
      ## untar(destfile, files = c(inv.file,data.file), list = FALSE,exdir=file.path(path,ver))	
      ## Cleaning the file from comment "#" character
      inv.path <- file.path(path,paste(src,upd.ver,sep="."),fsep= .Platform$file.sep)
      inv.file <- file.path(inv.path,paste(src,param,upd.ver,adj,"inv",sep="."),fsep= .Platform$file.sep)           
      text <- readLines(inv.file,encoding="US-ASCII")
      text1 <- gsub("#",replacement="X",x=text)
      writeLines(text1,inv.file)
      data.path <- file.path(path,paste(src,upd.ver,sep="."),fsep= .Platform$file.sep)
      data.file <- file.path(data.path,paste(src,param,upd.ver,adj,"dat",sep="."),fsep= .Platform$file.sep)
  #}
  if (verbose) print(paste("Version",upd.ver))
  if (verbose) print(paste("Adjusted", adj))
  ## Reading data as text ...
  ## setwd(newpath)		
  ## fdata <- paste(src,param,upd.version,adj,"dat",sep=".")  	
  ## fdata <- "ghcnm.tavg.v3.2.0.20130120.qca.dat"
  datatext = readLines(data.file)
  if (!is.null(stid)) datatext <- datatext[grep(stid,datatext)]       
  writeLines(datatext,file.path(path,"station.tmp",fsep= .Platform$file.sep))
  ghcnm.data <- read.fwf(file.path(path,"station.tmp",fsep= .Platform$file.sep),widths=c(3,8,4,4,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3), col.names= c("COUNTRY.CODE","ID","year","element","JAN","FLAG1","FEB","FLAG2","MAR","FLAG3","APR","FLAG4","MAY","FLAG5","JUI","FLAG6","JUL","FLAG7","AUG","FLAG8","SEP","FLAG9","OCT","FLAG10","NOV","FLAG11","DEC","FLAG12"),sep = "\t",as.is=TRUE)   
  if (!flag) ghcnm.data <- subset(ghcnm.data,select = c(seq(-6,-28,-2)))	
  ## Replacement of -9999 by "NA"		
  ghcnm.data[ghcnm.data == -9999] <- NA	
  ## Add attributes
  attr(ghcnm.data, "version") <- ver
  attr(ghcnm.data,"inv_file") <- destfile
  attr(ghcnm.data,"url") <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/",ver,sep="/")
  
  ##ghcnm.data <- read.fwf(invfile,widths=c(3,8,9,10,8,30,38),col.names=c			("COUNTRY.CODE","STNR","LATITUDE","LONGITUDE","STNELEV","NAME","EXTRA"),sep = "\t",as.is=TRUE)
  ## save(ghcnm.data,file="ghcnm.data.rda")
                                        # Remove temporary files
  file.remove(file.path(path,"station.tmp",fsep= .Platform$file.sep))
  ## } else {
  ##  if (verbose) print("Reading data from R data file ...")
  ##  load("ghcnm.data.rda")
  ##  if (verbose) print("Done !")
  ##}
  ## reset to hmoe directory
  ## setwd(oldpath)
  return(ghcnm.data)
}
## END OF FUCNTION dataghcnm
