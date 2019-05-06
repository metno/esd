## Decription  : Main Script for GHCND data format handling
## Author      : Abdelkader Mezghani
## Created     : 01-02-2013
## Last update : 13-06-2013
## Functions   : Dailyformat (in progess) ; metaghcnd (completed) ; dataghcnd (completed)
## Comments    : Lines 70 to 75 must be updated with the full list of variables and elements


## Sub function "metaghcnd"
ghcnd.meta <- function(param=NULL, src="ghcnd", adj=TRUE, path="data.GHCND",
                       url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily",
                       save.file=FALSE, test = TRUE,verbose=TRUE , force = TRUE) {
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
  setwd(path)
  
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
      if (verbose) print(paste("Adjusted", adj))
      if (verbose) print("Reading metadata...") 	
      #
      meta1 <- read.fwf(finventory,widths=c(11,9,11,4,5,5),
                        col.names=c("stid","lat","lon","param","start","end"), sep = "\t",
                        as.is=TRUE,header=FALSE)
      ## meta1 <- read.fwf(finventory,widths=c(3,8,9,11,4,5,5),col.names=c("country.code","wmo.code","lat","lon","param","start","end"), sep = "\t",as.is=TRUE) 
      meta <- read.fwf("ghcnd-stations.txt",widths=c(11,9,10,7,4,30,4,4,6),
                       col.names=c("stid","lat","lon","alt","state","stnm","gsnflag","hcnflag","wmo_id"), 
                       sep = "\t",as.is=TRUE,header=FALSE)
     
      ## meta <- read.fwf("ghcnd-stations.txt",widths=c(2,1,8,9,10,7,4,30,4,4,6),col.names=c("cntr.abb","cntr.netw.c","stid","lat","lon","alt","state","stnm","gsnflag","hcnflag","wmo_id"), sep = "\t",as.is=TRUE
      
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
## SUB-FUNCTION "dataghcnd"
ghcnd.data <- function(param = "PRCP", stid = "ACW00011604" , src = "ghcnd" , adj = TRUE, path = "data.GHCND",
                       url="ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all",flag = FALSE, miss2na = TRUE, 
                       force = TRUE , verbose=TRUE , save.file = FALSE, rm.file =TRUE) {
   ## AM 15-11-2013 removed from the argument list ver ="v3.2.0.20130120"
 
  ## Get old path
  ## oldpath <- getwd()
  src <- tolower(src)
  ## set working directory to new path
  ##newpath <- paste(path,paste(src,ver,sep="."),"/",sep="") 
  
  ## if (verbose) print(paste("Setting Work Directory to ",path, sep=" -> "))
  ##  wser()
  if (!file.exists(path) & save.file) {
    ##
    test <- readline(paste("Directory,",path," does not exist ! Would you like to create it (y or n)",sep=""))
    if ((tolower(test) == "yes") | (tolower(test) == "ye") | (tolower(test) == "y")) dir.create(path)
  } 
  
  if (verbose) print(paste("PARAM", param ,sep=" <<-->> "))
  if (src == "ghcnd"){
    if (verbose) print(paste("Source",src))     
    ## set destination full path filename to destfile
    file <- paste(stid,"dly",sep=".")
    if (!save.file) path <- url
    destfile <- file.path(path,file,fsep = .Platform$file.sep)
    ##
    if (!file.exists(destfile) | (file.info(destfile)$size==0) | !save.file) {  
      if (verbose) print("Reading data from ftp.ncdc.noaa.gov")    	
      url = paste(url,file,sep="/")
      if (!save.file) {
        destfile <- url(url)
      } else {
        test <- try(eval(download.file(url, destfile, method = "wget", quiet = FALSE, mode = "w", cacheOK = TRUE, extra = getOption("download.file.extra"))))
        if (test>0) {return(NULL)} ##;  setwd(oldpath)}
      }
    }   	 
    ## 
    ## browser()
    if (save.file) 
      if (file.info(destfile)$size==0) {
        print(paste("Warning : File",destfile,"has null size",sep=" "))
        return(NULL)}

    ## Reading data as text ...
    ## setwd(newpath)	
    ##	
    ## fdata <- paste(stid,"dly",sep=".")  	
    ##fdata <- "ghcnd.tavg.v3.2.0.20130120.qca.dat"   
    ## if (!is.null(stid) & save.file) datatext = readLines(destfile) ##readLines(fdata)    
   
    ghcnd.data <- read.fwf(destfile,widths=c(3,8,4,2,4,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3,5,3),col.names=c("COUNTRY.CODE","ID","YEAR","MONTH","ELEMENT","DAY1","MQSFLAG1","DAY2","MQSFLAG2","DAY3","MQSFLAG3","DAY4","MQSFLAG4","DAY5","MQSFLAG5","DAY6","MQSFLAG6","DAY7","MQSFLAG7","DAY8","MQSFLAG8","DAY9","MQSFLAG9","DAY10","MQSFLAG10","DAY11","MQSFLAG11","DAY12","MQSFLAG12","DAY13","MQSFLAG13","DAY14","MQSFLAG14","DAY15","MQSFLAG15","DAY16","MQSFLAG16","DAY17","MQSFLAG17","DAY18","MQSFLAG18","DAY19","MQSFLAG19","DAY20","MQSFLAG20","DAY21","MQSFLAG21","DAY22","MQSFLAG22","DAY23","MQSFLAG23","DAY24","MQSFLAG24","DAY25","MQSFLAG25","DAY26","MQSFLAG26","DAY27","MQSFLAG27","DAY28","MQSFLAG28","DAY29","MQSFLAG29","DAY30","MQSFLAG30","DAY31","MQSFLAG31"),sep = "\t",as.is=TRUE)   

    if (dim(ghcnd.data)[1]==0) {
      print("Warning : No data found for that station")
      ## setwd(oldpath)
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

##Dailyformat <- function(stid = NULL, param = "t2m", src = "GHCND", ver = "v3.2.0.20130120", adj = TRUE, path = "/klimadata/work/abdelkaderm/data/", force = TRUE, flags = FALSE, test = TRUE , verbose=TRUE) {
##  
##  obs.meta <- metaghcnd(stid = stations[2])
##  invfile <- attr(obs.meta,"inv_filename")
##  version.por <- attr(obs.meta,"version_upd")
##  version.upd <- attr(obs.meta,"version_por") 
##  
##  obs.data <- dataghcnd(stid = stations[2])
##  obs.data = data.frame(subset(obs.data,select=c(-1:-3))) # remove duplicate meta info
##  
##  obsfinalformat <- list(PARAM = rownames(table(obs.data$ELEMENT)),ID = paste(obs.meta$COUNTRY.CODE,obs.meta$STN,sep=""),NAME = obs.meta$NAME, LAT =obs.meta$LATITUDE,LON =obs.meta$LONGITUDE,ELEV = obs.meta$STNELEV,EXTRA = obs.meta$EXTRA,VAL = obs,SOURCE = upper.case(source), OR_VERSION = version.por , UPD_VERSION = version.upd, ADJUSTED = upper.case(adj), INV_FILE = invfile, DATA_FILE = fdata, ftp_src_link = "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn")}
