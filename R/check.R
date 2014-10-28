                                        # Name		: check.ncdf4
                                        # Description	: check the netcdf file attributes and data and update attributes if necessary for further use within ESD package.
                                        # Author 	: A. Mezghani, METNO
                                        # contact 	: abdelkaderm@met.no
                                        # Last Update	: 11-04-2013 ; 17-03-2014
                                        # require	: ncdf4

                                        # Define check as a method
                                        # check <- function(ncid,...) UseMethod("check")

## source("esd/R/frequency.R")

check.ncdf4 <- function(ncid, param="auto",verbose = FALSE) { ## use.cdfcont = FALSE - AM 22-10-2013 not used any more ! 
    browser()
    ## Load library
    library(ncdf4)
    
    ## Checking : Number of variables and select only one from the netcdf file, get variable attributes in v1. The user should decide between the list of variables
    if (tolower(param) == "auto") {
        if (ncid$nvars > 1) {
            i <- grep(param, names(ncid$var))
            if (length(i) == 0) i <- as.integer(readline(paste("Choose variable ",paste(namevars,collapse="/") ,"(from 1 - ",length(namevars), "): ",sep = "")))
            if (!is.integer(i)) stop("You should introduce an integer value and at least select one variable") 
        } else i <- 1
        param <- names(ncid$var)[i] # ; rm(i)
        v1 <- ncid$var[[i]] 
    } else {
        v1 <- NULL
        v1 <- eval(parse(text=paste("ncid$var$",param,sep="")))
        if (is.null(v1)) stop(paste("Variable ",param," could not be found !",sep=""))
    }
    ## Checking : Variable dimensions ...
    ndims <- eval(parse(text=paste("ncid$var$",param,"$ndims",sep="")))
    dimnames <- rep(NA,ndims)
    if (ndims>0) {
        for (j in 1:ndims) dimnames[j] <- eval(parse(text=paste("ncid$var$",param,"$dim[[",j,"]]$name",sep="")))
        if (verbose) print("Checking Dimensions --> [ok]")
        if (verbose) print(paste(as.character(ndims), " dimension(s) has(have) been found :"))
        if (verbose) print(dimnames)
    } else {
        stop("Checking Dimensions --> [fail]")
        if (verbose) print("The variable has no dimensions. The file may be corrupted!")  
    }
    dimnames <- tolower(dimnames)
    ## Get all attributes in model, check and update
    model <- ncatt_get(ncid,0)
    ## Update CMIP3 attributes to match those of CMIP5 
    mnames <- names(model)
    history <- ncatt_get(ncid,0,"history")
    ## browser()
    if (ncatt_get(ncid,0,"project_id")$hasatt) {
        if (model$project_id=="IPCC Fourth Assessment") {
            model$project_id <- "CMIP3"
            if (verbose) print("project_id IPCC Fourth Assessment set to CMIP3")
        }
    } else if (length(grep("sres",tolower(history$value)))>0) model$project_id <- "CMIP3"
    else if (length(grep("rcp",tolower(history$value)))>0) model$project_id <- "CMIP5"
    else {if (verbose) print("project_id is missing from file attributes")}
    
    if (!ncatt_get(ncid,0,"model_id")$hasatt & ncatt_get(ncid,0,"project_id")$hasatt) {   
        hist2 <- unlist(strsplit(history$value,split=c(" ")))
        ih <- grep("tas",hist2)
        if (model$project_id=="CMIP3") {
            txt <- hist2[ih[1]] 
            model$model_id <- cmip3.model_id(txt)
        }
        else if (model$project_id=="CMIP5") {
            txt <- hist2[ih[2]]
            model$model_id <- cmip5.model_id(txt)
        }
    }
    ## print(ncatt_get(ncid,0,"project_id")$hasatt)
    if (ncatt_get(ncid,0,"experiment_id")$hasatt & ncatt_get(ncid,0,"project_id")$hasatt) {
        if (tolower(model$project_id)=="cmip3") {
            txt <- unlist(strsplit(tolower(ncatt_get(ncid,0,"history")$value),split="/"))
            model$experiment_id <- paste(txt[grep("20c3m",txt)],txt[grep("sres",txt)],sep="-")
        }
    }
    ## browser()
    if (ncatt_get(ncid,0,"title")$hasatt) {
        title <- ncatt_get(ncid,0,"title")$value
        modelid <- unlist(strsplit(title,split=c(" ")))
        model$project_id <-modelid[1]
        model$experiment_id <-modelid[3]
        model$type <- modelid[2] 
    }
    
    ## END CMIP3 MODEL NAMES UPDATE
    ## browser()
    ## Checking : Time unit and origin
    ## Get system info
    a <- Sys.info()
    ## Get time dimension / val + attributes
    itime <- grep("tim", dimnames)
    if (length(itime) == 0)
        itime <- NULL
    else if (length(itime) > 1) stop("Error in time dim")
    else if (length(itime)==1) time <- eval(parse(text=paste("v1$dim[[",as.character(itime),"]]",sep="")))

    ## Get time unit and origin
    tatt <- tolower(names(time))
    itunit <- grep(c("unit"),tatt)
    itorigin <- grep(c("orig"),tatt)
    if (length(itunit)>0) {   
        tunit <- eval(parse(text = paste("time$",tatt[itunit],sep="")))
        if (verbose) print(paste("Time unit has been found in time$unit attribute (",tunit,")",sep=""))
    } else tunit <- NULL
    if (!is.null(tunit)) {
        if (verbose) print("Checking Time Unit --> [ok]")}
    else if (verbose) print("Checking Time Unit --> [fail]")
    if (!is.null(tunit) & (!is.null(grep("since",tunit)))) {
        if (verbose) print("Time unit and origin detected in time$unit attribute")
        tunit <- time$units
        torigin <- time$origin <- paste(unlist(strsplit(tunit,split=" "))[3:4],collapse=" ")
        tunit <- time$units <- unlist(strsplit(tunit,split=" "))[1]
        if (verbose) print(paste("Updating time$unit (",time$unit,") and creating time$origin (",time$origin,") attribute",sep= ""))
    } else if (length(itorigin)>0) {   
        torigin <- eval(parse(text = paste("time$",tatt[itorigin],sep="")))
        if (verbose) print(paste("Time origin has been found in time origin attribute and set to:",torigin,sep=" "))
    } else torigin <- NULL
    ##if (is.null(torigin) & is.null(tunit)) {
    ##   if (verbose) print("Attributes time unit and origin have not been found -> Sys.time() will be used !")
    ##   if ((tolower(a[1]) == "linux") & (use.cdfcont)) {
    ##      if (verbose) print("Linux & use.cdfcont : call cdfcont()")
    ##      time$origin <- torigin <- cdfcont(ncid$filename)$time.origin
    ##      time$units <- tunit <- cdfcont(ncid$filename)$time.unit
    ##   }
                                        #}
    if (is.null(torigin)) {
        if (verbose) print(paste("Time units:", tunit, " l=", min(tim[is.finite(time$vals)]),"-", max(tim[is.finite(time$vals)])))
        if (verbose) warning("Cannot determine the time origin!")
        if (verbose) warning("Example format: '15-Dec-1949'")
        if (verbose) warning("NCEP reanalysis typically: 01-01-01")
        if (verbose) warning("ERA-40 typically: 1900-01-01")
        torigin <- readline("Please enter a valid time origin: ")
    }
    
    if (!is.null(torigin)) {
        if (torigin == "1-01-01 00:00:00") {
            if (verbose) print("bad time origin")
            torigin <- "0001-01-01 00:00:00"
            if (verbose) print(paste("Re-setting time origin (",torigin,")",sep=""))
        }
    } else {
        torigin <- readline("Give me the time origin (format='YYYY-MM-DD' as '1949-22-01'):")
        if (!is.null(torigin)) {
            if (verbose) print(paste("Time origin set to =", torigin))
            else stop("Could not determine the time origin. The processing has been stopped !")
        }
    }

    if (!is.null(torigin)) {
        yorigin <- format.Date(as.Date(torigin),format="%Y")
        morigin <- format.Date(as.Date(torigin),format="%m")
        dorigin <- format.Date(as.Date(torigin),format="%d")
        if (as.numeric(yorigin) == 0) {
            if (verbose) warning("There is no year zero (Press et al., Numerical recipies)")
            yorigin <- 0
            if (verbose) print(paste("Warning : Year origin has been set to:",as.character(1900),sep="->"))     
        }
        if (is.na(dorigin)) {
            if (verbose) warning("Warning : Day origin is missing !")
            dorigin <- 1
            if (verbose) warning("Warning : Day origin has been set to:",dorigin)
        }  
        if (is.na(dorigin)) {
            if (verbose) warning("Warning : Month origin is missing !")
            morigin <- 1
            if (verbose) warning("Warning : Month origin has been set to:",morigin)
        }
        torigin1 <- paste(yorigin,morigin,dorigin,sep="-")
        torigin <- paste(torigin1,unlist(strsplit(torigin,split=" "))[2],sep=" ") 
    }

    if (!is.null(torigin)) {if (verbose) print("Checking Time Origin --> [ok]")} else if (verbose) print("Checking Time Origin --> [fail]")
    ##browser()
    ## Checking : Frequency
    type <- c("year","season","months","Days","Hours","minutes","seconds")
    type.abb <- substr(tolower(type),1,3)
    ## Initialize
    freq.att <- NULL
    ifreq <- grep("freq",names(model))
    if (length(ifreq)>0) {  
        itype <- grep(tolower(eval(parse(text=paste("model$",names(model)[ifreq],sep="")))),tolower(type))
        if (length(itype>0)) {
            if (verbose) print(paste("Frequency has been found in model$frequency attribute (",type[itype],")",sep="")) 
            freq.att <- frequency.name[grep(model$frequency,frequency.abb)]
        }
        if (verbose) print("Checking Frequency from attribute --> [ok]")
    } else {
        if (verbose) print("Checking Frequency from attribute --> [fail]")
        if (verbose) print("Frequency has not been found in the attributes") 
    }
    ## browser()
    ## Checking frequency from data
    frequency <- freq.data <- NULL
    freq.data <- frequency.data(data=as.vector(time$vals),unit=tunit,verbose=FALSE)
    if (!is.null(freq.data)) {if (verbose) print("Checking Frequency from the data --> [ok]")} else if (verbose) print("Checking Frequency from the data --> [fail]")
    ## Checking Calendar attribute if any, otherwise set to "ordinary"  # Possible values for CMIP5 files are : "365_day" , "standard" , "proleptic_gregorian" , "360_day"
    ical <- grep(c("calend"),tatt)
    if (length(ical)>0) {   
        calendar.att <- eval(parse(text = paste("time$",tatt[ical],sep="")))
        if (verbose) print("Checking Calendar from time attribute --> [ok]") 
        if (verbose) print(paste("Calendar attribute has been found in time$calendar (",time$calendar,")",sep =""))
    } else {
        if (verbose) print("Checking Calendar from time attribute --> [fail]")
        calendar.att <- NULL
        warnings("Calendar attribute has not been found in the meta data")
    }
    ## Identifying starting and ending dates for the data if possible
    if (!is.null(torigin)) {
        yorigin <- as.numeric(format.Date(torigin,format="%Y"))
        morigin <- as.numeric(format.Date(torigin,format="%m"))
        dorigin <- as.numeric(format.Date(torigin,format="%d"))
    }
    ## Get calendar from attribute if any and create vector of dates vdate
    browser()
    if (!is.null(calendar.att)) {
        if (grepl("gregorian",calendar.att) | grepl("standard",calendar.att)) {
            if (grepl("min",tunit)) time$vdate <- as.Date((time$vals/(24*60)),origin=as.Date(torigin))
            if (grepl("hou",tunit)) time$vdate <- as.Date((time$vals/24),origin=as.Date(torigin))
            if (grepl("day",tunit)) {
                time$vdate <- as.Date((time$vals),origin=as.Date(torigin))   
            }
            if (grepl("mon",tunit)) {
                if (sum(round(diff(time$vals)) > 1) < 1) {
                    year1 <- time$vals[1]%/%12 + yorigin
                    month1 <- morigin
                    time$vdate <- seq(as.Date(paste(as.character(year1),month1,"01",sep="-")), by = "month",length.out=length(time$vals))
                } else warnings("Monthly data are Mangeled") 
            } 
        } else if (!is.na(strtoi(substr(calendar.att, 1, 3)))) {
            if (verbose) print(paste(substr(calendar.att,1, 3), "-days' model year found in calendar attribute"))
            time$daysayear <- as.numeric(substr(calendar.att, 1, 3))
            if (!is.null(time$daysayear)) if (verbose) print(paste("Creating time$daysayear attribute and setting attribute to ", time$daysayear, sep=" "))
        }
        if (!is.null(time$daysayear)) {
            if (time$daysayear==365) 
                mndays <- c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in each month
            else if (time$daysayear==360)
                mndays <- rep(30,12) # Number of days in each month
            if (!is.null(time$daysayear) & !is.null(mndays)) {
                year1 <- time$vals[1]%/%time$daysayear + yorigin
                month1 <- morigin
                
                if (sum(diff(time$vals%/%time$daysayear) > 1) & (verbose)) warnings("Jumps of years has been found in the time series ")
                if (time$vals[1]%%time$daysayear > 27) {
                    year1 <- year1 + 1
                    month1 <- month1 + 1
                } 
                if (month1>12) month1 <- month1 - 12 
                                        # construct vdate
                months <- ((time$vals%%time$daysayear)%/%round(mean(mndays))) + 1
                years <- time$vals%/%time$daysayear + yorigin
                                        #shifting mndays by month1 to start with different initial months than january (1)
                mndays <- c(0,mndays[month1:length(mndays)-1],mndays[1:month1-1])
                days <- time$vals%%time$daysayear - rep(cumsum(mndays),time$len/12)
                if ((sum(diff(months) > 1) > 1) | (sum(diff(years) > 1) > 1) | (sum(round(abs(diff(days)))>2)) > 1) {
                    warnings("Warning: Jumps in data have been found !")
                    warnings("Warning: Trust the first date and force a continuous vector of dates !")
                    time$vdate <- seq(as.Date(paste(as.character(year1),month1,"01",sep="-")), by = "month",length.out=time$len)
                } else time$vdate <- as.Date(paste(years,months,"01",sep="-")) #round (days)                  
            }  
        } 
        if (verbose) print(paste("Starting date : ",time$vdate[1],"Ending date : ",time$vdate[length(time$vdate)], sep = " "))
    } else {
        if (verbose) print("warnings : Automatic detection of the calendar")
        calendar.detect <- "auto"
        ## browser()                                    # NOT COMPLETE ...
        if (grepl("sec",tunit)) time$vdate <- as.Date((time$vals/(24*3600)),origin=as.Date(torigin))
        if (grepl("hou",tunit)) time$vdate <- as.Date((time$vals/24),origin=as.Date(torigin))
        if (grepl("day",tunit)) time$vdate <- as.Date((time$vals),origin=as.Date(torigin))   
        if (grepl("mon",tunit)) {
            if (sum(diff(time$vals>1)) < 1) {
                year1 <- time$vals[1]%/%12 + yorigin
                month1 <- morigin
                time$vdate <- seq(as.Date(paste(as.character(year1),month1,"15",sep="-")), by = "month",length.out=length(time$vals))
            } else warnings("Monthly data are Mangeled") 
        } 
    }
    if ((length(time$vdate)>0) & (sum(diff(as.numeric(format.Date(time$vdate,"%m")))>1)) & (verbose)) stop("Vector date is mangeled ! Need extra check !")
    ## Checking the data / Extra checks / Automatic calendar detection / etc.
    ## Check 1 # Regular frequency
    ## browser()
    if (!is.null(time$vdate)) dt <- as.numeric(rownames(table(diff(time$vdate)))) else dt <- NULL
    if (!is.null(time$vdate)) {
        if (verbose) print("Vector of date is in the form :")
        if (verbose) print(str(time$vdate))
        if (verbose) print(diff(time$vdate))
    } else {
        if (grepl("min",tunit)) dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/(24*60)))))
        if (grepl("day",tunit)) dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
        if (grepl("hou",tunit)) dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/24))))
        if (grepl("mon",tunit)) dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
        if (length(dt)==1) {
            if (verbose) print("Regular frequency has been detected from the data")
        } else if (verbose) print("Irregular frequency has been detected from the data")
        if ((length(dt)==3) & grepl("day",tunit)) {
            if (verbose) print(paste("Calendar is likely to be a 365-",tunit," with: ",as.character(length(dt))," irregular frequencies",sep = ""))
            dt <- c(28,30,31)
            if (verbose) print(paste(as.character(dt),tunit,sep="-"))
        }
        if ((length(dt)==4) & grepl("day",tunit)) {
            if (verbose) print(paste("Calendar is likely to be a Gregorian 365/366-",tunit," with: ",as.character(length(dt))," irregular frequencies : ",sep=""))
            dt <- c(28,29,30,31)
            if (verbose) print(paste(as.character(dt),tunit,sep="-")) 
        }
        if (!is.null(time$daysayear)) {
            if ((length(dt)==3) & (time$daysayear != 365) & grepl("day",tunit))
                warning("Calendar does not match with detected frequencies")
        }
        if (length(ical)>0)
            if ((length(dt)!=4) & (grepl("gregorian",calendar.att) | grepl("standard",calendar.att)) & grepl("day",tunit))
                warning("Calendar does not match with detected frequencies")
        if ((length(dt)==2) | (length(dt)>4)) {
            if (verbose) print(paste("Warning : Irregular frequencies have been detected - The data might be corrupted and needs extra Checking !"))   
            if (verbose) print(paste(as.character(dt),tunit,sep=" "))
        }
    }
    ##browser()
    ## End check 1
    ## Begin check 2 if freq.att matches freq.data
    if (!is.null(freq.att)) {
        model$frequency <- freq.att 
        if (!is.null(freq.data)) {
            if (match(freq.att,freq.data)) {
                if (verbose) print("Frequency found in the attribute matches the frequency detected in data")
                model$frequency <- freq.data <- freq.att 
            } else warnings("Frequency found in the attribute does not match the frequency detected in data")
        } 
    } else if (!is.null(freq.data)) model$frequency <- freq.data
    else stop("Frequency could not be found, neither detected, the data might be corrupted !")
    
    if (!is.null(model$frequency)) {
        if (verbose) print(paste("Frequency set to ",model$frequency,sep=""))
        if (model$frequency=="month") {
            yr <-year(time$vdate)
            mo <- month(time$vdate)
            dy <- "01"
            time$vdate <- as.Date(paste(yr,mo,dy,sep="-"))     
        }
    }
    ## End check 2
    if (verbose) print("Checking --> [Done!]")
    ## use zoo library to format the data
    
    ## Extra Checking
    ##y.test <- data.e <- ncvar_get(ncid, v1$name, start = c(1,1,1), count = c(1,1,-1))
    ##ac.gcm <- data.frame(y = y.test, x1 = as.vector(cos(2 * pi * tim/daysayear)), x2 = as.vector(sin(2 * pi * tim/daysayear)))
    result <- list(model=model,time=time)
    invisible(result)
}

check.ncdf <- function(ncid, param="auto",verbose = FALSE) { ## use.cdfcont = FALSE - AM 22-10-2013 not used any more ! 
    
    ## Load library
    library(ncdf)
    
    ## Checking : Number of variables and select only one from the netcdf file, get variable attributes in v1. The user should decide between the list of variables
    if (tolower(param) == "auto") {
        if (ncid$nvars > 1) {
            i <- grep(param, names(ncid$var))
            if (length(i) == 0) i <- as.integer(readline(paste("Choose variable ",paste(namevars,collapse="/") ,"(from 1 - ",length(namevars), "): ",sep = "")))
            if (!is.integer(i)) stop("You should introduce an integer value and at least select one variable") 
        } else i <- 1
        param <- names(ncid$var)[i] # ; rm(i)
        v1 <- ncid$var[[i]] 
    } else {
        v1 <- NULL
        v1 <- eval(parse(text=paste("ncid$var$",param,sep="")))
        if (is.null(v1)) stop(paste("Variable ",param," could not be found !",sep=""))
    }
    ## Checking : Variable dimensions ...
    ndims <- eval(parse(text=paste("ncid$var$",param,"$ndims",sep="")))
    dimnames <- rep(NA,ndims)
    if (ndims>0) {
        for (j in 1:ndims) dimnames[j] <- eval(parse(text=paste("ncid$var$",param,"$dim[[",j,"]]$name",sep="")))
        if (verbose) print("Checking Dimensions --> [ok]")
        if (verbose) print(paste(as.character(ndims), " dimension(s) has(have) been found :"))
        if (verbose) print(dimnames)
    } else {
        stop("Checking Dimensions --> [fail]")
        if (verbose) print("The variable has no dimensions. The file may be corrupted!")  
    }
    dimnames <- tolower(dimnames)
    ## Get all attributes in model, check and update
    model <- att.get.ncdf(ncid,0,"global variable")
    ## Update CMIP3 attributes to match those of CMIP5 
    mnames <- names(model)
    history <- att.get.ncdf(ncid,0,"history")
    ## browser()
    if (att.get.ncdf(ncid,0,"project_id")$hasatt) {
        model$project_id <- att.get.ncdf(ncid,0,"project_id")$value
        if (model$project_id=="IPCC Fourth Assessment") {
            model$project_id <- "CMIP3"
            if (verbose) print("project_id IPCC Fourth Assessment set to CMIP3")
        }
    } else if (length(grep("sres",tolower(history$value)))>0) model$project_id <- "CMIP3"
    else if (length(grep("rcp",tolower(history$value)))>0) model$project_id <- "CMIP5"
    else {if (verbose) print("project_id is missing from file attributes")}
    
    if (!att.get.ncdf(ncid,0,"model_id")$hasatt & att.get.ncdf(ncid,0,"project_id")$hasatt) {   
        hist2 <- unlist(strsplit(history$value,split=c(" ")))
        ih <- grep("tas",hist2)
        if (model$project_id=="CMIP3") {
            txt <- hist2[ih[1]] 
            model$model_id <- cmip3.model_id(txt)
        }
        else if (model$project_id=="CMIP5") {
            txt <- hist2[ih[2]]
            model$model_id <- cmip5.model_id(txt)
        }
    }
    ## print(ncatt_get(ncid,0,"project_id")$hasatt)
    if (att.get.ncdf(ncid,0,"experiment_id")$hasatt & att.get.ncdf(ncid,0,"project_id")$hasatt) {
        if (tolower(model$project_id)=="cmip3") {
            txt <- unlist(strsplit(tolower(att.get.ncdf(ncid,0,"history")$value),split="/"))
            model$experiment_id <- paste(txt[grep("20c3m",txt)],txt[grep("sres",txt)],sep="-")
        }
    }
    ## browser()
    if (att.get.ncdf(ncid,0,"title")$hasatt) {
        title <- att.get.ncdf(ncid,0,"title")$value
        modelid <- unlist(strsplit(title,split=c(" ")))
        model$model_id <-modelid[1]
        model$experiment_id <-modelid[grep('rcp|sres',tolower(modelid))]
        model$project_id <-modelid[grep('cmip',tolower(modelid))]
        model$type <- modelid[2] 
    }
    
    ## END CMIP3 MODEL NAMES UPDATE
    ## browser()
    ## Checking : Time unit and origin
    ## Get system info
    a <- Sys.info()
    ## Get time dimension / val + attributes
    itime <- grep("tim", dimnames)
    if (length(itime) == 0)
        itime <- NULL
    else if (length(itime) > 1) stop("Error in time dim")
    else if (length(itime)==1) time <- eval(parse(text=paste("v1$dim[[",as.character(itime),"]]",sep="")))

    ## Get time unit and origin
    tatt <- tolower(names(time))
    itunit <- grep(c("unit"),tatt)
    itorigin <- grep(c("orig"),tatt)
    if (length(itunit)>0) {   
        tunit <- eval(parse(text = paste("time$",tatt[itunit],sep="")))
        if (verbose) print(paste("Time unit has been found in time$unit attribute (",tunit,")",sep=""))
    } else tunit <- NULL
    if (!is.null(tunit)) {
        if (verbose) print("Checking Time Unit --> [ok]")}
    else if (verbose) print("Checking Time Unit --> [fail]")
    if (!is.null(tunit) & (!is.null(grep("since",tunit)))) {
        if (verbose) print("Time unit and origin detected in time$unit attribute")
        tunit <- time$units
        torigin <- time$origin <- paste(unlist(strsplit(tunit,split=" "))[3:4],collapse=" ")
        tunit <- time$units <- unlist(strsplit(tunit,split=" "))[1]
        if (verbose) print(paste("Updating time$unit (",time$unit,") and creating time$origin (",time$origin,") attribute",sep= ""))
    } else if (length(itorigin)>0) {   
        torigin <- eval(parse(text = paste("time$",tatt[itorigin],sep="")))
        if (verbose) print(paste("Time origin has been found in time origin attribute and set to:",torigin,sep=" "))
    } else torigin <- NULL
    ##if (is.null(torigin) & is.null(tunit)) {
    ##   if (verbose) print("Attributes time unit and origin have not been found -> Sys.time() will be used !")
    ##   if ((tolower(a[1]) == "linux") & (use.cdfcont)) {
    ##      if (verbose) print("Linux & use.cdfcont : call cdfcont()")
    ##      time$origin <- torigin <- cdfcont(ncid$filename)$time.origin
    ##      time$units <- tunit <- cdfcont(ncid$filename)$time.unit
    ##   }
                                        #}
    if (is.null(torigin)) {
        if (verbose) print(paste("Time units:", tunit, " l=", min(tim[is.finite(time$vals)]),"-", max(tim[is.finite(time$vals)])))
        if (verbose) warning("Cannot determine the time origin!")
        if (verbose) warning("Example format: '15-Dec-1949'")
        if (verbose) warning("NCEP reanalysis typically: 01-01-01")
        if (verbose) warning("ERA-40 typically: 1900-01-01")
        torigin <- readline("Please enter a valid time origin: ")
    }
    
    if (!is.null(torigin)) {
        if (torigin == "1-01-01 00:00:00") {
            if (verbose) print("bad time origin")
            torigin <- "0001-01-01 00:00:00"
            if (verbose) print(paste("Re-setting time origin (",torigin,")",sep=""))
        }
    } else {
        torigin <- readline("Give me the time origin (format='YYYY-MM-DD' as '1949-22-01'):")
        if (!is.null(torigin)) {
            if (verbose) print(paste("Time origin set to =", torigin))
            else stop("Could not determine the time origin. The processing has been stopped !")
        }
    }

    if (!is.null(torigin)) {
        yorigin <- format.Date(as.Date(torigin),format="%Y")
        morigin <- format.Date(as.Date(torigin),format="%m")
        dorigin <- format.Date(as.Date(torigin),format="%d")
        if (as.numeric(yorigin) == 0) {
            if (verbose) warning("There is no year zero (Press et al., Numerical recipies)")
            yorigin <- 0
            if (verbose) print(paste("Warning : Year origin has been set to:",as.character(1900),sep="->"))     
        }
        if (is.na(dorigin)) {
            if (verbose) warning("Warning : Day origin is missing !")
            dorigin <- 1
            if (verbose) warning("Warning : Day origin has been set to:",dorigin)
        }  
        if (is.na(dorigin)) {
            if (verbose) warning("Warning : Month origin is missing !")
            morigin <- 1
            if (verbose) warning("Warning : Month origin has been set to:",morigin)
        }
        torigin1 <- paste(yorigin,morigin,dorigin,sep="-")
        torigin <- paste(torigin1,unlist(strsplit(torigin,split=" "))[2],sep=" ") 
    }

    if (!is.null(torigin)) {if (verbose) print("Checking Time Origin --> [ok]")} else if (verbose) print("Checking Time Origin --> [fail]")
    ##browser()
    ## Checking : Frequency
    type <- c("year","season","months","Days","Hours","minutes","seconds")
    type.abb <- substr(tolower(type),1,3)
    ## Initialize
    freq.att <- NULL
    ifreq <- grep("freq",names(model))
    if (length(ifreq)>0) {  
        itype <- grep(tolower(eval(parse(text=paste("model$",names(model)[ifreq],sep="")))),tolower(type))
        if (length(itype>0)) {
            if (verbose) print(paste("Frequency has been found in model$frequency attribute (",type[itype],")",sep="")) 
            freq.att <- frequency.name[grep(model$frequency,frequency.abb)]
        }
        if (verbose) print("Checking Frequency from attribute --> [ok]")
    } else {
        if (verbose) print("Checking Frequency from attribute --> [fail]")
        if (verbose) print("Frequency has not been found in the attributes") 
    }
    ## browser()
    ## Checking frequency from data
    frequency <- freq.data <- NULL
    freq.data <- frequency.data(data=as.vector(time$vals),unit=tunit,verbose=FALSE)
    if (!is.null(freq.data)) {if (verbose) print("Checking Frequency from the data --> [ok]")} else if (verbose) print("Checking Frequency from the data --> [fail]")
    ## Checking Calendar attribute if any, otherwise set to "ordinary"  # Possible values for CMIP5 files are : "365_day" , "standard" , "proleptic_gregorian" , "360_day"
    ical <- grep(c("calend"),tatt)
    if (length(ical)>0) {   
        calendar.att <- eval(parse(text = paste("time$",tatt[ical],sep="")))
        if (verbose) print("Checking Calendar from time attribute --> [ok]") 
        if (verbose) print(paste("Calendar attribute has been found in time$calendar (",time$calendar,")",sep =""))
    } else {
        if (verbose) print("Checking Calendar from time attribute --> [fail]")
        calendar.att <- NULL
        warnings("Calendar attribute has not been found in the meta data")
    }
    ## Identifying starting and ending dates for the data if possible
    if (!is.null(torigin)) {
        yorigin <- as.numeric(format.Date(torigin,format="%Y"))
        morigin <- as.numeric(format.Date(torigin,format="%m"))
        dorigin <- as.numeric(format.Date(torigin,format="%d"))
    }
    ## Get calendar from attribute if any and create vector of dates vdate
    ## browser()
    if (!is.null(calendar.att)) {
        if (grepl("gregorian",calendar.att) | grepl("standard",calendar.att)) {
            if (grepl("min",tunit)) time$vdate <- as.Date((time$vals/(24*60)),origin=as.Date(torigin))
            if (grepl("hou",tunit)) time$vdate <- as.Date((time$vals/24),origin=as.Date(torigin))
            if (grepl("day",tunit)) {
                time$vdate <- as.Date((time$vals),origin=as.Date(torigin))   
            }
            if (grepl("mon",tunit)) {
                if (sum(round(diff(time$vals)) > 1) < 1) {
                    year1 <- time$vals[1]%/%12 + yorigin
                    month1 <- morigin
                    time$vdate <- seq(as.Date(paste(as.character(year1),month1,"01",sep="-")), by = "month",length.out=length(time$vals))
                } else warnings("Monthly data are Mangeled") 
            } 
        } else if (!is.na(strtoi(substr(calendar.att, 1, 3)))) {
            if (verbose) print(paste(substr(calendar.att,1, 3), "-days' model year found in calendar attribute"))
            time$daysayear <- as.numeric(substr(calendar.att, 1, 3))
            if (!is.null(time$daysayear)) if (verbose) print(paste("Creating time$daysayear attribute and setting attribute to ", time$daysayear, sep=" "))
        }
        if (!is.null(time$daysayear)) {
            if (time$daysayear==365) 
                mndays <- c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in each month
            else if (time$daysayear==360)
                mndays <- rep(30,12) # Number of days in each month
            if (!is.null(time$daysayear) & !is.null(mndays)) {
                year1 <- time$vals[1]%/%time$daysayear + yorigin
                month1 <- morigin
                
                if (sum(diff(time$vals%/%time$daysayear) > 1) & (verbose)) warnings("Jumps of years has been found in the time series ")
                if (time$vals[1]%%time$daysayear > 27) {
                    year1 <- year1 + 1
                    month1 <- month1 + 1
                } 
                if (month1>12) month1 <- month1 - 12 
                                        # construct vdate
                months <- ((time$vals%%time$daysayear)%/%round(mean(mndays))) + 1
                years <- time$vals%/%time$daysayear + yorigin
                                        #shifting mndays by month1 to start with different initial months than january (1)
                mndays <- c(0,mndays[month1:length(mndays)-1],mndays[1:month1-1])
                days <- time$vals%%time$daysayear - rep(cumsum(mndays),time$len/12)
                if ((sum(diff(months) > 1) > 1) | (sum(diff(years) > 1) > 1) | (sum(round(abs(diff(days)))>2)) > 1) {
                    warnings("Warning: Jumps in data have been found !")
                    warnings("Warning: Trust the first date and force a continuous vector of dates !")
                    time$vdate <- seq(as.Date(paste(as.character(year1),month1,"01",sep="-")), by = "month",length.out=time$len)
                } else time$vdate <- as.Date(paste(years,months,"01",sep="-")) #round (days)                  
            }  
        } 
        if (verbose) print(paste("Starting date : ",time$vdate[1],"Ending date : ",time$vdate[length(time$vdate)], sep = " "))
    } else {
        if (verbose) print("warnings : Automatic detection of the calendar")
        calendar.detect <- "auto"
        ## browser()                                    # NOT COMPLETE ...
        if (grepl("hou",tunit)) time$vdate <- as.Date((time$vals/24),origin=as.Date(torigin))
        if (grepl("day",tunit)) time$vdate <- as.Date((time$vals),origin=as.Date(torigin))   
        if (grepl("mon",tunit)) {
            if (sum(diff(time$vals>1)) < 1) {
                year1 <- time$vals[1]%/%12 + yorigin
                month1 <- morigin
                time$vdate <- seq(as.Date(paste(as.character(year1),month1,"15",sep="-")), by = "month",length.out=length(time$vals))
            } else warnings("Monthly data are Mangeled") 
        } 
    }
    if ((length(time$vdate)>0) & (sum(diff(as.numeric(format.Date(time$vdate,"%m")))>1)) & (verbose)) stop("Vector date is mangeled ! Need extra check !")
    ## Checking the data / Extra checks / Automatic calendar detection / etc.
    ## Check 1 # Regular frequency
    ## browser()
    if (!is.null(time$vdate)) dt <- as.numeric(rownames(table(diff(time$vdate)))) else dt <- NULL
    if (!is.null(time$vdate)) {
        if (verbose) print("Vector of date is in the form :")
        if (verbose) print(str(time$vdate))
        if (verbose) print(diff(time$vdate))
    } else {
        if (grepl("min",tunit)) dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/(24*60)))))
        if (grepl("day",tunit)) dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
        if (grepl("hou",tunit)) dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/24))))
        if (grepl("mon",tunit)) dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
        if (length(dt)==1) {
            if (verbose) print("Regular frequency has been detected from the data")
        } else if (verbose) print("Irregular frequency has been detected from the data")
        if ((length(dt)==3) & grepl("day",tunit)) {
            if (verbose) print(paste("Calendar is likely to be a 365-",tunit," with: ",as.character(length(dt))," irregular frequencies",sep = ""))
            dt <- c(28,30,31)
            if (verbose) print(paste(as.character(dt),tunit,sep="-"))
        }
        if ((length(dt)==4) & grepl("day",tunit)) {
            if (verbose) print(paste("Calendar is likely to be a Gregorian 365/366-",tunit," with: ",as.character(length(dt))," irregular frequencies : ",sep=""))
            dt <- c(28,29,30,31)
            if (verbose) print(paste(as.character(dt),tunit,sep="-")) 
        }
        if (!is.null(time$daysayear)) {
            if ((length(dt)==3) & (time$daysayear != 365) & grepl("day",tunit))
                warning("Calendar does not match with detected frequencies")
        }
        if (length(ical)>0)
            if ((length(dt)!=4) & (grepl("gregorian",calendar.att) | grepl("standard",calendar.att)) & grepl("day",tunit))
                warning("Calendar does not match with detected frequencies")
        if ((length(dt)==2) | (length(dt)>4)) {
            if (verbose) print(paste("Warning : Irregular frequencies have been detected - The data might be corrupted and needs extra Checking !"))   
            if (verbose) print(paste(as.character(dt),tunit,sep=" "))
        }
    }
    ##browser()
    ## End check 1
    ## Begin check 2 if freq.att matches freq.data
    if (!is.null(freq.att)) {
        model$frequency <- freq.att 
        if (!is.null(freq.data)) {
            if (match(freq.att,freq.data)) {
                if (verbose) print("Frequency found in the attribute matches the frequency detected in data")
                model$frequency <- freq.data <- freq.att 
            } else warnings("Frequency found in the attribute does not match the frequency detected in data")
        } 
    } else if (!is.null(freq.data)) model$frequency <- freq.data
    else stop("Frequency could not be found, neither detected, the data might be corrupted !")
    
    if (!is.null(model$frequency)) {
        if (verbose) print(paste("Frequency set to ",model$frequency,sep=""))
        if (model$frequency=="month") {
            yr <-year(time$vdate)
            mo <- month(time$vdate)
            dy <- "01"
            time$vdate <- as.Date(paste(yr,mo,dy,sep="-"))     
        }
    }
    ## End check 2
    if (verbose) print("Checking --> [Done!]")
    ## use zoo library to format the data
    
    ## Extra Checking
    ##y.test <- data.e <- ncvar_get(ncid, v1$name, start = c(1,1,1), count = c(1,1,-1))
    ##ac.gcm <- data.frame(y = y.test, x1 = as.vector(cos(2 * pi * tim/daysayear)), x2 = as.vector(sin(2 * pi * tim/daysayear)))
    result <- list(model=model,time=time)
    invisible(result)
}
