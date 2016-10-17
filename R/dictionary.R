## Description : elements and variables dictionary for different data sources including conversion tools from elements to variables, and vice versa. 
## author  : Abdelkader Mezghani
## created : 12-06-2013
## updated : 12-06-2013

test.ele2param <- function() {

# RUN the following lines
ele2param(ele=NULL,src=NULL)
ele2param(ele="101",src=NULL)
ele2param(ele="101",src="GHCND")
ele2param(ele="101",src="GHCNM")
ele2param(ele="101",src="NORDKLIM")
ele2param(ele="601",src=NULL)
ele2param(ele="601",src="GHCND")
ele2param(ele="601",src="GHCNM")
}

esd2ele <- function(param = NULL) {
  if (!is.null(param)) ele <- switch(tolower(param),
                                     't2m' = "101",
                                     'tg' = "101",'tmean'="101",'tas'="101",'mean'="101", #REB 2016-07-25: more flexibility
                                     'rr' = "601",
                                     'slp' = "401",
                                     'cc' = "801",
                                     't2' = "101",
                                     'precip' = "601",
                                     '101' = "t2m",
                                     '401' = "slp",
                                     '601' = "precip",
                                     '801' = "cc",
                                     'tmin'="121",'tn'="121",                             #REB 2016-07-25                 
                                     'tmax'="111",'tx'="111",                             #REB 2016-07-25
                                     '121' ="tmin",
                                     '111' = "tmax",
                                     '901'  = "sd",
                                     'sd' = '901',
                                     'dd' = '502',
                                     '502' = 'dd',
                                     'fg' = '501', # Wind speed
                                     '501' = 'fg',
                                     'fx' = '503',
                                     '503' = 'fx',
                                     '201' = 'hu',
                                     'hu' = '201',
                                     '301' = 'ss',
                                     'ss' = '301',
                                     'sd'='901',
                                     '999'='sf',
                                     '122'='tl',
                                     '123'='tld',
                                     'tl'='122',
                                     'tld'='123',
                                     'dsc'='701',
                                     '701'='dsc',
                                     'th'='112',
                                     '112'='th',
                                     '113'='thd',
                                     'thd'='113',
                                     '602'='rx',
                                     'rx'='602')
 else ele <- 'NA'
return(ele)
}

param2ele <- function(param = NULL , src = NULL , verbose = TRUE) {
  ## convert arguments
  src <- toupper(src)
  param <- as.character(param)
  ## get elements and variable metadata info
                                        #browser()
  if (length(param)==0 & length(src)==0) {
    if (verbose) print("No selected variable")
    ele <- ele2param()
    print(ele)
  } else if (length(param) == 0 & length(src) >0) {  
    x <- ele2param(ele=NULL,src=src)
    ele <- as.character(x[,1])
    if (verbose) print("No conversion between variables and elements")
    ele <- NA
  } else if (length(param)> 0 & length(src) >0) {
    x <- ele2param(ele=NULL,src=src)
    ele <- as.character(x[x[,5]==param,1])
    if (length(ele)==0) ele <- NA
  } else ele <- NA 
  return(ele)
}
  
ele2param <- function(ele = NULL , src = NULL) {
# convert arguments
  src <- toupper(src)
  ele <- as.character(ele)
  ## combine all sources
  x <- merge(nordklim.ele(),nacd.ele(),all=TRUE)
  x <- merge(x,narp.ele(),all=TRUE)
  x <- merge(x,ecad.ele(),all=TRUE)
  x <- merge(x,ghcnm.ele(),all=TRUE)
  x <- merge(x,ghcnd.ele(),all=TRUE)
  x <- merge(x,metno.ele(),all=TRUE)
 
  if (length(src)>0) x <- subset(x, x[,6] == src)
  if (length(ele)>0) x <- subset(x, x[,1] == ele)
  ## if ((length(src)==0) & (length(ele)==0)) 
  ##df <- as.data.frame(x,stringsAsFactors=FALSE)
  if (length(x)==0) print(paste("Selected element does not exist in the",src,"database",sep=" "))
  invisible(x)
}
# Selected elements from NORDKLIM database 
nordklim.ele <- function() {
x <- rbind(c("101" , "Mean temperature"		 	, "0.1"   		, "degree*C" 	, "T"),
	   c("111" , "Mean maximum temperature"	 	, "0.1"   		, "degree*C" 	, "Tx"), 
 	   c("112" , "Highest maximum temperature"  	, "0.1"   		, "degree*C" 	, "Th"),
	   c("113" , "Day of Th" 			, "1"     		, "date", "Thd"),
	   c("121" , "Mean minimum temperature" 	, "0.1"   		, "degree*C" 	, "Tn"),
	   c("122" , "Lowest minimum temperature" 	, "0.1"   		, "degree*C"	, "Tl"),
	   c("123" , "Day of Tl"			, "1"	  		, "date", "Tld"),
	   c("401" , "Mean Pressure" 		 	, "0.1"   		, "hPa" , "P"),
	   c("601" , "Precipitation Sum"		, "0.1"   		, "mm"	, "R"),
	   c("602" , "Maximum 1-day precipitation"	, "0.1"	  		, "mm"	, "Rx"),
	   c("701" , "Number of days with snow cover (> 50% covered)" , "1"	, "days", "dsc"),
	   c("801" , "Mean cloud cover"	, "1" 			, "%"	, "N"))
y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "NORDKLIM",stringsAsFactors=FALSE)
return(y)
}
# Selected elements from NACD database / Same elements as nordklim database
nacd.ele <- function() {
x <- rbind(c("101" , "Mean temperature"			, "0.1"   		, "degree*C" 	, "T"),
	   c("111" , "Mean maximum temperature"		, "0.1"   		, "degree*C" 	, "Tx"), 
 	   c("112" , "Highest maximum temperature"  	, "0.1"   		, "degree*C" 	, "Th"),
	   c("113" , "Day of Th" 			, "1"     		, "date", "Thd"),
	   c("121" , "Mean minimum temperature" 	, "0.1"   		, "degree*C" 	, "Tn"),
	   c("122" , "Lowest minimum temperature" 	, "0.1"   		, "degree*C"	, "Tl"),
	   c("123" , "Day of Tl"			, "1"	  		, "date", "Tld"),
	   c("401" , "Mean Pressure" 		 	, "0.1"   		, "hPa" , "slp"),
	   c("601" , "Precipitation Sum"		, "0.1"   		, "mm"	, "R"),
	   c("602" , "Maximum 1-day precipitation"	, "0.1"	  		, "mm"	, "Rx"),
	   c("701" , "Number of days with snow cover (> 50% covered)" , "1"	, "days", "dsc"),
	   c("801" , "Mean cloud cover"		, "1" 			, "%"	, "Cloud")) 
y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "NACD",stringsAsFactors=FALSE)
return(y)
}
# Selected elements from NARP database / Same elements as nordklim database
narp.ele <- function() {
x <- rbind(c("101" , "Mean temperature"		 	, "0.1"   		, "degree*C" 	, "T"),
	   c("111" , "Mean maximum temperature"	 	, "0.1"   		, "degree*C" 	, "Tx"), 
 	   c("112" , "Highest maximum temperature"   	, "0.1"   		, "degree*C" 	, "Th"),
	   c("113" , "Day of Th" 			, "1"     		, "date", "Thd"),
	   c("121" , "Mean minimum temperature" 	, "0.1"   		, "degree*C" 	, "Tn"),
	   c("122" , "Lowest minimum temperature" 	, "0.1"   		, "degree*C"	, "Tl"),
	   c("123" , "Day of Tl"			, "1"	  		, "date", "Tld"),
	   c("401" , "Mean Pressure" 		 	, "0.1"   		, "hPa" , "P"),
	   c("601" , "Precipitation Sum"		, "0.1"   		, "mm"	, "R"),
	   c("602" , "Maximum 1-day precipitation"	, "0.1"	  		, "mm"	, "Rx"),
	   c("701" , "Number of days with snow cover (> 50% covered)" , "1"	, "days", "dsc"),
	   c("801" , "Mean cloud cover"		 	, "1" 			, "%"	, "N")) 
y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "NARP",stringsAsFactors=FALSE)
return(y)
}
# Selected elements from ECAD database
ecad.ele <- function(file='~/ecad_element.csv') {
    ## browser()
    ## if (file.exists(file)) {
    ##     x <- read.csv('~/ecad_element.csv',sep=',',stringsAsFactor=FALSE)
    ##     strsp <- strsplit(x$X.2,split=' ')
    ##     scale <- as.numeric(apply(as.matrix(x$X.2),1, function(x) strsplit(x,split=' ')[[1]][1]))
    ##     scale[is.na(scale)] <- 1
    ##     unit <- gsub("[[:digit:]]|\\.", "", x$X.2)
    ##     unit <- gsub(" ","",unit)
    ##     ## clean unit and scale
    ##     ok <- !(unit=="")
    ##     y <- data.frame(element=x$X[ok] , longname = x$X.1[ok] ,
    ##                     scale_factor = scale[ok],
    ##                     unit = unit[ok] ,
    ##                     param = x$X[ok] , source = "ECAD",stringsAsFactors=FALSE)
    ## } else {
        x <- rbind(c("601" , "Precipitation amount"		, "0.1"	  		, "mm"	        , "RR"),
                   c("401" , "Sea level pressure" 		, "0.1"	  		, "hPa"	        , "PP"),
                   c("901" , "Snow depth" 			, "1"	  		, "cm" 	        , "SD"),
                   c("111" , "Maximum temperature" 	 	, "0.1"   		, "degree*C" 	, "TX"),
                   c("121" , "Minimum temperature" 		, "0.1"	  		, "degree*C" 	, "TN"),
                   c("101" , "Mean temperature" 		, "0.1"   		, "degree*C" 	, "TG"),
                   c("501" , "Wind speed" 		        , "0.1"   		, "m/s" 	, "FG"),
                   c("502" , "Wind direction"         		, "1"   		, "degrees" 	, "DD"),
                   c("503" , "Wind Gust" 		        , "0.1"   		, "m/s" 	, "FX"),
                   c("201" , "Relative humidity" 		, "1"   		, "%   " 	, "HU"),
                   c("801" , "Cloud cover"      		, "1"   		, "oktas" 	, "CC"),
                   c("301" , "Sunshine"         		, "0.1"   		, "Hours" 	, "SS"))
        y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "ECAD",stringsAsFactors=FALSE)
    ## }
    return(y)   
}
# Selected elements from GHCNM database
ghcnm.ele <- function() {
x <- rbind(c("101" , "monthly mean temperature"	 	, "0.01" 		, "degree*C" , "TAVG"), 
	   c("111" , "monthly maximum temperature"  	, "0.01" 		, "degree*C" , "TMAX"),
	   c("121" , "monthly minimum temperature"  	, "0.01" 		, "degree*C" , "TMIN"))
y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "GHCNM",stringsAsFactors=FALSE)
return(y)
}
ghcnd.ele <- function() {
# Selected elements from GHCND database
x <- rbind(c("601" , "Precipitation"		 	, "0.1"	  		, "mm"	, "PRCP"),
	   c("999" , "Snowfall"			 	, "1" 	  		, "mm" 	, "SNOW"),
	   c("901" , "Snow depth"			, "1" 	  		, "mm" 	, "SNWD"),
	   c("111" , "Maximum temperature" 	 	, "0.1"   		, "degree*C"	, "TMAX"),
	   c("121" , "Minimum temperature" 	 	, "0.1"	  		, "degree*C" 	, "TMIN"))
y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "GHCND",stringsAsFactors=FALSE)
return(y)
}
metno.ele <- function() { ## must be updated - AM 2014-02-21
  ## Selected elements from GHCND database
    
    if (!file.exists('metno_element.txt')) {
      x <- rbind(c("601" , "Precipitation"		 	, "1"	  		, "mm"	, "RR"),
                 c("999" , "Snowfall"			, "1" 	  		, "mm" 	, "SNOW"),
                 c("901" , "Snow depth"			, "1" 	  		, "mm" 	, "SNWD"),
                 c("101" , "Mean temperature"	        , "1" 		        , "degree*C"  , "TAM"),
                 c("111" , "Maximum temperature" 	 	, "1"   		, "degree*C"	, "TAX"),
                 c("121" , "Minimum temperature" 	 	, "1"	  		, "degree*C" 	, "TAN"))
      y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "METNO",stringsAsFactors=FALSE)
      
      } else {
        x <- read.csv('metno_element.txt')

        y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = 1 , unit = x[,5] , param = x[,1] , source = "METNO",stringsAsFactors=FALSE)
      }
  return(y)
}

wmo.ele.code <- function(pattern='^wind.*.speed.*.10') {
    x <- read.csv('/disk1/downloads/fnmoc.B2L-058-001-B.txt',sep='\t',stringsAsFactor=FALSE,skip=5,header=TRUE)
    ele.sel <- subset(x,subset= grepl(pattern,Element.Name,ignore.case=TRUE))
    ele <- as.numeric(paste(ele.sel$X,ele.sel$Y,sep=''))
    if (dim(ele.sel)[1] >1) {
        print(ele.sel)
        print('please refine your selection')
    }
    invisible(ele.sel)
}
