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
                                     'tg' = "101",
                                     'rr' = "601",
                                     'slp' = "401",
                                     'cloud' = "801",
                                     't2' = "101",
                                     'precip' = "601",
                                     `101` = "t2m",
                                     `401` = "slp",
                                     `601` = "precip",
                                     `801` = "801",
                                     'tmin'="121",
                                     'tmax'="111",
                                     '121' ="tmin",
                                     '111' = "tmax",
                                     '901'  = "sd",
                                     'sd' = '901')
 else ele <- NULL
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
x <- rbind(c("101" , "Mean temperature"		 	, "0.1"   		, "°C" 	, "T"),
	   c("111" , "Mean maximum temperature"	 	, "0.1"   		, "°C" 	, "Tx"), 
 	   c("112" , "Highest maximum temperature"  	, "0.1"   		, "°C" 	, "Th"),
	   c("113" , "Day of Th" 			, "1"     		, "date", "Thd"),
	   c("121" , "Mean minimum temperature" 	, "0.1"   		, "°C" 	, "Tn"),
	   c("122" , "Lowest minimum temperature" 	, "0.1"   		, "°C"	, "Tl"),
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
x <- rbind(c("101" , "Mean temperature"			, "0.1"   		, "°C" 	, "T"),
	   c("111" , "Mean maximum temperature"		, "0.1"   		, "°C" 	, "Tx"), 
 	   c("112" , "Highest maximum temperature"  	, "0.1"   		, "°C" 	, "Th"),
	   c("113" , "Day of Th" 			, "1"     		, "date", "Thd"),
	   c("121" , "Mean minimum temperature" 	, "0.1"   		, "°C" 	, "Tn"),
	   c("122" , "Lowest minimum temperature" 	, "0.1"   		, "°C"	, "Tl"),
	   c("123" , "Day of Tl"			, "1"	  		, "date", "Tld"),
	   c("401" , "Mean Pressure" 		 	, "0.1"   		, "hPa" , "P"),
	   c("601" , "Precipitation Sum"		, "0.1"   		, "mm"	, "R"),
	   c("602" , "Maximum 1-day precipitation"	, "0.1"	  		, "mm"	, "Rx"),
	   c("701" , "Number of days with snow cover (> 50% covered)" , "1"	, "days", "dsc"),
	   c("801" , "Mean cloud cover"		, "1" 			, "%"	, "N")) 
y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "NACD",stringsAsFactors=FALSE)
return(y)
}
# Selected elements from NARP database / Same elements as nordklim database
narp.ele <- function() {
x <- rbind(c("101" , "Mean temperature"		 	, "0.1"   		, "°C" 	, "T"),
	   c("111" , "Mean maximum temperature"	 	, "0.1"   		, "°C" 	, "Tx"), 
 	   c("112" , "Highest maximum temperature"   	, "0.1"   		, "°C" 	, "Th"),
	   c("113" , "Day of Th" 			, "1"     		, "date", "Thd"),
	   c("121" , "Mean minimum temperature" 	, "0.1"   		, "°C" 	, "Tn"),
	   c("122" , "Lowest minimum temperature" 	, "0.1"   		, "°C"	, "Tl"),
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
ecad.ele <- function() {
x <- rbind(c("601" , "Precipitation amount"		, "0.1"	  		, "mm"	, "RR"),
	   c("401" , "Sea level pressure" 		, "0.1"	  		, "hPa"	, "PP"),
	   c("901" , "Snow depth" 			, "1"	  		, "cm" 	, "SD"),
	   c("111" , "Maximum temperature" 	 	, "0.1"   		, "°C" 	, "TX"),
	   c("121" , "Minimum temperature" 		, "0.1"	  		, "°C" 	, "TN"),
	   c("101" , "Mean temperature" 		, "0.1"   		, "°C" 	, "TG"))
y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "ECAD",stringsAsFactors=FALSE)
return(y)
}
# Selected elements from GHCNM database
ghcnm.ele <- function() {
x <- rbind(c("101" , "monthly mean temperature"	 	, "0.01" 		, "°C" , "TAVG"), 
	   c("111" , "monthly maximum temperature"  	, "0.01" 		, "°C" , "TMAX"),
	   c("121" , "monthly minimum temperature"  	, "0.01" 		, "°C" , "TMIN"))
y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "GHCNM",stringsAsFactors=FALSE)
return(y)
}
ghcnd.ele <- function() {
# Selected elements from GHCND database
x <- rbind(c("601" , "Precipitation"		 	, "0.1"	  		, "mm"	, "PRCP"),
	   c("999" , "Snowfall"			 	, "1" 	  		, "mm" 	, "SNOW"),
	   c("901" , "Snow depth"			, "1" 	  		, "mm" 	, "SNWD"),
	   c("111" , "Maximum temperature" 	 	, "0.1"   		, "°C"	, "TMAX"),
	   c("121" , "Minimum temperature" 	 	, "0.1"	  		, "°C" 	, "TMIN"))
y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "GHCND",stringsAsFactors=FALSE)
return(y)
}
metno.ele <- function() { ## must be updated - AM 2014-02-21
  ## Selected elements from GHCND database
  x <- rbind(c("601" , "Precipitation"		 	, "1"	  		, "mm"	, "RR"),
             c("999" , "Snowfall"			, "1" 	  		, "mm" 	, "SNOW"),
             c("901" , "Snow depth"			, "1" 	  		, "mm" 	, "SNWD"),
             c("101" , "Mean temperature"	        , "1" 		        , "°C"  , "TAM"),
             c("111" , "Maximum temperature" 	 	, "1"   		, "°C"	, "TAX"),
             c("121" , "Minimum temperature" 	 	, "1"	  		, "°C" 	, "TAN"))
  y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "METNO",stringsAsFactors=FALSE)
  return(y)
}
