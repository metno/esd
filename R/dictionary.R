test.ele2param <- function() {
  ele2param(ele=NULL,src=NULL)
  ele2param(ele="101",src=NULL)
  ele2param(ele="101",src="GHCND")
  ele2param(ele="101",src="GHCNM")
  ele2param(ele="101",src="NORDKLIM")
  ele2param(ele="601",src=NULL)
  ele2param(ele="601",src="GHCND")
  ele2param(ele="601",src="GHCNM")
}

#' Dictionary and conversion tools between esd element identifier and variables
#' names and specifications.
#' 
#' Converts between esd element/parameter identifier and variable names from
#' different data sources.
#' 
#' @aliases ele2param esd2ele param2ele metno.frost.ele
#'
#' @param param parameter identifier
#' @param ele element identifier
#' @param src a character string for the acronym of the data source
#'
#' @return A meta data matrix object with the glossary of the different
#' variables or element identifiers as originally defined by each data source
#'
#' @keywords parameter,element
#'
#' @examples
#' # Display the glossary of paramerters or element identifiers for 'GHCND' data source.
#' print(ele2param(ele=NULL,src='GHCND'))
#' # Display the glossary of parameters or element identifiers for all data sources. 
#' print(ele2param())
#' # Convert mean temperature parameter (param) to esd element (ele).
#' ele <- esd2ele(param='t2m')
#' print(ele)
#' 
#' 
#' @export ele2param
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
  x <- merge(x,metno.frost.ele(),all=TRUE)
  ## KMP 2021-09-21: subset is broken and doesn't do what it is supposed to.
  ## Replacing the following two lines. 
  #if (length(src)>0) x <- subset(x, is=toupper(x[,6])==toupper(src))
  #if (length(ele)>0) x <- subset(x, is=x[,1]==ele)
  is <- rep(TRUE, nrow(x))
  if (length(src)>0) is <- is & toupper(x[,6]) %in% toupper(src)
  if (length(ele)>0) is <- is & x[,1] %in% ele
  x <- x[is, ]
  ## if ((length(src)==0) & (length(ele)==0)) 
  ##df <- as.data.frame(x,stringsAsFactors=FALSE)
  if (length(x)==0) print(paste("Selected element does not exist in the",src,"database",sep=" "))
  invisible(x)
}

#' @export
esd2ele <- function(param = NULL) {
  if (!is.null(param)) {
    ele <- switch(tolower(param),
                  't2m' = "101",
                  'tg' = "101",'tmean'="101",'tas'="101",'mean'="101",
                  'rr' = "601",
                  'slp' = "401",
                  'pon' = "402",
                  'pom' = "401",
                  'pox' = "403",
                  'prn' = "402",
                  'prm' = "401",
                  'prx' = "403",
                  'pp'  = "401",
                  'cc' = "801",
                  't2' = "101",
                  'precip' = "601",
                  '101' = "t2m",
                  '401' = "slp",
                  '402' = "pon",
                  '403' = "pox",
                  '601' = "precip",
                  '801' = "cc",
                  'tmin'="121",'tn'="121",
                  'tmax'="111",'tx'="111",
                  '121' ="tmin",
                  '111' = "tmax",
                  '901'  = "sd",
                  'sd' = '901',
                  '901' = 'sd',
                  'dd' = '502',
                  '502' = 'dd',
                  'fg' = '501', # Wind speed
                  '501' = 'fg',
                  'fx' = '503',
                  '503' = 'fx',
                  '504' = 'dd06',
                  '505' = 'dd12',
                  '506' = 'dd18',
                  'ffm' = '501',
                  'ffx' = '503',
                  'dd06' = '504',
                  'dd12' = '505',
                  'dd18' = '506',
                  '201' = 'hu',
                  'hu' = '201',
                  '301' = 'ss',
                  'ss' = '301',
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
  } else {
    ele <- 'NA'
  }
  return(ele)
}

#' @export
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


# Selected elements from NORDKLIM database 
nordklim.ele <- function() {
  x <- rbind(c("101" , "Mean temperature"		, "0.1"   		, "degree*C" 	, "T"),
  	     c("111" , "Mean maximum temperature"	, "0.1"   		, "degree*C" 	, "Tx"), 
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
    
    if (!file.exists('metno_element.txt')) { ## AM :  Need to update the file name, the problem is that the meta file does not contain elements !!
        x <- rbind(c("601" , "Precipitation"		, "1"	, "mm"		, "RR"),
                   c("401" , "Sea level pressure"	, "1"	, "hPa"		, "POM"),
                   c("402" , "Sea level pressure"	, "1"	, "hPa"		, "PON"),
                   c("403" , "Sea level pressure"	, "1"	, "hPa"		, "POX"),
                   c("999" , "Snowfall"			, "1"	, "mm"		, "SNOW"),
                   c("901" , "Snow depth"		, "1"	, "cm"		, "SA"),
                   c("101" , "Mean temperature"		, "1"	, "degree*C"	, "TAM"),
                   c("111" , "Maximum temperature"	, "1"	, "degree*C"	, "TAX"),
                   c("121" , "Minimum temperature"	, "1"	, "degree*C"	, "TAN"),
                   c("501" , "Wind speed"		, "1"	, "m/s"		, "FFM"),
                   c("502" , "Wind direction"		, "1"	, "degrees"	, "DD"),
                   c("504" , "Wind direction at 06 UTC"	, "1"	, "degrees"	, "DD06"),
                   c("505" , "Wind direction at 12 UTC"	, "1"	, "degrees"	, "DD12"),
                   c("506" , "Wind direction at 18 UTC"	, "1"	, "degrees"	, "DD18"),
                   c("503" , "Wind Gust"		, "1"	, "m/s"		, "FGX"))
        y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "METNO",stringsAsFactors=FALSE)
        
    } else {
      x <- read.csv(file = 'metno_element.csv',sep = ';',header = TRUE,stringsAsFactors = FALSE)
      y <- data.frame(element= x$Elem_codes, longname = x$NAME , scale_factor = rep(1,dim(x)[1]) , unit = x$UNIT , param = x$Elem_codes , source = "METNO")
    }
  return(y)
}

#' @export
metno.frost.ele <- function() {
  # TODO: wind direction is here so we can get DD06, DD12 and DD18 and perform an average.. how?
  x <- rbind(
    c("601" , "Precipitation"		, "1"	, "mm"		, "sum(precipitation_amount *)"),
    c("401" , "Sea level pressure"	, "1"	, "hPa"		, "mean(surface_air_pressure *)"),
    c("402" , "Sea level pressure"	, "1"	, "hPa"		, "min(surface_air_pressure *)"),
    c("403" , "Sea level pressure"	, "1"	, "hPa"		, "max(surface_air_pressure *)"),
    c("901" , "Snow depth"		, "1"	, "cm"		, "surface_snow_thickness"),
    c("101" , "Mean temperature"	, "1"	, "degree*C"	, "mean(air_temperature *)"),
    c("111" , "Maximum temperature"	, "1"	, "degree*C"	, "max(air_temperature *)"),
    c("121" , "Minimum temperature"	, "1"	, "degree*C"	, "min(air_temperature *)"),
    c("501" , "Wind speed"		, "1"	, "m/s"		, "mean(wind_speed *)"),
    c("502" , "Wind direction"		, "1"	, "degrees"	, "wind_from_direction"),
    c("503" , "Wind Gust"		, "1"	, "m/s"		, "max(wind_speed_of_gust *)")
  )
  y <- data.frame(element=x[,1] , longname = x[,2] , scale_factor = x[,3] , unit = x[,4] , param = x[,5] , source = "METNO.FROST" , stringsAsFactors = FALSE)
  return(y)
}

## KMP 2018-11-08: This function doesn't work. Element.Name is not defined and it contains an absolute path. 
## Also it is not in the NAMESPACE and is not used in any other function.
#wmo.ele.code <- function(pattern='^wind.*.speed.*.10') {
#    x <- read.csv('/disk1/downloads/fnmoc.B2L-058-001-B.txt',sep='\t',stringsAsFactor=FALSE,skip=5,header=TRUE)
#    ele.sel <- subset(x,subset= grepl(pattern,Element.Name,ignore.case=TRUE))
#    ele <- as.numeric(paste(ele.sel$X,ele.sel$Y,sep=''))
#    if (dim(ele.sel)[1] >1) {
#        print(ele.sel)
#        print('please refine your selection')
#    }
#    invisible(ele.sel)
#}
