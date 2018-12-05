## Author : A. Mezghani
## used in check.ncdf4()

## Full time frequency names
frequency.name <- c("year","season","month","day","hours","minutes","seconds")

## Abbreviated time frequency names 
frequency.abb  <-  substr(tolower(frequency.name),1,3)

frequency.data <- function(data=NULL,unit=NULL,verbose=FALSE) {
  if (is.null(data) | is.null(unit)) stop("Both data and unit are mandatory inputs !")
  # Initialize
  freq <- NULL
  # compute interval
  # Automatic detection 
  if (!is.null(unit) & length(data)>1) {
    dt <- round(median(diff(data),na.rm=TRUE))
    # KMP 2018-10-23: for subdaily data with unit=day, e.g. dt=0.25 for 6-hourly data
    if(dt==0) dt <- median(diff(data),na.rm=TRUE)
    #print(dt)
    # KMP 2018-10-23: in case the unit is days or hours instead of day or hour
    unit <- substr(tolower(unit),1,3)
    if ((((dt==31) | (dt==30)) & grepl("day",unit)) | ((dt==1) & grepl("mon",unit)) | (((dt==744) | (dt==728) | (dt==720)) & grepl("hou",unit) | (((dt==31) | (dt==1440) | (dt==44640)) & grepl("min",unit))))
      freq <- "month"
    if (((dt==1) & grepl("day",unit)) | ((dt==24) & grepl("hou",unit)))
      freq <- "day" 
    if (dt==1 & grepl("hou",unit)) 
      freq <- "hour"
    if ((dt==7) & grepl("day",unit))
      freq <- "week"
    if ((dt==14) & grepl("day",unit))
      freq <- "2weeks"
    if ((dt>=360) & grepl("day",unit))
      freq <- "year"
    if (((dt==1) & grepl("hou",unit)) | ((dt==3600) & grepl("hou",unit)))
      freq <- "hour"
    if((dt>1) & grepl("hou",unit))
      freq <- paste(dt,"hour",sep="")
    #if ((dt==6) & grepl("hou",unit))
    #  freq <- "6hour"
    #if ((dt==12) & grepl("hou",unit))
    #  freq <- "12hour"
    # KMP 2018-10-23: for subdaily data with unit=day and e.g. dt=0.25 (6-hourly)
    if(dt<1 & grepl("day",unit))
      freq <- paste(round(dt*24),"hour",sep="")
    if (dt >= 672 & dt <= 744 & grepl("hou",unit)) 
      freq <- "month"
    if ((dt==12) & grepl("mon",unit))
      freq <- "year"
    #HBE 11/04/18 added check for daily and hourly seasonal data
    if (((dt==3) & grepl("mon",unit)) |  ((dt<93) & (dt>88) & grepl("day",unit)) |  ((dt<2209) & (dt>2159) & grepl("hou",unit))) 
      freq <- "season"
  } 
  if (is.null(freq)) {
    # User entry
    if (verbose) print("Frequency could not be set automatically !")
    #HBE 11/04/18 set length options integer equal length options 
    print(paste(as.character(seq(1,length(frequency.name),1)),frequency.name,sep=":"))
    ifreq <- as.integer(readline("Please select a frequency number from the list before continue and press Enter:"))   
    if (!is.na(ifreq)) {
      freq <- frequency.abb[ifreq]
    } else {
      stop("You must provide an integer from the list !")
    }
    if (is.null(freq)) stop("Process is stopped: User should provide the frequency before continuing !")
  }
  return(freq)
  #invisible(freq) # Why both return and invisible?
}
