## Author : A. Mezghani
## used in check.ncdf4()

#' Full time frequency names
#' @export
frequency.name <- c("year","season","month","day","hours","minutes","seconds")

#' Abbreviated time frequency names
#' @export
frequency.abb  <-  substr(tolower(frequency.name),1,3)

#' Calculate the frequency
#'
#' @param data input object containing a time index
#' @param unit unit of time index in 'data'
#'
#' @export
datafrequency <- function(data=NULL,unit=NULL,verbose=FALSE) {
  if (is.null(data) | is.null(unit)) stop("Both data and unit are mandatory inputs !")
  # Initialize
  freq <- NULL
  # compute interval
  # Automatic detection 
  if (!is.null(unit) & length(data)>1) {
    dt <- round(median(diff(data),na.rm=TRUE))
    if(dt==0) dt <- median(diff(data),na.rm=TRUE)
    unit <- substr(tolower(unit),1,3)
    if ( ((dt>=360) & grepl("day",unit)) | (dt==12) & grepl("mon",unit)) {
      freq <- "year"
    } else if (((dt==3) & grepl("mon",unit)) |  ((dt<93) & (dt>88) & grepl("day",unit)) |
               ((dt<2209) & (dt>2159) & grepl("hou",unit))) {
      freq <- "season"
    } else if ((((dt==31) | (dt==30)) & grepl("day",unit)) |
        ((dt==1) & grepl("mon",unit)) |
        ((dt>=672 & dt<=744) & grepl("hou",unit) |
	(((dt==31) | (dt==1440) | (dt==44640)) & grepl("min",unit)))) {
      freq <- "month"
    } else if ((dt==14) & grepl("day",unit)) {
      freq <- "2weeks"
    } else if ((dt==7) & grepl("day",unit)) {
      freq <- "week"
    } else if (((dt==1) & grepl("day",unit)) | ((dt==24) & grepl("hou",unit))) {
      freq <- "day"
    } else if (((dt==1) & grepl("hou",unit)) | ((dt==3600) & grepl("hou",unit))) {
      freq <- "hour"
    } else if((dt>1) & grepl("hou",unit)) {
      freq <- paste(dt,"hour",sep="")
    } else if(dt<1 & grepl("day",unit)) {
      freq <- paste(round(dt*24),"hour",sep="")
    }
  } 
  if (is.null(freq)) {
    # User entry
    if (verbose) print("Frequency could not be set automatically !")
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
}
