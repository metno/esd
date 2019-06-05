#' Conversion to esd objects.
#' 
#' \code{year}, \code{month}, \code{day}, \code{season} return the years, months, days, and seasons
#' associated with the data.
#' 
#' @aliases month day
#' @seealso season season.default
#'
#' @param x an object of, e.g., class 'station', 'field', or 'zoo', or a date
#'
#' @return a numeric for \code{year}, \code{month}, and \code{day}; A numeric or character for \code{season}
#'
#' @import zoo
#'
#' @keywords utilities
#'
#' @examples
#' data(bjornholt)
#' year(bjornholt)
#' month(bjornholt)
#' day(bjornholt)
#' season(bjornholt)
#' season(bjornholt, format="numeric")
#' 
#' @export
year <- function(x) {
  #str(x); print(class(x)); print(index(x))
  if (inherits(x,'integer')) x <- as.numeric(x)
  
  if ( (inherits(x,'numeric')) & (min(x,na.rm=TRUE) > 0) &
      (max(x,na.rm=TRUE) < 3000) )
    return(x)
  
  if (inherits(x,c('station','field','zoo'))) {
    y <- year(index(x))
    return(y)
  }
  if (inherits(x,'trajectory')) {
    y <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")$year + 1900
    return(y)
  }
  if (inherits(x,'events')) {
    y <- strptime(x$date,format="%Y%m%d")$year + 1900
    return(y)
  }
  if (inherits(x,c("POSIXt","PCICt"))) {
    y <- as.numeric(format(x, '%Y'))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==10)) {
    y <- year(as.Date(x))
    return(y)
  }
  if (class(index(x))=="numeric") {
    y <- trunc(index(x))
    return(y)
  }
  #print("here"); print(index(x))
  if (class(x)[1]=="Date")
    y <- as.numeric(format(x, '%Y')) else
  if (class(x)[1]=="yearmon") y <- trunc(as.numeric(x)) else
  if (class(x)[1]=="yearqtr") y <- trunc(as.numeric(x)) else
  if (class(x)[1]=="character") y <- trunc(as.numeric(x)) else
  if (class(x)[1]=="numeric") y <- trunc(x) else
  if (class(x)[1]=="season") {
    # If season, then the first month is really the December month of the previous year
    month <- round(12*(as.numeric(index(x)) - trunc(as.numeric(index(x)))) + 1)
    y[is.element(month,1)] <- y[is.element(month,1)] - 1
  } else {
    print(paste(class(x)[1],' confused...'))
  }
  return(y)
}

#' @export
month <- function(x) {
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) & (max(x,na.rm=TRUE) < 13) )
    return(x) 
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) ) y <- rep(1,length(x))
  if (inherits(x,c('station','field','zoo'))) {
    #
    y <- month(index(x))
    return(y)
  }
  if (inherits(x,'trajectory')) {
    y <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")$mon + 1
    return(y)
  }
  if (inherits(x,'events')) {
    y <- strptime(x$date,format="%Y%m%d")$mon + 1
    return(y)
  }
  if (inherits(x,c("POSIXt","PCICt"))) {
    y <- as.numeric(format(x, '%m'))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==10)) {
    y <- month(as.Date(x))
    return(y)
  }
  #print(class(index(x)))
  if (class(x)[1]=="Date")
    y <- as.numeric(format(x, '%m')) else
  if (class(x)[1]=="yearmon")
    y <- round(12*(as.numeric(x) - trunc(as.numeric(x))) + 1) else
  if (class(x)[1]=="yearqtr") y <- round(12*(as.numeric(x) - trunc(as.numeric(x))) + 1)
  if (class(x)[1]=="season") {
    # If season, then the first month is really the months are DJF, MAM, JJA, and OND:
    y <- y - 1
    y[y==-1] <- 12
  }
  return(y)
}

#' @export
day <- function(x) {
  if (inherits(x,c('station','field','zoo'))) {
    y <- day(index(x))
    return(y)
  }
  if (inherits(x,'trajectory')) {
    y <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")$mday
    return(y)
  }
  if (inherits(x,'events')) {
    y <- strptime(x$date,format="%Y%m%d")$mday
    return(y)
  }
  if (inherits(x,c('numeric','integer'))) x <- as.numeric(x)
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) & (max(x,na.rm=TRUE) < 32) )
    return(x)
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) ) y <- rep(1,length(x))  
  if (inherits(x,c("POSIXt","PCICt"))) {
    y <- as.numeric(format(x, '%d'))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==10) ) {
    y <- day(as.Date(x))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==4) ) {
    y <- rep(1,length(x)) 
    return(y)
  }
  if (class(x)[1]=="Date") y <- as.numeric(format(x, '%d'))
  if (class(x)[1] %in% c("yearmon","yearqtr","season")) {
    y <- rep(1,length(x))
  }
  return(y)
}

#' @export
season <- function(x,format="character",verbose=FALSE) UseMethod("season")

#' Conversion to esd objects.
#'
#' Used to estimate Dec-Feb, Mar-May, Jun-Aug, and Sep-Nov statistics
#' Manipulate the zoo-object by shifting the year/chonology so that
#' zoo thinks the year defined as December-November is January-December.
#' 
#' \code{season} return the seasons associated with the data.
#' 
#' @aliases season
#' @seealso year month 
#'
#' @param x an object of, e.g., class 'station', 'field', or 'zoo', or a date
#' @param format for season, set the format of the output 'character' or 'numeric'
#'
#' @return a numeric or character
#'
#' @keywords utilities
#'
#' @examples
#' data(bjornholt)
#' year(bjornholt)
#' month(bjornholt)
#' day(bjornholt)
#' season(bjornholt)
#' season(bjornholt, format="numeric")
#' 
#' @export
season.default <- function(x,format="character",verbose=FALSE) {
  if(verbose) print("season.default")
  nt <- length(index(x))
  season <- rep('',nt)
  m <- month(x)
  if ( (inherits(x,'zoo')) & (format=="character") ) {
    for (i in 1:nt)  season[i] <- switch(m[i],
                                        '1'='djf','2'='djf','12'='djf',
                                         '3'='mam','4'='mam','5'='mam',
                                         '6'='jja','7'='jja','8'='jja',
                                         '9'='son','10'='son','11'='son')
  } else if ( (inherits(x,'zoo')) & (format=="numeric") ){
    for (i in 1:nt)  season[i] <- switch(m[i],'1'=1,'2'=1,'12'=1,
                                         '3'=2,'4'=2,'5'=2,
                                         '6'=3,'7'=3,'8'=3,
                                         '9'=4,'10'=4,'11'=4)
    season <- as.numeric(season)
  } else {
    season <- paste(substr(month.abb[as.numeric(rownames(table(month(x))))],1,1),sep='')
  }
  return(season)
}

# do not export
seasonal.yearmon <- function(x) {
  attr(season,'history') <- history.stamp(x)
  season
}