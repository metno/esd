#' Count statistics
#' 
#' Returns a count of cases with values over s threshold value.
#' \code{exceedance} estimates statistics (default is mean) for peak-over-threshold, e.g. mean-over-threshold. 
#' Similar functions nclude: 
#' \code{nevents} returns the number of events with exceeding values (e.g. the
#' number of rainy days X > 1 mm/day).
#' \code{wetfreq} returns the wet-day frequency (a fraction) and \code{wetmean} returns the wet-day mean.
#' \code{qqgeom} produces a quantile-quantile plot of streak statistics comparing the
#' empirical quantiles with the distribution function quantiles (see \code{\link[stats]{qgeom}}). 
#'
#' @aliases count wetfreq wetmean nevents exceedance exceedance.default exceedance.station exceedance.field
#' @seealso spell hotsummerdays coldwinterdays nwetdays plot
#'
#' @importFrom stats glm qqline
#'
#' @param x station or field object
#' @param threshold threshold value
#' @param upper upper limit for maximum length - ignore any above this because they are likely erronoeus
#' @param fraction TRUE: divide the number of counts by number of samples
#' @param FUN function
#'
#' @return Station or field objects
#'
#' @keywords utilities
#' @examples
#' 
#' number of warm days:
#' data(ferder)
#' plot(as.seasons(ferder,FUN='count',threshold=20), new=FALSE)
#' 
#' Mild winter days - number of days in the winter season with
#' above freezing temperatures
#' try(coldwinterdays(ferder))
#'
#' Quantile-quantile plot
#' qqgeom(ferder, treshold=1, pois=TRUE)
#'
#'@export
count <- function(x,threshold=1,fraction=FALSE,...) {
  count <- sum(x > threshold,na.rm=TRUE)
  if (fraction) count <- count/sum(is.finite(x))
  return(count)
}

#' @export
wetfreq <- function(x,threshold=1,...) {
  if (inherits(x,'zoo')) x <- coredata(x)
  x <- x[is.finite(x)]
  x[x < threshold] <- NA
  y <- sum(is.finite(x))/length(x)
  return(y)
}

#' @export
nevents <- function(x,threshold=1,...) {
  y <- exceedance.default(x,threshold=threshold,FUN="count")
  return(y)
}

#' @export
wetmean <- function(x,threshold=1,...) {
  ## REB 2015-03-23 - slow
  #   y <- exceedance.default(x,threshold=threshold,FUN="mean")
  ## REB 2015-03-23 - faster
  ## Also add the standard error estimate based on the sample size
  ## and assuming an exponential distribtion for daily data
  ## (sigma = mu)
  if (inherits(x,'zoo')) x <- coredata(x)
  x[x < threshold] <- NA
  y <- mean(x,na.rm=TRUE)
  ##error <- sd(x,na.rm=TRUE)/sqrt(sum(is.finite(x))-1)
  return(y)
}

# Exceedance is a function that
#' @export
exceedance <- function(x,threshold=1,FUN='mean',...) UseMethod("exceedance")

#' @exportS3Method
#' @export exceedance.default
exceedance.default <- function(x,threshold=1,FUN='mean',na.rm=TRUE,...) {
  #print("HERE");  str(x)
  yrs <- year(x); d <- dim(x)
  X <- x; X[X <= threshold] <- NA
  # ns = number of stations
  if (is.null(d)) ns <- 1 else ns <- d[2]
  if ((FUN!="count") & (FUN!="freq")) {
    if ( (sum(is.element(names(formals(FUN)),'na.rm')==1)) |
         (sum(is.element(FUN,c('mean','min','max','sum','quantile')))>0 ) )
      y <- apply(matrix(X,length(X),ns),2,FUN,na.rm=na.rm, ...) else
        y <- apply(matrix(X,length(X),ns),2,FUN, ...)
      attr(y,'unit') <- attr(x,'unit')
  } else if (FUN=="count")  {
    #print("Wet-day count")
    y <- sum(is.finite(X))
    attr(y,'unit') <- paste("counts | X >",threshold,attr(y,'unit'))
  } else if (FUN=="freq") {
    #print(paste("Wet-day frequency",sum(is.finite(X)),sum(is.finite(x)),
    #            length(x),length(X),sum(is.finite(X))/sum(is.finite(x))))
    y <- sum(is.finite(X))/sum(is.finite(x))
    attr(y,'unit') <- paste("frequency | X >",threshold,attr(y,'unit'))
  }
  #str(y)
  #y <- attrcp(x,y)
  attr(y,'variable') <- paste(attr(x,'variable'),": exceedance above",threshold,
                              "-",FUN)
  attr(y,'history') <- history.stamp(x)
  return(y)
}

#' @exportS3Method
#' @export exceedance.station
exceedance.station <- function(x,threshold=1,FUN='mean',...) {
  y <- exceedance.default(x,threshold=threshold,FUN=FUN,...)
  return(y)
}

#' @exportS3Method
#' @export exceedance.field
exceedance.field <- function(x,threshold=1,FUN='mean',...) {
  y <- exceedance.default(x,threshold=threshold,FUN=FUN,...)
  #dimensions...
  return(y)
}


