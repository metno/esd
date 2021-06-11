#' Spell statistics
#' 
#' Statistics of spell durations (consecutive wet and dry days), e.g. dry and
#' wet periods or duration of extremes.
#' 
#' \code{exceedance} estimates statistics for peak-over-treshold, and
#' \code{nevents} returns the number of events with exceeding values (e.g. the
#' number of rainy days X > 1 mm/day).
#'
#' \code{wetfreq} returns the wet-day frequency (a fraction) and \code{wetmean} returns the wet-day mean.
#' 
#' \code{CDD}: Cooling degree day
#' \code{GDD}: Growing degree days (\url{http://en.wikipedia.org/wiki/Growing_degree-day})
#' \code{HDD}: Heating degree day
#'
#' \code{qqgeom} produces a quantile-quantile plot of streak statistics comparing the
#' empirical quantiles with the distribution function quantiles (see \code{\link[stats]{qgeom}}). 
#'
#' @aliases spell spell.default spell.station hist.spell count wetfreq wetmean
#' nevents exceedance exceedance.default exceedance.station exceedance.field
#' HDD CDD GDD qqgeom
#' @seealso hotsummerdays coldwinterdays coldspells heatwavespells nwetdays plot
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
#' # Example 1 : 
#' data(bjornholt)
#' x <- spell(bjornholt, threshold=.1)
#' x.ann <- annual(x, FUN="max")
#' plot(x.ann, plot.type="multiple",new=FALSE)
#' # Example 2 :
#' plot(x, new=FALSE)
#' 
#' # Growing degree days:
#' data(ferder)
#' plot(as.seasons(ferder,FUN='GDD'), new=FALSE)
#' 
#' # Mild winter days - number of days in the winter season with
#' # above freezing temperatures
#' try(coldwinterdays(ferder))
#'
#' # Quantile-quantile plot
#' qqgeom(ferder, treshold=1, pois=TRUE)
#'
#' @export spell
spell <- function(x,threshold,...) UseMethod("spell")

#' @exportS3Method
#' @export spell.default
spell.default <- function(x,threshold,upper=NULL,verbose=FALSE,...) {
  
  if (verbose) print('spell.default')
  
  ## Deal with missing data
  missing <- (1:length(x))[!is.finite(x)]
  if (min(missing)==1) {
    if (verbose) {
      print('strip away missing data at the beginning')
      print(paste('Data series length=',length(x),'number of missing data=',sum(missing)))
    }
    remove <- is.element(1:length(x),missing)
    # Make sure that the series does not start with missing data
    y <- zoo(x[!remove])
    class(y) <- class(x)
    y <- attrcp(x,y)
    y -> x ; rm('y')
  }
  ## Highligh the times when the values is above and below the given
  ## threshold:
  z <- coredata(x)
  mdate <- index(x)[!is.finite(x)]
  above <- z > threshold
  below <- z <= threshold
  above[!is.finite(above)] <- 0
  below[!is.finite(below)] <- 0
  if (verbose) print(paste('above=',sum(above,na.rm=TRUE),
                           'below=',sum(below,na.rm=TRUE)))
  
  ## Check if threshold is outside the range of data:
  if (sum(above,na.rm=TRUE)*sum(below,na.rm=TRUE)==0) {
    print(paste('The threshold',threshold,'is outside the range',
                paste(range(coredata(x,na.rm=TRUE)),collapse=' - ')))
    print("returning NULL")
    return(NULL)
  }
  
  ## Estimate the length of the streaks of above/below
  n <- length(index(x))
  t <- 1:n
  
  ## Use cumsum for calculating the number of days
  cath <- cumsum(above)
  cbth <- cumsum(below)
  ## Use diff to find the times when the streaks start and stop
  dt <- c(0,diff(above))
  if (verbose) print(table(dt))
  
  ## A streak starts the fist day when the value is above
  start <- t[dt > 0]
  ## A strak ends the day before the value is below 
  end <- t[dt < 0]-1
  
  ## Remove NA's:
  start <- start[is.finite(start)]
  end <- end[is.finite(end)]
  
  if (verbose) print(c(length(start),length(end)))
  
  ## Always start with a fresh spell
  if (start[1] > end[1]) {
    end <- end[-1]
  }
  ## Make sure that there are as many starts as endings
  if (length(start) > length(end)) {
    start <- start[1:length(end)]#-length(start)]
  }
  if (length(start) < length(end)) {
    end <- end[1:length(start)]#-length(end)]
  }
  
  chksum <- sum( (start - end) < 0)    
  
  if (verbose) {
    print(paste('Check sum=',chksum))
    print(c(length(start),length(end)))
    col <- c("blue","grey","red")
    dev.new()
    plot(t,cath,pch=19,cex=0.3,col=col[c(1+above-below)],...)
    points(t[start],cath[start],col='red',lwd=2,cex=0.7)
    points(t[end],cath[end],col='blue',lwd=2)
    dev.new(); hist(diff(start),col="blue")
    
    dev.new(); hist(diff(end),col="red")
    
    dev.new(); plot(start,end,pch=19,cex=0.3)
    lines(c(0,100000),c(0,100000),col=rgb(0.7,0.7,0.7,0.3))
  }
  
  ## Estimate the streaks:
  #  low <- t[start[-1]] - t[end[-length(end)]]
  #  high <- t[end] - t[start]
  high <- diff(cath[start])
  low <- diff(cbth[end])
  if (verbose) {print("high:"); print(summary(high))}  
  if (verbose) {print("low:"); print(summary(low))}
  
  ## If an upper limit is provided then ignore the long spells
  if (!is.null(upper)) {
    ignoreh <- high > upper 
    ignorel <- abs(low) > upper 
    high[ignoreh] <- NA
    low[ignorel] <- NA
  } else {
    ignoreh <- rep(FALSE,length(high))
    ignorel <- rep(FALSE,length(low))
  }
  
  Above <- zoo(high,order.by=index(x)[start])
  Below <- zoo(low,order.by=index(x)[end])
  
  ## Remove the spells with missing data:
  if (verbose) {
    print('remove spells when there are missing data')
    print(paste(sum(missing),'missing data')); print(summary(high)); print(summary(below))
  }
  if (verbose) print(paste('Above:',length(Above)))
  for (ii in 1:length(Above)) {
    fromd <- index(Above)[ii]; tod <- index(Above)[ii]+coredata(Above)[ii]
    if (!is.finite(tod)) tod <- fromd
    if (sum(is.element(seq(fromd,tod,by="1 day"),mdate)>0)) Above[ii] <- NA
  }
  if (verbose) print(paste('Below:',length(Below)))
  for (ii in 1:length(Below)) {
    fromd <- index(Below)[ii]; tod <- index(Below)[ii]+coredata(Below)[ii]
    if (!is.finite(tod)) tod <- fromd
    if (sum(is.element(seq(fromd,tod,by="1 day"),mdate)>0)) Below[ii] <- NA
  }
  
  #browser()
  y <- merge(Above,Below,all=TRUE)
  if (is.null(attr(x,'unit'))) attr(x,'unit') <- 'NA'
  
  if (verbose) {
    dev.new(); plot(y,main=paste(varid(x),'above/below',
                                 threshold,esd::unit(x),'at',loc(x)))
  }
  
  if (is.T(x)) {
    attr(y,'variable') <-  c("warm","cold") 
    attr(y,'longname') <-  c("duration of warm spells","duration of cold spells") 
  } else {
    attr(y,'variable') <-  c("wet","dry")
    attr(y,'longname') <-  c("duration of wet spells","duration of dry spells") 
  }
  attr(y,'location') <- rep(loc(x),2)
  attr(y,'station_id') <- rep(stid(x),2)
  attr(y,'longitude') <- rep(lon(x),2)
  attr(y,'latitude') <- rep(lat(y),2)
  attr(y,'altitude') <- rep(alt(x),2)
  attr(y,'unit') <- rep("days",2)
  attr(y,'threshold') <- rep(threshold,2)
  attr(y,'threshold.unit') <- rep(attr(x,'unit'),2)
  attr(y,'chksum') <- rep(chksum,2)
  attr(y,'uncredibly.high') <- t[ignoreh]
  attr(y,'uncredibly.low') <- t[ignorel]
  attr(y,'p.above') <- rep(sum(above)/length(above),2)
  attr(y,'interpolated.missing') <- index(x)[missing]
  class(y) <- c("spell",class(x))
  invisible(y)
}

#' @exportS3Method
#' @export spell.station
spell.station <-  function(x,threshold,upper=150,verbose=FALSE,...) {
  if (verbose) print('spell.station')
  if (!is.null(dim(x))) {
    if (verbose) print('group of stations')
    for (is in 1:dim(x)[2]) {
      y1 <- spell.default(subset(x,is=is),threshold=threshold,upper=upper,verbose=verbose,...)
      if (!inherits(y1,'spell')) {
        ## When failure, then create two series with NA
        y1 <- zoo(as.matrix(rep(NA,4),2,2),order.by=range(index(x)))
        if (is.T(x)) {
          attr(y1,'variable') <-  c("warm","cold") 
          attr(y1,'longname') <-  c("duration of warm spells","duration of cold spells") 
        } else {
          attr(y1,'variable') <-  c("wet","dry")
          attr(y1,'longname') <-  c("duration of wet spells","duration of dry spells") 
        }
        attr(y1,'location') <- loc(x)[is]
        attr(y1,'station_id') <- stid(x)[is]
        attr(y1,'longitude') <- lon(x)[is]
        attr(y1,'latitude') <- lat(y)[is]
        attr(y1,'altitude') <- alt(x)[is]
        attr(y1,'unit') <- "days"
        attr(y1,'threshold') <- rep(threshold,2)
        attr(y1,'threshold.unit') <- rep(attr(x,'unit'),2)
        attr(y1,'chksum') <- NA
        attr(y1,'uncredibly.high') <- NA
        attr(y1,'uncredibly.low') <- NA
        attr(y1,'p.above') <- NA
        attr(y1,'interpolated.missing') <- NA
        class(y1) <- c("spell",class(x))
      } 
      if (is==1) y <- y1 else y <- combine.stations(y,y1)
    }
  } else {
    ## Single station
    #
    missing <- (1:length(x))[!is.finite(x)]
    if (min(missing)==1) {
      if (verbose) {
        print('strip away missing data at the beginning')
        print(paste('Data series length=',length(x),'number of missing data=',sum(missing)))
      }
      remove <- is.element(1:length(x),missing)
      # Make sure that the series does not start with missing data
      y <- zoo(x[!remove])
      class(y) <- class(x)
      y <- attrcp(x,y)
      y -> x ; rm('y')
    }
    if (verbose) {print('Analyse spells'); str(x)}
    y <- spell.default(x,threshold=threshold,upper=upper,verbose=verbose,...)
    if (is.null(y)) return(y)
    y <- attrcp(x,y,ignore=c("variable","unit"))
    natr <- names(attributes(y))
    for (i in 1:length(natr)) 
      if (length(attr(y,natr[i]))==1) attr(y,natr[i]) <- rep(attr(y,natr[i]),2)
  }
  invisible(y)
}

#' @export
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

#' @exportS3Method
#' @export hist.spell
hist.spell <- function(x,family='geom',...) {
  n <- seq(0,ceiling(max(c(abs(x)),na.rm=TRUE))+1,by=1)
  hh <- hist(x[,1],breaks=n,plot=FALSE)
  hl <- hist(abs(x[,2]),breaks=n,plot=FALSE)
  
  #  dh <- dgeom(n,attr(x,'p.above')[1])
  #  dl <- dgeom(n,1-attr(x,'p.above')[2])
  if (substr(family,1,4)=='pois') {
    dh <- dpois(n,lambda=mean(x[,1],na.rm=TRUE))
    dl <- dpois(n,lambda=mean(abs(x[,2]),na.rm=TRUE))
  } else {
    dh <- dgeom(n,prob=1/(mean(x[,1],na.rm=TRUE)))
    dl <- dgeom(n,prob=1/(mean(abs(x[,2]),na.rm=TRUE)))
  }
  col <- c('red','blue')
  runs <- c('hot','cold')
  spelltype <- 'hot and cold'
  if (sum(is.element(attr(x,'variable'),c('wet','dry')))>0) {
    col <- c('darkblue','brown')
    runs <- c('wet','dry')
    spelltype <- 'wet and dry' 
  }
  
  main=paste(attr(x,'location')[1],spelltype,'spell duration')
  par(bty="n")
  plot(hh$mids,hh$density,type="s",col=col[1],lwd=3,
       ylab="",xlab="days",main=main,
       sub=paste("threshold=",attr(x,'threshold'),attr(x,'threshold.unit')))
  lines(hl$mids,hl$density,type="s",col=col[2],lwd=3)
  lines(n,dh,col=col[1],lty=2)
  lines(n,dl,col=col[2],lty=2)
  
  par(xaxt="n",yaxt="n",bty="n",fig=c(0.5,0.95,0.5,0.95),new=TRUE)
  plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
  legend(0.1,0.9,runs,col=col,lty=1,lwd=3,bty="n")
  
}

#' @export
qqgeom <- function(x,treshold=1,pois=FALSE,...) {
  s <- spell(x,threshold=treshold)
  x1 <- qgeom(seq(0,1,length=101),prob=1/(mean(coredata(s[,1]),na.rm=TRUE)))
  y1 <- quantile(as.numeric(s[,1]),probs=seq(0,1,length=101),na.rm=TRUE)
  x2 <- qgeom(seq(0,1,length=101),prob=1/(mean(coredata(s[,2]),na.rm=TRUE)))
  y2 <- quantile(as.numeric(s[,2]),probs=seq(0,1,length=101),na.rm=TRUE)
  xp1 <- qpois(seq(0,1,length=101),lambda=mean(coredata(s[,1]),na.rm=TRUE))
  xp2 <- qpois(seq(0,1,length=101),lambda=mean(coredata(s[,2]),na.rm=TRUE))
  
  par(bty='n')
  xy <- c(x1,x2,xp1,xp2,y1,y2,y1,y2); ok <- is.finite(xy)
  plot(range(xy[ok]),range(xy[ok]),
       type='l',main='q-q geometric of streak statistics',
       xlab='distribution function',ylab=expression(q[p]))
  if (pois) {
    points(xp1,y1,pch=15,col=rgb(0.7,0.7,1))
    points(xp2,y2,pch=15,col=rgb(1,0.7,0.7))
  }
  points(x1,y1,pch=19,col=rgb(0.2,0.2,1))
  points(x2,y2,pch=19,col=rgb(1,0.2,0.2))
  grid()
  
  if (pois) legend(min(xy[ok]),max(xy[ok]),
                   c(paste('geom.',varid(s)[1]),paste('geom.',varid(s)[2]),
                     paste('pois.',varid(s)[1]),paste('pois.',varid(s)[2])),
                   col=c(rgb(0.2,0.2,1),rgb(1,0.2,0.2),rgb(0.7,0.7,1),rgb(1,0.7,0.7)),
                   pch=c(rep(19,2),rep(15,2)),bty='n') else
                     legend(min(xy[ok]),max(xy[ok]),
                            c(paste('geom.',varid(s)[1]),paste('geom.',varid(s)[2])),
                            col=c(rgb(0.2,0.2,1),rgb(1,0.2,0.2)),
                            pch=rep(19,2),bty='n')
}

#' @export
HDD <- function(x,x0=18,na.rm=TRUE) {
  cold <- x < x0
  hdd <- sum(x0 - x[cold],na.rm=na.rm)
  return(hdd)
}


#' @export
CDD <- function(x,x0=22,na.rm=TRUE) {
  warm <- x > x0
  cdd <- sum(x[warm] - x0,na.rm=na.rm)
  return(cdd)
}

#' @export
GDD <- function(x,x0=10,na.rm=TRUE) {
  gdd <- CDD(x,x0=x0)
  attr(gdd,'url') <- 'http://en.wikipedia.org/wiki/Growing_degree-day'
  return(gdd)
}
