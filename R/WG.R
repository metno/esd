# This script provides a demonstration of how a weather generator can be designed:
# ESD is used to predict thanges in the pdfs for temperature and precipitation.
# Temperature: mean, stdv.
# Precipitation: mu, fw (assume similar spell statistics as in present time).
#
# (more advanced method may perhaps predict changes to the spell-statistics).
#
# Rasmus Benestad

#' Weather generators for conditioned on simulated climate aggregated
#' statistics.
#' 
#' Weather generators for conditional simulation of daily temperature and/or
#' precipitation, given mean and/or standard deviation. The family of WG
#' functions produce stochastic time series with similar characteristics as the
#' station series provided (if none if provided, it will use either ferder or
#' bjornholt provided by the esd-package). Here characteristics means similar
#' mean value, standard deviation, (frequency and wet-day mean precipitation for 
#' precipitation), and spectral properties. \code{FTscramble}
#' takes the Fourier components (doing a Fourier Transform - FT) of a series 
#' and reassigns random phase to each frequency and then returns a new series 
#' through an inverse FT. The FT scrambling is used for temperature, but not for 
#' daily precipitation that is non-Gaussian and involves sporadic events with rain. 
#' Instead, \code{FTscramble} is used for randomly scrambling annual wet-day frequency 
#' and wet-day mean precipitation, however, and for precipitation, the annual wet-day
#' frequency and the annual wet-day mean precipitation are used to randomly generate
#' exponentially distributed numbers to provide similar aggregated annual
#' statistics as the station or predicted though downscaling. The precipitation
#' WG can also take into account the number of consecutive number-of-wet-days
#' statistics using a geometric distribution.
#' 
#' The weather generator produces a series with similar length as the provided
#' sample data, but with shifted dates according to specified scenarios for
#' annual mean mean/standard deviation/wet-day mean/wet-day frequency.
#' 
#' \code{WG.FT.day.t2m} generates daily temperature from seasonal means and
#' standard deviations. It is given a sample station series, and uses
#' \code{FTscramble} is based on a Fourier Transform which generates a new series 
#' with random phase but nevertheless similar (or predicted - in the future) spectral 
#' characteristics as the original series. It then uses a quantile
#' transform to prescribe predicted mean and standard deviation, assuming the
#' distributions are normal, which usually is OK for seasonall/annually aggregated statistics 
#' (e.g. annual mean, annual wet-day frequency, or annual wet-day mean precipitation). 
#' The temporal structure (power spectrum) of the random series is therefore similar 
#' as the sample provided.
#' 
#' \code{WG.fwmu.day.precip} has been designed to be used with downscaled results for 
#' annual wet-day frequency and annual wet-day mean precipitation. It also tries to
#' simulate the wet-spell duration statistics (number of consecutive wet days) based on
#' provided sample data (argument x). The process can take annual wet-day mean precipitation 
#' and the wet-day frequency as input when used to make projections for the future, together 
#' with a sample station of daily values, to simulate stochastic numbers of consecutive 
#' wet days, based on its annual mean number of consecutive wet days. It also uses 
#' the mean annual cycle of wet-day frequency as well as the wet-day mean precipitation 
#' to guide the seasonal timing of wet days and amounts, and hence tries to mimic rain seasons. 
#' If not specified, it is taken from the sample data after being phase scrambled/shuffled 
#' (\code{FTscramble} - a bit like a deck of cards). If not specified, the annual wet-day 
#' frequency is a phase-scrambled version of annual aggregates from the sample data. The daily 
#' amount is taken from stochastic values generated with \code{\link{rexp}} scaled 
#' for the tail according to alpha in (described in DOI: 10.1088/1748-9326/ab2bb2) 
#' as in \code{\link{day2IDF}}. The number of consecutive wet days is approximated 
#' by a geometric distribution (\code{\link{rgeom}}), and the annual number of wet days 
#' is either given as input or estimated from the sample series. 
#' \code{test.WG.fwmu.day.precip} presents diagnostics of tests of \code{WG.fwmu.day.precip}.
#' 
#' @aliases WG WG.station WG.fwmu.day.precip WG.FT.day.t2m FTscramble
#'
#' @importFrom stats start end approx pnorm qnorm qqnorm sd dgeom rgeom rexp qexp pexp dpois
#' fft runif
#' @importFrom graphics hist
#'
#' @param x station object
#' @param option Define the type of WG
#' @param amean annual mean values. If NULL, use those estimated from x; if NA,
#' estimate using \code{\link{DSensemble.t2m}}, or if provided, assume a
#' 'dsensemble' object.
#' @param asd annual standard deviation. If NULL, use those estimated from x;
#' if NA, estimate using \code{\link{DSensemble.t2m}}, or if provided, assume a
#' 'dsensemble' object.
#' @param t Time axis. If null, use the same as x or the last interval of same
#' length as x from downscaled results.
#' @param ip passed on to \code{\link{DSensemble.t2m}}
#' @param select passed on to \code{\link{DSensemble.t2m}}
#' @param lon passed on to \code{\link{DSensemble.t2m}}
#' @param lat passed on to \code{\link{DSensemble.t2m}}
#' @param plot if TRUE, plot results
#' @param biascorrect passed on to \code{\link{DSensemble.t2m}}
#' @param verbose passed on to \code{\link{DSensemble.t2m}}
#' @param mu annual wet-mean values. If NULL, use those estimated from x; if
#' NA, estimate using \code{\link{DSensemble.t2m}}, or if provided, assume a
#' 'dsensemble' object.
#' @param fw annual wet-day frequency. If NULL, use those estimated from x; if
#' NA, estimate using \code{\link{DSensemble.t2m}}, or if provided, assume a
#' 'dsensemble' object.
#' @param ndd annual mean dry spell length. If NULL, use those estimated from
#' x; if NA, estimate using \code{\link{DSensemble.t2m}}, or if provided,
#' assume a 'dsensemble' object.
#' @param threshold Definition of a rainy day.
#' @param method Assume a gemoetric or a poisson distribution. Can also define
#' ownth methods.
#' @param t2m station object with temperature
#' @param precip station object with precipitation.
#' @param ndbr Number of 
#' @param n.spells.year = c('fw','spell') if 'fw' then estimate number of spells according to 365.25 else estimate number of events from \code{\link{spell}}.
#' @param alpha.scaling TRUE scale the low-probability events according to alpha in DOI:10.1088/1748-9326/ab2bb2
#' @param alpha values for alpha-scaling 
#' @param ensure.fw TRUE then WG tries to ensure that fw of simulations match those of observations or prescribed by adding or subtracting wet days.
#' @param \dots additional arguments
#' @author R.E. Benestad
#' @keywords manip
#' @examples
#' ## Temperature
#' data(ferder)
#' x <- WG(ferder)
#' ## Plot the results
#' plot(merge(ferder,x), xlab='', ylab=c('Obs T2m','WG T2m'), col='blue', main=paste(loc(x),' Obs/WG'))
#' 
#' ## Daily precipitation
#' data(bjornholt)
#' z <- WG(bjornholt)
#' ## Plot the results
#' plot(merge(bjornholt,z), xlab='', ylab=c('Obs precip','WG precip'), col='blue', main=paste(loc(z),' Obs/WG'))
#' sz <- sort(coredata(z)[index(z) %in% index(bjornholt)])
#' sy <- sort(coredata(bjornholt)[index(bjornholt) %in% index(z)])
#' ## Use WG to 'simulate' climate change
#' z2 <- WG(bjornholt, mu=annual(bjornholt, FUN='wetmean') + 2)
#' sz2 <- sort(coredata(z2)[index(z2) %in% index(bjornholt)])
#'
#' ## Plot the comparison of quantiles for fitted WG and simulated climate change
#' plot(sy, sz, pch=19, cex=0.7, main='QQ-plot', xlab='Observations', ylab='WG')
#' grid()
#' lines(c(0, max(sy,sz,na.rm=TRUE)), c(0,max(sy,sz,na.rm=TRUE)), lty=2, col='red')
#' points(sy, sz2, col='blue', cex=0.7)
#' legend('topleft',c('fitted WG','simulated change'),col=c('black','blue'),pch=21,bty='n')
#' 
#' 
#' ## Simple simulation of continued trends in wet-day mean precipitation and frequency
#' mu <- annual(bjornholt,FUN='wetmean',nmin=270) # Avoid missing values (NA)
#' fw <- annual(bjornholt,FUN='wetfreq',nmin=270) # Avoid missing values (NA)
#' mu.trend <- trend(mu)
#' fw.trend <- trend(fw)
#' ## Construct precipitation statistics for input to WG
#' mu2 <- c(mu,zoo(coredata(mu)+coredata(max(mu.trend)-min(mu.trend)),order.by=max(index(mu))+1:length(mu)))
#' fw2 <- c(fw,zoo(coredata(fw)+coredata(max(fw.trend)-min(fw.trend)),order.by=max(index(fw))+1:length(fw)))
#' z <- WG(bjornholt,mu=mu2,fw=fw2,verbose=TRUE)
#' plot(z)
#' 
#' #' ## Test the WG
#' z2 <- WG(bjornholt,plot=TRUE,verbose=TRUE)
#' plot(aggregate(z2,by=month,FUN='wetmean')); lines(aggregate(bjornholt,by=month,FUN='wetmean'))
#' z3 <- WG(bjornholt,plot=TRUE,verbose=TRUE)
#' plot(aggregate(z3,by=month,FUN='wetfreq')); lines(aggregate(bjornholt,by=month,FUN='wetfreq'))
#' 
#' ## Test-routine for WG
#' test.WG.fwmu.day.precip()
#' 
#' ## Demonstrate bi-variate
#' precip <- station(param='precip',stid=18700,src='metnod.thredds')
#' t2m <- station(param='t2m',stid=18700,src='metnod.thredds')
#' hist.2D <- bivariate.hist(precip,t2m,plot=TRUE)
#' @export WG
WG <- function(x,...) UseMethod("WG")

#' @exportS3Method
#' @export WG.station
WG.station <- function(x,...,option='default') {
  if (inherits(x,'day')) {
    if (length(varid(x))==1) {
      if (varid(x)=='t2m') y <- WG.FT.day.t2m(x,...) else
        if (varid(x)=='precip') y <- WG.fwmu.day.precip(x,...)
    }
  }
  return(y)
}

#' @exportS3Method
#' @export WG.FT.day.t2m
WG.FT.day.t2m <- function(x=NULL,...,amean=NULL,asd=NULL,t=NULL,
                          plot=FALSE,verbose=FALSE) {
  if (verbose) print('WG.FT.day.t2m')
  ## Single function for just temperature.
  ## The arguments amean and asd are time series predicted through ESD or
  ## adopted from a zoo or station object (x). Typically annually or seasonally
  ## aggregated data
  if (is.null(x)) {
    ## If no stations objects is given, use default 
    if (verbose) print("use default: Ferder, Norway")
    ferder <- NULL
    data("ferder",envir=environment())
    x <- ferder
    rm('ferder')
  }
  
  ## Get the daily anomalies and the climatology
  xa <- anomaly(x); clim <- x - xa
  
  ## KMP 2024-05-31: If amean or asd are NULL, the function fails! Is amean and asd supposed to be the
  ## annual mean and standard deviation of x? Adding a check to see if amean and asd exist and if not calculate them based on x.
  if(is.null(amean)) amean <- as.annual(x, FUN="mean", na.rm=TRUE)
  if(is.null(asd)) asd <- as.annual(x, FUN="sd", na.rm=TRUE)
  
  ## Define time axis for projection based on the annual mean data either from station or
  ## downscaled projections
  if (is.null(t)) {
    if (verbose) print("set the time index")
    ly <- max(year(amean)); ny <- length(rownames(table(year(amean)))) 
    interval <- c(ly-ny+1,ly)
    if(verbose) print(interval)
    t <- seq.Date(as.Date(paste(interval[1],substr(start(x),5,10),sep='')),
                  as.Date(paste(interval[2],substr(end(x),5,10),sep='')),
                  by="day")
    #browser()
    #str(t); print(paste(interval[1],month(x)[1],day(x)[1],sep='-'))
    #t <- t - julian(t[1]) +
    #  julian(as.Date(paste(interval[1],month(x)[1],day(x)[1],sep='-')))
  }
  
  ## Estimate a smooth curve for the annual mean and standard deviation that has a daily resolution
  if (verbose) print("Estimate smooth day-by-day changes in mean and sd:")
  ym <- approx(julian(as.Date(index(amean))),coredata(amean),xout=julian(as.Date(t)),rule=2)$y
  #print(summary(ym))
  ys <- approx(julian(as.Date(index(asd))),coredata(asd),xout=julian(as.Date(t)),rule=2)$y
  
  ## New object y that contains random variable as original data but with same spectral 
  ## characteristics and same climatology
  if (verbose) print("Construct a station object with random timing but original time structure:")
  y <- zoo(FTscramble(xa,t),order.by=t)
  ## Add climatology to get original type data
  if (verbose) print("add climatology")
  y <- y + matchdate(clim,y)
  
  if (plot) {
    dev.new()
    plot(merge(zoo(xa),zoo(anomaly(y))),plot.type='single',lwd=c(2,1),
         col=c('black','grey'))
  }
  
  ## qq-transform to transform the temperature distribution from present shape to future shape
  ## assuming a normal distribution: ~N(m1,s1) -> ~N(m2,s2). Estimate probabilities based on the 
  ## scrambeled series y and use these probabilities to derive new quantiles based on the shifted 
  ## pdf.
  cdf <- pnorm(q=y,mean=mean(y,na.rm=TRUE),sd=sd(y,na.rm=TRUE))
  q2 <- qnorm(cdf,mean=ym,sd=ys)
  #print(summary(cdf)); print(summary(q2))
  #hist(cdf); browser()
  z <- zoo(q2,order.by=t)
  #print(summary(z))
  if (verbose) print("Attach attributes")
  z <- attrcp(x,z)
  class(z) <- c('WG',class(x))
  attr(z,'mean') <- ym
  attr(z,'sd') <- ys
  attr(z,'aspect') <- paste(attr(z,'aspect'),'weather_generator',sep=', ')
  attr(z,'history') <- history.stamp(x)
  return(z)
}



## Use a simpler approach for wet-day mean precipitation mu - 
## Estimate daily climatology and smooth daily weights against the annual mean
## Then use these climatological weights after the amounts have been dealt out
## to ensure that mu has the same daily climatology.
## --- Precipitation 
#' @exportS3Method
#' @export WG.fwmu.day.precip
#' @exportS3Method
#' @export WG.fwmu.day.precip
WG.fwmu.day.precip <- function(x=NULL,...) {
  ## Argument x is a station object with daily data
  ## Collect the arguments passed on with ...
  args <- list(...)
  plot <- args$plot; if (is.null(plot)) plot <- FALSE
  verbose <- args$verbose;  if (is.null(verbose)) verbose <- FALSE
  mu=args$mu
  fw=args$fw
  t=args$t
  start <- args$start
  threshold <- args$threshold; if (is.null(threshold)) threshold <- 1
  ## Use alpha scaling estimates from DOI:10.1088/1748-9326/abd4ab - same as in ERL::IDF()
  ## The scaling for return values: alpha.scaling = 1.26 + 0.06 log(year)
  ## Supporting material in DOI:10.1088/1748-9326/ab2bb2
  ## 1 year return-period corresponds to a quantile of q(p) = q(1 - 1/365)
  ## estimate probabilities for the amounts and multiply by appropriate scaling factor by translating
  ## the above expression to the probability 
  
  ## tau return.year <- Pr = 1/(tau*365.25)
  ## Use pexp(amount) to estimate return.year and then the lin-log relation to estimate scaling factor 
  ## for each rainy day. Minimum return.year is set to 1/365.25 as a minimum possible return year 
  ## (every day) for which the formula will give 0.88
  x.tau <- seq(0,50,by=0.01)
  scaling <- function(x) { ## Table 1 from DOI:10.1088/1748-9326/ab2bb2
    y <- 1.256 + 0.064*log(x) 
    y[!is.finite(y)] <- 1
    return(y)
  } 
  
  if (verbose) print('WG.fwmu.day.precip')
  # Single function for just precipitation
  if (is.null(x)) {
    if (verbose) print('Use sample data from esd')
    data("bjornholt",envir=environment())
    x <- bjornholt
    rm('bjornholt')
  } else if (verbose) print(paste('Use data provided:',loc(x)))
  
  ## Estimate the climatologies in fw and mu
  fw.clim <- aggregate(x,by=month,FUN='wetfreq',threshold=threshold)
  mu.clim  <- aggregate(x,by=month,FUN='wetmean',threshold=threshold)
  ## Define the annual fw and mu
  x.fw <- annual(x,FUN='wetfreq',threshold=threshold,start=start)
  x.mu <- annual(x,FUN='wetmean',threshold=threshold,start=start)
  
  if (is.null(fw)) 
    fw <- zoo(FTscramble(x.fw),order.by=index(x.fw)) else
      ## fw is introduced as a change factor
      if (length(fw)==1) {
        fw <- fw + zoo(FTscramble(x.fw),order.by=index(x.mu))
      }
  if (is.null(mu)) mu <- zoo(FTscramble(x.mu),order.by=index(x.mu)) else
    ## mu is introduced as a change factor
    if (length(mu)==1) {
      mu <- mu + zoo(FTscramble(x.mu),order.by=index(x.mu))
    } 
  rm('x.fw','x.mu')
  
  ## Define the time axis
  if (is.null(t)) t <- index(x)
  ## Define the years
  yrs <- rownames(table(year(t)))
  if (verbose) cat('The WG simulates the years: ',range(yrs),' \n')
  ## For each year - loop
  Z <- NULL
  for (it in 1:length(yrs)) { 
    ## fw gives number of wet days per year
    n.wet <- round(fw[it] * 365.25)
    ## mu gives the amounts
    ## since we sort the data according to magnitude, we need to draw just enough numbers
    ## to correspond to the exponential distribution
    amount <- sort(rexp(n.wet,rate=1/mu[it]),decreasing=TRUE)
    p.amount <- 1 - pexp(amount,rate=1/mu[it])
    ## The scaling was defined for return values x_tau = scaling * mu * ln(tau*fw)
    ## where tau is the return period (years) and is 1/(n_tau * 365.26)
    tau.amount <- 1/(coredata(fw[it])*p.amount*365.25)
    ## Interpolate the scaling from the linear-log expression for scaling
    scaling.amount <- approx(x=x.tau,y=scaling(x.tau),xout=tau.amount)$y
    amount <- amount * scaling.amount
    ## Use fw.clim to estimate the probability of a wet day and normalise so that it theoretically gives 
    ## N.wet wet days in a year. Sort from highest to lowest probability/frequency, keeping track of
    ## the Julian days to avoid a gap at the end of the year
    ndaysthisyear <- as.numeric(as.Date(paste0(yrs[it],'-12-31'))-as.Date(paste0(yrs[it],'-01-01')))+1
    jdays <- seq(1,ndaysthisyear,by=1)
    ## Use mu.clim to deal out appropriate amounts on the wet days, drawing from amount. Base it on
    wPr <- approx(x=seq(1,ndaysthisyear,length=12),y=coredata(fw.clim),xout=jdays)$y
    ## Weight  the wet-probability so that it matches that of the annual wet-day frequency
    wPr <- wPr* coredata(fw)[it]/mean(wPr)
    ## First guess: 
    
    wet <- runif(ndaysthisyear) < wPr
    
    ## Add or remove wet days so that the total number of wet days corresponds with fw
    while (sum(wet) != n.wet) {
      #if (verbose) cat(n.wet,' == ',sum(wet),'? \n')
      if (sum(wet) > n.wet) {
        ## Find a random wet day and flip it to dry
        w2d <- jdays[wet]; pick <- w2d[order(rnorm(length(w2d)))][1]
        wet[pick] <- FALSE
      }
      if (sum(wet) < n.wet) {
        ## Find a random wet day and flip it to dry
        d2w <- jdays[!wet]; pick <- d2w[order(rnorm(length(d2w)))][1]
        wet[pick] <- TRUE
      }
    } 
    ## Rank typical intensity variation with the season based on mu.clim and the Julian day.
    mu.jday <- approx(x=seq(1,ndaysthisyear,length=12),y=coredata(mu.clim),xout=jdays)$y
    ## For wet days, deal out the amount according to mu.clim
    ## Daily precipitation this year
    z <- rep(0,ndaysthisyear)
    itmu <- order(mu.jday[wet] + 2*sd(mu.clim)*rnorm(sum(wet)),decreasing=TRUE)
    z[wet] <- amount[order(itmu)]
    ## Ensure that the year has the same wet-day mean as prescribed
    z[wet] <- round(z[wet]*coredata(mu)[it]/mean(z[wet],na.rm=TRUE),1)
    z <- zoo(z,order.by=seq(as.Date(paste0(yrs[it],'-01-01')),as.Date(paste0(yrs[it],'-12-31')),by=1))
    if (is.null(Z)) Z <- z else Z <- c(Z,z)
    if (verbose) cat(it, yrs[it], 'fw=',fw[it], 'n.wet=',n.wet,'=',sum(wet),
                     'mu=',round(mu[it],1),'=', round(mean(z[wet]),1),'\n')
  }  ## end of loop over the years
  class(Z) <- class(x)
  Z <- attrcp(x,Z)
  attr(Z,'info') <- 'Weather generator'
  attr(Z,'history') <- history.stamp()
  invisible(Z)
}

#' This function tests the WG for precipitation:
#' Quantile-quantile plots of wet-day amounts
#' Number of wet days
#' @exportS3Method
#' @export test.WG.fwmu.day.precip
test.WG.fwmu.day.precip <- function(x=NULL,verbose=TRUE) {
  if (verbose) print('test.WG.fwmu.day.precip')
  if (is.null(x)) {data('bjornholt',envir=environment()); x <- bjornholt; rm('bjornholt')}
  print(paste('test.WG.fwmu.day.precip for',loc(x)))
  z <- WG(x,verbose=verbose)
  ## sort magnitudes to plot quantile-quantile plots
  if (verbose) print('Sort precipitation magnitudes')
  xw <- sort(coredata(x)[x > 1])
  zw <- sort(coredata(z)[z > 1])
  ## There may be different number of wet days - pad the shortest series with 0s
  nx <- length(xw)
  nz <- length(zw)
  print(paste(nx,'observed wet days and',nz,'simulated wet days'))
  print('Obs:');print(summary(x)); print('WG:');print(summary(z))
  if (nx > nz) zw <- c(rep(0,nx-nz),zw) else 
    if (nz > nx) xw <- c(rep(0,nz-nx),xw)
  nx <- length(xw); nz <- length(zw)

  xylim <- c(0,max(c(xw,zw)))
  par(mfcol=c(2,2))
  plot(xw,zw,ylim=xylim,xlim=xylim,xlab='Observed amount (mm/day)',
       ylab='WG amount (mm/day)',main=paste(loc(x),'wet-day amounts'))
  grid()
  lines(xylim,xylim,lty=2,col='red')
  ## compare the number of wet days
  xnw <- zoo(annual(x,FUN='count',1))
  znw <- zoo(annual(z,FUN='count',1))
  plot(merge(xnw,znw),plot.type='single',col=c('black','red'),lty=c(1,2),
       main='Number of annual wet days',ylab='days',xlab='')
  legend('topleft',c('Original','WG'),col=c('black','red'),lty=c(1,2),bty='n',
         cex=0.6)
  grid()
  ## Compare the spell duration statistics
  sx <- spell(x,1)
  sz <- spell(z,1)
  dryx <- coredata(sx[,1]); dryx <- sort(dryx[is.finite(dryx)])
  dryz <- coredata(sz[,1]); dryz <- sort(dryz[is.finite(dryz)])
  wetx <- coredata(sx[,2]); wetx <- sort(wetx[is.finite(wetx)])
  wetz <- coredata(sz[,2]); wetz <- sort(wetz[is.finite(wetz)])
  if (length(dryx) < length(dryz)) dryx <- c(rep(0,length(dryz)-length(dryx)),dryx) else 
    if (length(dryx) > length(dryz)) dryz <- c(rep(0,length(dryx)-length(dryz)),dryz)
  xylim <- c(1,max(c(dryx,dryz),na.rm=TRUE))
  ## Very long spells are few and more influenced by randomness
  #dryx[dryx>30] <- NA; dryz[dryz>30] <- NA
  plot(dryx,dryz,main='Dry spell durations',xlab='obs',ylab='WG',
       xlim=xylim,ylim=xylim)
  grid()
  lines(xylim,xylim,lty=2,col='red')
  if (length(wetx) < length(wetz)) wetx <- c(rep(0,length(wetz)-length(wetx)),wetx) else 
    if (length(wetx) > length(wetz)) wetz <- c(rep(0,length(wetx)-length(wetz)),wetz)
  xylim <- c(1,max(c(wetx,wetz),na.rm=TRUE))
  ## Very long spells are few and more influenced by randomness
  #wetx[wetx>30] <- NA; wetz[wetz>30] <- NA
  plot(wetx,wetz,main='wet spell durations',xlab='obs',ylab='WG',
       xlim=xylim,ylim=xylim)
  grid()
  lines(xylim,xylim,lty=2,col='red')
  
  ## Compare the mean annual cycle of observations and simulations
  zx <- combine.stations(x,z)
  col <- c('black','red')
  ## Compare annual statistics
  plot(zoo(annual(zx,FUN='sum')),main='Annual total precipitation',col=col,
       plot.type='single',ylab=expression(sum(x)*phantom(0)*(mm/day)),lty=c(1,2))
  legend('topleft',c('Original','WG'),col=c('black','red'),lty=c(1,2),bty='n',
         cex=0.6)
  grid()
  ## Compare mean seasonal 
  plot(zoo(aggregate(zx,by=month,FUN='sum')),main='Seasonal total precipitation',
       ylab=expression(sum(x)*phantom(0)*(mm/day)),col=col,plot.type='single',lty=c(1,2)); grid()
  legend('topleft',c('Original','WG'),col=c('black','red'),lty=c(1,2),bty='n',
         cex=0.6)
  plot(zoo(aggregate(zx,by=month,FUN='wetfreq')),main='Seasonal wet-day frequency',
       ylab=expression(f[w]),col=col,plot.type='single',lty=c(1,2)); grid()
  legend('topleft',c('Original','WG'),col=c('black','red'),lty=c(1,2),bty='n',
         cex=0.6)
  plot(zoo(aggregate(zx,by=month,FUN='wetmean')),main='Seasonal wet-day mean',
       ylab=expression(mu*phantom(0)*(mm/day)),col=col,plot.type='single',lty=c(1,2)); grid()
  legend('topleft',c('Original','WG'),col=c('black','red'),lty=c(1,2),bty='n',
         cex=0.6)
  invisible(merge(x,z))
}


#' @export
FTscramble <- function(x,t=NULL,interval=NULL,spell.stats=FALSE,
                       wetfreq.pred=FALSE) {
  attributes(x) <- NULL
  n <- length(x)
  
  # This function scrambles the phase information of the FT components of a
  # time series, maintaining the same spectral and time structure
  
  if (sum(is.na(x))>0) {
    ## If there are some missing data, fill in with interpolated values
    ok <- is.finite(x)
    if (sum(ok)>2) y <- approx((1:n)[ok],x[ok],xout=1:n,rule=2)$y else {
      stop(paste('FTscramble - only',sum(ok),'valid data points'))
    }
    x <- y
    rm('y')
  }
  
  # Fourier transform (FT) to obtain power and phase information
  X <- fft(x)
  #print(summary(Re(X))); print(summary(Im(X)))
  
  # Z contains the phase information
  Z <- Mod(X)
  #print(summary(Z))
  ReX <- Re(X)
  ImX <- Im(X)
  phiX <- Arg(X)
  #plot(phiX)
  
  # Set new phase information to random
  phiY <- runif(n,min=-pi,max=pi)
  ReY <- Z*cos(phiY)
  ImY <- Z*sin(phiY)
  ReY[1] <- ReX[1]; ImY[1] <- ImX[1]
  #  ReY[n] <- ReX[n]; ImY[n] <- ImX[n]
  Y <- complex(real=ReY, imaginary=ImY)
  
  # Inverse FT to generate new time series:
  y <- Re(fft(Y,inverse=TRUE))/n
  # Make sure that the new scrambled series has the same mean and standard deviation
  # as the original data:
  y <- sd(x,na.rm=TRUE)*(y - mean(y,na.rm=TRUE))/sd(y,na.rm=TRUE) + mean(x,na.rm=TRUE)
  
  if (!is.null(t)) {
    if (length(t) <= length(y)) y <- y[1:length(t)] else
      y <- c(y,rep(NA,length(t)-length(y)))
  }
  invisible(y)
}

## This weather generator assumes daily precipitation and temperature and samples from a
## joint statistical distribution. Since precipitation involves dry and wet days, it uses
## the quantiles from precipitation to slice the 2D-distribution function to re-sample
## temperature (temperatures are provided for all the time whereas it only rains in
## occasions). Use mean and standard deviation to standardise temperature and apply 
## new mean and standard deviation to simulate changed temperature. The function returns 
## the original daily precipitation series. 
#' @exportS3Method
#' @export WG.bivariate
WG.bivariate <- function(x,y,...) {
  z.x <- WG(x,...)
  z.y <- WG(y,...)
  Z.xy <- bivariate.hist(x,y) 
  sx <- stand(x)
}

#' @export
bivariate.hist <- function(x,y,plot=TRUE) {
  ## Standardising
  stand <- function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  
  ## Keep a copy
  pr <- x; tas <- y
  ## Ensure rigth format: zoo-objects
  x <- zoo(x)
  y <- stand(zoo(y))
  
  ## Make 2D bins and count events in each 
  nx = 30; ny <- nx
  
  seqx <- seq(floor(min(x)),ceiling(max(x)),length=nx+1)
  seqy <- seq(floor(min(y)),ceiling(max(y)),length=ny+1)
  seqX <- seq(floor(min(x)),ceiling(max(x)),length=10*nx)
  seqY <- seq(floor(min(y)),ceiling(max(y)),length=10*ny)
  Z <- matrix(rep(NA,nx*ny),nx,ny)
  for (i in 1:nx)
    for (j in 1:ny) Z[i,j] <- sum( (x >= seqx[i]) & (x < seqx[i+1]) &
                                     (y >= seqy[j]) & (y < seqy[j+1]) )
  ## Make a smoother 2D surface with higher resolution
  ## Start with x dimension
  Zx <- matrix(rep(NA,10*nx*ny),length(seqX),ny)
  for (j in 1:ny) Zx[,j] <- approx(1:nx,Z[,j],xout=seq(1,nx,length=10*nx))$y
  Zxy <- matrix(rep(NA,100*nx*ny),length(seqX),length(seqY))
  for (i in seq(1,10*nx,by=1)) Zxy[i,] <- approx(1:ny,Zx[i,],xout=seq(1,ny,length=10*ny))$y
  
  seqY <- seqY*sd(y,na.rm=TRUE) + mean(y,na.rm=TRUE)
  if (plot) { 
    image(seqX,seqY,log(Zxy),main='Bivariate statistical distribution',
          xlab='Precipitation (mm/day)',ylab='Temperature (C)')
    contour(seqX,seqY,log(Zxy),add=TRUE)
    points(x,y,pch=19,cex=0.25,col=rgb(0,0,0,0.1))
    grid()
  }
  attr(Zxy,'x') <- seqX
  attr(Zxy,'y') <- seqY
  invisible(Zxy)
}
