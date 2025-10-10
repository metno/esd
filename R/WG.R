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
#' @aliases WG WG.station WG.fwmu.day.precip WG.FT.day.t2m
#' WG.pca.day.t2m.precip FTscramble
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
#' ## Plot the comparison of quantiles
#' plot(sy, sz, pch=19, cex=0.7, main='QQ-plot', xlab='Observations', ylab='WG')
#' grid()
#' lines(c(0, max(sy,sz,na.rm=TRUE)), c(0,max(sy,sz,na.rm=TRUE)), lty=2, col='red')
#' points(sy, sz2, col='blue', cex=0.7)
#' 
#' 
#' ## Simple simulation of contnued trends in wet-day mean precipitation and frequency
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
WG.fwmu.day.precip <- function(x=NULL,mu=NULL,fw=NULL,t=NULL,...) {
  ## Argument x is a station object with daily data
  ## Collect the arguments passed on with ...
  args <- list(...)
  plot <- args$plot; if (is.null(plot)) plot <- FALSE
  verbose <- args$verbose;  if (is.null(verbose)) verbose <- FALSE
  # mu=args$mu
  # fw=args$fw
  # t=args$t
  threshold <- args$threshold; if (is.null(threshold)) threshold <- 1
  alpha.scaling <-args$alpha.scaling
  if (is.null(alpha.scaling)) alpha.scaling <- TRUE
  ## Use alpha scaling estimates from DOI:10.1088/1748-9326/abd4ab - same as in ERL::IDF()
  alpha <-args$alpha; if (is.null(alpha)) alpha=c(1.256,0.064)
  ## Weighting function to determine the degree which the mean seasonal cycle determines the results 
  ncd.max <- args$ncd.max; if (is.null(ncd.max)) ncd.max <- 100
  
  if (verbose) print('WG.fwmu.day.precip')
  # Single function for just precipitation
  if (is.null(x)) {
    if (verbose) print('Use sample data from esd')
    bjornholt <- NULL
    data("bjornholt",envir=environment())
    x <- bjornholt
    rm('bjornholt')
  } else if (verbose) print(paste('Use data provided:',loc(x)))
  
  ## Estimate climatology for mean seasonal cycle in total precipitation. Use this information
  ## as a guide for which months to add wet days to ensure correct wet-day frequency fw - 
  ## this is important for locations with a rainy season
  #if (verbose) print('Get the seasonal cycle')
  fw.ac <- aggregate(x,month,FUN='wetfreq',threshold=1,na.rm=TRUE)
  
  ## Also find the climatology for the wet-day mean precipitation mu
  mu.ac <- aggregate(x,month,FUN='wetmean',threshold=1,na.rm=TRUE)
  ## If mu is not defined because there is no wet days, set the weights to zero 
  mu.ac[!is.finite(mu.ac)] <- 0
  if (verbose) print(cbind(fw.ac,mu.ac))
  # use fw to estimate the number of rainy days per year:
  x.fw <-  annual(x,'wetfreq',threshold=threshold,nmin=30)
  
  # Use predicted mu to generate exponentially distributed data:
  x.mu <- annual(x,'wetmean',threshold=threshold,nmin=30)
  
  # Number of consecutive wet/dry days
  if (verbose) print('Get the spell duration statistics')
  ncd <- subset(spell(x,threshold=threshold),is=1)
  good <- !is.na(index(ncd))
  ncd <- ncd[good]
  ncd[ncd > ncd.max] <- NA
  ## Annual mean number of consecutive wet days
  amncwd <- subset(annual(ncd, nmin=1), is=1)
  if (sum(is.finite(amncwd))==0) browser()
  ismissing <- !is.finite(amncwd)
  ## If there are missing data, use the mean value
  if (sum(ismissing)>0) amncwd[ismissing] <- mean(amncwd,na.rm=TRUE)
  # extract the time interval between the start of each dry spell
  dt1 <- diff(julian(as.Date(index(ncd[is.finite(ncd[,1]),1]))))
  
  if (plot) {
    dev.new()
    par(mfrow=c(2,2),cex.main=0.7)
    ## Exponential distribution
    x.wet <- coredata(x); x.wet <- x.wet[x.wet > threshold]
    f.x <- dexp(0:ceiling(max(x.wet)),rate=1/wetmean(x))
    hist(x.wet,freq=FALSE,col="grey",xlab="mm",
         main="The precipitation depth",
         sub="Test: Fitted exponential distribution")
    lines(0:ceiling(max(x.wet)),f.x,lwd=5,col="red")
    grid()
    
    ## Timing between each precipitation event
    f.k <- dgeom(0:max(dt1), prob=1/(mean(dt1)+1))
    hist(dt1,freq=FALSE,col="grey",xlab="days",
         main="Time between the start of each rainfall",
         sub="Test: Fitted geometric distribution")
    lines(0:max(dt1),f.k,lwd=5,col="red")
    grid()
  }
  
  # Wet-day mean: from DS or from observations
  #if (verbose) print('wet-day mean')
  if (is.null(mu))
    mu <- zoo(FTscramble(x.mu),order.by=index(x.mu)) else
      ## mu is introduced as a change factor
      if (length(mu)==1) {
        mu <- mu + zoo(FTscramble(x.mu),order.by=index(x.mu))
      }
  rm('x.mu')
  coredata(mu)[mu<=0] <- NA
  coredata(mu)[!is.finite(coredata(mu))] <- mean(mu,na.rm=TRUE)
  
  # Wet-day frequency: from DS or from observations
  #if (verbose) print('wet-day frequency')
  if (is.null(fw)) 
    fw <- zoo(FTscramble(x.fw),order.by=index(x.fw)) else
      ## fw is introduced as a change factor
      if (length(fw)==1) {
        fw <- fw + zoo(FTscramble(x.mu),order.by=index(x.mu))
      }
  rm('x.fw')
  coredata(fw)[fw==0] <- mean(fw,na.rm=TRUE)
  
  if (plot) {
    ## Number of events per year
    hist(coredata(ncd),breaks=seq(0,40,by=2),freq=FALSE,col="grey",
         main="Duration of wet spells",xlab="days",
         sub="Test: Ftted geometric distribution")
    lines(seq(0,40,by=1),dgeom(seq(0,40,by=1),prob=1/mean(coredata(amncwd))),
          col="red",lwd=3)
    grid()
  }
  
  ## Time axis for projection:
  if (verbose) print('Time axis for projection')
  if (is.null(t)) {
    nxy <- range(year(mu))
    t <- seq(as.Date(paste0(nxy[1],'-01-01')),as.Date(paste0(nxy[2],'-12-31')),by=1)
    ## Number of years
    ny <- length(rownames(table(year(mu)))) 
    if (verbose) print(range(t))
  }
  ## Number of days
  nd <- length(t)
  yrs <- as.numeric(rownames(table(year(t))))
  
  # Estimate the annual number of rainy days:
  if (verbose) print('Number of wet days each year:')
  anwd <- round( ( julian(as.Date(paste(year(fw),'12-31',sep='-'))) -
                     julian(as.Date(paste(year(fw),'01-01',sep='-'))) + 1) *
                   coredata(fw) )
  
  # Error bars for mu: var = mu**2 for exponential distribution:
  mu.err <- mu/sqrt(anwd - 1)
  
  # set up a record with no rain:
  z <- zoo(rep(NA,nd),order.by=t)
  
  # add rain events:
  if (verbose) print(paste('loop over year:',1,'-',ny,'number of days=',nd,length(z),
                           'length(mu)=',length(mu),'length(fw)=',length(fw),length(anwd)))
  for ( i in 1:ny ) {
    #if (verbose) print(i)
    ii <- is.element(year(t),yrs[i])
    ndaysyr <- sum(ii)
    ## Duration of wet events
    if (i <= length(amncwd)) ncwd <- rgeom(ndaysyr,prob=1/amncwd[i]) + 1 else
      ncwd <- rgeom(ndaysyr,prob=1/mean(amncwd,na.rm=TRUE)) + 1
    
    ## White noise to introduce stochastic weather and find most suitable times of 
    ## the year with stochastic influence
    fw.ac.wn <- approx(1:12,fw.ac,xout = seq(1,12,length=ndaysyr))$y + 0.01*rnorm(ndaysyr)
    ## Use ij as index for timing wet events
    ij <- order(fw.ac.wn,decreasing=TRUE)
    
    ## Go through each event and place according to climatology and stochastic weather
    #if (verbose) print('Select wet and dry days')
    dry <- c(); wet <- c(); nes <- 1
    while ( (length(wet) < anwd[i]) & (nes <= ndaysyr) ) {
      ## Check whether the selected days are available: start with the first julian day in the year
      idy1 <- 1
      ## TRUE if found available sequence of wet days
      d.available <- FALSE
      ## We search the days in the year for sequences that include the wet spell duration 
      ## padded by dry days. iseq is an index of sequences and ncwd is the number of consecutive
      ## wet days.
      while( (!d.available) & (idy1 <= ndaysyr) ) { 
        ## sequence of days: wet spell padded by dry days
        if (is.finite(ncwd[1])) iseq <- ij[idy1] + seq(-1,ncwd[1]+1,by=1) else
          iseq <- ij[idy1] + seq(-1,2,by=1)
        ## Test if all elements in iseq also are in ij
        if (length(intersect(iseq,ij))==length(iseq)) d.available <- TRUE else
          idy1 <- idy1 + 1
      }
      ## If no available sequence of days was found, then pick just random individual days available 
      ## from the pool of remaining days
      if (!d.available) {
        iseq <- ij[seq(1,length(ncwd[1])+2,by=1)]
      }
      ## Check that dry and wet contain valid Julian days from ij also, if there are elements
      ## out of sample, then add new random elements from ij. nseq is number of sequences
      nseq <- length(iseq)
      iseq <- intersect(iseq,ij)
      dseq <- setdiff(iseq,ij)
      ## Test if all indices stored in iseq also are listed in ij
      diffseq <- nseq - length(iseq)
      if (diffseq > 0) iseq <- c(iseq,dseq[sort(rnorm(length(dseq)))][1:diffseq])
      ## Once a suitable sequence of days have been located, use it to define wet spell padded
      ## with dry days
      
      dry <- c( dry, iseq[c(1,length(iseq))] )
      wet <- c( wet, iseq[2:(length(iseq)-1)] )
      ## Remove duplicated (overlapping) days in wet and dry
      wet <- wet[!duplicated(wet)]
      dry <- dry[!duplicated(dry)]
      
      ## Remove duplicates - for some reason, there are some of them...
      ndupl <- sum(duplicated(c(dry,wet)))
      #if (verbose) print(paste(ndupl,'duplicated in wet and dry'))
      
      ## Remove used indices and used wet-spell duration
      ij <- ij[!is.element(ij,intersect(c(dry,wet),ij))]; ncwd <- ncwd[-1]
      
      ## Increment number of events
      nes <- nes + 1
    }
    
    ## Finish sorting all the days into wet and dry  
    dry <- sort(c(dry,ij)); wet <- sort(wet)
    ## This should not happen, but ...
    if (sum(duplicated(dry)>0)| sum(duplicated(wet))>0) browser("WG.fwmu: should not happen")
    #dry <- dry[!duplicated(dry)]; wet <- wet[!duplicated(wet)]
    ## deal with cases where days are classified as both dry and wet
    inboth <- intersect(wet,dry)
    if (length(inboth)>0) browser("WG.fwmu: same days in in bot wet and dry")
    #dry <- dry[!is.element(dry,inboth)]
    
    ## Quality control: If there are too few wet days, add random wet days to or
    ## if there are too many, then replace the excess with dry days
    ## anwd[i] is the number of wet days for year i based on fw and the number of days per year
    nwdd <- length(wet) - anwd[i]
    #if (verbose) print(paste(nwdd,'too many wet days'))
    if (nwdd < 0) {
      swap <- order(rnorm(length(dry)))[1:abs(nwdd)]
      wet <- sort(c(wet,dry[swap])); dry <- sort(dry[-swap]) 
    } else if (nwdd > 0) {
      swap <- order(rnorm(length(wet)))[1:nwdd]
      dry <- sort(c(dry,wet[swap])); wet <- sort(wet[-swap]) 
    }
    
    if (i > length(mu)) browser("WG.fwmu: should not happen")
    ## The wet-day mean precipitation amount
    if (!is.finite(mu[i])) mu[i] <- mean(mu,na.rm=TRUE)
    
    ## The daily amounts for wet days - first sort the data according to magnitude
    ## then shuffle them according to a mix of chance and mu climatology
    ## add rain to the appropriate year:
    rain <- rexp(ndaysyr,rate=1/coredata(mu[i]))
    ## amounts less then thresholds have been set to dry days - reset these by repeat throwing the dice
    iybt <- rain < threshold
    while (sum(iybt)>0) {
      rain[iybt] <- rexp(sum(iybt),rate=1/coredata(mu[i]))
      iybt <- rain < threshold
    }
    
    if (alpha.scaling) {
      ## REB 2024-05-13
      ## Scale the amounts according to return-period according to 
      ## DOI:https://doi.org/10.1088/1748-9326/ab2bb2 see day2IDF
      ## tau - return-interval in years 
      if (verbose & (i==1)) print('Scale by alpha according to return-interval')
      ## Estimate the return period tau =  1/Pr(X>x)
      if (is.finite(mu[i])) tau <- 1/( 1-pexp(rain,rate=1/coredata(mu[i])) ) else
        tau <- 1/( 1- pexp(rain,rate=1/coredata(mean(mu,na.rm=TRUE))) )
      ## Take into account the fraction of wet days and express return interval in years
      if (fw[i] > 0) tau <- tau/(365.25*coredata(fw[i])) else
        tau <- tau/(365.25*coredata(mean(fw,na.rm=TRUE)))
      alphas <- alpha[1] + alpha[2]*as.numeric(log(tau))
      alphas[alphas < 1] <- 1
      rain <- rain * alphas
      #if (verbose) {print(summary(tau)); print(summary(alphas)); print(length(rain))}
    }
    
    ## Set zero rainfall on dry days
    dry <- dry[dry <= ndaysyr]
    rain[dry] <- 0
    ## ensure that the rain has correct mu climatology 
    #if (verbose) print('Apply daily climatological scaling for mu')
    rain <- rain*approx(1:12,mu.ac,xout = seq(1,12,length=length(rain)))$y
    ## ensure that wet-day mean mu in rain matches mu[i]
    if (is.finite(mu[i])) rain <- rain * coredata(mu[i])/wetmean(rain) else 
      rain <- rain * coredata(mean(mu,na.rm=TRUE))/wetmean(rain)
    
    if (verbose) print(paste('i=',i,'year=',yrs[i],'days=',length(rain),'=',ndaysyr,
                             'tot rain',round(sum(rain,na.rm=TRUE)),
                             'mm/year, #wet days=',length(wet),'=','n*fw[i]=',anwd[i],
                             'mu[i]=',round(mu[i],1),'=',round(wetmean(rain),1),
                             ' #events=',nes,' days done:',sum(is.finite(z)),
                             ' [',min((1:nd)[ii]),',',max((1:nd)[ii]),']'))
    z[ii] <- rain
  }
  z <- zoo(z,order.by=t)
  if (plot) {
    zx <- coredata(matchdate(x,z))
    zz <- coredata(matchdate(z,x))
    if (length(zx)!=length(zz)) browser("different data lengths")
    plot(sort(zx,decreasing=TRUE),sort(zz,decreasing=TRUE),pch=19,
         xlim=c(0,max(zx,zz)),ylim=c(0,max(zx,zz)),
         main='QQ-plot: Wet-day amounts (mm)',xlab='Observed',ylab='WG')
    grid()
    maxzy <- max(z,y,na.rm=TRUE)
    lines(c(0,maxzy),c(0,maxzy),lty=2,col='red')
  }
  class(z) <- class(x)
  z <- attrcp(x,z)
  attr(z,'original_fw') <- fw
  attr(z,'original_mu') <- mu
  attr(z,'alpha_scaling') <- alpha.scaling
  attr(z,'mu_error') <- mu.err
  attr(z,'aspect') <- paste(attr(z,'aspect'),'weather_generator',sep=', ')
  attr(z,'history') <- history.stamp(x)
  #if (verbose) print('EXIT WG')
  return(z)
}


#' This function tests the WG for precipitation:
#' Quantile-quantile plots of wet-day amounts
#' Number of wet days
#' @exportS3Method
#' @export test.WG.fwmu.day.precip
test.WG.fwmu.day.precip <- function(x=NULL,verbose=FALSE) {
  if (verbose) print('test.WG.fwmu.day.precip')
  if (is.null(x)) {data('bjornholt'); x <- bjornholt; rm('bjornholt')}
  print(paste('test.WG.fwmu.day.precip for',loc(x)))
  z <- WG(x,verbose=TRUE)
  z0 <- WG(x,alpha.scaling=FALSE)
  ## sort magnitudes to plot quantile-quantile plots
  if (verbose) print('Sort precipitation magnitudes')
  xw <- sort(coredata(x)[x > 1])
  zw <- sort(coredata(z)[z > 1])
  zw0 <- sort(coredata(z0)[z0 > 1])
  ## There may be different number of wet days - pad the shortest series with 0s
  nx <- length(xw)
  nz <- length(zw)
  nz0 <- length(zw0)
  print(paste(nx,'observed wet days and',nz,'simulated wet days - without scaling, there were ',nz0,'days'))
  print('Obs:');print(summary(x)); print('WG:');print(summary(z))
  if (nx > nz) zw <- c(rep(0,nx-nz),zw) else 
    if (nz > nx) xw <- c(rep(0,nz-nx),xw)
  nx <- length(xw); nz <- length(zw)
  if (nx > nz0) zw0 <- c(rep(0,nx-nz0),zw0) else 
    if (nz0 > nx) zw0 <- zw0[1:nx] 
  nz0 <- length(zw0)
  
  xylim <- c(0,max(c(xw,zw)))
  par(mfcol=c(2,2))
  plot(xw,zw,ylim=xylim,xlim=xylim,xlab='Observed amount (mm/day)',
       ylab='WG amount (mm/day)',main=paste(loc(x),'wet-day amounts'))
  points(xw,zw0,col=rgb(0.5,0.5,0.5,0.5),cex=0.5)
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
  
  segY <- segY*sd(y,na.rm=TRUE) + mean(y,na.rm=TRUE)
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
