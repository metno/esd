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
#' functions procude stochastic time series with similar characteristics as the
#' station series provided (if none if provided, it will use either ferder or
#' bjornholt provided by the esd-package). Here characteristics means similar
#' mean value, standard deviation, and spectral properties. \code{FTscramble}
#' takes the Fourier components (doing a Fourier Transform - FT) of a series
#' and reassigns random phase to each frequency and then returns a new series
#' through an inverse FT. The FT scrambling is used for temperature, but not
#' for precipitation that is non-Gaussian and involves sporadic events with
#' rain. For precipitation, a different approach is used, taking the wet-day
#' frequency of each year and using the wet-day mean and ranomly generated
#' exponentially distributed numbers to provide similar aggregated annual
#' statistics as the station or predicted though downscaling. The precipitation
#' WG can also take into account the number of consequtive number-of-dry-days
#' statistics, using either a Poisson or a gemoetric distribution.
#' 
#' The weather generater produces a series with similar length as the provided
#' sample data, but with shifted dates according to specified scenarios for
#' annual mean mean/standard deviation/wet-day mean/wet-day frequency.
#' 
#' \code{WG.FT.day.t2m} generates daily temperature from seasonal means and
#' standard deviations. It is given a sample station series, and uses
#' \code{FTscramble} to generate a series with random phase but similar (or
#' predicted - in the future) spectral characteristics. It then uses a quantile
#' transform to prescribe predicted mean and standard deviation, assuming the
#' distributions are normal. The temperal structure (power spectrum) is
#' therefore similar as the sample provided.
#' 
#' \code{WG.fw.day.precip} uses the annual wet-day mean and the wet-day
#' frequency as input, and takes a sample station of daily values to
#' stochastically simulate number consequtive wet days based on its annual mean
#' number. If not specified, it is taken from the sample data after being phase
#' scrambeled (\code{FTscramble}) The number of wet-days per year is estimated
#' from the wed-day frequency, it too taken to be phase scrambled estimates
#' from the sample data unless specifically specified. The daily amount is
#' taken from stochastic values generated with \code{\link{rexp}}. The number
#' of consequtive wet days can be approximated by a geometric distribution
#' (\code{\link{rgeom}}), and the annual mean number was estimated from the
#' sample series.
#' 
#' @aliases WG WG.station WG.fw.day.precip WG.FT.day.t2m
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
#' @param \dots additional arguments
#' @author R.E. Benestad
#' @keywords manip
#' @examples
#' 
#' data(ferder)
#' t2m <- WG(ferder)
#' data(bjornholt)
#' pr <- WG(bjornholt)
#' 
#' @export WG
WG <- function(x,...) UseMethod("WG")

#' @exportS3Method
#' @export WG.station
WG.station <- function(x,...,option='default') {
  if (inherits(x,'day')) {
    if (length(varid(x))==1) {
      if (varid(x)=='t2m') y <- WG.FT.day.t2m(x,...) else
        if (varid(x)=='precip') y <- WG.fw.day.precip(x,...)
    }
  }
  return(y)
}

#' @exportS3Method
#' @export WG.FT.day.t2m
WG.FT.day.t2m <- function(x=NULL,...,amean=NULL,asd=NULL,t=NULL,ip=1:4,
                          select=NULL,lon=c(-20,20),lat=c(-20,20),
                          plot=FALSE,biascorrect=TRUE,verbose=FALSE) {
  if (verbose) print('WG.FT.day.t2m')
  ## Single function for just temperature.
  ## The arguments mean and sd are time series predicted through ESD or
  ## adopted from a zoo or station object (x). 
  if (is.null(x)) {
    ## If no stations objects is given, use default 
    if (verbose) print("use default: Ferder, Norway")
    ferder <- NULL
    data("ferder",envir=environment())
    x <- ferder
    rm('ferder')
  }
  
  ## Different options for annual mean temperature. Default - estimate from the station
  if (is.null(amean)) amean <- annual(x) else
    ## If NA, then compute using DSensemble
    if (is.na(amean)) {
      if (verbose) print('Estimate mean change')
      T2M <- retrieve('~/data/ERAINT/ERAINT_t2m_mon.nc',
                      lon=lon(x) + lon,lat=lat(x) + lat)
      ztm <- DSensemble.t2m(x,predictor=T2M,biascorrect=biascorrect,
                            plot=plot,lon=lon,lat=lat,ip=ip,
                            select=select,verbose=verbose)
      amean <- zoo(rowMeans(ztm,na.rm=TRUE) - mean(ztm,na.rm=TRUE),
                   order.by=index(ztm))
    } else if (inherits(amean,'dsensemble'))
      ## Or use prescribed projections
      amean <- rowMeans(amean,na.rm=TRUE) - mean(amean,na.rm=TRUE)
    if(verbose) print(paste('mean(amean)=',mean(amean)))
    
    ## Also select annual standard deviations estimated from daly anomalies -
    ## repeat the same procedure as for the mean.
    if (is.null(asd)) asd <- annual(anomaly(x,verbose=verbose),FUN='sd') else
      if (is.na(asd)) {
        if (verbose) print('Estimate standard deviation change')
        SLP <- retrieve('~/data/ERAINT/ERAINT_slp_mon.nc',
                        lon=lon(x) + lon,lat=lat(x) + lat)
        coredata(SLP) <- 100*coredata(SLP)  # The CMIP5 units are in Pa!
        if (plot) dev.new()
        #    zts <- DSensemble.t2m(x,predictor=SLP,biascorrect=biascorrect,
        #                          FUN='sd',plot=plot,lon=lon,lat=lat,ip=ip,
        #                          path='data/CMIP5.mslp/',pattern='psl_Amon_ens',
        #                          select=select,verbose=verbose)
        zts <- DSensemble.t2m(x,predictor=T2M,biascorrect=biascorrect,
                              FUN='sd',plot=plot,lon=lon,lat=lat,ip=ip,
                              FUNX='sd',select=select,verbose=verbose)
        asd <- zoo(rowMeans(zts,na.rm=TRUE) - mean(zts,na.rm=TRUE),
                   order.by=index(zts))
      } else if (inherits(asd,'dsensemble'))
        asd <- rowMeans(asd,na.rm=TRUE) - mean(asd,na.rm=TRUE) 
    
    ## Get the daily anomalies and the climatology
    xa <- anomaly(x); clim <- x - xa
    
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
    attr(z,'mean') <- ym
    attr(z,'sd') <- ys
    attr(z,'aspect') <- paste(attr(z,'aspect'),'weather_generator',sep=', ')
    attr(z,'history') <- history.stamp(x)
    return(z)
}

## Fuure considerations -lso allow for estimating the AR(1) coefficient of the Hurst coefficient?
## Fractional Gaussian noise...?

## --- Precipitation  
#' @exportS3Method
#' @export WG.fw.day.precip
WG.fw.day.precip <- function(x=NULL,...,mu=NULL,fw=NULL,
                             ncwd=NULL,ndbr=NULL,t=NULL,
                             threshold=1,select=NULL,
                             ip=1:6,lon=c(-10,10),lat=c(-10,10),
                             plot=FALSE,biascorrect=TRUE,
                             verbose=TRUE) {
  
  if (verbose) print('WG.fw.day.precip')
  # Single function for just precipitation
  if (is.null(x)) {
    bjornholt <- NULL
    data("bjornholt",envir=environment())
    x <- bjornholt
    rm('bjornholt')
  }
  
  # use fw to estimate the number of rainy days per year:
  x.fw <-  as.annual(x,'wetfreq',threshold=threshold)
  
  # Use predicted mu to generate exponentially distributed data:
  x.mu <- as.annual(x,'exceedance',threshold=threshold)
  
  # according to a geometric (default) or Poisson distribution
  ncdd.cwd <- spell(x,threshold=threshold)
  x.nd <- subset(annual(ncdd.cwd),is=1)
  # extract the time interval between the start of each dry spell
  ndbr <- diff(julian(as.Date(index(ncdd.cwd[is.finite(ncdd.cwd[,1]),1]))))
  if (plot) {
    dev.new()
    f.k <- dgeom(0:max(ndbr), prob=1/(mean(ndbr)+1))
    hist(ndbr,freq=FALSE,col="grey",xlab="days",
         main="The time between the start of each precipitation event",
         sub="Test: Red curve is the fitted geometric distribution")
    lines(0:max(ndbr),f.k,lwd=5,col="red")
  }
  
  # Aggregate the number of days between start of each rain event
  # to annual mean
  ndbram <- annual(zoo(x=ndbr,order.by=index(ncdd.cwd)[-1]))
  
  # Estimate number of wet events each year:
  wet <-   subset(ncdd.cwd,is=1)
  nawe <- aggregate(wet,by=year(wet),FUN="nv")
  attributes(nawe) <- NULL
  if (verbose) print(coredata(nawe))
  
  if (plot) {
    dev.new()
    hist(coredata(nawe),breaks=seq(0,100,by=5),freq=FALSE,col="grey",
         main="Number of wet events per year",xlab="days",
         sub="Test: Red curve is the fitted Poisson distribution")
    lines(seq(0,100,by=1),dpois(seq(0,100,by=1),lambda=mean(coredata(nawe))),
          col="red",lwd=3)
  }
  
  # Wet-day mean: from DS or from observations
  if (verbose) print('wet-day mean')
  if (is.null(mu))
    mu <- zoo(FTscramble(x.mu),order.by=index(x.mu)) else
      if (is.na(mu)) {
        if (verbose) print('Estimate mean change')
        PRE <- retrieve('~/data/ERAINT/ERAINT_precip_mon.nc',
                        lon=lon(x) + lon,lat=lat(x) + lat)
        zmu <- DSensemble.precip(x,predictor=PRE,biascorrect=biascorrect,
                                 plot=plot,lon=lon,lat=lat,ip=ip,
                                 treshold=threshold,
                                 select=select,verbose=verbose)
        mu <- rowMeans(zmu,na.rm=TRUE) - mean(zmu,na.rm=TRUE)
      } else if (inherits(mu,'dsensemble'))
        mu <- rowMeans(mu,na.rm=TRUE) - mean(mu,na.rm=TRUE) 
  
  # Wet-day frequency: from DS or from observations
  if (verbose) print('wet-day frequency')
  if (is.null(fw))
    fw <- zoo(FTscramble(x.fw),order.by=index(x.fw)) else
      if (is.na(fw)) {
        SLP <- retrieve('~/data/ERAINT/ERAINT_slp_mon.nc',
                        lon=lon(x) + lon,lat=lat(x) + lat)
        coredata(SLP) <- 100*coredata(SLP)  # The CMIP5 units are in Pa!
        if (plot) dev.new()
        zfw <- DSensemble.precip(x,predictor=SLP,biascorrect=biascorrect,
                                 FUN='wetfreq',threshold=threshold,
                                 plot=plot,lon=lon,lat=lat,ip=ip,
                                 path='data/CMIP5.mslp/',pattern='psl_Amon_ens',
                                 select=select,verbose=verbose)
        fw <- rowMeans(zfw,na.rm=TRUE) - mean(zfw,na.rm=TRUE)
      } else if (inherits(fw,'dsensemble'))
        fw <- rowMeans(fw,na.rm=TRUE) - mean(fw,na.rm=TRUE)
  
  # Number of consequtive wet days: from DS or from observations
  if (verbose) print('random annual mean number of n_cwd:')
  rnd <- rnorm(length(mu),mean=mean(coredata(x.nd),na.rm=TRUE),
               sd=sd(coredata(x.nd),na.rm=TRUE))
  rnd[rnd < 1] <- 1;
  if (verbose) print(rnd)
  prob <- 1/rnd
  if (verbose) print('the annual mean probability of successive wet days: prob')
  if (verbose) print(prob)
  if (verbose) print('the annual mean number of consecutive wet days: ncwd')
  if (is.null(ncwd)) ncwd <- rgeom(length(mu),prob=prob)+1 else
    if (is.na(ncwd)) {
      if (verbose) print('estimate from ERAINT')
      SLP <- retrieve('~/data/ERAINT/ERAINT_slp_mon.nc',
                      lon=lon(x) + lon,lat=lat(x) + lat)
      coredata(SLP) <- 100*coredata(SLP)  # The CMIP5 units are in Pa!
      if (plot) dev.new()
      y.ncwd <- annual(subset(spell(x,threshold=threshold),is=1))
      znd <- DSensemble.precip(y.ncwd,predictor=SLP,biascorrect=biascorrect,
                               plot=plot,lon=lon,lat=lat,ip=ip,
                               path='data/CMIP5.mslp/',pattern='psl_Amon_ens',
                               select=select,verbose=verbose)
      ncwd <- rowMeans(znd,na.rm=TRUE) - mean(znd,na.rm=TRUE)
    } else {
      if (verbose) print("inherits(ncwd,'dsensemble')")
      if (inherits(ncwd,'dsensemble'))
        ncwd <- rowMeans(ncwd,na.rm=TRUE) - mean(ncwd,na.rm=TRUE)
      if (verbose) print(ncwd)
    }
  
  # Time axsis for projection:
  if (verbose) print('Time axis for projection')
  if (is.null(t)) {
    ly <- max(year(mu))
    ny <- length(rownames(table(year(mu)))) 
    interval <- c(ly-ny+1,ny)
    t <- index(x)
    t <- t - julian(as.Date(t[1])) +
      julian(as.Date(paste(interval[1],month(x)[1],day(x)[1],sep='-')))
    if (verbose) print(interval)
  }
  n <- length(t)
  yrs <- rownames(table(year(t)))
  
  # One alternative: qq-transform: precip(exp1 -> exp2)
  #  pr.x.wet <- qexp(pexp(q=round(365.25*coredata(x.fw)),
  #                        rate=1/coredata(x.mu),na.rm=TRUE)),
  #                   rate=1/coredata(mu))
  
  # Estimate the number of rainy days for each year
  if (verbose) print('Number of wet days each year:')
  nwet <- round( ( julian(as.Date(paste(year(fw),'12-31',sep='-'))) -
                     julian(as.Date(paste(year(fw),'01-01',sep='-'))) + 1) *
                   coredata(fw) )
  #print(rbind(nwet,coredata(mu)))
  
  # Errorbars for mu: var = mu**2 for exponential distrib:
  mu.err <- mu/sqrt(nwet -1)
  
  # set up a record with no rain:
  z <- zoo(rep(0,length(t)),order.by=t)
  
  #print(c(length(annual(x)),length(x.fw),length(x.mu),length(fw),length(mu)))
  #print(start(annual(x)));print(end(annual(x)))
  #print(start(x.mu));print(end(x.mu))
  #print(start(x.fw));print(end(x.fw))
  #print(start(mu));print(end(mu))
  #print(start(fw));print(end(fw))
  
  #browser()
  
  # add rain events:
  if (verbose) print(paste('loop over year:',1,'-',ny))
  for ( i in 1:ny ) {
    # Simulate precipitation amount: reduce to one decimal point and add
    # threshol to mimic the original data...
    
    if (!is.finite(mu[i])) mu[i] <- mean(mu,na.rm=TRUE)
    if (!is.finite(ndbram[i])) ndbram[i] <- mean(ndbram,na.rm=TRUE)
    if (!is.finite(nwet[i])) nwet[i] <- mean(nwet,na.rm=TRUE)
    
    y <- round(rexp(nwet[i],1/coredata(mu[i])),1) + threshold # Check!
    
    # Simulate the start of each rain event: 
    t0 <- cumsum(rgeom(max(coredata(nawe),na.rm=TRUE),
                       prob=1/(coredata(ndbram[i]))))
    #simulate the duration of wet events:
    #str(t0); print(prob)
    nwd <- rgeom(length(t0),prob=prob[i])+1 
    nwd <- nwd[cumsum(nwd <= nwet[i])]
    t0 <- t0[1:length(nwd)]
    
    # add rain to the appropriate year:
    ii <- is.element(year(t),yrs[i])
    rain <- rep(0,sum(ii)); iii <- 0
    
    #browser()
    
    # simulate the rain-eve-t occurrances which start at t0 and vary in
    # duration: nwd
    for (iv in 1:length(nwd)) {
      iii <- max(iii) + (1:nwd[iv])
      rain[t0[iv] + 0:(nwd[iv]-1)] <-  y[iii]
    }
    
    if (verbose) print(paste(i,yrs[i],' total rain=',sum(rain,na.rm=TRUE),
                             ' #wet days=',sum(nwd),
                             'nwet[i]=',nwet[i],' #events=',length(nwd)))
    z[ii] <- rain
  }
  z <- attrcp(x,z)
  attr(z,'original fw') <- fw
  attr(z,'ncc') <- nwd
  attr(z,'original mu') <- mu
  attr(z,'mu.err') <- mu.err
  attr(z,'aspect') <- paste(attr(z,'aspect'),'weather_generator',sep=', ')
  attr(z,'history') <- history.stamp(x)
  return(z)
}


# This weather generator assumes that the past covariate structure between
# temperature and precipitation is constant and doesn't change in the future.
# Moreover, the method also assumes that the spell-statistics will stay the same.
#' @exportS3Method
#' @export WG.pca.day.t2m.precip
WG.pca.day.t2m.precip <- function(x=NULL,...,precip=NULL,threshold=1,select=NULL,
                                  wetfreq.pred=FALSE,spell.stats=FALSE,
                                  verbose=FALSE) {
  if(verbose) print("WG.pca.day.t2m.precip")			
  t2m <- x
  if (is.null(t2m)) {
    ferder <- NULL
    data("ferder",envir=environment())
    t2m <- ferder
    rm('ferder')
  }
  if (is.null(precip)) {
    bjornholt <- NULL
    data("bjornholt",envir=environment())
    pr <- bjornholt
    rm('bjornholt')
  }
  
  lon <- lon(t2m) + c(-10,10)
  lat <- lat(t2m) + c(-10,10)
  
  SLP <- retrieve('~/data/ERAINT/ERAINT_slp_mon.nc',lon=lon,lat=lat)
  coredata(SLP) <- 100*coredata(SLP)  # The CMIP5 units are in Pa!
  T2M <- retrieve('~/data/ERAINT/ERAINT_t2m_mon.nc',lon=lon,lat=lat)
  PRE <- retrieve('~/data/ERAINT/ERAINT_pr_mon.nc',lon=lon,lat=lat)
  
  ztm = DSensemble.t2m(t2m,predictor=T2M,biascorrect=TRUE,plot=FALSE,
                       lon=c(-10,10),lat=c(-10,10),select=select)
  zts = DSensemble.t2m(t2m,FUN='sd',predictor=SLP,biascorrect=TRUE,plot=FALSE,
                       path='data/CMIP5.mslp/',pattern='psl_Amon_ens',
                       lon=c(-10,10),lat=c(-10,10),select=select)
  
  zpm = DSensemble.precip(pr,FUN='wetmean',predictor=PRE,biascorrect=TRUE,
                          plot=FALSE,select=select)
  if (wetfreq.pred) {
    zpf = DSensemble.precip(pr,FUN='wetfreq',predictor=SLP,biascorrect=TRUE,plot=FALSE,
                            path='data/CMIP5.mslp/',pattern='psl_Amon_ens',select=select)
  }
  # select an equal number of years at the end of the downscaled results if none is specified
  if (is.null(interval)) {
    interval <- c(max(year(ztm))-ny+1,max(year(ztm)))
  }
  ztm  <- subset(ztm,it=interval)
  zts  <- subset(zts,it=interval)
  zpm  <- subset(zpm,it=interval)
  
  # Generate the daily variability from the observations:
  t2ma <- anomaly(t2m); t2mc <- t2m - t2ma
  xy <- merge(zoo(t2m),zoo(pr),zoo(t2mc),all=FALSE)
  years <- as.numeric(rownames(table(year(xy))))
  ny <- length(years)
  
  wet <- (xy[,2] > threshold) & (is.finite(xy[,1])) & (is.finite(xy[,2]))
  dry <- (xy[,2] <= threshold) & (is.finite(xy[,1])) & (is.finite(xy[,2]))
  if (spell.stats) L <- spell(pr,threshold=threshold)
  
  # Rainy days:
  # Find the co-variate structures in daily temperature and precipitation
  X <- coredata(xy[wet,1:2])
  X[,2] <- log(X[,2])
  qqnorm(X[,2])
  pca <- svd(X)
  # scramble the principal components:
  z1 <- FTscramble(pca$u[,1])
  z2 <- FTscramble(pca$u[,2])
  pca$u <- cbind(z1,z2)
  X <- pca$u %*% diag(pca$d) %*% t(pca$v)
  X[,2] <- exp(X[,2])
  
  # All days: set up date string
  t <- julian(as.Date(paste(year(xy)-year(xy)[1]+interval[1],
                            month(xy),day(xy),sep='-')))
  x <- zoo(coredata(xy[,1]),order.by=t)
  attributes(x) <- NULL
  
  # initialise temperature and precip - for t2m, add observed daily climatology
  n <- length(x)
  pr.x <- rep(0,n) 
  
  # shoehorn the rainy day temperature and precipitation into series
  t2m.x[wet] <-  X[,1]
  
  # generate time series for predicted mean and sd with same length as the scrambled
  # daily series: use the seasonal estimates
  tensm <- rowMeans(ztm,na.rm=TRUE) - mean(ztm,na.rm=TRUE)
  tenss <- rowMeans(zts,na.rm=TRUE) - mean(zts,na.rm=TRUE)
  penss <- rowMeans(zpm,na.rm=TRUE) - mean(zpm,na.rm=TRUE)
  
  ytm <- approx(julian(as.Date(index(ztm))),tensm,xout=julian(as.Date(t)),rule=2)$y
  yts <- approx(julian(as.Date(index(zts))),tensm,xout=julian(as.Date(t)),rule=2)$y
  ypm <- approx(julian(as.Date(index(ztm))),tensm,xout=julian(as.Date(t)),rule=2)$y
  
  #qq-transform: temp(N1 -> N2) - year by year or for a given interval?
  t2m.x <- qnorm(pnorm(q=t2m.x,mean=mean(t2m.x,na.rm=TRUE),sd=sd(t2m.x,na.rm=TRUE)),
                 mean=ytm,sd=yts)
  
  #qq-transform: precip(exp1 -> exp2)
  pr.x.wet <- qexp(pexp(q=X[,2],rate=1/mean(xy[,2],na.rm=TRUE)),
                   rate=1/ypm)
  
  #empirical adjustment to precipitation?
  
  # shoehorn the rainy day temperature and precipitation into series
  pr.x[wet] <- pr.x.wet 
  t2m.x <- zoo(t2m.x,order.by=t)
  t2m.x <- attrcp(t2m)
  attr(t2m.x,'method') <- 'DShydro'
  attr(t2m.x,'type') <- 'result: ESD + WG'
  attr(t2m.x,'activity') <- 'CMIP5'
  attr(t2m.x,'experiment') <- 'RCP 4.5'
  pr.x <- zoo(pr.x,order.by=t)
  pr.x <- attrcp(pr)
  attr(pr.x,'method') <- 'DShydro'
  attr(pr.x,'type') <- 'result: ESD + WG'
  attr(pr.x,'activity') <- 'CMIP5'
  attr(pr.x,'experiment') <- 'RCP 4.5'
  
  y <- list(t2m=t2m.x,pr=pr.x)
  attr(y,'history') <- history.stamp(t2m)
  class(y) <- class(t2m)
  return(y)
}

#' @export
FTscramble <- function(x,t=NULL,interval=NULL,spell.stats=FALSE,
                       wetfreq.pred=FALSE) {
  attributes(x) <- NULL
  n <- length(x)
  
  # This function scramles the phase information of the FT components of a
  # time series, maintaining the same spectral and time structure
  
  if (sum(is.na(x))>0) {
    ok <- is.finite(x)
    y <- approx((1:n)[ok],x[ok],xout=1:n,rule=2)$y
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