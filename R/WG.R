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
#' \code{WG.fw.day.precip} has been designed to use with downscaled results for 
#' annual wet-day frequency and annual wet-day mean precipitation. 
#' It uses the annual wet-day mean and the wet-day
#' frequency as input, and takes a sample station of daily values to
#' simulate stochastic numbers of consecutive wet days based on its annual mean
#' number of consecutive wet days. If not specified, it is taken from the sample 
#' data after being phase scrambled (\code{FTscramble}). 
#' The number of wet-spells per year is estimated from the wet-day 
#' frequency divided by annual mean wet-spell duration. The annual wet-day frequency  
#' is also taken to be phase scrambled estimates from the sample data unless 
#' specifically specified. The algorithm uses the mean duration between the 
#' start of each precipitation event (Dpe) in addition to the annual number of 
#' precipitation events (Npe). The daily amount is taken from stochastic values 
#' generated with \code{\link{rexp}} scaled for the tail according to alpha in 
#' (described in DOI: 10.1088/1748-9326/ab2bb2) as in \code{\link{day2IDF}}. 
#' The number of consecutive wet days and duration between start of each event 
#' can be approximated by a geometric distribution (\code{\link{rgeom}}), and the 
#' annual mean number was estimated from the sample series.
#' \code{test.WG.fw.day.precip} presents diagnostics of tests of \code{WG.fw.day.precip}.
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
#' @param ndbr Number of 
#' @param n.spells.year = c('fw','spell') if 'fw' then estimate number of spells according to 365.25 else estimate number of events from \code{\link{spell}}.
#' @param alpha.scaling TRUE scale the low-probability events according to alpha in DOI:10.1088/1748-9326/ab2bb2
#' @param alpha values for alpha-scaling 
#' @param ensure.fw TRUE then WG tries to ensure that fw of simulations match those of observations or prescribed by adding or substracting wet days.
#' @param \dots additional arguments
#' @author R.E. Benestad
#' @keywords manip
#' @examples
#' ## Temperature
#' data(ferder)
#' x <- WG(ferder)
#' ## Plot the results
#' plot(merge(ferder,x),xlab='',ylab=c('Obs T2m','WG T2m'), col='blue',main=paste(loc(y),' Obs/WG'))
#' 
#' ## Daily precipitation
#' data(bjornholt)
#' z <- WG(bjornholt)
#' ## Plot the results
#' plot(merge(bjornholt,z),xlab='',ylab=c('Obs precip','WG precip'), col='blue',main=paste(loc(y),' Obs/WG'))
#' sz <- sort(coredata(z))
#' sy <- sort(coredata(bjornholt))
#' ## Use WG to 'simulate' climate change
#' z2 <- WG(bjornholt,mu=annual(bjornholt,FUN='wetmean')+2)
#' sz2 <- sort(coredata(z2))
#' ## Plot the comparison of quantiles
#' plot(sy,sz,pch=19,cex=0.7,main='QQ-plot',xlab='Observations',ylab='WG')
#' grid()
#' lines(c(0,max(sy,sz,na.rm=TRUE)),c(0,max(sy,sz,na.rm=TRUE)),lty=2,col='red')
#' points(sy,sz2,col='blue',cex=0.7)
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

## Future considerations -also allow for estimating the AR(1) coefficient of the Hurst coefficient?
## Fractional Gaussian noise...?
## N wet days from fw
## Timing from spell
## Duration from spell.
## --- Precipitation  
#' @exportS3Method
#' @export WG.fw.day.precip
WG.fw.day.precip <- function(x=NULL,...,mu=NULL,fw=NULL,ndbr=NULL,t=NULL,
                             threshold=1,alpha.scaling=TRUE,alpha=c(1.256,0.064),
                             plot=FALSE,verbose=FALSE) {
  
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
  x.mu <- as.annual(x,'wetmean',threshold=threshold)
  
  # Number of consecutive wet/dry days
  ncd <- spell(x,threshold=threshold)
  ## Annual mean number of consecutive wet days
  amncwd <- subset(annual(ncd),is=1)
  # extract the time interval between the start of each dry spell
  dt1 <- diff(julian(as.Date(index(ncd[is.finite(ncd[,1]),1]))))
  
  if (plot) {
    ## Timing between each precipitation event
    dev.new()
    par(mfrow=c(2,2),cex.main=0.7)
    f.k <- dgeom(0:max(dt1), prob=1/(mean(dt1)+1))
    hist(dt1,freq=FALSE,col="grey",xlab="days",
         main="The time between the start of each precipitation event",
         sub="Test: Red curve is the fitted geometric distribution")
    lines(0:max(dt1),f.k,lwd=5,col="red")
    grid()
  }
  
  ## Annual mean number of days between start of each rain event
  ## Remove first and last elements to avoid cut-off problems at start and
  ## end of the time series
  amndse <- annual(zoo(x=dt1,order.by=index(dt1)))[-c(1,length(dt1))]
  
  ## Wet-day spell duration statistics:
  wetsd <-   subset(ncd,is=1)
  ## Remove the first and last estimate to avoid cut-off problems
  wetsd <-  subset(wetsd,it=c(FALSE,rep(TRUE,length(wetsd)-2),FALSE))
  amwetsd <- annual(wetsd,FUN='mean',nmin=1)
  ## Annual number of wet events 
  nwes <- aggregate(wetsd,by=year(wetsd),FUN="nv")
  if (plot) {
    ## Number of events per year
    hist(coredata(nwes),breaks=seq(0,100,by=5),freq=FALSE,col="grey",
         main="Number of wet events per year",xlab="days",
         sub="Test: Red curve is the fitted Poisson distribution")
    lines(seq(0,100,by=1),dpois(seq(0,100,by=1),lambda=mean(coredata(nwes))),
          col="red",lwd=3)
    grid()
  }
  
  ## Estimate climatology for mean seasonal cycle in total precipitation. Use this information
  ## as a guide for which months to add wet days to ensure correct wet-day frequency fw
  pt.ac <- aggregate(y,month,FUN='mean',na.rm=TRUE)
  
  # Wet-day mean: from DS or from observations
  if (verbose) print('wet-day mean')
  if (is.null(mu))
    mu <- zoo(FTscramble(x.mu),order.by=index(x.mu)) 
  rm('x.mu')
  
  # Wet-day frequency: from DS or from observations
  if (verbose) print('wet-day frequency')
  if (is.null(fw)) fw <- zoo(FTscramble(x.fw),order.by=index(x.fw)) 
  rm('x.fw')
  
  ## Use the time between start of each wet spell and statistics of the wet spell duration 
  ## to populate the year with wet days
  if (verbose) print('random annual mean number of n_cwd:')
  ## Stochastic annual mean wet-spell duration from number of consecutive wet day x.nd:
  amwsd <- rnorm(length(amwetsd),mean=mean(coredata(amwetsd),na.rm=TRUE),
               sd=sd(coredata(amwetsd),na.rm=TRUE))
  if (plot) {
    ## Number of events per year
    hist(coredata(wetsd),breaks=seq(0,40,by=2),freq=FALSE,col="grey",
         main="Duration of wet spells",xlab="days",
         sub="Test: Red curve is the fitted geometric distribution")
    lines(seq(0,40,by=1),dgeom(seq(0,40,by=1),prob=1/mean(coredata(wetsd))),
          col="red",lwd=3)
    grid()
  }
  ## Constraint: at least one wet event per year
  amwsd[amwsd < 1] <- 1;
  if (verbose) print(amwsd)

  ## Time axis for projection:
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
  
  # Estimate the annual number of rainy days:
  if (verbose) print('Number of wet days each year:')
  anwd <- round( ( julian(as.Date(paste(year(fw),'12-31',sep='-'))) -
                   julian(as.Date(paste(year(fw),'01-01',sep='-'))) + 1) *
                   coredata(fw) )
  
  # Error bars for mu: var = mu**2 for exponential distribution:
  mu.err <- mu/sqrt(anwd - 1)
  
  # set up a record with no rain:
  z <- zoo(rep(0,length(t)),order.by=t)
  j <- 1:366
  
  # add rain events:
  if (verbose) print(paste('loop over year:',1,'-',ny))
  for ( i in 1:ny ) {
    ## Simulate daily: first sort days into wet and dry days, first round based on the statistics of the timing between the 
    ## start of each wet event and its duration
    
    ## Time between the start of each event
    tbsee <- rgeom(366,prob=1/(amndse[i])) + 1
    ## Duration of wet events
    ncwd <- rgeom(366,prob=1/(amwsd[i])) + 1
    
    ## Find number of events needed to get right number of wet days
    i2 <- cumsum(ncwd) >= anwd[i]
    #if (verbose) print(table(i2))
    nes <- min(j[i2],na.rm=TRUE)
    #if (verbose) print(paste(anwd[i],'number of wet days in',nes,'wet events'))
    ## Keep to number of events: timing between them and duration
    tbsee <- tbsee[1:nes]; ncwd <- ncwd[1:nes]
    ## Go through each event:
    dry <- c(); wet <- c(); t2 <- 1
    for (ie in 1:nes) {
      dry <- c( dry, t2 + 0:(tbsee[ie]-1) )
      wet <- c( wet, t2 + tbsee[ie] + 0:(ncwd[ie]-1) )
      t2 <- t2 + tbsee[ie] + ncwd[ie]
    }
    # if (verbose) {
    #   print(paste('Number of wet days in',year(anwd[i]),anwd[i],length(wet)))
    #   print('Wet days'); print(wet); print('Dry days'); print(dry)
    #   print('Missing assignements'); print(setdiff(j,sort(c(dry,wet))))
    #   print('-----------------')
    # }
    ## Need to pad up so the total number of wet days matches the wet-day frequency fw
    if (anwd[i] > length(wet)) {
      nswap <- anwd[i] - length(wet)
      ## Pick random dry days
      swap <- dry[order(rnorm(length(dry)))][1:nswap]
      dry <- dry[-swap]
      wet <- sort(c(wet,swap))
    } else if (anwd[i] < length(wet)) {
      ## Set excess random rainy days to dry days
      nswap <- length(wet) - anwd[i]
      swap <- wet[order(rnorm(length(wet)))][1:nswap]
      wet <- wet[-swap]
      dry <- sort(c(dry,swap))
    } 
    
    ## The wet-day mean precipitation amount
    if (!is.finite(mu[i])) mu[i] <- mean(mu,na.rm=TRUE)
    ## The daily amounts for wet days
    y <- round(rexp(366,rate=1/coredata(mu[i])),1)
    if (alpha.scaling) {
      ## REB 2024-05-13
      ## Scale the amounts according to return-period according to 
      ## DOI:https://doi.org/10.1088/1748-9326/ab2bb2 see day2IDF
      ## tau - return-interval in years 
      #if (verbose) print('Scale by alpha according to return-interval')
      tau <- 1/(365.25*(1- pexp(y,rate=1/coredata(mu[i]))))
      ## Take into account the fraction of wet days
      tau <- tau/coredata(fw[i])
      alphas <- alpha[1] + alpha[2]*log(tau)
      alphas[alphas < 1] <- 1
      y <- y * alphas
      #if (verbose) {print(summary(tau)); print(summary(alphas))}
    }
    if (plot & (i==1)) {
      z <- coredata(subset(x,it=rep(year(x)[1],2)))
      z <- z[z >= 1]
      plot(sort(z),sort(y[1:length(z)]),main=paste('Wet-day amounts (mm) for',year(x)[1]),
           xlab='Observed',ylab='WG')
      grid()
      maxzy <- max(z,y,na.rm=TRUE)
      lines(c(0,maxzy),c(0,maxzy),lty=2,col='red')
    }
    # add rain to the appropriate year:
    ii <- is.element(year(t),yrs[i])
    rain <- rep(0,sum(ii)); iii <- 0
    rain[wet] <- y[1:length(wet)]
    ## Make it a zoo object to assign months
    #if (verbose) print(range(as.Date(paste0(year(fw[i])-1,'-12-31'))+1:length(rain)))
    #rain <- zoo(rain,order.by=as.Date(paste0(year(fw[i])-1,'-12-31'))+1:length(rain))
    
    if (verbose) print(paste(yrs[i],'tot rain',round(sum(rain,na.rm=TRUE)),
                             '#wet days=',length(wet),'n*fw[i]=',round(365.25*fw[i]),
                             'max(wet)=',max(wet),'mu[i]=',round(mu[i],1),'#events=',nes,
                             'length(dry)=',length(dry),'max(dry)=',max(dry)))
    z[ii] <- rain
  }
  z <- zoo(z,order.by=t)
  class(z) <- class(x)
  z <- attrcp(x,z)
  attr(z,'original_fw') <- fw
  attr(z,'n_cwd') <- amwetsd
  attr(z,'t_between_events') <- amndse
  attr(z,'original_mu') <- mu
  attr(z,'mu_error') <- mu.err
  attr(z,'aspect') <- paste(attr(z,'aspect'),'weather_generator',sep=', ')
  attr(z,'history') <- history.stamp(x)
  return(z)
}

#' This function tests the WG for precipitation:
#' Quantile-quantile plots of wet-day amounts
#' Number of wet days
#' @exportS3Method
#' @export test.WG.fw.day.precip
test.WG.fw.day.precip <- function(x) {
  z <- WG(x)
  ## sort magnitudes to plot quantile-quantile plots
  xw <- sort(coredata(x[x > 1]))
  zw <- sort(coredata(z[z > 1]))
  ## There may be different number of wet days - pad the shortest series with 0s
  nx <- lenth(xw)
  nz <- length(zw)
  if (nx > nz) zw <- c(rep(0,nx-nz),zw) else 
    if (nz > nx) xw <- c(rep(0,nz-nx),xw)
  xylim <- c(0,max(c(xw,zw)))
  par(mfcol=c(2,1))
  plot(xw,zw,ylim=xylim,xlim=xylim,xlab='Observed amount (mm/day)',
       ylab='WG amount (mm/day)',main=paste(loc(x),'wet-day amounts'))
  grid()
  lines(xylim,xylim,lty=2,col='red')
  ## compare the number of wet days
  xnw <- zoo(annual(x,FUN='count',1))
  znw <- zoo(annual(z,FUN='count',1))
  plot(merge(xnw,znw),plot.type='single',col=c('black','red'),lty=c(1,2),
       main='Number of annual wet days')
  grid()
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