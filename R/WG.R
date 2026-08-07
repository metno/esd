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
#' \code{cluster_spells} can adjust the timing of the dry and wet days to increase 
#' the spell lengths. an example of usage:
#' y <- WG(station_data)
#' y_clustered <- cluster_spells(y, prob_shift = 0.4)
#' 
#' @aliases WG WG.station WG.fwmu.day.precip WG.FT.day.t2m FTscramble bivariate.hist
#'
#' @importFrom stats start end approx pnorm qnorm qqnorm sd dgeom rgeom rexp qexp pexp dpois
#' fft runif
#' @importFrom graphics hist
#'
#' @usage WG(x, ...)
#' @usage test.WG.fwmu.day.precip(x,n=30,verbose)
#' @usage FTscramble(x,t,interval,spell.stats,wetfreq.pred)
#' @usage bivariate.hist(x,y,plot,verbose, nx)
#' @usage cluster_spells(x, prob_shift = 0.3)
#' 
#' @param x (and y) station object
#' @param n is the number of simulations for the test of WG
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
#' @param threshold Definition of a rainy day.
#' @param method Assume a gemoetric or a poisson distribution. Can also define
#' ownth methods.
#' @param t2m station object with temperature
#' @param nx number of bins for 2D histogram
#' @param precip station object with precipitation.
#' @param mu.smudge a factor that adds noise to the ranking of mu according to climatology (\code{WG.fwmu.day.precip}). It smudges/smears out the simulated mu climatology (default: mu.smudge = 0.5*mean(mu.clim)).
#' @param start start of the year - se \code{\link{annual}} (default January 1st)
#' @param cluster an argument for \code{\link{WG.fwmu.day.precip}} to improve the spell duration statistics. 
#' @param prob_shift adjusts prob_shift (e.g., 0.2 for minor adjustments up to 0.5 for strong clustering) to match the observed mean spell length of your station.
#' @param \dots additional arguments
#' 
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
#' z2 <- WG(bjornholt, fw=0.1, mu=1)
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
#' w <- WG(bjornholt,plot=TRUE,verbose=TRUE)
#' plot(aggregate(w,by=month,FUN='wetmean')); lines(aggregate(bjornholt,by=month,FUN='wetmean'))
#' plot(aggregate(w,by=month,FUN='wetfreq')); lines(aggregate(bjornholt,by=month,FUN='wetfreq'))
#' ## Simple example of using WG to simulate climate change
#' ## We take an estimate of changed rainfall statistics from Table 1 in 
#' ## https://doi.org/10.5194/hess-29-45-2025
#' ## There is a difference in precipitation measured at Oslo-Blindern and Bjørnholt
#' ## > c(wetfreq(Oslo.Blindern),wetfreq(bjornholt)):  0.3139275 0.3717596; diff= 0.06
#' ## > c(wetmean(Oslo.Blindern),wetmean(bjornholt)): 6.804407 8.778730; diff= 1.97
#' q <- WG(bjornholt,fw=0.32 - wetfreq(bjornholt)+0.06,mu=7.24 - wetmean(bjornholt)+1.97)
#' ## Give it another time interval
#' index(q) <- index(q) - index(q)[1] + as.Date('2071-01-01') 
#' plot(zoo(combine.stations(bjornholt,q)),plot.type='single',col=c('black','grey'),
#'      main='Precipitation at Bjørnholt',xlab='',ylab='mm/day')
#' grid()
#' legend('topleft',c('Rain gauge measurement','Weather generator simulation'),
#' col=c('black','grey'),lty=1,bty='n')
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
      if (is.T(x)) y <- WG.FT.day.t2m(x,...) else
        if (is.precip(x)) y <- WG.fwmu.day.precip(x,...)
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
WG.fwmu.day.precip <- function(x=NULL,...) {
  ## Argument x is a station object with daily data
  ## Collect the arguments passed on with ...
  args <- list(...)
  plot <- args$plot; if (is.null(plot)) plot <- FALSE
  verbose <- args$verbose;if (is.null(verbose)) verbose <- FALSE
  cluster <- args$cluster; if (is.null(cluster)) cluster <- FALSE
  prob_shift <- args$prob_shift; if (is.null(prob_shift)) prob_shift <- 0.25
  mu=args$mu
  fw=args$fw
  t=args$t
  mu.smudge <- args$mu.smudge 
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
  if (is.null(mu.smudge)) mu.smudge <- 0.5*mean(mu.clim,na.rm=TRUE)
  ## If some values are undefined (it never rained), set mu.clim to zero
  coredata(mu.clim)[!is.finite(mu.clim)] <- 0
  ## Also set mu.clim to zero if there are very few events (low wet-day frequency)
  coredata(mu.clim)[fw.clim < 0.01] <- 0
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
  min.year <- max(c(min(year(fw)),min(year(mu))))
  max.year <- min(c(max(year(fw)),max(year(mu))))
  t1 <- as.Date(paste0(min.year,'-01-01'))
  t2 <- as.Date(paste0(max.year,'-01-01'))
  if (is.null(t)) t <- seq(t1,t2,by=1)
  if (verbose) {print(summary(c(fw))); print(summary(c(mu)))}
  ## Define the years
  yrs <- rownames(table(year(t)))
  if (verbose) cat('The WG simulates the years: ',range(yrs),' \n')
  ## For each year - loop
  Z <- NULL
  for (it in 1:length(yrs)) { 
    ## fw gives number of wet days per year
    if ( is.finite(fw[it]) & is.finite(mu[it]) ) { 
      n.wet <- round(fw[it] * 365.25)
      ## mu gives the amounts
      ## since we sort the data according to magnitude, we need to draw just enough numbers
      ## to correspond to the exponential distribution
      amount <- sort(rexp(n.wet,rate=1/mu[it]),decreasing=TRUE) + threshold
      p.amount <- 1 - pexp(amount,rate=1/mu[it])
      ## The scaling was defined for return values x_tau = scaling * mu * ln(tau*fw)
      ## where tau is the return period (years) and is 1/(n_tau * 365.26)
      tau.amount <- 1/(coredata(fw[it])*p.amount*365.25)
      ## Interpolate the scaling from the linear-log expression for scaling
      scaling.amount <- approx(x=x.tau,y=scaling(x.tau),xout=tau.amount)$y
      scaling.amount[!is.finite(scaling.amount)] <- 1
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
      itmu <- order(mu.jday[wet] + mu.smudge*sd(mu.clim,na.rm=TRUE)*rnorm(sum(wet)),decreasing=TRUE)
      z[wet] <- amount[order(itmu)]
      ## Ensure that the year has the same wet-day mean as prescribed
      z[wet] <- round(z[wet]*coredata(mu)[it]/mean(z[wet],na.rm=TRUE),1)
      z <- zoo(z,order.by=seq(as.Date(paste0(yrs[it],'-01-01')),as.Date(paste0(yrs[it],'-12-31')),by=1))
      if (is.null(Z)) Z <- z else Z <- c(Z,z)
      if (sum(!is.finite(z))) {
        cat('!!! NA warning !!! it=',it,'length(wet)=',length(wet),
            'length(amount)=',length(amount),'\n')
        print(summary(wet))
        print(summary(amount))
      }
      if (verbose) cat(it, yrs[it], 'fw=',fw[it], 'n.wet=',n.wet,'=',sum(wet),
                       'mu=',round(mu[it],1),'=', round(mean(z[wet]),1),'\n')
    } ## end of finite-value check
  }  ## end of loop over the years
  class(Z) <- class(x)
  Z <- attrcp(x,Z)
  if (cluster) Z <- cluster_spells(Z, prob_shift = prob_shift)
  attr(Z,'info') <- 'Weather generator'
  attr(Z,'history') <- history.stamp()
  invisible(Z)
}

#' This function tests the WG for precipitation:
#' Quantile-quantile plots of wet-day amounts
#' Number of wet days
#' @exportS3Method
#' @export test.WG.fwmu.day.precip
test.WG.fwmu.day.precip <- function(x=NULL,n=30,verbose=TRUE,
                                    cluster=FALSE,prob_shift=0.25) {
  if (verbose) print('test.WG.fwmu.day.precip')
  if (is.null(x)) {data('bjornholt',envir=environment()); x <- bjornholt; rm('bjornholt')}
  cat('test.WG.fwmu.day.precip for',loc(x),' and ',n,' runs \n')
  ## Compare the spell duration statistics for observations
  sx <- spell(x,1)
  probs_spells <- seq(0, 1, length.out = 100)
  dryx <- quantile(coredata(sx[,1]), probs = probs_spells, na.rm = TRUE)
  wetx <- quantile(coredata(sx[,2]), probs = probs_spells, na.rm = TRUE)
  probs <- seq(0, 1, length.out = 1000)
  xw_q  <- quantile(coredata(x)[x > 1], probs = probs, na.rm = TRUE)
  cat('length=',length(x),', n quantiles =',length(dryx),' = ',length(wetx),'\n')
  
  Z <- list()
  
  for (i in 1:n) {
    z <- WG(x, verbose = verbose)
    zw_q  <- quantile(coredata(z)[z > 1], probs = probs, na.rm = TRUE)
    
    # Do the same for dry/wet spells padding with NA
    sz <- spell(z, 1)
    dryz <- quantile(coredata(sz[,1]), probs = probs_spells, na.rm = TRUE)
    wetz <- quantile(coredata(sz[,2]), probs = probs_spells, na.rm = TRUE)
    
    # Inside the loop:
    Z[[paste0('run.', i, '.zw_q')]]   <- zw_q
    Z[[paste0('run.', i, '.dryz')]] <- dryz
    Z[[paste0('run.',i,'.wetz')]] <- wetz
    if (i==1) zz <- z else zz <- combine.stations(zz,z)
  }
  save(Z,zz,file='Z.test.WG.fwmu.rda')
  cat('Simulations finished - collect and inspect the results: ')
  # Extract matrix per variable (rows = days/spells, cols = runs)
  zw_mat   <- do.call(cbind, Z[grep("\\.zw_q$", names(Z))])
  dryz_mat <- do.call(cbind, Z[grep("\\.dryz$", names(Z))])
  wetz_mat <- do.call(cbind, Z[grep("\\.wetz$", names(Z))])
  
  # Compute ensemble statistics across simulation runs
  cat(' ensemble statistics')
  zw_mean   <- rowMeans(zw_mat, na.rm = TRUE)
  zw_bounds <- apply(zw_mat, 1, quantile, probs = c(0.05, 0.95), na.rm = TRUE)
  
  dryz_mean <- rowMeans(dryz_mat, na.rm = TRUE)
  dry_bounds <- apply(dryz_mat, 1, quantile, probs = c(0.05, 0.95), na.rm = TRUE)
  wetz_mean <- rowMeans(wetz_mat, na.rm = TRUE)
  wet_bounds <- apply(wetz_mat, 1, quantile, probs = c(0.05, 0.95), na.rm = TRUE)
  
  cat(', plot \n')
  par(mfrow = c(2, 2), mar = c(4, 5, 2.5, 1), 
      cex.axis=1.25,cex.lab=1.25)
  
  # --- Panel 1: Wet-day Amounts (Q-Q Plot) ---
  lim_w <- range(c(xw_q, zw_mat), na.rm = TRUE)
  plot(xw_q, zw_mean, type = "n", xlim = lim_w, ylim = lim_w,
       xlab = "Observed amount (mm/day)", ylab = "WG amount (mm/day)",
       main = paste(loc(x), "Wet-day Amounts"))
  grid()
  
  # 95% Ensemble range polygon & individual run lines
  polygon(c(xw_q, rev(xw_q)), c(zw_bounds[1, ], rev(zw_bounds[2, ])), 
          col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
  matlines(xw_q, zw_mat, col = adjustcolor("gray60", alpha.f = 0.3), lty = 1)
  
  # Ensemble mean & 1:1 reference line
  lines(xw_q, zw_mean, col = "blue", lwd = 2)
  abline(0, 1, lty = 2, col = "red")
  legend("topleft", legend = c("Ensemble Mean", "5-95% Range", "1:1 Line"),
         col = c("blue", "skyblue", "red"), lty = c(1, 1, 2), lwd = c(2, 6, 1), bty = "n")
  
  # --- Panel 2: Dry Spell Duration Q-Q Plot ---
  lim_d <- range(c(dryx, dryz_mat), na.rm = TRUE)
  plot(dryx, dryz_mean, type = "n", xlim = lim_d, ylim = lim_d,
       xlab = "Observed dry spell (days)", ylab = "WG dry spell (days)",
       main = "Dry Spell Duration")
  grid()
  polygon(c(dryx, rev(dryx)), c(dry_bounds[1, ], rev(dry_bounds[2, ])), 
          col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
  matlines(dryx, dryz_mat, col = adjustcolor("gray60", alpha.f = 0.3), lty = 1)
  lines(dryx, dryz_mean, col = "blue", lwd = 2)
  abline(0, 1, lty = 2, col = "red")
  
  # --- Panel 3: Wet Spell Duration Q-Q Plot ---
  lim_w <- range(c(wetz_mat), na.rm = TRUE) # add wetx if defined in your workspace
  plot(wetx, wetz_mean, type = "n", ylim = lim_w,
       xlab = "Observed wet spell (days)", ylab = "WG wet spell (days)",
       main = "Wet Spell Duration")
  grid()
  polygon(c(wetx, rev(wetx)), c(wet_bounds[1, ], rev(wet_bounds[2, ])), 
          col = adjustcolor("skyblue", alpha.f = 0.4), border = NA)
  matlines(wetx,wetz_mat, col = adjustcolor("gray60", alpha.f = 0.3), lty = 1)
  lines(wetx, wetz_mean, col = "blue", lwd = 2)
  abline(0, 1, lty = 2, col = "red")
  
  # --- Panel 4: Dry Spell Duration Q-Q Plot ---
  xz <- combine.stations(x,zz)
  col <- c('black',rep(rgb(1,0,0,0.5),n))
  nwxz <- zoo(annual(xz, FUN = "count", threshold = 1))
  plot(nwxz,main='Annual number of wet days',col=col,
       plot.type='single',ylab='days',lty=c(1,rep(2,n)))
  legend('topleft',c('Original','WG'),col=c('black',rgb(1,0,0,0.5)),
         lty=c(1,2),bty='n',cex=0.6)
  grid()
  lines(subset(nwxz,is=1))
  
  ## Compare the mean annual cycle of observations and simulations
  
  ## Compare annual statistics
  axz <- zoo(annual(xz,FUN='sum'))
  nyrs <- length(index(axz))
  plot(axz,main='Annual total precipitation',col=col,
       plot.type='single',ylab=expression(sum(x)*phantom(0)*(mm/day)),
       lty=c(1,rep(2,n)))
  legend('topleft',c('Original','WG'),col=c('black',rgb(1,0,0,0.5)),
                                            lty=c(1,2),bty='n',cex=0.6)
  grid()
  lines(subset(axz,is=1))
  
  ## Compare mean seasonal
  plot(zoo(aggregate(xz,by=month,FUN='sum'))/nyrs,main='Seasonal total precipitation',
       ylab=expression(sum(x)*phantom(0)*(mm/day)),col=col,plot.type='single',
       lty=c(1,rep(2,n))); grid()
  legend('topleft',c('Original','WG'),col=c('black',rgb(1,0,0,0.5)),
                                            lty=c(1,2),bty='n',cex=0.6)
  plot(zoo(aggregate(xz,by=month,FUN='wetfreq')),main='Seasonal wet-day frequency',
       ylab=expression(f[w]),col=col,plot.type='single',lty=c(1,rep(2,n)),
       ylim=c(0,1)); grid()
  legend('topleft',c('Original','WG'),col=c('black',rep(rgb(1,0,0,0.5)),
                                            lty=c(1,2)),bty='n',cex=0.6)
  plot(zoo(aggregate(xz,by=month,FUN='wetmean')),main='Seasonal wet-day mean',
       ylab=expression(mu*phantom(0)*(mm/day)),col=col,plot.type='single',
       lty=c(1,rep(2,n))); grid()
  legend('topleft',c('Original','WG'),col=c('black',rgb(1,0,0,0.5)),
                                            lty=c(1,2),bty='n',cex=0.6)
  
  ## Quantitative statistics for comparing the data
  # 1. Aggregate to monthly totals
  
  clim_sum <- aggregate(as.monthly(xz,FUN='sum'), by = month, FUN = "mean")
  clim_mu <- aggregate(xz, by = month, FUN = "wetmean")
  clim_fw <- aggregate(xz, by = month, FUN = "wetfreq")
  
  # 2. Separate observed (1st column) and ensemble simulations (rest of columns)
  obs_sum  <- clim_sum[, 1];  sims_sum <- clim_sum[, -1]
  obs_mu   <- clim_mu[, 1];   sims_mu  <- clim_mu[, -1]
  obs_fw   <- clim_fw[, 1];   sims_fw  <- clim_fw[, -1]
  
  # 3. Calculate RMSE across each ensemble member (vectorized across columns)
  rmse_q <- sqrt(colMeans((zw_mat - xw_q)^2, na.rm = TRUE))
  rmse_wet <- sqrt(colMeans((wetz_mat - wetx)^2, na.rm = TRUE))
  rmse_dry <- sqrt(colMeans((dryz_mat - dryx)^2, na.rm = TRUE))
  rmse_sum <- sqrt(colMeans((sims_sum - obs_sum)^2, na.rm = TRUE))
  rmse_mu <- sqrt(colMeans((sims_mu - obs_mu)^2, na.rm = TRUE))
  rmse_fw <- sqrt(colMeans((sims_fw - obs_fw)^2, na.rm = TRUE))
  
  # 4. Extract statistical summary across the ensemble
  cat('RMSE for daily wet-day quantiles: \n')
  print(summary(rmse_q))
  cat('RMSE for wet-spell duration quantiles: \n')
  print(summary(rmse_wet))
  cat('RMSE for dry-spell duration quantiles: \n')
  print(summary(rmse_dry))
  cat('RMSE for climatology in total precipitation: \n')
  print(summary(rmse_sum))
  cat('RMSE for climatology in wet-day mean precipitation: \n')
  print(summary(rmse_mu))
  cat('RMSE for climatology in wet-day frequency: \n')
  print(summary(rmse_fw))
  attr(xz,'location') <- paste(loc(xz),c('obs',
                                         paste(rep('WG',n),1:n)))
  invisible(xz)
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
bivariate.hist <- function(x,y,plot=TRUE,verbose=FALSE, nx=30) {
  ## Standardising
  if (verbose) cat('bivariate.hist: \n')
  stand <- function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  
  ## Keep a copy
  pr <- x; tas <- y
  ## Ensure rigth format: zoo-objects
  x <- zoo(x)
  y <- stand(zoo(y))
  
  ## Make 2D bins and count events in each 
  ny <- nx
  seqx <- seq(floor(min(x,na.rm=TRUE)),ceiling(max(x,na.rm=TRUE)),length=nx+1)
  seqy <- seq(floor(min(y,na.rm=TRUE)),ceiling(max(y,na.rm=TRUE)),length=ny+1)
  seqX <- seq(floor(min(x,na.rm=TRUE)),ceiling(max(x,na.rm=TRUE)),length=10*nx)
  seqY <- seq(floor(min(y,na.rm=TRUE)),ceiling(max(y,na.rm=TRUE)),length=10*ny)
  Z <- matrix(rep(NA,nx*ny),nx,ny)
  for (i in 1:nx)
    for (j in 1:ny) Z[i,j] <- sum( (x >= seqx[i]) & (x < seqx[i+1]) &
                                     (y >= seqy[j]) & (y < seqy[j+1]), na.rm=TRUE )
  if (verbose) print(summary(c(Z)))
  
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

# Function to increase wet and dry spell clustering
#' @export
cluster_spells <- function(x, prob_shift = 0.3) {
  # x: zoo/station object from esd::WG or numeric vector
  # prob_shift: fraction of isolated wet days (1-day spells) to move (0.0 to 1.0)
  
  vals <- coredata(x)
  n <- length(vals)
  is_wet <- vals > 0
  
  # Identify isolated wet days (surrounded by dry days on both sides)
  isolated <- which(is_wet & 
                      c(FALSE, !is_wet[-n]) & 
                      c(!is_wet[-1], FALSE))
  
  if (length(isolated) == 0) return(x)
  
  # Select a random subset of isolated days to relocate
  n_move <- round(prob_shift * length(isolated))
  to_move <- sample(isolated, size = n_move)
  
  # Move precipitation from isolated days to dry days adjacent to existing wet days
  for (idx in to_move) {
    # Find dry days that are immediately adjacent to at least one wet day
    adj_dry <- which(!is_wet & (c(FALSE, is_wet[-n]) | c(is_wet[-1], FALSE)))
    if (length(adj_dry) == 0) break
    
    # Pick a random adjacent dry day and transfer the daily amount
    target <- sample(adj_dry, 1)
    vals[target] <- vals[idx]
    vals[idx] <- 0
    
    # Update state vector
    is_wet[target] <- TRUE
    is_wet[idx] <- FALSE
  }
  
  # Reassign modified values while preserving esd object attributes and dates
  x[] <- vals
  return(x)
}

