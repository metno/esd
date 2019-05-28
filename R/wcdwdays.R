' Projection of hot and cold day statistics
#' 
#' The functions \code{hotsummerdays}, \code{heatwavespells},
#' \code{coldwinterdays}, and \code{coldspells} estimate statistics for
#' heatwaves/hot days or cold spells based on seasonal mean temperatures. The
#' estimations are based on a regression analysis (GLM) between observed number
#' of events or spell lengths and seasonal mean from station data.
#' \code{nwetdays} estimates the number of days per year with precipitation
#' amount exceeding a threshold values.
#' 
#' The estimation of these statistics makes use of general linear models (GLMs)
#' and take the counts to follow the 'Poisson family' whereas the spall lengths
#' belong to the geometric distribution. The seasonal mean temperature or
#' annual wet-mean precipitation are used as independent variable.
#' 
#' @aliases hotsummerdays coldwinterdays coldspells heatwavespells nwetdays
#' @param x station object, e.g. the temperature. Matches the element used in
#' the dsensemble object 'dse'
#' @param y station object which may be some statistics with dependency to x,
#' e.g. snow depth.
#' @param dse a dsensembel object. If NULL, then run DSensemble
#' @param it Default season set for northern hemisphere. %% ~~Describe
#' \code{it} here~~
#' @param threshold Temperature threshold
#' @param verbose TRUE for trouble shooting, debugging etc.
#' @param plot TRUE - produce graphics
#' @author R.E. Benestad
#' @examples
#' 
#' data(ferder)
#' data(dse.ferder)
#' hw <- hotsummerdays(ferder,dse.ferder,threshold=20)
#' 
#' @export hotsummerdays
hotsummerdays <- function(x,y=NULL,dse=NULL,it='jja',threshold=30,
                          verbose=FALSE,plot=TRUE,nmin=90,new=TRUE,...) {
    # Estimate number of days with low temperatures
  if (verbose) print('hotsummerdays')
  if ( (inherits(y,'dsensemble')) & is.null(dse)) {
    ## Swap y & dse.
    print('use y  to set dse')
    dse <- y; y <- NULL
  }
  stopifnot(inherits(x,'station'))
  if (is.null(y)) y <- x
  djf <- subset(x,it=it)      # default: summer
  djfy <- subset(y,it=it)     # default: summer
  nwd1 <- annual(djfy,FUN='count',threshold=threshold,nmin=nmin)
  mwd1 <- annual(djfy,FUN='mean',nmin=nmin)
  ## REB 2016-11-17: exclude temperatures far from the threshold:
  xcld1 <- (mwd1 < threshold - 15) | (mwd1 > threshold + 20)
  xcld1 <- xcld1[is.finite(xcld1)]
  if (verbose) print(paste('Exclude',sum(xcld1),'outliers'))
  if (sum(xcld1)>0) mwd1[xcld1] <- NA
  
  if (verbose) print(c(length(mwd1),length(nwd1)))
  cal <- data.frame(x=c(coredata(mwd1)),
                    y=c(coredata(nwd1)))
  ## Use linear fit rather than polynomial as this analysis only
  ## involves one season and to a greater degree extrapolation.
  dfit <- glm(y ~ x,family='poisson',data=cal)

  if (plot) {
    if (new) dev.new()
    par(bty='n')
    plot(cal,pch=19,ylim=c(0,90),
         xlab=expression(paste('mean temperature ',(degree*C))),
         ylab='number of hot days',main=loc(x))
    pre <- data.frame(x=seq(min(cal$x,na.rm=TRUE)-1,
                            max(cal$x,na.rm=TRUE)+5,by=0.1))
    lines(pre$x,exp(predict(dfit,newdata=pre)),col=rgb(1,0,0,0.3),lwd=3)

    if (new) dev.new()
    ## Use the distribution about the mean
    djf.sd <- sd(coredata(djf),na.rm=TRUE)
    ## Check it it is normally distributed
    qqnorm(coredata(djf))
    qqline(coredata(coredata(djf)),col='red')
    grid()
  }

  if (is.null(dse)) dse <-  DSensemble.t2m(x,biascorrect=TRUE,
                                           verbose=verbose,plot=plot)
  if (length(table(month(dse)))==4) djf.dse <- subset(dse,it='jja') else djf.dse <- dse
  index(djf.dse) <- year(djf.dse)
  ovl <- window(djf.dse,start=year(start(x)),end=year(end(x)))
  djf.dse <- djf.dse - mean(coredata(ovl),na.rm=TRUE) +
                       mean(coredata(mwd1),na.rm=TRUE)

  q1 <- data.frame(x=apply(coredata(djf.dse),1,quantile,probs=0.05,na.rm=TRUE))
  q2 <- data.frame(x=apply(coredata(djf.dse),1,quantile,probs=0.95,na.rm=TRUE))
  qm <- data.frame(x=apply(coredata(djf.dse),1,mean,na.rm=TRUE))
  obs <- data.frame(x=coredata(mwd1))

  t <- year(index(djf.dse))
  preq1 <- exp(predict(dfit,newdata=q1))
  preq1[preq1 > 92] <- 92  # maximum length of the summer season
  tr1 <- predict(lm(preq1 ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
  tr1[!is.finite(preq1)] <- NA
  preq2 <- exp(predict(dfit,newdata=q2))
  preq2[preq2 > 92] <- 92  # maximum length of the summer season
  tr2 <- predict(lm(preq2 ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
  tr2[!is.finite(preq2)] <- NA
  prem  <- exp(predict(dfit,newdata=qm))
  prem[prem > 92] <- 92  # maximum length of the summer season 
  tr3 <- predict(lm(prem ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
  tr3[!is.finite(prem)] <- NA
  Nwd <- zoo(cbind(preq1,preq2,prem,tr1,tr2,tr3),order.by=t)
  nwd.pre <- zoo(exp(predict(dfit,newdata=obs)),order.by=year(mwd1))
  
  
  if (plot) {
    if (verbose) print('plots')
    if (new) dev.new()
    par(bty='n')
    plot(zoo(djf.dse,order.by=year(djf.dse)),
         plot.type='single',col=rgb(0.5,0.5,0.5,0.2),
         ylab=expression(paste('mean temperature',(degree*C))),
         xlab='',main=loc(x))
    points(mwd1,pch=19)
    grid()
    
    if (new) dev.new()
    par(bty='n')
    plot(Nwd,plot.type='single',lwd=5,main=loc(x),ylim=c(0,90),
         xlab="",ylab=paste('number of hot days: T(2m) > ',threshold,unit(x)),
         col=c(rgb(0.5,0.5,0.7,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.8,0.5,0.8,0.5),
               rgb(0.3,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.6,0.3,0.6,0.5)),
         ...)
    grid()
    points(nwd1,pch=19)
    lines(nwd.pre,col=rgb(0.5,0.5,0.5,0.5))
  }

  Nwd <- attrcp(x,Nwd)
  attr(Nwd,'unit') <- 'days'
  attr(Nwd,'info') <- paste('number of hot days: t2m > ',threshold)
  attr(Nwd,'observation') <- nwd1
  attr(Nwd,'nwd.pre') <- nwd.pre
  index(Nwd) <- t
  class(Nwd) <- c('nevents','zoo')
  if (verbose) print('exit hotsummerdays')
  invisible(Nwd)
}

#' @export
coldwinterdays <- function(x,y=NULL,dse=NULL,it='djf',threshold=0,
                           verbose=FALSE,plot=TRUE,nmin=90,new=TRUE,...) {
  # Estimate number of days with low temperatures or number of days with y < threshold
  if (verbose) print('coldwinterdays')
  stopifnot(inherits(x,'station'))
  if (is.null(y)) y <- x
  djf <- subset(x,it=it)     # Winter
  djfy <- subset(y,it=it)     # Winter
  mam <- subset(x,it='mam')  # Spring/autumn
  mamy <- subset(y,it='mam')  # Spring/autumn
  nwd1 <- annual(-djfy,FUN='count',threshold=threshold,nmin=nmin)
  mwd1 <- annual(djf,FUN='mean',nmin=nmin)
  nwd2 <- annual(-mamy,FUN='count',threshold=threshold,nmin=nmin)
  mwd2 <- annual(mam,FUN='mean',nmin=nmin)

  ## REB 2016-11-17: exclude temperatures far from the threshold:
  xcld1 <- (mwd1 < threshold - 15) | (mwd1 > threshold + 20)
  if (verbose) print(paste('Exclude',sum(xcld1),'outliers'))
  mwd1[xcld1] <- NA
  xcld2 <- (mwd2 < threshold - 15) | (mwd2 > threshold + 20)
  if (verbose) print(paste('Exclude',sum(xcld2),'outliers'))
  mwd1[xcld2] <- NA
  
  cal <- data.frame(x=c(coredata(mwd1),coredata(mwd2)),
                    y=c(coredata(nwd1),coredata(nwd2)))
  ## Use polymomial as two different seasons are involved and the
  ## analysis involves a interpolation more than an extrapolation
  ## since winter is expected to become more similar to spring/autumn
  dfit <- glm(y ~ x + I(x^2) + I(x^3),family='poisson',data=cal)
  
  if (plot) {
    if (verbose) print('plot calibration results')
    if (new) dev.new()
    par(bty='n')
    plot(cal,ylim=c(0,90),xlim=range(c(coredata(mwd1),coredata(mwd2)),na.rm=TRUE),pch=19,
         xlab=expression(paste('mean temperature ',(degree*C))),
         ylab='number of cold days',main=loc(x))
    pre <- data.frame(x=seq(min(c(coredata(mwd1),coredata(mwd2)),na.rm=TRUE),
                            max(c(coredata(mwd1),coredata(mwd2)),na.rm=TRUE),by=0.1))
    lines(pre$x,exp(predict(dfit,newdata=pre)),col=rgb(1,0,0,0.3),lwd=3)

    if (new) dev.new()
    ## Use the distribution about the mean
    djf.sd <- sd(coredata(djf),na.rm=TRUE)
    ## Check it it is normally distributed
    qqnorm(coredata(djf))
    qqline(coredata(coredata(djf)),col='red')
    grid()
  }

  if (verbose) print(class(dse))
  if (is.null(dse)) dse <-  DSensemble.t2m(x,biascorrect=TRUE,
                                           verbose=verbose,plot=plot)
  if (verbose) {print('Seasons/annual?'); print(table(month(dse))); print(class(dse))}
  if (length(table(month(dse)))==4) djf.dse <- subset(dse,it='djf') else djf.dse <- dse
  index(djf.dse) <- year(djf.dse)
  ovl <- window(djf.dse,start=year(start(x)),end=year(end(x)))
  djf.dse <- djf.dse - mean(coredata(ovl),na.rm=TRUE) +
                       mean(coredata(mwd1),na.rm=TRUE)
  if (verbose) {str(ovl); print(mean(coredata(ovl),na.rm=TRUE)); print(mean(coredata(mwd1),na.rm=TRUE))}

  t <- year(index(djf.dse))
  q1 <- data.frame(x=apply(coredata(djf.dse),1,quantile,probs=0.05,na.rm=TRUE))
  q2 <- data.frame(x=apply(coredata(djf.dse),1,quantile,probs=0.95,na.rm=TRUE))
  qm <- data.frame(x=apply(coredata(djf.dse),1,mean,na.rm=TRUE))
  obs <- data.frame(x=coredata(mwd1))
  
  ## If the values are far from zero, set to NA
  #q1$x[(q1$x < threshold - 15) | (q1$x > threshold + 20)] <- NA
  #q2$x[(q2$x < threshold - 15) | (q2$x > threshold + 20)] <- NA
  #qm$x[(qm$x < threshold - 15) | (qm$x > threshold + 20)] <- NA
  
  if (verbose) str(qm)
  
  if (verbose) {print('Fit trend cubic models'); range(t)}
  if (verbose) print('ensemble q05')
  preq1 <- exp(predict(dfit,newdata=q1))
  preq1[preq1 > 92] <- 92  # maximum length of the winter season
  tr1 <- predict(lm(preq1 ~ t + I(t^2) + I(t^3)))
  tr1[!is.finite(preq1)] <- NA
  if (verbose) print('ensemble q95')
  preq2 <- exp(predict(dfit,newdata=q2))
  preq2[preq2 > 92] <- 92  # maximum length of the winter season
  tr2 <- predict(lm(preq2 ~ t + I(t^2) + I(t^3)))
  tr2[!is.finite(preq2)] <- NA
  if (verbose) print('ensemble mean')
  prem  <- exp(predict(dfit,newdata=qm))
  prem[prem > 92] <- 92  # maximum length of the winter season
  tr3 <- predict(lm(prem ~ t + I(t^2) + I(t^3)))
  tr1[!is.finite(prem)] <- NA
  
  if (verbose) print('combine predictions/trends into one object respectively')
  Nwd <- zoo(cbind(preq1,preq2,prem,tr1,tr2,tr3),order.by=t)
  nwd.pre <- zoo(exp(predict(dfit,newdata=obs)),order.by=year(mwd1))
  
  if (plot) {
    if (verbose) print('plots')
    if (new) dev.new()
    par(bty='n')
    plot(zoo(djf.dse,order.by=year(djf.dse)),
         plot.type='single',col=rgb(0.5,0.5,0.5,0.2),
         ylab=expression(paste('mean temperature - winter',(degree*C))),
         xlab='',main=loc(x))
    points(mwd1,pch=19)
    grid()
    
    if (new) dev.new()
    par(bty='n')
    plot(Nwd,plot.type='single',lwd=5,main=loc(x),ylim=c(0,100),
         xlab="",ylab=paste('number of cold days: T(2m) < ',threshold,unit(x)),
         col=c(rgb(0.5,0.5,0.7,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.8,0.5,0.8,0.5),
               rgb(0.3,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.6,0.3,0.6,0.5)),
         ...)
    grid()
    points(nwd1,pch=19)
    lines(nwd.pre,col=rgb(0.5,0.5,0.5,0.5))
  }

  Nwd <- attrcp(x,Nwd)
  attr(Nwd,'unit') <- 'days'
  attr(Nwd,'info') <- paste('number of cold days: t2m < ',threshold)
  attr(Nwd,'observation') <- nwd1
  attr(Nwd,'nwd.pre') <- nwd.pre
  index(Nwd) <- t
  class(Nwd) <- c('nevents','zoo')
  if (verbose) print('exit coldwinterdays')
  invisible(Nwd)
}

#' @export
heatwavespells <- function(x,y=NULL,dse=NULL,it='jja',threshold=30,
                           verbose=FALSE,plot=TRUE,ylab=NULL,is=1,new=TRUE,...) {
  ## Use the downscaled temperatures from ensembles to estimate the
  ## mean length og heatwaves
  ## Or use Warm Spell Duration Index (WSDI)?
  ## http://www.clipc.eu/media/clipc/org/documents/Deliverables/clipc_del7%201_final.pdf
  
  if (verbose) print('heatwaves')
  stopifnot(inherits(x,'station'))
  if (is.null(y)) y <- x
  ## Annual number of consequtive warm days
  ncwd <- aggregate(subset(spell(y,threshold=threshold),is=is),
                    year,FUN='mean')
  
  if (is.null(dse)) dse <-  DSensemble.t2m(x,biascorrect=TRUE,
                                           verbose=verbose,plot=plot)
  ## Warm season of dse
  wdse <- subset(dse,it=it)
  ## Warm season of x
  xws <- annual(subset(x,it=it),FUN='mean',nmin=90)
  ## Same period as in x (synchronise)
  index(wdse) <- year(wdse)
  ovl <- window(wdse,start=year(start(x)),end=year(end(x)))
  ## Bias correction: ensure same mean
  wdse <- wdse - mean(coredata(ovl),na.rm=TRUE) +
                       mean(coredata(xws),na.rm=TRUE)
  ## Predictor data
  irm <- setdiff(index(xws),index(ncwd))
  cal <- data.frame(y=coredata(ncwd)[!is.element(index(ncwd),irm)],
                    x=coredata(xws)[!is.element(index(xws),irm)])
  obs <- data.frame(x=coredata(xws))
  model <- lm(y ~ x, data=cal)
  
  q1 <- data.frame(x=apply(coredata(wdse),1,quantile,probs=0.05,na.rm=TRUE))
  q2 <- data.frame(x=apply(coredata(wdse),1,quantile,probs=0.95,na.rm=TRUE))
  qm <- data.frame(x=apply(coredata(wdse),1,mean,na.rm=TRUE))
  
  t <- year(index(wdse))
  preq1 <- predict(model,newdata=q1)
  tr1 <- predict(lm(preq1 ~ t + I(t^2) + I(t^3)))
  preq2 <- predict(model,newdata=q2)
  tr2 <- predict(lm(preq2 ~ t + I(t^2) + I(t^3)))
  prem  <- predict(model,newdata=qm)
  tr3 <- predict(lm(prem ~ t + I(t^2) + I(t^3)))
  Nwd <- zoo(cbind(preq1,preq2,prem,tr1,tr2,tr3),order.by=t)
  nwd.pre <- zoo(predict(model,newdata=obs),order.by=year(xws))

  if (is.null(ylab)) ylab <- paste('mean spell duration in days: ',varid(x),
                                   '> ',threshold,unit(x))
  if (plot) {
    if (new) dev.new()
    qqnorm(ncwd); qqline(ncwd)
    
    dev.new()
    par(bty='n')
    plot(zoo(wdse,order.by=year(wdse)),
         plot.type='single',col=rgb(0.5,0.5,0.5,0.2),
         ylab=expression(paste('mean temperature - ',toupper(it),(degree*C))),
         xlab='',main=loc(x))
    points(xws,pch=19)
    grid()
    
    if (new) dev.new()
    par(bty='n')
    plot(Nwd,plot.type='single',lwd=5,main=loc(x),
         xlab="",ylab=ylab,
         col=c(rgb(0.5,0.5,0.7,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.8,0.5,0.8,0.5),
               rgb(0.3,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.6,0.3,0.6,0.5)),
         ...)
    grid()
    points(ncwd,pch=19)
    lines(nwd.pre,col=rgb(0.5,0.5,0.5,0.5))
  }

  Nwd <- attrcp(x,Nwd)
  attr(Nwd,'unit') <- 'days'
  attr(Nwd,'info') <- ylab
  attr(Nwd,'observation') <- ncwd
  attr(Nwd,'nwd.pre') <- nwd.pre
  index(Nwd) <- t
  class(Nwd) <- c('nevents','zoo')
  return(Nwd)
}

#' @export
coldspells <- function(x,y=NULL,dse=NULL,it='djf',threshold=0,
                       verbose=FALSE,plot=TRUE,new=TRUE,...) {

  ylab <- paste('mean spell duration in days: ',varid(x),
                            '< ',threshold,unit(x))
  y <- heatwavespells(x,y=y,dse=dse,it=it,threshold=threshold,
                      verbose=verbose,plot=plot,ylab=ylab,is=2,new=new,...)

  invisible(y)
}

#' @export
nwetdays <- function(x,y=NULL,dse=NULL,threshold=10,
                     verbose=FALSE,plot=TRUE,new=TRUE) {
  if (is.null(y)) y <- x
  nw <- annual(y,FUN='count',threshold = threshold)
  mu <- annual(x,FUN='wetmean')
  cal <- data.frame(x=mu,y=nw)
  dfit <- glm(y ~ x,family='poisson',data=cal)

  if (plot) {
    if (new) dev.new()
    par(bty='n')
    plot(cal,pch=19,
         xlab=paste(varid(mu),unit(mu)),
         ylab='count (events/year)',main=loc(x))
    pre <- data.frame(x=seq(min(cal$x,na.rm=TRUE)-1,
                        max(cal$x,na.rm=TRUE)+5,by=0.1))
    lines(pre$x,exp(predict(dfit,newdata=pre)),col=rgb(1,0,0,0.3),lwd=3)
  }

   if (is.null(dse)) dse <-  DSensemble.precip(x,biascorrect=TRUE,
                                               predictor='air.mon.mean.nc',
                                               pattern="tas_Amon_ens_",
                                               verbose=verbose,plot=plot)

  dse <- subset(dse,it=0)
  index(dse) <- year(dse)
  ovl <- window(dse,start=year(start(x)),end=year(end(x)))
  dse <- dse - mean(coredata(ovl),na.rm=TRUE) +
               mean(coredata(mu),na.rm=TRUE)

  q1 <- data.frame(x=apply(coredata(dse),1,quantile,probs=0.05,na.rm=TRUE))
  q2 <- data.frame(x=apply(coredata(dse),1,quantile,probs=0.95,na.rm=TRUE))
  qm <- data.frame(x=apply(coredata(dse),1,mean,na.rm=TRUE))
  obs <- data.frame(x=coredata(mu))

  t <- year(index(dse))
  preq1 <- exp(predict(dfit,newdata=q1))
  preq1[preq1 > 360] <- NA  # maximum length of the summer season
  tr1 <- predict(lm(preq1 ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
  tr1[!is.finite(preq1)] <- NA
  preq2 <- exp(predict(dfit,newdata=q2))
  preq2[preq2 > 360] <- NA  # maximum length of the summer season
  tr2 <- predict(lm(preq2 ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
  tr2[!is.finite(preq2)] <- NA
  prem  <- exp(predict(dfit,newdata=qm))
  prem[prem > 360] <- NA  # maximum length of the summer season 
  tr3 <- predict(lm(prem ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
  tr3[!is.finite(prem)] <- NA
  Nwd <- zoo(cbind(preq1,preq2,prem,tr1,tr2,tr3),order.by=t)
  nwd.pre <- zoo(exp(predict(dfit,newdata=obs)),order.by=year(mu))

  if (plot) {
    if (new) dev.new()
    par(bty='n')
    plot(zoo(dse,order.by=year(dse)),
         plot.type='single',col=rgb(0.5,0.5,0.5,0.2),
         ylab=paste('annual mean',varid(x),' (',unit(x),')',sep=''),
         xlab='',main=loc(x))
    points(mu,pch=19)
    grid()
    
    if (new) dev.new()
    par(bty='n')
    plot(Nwd,plot.type='single',lwd=5,main=loc(x),
         xlab="",ylab=paste('number of days per year with P > ',
                   threshold,'mm/day'),
         col=c(rgb(0.5,0.5,0.7,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.8,0.5,0.8,0.5),
               rgb(0.3,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.6,0.3,0.6,0.5)))
    grid()
    points(nw,pch=19)
    lines(nwd.pre,col=rgb(0.5,0.5,0.5,0.5))
  }

  Nwd <- attrcp(x,Nwd)
  attr(Nwd,'unit') <- 'days'
  attr(Nwd,'info') <- paste('number of cold days: t2m < ',threshold)
  attr(Nwd,'observation') <- nw
  attr(Nwd,'nwd.pre') <- nwd.pre
  index(Nwd) <- index(dse)
  class(Nwd) <- c('nevents','zoo')
  invisible(Nwd)
}

## Catalog of potential impact and climate change indicators for CLIPC.
## Tier-1; Tier-2; Tier-3.
#    Arctic and Baltic Sea ice extent  (Tier-1)
#    River flow change (Tier-2)
#    Bathing water quality (Tier-1)
#    100 years flood return level (Tier-2)
#    Chlorophyll-a concentration (Tier-1)
#    Water-limited crop yield (Tier-2)
# Cold days (Tier-1)
#    River flood occurrence (Tier-2)
# Cold nights (Tier-1)
#    River flow (Tier-2)
# Cold spell duration index (Tier-1)
#    Water scarcity (Tier-2)
# Consecutive dry days (Tier-1)
#    Water temperature (Tier-2)
# Consecutive wet days (Tier-1)
#    Water-limited crop productivity (Tier-2)
# Diurnal temperature range (Tier-1)
# Intensity of urban heat island with city size (Tier-2)
# Frost days Heating degree-days (Tier-1)
#    Mass balance of glaciers (Tier-2)
# Rainfall Deciles (Tier-2)
#    Sea level change (Tier-1)
#    Reconnaissance Drought Index (Tier-2)
#    Greenland ice sheet mass balance (Tier-1)
# Growing Degree Days (Tier-2)
#    Grow season length of vegetation (Tier-1)
#    Chilling Units (Tier-2)
#    Hazardous substances in marine organisms (Tier-1)
#    Climatic favorability of tree species (Tier-2)
# Heavy precipitation days (Tier-1)
#    Distribution of marine species (Tier-2)
#    Ice days (Tier-1)
#    Freshwater biodiversity and water quality (Tier-2)
#    Lake and river ice cover duration (Tier-1)
#    Growing season for agriculture (Tier-2)
#    Lake and river ice phenology (Tier-1)
#    Land-cover extension below projected sea-level (Tier-2)
#    Lake Ice extension (Tier-1)
#    Moth Phenology Index (Tier-2)
# Max 1 day precipitation (Tier-1)
#    Coastal flood damage and adaptation costs (Tier-3)
