# number of wet, cold, dry, or wet days
coldwinterdays <- function(x,dse=NULL,it='djf',threshold=0,
                           verbose=FALSE,plot=TRUE,...) {
  # Estimate number of days with low temperatures
  if (verbose) print('mildwinterdays')
  stopifnot(inherits(x,'station'))
  djf <- subset(x,it=it)     # Winter
  mam <- subset(x,it='mam')  # Spring/autumn
  nwd1 <- annual(-djf,FUN='count',threshold=threshold,nmin=90)
  mwd1 <- annual(djf,FUN='mean',nmin=90)
  nwd2 <- annual(-mam,FUN='count',threshold=threshold,nmin=90)
  mwd2 <- annual(mam,FUN='mean',nmin=90)

  cal <- data.frame(x=c(coredata(mwd1),coredata(mwd2)),
                    y=c(coredata(nwd1),coredata(nwd2)))
  ## Use polymomial as two different seasons are involved and the
  ## analysis involves a interpolation more than an extrapolation
  ## since winter is expected to become more similar to spring/autumn
  dfit <- glm(y ~ x + I(x^2) + I(x^3),family='poisson',data=cal)

  if (plot) {
    dev.new()
    par(bty='n')
    plot(cal,ylim=c(0,90),xlim=c(-10,10),pch=19,
         xlab=expression(paste('mean temperature ',(degree*C))),
         ylab='number of cold days',main=loc(x))
    pre <- data.frame(x=-10:10)
    lines(pre$x,exp(predict(dfit,newdata=pre)),col=rgb(1,0,0,0.3),lwd=3)

    dev.new()
    ## Use the distribution about the mean
    djf.sd <- sd(coredata(djf),na.rm=TRUE)
    ## Check it it is normally distributed
    qqnorm(coredata(djf))
    qqline(coredata(coredata(djf)),col='red')
    grid()
  }

  if (is.null(dse)) dse <-  DSensemble.t2m(x,biascorrect=TRUE,
                                           verbose=verbose,plot=plot)
  djf.dse <- subset(dse,it='djf')
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
  tr1 <- predict(lm(preq1 ~ t + I(t^2) + I(t^3)))
  preq2 <- exp(predict(dfit,newdata=q2))
  tr2 <- predict(lm(preq2 ~ t + I(t^2) + I(t^3)))
  prem  <- exp(predict(dfit,newdata=qm))
  tr3 <- predict(lm(prem ~ t + I(t^2) + I(t^3)))
  Nwd <- zoo(cbind(preq1,preq2,prem,tr1,tr2,tr3),order.by=t)
  nwd.pre <- zoo(exp(predict(dfit,newdata=obs)),order.by=year(mwd1))

  if (plot) {
    dev.new()
    par(bty='n')
    plot(zoo(djf.dse,order.by=year(djf.dse)),
         plot.type='single',col=rgb(0.5,0.5,0.5,0.2),
         ylab=expression(paste('mean temperature - winter',(degree*C))),
         xlab='',main=loc(x))
    points(mwd1,pch=19)
    grid()
    
    dev.new()
    par(bty='n')
    plot(Nwd,plot.type='single',lwd=5,main=loc(x),ylim=c(0,100),
         xlab="",ylab=paste('number of cold days: T(2m) < ',threshold),
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
  invisible(Nwd)
}


hotsummerdays <- function(x,dse=NULL,it='jja',threshold=30,
                          verbose=FALSE,plot=TRUE,...) {
    # Estimate number of days with low temperatures
  if (verbose) print('mildwinterdays')
  stopifnot(inherits(x,'station'))
  djf <- subset(x,it=it)     # summer
  nwd1 <- annual(djf,FUN='count',threshold=threshold,nmin=90)
  mwd1 <- annual(djf,FUN='mean',nmin=90)

  cal <- data.frame(x=c(coredata(mwd1)),
                    y=c(coredata(nwd1)))
  ## Use linear fit rather than polynomial as this analysis only
  ## involves one season and to a greater degree extrapolation.
  dfit <- glm(y ~ x,family='poisson',data=cal)

  if (plot) {
    dev.new()
    par(bty='n')
    plot(cal,pch=19,ylim=c(0,90),
         xlab=expression(paste('mean temperature ',(degree*C))),
         ylab='number of hot days',main=loc(x))
    pre <- data.frame(x=seq(min(cal$x,na.rm=TRUE)-1,
                            max(cal$x,na.rm=TRUE)+5,by=0.1))
    lines(pre$x,exp(predict(dfit,newdata=pre)),col=rgb(1,0,0,0.3),lwd=3)

    dev.new()
    ## Use the distribution about the mean
    djf.sd <- sd(coredata(djf),na.rm=TRUE)
    ## Check it it is normally distributed
    qqnorm(coredata(djf))
    qqline(coredata(coredata(djf)),col='red')
    grid()
  }

  if (is.null(dse)) dse <-  DSensemble.t2m(x,biascorrect=TRUE,
                                           verbose=verbose,plot=plot)
  djf.dse <- subset(dse,it='djf')
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
  preq1[preq1 > 90] <- NA  # maximum length of the summer season
  tr1 <- predict(lm(preq1 ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
  tr1[!is.finite(preq1)] <- NA
  preq2 <- exp(predict(dfit,newdata=q2))
  preq2[preq2 > 90] <- NA  # maximum length of the summer season
  tr2 <- predict(lm(preq2 ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
  tr2[!is.finite(preq2)] <- NA
  prem  <- exp(predict(dfit,newdata=qm))
  prem[prem > 90] <- NA  # maximum length of the summer season 
  tr3 <- predict(lm(prem ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
  tr3[!is.finite(prem)] <- NA
  Nwd <- zoo(cbind(preq1,preq2,prem,tr1,tr2,tr3),order.by=t)
  nwd.pre <- zoo(exp(predict(dfit,newdata=obs)),order.by=year(mwd1))

  if (plot) {
    dev.new()
    par(bty='n')
    plot(zoo(djf.dse,order.by=year(djf.dse)),
         plot.type='single',col=rgb(0.5,0.5,0.5,0.2),
         ylab=expression(paste('mean temperature',(degree*C))),
         xlab='',main=loc(x))
    points(mwd1,pch=19)
    grid()
    
    dev.new()
    par(bty='n')
    plot(Nwd,plot.type='single',lwd=5,main=loc(x),ylim=c(0,90),
         xlab="",ylab=paste('number of hot days: T(2m) > ',threshold),
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
  invisible(Nwd)
}

heatwavespells <- function(x,dse=NULL,it='jja',threshold=30,
                           verbose=FALSE,plot=TRUE,ylab=NULL,is=1,...) {
  ## Use the 
  if (verbose) print('heatwaves')
  stopifnot(inherits(x,'station'))
  ## Annual number of consequtive warm days
  ncwd <- aggregate(subset(spell(x,threshold=threshold),is=is),
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
    dev.new()
    qqnorm(ncwd); qqline(ncwd)
    
    dev.new()
    par(bty='n')
    plot(zoo(wdse,order.by=year(wdse)),
         plot.type='single',col=rgb(0.5,0.5,0.5,0.2),
         ylab=expression(paste('mean temperature - ',toupper(it),(degree*C))),
         xlab='',main=loc(x))
    points(xws,pch=19)
    grid()
    
    dev.new()
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

coldspells <- function(x,dse=NULL,it='djf',threshold=0,
                       verbose=FALSE,plot=TRUE,...) {

  ylab <- paste('mean spell duration in days: ',varid(x),
                            '< ',threshold,unit(x))
  y <- heatwavespells(x,dse=dse,it=it,threshold=threshold,
                      verbose=verbose,plot=plot,ylab=ylab,is=2,...)

  invisible(y)
}


nwetdays <- function(x,dse=NULL,threshold=10,
                     verbose=FALSE,plot=TRUE) {
  nw <- annual(x,FUN='count',threshold = threshold)
  mu <- annual(x,FUN='wetmean')
  cal <- data.frame(x=mu,y=nw)
  dfit <- glm(y ~ x,family='poisson',data=cal)

  if (plot) {
    dev.new()
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
    dev.new()
    par(bty='n')
    plot(zoo(dse,order.by=year(dse)),
         plot.type='single',col=rgb(0.5,0.5,0.5,0.2),
         ylab=paste('annual mean',varid(x),' (',unit(x),')',sep=''),
         xlab='',main=loc(x))
    points(mu,pch=19)
    grid()
    
    dev.new()
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
  
