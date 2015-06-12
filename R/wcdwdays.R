# number of wet, cold, dry, or wet days
coldwinterdays <- function(x,dse=NULL,it='djf',threshold=0,verbose=FALSE,plot=TRUE) {
  # Estimate
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

  if (plot) {
    dev.new()
    par(bty='n')
    plot(cal,ylim=c(0,90),xlim=c(-10,10),pch=19,
         xlab=expression(paste('mean temperature ',(degree*C))),
         ylab='number of cold days',main=loc(x))
    dfit <- glm(y ~ x + I(x^2) + I(x^3),family='poisson',data=cal)
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
  ovl <- window(djf.dse,start=start(y),end=end(y))
  djf.dse <- djf.dse - mean(coredata(ovl),na.rm=TRUE) +
                       mean(coredata(mwd1),na.rm=TRUE)

  q1 <- data.frame(x=apply(coredata(djf.dse),1,quantile,probs=0.05))
  q2 <- data.frame(x=apply(coredata(djf.dse),1,quantile,probs=0.95))
  qm <- data.frame(x=apply(coredata(djf.dse),1,mean))
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
               rgb(0.3,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.6,0.3,0.6,0.5)))
    grid()
    points(nwd1,pch=19)
    lines(nwd.pre,col=rgb(0.5,0.5,0.5,0.5))
  }

  Nwd <- attrcp(x,Nwd)
  attr(Nwd,'unit') <- 'days'
  attr(Nwd,'info') <- paste('number of cold days: t2m < ',threshold)
  invisible(Nwd)
}

heatwaves <- function(x,dse=NULL,it='jja',threshold=30,verbose=FALSE,plot=TRUE) {
  ## Use the 
  if (verbose) print('heatwaves')
  stopifnot(inherits(x,'station'))
  ## Annual number of consequtive warm days
  ncwd <- annual(subset(spell(x,threshold=treshold),is=1))
  
  if (is.null(dse)) dse <-  DSensemble.t2m(x,biascorrect=TRUE,
                                           verbose=verbose,plot=plot)
  ## Warm season of dse
  wdse <- subset(dse,it=it)
  ## Warm season of x
  xws <- annual(subset(x,it=it),FUN='mean',nmin=90)
  ## Same period as in x (synchronise)
  ovl <- window(wdse,start=start(y),end=end(y))
  ## Bias correction: ensure same mean
  wdse <- wdse - mean(coredata(ovl),na.rm=TRUE) +
                       mean(coredata(xws),na.rm=TRUE)
  ## Predictor data
  cal <- data.frame(y=coredata(ncwd), x=coredata(xws))
  model <- lm(y ~ x, data=cal)
  
  q1 <- data.frame(x=apply(coredata(wdse),1,quantile,probs=0.05))
  q2 <- data.frame(x=apply(coredata(wdse),1,quantile,probs=0.95))
  qm <- data.frame(x=apply(coredata(wdse),1,mean))
  
  t <- year(index(wdse))
  preq1 <- predict(model,newdata=q1)
  tr1 <- lm(preq1 ~ t + I(t^2) + I(t^3))
  preq2 <- predict(model,newdata=q2)
  tr2 <- lm(preq2 ~ t + I(t^2) + I(t^3))
  prem  <- predict(model,newdata=qm)
  tr3 <- lm(prem ~ t + I(t^2) + I(t^3))
  Nwd <- zoo(cbind(preq1,preq2,prem,tr1,tr2,tr3),order.by=t)

  if (plot) {
    dev.new()
    par(bty='n')
    plot(zoo(wdse,order.by=year(wdse)),
         plot.type='single',col=rgb(0.5,0.5,0.5,0.2),
         ylab=expression(paste('mean temperature - ',toupper(it),(degree*C))),
         xlab='',main=loc(x))
    points(mwd1,pch=19)
    grid()
    
    dev.new()
    par(bty='n')
    plot(Nwd,plot.type='single',lwd=5,main=loc(x),ylim=c(0,100),
         xlab="",ylab=paste('number of hot days: T(2m) > ',threshold),
         col=c(rgb(0.5,0.5,0.7,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.8,0.5,0.8,0.5),
               rgb(0.3,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.6,0.3,0.6,0.5)))
    grid()
    points(nwd1,pch=19)
    lines(nwd.pre,col=rgb(0.5,0.5,0.5,0.5))
  }

  Nwd <- attrcp(x,Nwd)
  attr(Nwd,'unit') <- 'days'
  attr(Nwd,'info')
  return(Nwd)
}
