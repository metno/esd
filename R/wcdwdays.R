# number of wet, cold, dry, or wet days
mildwinterdays <- function(x,dse=NULL,threshold=0,verbose=FALSE,plot=TRUE) {
  # Estimate
  if (verbose) print('mildwinterdays')
  stopifnot(inherits(x,'station'))
  djf <- subset(x,it='djf') 
  mam <- subset(x,it='mam')
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
         xlab="",ylab=paste('number of cold days: t2m < ',threshold),
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
