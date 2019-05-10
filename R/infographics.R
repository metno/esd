# Show the cumulative sum of station value from January 1st. Use
# different colours for different year.
cumugram <- function(x,it=NULL,start='-01-01',prog=FALSE,verbose=FALSE,FUN='mean',main=NULL,...) {
  stopifnot(!missing(x),inherits(x,"station"))
  
  #print("cumugram")
  yrs <- as.numeric(rownames(table(year(x))))
  today <- Sys.Date(); yesterday <- seq(today, length.out=2, by=-1)[2]
  
  #print(yrs)
  ny <- length(yrs)
  j <- 1:ny
  col <- rgb(j/ny,abs(sin(pi*j/ny)),(1-j/ny),0.3)
  class(x) <- "zoo"

  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") )
      unit <- expression(degree*C) else
      unit <- attr(x,'unit')
  titletext <- paste('Running cumulative',FUN,'of')
  if (is.null(main)) 
    eval(parse(text=paste("main <- paste('",titletext,"',
                          tolower(attr(x,'longname')),sep=' ')")))
  dev.new()
  par(bty="n")
  z <- coredata(x)
  ylim <- c(NA,NA)

  #print('Find the y-range')
  y.rest <- rep(NA,ny); y2n <- y.rest
  ylim <- max(coredata(x),na.rm=TRUE) # to avoid getting warnings with empty vectors.
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],start,sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    y.rest[i] <- mean(coredata(window(x,start=as.Date(paste(yrs[i],format(Sys.Date(),'-%m-%d'),sep='')),
                                      end=as.Date(paste(yrs[i],'-12-31',sep='')))))
    y2n[i] <- mean(coredata(window(x,end=as.Date(paste(yrs[i],format(Sys.Date()-1,'-%m-%d'),sep='')),
                                     start=as.Date(paste(yrs[i],'-01-01',sep='')))))                                  
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    if (FUN=='mean') z <- cumsum(coredata(y))/1:length(y) else
    if (FUN=='sum') z <- cumsum(coredata(y))
    ok <- is.finite(z)
    #rint(c(i,yrs[i],range(z[ok],na.rm=TRUE),ylim))
    ylim[!is.finite(ylim)] <- NA
    ylim[1] <- min(c(ylim,y[ok]),na.rm=TRUE)
    ylim[2] <- max(c(ylim,y[ok]),na.rm=TRUE)
  }
  #print(ylim)
  names(y2n) <- yrs
  y2n <- round(sort(y2n,decreasing=TRUE),2)
  
  plot(c(0,length(y)),ylim,
       type="n",xlab="",
       main=main,sub=attr(x,'location'),ylab=ylab(x),...)
  grid()

  cm <- rep(NA,ny)
  
  #browser()

  mm <- format(yesterday, "%m")
  dd <- format(yesterday, "%d")
  period <- paste('YYYY',start,' to YYYY-',paste(mm,dd,sep='-'),sep='')
  if (verbose) {print(yesterday); print(mm); print(dd); print(period)}
  
  if (verbose) print('No. year min max ylim[1] ylim[2]')
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],start,sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],start,sep='')))
    if (FUN=='mean') z <- cumsum(coredata(y))/1:length(y) else
    if (FUN=='sum') z <- cumsum(coredata(y))
    
    if (FUN=='mean') cm[i] <- mean(coredata(window(x,
                                   start=as.Date(paste(yrs[i],start,sep='')),
                                   end=as.Date(paste(yrs[i],mm,dd,sep='-'))))) else 
                     cm[i] <- sum(coredata(window(x,
                                    start=as.Date(paste(yrs[i],start,sep='')),
                                    end=as.Date(paste(yrs[i],mm,dd,sep='-')))))
    lines(t,z,lwd=2,col=col[i])
    if (verbose) print(c(i,yrs[i],cm[i],range(z[ok],na.rm=TRUE),ylim))
  }
  if (is.null(it)) {
    lines(t,z,lwd=5,col="black")
    lines(t,z,lwd=2,col=col[i])
  } else {
    y <- window(x,start=as.Date(paste(it,start,sep='')),
                    end=as.Date(paste(it,'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(it,start,sep='')))
    if (FUN=='mean') z <- cumsum(coredata(y))/1:length(y)  else
    if (FUN=='sum') z <- cumsum(coredata(y))  
    lines(t,z,lwd=5,col="black")
    lines(t,z,lwd=2,col=col[i])
  }
  tn <- t[length(t)]; 

  ## Here is some difference between monthly and daily data:
  if (!is.na(coredata(z[length(z)]))) zn <- coredata(z[length(z)]) else
                                      zn <- coredata(z[length(z)-1])
  n <- max(table(year(x)))
  if (n>=365) n <- as.numeric(diff(as.Date(c(paste(yrs[ny],start,sep=''),
                                             paste(yrs[ny],'-12-31',sep=''))))+1)
  if (n>=365) tm <- julian(as.Date('1900-12-31')) - julian(as.Date('1900-01-01')) else
              tm <- julian(as.Date('1900-12-01')) - julian(as.Date('1900-01-01'))
  #browser()
  zp <- length(z)/n * zn + (n-length(z))/n * quantile(y.rest,0.95,na.rm=TRUE)
  zm <- length(z)/n * zn + (n-length(z))/n * quantile(y.rest,0.05,na.rm=TRUE)
  zz <- length(z)/n * zn + (n-length(z))/n * mean(y.rest,na.rm=TRUE)
  if (prog) {
    polygon(c(tn,rep(tm,2),tn),c(zn,zp,zm,zn),
            col=rgb(0.5,0.5,0.5,0.1),border=rgb(0.5,0.5,0.5,0.2),lwd=2)
    lines(c(tn,tm),c(zn,zz),col=rgb(0.3,0.3,0.3,0.1),lwd=3)
    text(tm,zp,round(zp,1),pos=3,cex=0.5,col='grey40')
    text(tm,zm,round(zm,1),pos=1,cex=0.5,col='grey40')
    text(tm,zz,round(zz,1),pos=4,cex=0.75)
    print(paste('Prognosis for end-of-year: ',round(zz,1),' (',round(zm,1),',',round(zp,1),')',sep=''))
  }

  if (!is.precip(x))
    par(new=TRUE,fig=c(0.70,0.85,0.20,0.35),mar=c(0,3,0,0),
        cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  else
    par(new=TRUE,fig=c(0.70,0.85,0.70,0.85),mar=c(0,3,0,0),
        cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:ny,1:ny)
  image(1:2,yrs,colbar,col=col)

  srt <- order(cm,decreasing=TRUE)
  if (verbose) print(y2n)
  result <- cbind(yrs[srt],cm[srt])
  if (verbose) print(round(t(result)))
  colnames(result) <- c('year','cumulated')
  attr(result,'period')  <- period
  invisible(result)
  
}

# Estimate how the variance varies with season 
# sd from inter-annual variability of daily values
climvar <- function(x,FUN='sd',plot=TRUE,...) {
  yrs <- as.numeric(rownames(table(year(x))))
  #print(yrs)
  ny <- length(yrs)
  X <- x; class(X) <- "zoo"
  
  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") ) {
    unit <- expression(degree*C) 
  } else {
    unit <- attr(x,'unit')
  }
  eval(parse(text=paste("main <- expression(paste('seasonal ",
#               deparse(substitute(FUN))," of ',",
               FUN," of ',",attr(x,'variable'),"))")))
  Z <- matrix(rep(NA,ny*365),ny,365)
  
  for (i in 1:ny) {
    y <- window(X,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    Z[i,] <- approx(t,y,1:365)$y
  }
  z <- apply(Z,2,FUN,na.rm=TRUE,...)
  wt <- 2*pi*(1:365)/365
  s1 <- sin(wt); c1 <- cos(wt); s2 <- sin(2*wt); c2 <- cos(2*wt)
  s3 <- sin(3*wt); c3 <- cos(3*wt); s4 <- sin(4*wt); c4 <- cos(4*wt)
  acfit <- predict(lm(z ~s1 + c1 + s2 + c2 + s3 + c3 + c4 + s4))
    
  if (plot) {
    dev.new()
    par(bty="n")
    plot(c(0,365),range(z,na.rm=TRUE),
         type="n",xlab="",
         main=main,
        sub=attr(x,'location'),ylab=ylab(x))
    grid()
    lines(z,lwd=5)
    lines(z,lwd=3,col="grey")
    lines(acfit,lwd=5)
    lines(acfit,lwd=3,col="red")

    par(new=TRUE,fig=c(0.15,0.35,0.70,0.90),mar=c(0,0,0,0),
        yaxt="n",xaxt="n",las=1)
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    legend(0,1,c("raw data","harmonic fit"),lwd=3,col=c("grey","red"),bty="n",cex=0.6)  
  }
  
  acfit <- attrcp(x,acfit)
  attr(acfit,'raw_data') <- z
  attr(acfit,'history') <- history.stamp(x)
  return(z)
}



histwet <- function(x,breaks=NULL,threshold=1) {
  if (is.null(breaks)) breaks=seq(0,1.1*max(x,na.rm=TRUE),by=5)
  h <- hist(x[x > threshold],breaks=breaks,plot=FALSE)
  return(h$density)
}

## Some additional infographics: show p.d.f. for each year    
visprob <- function(x,...) UseMethod("visprob")

visprob.default <- function(x,...) {

}

visprob.station <- function(x,...,y=NULL,is=1,dy=0.01,verbose=FALSE) {
  if (is.precip(x)) visprob.station.precip(x,y=y,is=is,
                                           dy=dy,verbose=verbose,...) 
}

visprob.station.precip <- function(x,...,y=NULL,is=1,threshold=1,dy=0.005,
                                   breaks=NULL,pdf=FALSE,verbose=FALSE) {
## Plot the histogram for each year in different colours, depending on y. Iy y
## is NULL, use the year to set the colour
  if (verbose) print('visprob.station.precip')
  ## If y is provided, synchronise the two time series:
  if (!is.null(y)) {
    y <- subset(y,it=c(start(x),end(x)))
    x <- subset(x,it=c(start(y),end(y)))
    y <- annual(y)
  } else y <- year(annual(x))
  mu <- aggregate(x,year,FUN='wetmean',threshold=threshold)
  fw <- aggregate(x,year,FUN='wetfreq',threshold=threshold)
  col <- colscal(n=length(y),alpha=coredata(fw))
  srtc <- order(y)
  col <- col[srtc]
  if (is.null(breaks))
    breaks <- seq(floor(min(x,na.rm=TRUE)),ceiling(max(x,na.rm=TRUE))+5,by=5)
  z <- aggregate(x,year,FUN='histwet',breaks=breaks,threshold=threshold)
  if (verbose) print(c(dim(z),length(y)))
  dy <- abs(max(z,na.rm=TRUE)*dy)
  mids <- 0.5*(breaks[-1] + breaks[-length(breaks)])
  par(bty='n',yaxt='n')
  plot(range(breaks),c(0,max(z,na.rm=TRUE) + length(y)*dy),type='n',
  ylab='f(x)',xlab=paste(varid(x),unit(x)),
       main=paste('Statistical distribution for',loc(x)),...)
  for (i in seq(length(y),1,by=-1)) {
    lines(mids,z[i,]+dy*i,col="grey",lwd=5)
    lines(mids,z[i,]+dy*i,col=col[i],lwd=4)
  }
  if (pdf) {
    if (verbose) print('add pdfs')
    for (i in 1:length(y)) {
      lines(mids,dy*i + exp(-mids/coredata(mu[i]))/coredata(mu[i]),
            col='black',lty=2)
    }
  }
  if (!is.null(loc(y)))
    text(par()$xaxp[2],par()$yaxp[2],loc(y),pos=2)
}

graph <- function(x,...) UseMethod("graph")

graph.default <- function(x,img=NULL,pch='fancy',it=NULL,col=rgb(0.5,0.5,0.5,0.5),lwd=5,
                          xlim=NULL,ylim=NULL,new=TRUE,col.obs='black',...) {
    print('graph.default')
    ## Produce the graphics:
    
    if (new) dev.new()
    if (!is.null(img)) {
      par0 <- par()
      par(mar=rep(0,4))
      plot(c(0,1),c(0,1),type='n')
      rasterImage(img, -0.05, -0.05, 1.05, 1.05)
      par(new=TRUE,col.axis='white',col.lab='white',xaxt='n',yaxt='n',
          mar=par0$mar,bty='n',col.sub='white')
    }
    if (!is.null(it)) y <- subset(x,it=it) else y <- x
    plot.zoo(y,lwd=lwd,col=col,ylim=ylim,xlim=xlim,
             ylab=ylab(y),sub=loc(y))
    #balls(y)
    if (!is.null(pch)) if (pch=='fancy') balls(y,col=col.obs) else points(zoo(y),pch=pch,col=col.obs)
    par(xaxt='s',yaxt='s')
    axis(1,col='white')
    axis(2,col='white')
}

graph.dsensemble <- function(x,img=NULL,pch='fancy',it=NULL,col=rgb(1,0.7,0.7,0.1),
                             lwd=5,xlim=NULL,ylim=NULL,add=FALSE,new=TRUE,ensmean=FALSE,col.obs='black') {
    #print('graph.dsensemble')
    ## Produce the graphics:
    if ((!add) & (new)) dev.new()
    if (!is.null(img)) {
      par0 <- par()
      par(mar=rep(0,4))
      plot(c(0,1),c(0,1),type='n')
      rasterImage(img, -0.05, -0.05, 1.05, 1.05)
      par(new=TRUE,col.axis='white',col.lab='white',xaxt='n',yaxt='n',
          mar=par0$mar,bty='n',col.sub='white')
    }
    if (!is.null(it)) y <- subset(x,it=it) else y <- x
    index(y) <- year(y)
    index(attr(y,'station')) <- year(attr(y,'station'))
    if (is.null(xlim)) xlim <- range(index(y))
    if (is.null(ylim)) ylim <- range(coredata(y),na.rm=TRUE)
    
    if (!add) plot.zoo(attr(y,'station'),lwd=lwd,col=rgb(0.5,0.5,0.5,0.5),
                       ylim=ylim,xlim=xlim,ylab=ylab(attr(y,'station')),
                       sub=loc(x),plot.type='single',xlab='')
    for (i in 1:dim(x)[2]) lines(y[,i],lwd=7,col=col)
    if (ensmean) lines(index(x),apply(coredata(x),1,'mean',na.rm=TRUE),
                                 lwd=3,col='red')

    #balls(attr(y,'station'))
    if (!is.null(pch)) if (pch=='fancy') balls(attr(y,'station'),col=col.obs) else points(zoo(attr(y,'station')),pch=pch,col=col.obs)
    par(xaxt='s',yaxt='s')
    if (!is.null(img)) col.axis <- 'white' else col.axis <- 'black'
    axis(1,col=col.axis)
    axis(2,col=col.axis)
}

graph.list <- function(x,img=NULL,pch='fancy',it=NULL,
                       col=c(rgb(1,1,0.5,0.05),rgb(1,0.5,0.5,0.05),rgb(0.5,1,0.5,0.05),
                             rgb(0.5,0.5,0.5,0.05) ),
                       lwd=5,xlim=NULL,ylim=NULL,add=FALSE,new=TRUE,ensmean=FALSE,col.obs='black') {
  if ((!is.null(it)) & (inherits(x[[1]],'dsensemble')))
    y <- subset(x[[1]],it=it) else y <- x[[1]]
    index(y) <- year(y)
  graph(y,img=img,pch=pch,col=col[1],lwd=lwd,xlim=xlim,ylim=ylim,add=add,new=new,col.obs=col.obs)
  if (!is.null(attr(x,'obs')) & is.null(attr(y,'dsensemble'))) obs <- attr(x,'obs') else
                                                               obs <- attr(y,'station')
  index(obs) <- year(obs)
  
  for (j in c(2:length(x),1)) {
    if ((!is.null(it)) & (inherits(x[[j]],'dsensemble')))
         y <- subset(x[[j]],it=it) else y <- x[[j]]
         index(y) <- year(y)
    for (i in 1:dim(y)[2]) lines(y[,i],lwd=7,col=col[j])
  }
  lines(obs,lwd=3,col=rgb(0.5,0.5,0.5,0.25))

  if (ensmean) {
    emcol <- c('wheat','red','green','grey')    
    for (i in 1:length(x)) lines(index(x[[i]]),apply(coredata(x[[i]]),1,'mean',na.rm=TRUE),
                                 lwd=3,col=emcol[i])
    legend(index(y)[1],max(coredata(y)),names(x),col=emcol,lty=1,lwd=3,bty='n')
  }
  if (!is.null(pch)) if (pch=='fancy') balls(obs,col=col.obs) else points(obs,pch=pch,col=col.obs)
}


graph.zoo <- function(x,img=NULL,it=NULL,col=rgb(1,0.7,0.7,0.1),pch=1,
                      lwd=5,xlim=NULL,ylim=NULL,xlab='',ylab='',add=FALSE,
                      new=TRUE,ensmean=FALSE,col.obs='black') {
  #print('graph.zoo')
    ## Produce the graphics:
    if ((!add) & (new)) dev.new()
    if (!is.null(img)) {
      par0 <- par()
      par(mar=rep(0,4))
      plot(c(0,1),c(0,1),type='n')
      rasterImage(img, -0.05, -0.05, 1.05, 1.05)
      par(new=TRUE,col.axis='white',col.lab='white',xaxt='n',yaxt='n',
          mar=par0$mar,bty='n',col.sub='white')
    }
    if (!is.null(it)) y <- subset(x,it=it) else y <- x
    index(y) <- year(y)
    if (is.null(xlim)) xlim <- range(index(y))
    if (is.null(ylim)) ylim <- range(coredata(y),na.rm=TRUE)
    if (!add) plot.zoo(y[,1],lwd=lwd,col=col,ylim=ylim,xlim=xlim,
                       ylab=ylab,xlab=xlab,plot.type='single')
    grid()
    for (i in 1:dim(x)[2]) lines(y[,i],lwd=7,col=col)

    if (!is.null(pch)) {
      if (pch=='fancy') {
        balls(attr(y,'station'),col=col.obs) 
      } else {
        points(zoo(attr(y,'station')),pch=pch,col=col.obs)
      }
    }
    par(xaxt='s',yaxt='s')
    if (!is.null(img)) col.axis <- 'white' else col.axis <- 'black'
    axis(1,col=col.axis)
    axis(2,col=col.axis)
}

qp.test <- function(x,...) UseMethod("qp.test")

qp.test.station <- function(x,...) {
  if (is.precip(x)) qp.test.precip(x,...) else
  if (is.T(x)) qp.test.t2m(x,...)
}


qp.test.precip <- function(x,p=c(seq(0.1,0.95,0.05),0.97,0.98,0.99),threshold=1,...) {
  # From qqplotter:
  x[x < threshold] <- NA
  if (is.null(dim(x))) {
    qp <- quantile(x,prob=p,na.rm=TRUE)
    qmu <- -log(1-p)*mean(x,na.rm=TRUE)
  } else {
    qp <- apply(coredata(x),2,quantile,prob=p,na.rm=TRUE)
    qmu <- -log(1-p)%o%apply(coredata(x),2,wetmean,na.rm=TRUE)
    #fw <- round(100*apply(coredata(x),2,wetfreq))
  }
  plot(qp,qmu,pch=19,col=rgb(0,0,1,0.2),
       xlab=expression(q[p]),ylab=expression(-log(1-p)*mu))
  lines(range(qp,qmu),range(qp,qmu))
  grid()
}

qp.test.t2m <- function(x,p=c(0.01,0.02,0.03,0.04,seq(0.1,0.95,0.05),
                                0.97,0.98,0.99),...) {
  d <- dim(x); if (is.null(d)) d <- c(length(x),1)
  if (is.null(dim(x))) {
    qp <- quantile(x,prob=p,na.rm=TRUE)
    qmu <-  qnorm(p=p,mean=mean(coredata(x),na.rm=TRUE),sd=sd(coredata(x),na.rm=TRUE))
  } else {
    qp <- apply(coredata(x),2,quantile,prob=p,na.rm=TRUE)
    qmu <- qp + NA
    for (i in 1:length(p))
      qmu[i,] <- qnorm(p=p[i],mean=apply(coredata(x),2,mean,na.rm=TRUE),
                       sd=apply(coredata(x),2,sd,na.rm=TRUE))
  }
  plot(qp,qmu,pch=19,col=rgb(1,0,0,0.2),
       xlab=expression(q[p]),ylab='qnorm(p,mean(x),sd(x))')
  lines(range(qp,qmu),range(qp,qmu))
  grid() 
}
