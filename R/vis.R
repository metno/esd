# Visualise - different type of plotting... 'Infographics' type.

vis <- function(x,...) UseMethod("vis")

vis.station <- function(x,...) {
  if (is.precip(x)) vis.station.precip(x,...) else
  if (is.T(x)) vis.station.t2m(x,...)
}

vis.station.precip <- function(x,p=c(seq(0.1,0.95,0.05),0.97,0.98,0.99),...) {
  # From qqplotter:
  qp <- apply(coredata(x),2,quantile,prob=p,na.rm=TRUE)
  qmu <- -log(1-p)%o%apply(coredata(x),2,wetmean)
  fw <- round(100*apply(coredata(x),2,wetfreq))
  plot(qp,qmu,pch=19,col=rgb(0,0,1,0.2))
  lines(range(qp,qmu),range(qp,qmu))
  grid()
}

vis.station.t2m <- function(x,p=c(0.01,0.02,0.03,0.04,seq(0.1,0.95,0.05),
                                0.97,0.98,0.99),...) {
  qp <- apply(coredata(x),2,quantile,prob=p,na.rm=TRUE)
  qmu <- qp + NA
  for (i in 1:dim(qmu)[1])
    qmu[i,] <- qnorm(p=p[i],mean=apply(coredata(x),2,mean,na.rm=TRUE),
                     sd=apply(coredata(x),2,sd,na.rm=TRUE))
  plot(qp,qmu,pch=19,col=rgb(1,0,0,0.2))
  lines(range(qp,qmu),range(qp,qmu))
  grid() 
}

vis.field <- function(x,...) {
}

vis.eof <- function(x,...) {
}

vis.spell <- function(x,...) {
}

vis.cca <- function(x,...) {
}

vis.mvr <- function(x,...) {
}

vis.dsensemble <- function(x,...) {
}

vis.ds <- function(x,...) {
}


diagram <- function(x,...) UseMethod("diagram")

diagram.dsensemble <- function(x,it=0,...) {
  stopifnot(inherits(x,'dsensemble'))
  #print("subset") 
  if (!inherits(attr(x,'station'),'annual')) z <- subset(x,it=it) else
                                             z <- x
  y <- attr(z,'station')
  #print("diagnose")
  pscl <- c(0.9,1.1)
  if (max(coredata(z),na.rm=TRUE) < 0) pscl <- rev(pscl)
  #print("...")
  plot(y,type="b",pch=19,
       xlim=range(year(z)),
       ylim=pscl*range(coredata(z),na.rm=TRUE))
  grid()
  usr <- par()$usr; mar <- par()$mar; fig <- par()$fig
  t <- year(z); n <- dim(z)[2]
  col <- rgb(seq(1,0,length=n)^2,sin(seq(0,pi,length=n))^2,seq(0,1,length=n)^2,0.2)
  for (i in 1:n) lines(t,z[,i],col=col[i],lwd=2)
  points(y,pch=19,lty=1)
}

diagram.ds <- function(x,...) {
}


colscal <- function(n=30,col="bwr",test=FALSE) {

  test.col <- function(r,g,b) {
    dev.new()
    par(bty="n")
    plot(r,col="red")
    points(b,col="blue")
    points(g,col="green")
  }
  
  # Set up colour-palette
  x <- 1:n
  r0 <- round(n*0.55)
  g0 <- round(n*0.5)
  b0 <- round(n*0.45)
  s <- -0.1/n
  if (n < 30) sg <- s*2.5 else sg <- s
  n1 <- g0; n2 <- n-n1

  if (col=="bwr") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
    g <- exp(sg*(x - g0)^2)^2
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
    col <- rgb(r,g,b)
  } else if (col=="rwb") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
    g <- exp(sg*(x - g0)^2)^2
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
    col <- rgb(b,g,r)
  } else if (col=="faint.bwr") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
    g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
    col <- rgb(r,g,b)
  } else if (col=="faint.rwb") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
    g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
    col <- rgb(b,g,r)
  }

  if (test) test.col(r,g,b)
  return(col)
}

colbar <- function(scale,col,fig=c(0.15,0.2,0.15,0.3)) {
  par(xaxt="n",yaxt="s",fig=fig,mar=c(0,1,0,0),new=TRUE,las=1,cex.axis=0.6)
  image(1:2,scale,rbind(scale,scale),col=col)
}


# Show the cumulative sum of station value from January 1st. Use
# different colours for different year.
cumugram <- function(x,it=NULL,prog=FALSE,verbose=FALSE,...) {
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
  eval(parse(text=paste("main <- expression(paste('Running cumulative mean of ',",
               attr(x,'variable'),"))")))
  dev.new()
  par(bty="n")
  z <- coredata(x)
  ylim <- c(NA,NA)

  #print('Find the y-range')
  y.rest <- rep(NA,ny); y2n <- y.rest
  ylim <- max(coredata(x),na.rm=TRUE) # to avoid getting warnings with empty vectors.
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    y.rest[i] <- mean(coredata(window(x,start=as.Date(paste(yrs[i],format(Sys.Date(),'-%m-%d'),sep='')),
                                      end=as.Date(paste(yrs[i],'-12-31',sep='')))))
    y2n[i] <- mean(coredata(window(x,end=as.Date(paste(yrs[i],format(Sys.Date()-1,'-%m-%d'),sep='')),
                                     start=as.Date(paste(yrs[i],'-01-01',sep='')))))                                  
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    z <- cumsum(coredata(y))/1:length(y)
    ok <- is.finite(z)
    #rint(c(i,yrs[i],range(z[ok],na.rm=TRUE),ylim))
    ylim[!is.finite(ylim)] <- NA
    ylim[1] <- min(c(ylim,z[ok]),na.rm=TRUE)
    ylim[2] <- max(c(ylim,z[ok]),na.rm=TRUE)
  }
  #print(ylim)
  names(y2n) <- yrs
  y2n <- round(sort(y2n,decreasing=TRUE),2)
  
  plot(c(0,365),ylim,
       type="n",xlab="",
       main=main,sub=attr(x,'location'),ylab=unit,...)
  grid()

  cm <- rep(NA,ny)
  
  #browser()
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    z <- cumsum(coredata(y))/1:length(y)

    mm <- format(yesterday, "%m")
    dd <- as.numeric(yesterday, "%d")
    
    cm[i] <- mean(coredata(window(x,
            start=as.Date(paste(yrs[i],'-01-01',sep='')),
            end=as.Date(paste(yrs[i],mm,dd,sep='-')))))
    lines(t,z,lwd=2,col=col[i])
    print(c(i,yrs[i],range(z[ok],na.rm=TRUE),ylim))
  }
  if (is.null(it)) {
    lines(t,z,lwd=5,col="black")
    lines(t,z,lwd=2,col=col[i])
  } else {
    y <- window(x,start=as.Date(paste(it,'-01-01',sep='')),
                    end=as.Date(paste(it,'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(it,'-01-01',sep='')))
    z <- cumsum(coredata(y))/1:length(y)   
    lines(t,z,lwd=5,col="black")
    lines(t,z,lwd=2,col=col[i])
  }
  tn <- t[length(t)]; 
  tm <- julian(as.Date('1900-12-31')) - julian(as.Date('1900-01-01'))
  zn <- coredata(z[length(z)-1])
  n <- 365
  #browser()
  zp <- length(z)/n * zn + (1-length(z)/n) * quantile(y.rest,0.95,na.rm=TRUE)
  zm <- length(z)/n * zn + (1-length(z)/n) * quantile(y.rest,0.05,na.rm=TRUE)
  zz <- length(z)/n * zn + (1-length(z)/n) * mean(y.rest,na.rm=TRUE)
  if (prog) {
    polygon(c(tn,rep(tm,2),tn),c(zn,zp,zm,zn),
            col=rgb(0.5,0.5,0.5,0.1),border=rgb(0.5,0.5,0.5,0.2),lwd=2)
    lines(c(tn,tm),c(zn,zz),col=rgb(0.3,0.3,0.3,0.1),lwd=3)
    text(tm,zp,round(zp,1),pos=3,cex=0.5,col='grey40')
    text(tm,zm,round(zm,1),pos=1,cex=0.5,col='grey40')
    text(tm,zz,round(zz,1),pos=4,cex=0.75)
    print(paste('Prognosis for end-of-year: ',round(zz,1),' (',round(zm,1),',',round(zp,1),')',sep=''))
  }

  if (varid(x)!='precip') 
    par(new=TRUE,fig=c(0.70,0.85,0.20,0.35),mar=c(0,3,0,0),
        cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  else
    par(new=TRUE,fig=c(0.70,0.85,0.70,0.85),mar=c(0,3,0,0),
        cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:ny,1:ny)
  image(1:2,yrs,colbar,col=col)

  srt <- order(cm,decreasing=TRUE)
  invisible(cbind(yrs[srt],cm[srt]))
  if (verbose) print(y2n)
  invisible(y2n)
}

# Estimate how the variance varies with season 
# sd from inter-annual variability of daily values

climvar <- function(x,FUN='sd',plot=TRUE,...) {
  yrs <- as.numeric(rownames(table(year(x))))
  #print(yrs)
  ny <- length(yrs)
  X <- x; class(X) <- "zoo"
  
  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") )
      unit <- expression(degree*C) else
      unit <- attr(x,'unit')
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
        sub=attr(x,'location'),ylab=unit)
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


