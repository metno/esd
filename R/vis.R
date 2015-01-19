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

# Show the temperatures against the day of the year. Use
# different colours for different year.
diagram.station <- function(x,it=NULL,...) {
  yrs <- as.numeric(rownames(table(year(x))))
  d <- dim(x)
  #print(yrs)
  ny <- length(yrs)
  j <- 1:ny
  col <- rgb(j/ny,abs(sin(pi*j/ny)),(1-j/ny),0.2)
  class(x) <- "zoo"

  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") )
      unit <- expression(degree*C) else
      unit <- attr(x,'unit')
  eval(parse(text=paste("main <- expression(paste('Seasonal evaution: ',",
               attr(x,'variable'),"))")))
  dev.new()
  par(bty="n")
  z <- coredata(x)
  plot(c(0,365),1.25*range(z,na.rm=TRUE),
       type="n",xlab="",
       main=main,
       sub=attr(x,'location'),ylab=unit)
  grid()
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    if (is.null(d)) points(t,coredata(y),lwd=2,col=col[i],pch=19,cex=0.5) else
                    points(rep(t,d[2]),coredata(y),lwd=2,col=col[i],pch=19,cex=0.5)
  }
  if (!is.null(it)) {
    y <- window(x,start=as.Date(paste(it,'-01-01',sep='')),
                    end=as.Date(paste(it,'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(it,'-01-01',sep='')))
  }
  if (is.null(d)) points(t,coredata(y),col="black",cex=0.7) else
                  points(rep(t,d[2]),coredata(y),col="black",cex=0.7)

  par(new=TRUE,fig=c(0.70,0.85,0.70,0.85),mar=c(0,3,0,0),
      cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:ny,1:ny)
  image(1:2,yrs,colbar,col=col)
}

#seNorge
#nodata	10000
#nobits	16
#header1	*Temperatur
#header2	*Siste døgn (18-18 UTC)
#legend	*Grader Celsius
#From	To	R	G	B	Forklaring
#2931	10000	204	0	0	*Over 20 
#2881	2931	255	25	0	*15 - 20
#2831	2881	255	102	0	*10 - 15
#2781	2831	255	179	0	*5 - 10
#2761	2781	255	230	77	*3 - 5
#2741	2761	255	255	64	*1 - 3
#2731	2741	255	255	190	*0 - 1
#2721	2731	217	255	255	*÷1 - 0
#2701	2721	179	255	255	*÷3 - ÷1
#2681	2701	128	235	255	*÷5 - ÷3
#2631	2681	64	204	255	*÷10 - ÷5
#2581	2631	0 	153	255	*÷15 - ÷10
#2531	2581	0 	25	255	*÷20 - ÷15
#0	2531	0	0	153	*Under ÷20


#nodata	10000
#nobits	16
#header1	*Nedbør
#header2	*Siste døgn (06-06 UTC)
#legend	*mm
#From	To	R	G	B	Forklaring
#1500	10000	0	0	153	*Over 150 
#750	1500	0 	25	255	*75 - 150 
#500	750	0 	153	255	*50 - 75
#300	500	64	204	255	*30 - 50
#200	300	128	235	255	*20 - 30
#100	200	179	255	255	*10 - 20
#1	100	217	255	255	*Under 10
#0	1	229	229	229	*Ikke nedbør

colscal <- function(n=14,col="bwr",test=FALSE) {

  test.col <- function(r,g,b) {
    dev.new()
    par(bty="n")
    plot(r,col="red")
    points(b,col="blue")
    points(g,col="green")
  }

#R	G	B
  seNorgeT <- c(204,  0,    0,	
               255, 25,    0,	
               255, 102,   0,	
               255, 179,   0,	
               255, 230,  77,	
               255, 255,  64,	
               255, 255, 190,	
               217, 255, 255,	
               179, 255, 255,	
               128, 235, 255,	
               64, 204, 255,	
               0, 153, 255,	
               0,  25, 255,	
               0,   0, 153)	
  dim(seNorgeT) <- c(3,14)

  seNorgeP <- c(0, 0, 153,
                0, 25, 255,
                0, 153, 255,
                64, 204, 255,
                128, 235, 255,
                179, 255, 255,
                217, 255, 255,
                229, 229, 229)
  
  dim(seNorgeP) <- c(3,8)
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
  } else if (col=="t2m") {
    r <- round(0.01*approx(seNorgeT[1,],n=n)$x)
    g <- round(0.01*approx(seNorgeT[2,],n=n)$x)
    b <- round(0.01*approx(seNorgeT[3,],n=n)$x)
    col <- rgb(b,g,r)
  } else if (col=="precip") {
    r <- round(0.01*approx(seNorgeP[1,],n=n)$x)
    g <- round(0.01*approx(seNorgeP[2,],n=n)$x)
    b <- round(0.01*approx(seNorgeP[3,],n=n)$x)
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


diagnose <-function(x,...) UseMethod("diagnose")

diagnose.default <- function(x,...) {
}

diagnose.comb <- function(x) {
  n.app <- attr(x,'n.apps')
  cols <- c("black","red","blue","darkgreen","darkred","darblue",
            "grey","green","mangenta","cyan")
  par(bty="n")
  plot(colMeans(x),cex=0.5,pch=19,
       main="Grid box mean value for combined fields")
  for (i in 1:n.app) {
    col <- cols[i%%10+1]
    y <- eval(parse(text=paste("attr(x,'appendix.",i,"')",sep="")))
    points(colMeans(y,na.rm=TRUE),col=col,cex=0.3)
  }
}


diagnose.eof <- function(x) {
  if (inherits(x,'comb')) y <- diagnose.comb.eof(x) else
                          y <- x
  return(y)
}

diagnose.comb.eof <- function(x) {

  ACF <- function(x) acf(x,plot=FALSE,na.action=na.omit)$acf[2]
  sign <- function(x,y) {z<-x*y; z[z<0] <- -1; z[z>0] <- 1; z}
  
  stopifnot(!missing(x), inherits(x,"eof"),inherits(x,"comb"))
  #print("diagnose.comb.eof")

  # The original field, e.g. reanalyses
  Y <- zoo(coredata(x),order.by=index(x))
  n <- attr(x,'n.apps')
  m <- length(attr(x,'eigenvalues'))
  dm <- rep(NA,n*m); dim(dm) <- c(n,m)
  sr <- dm; ar <- sr

  # The appended fields, e.g. GCM results
  rowname <- rep("GCM",n)
  for ( i in 1:n ) {
    eval(parse(text=paste("z <- attr(x,'appendix.",i,"')",sep="")))
    y <- zoo(coredata(z),order.by=index(z))

    # Extract a comon period:
    X <- merge(Y,y,all=FALSE)
    Ym <- apply(coredata(X),2,mean,na.rm=TRUE)
    #print(Ym)
    #plot(Ym)
    Ys <- apply(coredata(X),2,sd,na.rm=TRUE)
    AR <- apply(coredata(X),2,ACF)
    dm[i,] <- Ym[1:m] - Ym[(m+1):(2*m)]
    # ratio: GCM/original
    # The problem is when the denominator is close to zero...
    sr[i,] <- (0.01 + Ys[(m+1):(2*m)])/(0.01 + Ys[1:m])
    ar[i,] <- (0.01 + AR[(m+1):(2*m)])/(0.01 + AR[1:m])*sign(AR[(m+1):(2*m)],AR[1:m])
    if (!is.null(attr(z,'source'))) rowname[i] <- attr(z,'source') else
    if (!is.null(attr(z,'model_id'))) rowname[i] <- attr(z,'model_id') 
  }
  rownames(dm) <- rowname
  rownames(sr) <- rowname
  rownames(ar) <- rowname
  diag <- list(mean.diff=dm,sd.ratio=sr,autocorr.ratio=ar,
               common.period=range(index(Y)),sd0=Ys,
               calibrationdata=attr(x,'source'))
  attr(diag,'variable') <- attr(x,'variable')
  #print(summary(diag))
  attr(diag,'history') <- history.stamp(x)
  class(diag) <- c("diagnose","comb","eof","list")
  invisible(diag)
}


diagnose.mvr <- function(x) {
  print("Not finished")
}

diagnose.cca <- function(x) {
  par(bty="n")
  plot(x$r,pch=19,cex=1.5,main="Canonical correlations",
       ylim=c(-1,1),ylab="correlation",xlab="pattern number")
  lines(c(0,length(x$r)),rep(0,2),col="grey")
  grid()
}

# Display cross-validation and statistics on the residual
diagnose.ds <- function(x,plot=FALSE) {

  # the attribute 'evaluation' contains cross-validation
  if (!is.null(attr(x,'evaluation'))) xval <- attr(x,'evaluation') else
                                      xval <- crossval(x)
  y <- as.residual(x)
  z <- as.original.data(x)
  anova <- summary(attr(x,'model'))
  eof <- attr(x,'eof')
  if (inherits(eof,'comb')) bias.diag <- diagnose(eof) else
                            bias.diag <- NULL

  spectrum(coredata(y),plot=FALSE) -> s
  sp <- data.frame(y=log(s$spec),x=log(s$freq))
  if (length(dim(y))==0) {
    beta <- -summary(lm(y ~ x, data=sp))$coefficient[2]
    beta.error <- summary(lm(y ~ x, data=sp))$coefficient[4]
    ar1 <- acf(y,plot=FALSE)$acf[2]
  } else {beta <- NA; beta.error <- NA; ar1 <- NA}
  
  if (plot) {
    plot(xval)

    dev.new()
    par(bty="n",mfcol=c(3,2))
    plot(y)
    
    acf(y)

    plot(z,y)
    spectrum(y)

    qqnorm(y)
    qqline(y)

    if  (!is.null(attr(x,'diagnose'))) 
      plot(attr(x,'diagnose'))
  }
  
  diagnostics <- list(residual=y,anova=anova,xval=xval,bias.diag=bias.diag,
                      ar1=ar1,beta=beta, H=(beta+1)/2, beta.error=beta.error)
  return(diagnostics)
}


diagnose.dsensemble <- function(x,plot=TRUE,type='target',...) {
  # Trend-evaluation: rank
  # Counts outside 90% confidence: binomial distrib. & prob.
  stopifnot(!missing(x),inherits(x,"dsensemble"))
  z <- x
  # Remove the results with no valid data:
  n <- apply(z,2,FUN=nv)
  z <- subset(z,is=(1:length(n))[n > 0])
  
  d <- dim(z)
  t <- index(z)
  y <- attr(x,'station')
  
  # statistics: past trends
  #browser()
  i1 <- is.element(year(y)*100 + month(y),year(z)*100 + month(z))
  i2 <- is.element(year(z)*100 + month(z),year(y)*100 + month(y))
  obs <- data.frame(y=y[i1],t=year(y)[i1])
  #print(summary(obs)); print(sum(i1)); print(sum(i2)); browser()
  deltaobs <- lm(y ~ t,data=obs)$coefficients[2]*10  # deg C/decade
  deltagcm <- rep(NA,d[2])
  for (j in 1:d[2]) {
    gcm <- data.frame(y=z[i2,j],t=year(z)[i2])
    deltagcm[j] <- lm(y ~ t,data=gcm)$coefficients[2]*10  # deg C/decade
  }
  robs <- round(100*sum(deltaobs < deltagcm)/d[2])
  #print(deltaobs); print(deltagcm); print(order(c(deltaobs,deltagcm))[1])

  # apply to extract mean and sd from the selected objects:
  mu <- apply(coredata(z),1,mean,na.rm=TRUE)
  si <- apply(coredata(z),1,sd,na.rm=TRUE)
  q05 <- qnorm(0.05,mean=mu,sd=si)
  q95 <- qnorm(0.95,mean=mu,sd=si)
  # number of points outside conf. int. (binom)
  above <- y[i1] > q95[i2]
  below <- y[i1] < q05[i2]
  #browser()
  outside <- sum(above) + sum(below)
  N <- sum(i1)
  
  if (plot) {
    x <- -round(200*(0.5-pbinom(outside,size=N,prob=0.1)),2)
    y <- -round(200*(0.5-pnorm(deltaobs,mean=mean(deltagcm),sd=sd(deltagcm))),2)
    #print(c(x,y))
    
    par(bty="n",xaxt="n",yaxt="n")
    plot(c(-100,100),c(-100,100),type="n",ylab="magnitude",xlab="trend")
    
    bcol=c("grey95","grey40")
    for (i in 1:10) {
      r <- (11-i)*10
      polygon(r*cos(pi*seq(0,2,length=360)),
              r*sin(pi*seq(0,2,length=360)),
              col=bcol[i %% 2 + 1],border="grey15")
    }
    for (i in seq(0,90,by=1))
      points(x,y,pch=19,cex=2 - i/50,col=rgb(i/90,0,0))
  }
  diag <- list(robs=robs,deltaobs=deltaobs,deltagcm=deltagcm,
               outside=outside,above=above,below=below,
               y=y[i1],N=N,i1=i1,
               mu=zoo(mu,order.by=index(x)),
               si=zoo(si,order.by=index(x)),
               q05=zoo(q05,order.by=index(x)),
               q95=zoo(q95,order.by=index(x)))
  attr(diag,'history') <- history.stamp(x)

  invisible(diag)
}

diagnose.station <- function(x,main='Data availability',
                            xlab='',ylab='station',
                            sub=src(x),...) {
  d <- dim(x)
  if (is.null(d)) return('Need more than one station')
  par(mar=c(5, 4, 4, 5),las=1,xpd=TRUE,cex.lab=0.5,cex.axis=0.5)
  
  image(index(x),1:d[2],coredata(x),
        main=main,xlab=xlab,ylab=ylab,
        sub=sub,...)
  axis(4,at=1:d[2],labels=substr(loc(x),1,6),cex.lab=0.5,col='grey')
  par(xpd=FALSE)
  grid(nx=d[1],ny=d[2])
}
