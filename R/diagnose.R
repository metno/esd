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
  beta <- -summary(lm(y ~ x, data=sp))$coefficient[2]
  beta.error <- summary(lm(y ~ x, data=sp))$coefficient[4]
  
  if (plot) {
    plot(xval)

    dev.new()
    par(bty="n",mfcol=c(3,2))
    plot(y)
    
    acf(y) -> ar

    plot(z,y)
    spectrum(y)

    qqnorm(y)
    qqline(y)

    if  (!is.null(attr(x,'diagnose'))) 
      plot(attr(x,'diagnose'))
  }
  
  diagnostics <- list(residual=y,anova=anova,xval=xval,bias.diag=bias.diag,
                      ar1=ar$acf[2],beta=beta, H=(beta+1)/2, beta.error=beta.error)
  return(diagnostics)
}

# Show the temperatures against the day of the year. Use
# different colours for different year.
diagnose.station <- function(x,it=NULL,...) {
  yrs <- as.numeric(rownames(table(year(x))))
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
    points(t,coredata(y),lwd=2,col=col[i],pch=19,cex=0.5)
  }
  if (!is.null(it)) {
    y <- window(x,start=as.Date(paste(it,'-01-01',sep='')),
                    end=as.Date(paste(it,'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(it,'-01-01',sep='')))
  }
  points(t,coredata(y),col="black",cex=0.7)

  par(new=TRUE,fig=c(0.70,0.85,0.70,0.85),mar=c(0,3,0,0),
      cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:ny,1:ny)
  image(1:2,yrs,colbar,col=col)
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

