

iid.test <- function(x,...) UseMethod("iid.test")

iid.test.station <- function(x,verbose=TRUE,...) {
  # Re-orders the station data into parallel time series for each calendar
  # month into new matrix X. Then apply the iid.test to this matrix.

#  ts2mon <- function(x,verbose=TRUE) {
#    if (verbose) print('ts2mon')
#    yrs <- year(x); n <- length(rownames(table(yrs)))
#    d <- dim(x)
#  # Test for multiple series:
#    if (is.null(d)) m <- 1 else  # single series
#                    m <- d[2]    # multiple
#    X <- matrix(rep(NA,m*n*12),n,m*12)
#    dim(X) <- c(n,m,12)
#    if (verbose) print(dim(X))
#  
#    for (i in 1:12) {
#      y <- subset(x,it=month.abb[i],verbose=verbose)
#      y <- aggregate(y,year,FUN='max',na.rm=TRUE) # one estimate for each month/year
#      if (verbose) print(dim(y))
#      if (verbose) print(paste(month.abb[i],(1 + n-length(index(y))),
#                               length((1 + n-length(index(y))):n)))
#      if (verbose) print(length(index(y)))
#      X[(1 + n-length(index(y))):n,1:m,i] <- coredata(y)
#    }
#    if (verbose) print('set dimensions')
#    dim(X) <- c(n,m*12)
#    attr(X,'description') <- 'data matrix re-orderd on month and location'
#    attr(X,'original_dimensions') <- c(n,m,12)
#    attr(X,'history') <- history.stamp(x)
#    return(X)
#  }

  
  if (verbose) print('iid.test.station')
  #X <- ts2mon(x,verbose=verbose)
  #if (verbose) print('weed out bad data')
  #good <- is.finite(rowMeans(X))
  X <- as.monthly(x,FUN='max')
  iid <- iid.test.default(X,verbose=verbose)
  invisible(iid)
}

iid.test.field <- function(x,verbose=TRUE,...) {
  # Uses EOFs to account for spatial co-variance, and test the PCs rather
  # than the grid points.
  # Re-orders the PCs into parallel time series for each calendar
  # month into new matrix X. Then apply the iid.test to this matrix. 

  if (verbose) print('iid.test.field')
  yrs <- year(x); n <- length(rownames(table(yrs)))
  X <- matrix(rep(NA,20*12*n),n,12*20)
  for (i in 1:12) {
    eof <- EOF(subset(x,it=month.abb[i]))
    # The PCs and EOFs have a random sign. ensure that these are
    # consistent with the original data:
    y <- spatial.avg.field(subset(x,it=month.abb[i]))
    #print(length(y)); dim(eof)
    r <- y %*% coredata(eof)
    #print(paste('inner-product: ',round(r,2)))
    eof[,r < 0] <- -eof[,r < 0]
    m <- length(index(eof))
    #str(eof)
    #print(c(n,m,(1 + n-m):n,1:20 + (i-1)*20));
    #print(dim(X[(1 + n-m):n,1:20 + (i-1)*20]))
    #print(dim(coredata(eof)))
    X[(1 + n-m):n,1:20 + (i-1)*20] <- coredata(eof)
  }
  good <- is.finite(rowMeans(X))
  iid <- iid.test.default(X[good,],verbose=verbose)
  invisible(iid)
}



iid.test.default <- function(x,plot=TRUE,Monte.Carlo=TRUE,
                             N.test=200,rev.plot.rev=TRUE,verbose=TRUE) {
  if (verbose) print('iid.test.default')
  Y <- as.matrix(x)
  Y[!is.finite(Y)] <- NA
  t.r <- dim(Y)
  events <- matrix(rep(FALSE,t.r[1]*t.r[2]),t.r[1],t.r[2])
  events.rev <- events
  N.records <- rep(NA,t.r[2])
# Use binomial distribution to look for suspicious clusters (dependencies)
  CI.95 <- rep(t.r[2],t.r[1]*2); dim(CI.95) <- c(t.r[1],2); CI.95.rev <- CI.95 
  p.val <- rep(NA,t.r[1]); i.cluster <- rep(FALSE,t.r[1])
  p.val.rev <- p.val; i.cluster.rev <- i.cluster

  if (plot) {
    par(col.axis="white")
    plot(c(1,t.r[1]),c(1,2*t.r[2]),type="n",main="iid-test",
         xlab="record length",ylab="location")
    par(col.axis="black")
    axis(1)
    lines(c(-5,t.r[1]+6),rep(t.r[2],2)+0.5,lwd=3)
    par.0 <- par(); par(srt=90)
    text(0,round(t.r[2]/2),"Forward",cex=1,vfont=c("sans serif","italic"))
    text(0,round(3*t.r[2]/2),"Backward",cex=1,vfont=c("sans serif","italic"))
    par(par.0)
  }

  
  for (ir in 1:t.r[2]) {
    #record.stats <- n.records(Y[,ir])
    record.stats <- n.records(subset(x,is=ir))
    if (verbose) str(record.stats)
    N.records[ir] <- record.stats$N
    events[record.stats$t,ir] <- TRUE
    events.rev[record.stats$t.rev,ir] <- TRUE

    if (plot) {
      # Timing index for record.
      t1 <- record.stats$t
      if (rev.plot.rev) t2 <- record.stats$t.rev  else 
                        t2 <- t.r[1] - record.stats$t.rev + 1

      #lines(c(1,t.r[1]),rep(ir,2),col="grey70")            
      points(t1,rep(ir,record.stats$N)+0.025,pch=20,cex=1.50,col="grey30")
      points(t1+0.05,rep(ir,record.stats$N)+0.050,pch=20,cex=0.70,col="grey50")
      points(t1+0.07,rep(ir,record.stats$N)+0.075,pch=20,cex=0.50,col="grey70")
      points(t1+0.1,rep(ir,record.stats$N)+0.100,pch=20,cex=0.30,col="white")   
   
      #lines(c(1,t.r[1]),rep(ir+t.r[2],2),col="grey70")
      points(t2,rep(ir,record.stats$N.rev)+0.025+t.r[2],pch=20,
             cex=1.50,col="grey30")
      points(t2+0.05,rep(ir,record.stats$N.rev)+0.050+t.r[2],pch=20,
             cex=0.70,col="grey50")
      points(t2+0.07,rep(ir,record.stats$N.rev)+0.075+t.r[2],pch=20,
             cex=0.50,col="grey70")
      points(t2+0.1,rep(ir,record.stats$N.rev)+0.100+t.r[2],pch=20,
             cex=0.30,col="white")
      text(loc(x)[ir],0,ir,pos=3)
    }
  }

  if (plot) {
    for (it in 2:t.r[1]) {
      CI.95[it,] <- qbinom(p=c(0.025,0.975),size=t.r[2],prob=1/it)
      p.val[it] <- pbinom(events[it,],size=t.r[2],prob=1/it)
      if ( (sum(events[it,]) < CI.95[it,1]) |
          (sum(events[it,]) > CI.95[it,2]) ) {
        i.cluster[it] <- TRUE
        lines(rep(it,2),c(0,t.r[2]),lwd=1,lty=2,col=rgb(1,0.5,0.5,0.3))
      }
      CI.95.rev[it,] <- qbinom(p=c(0.025,0.975),size=t.r[2],prob=1/it)
      p.val.rev[it] <- pbinom(events[it,],size=t.r[2],prob=1/it)
      if ( (sum(events.rev[t.r[1]-it+1,]) < CI.95.rev[it,1]) |
          (sum(events.rev[t.r[1]-it+1,]) > CI.95.rev[it,2]) ) {
        i.cluster.rev[it] <- TRUE
        if (rev.plot.rev)
          lines(rep(t.r[1]-it+1,2),c(t.r[2]+1,2*t.r[2]),lwd=1,
                lty=2,col=rgb(1,0.5,0.5,0.3)) else
          lines(rep(it,2),c(t.r[2]+1,2*t.r[2]),lwd=1,lty=2,col=rgb(1,0.5,0.5,0.3)) 
      }    
    }
  }

  events[!is.finite(Y)] <- NA
  events.rev[!is.finite(Y)] <- NA
  record.density <- rowMeans(events,na.rm=TRUE)
  record.density.rev <- rev(rowMeans(events.rev,na.rm=TRUE))
  N <- length(record.density)

  q025=rep(NA,N); q975=q025    
  if (Monte.Carlo) {
    print(paste("Please be patient -",N.test,"Monte Carlo runs in progress..."))
    record.mc <- rep(NA,2*N.test*N); dim(record.mc) <- c(N,N.test,2)
    for (ii in 1:N.test) {
      mc.stats <- test.iid.test(d=dim(events),plot=FALSE,Monte.Carlo=FALSE)  
      record.mc[,ii,1] <- cumsum(mc.stats$record.density)
      record.mc[,ii,2] <- cumsum(mc.stats$record.density.rev)
    } 

    for (i in 1:N) {
      q025[i] <- quantile(record.mc[i,,],0.025)
      q975[i] <- quantile(record.mc[i,,],0.955)
    }
    sub <- paste("Shaded region= 95% conf.int. from Monte-Carlo with N=",N.test)
  } else {
#    for (i in 1:N) {
#      q025[i] <- sum(CI.95[1:i,1])/t.r[2]
#      q975[i] <- sum(CI.95[1:i,2])/t.r[2]
#    }     
    sub <- ""
  }

  if (plot) {
    dev.new(); par(col.axis="white")
    Time <- 1:N
    plot(Time,exp( cumsum( 1/(1:N)) ),type="l",lwd=3,col="grey60",
           xlab="Time",ylab="exp( cumsum( 1/(1:n)) )",
         main="Observed & Expected number of record-events",
           sub=sub)
    par(col.axis="black")
    axis(1)
    axis(2,at=exp(1:(2*sum(record.density,na.rm=TRUE))),
         label=1:(2*sum(record.density,na.rm=TRUE)))
    legend(1,exp(sum(1/(1:N))),c("Theoretical","Forward","Backward"),
           pch=c(-1,19,21),lwd=c(3,1,1),lty=c(1,0,0),col=c("grey60",
                                                    rep("black",2)))

    if (Monte.Carlo) {
      polygon(c(Time,rev(Time)),c(exp(q025),rev(exp(q975))),
              col="grey90",border="grey85")
      lines(Time,exp( cumsum( 1/(1:N)) ),lwd=3,col="grey60")
    }
    grid()
    points(Time,exp(cumsum(record.density)),pch=19,cex=0.9)
    points(Time,exp(cumsum(record.density.rev)),pch=21,cex=0.9)

    
    dev.new()
    plot(Time,1/(1:N),type="n",xlab="Time",
         main="Observed & Expected record-occurence",
           sub="95% Confidence interval from binomial distribution",
           ylab="record-density")
    for (i in 1:N) {
      lines(rep(Time[i],2)-0.2,c(CI.95[i,1]/t.r[2],CI.95[i,2]/t.r[2]),lwd=2,
            col="grey40")
      lines(Time[i]+c(-0.2,0.2),rep(CI.95[i,1]/t.r[2],2),lwd=1,col="grey40")
      lines(Time[i]+c(-0.2,0.2),rep(CI.95[i,2]/t.r[2],2),lwd=1,col="grey40")
      lines(rep(Time[i],2)+0.2,c(CI.95.rev[i,1]/t.r[2],CI.95.rev[i,2]/t.r[2]),
            lwd=2,col="grey70")
      lines(Time[i]+c(-0.2,0.2),rep(CI.95.rev[i,1]/t.r[2],2),lwd=1,col="grey70")
      lines(Time[i]+c(-0.2,0.2),rep(CI.95.rev[i,2]/t.r[2],2),lwd=1,col="grey70")
    }
    points(Time,record.density,pch=20,cex=0.9,col="grey20")    
    points(Time,record.density.rev,pch=21,cex=0.9)
    lines(Time,1/(1:N),lwd=2,col="red")
  }
  
  results <- list(record.density=record.density,
                  record.density.rev=record.density.rev,
                  CI.95=CI.95,p.val=p.val,i.cluster=i.cluster,
                  CI.95.rev=CI.95.rev,p.val.rev=p.val.rev,
                  i.cluster.rev=i.cluster.rev)
  invisible(results)
}

  

test.iid.test <- function(distr="rnorm",d=c(70,50),plot=TRUE,
                          Monte.Carlo=TRUE) {
  rnd <- eval(parse(text=paste(distr,"(",d[1]*d[2],")",sep="")))
  dim(rnd) <- c(d[1],d[2])
  test.results <- iid.test(zoo(rnd,order.by=1:d[1]),
                           plot=plot,Monte.Carlo=Monte.Carlo)
  invisible(test.results)
}


n.records <- function(x,verbose=FALSE) {
  if (verbose) print('n.records')
  y <- x
  m <- length(y)
  y[!is.finite(y)] <- min(y,na.rm=TRUE)
  if (length(rownames(table(y))) < 0.99 * length(y)) {
    print("---Warning---Warning---Warning---Warning---Warning---Warning---")
    print("r.records (iid.test): Warning, the time series contains many similar values!")
    print("The test does not work for cases where ties are common")
    print("See Benestad (2004) 'Record-values, non-stationarity tests and extreme value distributions' Global and Planetary Change, 44, 11-26")
    print("http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6VF0-4D6373Y-2&_coverDate=12%2F01%2F2004&_alid=228212815&_rdoc=1&_fmt=&_orig=search&_qd=1&_cdi=5996&_sort=d&view=c&_acct=C000056508&_version=1&_urlVersion=0&_userid=2181464&md5=632559476e84eb8c48287cf8038690d2")
  }
  y.rev <- rev(y); index(y.rev) <- index(y) 
  if (verbose) {str(y); str(y.rev)}
  N <- 1; N.rev <- N
  t <- rep(1,m); t.rev <- rep(m,m)
  events <- rep(FALSE,m); events.rev <- events
  
  if (verbose) print('fast algorithm')
  events <- records(y,verbose=verbose)
  if (verbose) print('reverse series')
  events.rev <- records(y.rev,verbose=verbose)
    
  if (is.numeric(events)) { 
      if (verbose) print('single series')
      N <- sum(is.finite(events)) 
      t <- attr(events,'t')
      N.rev <- sum(is.finite(events.rev))
      t.rev <- attr(events.rev,'t')
  } else if (is.list(events)) {
    if (verbose) print('matrix')
      N <- lapply(events,function(x) sum(is.finite(x)))
      t <- lapply(events,function(x) attr(x,'t'))
      N <- lapply(events.rev,function(x) sum(is.finite(x)))
      t.rev <- lapply(events.rev,function(x) attr(x,'t'))
  } else stop(paste('n.records - naot programmed to handle',class(events)))

  if (verbose) print('organise into list object')
  records <- list(N=N,t=t,events=events,N.rev=N.rev, 
                  t.rev=t.rev, events.rev=events.rev)
  invisible(records)
} 

## This algorithm is faster than the older code that used for-loop
## Search 'back-ward' statring with the highest value:

records <- function(x,verbose=FALSE,diff=FALSE) {
  if (verbose) print('records')
  if (!is.null(dim(x))) {
    if (verbose) print('matrix: apply recursive routine')
    r <- apply(x,2,'records')
    return(r)
  }
  n <- length(x)
  r <- rep(NA,length(x)); t <- rep(NA,length(x))
  ii <- 1; i <- length(x)
  if (verbose) print(n)
  while(length(x) >= 1 & i > 1) {
    z <- max(x,na.rm=TRUE)
    r[ii] <- z
    i <- (1:n)[is.element(x,z)][1]
    t[ii] <- i; ii <- ii  + 1
    if (is.finite(i)) x <- x[1:(i-1)] else i <- 0
    if (verbose) print(c(z,ii,i,length(x)))
  }

  t <- t[is.finite(r)]
  r <- r[is.finite(r)]
  r <- rev(r)
  if (diff) r <- c(min(r),diff(r))
  attr(r,'t') <- rev(t)
  return(r)
}

test.records <- function(N=1000) {
  y <- rnorm(N) + seq(0,1,length=N)
  rbv <- records(y)
  plot(y,type='l')
  points(attr(rbv,'t'),rbv,col='red')
}

