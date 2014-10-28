



iid.test <- function(x,...) UseMethod("iid.test")

iid.test.station <- function(x,...) {
  # Re-orders the station data into parallel time series for each calendar
  # month into new matrix X. Then apply the iid.test to this matrix.

  ts2mon <- function(x) {
    print('ts2mon')
    yrs <- year(x); n <- length(rownames(table(yrs)))
    d <- dim(x)
  # Test for multiple series:
    if (is.null(d)) m <- 1 else  # single series
                    m <- d[2]    # multiple
    X <- matrix(rep(NA,m*n*12),n,m*12)
    dim(X) <- c(n,m,12)
    print(dim(X))
  
    for (i in 1:12) {
      y <- subset(x,it=month.abb[i])
    #print(dim(y))
      X[(1 + n-length(y)):n,1:m,i] <- coredata(y)
    }
    dim(X) <- c(n,m*12)
    attr(X,'description') <- 'data matrix re-orderd on month and location'
    attr(X,'original_dimensions') <- c(n,m,12)
    attr(X,'history') <- history.stamp(x)
    return(X)
  }

  
  print('iid.test.station')
  X <- ts2mon(x)
  good <- is.finite(rowMeans(X))
  iid <- iid.test.default(X[good,])
  invisible(iid)
}

iid.test.field <- function(x,...) {
  # Uses EOFs to account for spatial co-variance, and test the PCs rather
  # than the grid points.
  # Re-orders the PCs into parallel time series for each calendar
  # month into new matrix X. Then apply the iid.test to this matrix. 

  print('iid.test.field')
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
  iid <- iid.test.default(X[good,])
  invisible(iid)
}



iid.test.default <- function(x,plot=TRUE,Monte.Carlo=TRUE,
                             N.test=200,rev.plot.rev=TRUE) {
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
         xlab="time",ylab="location")
    par(col.axis="black"); axis(1)
    lines(c(-5,t.r[1]+6),rep(t.r[2],2)+0.5,lwd=3)
    par.0 <- par(); par(srt=90)
    text(0,round(t.r[2]/2),"Forward",cex=1,vfont=c("sans serif","italic"))
    text(0,round(3*t.r[2]/2),"Backward",cex=1,vfont=c("sans serif","italic"))
    par(par.0)
  }

  for (ir in 1:t.r[2]) {
    record.stats <- n.records(Y[,ir])
    N.records[ir] <- record.stats$N
    events[,ir] <- as.numeric(record.stats$events)
    events.rev[,ir] <- as.numeric(record.stats$events.rev)

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
    }
  }

  if (plot) {
    for (it in 2:t.r[1]) {
      CI.95[it,] <- qbinom(p=c(0.025,0.975),size=t.r[2],prob=1/it)
      p.val[it] <- pbinom(events[it,],size=t.r[2],prob=1/it)
      if ( (sum(events[it,]) < CI.95[it,1]) |
          (sum(events[it,]) > CI.95[it,2]) ) {
        i.cluster[it] <- TRUE
        lines(rep(it,2),c(0,t.r[2]),lwd=2,lty=2,col="pink")
      }
      CI.95.rev[it,] <- qbinom(p=c(0.025,0.975),size=t.r[2],prob=1/it)
      p.val.rev[it] <- pbinom(events[it,],size=t.r[2],prob=1/it)
      if ( (sum(events.rev[t.r[1]-it+1,]) < CI.95.rev[it,1]) |
          (sum(events.rev[t.r[1]-it+1,]) > CI.95.rev[it,2]) ) {
        i.cluster.rev[it] <- TRUE
        if (rev.plot.rev)
          lines(rep(t.r[1]-it+1,2),c(t.r[2]+1,2*t.r[2]),lwd=2,
                lty=2,col="pink") else
          lines(rep(it,2),c(t.r[2]+1,2*t.r[2]),lwd=2,lty=2,col="pink") 
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
  test.results <- iid.test(rnd,plot=plot,Monte.Carlo=Monte.Carlo)
  invisible(test.results)
}

daily.station.records <- function(obs,element="precip",subsample=5,tolerance=2,remove.zeroes=FALSE,rev.plot.rev=FALSE) {

  if (class(obs)[2] != "daily.station.record") 
     stop("Need a 'daily.station.record' object!")
  years <- as.numeric(rownames(table(obs$yy)))
  ny <- length(years)
  dat <- rep(NA,ny*366); dim(dat) <- c(ny,366)

  for (i in 1:ny) {
    iyear <- is.element(obs$yy,years[i])
    ii <- julday(obs$mm[iyear],obs$dd[iyear],obs$yy[iyear]) -
          julday(1,1,years[i])+1
    data.thisyear <- eval(parse(text=paste("obs$",element,"[iyear]",sep="")))
    print(c(years[i],sum(iyear),NA,range(ii),NA,length(data.thisyear)))
  
    if (sum(is.finite(data.thisyear)) > 2) plot(data.thisyear)
    dat[i,ii] <- data.thisyear
  }
  
  #plot(dat[1,],main=obs$location,ylab=element,xlab="Day in the year",
  #     pch=20,cex=0.8,col="grey50")
  #for (i in 1:ny) points(dat[i,],pch=20,cex=0.8,col="grey50")

  print(paste("sub-sample every",subsample,"points"))

  ykeep <- rep(TRUE,366)
  for (i in 1:length(ykeep)) {
    if (mod(i,subsample) != 1) ykeep[i] <- FALSE
    if (sum(!is.finite(dat[,i])) > tolerance) ykeep[i] <- FALSE
    if ( (remove.zeroes) & (sum(dat[,i]==0) > tolerance) ) ykeep[i] <- FALSE
  }
  dat <- dat[,ykeep] # extra days during leap year will bias the results 

  #image(log(dat)); x11()

  dev.new()
  iid.test(dat,rev.plot.rev=rev.plot.rev)
}

n.records <- function(x) {
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
  y.rev <- rev(y)
  N <- 1; N.rev <- N
  t <- rep(1,m); t.rev <- rep(m,m)
  events <- rep(FALSE,m); events.rev <- events
  events[1] <- TRUE; events.rev[m] <- TRUE
  for (i in 2:m) {
    if (y[i] > max(y[1:(i-1)],na.rm=TRUE)) {
      N <- N + 1
      t[N] <- i
      events[i] <- TRUE
    }
#    if (y[m-i+1] > max(y[(m-i+2):m],na.rm=TRUE)) {
#      N.rev <- N.rev + 1
#      t.rev[N.rev] <- length(y)-i
#      events.rev[length(y)-i] <- TRUE
#    }
# Changed 18.08.2006: REB - a more readable code... shouldn't make any difference.
    if (y.rev[i] > max(y.rev[1:(i-1)],na.rm=TRUE)) {
      N.rev <- N.rev + 1
      t.rev[N.rev] <- length(y)-i+1
      events.rev[length(y)-i+1] <- TRUE
    }
  }
  t <- t[1:N]; t.rev <- t.rev[1:N.rev]
  records <- list(N=N,t=t,events=events,N.rev=N.rev, 
                  t.rev=t.rev, events.rev=events.rev)
  invisible(records)
} 
