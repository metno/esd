# K Parding, 15.10.2015

track <- function(x,...) UseMethod("track")

track.events <- function(x,verbose=FALSE,...) {
  if(verbose) print("track.events")
  track.default(x,...)
}

track.default <- function(x,x0=NULL,it=NULL,is=NULL,dmax=1.2E6,amax=90,
                         nmax=124,nmin=3,dE=0.3,dN=0,dmin=5E5,
                         lplot=FALSE,progress=TRUE,verbose=FALSE) {
  if(verbose) print("track.default")
  x <- subset(x,it=!is.na(x["date"][[1]]))
  x <- subset(x,it=it,is=is)
  yrmn <- as.yearmon(strptime(x["date"][[1]],"%Y%m%d"))
  if (length(unique(yrmn))>1 & length(yrmn)>1000) {
    x.tracked <- NULL
    if (progress) pb <- txtProgressBar(style=3)
    for (i in 1:length(unique(yrmn))) {
      if(verbose) print(unique(yrmn)[i])
      x.y <- subset(x,it=(yrmn==unique(yrmn)[i]))
      if (is.null(x.tracked)) {
        x.t <- Track(x.y,x0=x0,lplot=lplot,x0cleanup=FALSE,dE=dE,dN=dN,
                     amax=amax,dmax=dmax,dmin=dmin,nmax=nmax,nmin=nmin,
                     progress=FALSE,verbose=verbose)
        x.tracked <- x.t$y
      } else {
        x.t <- Track(x.y,x0=x.tracked,lplot=lplot,dE=dE,dN=dN,
                     amax=amax,dmax=dmax,dmin=dmin,nmax=nmax,nmin=nmin,
                     progress=FALSE,verbose=verbose)
        x.tracked <- x.t$y0
        x.tracked <- merge(x.tracked,x.t$y,all=TRUE)
      }
      if (progress) setTxtProgressBar(pb,i/(length(unique(yrmn))))
    }
    y <- x.tracked
  } else {
    x.tracked <- Track(x,x0=NULL,lplot=lplot,dE=dE,dN=dN,
                       amax=amax,dmax=dmax,dmin=dmin,nmax=nmax,nmin=nmin,
                       progress=progress,verbose=verbose)
    y <- x.tracked$y
  }
  y <- attrcp(x,y)
  class(y) <- class(x)
  if(any(is.na(y$trajectory))) {
    nok <- is.na(y$trajectory)
    y$trajectory[nok] <- seq(sum(nok)) + max(y$trajectory,na.rm=TRUE)
  }
  invisible(y)
}

Track <- function(x,x0=NULL,it=NULL,is=NULL,dmax=1.2E6,amax=90,
                         nmax=124,nmin=3,dE=0.3,dN=0.2,dmin=1E5,
                         x0cleanup=TRUE,lplot=FALSE,
                         progress=TRUE,verbose=FALSE) {
  if (verbose) print("Track - cyclone tracking based on the distance and change in angle of direction between three subsequent time steps")
  options(digits=12)
  d <- sort(unique(strptime(paste(x$date,x$time),"%Y%m%d %H")))
  dh <- mean(as.numeric(difftime(d[2:length(d)],
                                 d[1:(length(d)-1)],units="hours")))
  #if (!("trajectory" %in% names(x))) {
  if(verbose) print("calculate trajectories")
  dates <- x$date
  times <- x$time
  lons <- x$lon
  lats <- x$lat
  num <- rep(NA,dim(x)[1])
  dx <- rep(NA,dim(x)[1])
  datetime <- strptime(paste(dates,times),"%Y%m%d %H")
  d <- sort(unique(datetime))
  if(!is.null(x0)) {
    d0 <- sort(unique(strptime(paste(x0$date,x0$time),"%Y%m%d %H")))
    dh0 <- as.numeric(difftime(min(d),max(d0),units="hours"))
  } else {
    dh0 <- 1E5
  }
  if (all(c("date","time","lon","lat","trajectory") %in% names(x0)) & dh0<dh*2) {
    dt0 <- x0$date*1E2 + x0$time
    num0 <- x0$trajectory
    x00 <- x0[dt0>=max(dt0),]
    dates <- c(x00$date,dates)
    times <- c(x00$time,times)
    lons <- c(x00$lon,lons)
    lats <- c(x00$lat,lats)
    num <- c(x00$trajectory,num)
    dx <- c(rep(0,dim(x00)[1]),dx)
    datetime <- strptime(paste(dates,times),"%Y%m%d %H")
    i.start <- length(x00["lon"][[1]])-1
    d <- unique(c(datetime[i.start:length(datetime)],d))
    n00 <- max(num0,na.rm=TRUE)
    nend0 <- unique(num0[dt0<max(dt0)])
  } else {
    x00 <- NULL
    num[datetime==d[1]] <- 1:sum(datetime==d[1])
    i.start <- 1
    n00 <- 0
    nend0 <- 0
  }
  t1 <- Sys.time()
  if (progress) pb <- txtProgressBar(style=3)
  nend <- nend0
  n0 <- n00
  for (i in 2:(length(d)-1)) {
    if (progress) setTxtProgressBar(pb,i/(length(d)-1))
    nn <- Track123(
        step1=list(lon=lons[datetime==d[i-1]],lat=lats[datetime==d[i-1]],
                   num=num[datetime==d[i-1]],dx=dx[datetime==d[i-1]]),
        step2=list(lon=lons[datetime==d[i]],lat=lats[datetime==d[i]],
                   num=num[datetime==d[i]],dx=dx[datetime==d[i]]),
        step3=list(lon=lons[datetime==d[i+1]],lat=lats[datetime==d[i+1]],
                   num=num[datetime==d[i+1]],dx=dx[datetime==d[i+1]]),
             dmax=dmax,n0=n0,amax=amax,nend=nend,dE=dE,dN=dN)
    num[datetime==d[i-1]] <- nn$step1$num
    num[datetime==d[i]] <- nn$step2$num
    num[datetime==d[i+1]] <- nn$step3$num
    dx[datetime==d[i]] <- nn$step2$dx
    dx[datetime==d[i+1]] <- nn$step3$dx
    n0 <- max(nn$n0,n0,na.rm=TRUE)
    nend <- nn$nend
  }
  t2 <- Sys.time()
  if (verbose) print(paste('Three step tracking took',
            round(as.numeric(t2-t1,units="secs")),'s'))

  tracklen <- data.frame(table(num))
  if (!is.null(nmax)) {
    while (any(tracklen$Freq>nmax)) {
      for (n in tracklen[tracklen$Freq>nmax,]$num) {
        i <- which(num==n)
        dn <- mapply(distAB,lons[i][2:length(i)],lats[i][2:length(i)],
                   lons[i][1:(length(i)-1)],lats[i][1:(length(i)-1)])
        dn[is.na(dn)] <- 0
        num[!is.na(num) & num>as.numeric(n)] <- num[!is.na(num) & num>as.numeric(n)] + 1
        num[i[which(dn==max(dn))+1:length(i)]] <- as.numeric(n) + 1
        dx[i[which(dn==max(dn))+1]]<- 0
      }
      tracklen <- data.frame(table(num))
    }
  }
  dx[is.na(dx)] <- 0
  i.start <- length(num)-length(x$date)+1
  x$trajectory <- num[i.start:length(num)]
  x$distance <- dx[i.start:length(num)]*1E-3
  if ( length(names(x))>length(attr(x,"unit")) ) {
    attr(x,"unit") <- c(attr(x,"unit"),c("numeration","km"))
  }

  if (verbose) print("calculate trajectory statistics")
  x <- Trackstats(x)
  if(!is.null(x00)) {
    x00$trajectory <- num[1:(i.start-1)]
    x00$distance <- dx[1:(i.start-1)]*1E-3
    x0[dt0>=max(dt0),] <- x00
    x01 <- rbind(x0,x)
    c01 <- attrcp(x,x01)
    class(x01) <- class(x)
    x01 <- Trackstats(x01)
    dnum <- x01$date*1E2 + x01$time
    if (is.null(nmin)) nmin <- 3
    if (is.null(dmin)) dmin <- 0
    starts <- dnum < as.numeric( strftime(strptime(min(dnum),"%Y%m%d%H") +
              (nmin-1)*dh*3600,"%Y%m%d%H") )
    ends <- dnum > as.numeric( strftime(strptime(max(dnum),"%Y%m%d%H") -
              (nmin-1)*dh*3600,"%Y%m%d%H") )
    ok <- x01$trackcount>=nmin | ends | starts
    ok <- ok & (x01$tracklength>=(dmin*1E-3) | ends | starts)
    if(!x0cleanup) ok[dnum<min(x$date*1E2 + x$time)]
    y01 <- subset(x01,it=ok)
    rnum <- y01$trajectory
    for(k in unique(x01$trajectory[!ok])) {
      rnum[y01$trajectory>k & !is.na(rnum)] <- rnum[y01$trajectory>k & !is.na(rnum)]-1
    }
    dy <- Displacement(y01)
    y01$distance <- dy
    dnum <- y01$date*1E2 + y01$time
    dnum1 <- x$date*1E2 + x$time
    y0 <- subset(y01,it=dnum<min(dnum1))
    y <- subset(y01,it=dnum>=min(dnum1))
  } else {
    dnum <- x["date"][[1]]*1E2 + x["time"][[1]]
    ok <- rep(TRUE,length(x$trackcount))
    ends <- dnum < as.numeric( strftime(strptime(min(dnum),"%Y%m%d%H") +
                   (nmin-1)*dh*3600,"%Y%m%d%H") ) |
            dnum > as.numeric( strftime(strptime(max(dnum),"%Y%m%d%H") -
                   (nmin-1)*dh*3600,"%Y%m%d%H") )
    ok <- x$trackcount>=nmin | ends
    ok <- ok & (x$tracklength>=(dmin*1E-3) | ends)
    y <- subset(x,it=ok,verbose=verbose)
    rnum <- Enumerate(y,verbose=verbose)
    rnum[y$trajectory<=n00 & !is.na(rnum)] <- y$trajectory[y$trajectory<=n00 & !is.na(rnum)]
    rnum[y$trajectory>n00 & !is.na(rnum)] <- n00 + rnum[y$trajectory>n00 & !is.na(rnum)]
    y$trajectory <- rnum
    y0 <- NULL
  }
 
  if(lplot) {
    lons <- y$lon
    lats <- y$lat
    num <- y$trajectory
    dates <- y$date
    times <- y$time
    dplot <- unique(dates)[min(5,length(unique(dates)))]
    tplot <- unique(times[dates==dplot])[min(3,
                    length(unique(times[dates==dplot])))]
    nvec <- unique(num[dates==dplot & times==tplot])
          #tracklen$num[tracklen$Freq>q995(tracklen$Freq)]
    cols <- rainbow(length(nvec))
    data(geoborders,envir=environment())
    plot(geoborders,type="l",col="grey20",lwd=0.5,
         xlim=c(-180,180),ylim=c(0,90),
         main=paste(dplot,tplot))
    for (i in seq_along(nvec)) {
      points(lons[num==nvec[i]],lats[num==nvec[i]],col=cols[i],pch=1,cex=1.5)
      points(lons[num==nvec[i]][1],lats[num==nvec[i]][1],
             col=cols[i],pch=".",cex=3)
      points(lons[num==nvec[i] & dates==dplot & times==tplot][1],
           lats[num==nvec[i] & dates==dplot & times==tplot][1],
           col=adjustcolor(cols[i],alpha=0.5),pch=19,cex=1.5)
    }
  }
  invisible(list(y=y,y0=y0))
}

Track123 <- function(step1,step2,step3,n0=0,amax=90,dmax=1.2E6,
                     dE=0.3,dN=0.2,dmax.s=2E5,nend=NA,lplot=FALSE,
                     verbose=FALSE) {
  if (verbose) print("Three step cyclone tracking")
  if (is.na(n0) & !all(is.na(step1$num))) {
    n0 <- max(step1$num,na.rm=TRUE)
  } else if (is.na(n0) & all(is.na(step1$num))) {
    n0 <- 0
  }
  if (any(is.na(step1$num))) {
    step1$num[is.na(step1$num)] <- seq(sum(is.na(step1$num))) + n0
  }
  n0 <- max(c(n0,step1$num),na.rm=TRUE)
  n1 <- length(step1$lon)
  n2 <- length(step2$lon)
  n3 <- length(step3$lon)
  d23 <- mapply(function(x,y) distAB(x,y,step3$lon,step3$lat),step2$lon,step2$lat)
  d23[is.na(d23)] <- 0
  a23 <- mapply(function(x,y) angle(x,y,step3$lon,step3$lat),step2$lon,step2$lat)
  dmax23 <- adjustdmax(a23,dmax=dmax,dE=dE,dN=dN)
  d12 <- mapply(function(x,y) distAB(x,y,step2$lon,step2$lat),step1$lon,step1$lat)
  d12[is.na(d12)] <- 0
  a12 <- mapply(function(x,y) angle(x,y,step2$lon,step2$lat),step1$lon,step1$lat)
  dmax12 <- adjustdmax(a12,dmax=dmax,dE=dE,dN=dN)
  da <- sapply(a12,function(x) abs(x-a23))
  da[da>180 & !is.na(da)] <- abs(da[da>180 & !is.na(da)]-360)
  dd <- sapply(d12,function(x) abs(x-d23))
  d123 <- sapply(d12,function(x) x+d23)
  ok.d <- sapply(d12<dmax12,function(x) sapply(d23<dmax23,function(y) y & x ))
  ok.d2 <- sapply(d12<dmax.s,function(x) sapply(d23<dmax.s,function(y) y & x ))
  ok.dd <- (dd/d123 < 0.5 | ok.d2)
  ok <- ok.d & (da <= amax | ok.d2) & ok.dd
  j1 <- as.vector(sapply(seq(n1),function(x) rep(x,n2)))
  j2 <- rep(seq(n2),n1)
  i2 <- as.vector(sapply(seq(n2),function(x) rep(x,n3)))
  i3 <- rep(seq(n3),n2)
  for(i in unique(i2)) ok[i==i2,i!=j2] <- FALSE
  if(any(ok)) {
    da[!ok] <- NA
    dim(da) <- c(n2*n3,n1*n2)
    dim(d12) <- c(n2,n1)
    dim(d23) <- c(n3,n2)
    nok <- matrix(rep(FALSE,length(da)),dim(da))
    if(any(!is.na(step2$num))) {
      for(k in unique(step2$num[!is.na(step2$num)])) {
        nok[,(j1==which(step1$num==k) & j2!=which(step2$num==k))] <- TRUE
        nok[,(j1!=which(step1$num==k) & j2==which(step2$num==k))] <- TRUE
        nok[(i2==which(step2$num==k)),j1!=which(step1$num==k)] <- TRUE    
      }
    }
    if(any(!is.na(nend))) {
      for(k in nend[!is.na(nend) & nend %in% step1$num]) {
        nok[,(j1==which(step1$num==k))] <- TRUE
      }
      for(k in nend[!is.na(nend) & nend %in% step2$num]) {
        nok[,(j2==which(step2$num==k))] <- TRUE
        nok[(i2==which(step2$num==k)),] <- TRUE
      }
    }
    da[nok] <- NA
    dd[is.na(da)] <- NA
    d123[is.na(da)] <- NA
    rank.da <- matrix(rank(da),dim(da))
    rank.dd <- matrix(rank(dd),dim(dd))
    rank.d123 <- matrix(rank(d123),dim(d123))
    rank.all <- rank.da + rank.dd + rank.d123
    rank.all[is.na(da)] <- NA
    while(any(!is.na(rank.all))) {
      ij <- which(rank.all==min(rank.all,na.rm=TRUE),arr.ind=TRUE)
      if(dim(ij)[1]>1) {
        k <- which.min(apply(ij,1,function(x) da[x[1],x[2]]))
        ij <- ij[k,]
      }
      step2$num[j2[ij[2]]] <- step1$num[j1[ij[2]]]
      step3$num[i3[ij[1]]] <- step1$num[j1[ij[2]]]
      step2$dx[j2[ij[2]]] <- d12[j2[ij[2]],j1[ij[2]]]
      step3$dx[i3[ij[1]]] <- d23[i3[ij[1]],j2[ij[2]]]
      rank.all[,j1==j1[ij[2]]] <- NA
      rank.all[,j2==j2[ij[2]]] <- NA
      rank.all[i2==i2[ij[1]],] <- NA
      rank.all[i3==i3[ij[1]],] <- NA
    }
  }
  if(all(is.na(step2$num))) {
    step2$num[is.na(step2$num)] <- seq(sum(is.na(step2$num))) + n0
  } else if (any(is.na(step2$num))) {
    step2$num[is.na(step2$num)] <- seq(sum(is.na(step2$num))) + max(c(n0,step2$num),na.rm=TRUE)    
  }
  nok <- (step2$num %in% step1$num) & !(step2$num %in% step3$num) & !is.na(step2$num)
  if(any(nok)) nend <- c(nend,step2$num[nok]); nend <- nend[!is.na(nend)]
  return(list(step1=step1,step2=step2,step3=step3,nend=nend,
         n0=max(c(step1$num,step2$num,step3$num,n0),na.rm=TRUE)))
}

angle <- function(lon1,lat1,lon2,lat2) {
  a <- 360 - (atan2(lat2-lat1,lon2-lon1)*(180/pi) + 360) %% 360
  #a[a>180] <- 360-a[a>180]
  return(a)
  #return(atan2(lat2-lat1,lon2-lon1)*180/pi+90)
}

## adjust maximum distance based on angle of direction: max eastward, min westward 
adjustdmax <- function(a,dmax=1.2E6,dE=0.3,dN=0,width=1,height=1,lplot=FALSE) {
  rad <- a*pi/180
  east <- rad > -pi/2 & rad < pi/2
  north <- rad < pi & rad > 0
  x <- width*height/sqrt(height^2 + width^2*(tan(rad)^2))
  x[!east] <- -x[!east]
  y <- height*sqrt(1-(x/width)^2)
  y[!north] <- -y[!north]
  #x <- (x + dE)/sqrt((cos(atan(dN/dE)) + dE)^2 + (sin(atan(dN/dE)) + dN)^2)
  #y <- (y + dN)/sqrt((cos(atan(dN/dE)) + dE)^2 + (sin(atan(dN/dE)) + dN)^2)
  x <- x + dE
  y <- y + dN
  d <- dmax*sqrt(x^2 + y^2)
  if(lplot) {
    plot(0,0,cex=2,pch=3,col="black",
         xlim=c(-1.5,1.5)*dmax*1E-3,ylim=c(-1.5,1.5)*dmax*1E-3,
         xlab="dmax (km)",ylab="dmax (km)",main="maximum displacement radius")
    points(dmax*x*1E-3,dmax*y*1E-3,pch=19,cex=0.5)
    grid()
  }
  return(d)
}

Trackstats <- function(x,verbose=FALSE) {
  if(verbose) print("Trackstats")
  if (!any("trajectory" %in% names(x))) x <- tracking.events(x,verbose=verbose)
  y <- x[order(x$trajectory),]
  y <- attrcp(x,y)
  class(y) <- class(x)
  lons <- y["lon"][[1]]
  lats <- y["lat"][[1]]
  if(verbose) print("Enumerate")
  rnum <- Enumerate(y)
  nummax <- max(rnum,na.rm=TRUE)
  rnum[is.na(rnum)] <- nummax+1
  if(verbose) print("displacement")
  dy <- Displacement(y)
  y$distance <- dy
  if(!"distance" %in% colnames(x)) {
    attr(y,"unit") <- c(attr(y,"unit"),"km")
  }
  if(verbose) print("trackcount")
  trackcount <- data.frame(table(rnum))
  if(verbose) print("timestep")
  ts <- unlist(sapply(unique(rnum),function(i) 1:trackcount$Freq[i]))
  timestep <- rep(NA,length(ts))
  timestep[order(rnum)] <- ts
  if(verbose) print("tracklength")
  tracklength <- as.numeric(by(y$distance,rnum,function(x) sum(x,na.rm=TRUE)))
  if (any(rnum>nummax)) {
    trackcount$Freq[nummax+1] <- 1
    timestep[rnum==(nummax+1)] <- 1
    tracklength[nummax+1] <- 0
  }
  y$trackcount <- trackcount$Freq[rnum]
  y$timestep <- timestep
  y$tracklength <- tracklength[rnum]
  if (length(attr(y,"unit"))<dim(y)[2]) {
    attr(y,"unit") <- c(attr(y,"unit"),c("count","numeration","km"))
  }
  invisible(y)
}

Displacement <- function(x,verbose=FALSE) {
  if(verbose) print("Calculate displacement")
  if(!"trajectory" %in% colnames(x)) {
    if(verbose) print("track.events")
    x <- track.events(x)
  }
  if("distance" %in% colnames(x)) {
    dx <- x$distance
  } else {
    dx <- rep(0,dim(x)[1])
  }
  for(k in unique(x$trajectory[x$trackcount>1 & x$timestep>1 & dx==0])) {
    ik <- x$trajectory==k & !is.na(x$trajectory)
    if(sum(ik)>1) {
      dk <- mapply(distAB,x$lon[ik][1:(sum(ik)-1)],x$lat[ik][1:(sum(ik)-1)],
                          x$lon[ik][2:sum(ik)],x$lat[ik][2:sum(ik)])
      dk[is.na(dk)] <- 0
      if(length(dk)<sum(ik)) dk <- c(0,dk)
      dx[ik] <- dk*1E-3
    }
  }
  invisible(dx)
}
    
Enumerate <- function(x,param="trajectory",verbose=FALSE) {
  if(verbose) print("Enumerate")
  stopifnot(inherits(x,"data.frame"))
  num <- x[param][[1]]
  if (identical(unique(num),seq_along(unique(num)))) {
    rnum <- num
  } else {
    if(verbose) print(paste("number of events:",length(num)))
    if(verbose) print(paste("number of trajectories:",length(unique(num))))
    rnum <- num[2:length(num)]-num[1:(length(num)-1)]
    rnum[is.na(rnum)] <- 0
    rnum[rnum>0] <- 2:(sum(rnum>0)+1)
    rnum <- c(1,rnum)
    while(any(rnum==0)) {
      rnum[rnum==0] <- rnum[which(rnum==0)-1]
    }
    rnum[is.na(num)] <- NA
  }
  invisible(rnum)
}



## NearestNeighbour <- function(lon1,lat1,lon2,lat2,dmax=1E6,
##                              lplot=FALSE,verbose=FALSE) {
##   if (verbose) print("NearestNeighbour")
##   distance <- mapply(function(x,y) distAB(x,y,lon2,lat2),lon1,lat1)
##   n <- length(lon1)
##   if (n==1) {
##     num <- as.numeric(which(rank(distance)==1))
##   } else {
##     num <- as.numeric(apply(distance,1,function(x) which(rank(x)==1)))
##   }
##   for (i in which(is.na(num))) {
##     num.i <- which(rank(distance[,i])==min(rank(distance[,i])))
##     if(!all(num.i %in% num[-i]) & !all(!num.i %in% num[-i])) {
##       num[i] <- num.i[!num.i %in% num[-i]][1]
##     } else {
##       num[i] <- num.i[1]
##     }
##   }
##   if(length(num)==1) {
##     d.num <- distance[num]
##   } else {
##     d.num <- sapply(1:length(num),function(x) distance[x,num[x]])
##   }
##   if (any(d.num>dmax)) {
##     num[d.num>dmax] <- NA
##   }
##   while (any(duplicated(num[!is.na(num)]))) {
##     for (i in num[duplicated(num) & !is.na(num)]) {
##       j <- which(num==i)
##       k <- j[d.num[j]==max(d.num[j])][1]
##       num.k <- which(rank(distance[,k],ties.method="first")==2)
##       if (!any(num==num.k,na.rm=TRUE) & distance[num.k,k]<=dmax) {
##         num[k] <- num.k
##       } else {
##         num[k] <- NA
##       }
##     }
##   }
##   d.num[is.na(num)] <- 0
##   d.num <- d.num*1E-3
##   if (lplot) {
##     plot(lon1,lat1,type="p",col="black",pch=1:length(lon1),lwd=1,
##        xlim=c(-180,180),ylim=c(0,90))
##     points(lon2,lat2,col="blue",pch=num,lwd=1)
##   }
##   invisible(cbind(num,d.num))
## }
