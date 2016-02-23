# K Parding, 15.10.2015

track.events <- function(x,x0=NULL,it=NULL,is=NULL,dmax=8E4,amax=90,
                         nmax=31*24,nmin=5,dmin=5E5,lplot=FALSE,
                         progress=TRUE,verbose=FALSE) {

  if(verbose) print("tracking.events")
  x <- subset(x,it=!is.na(x["date"][[1]]))
  x <- subset(x,it=it,is=is)
  yrmn <- as.yearmon(strptime(x["date"][[1]],"%Y%m%d"))
  if (length(unique(yrmn))>1) {
    x.tracked <- NULL
    if (progress) pb <- txtProgressBar(style=3)
    for (i in 1:length(unique(yrmn))) {
      if(verbose) print(unique(yrmn)[i])
      if (progress) setTxtProgressBar(pb,i/(length(unique(yrmn))))
      x.y <- subset(x,it=(yrmn==unique(yrmn)[i]))
      if (is.null(x.tracked)) {
        x.t <- Track(x.y,x0=x0,lplot=lplot,x0cleanup=FALSE,
                     amax=amax,dmax=dmax,dmin=dmin,nmax=nmax,nmin=nmin,
                     progress=FALSE,verbose=verbose)
        x.tracked <- x.t$y
      } else {
        x.t <- Track(x.y,x0=x.tracked,lplot=lplot,
                     amax=amax,dmax=dmax,dmin=dmin,nmax=nmax,nmin=nmin,
                     progress=FALSE,verbose=verbose)
        x.tracked <- x.t$y0
        x.tracked <- merge(x.tracked,x.t$y,all=TRUE)
      }
    }
    y <- x.tracked
  } else {
    x.tracked <- Track(x,x0=NULL,lplot=lplot,
                       amax=amax,dmax=dmax,nmax=nmax,nmin=nmin,
                       progress=progress,verbose=verbose)
    y <- x.tracked$y
  }
  y <- attrcp(x,y)
  class(y) <- class(x)
  if(any(is.na(y$trajectory))) {
    nmax <- max(y$trajectory,na.rm=TRUE)
    nok <- is.na(y$trajectory)
    y$trajectory[nok] <- seq(sum(nok)) + nmax
  }
  invisible(y)
}

Track <- function(x,x0=NULL,it=NULL,is=NULL,dmax=8E4,amax=map90,
                         nmax=31*24,nmin=5,dmin=5E5,
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
             dmax=dmax*dh,n0=n0,amax=amax,nend=nend)
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
  #} else {
  #  num <- x["trajectory"][[1]]
  #  if("distance" %in% colnames(x)) {
  #    dx <- x["distance"][[1]]
  #  } else {
  #    dx <- rep(NA,length(x["trajectory"][[1]]))
  #  }
  #  i.start <- 1
  #}

  tracklen <- data.frame(table(num))
  if (!is.null(nmax)) {
    while (any(tracklen$Freq>nmax/dh)) {
      for (n in tracklen[tracklen$Freq>nmax/dh,]$num) {
        i <- which(num==n)
        dn <- mapply(distAB,lons[i][2:length(i)],lats[i][2:length(i)],
                   lons[i][1:(length(i)-1)],lats[i][1:(length(i)-1)])
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

Track123 <- function(step1,step2,step3,n0=0,amax=90,dmin=0,dmax=1E6,
                             nend=NA,lplot=FALSE,verbose=FALSE) {
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
  dmax23 <- adjustdmax(dmax,a23)
  a23[d23>dmax23] <- NA
  a23[d23<dmin] <- NA
  d12 <- mapply(function(x,y) distAB(x,y,step2$lon,step2$lat),step1$lon,step1$lat)
  d12[is.na(d12)] <- 0
  a12 <- mapply(function(x,y) angle(x,y,step2$lon,step2$lat),step1$lon,step1$lat)
  dmax12 <- adjustdmax(dmax,a12)
  a12[d12>dmax12] <- NA
  a12[d12<dmin] <- NA
  if(!all(is.na(a12)) & !all(is.na(a23))) {
    da <- sapply(a12,function(x) abs(x-a23))
    da[da > amax & !is.na(da)] <- NA
    dim(da) <- c(n2*n3,n1*n2)
    dim(d12) <- c(n2,n1)
    dim(d23) <- c(n3,n2)
    j1 <- as.vector(sapply(seq(n1),function(x) rep(x,n2)))
    j2 <- rep(seq(n2),n1)
    i2 <- as.vector(sapply(seq(n2),function(x) rep(x,n3)))
    i3 <- rep(seq(n3),n2)
    if(any(!is.na(step2$num))) {
      kvec <- unique(c(step1$num,step2$num))
      kvec <- kvec[!is.na(kvec)]
      for(k in kvec) {
        if(k %in% step1$num  & k %in% step2$num) {
          nok <- (j1==which(step1$num==k) & !j2==which(step2$num==k)) |
                 (j1!=which(step1$num==k) & j2==which(step2$num==k))
        } else if (k %in% step1$num) {
          nok <- j1==which(step1$num==k)
        } else {
          nok <- j2==which(step2$num==k)
        }
        da[,nok] <- NA
      }
    }
    if(!is.na(nend)) {
      for(k in nend) {
        if (k %in% c(step1$num,step2$num)) {
          nok <- j1==which(step1$num==k | step2$num==k)
          da[,nok] <- NA
        }
      }
    }
    for(i in unique(i2)) da[i==i2,i!=j2] <- NA
    rank.j1 <- matrix(rep(NA,n1*n2*n2*n3),c(n2*n3,n1*n2))
    rank.j2 <- matrix(rep(NA,n1*n2*n2*n3),c(n2*n3,n1*n2))
    rank.i2 <- matrix(rep(NA,n1*n2*n2*n3),c(n2*n3,n1*n2))
    rank.i3 <- matrix(rep(NA,n1*n2*n2*n3),c(n2*n3,n1*n2))
    for(j in seq(n1)) rank.j1[,which(j1==j)] <- rank(da[,which(j1==j)])
    for(j in seq(n2)) rank.j2[,which(j2==j)] <- rank(da[,which(j2==j)])
    for(i in seq(n2)) rank.i2[which(i2==i),] <- rank(da[which(i2==i),])
    for(i in seq(n3)) rank.i3[which(i3==i),] <- rank(da[which(i3==i),])
    rank.j1[is.na(da)] <- NA
    rank.j2[is.na(da)] <- NA
    rank.i2[is.na(da)] <- NA
    rank.i3[is.na(da)] <- NA
    rank.all <- rank.j1 + rank.j2 + rank.i2 + rank.i3
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
  a[a>180] <- 360-a[a>180]
  return(a)
  #return(atan2(lat2-lat1,lon2-lon1)*180/pi+90)
}

## adjust maximum distance based on angle of direction: max eastward, min westward 
adjustdmax <- function(dmax,a) {
  a[a>180] <- 360-a
  a <- abs(a)
  return(dmax*(1.3-0.2*(a/180)-0.8*(a/180)**2))
}

Trackstats <- function(x,verbose=FALSE) {
  if(verbose) print("Trackstats")
  if (!any("trajectory" %in% names(x))) x <- tracking.events(x,verbose=verbose)
  y <- x[order(x$trajectory),]
  y <- attrcp(x,y)
  class(y) <- class(x)
  lons <- y["lon"][[1]]
  lats <- y["lat"][[1]]
  rnum <- Enumerate(y)
  nummax <- max(rnum,na.rm=TRUE)
  rnum[is.na(rnum)] <- nummax+1
  if(verbose) print("trackcount")
  trackcount <- data.frame(table(rnum))
  if(verbose) print("timestep")
  ts <- unlist(sapply(unique(rnum),function(i) 1:trackcount$Freq[i]))
  timestep <- rep(NA,length(ts))
  timestep[order(rnum)] <- ts
  if(verbose) print("tracklength")
  lon1 <- as.numeric(by(y$lon,rnum,function(x) x[1]))
  lat1 <- as.numeric(by(y$lat,rnum,function(x) x[1]))
  lon2 <- as.numeric(by(y$lon,rnum,function(x) x[length(x)]))
  lat2 <- as.numeric(by(y$lat,rnum,function(x) x[length(x)]))
  tracklength <- round(mapply(distAB,lon1,lat1,lon2,lat2)*1E-3)
  tracklength[lat1==lat2 & lon1==lon2] <- 0
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
  if(verbose) print("displacement")
  dy <- Displacement(y)
  y$distance <- dy
  if(!"distance" %in% colnames(x)) {
    attr(y,"unit") <- c(attr(y,"unit"),"km")
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
