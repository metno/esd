# K Parding, 15.10.2015

Track.events <- function(x,x0=NULL,it=NULL,is=NULL,dmax=8E4,nmax=31*24,
                         nmin=5,dmin=5E5,lplot=FALSE,progress=TRUE,
                         verbose=FALSE) {
  if (verbose) print("Track.events - nearest neighbour cyclone tracking")
  x <- subset(x,it=!is.na(x["date"][[1]]))
  x <- subset(x,it=it,is=is)
  yrmn <- as.yearmon(strptime(x["date"][[1]],"%Y%m%d"))
  d <- sort(unique(strptime(paste(x["date"][[1]],x["time"][[1]]),"%Y%m%d %H")))
  dh <- mean(as.numeric(difftime(d[2:length(d)],
                                 d[1:(length(d)-1)],units="hours")))
  if (length(unique(yrmn))>1) {
    x.tracked <- NULL
    if (progress) pb <- txtProgressBar(style=3)
    for (i in 1:length(unique(yrmn))) {
      if(verbose) print(unique(yrmn)[i])
      if (progress) setTxtProgressBar(pb,i/(length(unique(yrmn))))
      x.y <- subset(x,it=(yrmn==unique(yrmn)[i]))
      if (is.null(x.tracked)) {
        x.t <- Track.events(x.y,x0=x0,lplot=lplot,
                        progress=FALSE,verbose=verbose)
        x.tracked <- x.t
      } else {
        x.t <- Track.events(x.y,x0=x.tracked,lplot=lplot,
                        progress=FALSE,verbose=verbose)         
        x.tracked <- merge(x.tracked,x.t,all=TRUE)
      }
    }
    x <- x.tracked
    x <- attrcp(x.t,x)
    class(x) <- class(x.t)
    invisible(x)
  } else { 
    if (!("trajectory" %in% names(x))) {
      if(verbose) print("calculate trajectories")
      dates <- x["date"][[1]]
      times <- x["time"][[1]]
      lons <- x["lon"][[1]]
      lats <- x["lat"][[1]]
      num <- rep(NA,dim(x)[1])
      dx <- rep(NA,dim(x)[1])
      datetime <- strptime(paste(dates,times),"%Y%m%d %H")
      d <- sort(unique(datetime))
      if (all(c("date","time","lon","lat","trajectory") %in% names(x0))) {
        dt0 <- strptime(paste(x0["date"][[1]],x0["time"][[1]]),"%Y%m%d %H")
        num0 <- x0["trajectory"][[1]]
        x0 <- x0[num0 %in% num0[dt0==max(dt0[dt0<min(datetime)])],]
        dates <- c(x0["date"][[1]],dates)
        times <- c(x0["time"][[1]],times)
        lons <- c(x0["lon"][[1]],lons)
        lats <- c(x0["lat"][[1]],lats)
        num <- c(x0["trajectory"][[1]],num)
        dx <- c(rep(0,dim(x0)[1]),dx)
        datetime <- strptime(paste(dates,times),"%Y%m%d %H")
        i.start <- length(x0["lon"][[1]])+1
        d <- c(max(datetime[1:i.start-1]),d)
      } else {
        num[datetime==d[1]] <- 1:sum(datetime==d[1])
        i.start <- 1
      }
      t1 <- Sys.time()
      if (progress) pb <- txtProgressBar(style=3)
      nnmax <- 0
      for (i in 2:(length(d)-1)) {
        if (verbose) print(i)
        if (progress) setTxtProgressBar(pb,i/(length(d)-1))
        nn <- NearestNeighbour(
            step1=list(lon=lons[datetime==d[i-1]],lat=lats[datetime==d[i-1]],
                       num=num[datetime==d[i-1]],dx=dx[datetime==d[i-1]]),
            step2=list(lon=lons[datetime==d[i]],lat=lats[datetime==d[i]],
                       num=num[datetime==d[i]],dx=dx[datetime==d[i]]),
            step3=list(lon=lons[datetime==d[i+1]],lat=lats[datetime==d[i+1]],
                       num=num[datetime==d[i+1]],dx=dx[datetime==d[i+1]]),
                 dmax=dmax*dh,nnmax=nnmax)
        num[datetime==d[i-1]] <- nn$step1$num
        num[datetime==d[i]] <- nn$step2$num
        num[datetime==d[i+1]] <- nn$step3$num
        dx[datetime==d[i]] <- nn$step2$dx
        dx[datetime==d[i+1]] <- nn$step3$dx
        nnmax <- max(c(nnmax,nn$nnmax),na.rm=TRUE)
      }
      t2 <- Sys.time()
      if (verbose) print(paste('Nearest neighbour tracking took',
                round(as.numeric(t2-t1,units="secs")),'s'))
    } else {
      num <- x["trajectory"][[1]]
      if("distance" %in% colnames(x)) {
        dx <- x["distance"][[1]]
      } else {
        dx <- rep(NA,length(x["trajectory"][[1]]))
      }
      i.start <- 1
    }
    browser()
    tracklen <- data.frame(table(num))
    if (!is.null(nmax)) {
      while (any(tracklen$Freq>nmax/dh)) {
        for (n in tracklen[tracklen$Freq>nmax/dh,]$num) {
          i <- which(num==n)
          dn <- mapply(distAB,lons[i][2:length(i)],lats[i][2:length(i)],
                     lons[i][1:(length(i)-1)],lats[i][1:(length(i)-1)])
          num[num>as.numeric(n)] <- num[num>as.numeric(n)] + 1
          num[i[which(dn==max(dn))+1:length(i)]] <- as.numeric(n) + 1
          dx[i[which(dn==max(dn))+1]]<- 0
        }
        tracklen <- data.frame(table(num))
      }
    }
    dx[is.na(dx)] <- 0
    x$trajectory <- num[i.start:length(num)]
    x$distance <- dx[i.start:length(num)]
    if ( length(names(x))>length(attr(x,"unit")) ) {
      attr(x,"unit") <- c(attr(x,"unit"),c("numeration","km"))
    }
    if (verbose) print("calculate trajectory statistics")
    x <- Trackstats(x)
    num <- x["trajectory"][[1]]
    n <- x["trackcount"][[1]]
    d <- x["tracklength"][[1]]
    options(digits=12)
    dnum <- as.numeric(paste(x["date"][[1]],x["time"][[1]],sep="."))
    ok <- rep(TRUE,length(n))
    ends <- dnum==min(dnum) | dnum==max(dnum)
    if (!is.null(nmin)) ok <- n>=nmin | (num %in% num[ends])
    if (!is.null(dmin)) ok <- ok & (d>=(dmin*1E-3) | (num %in% num[ends]))
    y <- subset(x,it=ok)
    rnum <- Enumerate(y,verbose=verbose)
    if (i.start>1) rnum <- rnum + max(num[1:(i.start-1)])
    y["trajectory"] <- rnum
    if(lplot) {
      lons <- y["lon"][[1]]
      lats <- y["lat"][[1]]
      num <- y["trajectory"][[1]]
      dates <- y["date"][[1]]
      times <- y["time"][[1]]
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
      #lines(lons[num==i],lats[num==i],col="blue",lty=1,lwd=1)
      }
    }
  invisible(y)
  }
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

NearestNeighbour <- function(step1,step2,step3,nnmax=0,damax=90,dmax=1E6,lplot=FALSE,verbose=FALSE) {
  if (verbose) print("NearestNeighbour")
  if (is.na(nnmax) & !all(is.na(step1$num))) {
    nnmax <- max(step1$num,na.rm=TRUE)
  } else if (is.na(nnmax) & all(is.na(step1$num))) {
    nnmax <- 0
  }
  if (any(is.na(step1$num))) {
    step1$num[is.na(step1$num)] <- seq(sum(is.na(step1$num))) + nnmax
  }
  n1 <- length(step1$lon)
  n2 <- length(step2$lon)
  n3 <- length(step3$lon)
  d23 <- mapply(function(x,y) distAB(x,y,step3$lon,step3$lat),step2$lon,step2$lat)
  a23 <- mapply(function(x,y) angle(x,y,step3$lon,step3$lat),step2$lon,step2$lat)
  dmax23 <- adjustdmax(dmax,a23)
  a23[d23>dmax23] <- NA
  d12 <- mapply(function(x,y) distAB(x,y,step2$lon,step2$lat),step1$lon,step1$lat)
  a12 <- mapply(function(x,y) angle(x,y,step2$lon,step2$lat),step1$lon,step1$lat)
  dmax12 <- adjustdmax(dmax,a12)
  a12[d12>dmax12] <- NA
  if(!all(is.na(a12)) & !all(is.na(a23))) {
    da <- sapply(a12,function(x) abs(x-a23))
    da[da > damax & !is.na(da)] <- NA
    j1 <- ceiling(seq(dim(da)[2])/n2)
    j2 <- rep(seq(n2),n1)
    i3 <- rep(seq(n3),n2)
    if(any(!is.na(step2$num))) {
      for(k in unique(c(step1$num,step2$num))) {
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
    rank.j1 <- matrix(rep(NA,length(da)),dim(da))
    rank.j2 <- matrix(rep(NA,length(da)),dim(da))
    rank.i3 <- matrix(rep(NA,length(da)),dim(da))
    for(j in unique(j1)) rank.j1[,j1==j] <- rank(da[,j1==j])
    for(j in unique(j2)) rank.j2[,j2==j] <- rank(da[,j2==j])
    for(i in unique(i3)) rank.i3[i3==i,] <- rank(da[i3==i,])
    rank.j1[is.na(da)] <- NA
    rank.j2[is.na(da)] <- NA
    rank.i3[is.na(da)] <- NA
    rank.all <- rank.j1 + rank.j2 + rank.i3
    while(any(!is.na(rank.all))) {
      ij <- which(rank.all==min(rank.all,na.rm=TRUE),arr.ind=TRUE)
      step2$num[j2[ij[2]]] <- step1$num[j1[ij[2]]]
      step3$num[i3[ij[1]]] <- step1$num[j1[ij[2]]]
      step2$dx[j2[ij[2]]] <- d12[j2[ij[2]],j1[ij[2]]]
      step3$dx[i3[ij[1]]] <- d23[i3[ij[1]],j2[ij[2]]]
      rank.all[,j1==j1[ij[2]]] <- NA
      rank.all[,j2==j2[ij[2]]] <- NA
      rank.all[i3==i3[ij[1]],] <- NA
    }
  }
  if(all(is.na(step2$num))) {
    step2$num[is.na(step2$num)] <- seq(sum(is.na(step2$num))) + max(step1$num,na.rm=TRUE)
  } else if (any(is.na(step2$num))) {
    step2$num[is.na(step2$num)] <- seq(sum(is.na(step2$num))) + max(step2$num,na.rm=TRUE)    
  }
  return(list(step1=step1,step2=step2,step3=step3,
              nnmax=max(c(step1$num,step2$num,step3$num),na.rm=TRUE)))
}

NearestNeighbour_org <- function(lon1,lat1,lon2,lat2,dmax=1E6,
                             lplot=FALSE,verbose=FALSE) {
  if (verbose) print("NearestNeighbour")
  distance <- mapply(function(x,y) distAB(x,y,lon2,lat2),lon1,lat1)
  n <- length(lon1)
  if (n==1) {
    num <- as.numeric(which(rank(distance)==1))
  } else {
    num <- as.numeric(apply(distance,1,function(x) which(rank(x)==1)))
  }
  for (i in which(is.na(num))) {
    num.i <- which(rank(distance[,i])==min(rank(distance[,i])))
    if(!all(num.i %in% num[-i]) & !all(!num.i %in% num[-i])) {
      num[i] <- num.i[!num.i %in% num[-i]][1]
    } else {
      num[i] <- num.i[1]
    }
  }
  if(length(num)==1) {
    d.num <- distance[num]
  } else {
    d.num <- sapply(1:length(num),function(x) distance[x,num[x]])
  }
  if (any(d.num>dmax)) {
    num[d.num>dmax] <- NA
  }
  while (any(duplicated(num[!is.na(num)]))) {
    for (i in num[duplicated(num) & !is.na(num)]) {
      j <- which(num==i)
      k <- j[d.num[j]==max(d.num[j])][1]
      num.k <- which(rank(distance[,k],ties.method="first")==2)
      if (!any(num==num.k,na.rm=TRUE) & distance[num.k,k]<=dmax) {
        num[k] <- num.k
      } else {
        num[k] <- NA
      }
    }
  }
  d.num[is.na(num)] <- 0
  d.num <- d.num*1E-3
  if (lplot) {
    plot(lon1,lat1,type="p",col="black",pch=1:length(lon1),lwd=1,
       xlim=c(-180,180),ylim=c(0,90))
    points(lon2,lat2,col="blue",pch=num,lwd=1)
  }
  invisible(cbind(num,d.num))
}

Trackstats <- function(x,verbose=FALSE) {
  if(verbose) print("Trackstats")
  if (!any("trajectory" %in% names(x))) x <- Track.events(x,verbose=verbose)
  y <- x[order(x$trajectory),]
  y <- attrcp(x,y)
  class(y) <- class(x)
  lons <- y["lon"][[1]]
  lats <- y["lat"][[1]]
  rnum <- Enumerate(y)
  if(verbose) print("trackcount")
  trackcount <- data.frame(table(rnum))
  if(verbose) print("timestep")
  ts <- unlist(sapply(unique(rnum),function(i) 1:trackcount$Freq[i]))
  timestep <- rep(NA,length(ts))
  timestep[order(rnum)] <- ts
  if(verbose) print("tracklength")
  #dy <- y["distance"][[1]]
  #dy[is.na(dy)] <- 0
  #tracklength <- as.numeric(by(dy,rnum,sum))
  lon1 <- as.numeric(by(y$lon,rnum,function(x) x[1]))
  lat1 <- as.numeric(by(y$lat,rnum,function(x) x[1]))
  lon2 <- as.numeric(by(y$lon,rnum,function(x) x[length(x)]))
  lat2 <- as.numeric(by(y$lat,rnum,function(x) x[length(x)]))
  tracklength <- round(mapply(distAB,lon1,lat1,lon2,lat2)*1E-3)
  tracklength[lat1==lat2 & lon1==lon2] <- 0
  y$trackcount <-  trackcount$Freq[rnum]
  y$timestep <- timestep
  y$tracklength <- tracklength[rnum]
  if (length(attr(y,"unit"))<dim(y)[2]) {
    attr(y,"unit") <- c(attr(y,"unit"),c("count","numeration","km"))
  }
  invisible(y)
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
