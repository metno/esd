# K Parding, 15.10.2015

Track.events <- function(x,x0=NULL,it=NULL,is=NULL,dmax=8E4,nmax=31*24,
                         nmin=5,dmin=5E5,lplot=FALSE,progress=TRUE,
                         verbose=FALSE) {
  if (verbose) print("Track.events - nearest neighbour cyclone tracking")
  x <- subset(x,it=!is.na(x["date"][[1]]))
  x <- subset(x,it=it,is=is)
  yrmn <- as.yearmon(strptime(x["date"][[1]],"%Y%m%d"))
  d <- sort(unique(strptime(paste(x["date"][[1]],x["time"][[1]]),"%Y%m%d %H")))
  dh <- mean(as.numeric(diff(d,units="hours")))
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
      for (i in 1:(length(d)-1)) {
        #if (verbose) print(i)
        if (progress) setTxtProgressBar(pb,i/(length(d)-1))
        nn <- NearestNeighbour(lons[datetime==d[i]],lats[datetime==d[i]],
                 lons[datetime==d[i+1]],lats[datetime==d[i+1]],dmax=dmax*dh)
        num.i <- num[datetime==d[i]][nn[,1]]
        num.i[is.na(num.i)] <- 1:sum(is.na(num.i)) + max(num,na.rm=TRUE)
        num[datetime==d[i+1]] <- num.i
        dx[datetime==d[i+1]] <- nn[,2]
      }
      t2 <- Sys.time()
      if (verbose) print(paste('Nearest neighbour tracking took',
                round(as.numeric(t2-t1,units="secs")),'s'))
    } else {
      num <- x["trajectory"][[1]]
      dx <- x["distance"][[1]]
      i.start <- 1
    }
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

NearestNeighbour <- function(lon1,lat1,lon2,lat2,dmax=1E6,lplot=FALSE,
                             verbose=FALSE) {
  if (verbose) print("NearestNeighbour")
  distance <- mapply(function(x,y) distAB(x,y,lon1,lat1),lon2,lat2)
  num <- as.numeric(apply(distance,2,function(x) which(rank(x)==1)))
  for (i in which(is.na(num))) {
    num.i <- which(rank(distance[,i])==min(rank(distance[,i])))
    if(!all(num.i %in% num[-i]) & !all(!num.i %in% num[-i])) {
      num[i] <- num.i[!num.i %in% num[-i]][1]
    } else {
      num[i] <- num.i[1]
    }
  }
  d.num <- sapply(1:length(num),function(x) distance[num[x],x])
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
