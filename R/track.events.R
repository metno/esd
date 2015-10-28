# K Parding, 15.10.2015

Track.events <- function(x,x0=NULL,it=NULL,is=NULL,dmax=8E4,nmax=31*24,
                         nmin=5,dmin=1E5,lplot=FALSE,progress=TRUE,
                         verbose=FALSE) {
  if (verbose) print("Track.events - nearest neighbour cyclone tracking")
  x <- subset(x,it=!is.na(x["date"][[1]]))
  x <- subset(x,it=it,is=is)
  if (!("trajectory" %in% names(x))) {
    yrmn <- as.yearmon(strptime(x["date"][[1]],"%Y%m%d"))
    if (length(unique(yrmn))>1) {
      x.tracked <- NULL
      if (progress) pb <- txtProgressBar(style=3)
      for (i in 1:length(unique(yrmn))) {
        if (progress) setTxtProgressBar(pb,i/(length(unique(yrmn))))
        if(verbose) print(unique(yrmn)[i])
        x.y <- subset(x,it=(yrmn==unique(yrmn)[i]))
        x.t <- Track.events(x.y,x.tracked,lplot=lplot,
                          progress=FALSE,verbose=verbose)
        if (is.null(x.tracked)) {
          x.tracked <- x.t
        } else {
          x.tracked <- merge(x.tracked,x.t,all=TRUE)
        }
      }
      x <- x.tracked
      x <- attrcp(x.t,x)
      class(x) <- class(x.t)
      invisible(x)
    } 
    dates <- x["date"][[1]]
    times <- x["time"][[1]]
    lons <- x["lon"][[1]]
    lats <- x["lat"][[1]]
    num <- rep(NA,dim(x)[1])
    dx <- rep(NA,dim(x)[1])
    datetime <- strptime(paste(dates,times),"%Y%m%d %H")
    d <- sort(unique(datetime))
    dh <- mean(as.numeric(diff(d,units="hours")))
    dmax <- dmax*dh
    nmax <- nmax/dh
    if (all(c("date","time","lon","lat","trajectory") %in% names(x0))) {
      dt0 <- strptime(paste(x0["date"][[1]],x0["time"][[1]]),"%Y%m%d %H")
      x0 <- x0[dt0==max(dt0[dt0<min(datetime)]),]
      dates <- c(x0["date"][[1]],dates)
      times <- c(x0["time"][[1]],times)
      lons <- c(x0["lon"][[1]],lons)
      lats <- c(x0["lat"][[1]],lats)
      num <- c(x0["trajectory"][[1]],num)
      dx <- c(rep(0,dim(x0)[1]),dx)
      datetime <- strptime(paste(dates,times),"%Y%m%d %H")
      i.start <- length(x0["lon"][[1]])+1
      d <- c(datetime[i.start-1],d)
    } else {
      num[datetime==d[1]] <- 1:sum(datetime==d[1])
      i.start <- 1
    }
    t1 <- Sys.time()
    if (progress) pb <- txtProgressBar(style=3)
    for (i in 1:(length(d)-1)) {
      if (progress) setTxtProgressBar(pb,i/(length(d)-1))
      nn <- NearestNeighbour(lons[datetime==d[i]],lats[datetime==d[i]],
                 lons[datetime==d[i+1]],lats[datetime==d[i+1]],dmax=dmax)
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
    while (any(tracklen$Freq>nmax)) {
      for (n in tracklen[tracklen$Freq>nmax,]$num) {
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
  n <- x["trackcount"][[1]]
  d <- x["tracklength"][[1]]
  options(digits=12)
  dnum <- as.numeric(paste(x["date"][[1]],x["time"][[1]],sep="."))
  ok <- rep(TRUE,length(n))
  ends <- dnum==min(dnum) | dnum==max(dnum)
  if (!is.null(nmin)) ok <- n>=nmin | ends
  if (!is.null(dmin)) ok <- ok & (d>=dmin | ends)
  browser()
  y <- subset(x,it=ok)
  rnum <- Enumerate(y)
  if (i.start>1) rnum <- rnum + max(num[1:(i.start-1)])
  y["trajectory"] <- rnum
  if(lplot) {
    lons <- y["lon"][[1]]
    lats <- y["lat"][[1]]
    num <- y["trajectory"][[1]]
    data(geoborders,envir=environment())
    plot(geoborders,type="l",col="grey20",lwd=0.5,xlim=c(-180,180),ylim=c(0,90))
    for (i in tracklen$num[tracklen$Freq>q995(tracklen$Freq)]) {
      points(lons[num==i],lats[num==i],col="blue",pch=1,cex=0.75)
      points(lons[num==i][1],lats[num==i][1],col="blue",pch=19,cex=0.75)
      lines(lons[num==i],lats[num==i],col="blue",lty=1,lwd=1)
    }
  }
  
  invisible(x)
}

NearestNeighbour <- function(lon1,lat1,lon2,lat2,dmax=1E6,lplot=FALSE,
                             verbose=FALSE) {
  if (verbose) print("NearestNeighbour")
  distance <- mapply(function(x,y) distAB(x,y,lon1,lat1),lon2,lat2)
  num <- as.numeric(apply(distance,2,function(x) which(rank(x)==1)))
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
  if (!identical(unique(y$trajectory),1:length(unique(y$trajectory)))) {
    rnum <- Enumerate(y)
  } else {
    rnum <- y$trajectory
  }
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
  rnum <- num[2:length(num)]-num[1:(length(num)-1)]
  rnum[rnum>0] <- 2:(sum(rnum>0)+1)
  rnum <- c(1,rnum)    
  while(any(rnum==0)) {
    rnum[rnum==0] <- rnum[which(rnum==0)-1]
  }
  invisible(rnum)
}

