# K Parding, 15.10.2015

Track.events <- function(x,it=NULL,is=NULL,dmax=1E6,lplot=FALSE,verbose=FALSE) {
  if (verbose) print("Track.events - nearest neighbour cyclone tracking")
  stopifnot(inherits(x,"events"))
  x <- subset(x,it=!is.na(x["date"][[1]]))
  dates <- x["date"][[1]]
  times <- x["time"][[1]]
  lons <- x["lon"][[1]]
  lats <- x["lat"][[1]]
  num <- rep(NA,dim(x)[1])
  datetime <- strptime(paste(dates,times),"%Y%m%d %H")
  d <- sort(unique(datetime))
  num[datetime==d[1]] <- 1:sum(datetime==d[1])
  t1 <- Sys.time()
  pb <- txtProgressBar(style=3)
  for (i in 1:(length(d)-1)) {
    setTxtProgressBar(pb,i/(length(d)-1))
    nn <- NearestNeighbour(lons[datetime==d[i]],lats[datetime==d[i]],
                 lons[datetime==d[i+1]],lats[datetime==d[i+1]],dmax=dmax)
    num.i <- num[datetime==d[i]][nn]
    num.i[is.na(num.i)] <- 1:sum(is.na(num.i)) + max(num,na.rm=TRUE)
    num[datetime==d[i+1]] <- num.i
  }
  t2 <- Sys.time()
  if (verbose) print(paste('Nearest neighbour tracking took',
                round(as.numeric(t2-t1,units="secs")),'s'))
  if(lplot) {
    data(geoborders)
    plot(geoborders,type="l",col="grey20",lwd=0.5,xlim=c(-180,180),ylim=c(0,90))
    for (i in 1:20) {
      points(lons[num==i],lats[num==i],type="both",col="blue",pch=num[i],lwd=1)
    }
  }
  x$tracknumber <- num
  invisible(x)
}

NearestNeighbour <- function(lon1,lat1,lon2,lat2,dmax=1E6,lplot=FALSE,verbose=FALSE) {
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
  if (lplot) {
    plot(lon1,lat1,type="p",col="black",pch=1:length(lon1),lwd=1,
       xlim=c(-180,180),ylim=c(0,90))
    points(lon2,lat2,col="blue",pch=num,lwd=1)
  }
  invisible(num)
}
