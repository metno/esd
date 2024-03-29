#' 3-step cyclone tracking algorithm.
#' 
#' Applies a tracking algorithm to a set of cyclones (\code{\link{CCI}}).
#' 
#' The algorithm connects events in three subsequent time steps, chosing the
#' path that minimizes the total displacement as well as the change in angle
#' and displacement between them. The relative weight of these criteria can be
#' adjusted. The analysis can be applied to 'events' objects.
#' 
#' Note: The algorithm has been developed for tracking midlatitude cyclones in
#' the northern hemisphere and may not work as well for other regions or
#' 'events' of different types, e.g., anti-cyclones.
#' 
#' 
#' @aliases track track.default track.events
#' @param X An 'events' object containing temporal and spatial information
#' about a set of cyclones or anticyclones.
#' @param x0 A tracked 'events' object from previous time steps, used as a
#' starting point for the tracking of X so that trajectories can continue from
#' x0 to X.
#' @param it A list providing time index, e.g. month.
#' @param is A list providing space index, lon and/or lat.
#' @param dmax Maximum displacement of events between two time steps. Unit: m.
#' @param f.d Relative weight of the total displacement criterion in finding
#' the most probable trajectories.
#' @param f.dd Relative weight of the change in displacement as a criterion in
#' finding the most probable trajectories.
#' @param f.da Relative weight of the change in direction (angle) as a
#' criterion in finding the most probable trajectories.
#' @param nmax Maximum total lifetime of a trajectory. Unit: number of time
#' steps.
#' @param nmin Minimum total lifetime of a trajectory. Unit: number of time
#' steps.
#' @param dmin Minimum total length of a trajectory. Unit: m.
#' @param plot If TRUE, show plots of trajectories for selected time steps.
#' @param progress If TRUE, show progress bar.
#' @param verbose If TRUE, print out diagnosics.
#' @return An 'events' object containing the original information as well as
#' the trajectory number ('trajectory') of each event and statistical
#' properties of the trajectories ('trackcount' - number of events in path;
#' 'tracklen' - distance between start and end point of path').
#' @author K. Parding
#' @seealso CCI,as.trajectory
#' @keywords track
#' @examples
#' 
#' # Load sample data to use for example
#' # ERA5 6-hourly SLP data from the North Atlantic region, 2016-09-15 to 2016-10-15
#' data(slp.ERA5)
#' 
#' ## Cyclone identification
#' Cstorms <- CCI(slp.ERA5, m=20, label='ERA5', pmax=1000, verbose=TRUE, plot=FALSE)
#' 
#' ## Cyclone tracking
#' Ctracks <- track(Cstorms, plot=FALSE, verbose=TRUE)
#' 
#' ## Map with points and lines showing the cyclone centers and trajectories
#' map(Ctracks, type=c("trajectory","points"), col="blue", new=FALSE)
#' ## Map with only the trajectory and start and end points
#' map(Ctracks, type=c("trajectory","start","end"), col="red", new=FALSE)
#' ## Map showing the cyclone depth (slp) as a color scale (rd = red scale)
#' map(Ctracks, param="pcent", type=c('trajectory','start'), 
#'     colbar=list(pal="rd", rev=TRUE, breaks=seq(980,1010,5)), 
#'     alpha=0.9, new=FALSE)
#' 
#' ## Select only the long lasting trajectories...
#' Ct <- subset(Ctracks, ic=list(param='trackcount', pmin=12) )
#' map(Ct, new=FALSE)
#' ## ...or only the long distance ones...
#' Ct <- subset(Ctracks, ic=list(param='tracklength', pmin=3000) )
#' map(Ct, new=FALSE)
#' ## ...or only the deep cyclones
#' Ct <- subset(Ctracks, ic=list(param='pcent', pmax=980) )
#' map(Ct, new=FALSE)
#' 
#' ## Map of cyclone trajectories with the slp field in background
#' cb <- list(pal="budrd",breaks=seq(990,1040,5))
#' map(Ctracks, Y=slp.ERA5, it=as.POSIXct("2016-09-30 19:00"), colbar=cb, 
#'     verbose=TRUE, new=FALSE)
#' 
#' ## Transform the cyclones into a 'trajectory' object which takes up less space
#' Ctraj <- as.trajectory(Ctracks)
#' map(Ctraj, new=FALSE)
#' print(object.size(Ctracks), units="auto")
#' print(object.size(Ctraj), units="auto")
#' 
#' @export
track <- function(x,...) UseMethod("track")

#' @exportS3Method
#' @export track.events
track.events <- function(x,...,verbose=FALSE) {
  if(verbose) print("track.events")
  track.default(x,...,verbose=verbose)
}

#' @exportS3Method
#' @export track.default
track.default <- function(x,...,x0=NULL,it=NULL,is=NULL,dmax=1E6,nmax=200,nmin=3,dmin=1E5,
                          f.d=0.5,f.da=0.3,f.dd=0.2,f.dp=0,f.depth=0,dh=NULL,
		          fill.last=TRUE,greenwich=NULL,plot=FALSE,progress=TRUE,verbose=FALSE) {
  if(verbose) print("track.default")
  x <- subset(x,it=!is.na(x["date"][[1]]))
  x <- subset(x,it=it,is=is)
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  if (requireNamespace("PCICt", quietly = TRUE)) {
    d <- PCICt::as.PCICt(paste(x$date,x$time),format="%Y%m%d %H",cal=calendar)
  } else {
    d <- as.POSIXct(paste(x$date,x$time),format="%Y%m%d %H")
  }
  y0 <- NULL
  if(is.null(dh)) {
    dh <- min(as.numeric(diff(sort(unique(d)),units="hours")))
    ## Temporary fix for daylight saving time. Should try to find a more solid solution. 
    if(dh==23 & length(unique(x$time))==1) dh <- 24
    if(dh==11 & length(unique(x$time))==2) dh <- 12
    if(dh==5 & length(unique(x$time))==4) dh <- 6
    if(dh==2 & length(unique(x$time))==8) dh <- 3
  }
  if(!is.null(greenwich)) x <- g2dl(x,"greenwich")
  yrmn <- as.numeric(format(d,"%Y"))+as.numeric(format(d,"%m"))/12
  if (length(unique(yrmn))>1 & length(yrmn)>1000) {
    x.tracked <- NULL
    if (progress) pb <- txtProgressBar(style=3)
    for (i in 1:length(unique(yrmn))) {
      if(verbose) print(unique(yrmn)[i])
      x.y <- subset(x,it=(yrmn==unique(yrmn)[i]))
      if (is.null(x.tracked)) {
        x.t <- Track(x.y,x0=x0,plot=plot,cleanup.x0=FALSE,dh=dh,
                     dmax=dmax,dmin=dmin,nmax=nmax,nmin=nmin,
		     f.d=f.d,f.da=f.da,f.dd=f.dd,f.dp=f.dp,f.depth=f.depth,
                     progress=FALSE,verbose=verbose)
        if("trajectory" %in% names(x.t$y)) {
	  x.tracked <- x.t$y
	} else {
	  x0 <- merge(x0,x.t$y,all=TRUE)
	}
	if(!is.null(x0)) y0 <- subset(x.t$y0, it=(x.t$y0$date==max(x.t$y0$date)))
      } else {
        x.t <- Track(x.y,x0=x.tracked,plot=plot,dh=dh,
                     dmax=dmax,dmin=dmin,nmax=nmax,nmin=nmin,
		     f.d=f.d,f.da=f.da,f.dd=f.dd,f.dp=f.dp,f.depth=f.depth,
                     progress=FALSE,verbose=verbose)
        x.tracked <- x.t$y0
        x.tracked <- merge(x.tracked,x.t$y,all=TRUE)
      }
      if (progress) setTxtProgressBar(pb,i/(length(unique(yrmn))))
    }
    y <- x.tracked
  } else {
    x.tracked <- Track(x,x0=x0,plot=plot,dh=dh,#x0=NULL,
                       dmax=dmax,dmin=dmin,nmax=nmax,nmin=nmin,
		       f.d=f.d,f.da=f.da,f.dd=f.dd,f.dp=f.dp,f.depth=f.depth,
                       progress=progress,verbose=verbose)
    y <- x.tracked$y
    if(!is.null(x0)) y0 <- subset(x.tracked$y0, it=x.tracked$y0$date==max(x.tracked$y0$date))
  }
  y <- attrcp(x,y) 
  class(y) <- class(x)
  attr(y, "x0") <- y0
  if(any(is.na(y$trajectory)) & fill.last) {
    nok <- is.na(y$trajectory)
    y$trajectory[nok] <- seq(sum(nok)) + max(y$trajectory,na.rm=TRUE)
  }
  invisible(y)
}

# NOT EXPORTED
Track <- function(x,x0=NULL,it=NULL,is=NULL,dmax=1E6,nmax=124,nmin=3,dmin=1E5,
		  f.d=0.5,f.da=0.3,f.dd=0.2,f.dp=0,f.depth=0,dh=6,
		  cleanup.x0=TRUE,plot=FALSE,progress=TRUE,verbose=FALSE) {
  if (verbose) print("Track - cyclone tracking based on the distance and change in angle of direction between three subsequent time steps")
  options(digits=12)
  x <- x[!is.na(x[,1]), ]
  x <- x[order(x$date*1E2+x$time), ]
  x <- x[!duplicated(x), ]
  if(!is.null(x0)) {
    x0 <- x0[!is.na(x0[,1]), ]
    x0 <- x0[!duplicated(x0), ]
  }
  if(verbose) print("calculate trajectories")
  dates <- x$date
  times <- x$time
  lons <- x$lon
  lats <- x$lat
  pcent <- x$pcent
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  num <- rep(NA,dim(x)[1])
  dx <- rep(NA,dim(x)[1])
  if (requireNamespace("PCICt", quietly = TRUE)) {
    datetime <- PCICt::as.PCICt(paste(dates,times),format="%Y%m%d %H",cal=calendar)
  } else {
    datetime <- strptime(paste(dates,times),"%Y%m%d %H")  
  }
  d <- unique(sort(datetime))
  #dh <- min(as.numeric(diff(d,units="hours")))
  if(!is.null(x0)) {
    if (dim(x0)[1]>0) {
      if (requireNamespace("PCICt", quietly = TRUE)) {
        d0 <- sort(unique(PCICt::as.PCICt(paste(x0$date,x0$time),format="%Y%m%d %H",cal=calendar)))
      } else {
        d0 <- sort(unique(strptime(paste(x0$date,x0$time),"%Y%m%d %H")))
      }
      dh0 <- as.numeric((min(d)-max(d0))/(60*60))
    } else dh0 <- 1E5
  } else {
    dh0 <- 1E5
  }
  if (all(c("date","time","lon","lat","trajectory") %in% names(x0))) {
    num0 <- x0$trajectory
    n00 <- max(num0,na.rm=TRUE)
    dt0 <- x0$date*1E2 + x0$time
    nend0 <- unique(num0[dt0<max(dt0)])
    if (dh0<=dh*2) {
      dt0.max <- unique(sort(dt0))[max(1,length(unique(dt0)))]
      x00 <- x0[dt0>=dt0.max,]
      dates <- c(x00$date,dates)
      times <- c(x00$time,times)
      lons <- c(x00$lon,lons)
      lats <- c(x00$lat,lats)
      if(!is.null(pcent)) pcent <- c(x00$pcent,pcent)
      num <- c(x00$trajectory,num)
      if(!is.null(x00$dx)) dx <- c(x00$dx, dx) else dx <- c(rep(0,dim(x00)[1]),dx)
      nend0 <- nend0[!nend0 %in% num]
      i.start <- 1#length(x00$lon)-1
    } else {
      x00 <- NULL
      i.start <- 1
    }
    if (requireNamespace("PCICt", quietly = TRUE)) {
      datetime <- PCICt::as.PCICt(paste(dates,times),format="%Y%m%d %H",cal=calendar)
    } else {
      datetime <- strptime(paste(dates,times),"%Y%m%d %H")
    }
    d <- unique(c(datetime[i.start:length(datetime)],d))
    d <- d[!is.na(d)]
  } else {
    x00 <- NULL
    num[datetime==d[1]] <- 1:sum(datetime==d[1])
    i.start <- 1
    n00 <- 0
    nend0 <- 0
  }
  if(length(d)<3) {
    y <- x
    y0 <- x0
  } else {
  t1 <- Sys.time()
  if (progress) pb <- txtProgressBar(style=3)
  nend <- nend0
  n0 <- n00
  for (i in 2:(length(d)-1)) {
    if (progress) setTxtProgressBar(pb,i/(length(d)-1))
    dhi <- as.numeric((d[i:(i+1)]-d[(i-1):i])/(60*60))
    if(all(dhi<=dh*2)) {
      nn <- Track123(
        step1=list(lon=lons[datetime==d[i-1]],lat=lats[datetime==d[i-1]],
                   num=num[datetime==d[i-1]],dx=dx[datetime==d[i-1]],
                   pcent=pcent[datetime==d[i-1]]),
        step2=list(lon=lons[datetime==d[i]],lat=lats[datetime==d[i]],
                   num=num[datetime==d[i]],dx=dx[datetime==d[i]],
                   pcent=pcent[datetime==d[i]]),
        step3=list(lon=lons[datetime==d[i+1]],lat=lats[datetime==d[i+1]],
                   num=num[datetime==d[i+1]],dx=dx[datetime==d[i+1]],
                   pcent=pcent[datetime==d[i+1]]),
             f.d=f.d,f.da=f.da,f.dd=f.dd,f.dp=f.dp,f.depth=f.depth,
             dmax=dmax,n0=n0,nend=nend)
      num[datetime==d[i-1]] <- nn$step1$num
      num[datetime==d[i]] <- nn$step2$num
      num[datetime==d[i+1]] <- nn$step3$num
      dx[datetime==d[i]] <- nn$step2$dx
      dx[datetime==d[i+1]] <- nn$step3$dx
      n0 <- max(nn$n0,n0,na.rm=TRUE)
      nend <- nn$nend
    } else {
      if(any(is.na(num[datetime==d[i-1]]))) {
        i.solo <- datetime==d[i-1] & is.na(num)
	solo <- seq(sum(is.na(num[datetime==d[i-1]]))) +
            max(num[datetime %in% d[(i-1):(i+1)]],n0,na.rm=TRUE)
        num[i.solo] <- solo
        dx[i.solo] <- 0
        nend <- c(nend, solo)
      }
      if(any(is.na(num[datetime==d[i]]))) {
        i.solo <- datetime==d[i] & is.na(num)
	solo <- seq(sum(is.na(num[datetime==d[i]]))) +
            max(num[datetime %in% d[(i-1):(i+1)]],n0,na.rm=TRUE)
        num[i.solo] <- solo
        dx[i.solo] <- 0
        nend <- c(nend, solo)
      }
      if(any(is.na(num[datetime==d[i+1]]))) {
        i.solo <- datetime==d[i+1] & is.na(num)
	solo <- seq(sum(is.na(num[datetime==d[i+1]]))) +
            max(num[datetime %in% d[(i-1):(i+1)]],n0,na.rm=TRUE)
        num[i.solo] <- solo
        dx[i.solo] <- 0
      }
      n0 <- max(num[datetime %in% d[(i-1):(i+1)]],n0,na.rm=TRUE)
    }
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
  x$dx <- dx[i.start:length(num)]*1E-3
  if ( length(names(x))>length(attr(x,"unit")) ) {
    attr(x,"unit") <- c(attr(x,"unit"),c("numeration","km"))
  }
  if (verbose) print("calculate trajectory statistics")
  x <- trackstats(x,verbose=verbose)
  if(!is.null(x0)) {
    if(!is.null(x00)) {
      if (dim(x00)[1]>0) {
        x00$trajectory <- num[1:(i.start-1)]
        x00$dx <- dx[1:(i.start-1)]*1E-3
        x0[dt0>=max(dt0),] <- x00
      }
    }
    x01 <- rbind(x0,x)
    c01 <- attrcp(x,x01)
    class(x01) <- class(x)
    x01 <- trackstats(x01)
    dnum <- x01$date*1E2 + x01$time
    if (is.null(nmin)) nmin <- 3
    if (is.null(dmin)) dmin <- 0
    if (requireNamespace("PCICt", quietly = TRUE)) {
      starts <- dnum < as.numeric( format(PCICt::as.PCICt(as.character(min(dnum)),"%Y%m%d%H",cal=calendar) +
                                          (nmin-1)*dh*3600,"%Y%m%d%H") )
      ends <- dnum > as.numeric( format(PCICt::as.PCICt(as.character(max(dnum)),"%Y%m%d%H",cal=calendar) -
                                          (nmin-1)*dh*3600,"%Y%m%d%H") )
    } else {
      starts <- dnum < as.numeric( format(strptime(min(dnum),"%Y%m%d%H") +
                                          (nmin-1)*dh*3600,"%Y%m%d%H") )
      ends <- dnum > as.numeric( format(strptime(max(dnum),"%Y%m%d%H") -
                                          (nmin-1)*dh*3600,"%Y%m%d%H") )
    }
    ok <- x01$n>=nmin | ends | starts
    ok <- ok & (x01$distance>=(dmin*1E-3) | ends | starts)
    if(!cleanup.x0) ok[dnum<min(x$date*1E2 + x$time)]
    y01 <- subset(x01,it=ok)
    rnum <- y01$trajectory
    kvec <- unique(x01$trajectory[!ok])
    kvec <- kvec[!is.na(kvec)]
    for(k in kvec) {
      ik <- y01$trajectory>k & !is.na(rnum)
      if(any(ik)) rnum[ik] <- rnum[ik] - 1
    }
    dx <- Displacement(y01)
    y01$dx <- dx
    dnum <- y01$date*1E2 + y01$time
    dnum1 <- x$date*1E2 + x$time
    y0 <- subset(y01,it=dnum<min(dnum1))
    y <- subset(y01,it=dnum>=min(dnum1))
    attr(y0,"calendar") <- calendar
    attr(y,"calendar") <- calendar
  } else {
    dnum <- x["date"][[1]]*1E2 + x["time"][[1]]
    ok <- rep(TRUE,length(x$n))
    if (requireNamespace("PCICt", quietly = TRUE)) {
      ends <- dnum < as.numeric( format(PCICt::as.PCICt(as.character(min(dnum)),"%Y%m%d%H",cal=calendar) +
                                          (nmin-1)*dh*3600,"%Y%m%d%H") ) |
      dnum > as.numeric( format(PCICt::as.PCICt(as.character(max(dnum)),"%Y%m%d%H",cal=calendar) -
                                          (nmin-1)*dh*3600,"%Y%m%d%H") )
    } else {
      ends <- dnum < as.numeric( format(strptime(min(dnum),"%Y%m%d%H") +
                                          (nmin-1)*dh*3600,"%Y%m%d%H") ) |
      dnum > as.numeric( format(strptime(max(dnum),"%Y%m%d%H") -
                                          (nmin-1)*dh*3600,"%Y%m%d%H") )
    }
    ok <- x$n>=nmin | ends
    ok <- ok & (x$distance>=(dmin*1E-3) | ends)
    y <- subset(x,it=ok,verbose=verbose)
    rnum <- enumerate(y,verbose=verbose)
    rnum[y$trajectory<=n00 & !is.na(rnum)] <- y$trajectory[y$trajectory<=n00 & !is.na(rnum)]
    rnum[y$trajectory>n00 & !is.na(rnum)] <- n00 + rnum[y$trajectory>n00 & !is.na(rnum)]   
    y$trajectory <- rnum
    attr(y,"calendar") <- calendar
    y0 <- NULL
  }
  
  if(plot) {
    lons <- y$lon
    lats <- y$lat
    num <- y$trajectory
    dates <- y$date
    times <- y$time
    dplot <- unique(dates)[min(5,length(unique(dates)))]
    tplot <- unique(times[dates==dplot])[min(3,
                    length(unique(times[dates==dplot])))]
    nvec <- unique(num[dates==dplot & times==tplot])
    cols <- rainbow(length(nvec))
    if(max(lats)>0) ylim <- c(0,90) else ylim <- c(-90,0)
    if(attr(x,"greenwich")) xlim <- c(0,360) else xlim <- c(-180,180)
    geoborders <- NULL
    data("geoborders",envir=environment())
    if(!attr(x,"greenwich")) {
      plot(geoborders,type="l",col="grey20",lwd=0.5,
           xlim=xlim,ylim=ylim,main=paste(dplot,tplot))
    } else {
      lon.gb <- geoborders[,1] 
      lon.gb[lon.gb<0 & !is.na(lon.gb)] <-
          lon.gb[lon.gb<0 & !is.na(lon.gb)] + 360
      lon.w <- lon.gb
      lon.e <- lon.gb
      lon.w[lon.w>=180] <- NA
      lon.e[lon.e<180] <- NA
      plot(lon.w,geoborders[,2],type="l",col="grey20",lwd=0.5,
           xlim=xlim,ylim=ylim,main=paste(dplot,tplot))
      lines(lon.e,geoborders[,2],col="grey20",lwd=0.5) 
    }
    for (i in seq_along(nvec)) {
      points(lons[num==nvec[i]],lats[num==nvec[i]],col=cols[i],pch=1,cex=1.5)
      points(lons[num==nvec[i]][1],lats[num==nvec[i]][1],
             col=cols[i],pch=".",cex=3)
      points(lons[num==nvec[i] & dates==dplot & times==tplot][1],
           lats[num==nvec[i] & dates==dplot & times==tplot][1],
           col=adjustcolor(cols[i],alpha.f=0.5),pch=19,cex=1.5)
    }
  }
  }
  invisible(list(y=y,y0=y0))
}

# NOT EXPORTED
Track123 <- function(step1,step2,step3,n0=0,dmax=1E6,
                     f.d=0.5,f.da=0.3,f.dd=0.2,f.dp=0,f.depth=0,
		     nend=NA,plot=FALSE,verbose=FALSE) {
  if (verbose) print("Three step cyclone tracking")
  ## Set constants
  amax <- 90
  ## Normalize relative weights of the different criteria:
  if(!is.null(step3$pcent) & !is.null(step2$pcent) & !is.null(step1$pcent)) {
    f.dp <- 0
    f.depth <- 0
  }
  if( (f.d+f.da+f.dd+f.dp)!=1 ) {
    f.all <- f.d+f.da+f.dd+f.dp
    f.d <- f.d/f.all
    f.da <- f.da/f.all
    f.dd <- f.dd/f.all
    f.dp <- f.dp/f.all
  }
  
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
  ## Distance and direction between cyclones in timesteps 2 and 3
  d23 <- mapply(function(x,y) distAB(x,y,step2$lon,step2$lat),step3$lon,step3$lat)
  d23[is.na(d23)] <- 0
  a23 <- mapply(function(x,y) 180-angle(x,y,step2$lon,step2$lat),step3$lon,step3$lat)
  dim(d23) <- c(n2,n3)
  dim(a23) <- c(n2,n3)
  ## Distance and direction between cyclones in timesteps 2 and 3
  d12 <- mapply(function(x,y) distAB(x,y,step2$lon,step2$lat),step1$lon,step1$lat)
  d12[is.na(d12)] <- 0
  a12 <- mapply(function(x,y) angle(x,y,step2$lon,step2$lat),step1$lon,step1$lat)
  dim(d12) <- c(n2,n1)
  dim(a12) <- c(n2,n1)
  ## Calculate the change in direction for potential trajectorues
  da <- sapply(a12,function(x) abs(x-a23))
  da[da>180 & !is.na(da)] <- abs(da[da>180 & !is.na(da)]-360)
  dim(da) <- c(n2*n3,n1*n2)
  ## Displacement and change in displacement
  d <- sapply(d12,function(x) apply(d23,c(1,2),function(y) max(y,x)))
  ## KMP 2018-10-30: dd/d doesn't work if d=0, i.e., for cyclone that is stationary in steps 1,2,3
  d[d==0 & !is.na(d)] <- 0.1
  dd <- sapply(d12,function(x) apply(d23,c(1,2),function(y) abs(y-x)))
  ok.d <- sapply(d12<dmax,function(x) sapply(d23<dmax,function(y) y & x ))
  dim(dd) <- c(n2*n3,n1*n2)
  dim(ok.d) <- c(n2*n3,n1*n2)
  if(!is.null(step3$pcent) & !is.null(step2$pcent) & !is.null(step1$pcent)) {
    ## Mean sea level pressure at cyclone centers
    p23 <- sapply(step3$pcent, function(x) (step2$pcent+x)/2)
    p12 <- sapply(step1$pcent, function(x) (step2$pcent+x)/2)
    dim(p12) <- c(n2,n1)
    p <- sapply(p12,function(x) (x+p23)/2)
    dim(p) <- c(n2*n3,n1*n2)
    ## Change in sea level pressure at cyclone centers
    dp23 <- sapply(step3$pcent, function(x) step2$pcent-x)
    dp12 <- sapply(step1$pcent, function(x) step2$pcent-x)
    dim(dp12) <- c(n2,n1)
    dp <- sapply(dp12,function(x) (abs(x)+abs(dp23))/2)
    dim(dp) <- c(n2*n3,n1*n2)
  }
  
  ## Probability factors for three step trajectories...
  
  ## Indices for the probability factor (pf) matrix: Each point in pf represents the probability 
  ## of a possible trajectory. The index vectors j1 and j2 indicate which columns correspond to 
  ## which cyclones in time steps 1 and 2, and i2 and i3 indicate which which rows correspond
  ## to which cyclones in time steps 2 and 3.
  j1 <- as.vector(sapply(seq(n1),function(x) rep(x,n2)))
  j2 <- rep(seq(n2),n1)
  i2 <- rep(seq(n2),n3)
  i3 <- as.vector(sapply(seq(n3),function(x) rep(x,n2)))
  
  ## ...based on the pressure at the center of the cyclones
  if(!is.null(step3$pcent) & !is.null(step2$pcent) & !is.null(step1$pcent)) {
    p.range <-(max(p[!is.na(p)])-min(p[!is.na(p)]))
    pf.depth <- (max(p[!is.na(p)])-p)/p.range
    pf.dp <- dp/p.range#max(dp[!is.na(dp)])
    pf.depth[p.range==0] <- 0
    pf.dp[p.range==0] <- 0
  } else {
    pf.depth <- matrix(1,nrow(da),ncol(da))
    pf.dp <- matrix(1,nrow(da),ncol(da))
  }
  ## ...and the change in displacement, direction and pressure of cyclones
  pf <- 1 - f.d*d/dmax - f.da*da/amax - f.dd*dd/d - f.dp*pf.dp
  pf[pf<0] <- 0
  pf[!ok.d] <- 0
  pf <- pf + f.depth*pf.depth
  
  ## Set pf of impossible trajectories to NA...
  for(i in unique(i2)) pf[i==i2,i!=j2] <- NA
  ## ...and make sure not to replace already existing trajectories.
  if(any(!is.na(step2$num))) {
    for(sn in step2$num[!is.na(step2$num)]) {
      if(any(sn %in% step1$num)) {
        i <- which(step2$num==sn & !is.na(step2$num))
        j <- which(step1$num!=sn | is.na(step1$num))
        pf[i2 %in% i, j1 %in% j] <- NA
        j <- which(step1$num==sn & !is.na(step1$num))
        i <- which(step2$num!=sn | is.na(step2$num))
        pf[i2 %in% i, j1 %in% j] <- NA
      }
    }
  }

  ## Set pf of trajectories that have ended...
  if(any(step1$num %in% nend)) {
    i <- which(step1$num %in% nend & !is.na(step1$num))
    pf[, j1 %in% i] <- NA
  }

  # Probability factors for broken trajectories - disabled
  # f.break <- 0.2
  # if(!is.null(step3$pcent) & !is.null(step2$pcent) & !is.null(step1$pcent)) {
  #   pf.depth12 <- (max(p[!is.na(p)])-p12)/(max(p[!is.na(p)])-min(p[!is.na(p)]))
  #   pf.dp12 <- dp12/(max(p[!is.na(p)])-min(p[!is.na(p)]))
  #   pf.depth23 <- (max(p[!is.na(p)])-p23)/(max(p[!is.na(p)])-min(p[!is.na(p)]))
  #   pf.dp23 <- dp23/(max(p[!is.na(p)])-min(p[!is.na(p)]))
  # } else {
  #   pf.dp12 <- matrix(1,nrow(d12),ncol(d12))
  #   pf.dp23 <- matrix(1,nrow(d23),ncol(d23))
  #   pf.depth12 <- matrix(1,nrow(d12),ncol(d12))
  #   pf.depth23 <- matrix(1,nrow(d23),ncol(d23))
  # } 
  # pf.12 <- 1 - f.d*d12/dmax - f.dp*pf.dp12 - (f.da+f.dd) - f.break
  # pf.23 <- 1 - f.d*d23/dmax - f.dp*pf.dp23 - (f.da+f.dd) - f.break
  # pf.12[pf.12<0] <- 0
  # pf.23[pf.23<0] <- 0
  # pf.12[d12>=dmax] <- 0
  # pf.23[d23>=dmax] <- 0
  # pf.12 <- pf.12 + f.depth*pf.depth12
  # pf.23 <- pf.23 + f.depth*pf.depth23

  ## Put the probability factors into a common matrix.
  ## The indices in j1.all, j2.all, i2.all, i3.all keeps track
  ## of which row (i) and column (j) represents which cyclones in 
  ## timesteps 1, 2, and 3.
  # pf.all <- matrix(0,ncol=n1*n2+n2,nrow=n2*n3+n2)
  # j1.all <- c(j1,rep(NA,n2))
  # j2.all <- rep(seq(n2),n1+1)
  # i2.all <- rep(seq(n2),n3+1)
  # i3.all <- c(i3,rep(NA,n2))
  pf.all <- pf
  j1.all <- j1
  j2.all <- j2
  i2.all <- i2
  i3.all <- i3
  
  ## Put probability factors for 3-step trajectories into matrix
  pf.all[1:nrow(pf),1:ncol(pf)] <- pf #pf.d*pf.change
  ## Put probability factors for broken trajectories into matrix - disabled
  # for(k in unique(j2)) {
  #   j.k <- which(is.na(j1.all) & j2.all==k)
  #   i.k <- which(!is.na(i3.all) & i2.all==k)
  #   pf.all[i.k,j.k] <- pf.23[k,]
  #   j.k <- which(!is.na(j1.all) & j2.all==k)
  #   i.k <- which(is.na(i3.all) & i2.all==k)
  #   pf.all[i.k,j.k] <- pf.12[k,]
  # }
  
  ## Connect the most likely trajectories based on the probability factors
  if(any(pf.all>0 & !is.na(pf.all))) {
    rank.all <- matrix(rank(1-pf.all),dim(pf.all))
    if(plot) {
      dev.new()
      geoborders <- NULL
      data("geoborders",envir=environment())
      plot(geoborders,type="l",col="grey20",lwd=0.5,
           xlim=c(-90,90),ylim=c(30,90))      
      points(step1$lon,step1$lat,col="hotpink",pch=21,lwd=3)
      points(step2$lon,step2$lat,col="red",pch=21,lwd=3)
      points(step3$lon,step3$lat,col="orange",pch=21,lwd=3)
      for(k in sort(unique(rank.all[pf.all>0]))) {
        ij <- which(rank.all==k,arr.ind=TRUE)
        i1.k <- j1.all[ij[2]]
        i2.k <- i2.all[ij[1]]
        i3.k <- i3.all[ij[1]]
        lon.k <- c(step1$lon[i1.k],step2$lon[i2.k],step3$lon[i3.k])
        lat.k <- c(step1$lat[i1.k],step2$lat[i2.k],step3$lat[i3.k])
        lines(lon.k,lat.k,lwd=1,lty=1,
         col=adjustcolor("black",alpha.f=pf.all[rank.all==k]))
      }
    }
    
    #print(rank.all[i2.all==8&i3.all==10,j1.all==11&j2.all==8])
    rank.all[pf.all<=0 | is.na(pf.all)] <- NA
    while(any(!is.na(rank.all))) {
      ij <- which(rank.all==min(rank.all,na.rm=TRUE),arr.ind=TRUE)
      # If more than one trajectory with same ranking:
      if(dim(ij)[1]>1) {
	# If any three step trajectories, choose among those based on angle change
        is.123 <- ij[,1]<=dim(da)[1] & ij[,2]<=dim(da)[2]
        if(sum(is.123)>0 & sum(!is.123)>0) { 
	        ij <- ij[which(is.123),]
	        if(is.null(dim(ij))) dim(ij) <- c(1,length(ij))
	        is.123 <- ij[,1]<=dim(da)[1] & ij[,2]<=dim(da)[2]
	      }
        if(sum(is.123)>1) { 
          ij <- ij[which.min(apply(ij,1,function(x) da[x[1],x[2]])),]
	        if(is.null(dim(ij))) dim(ij) <- c(1,length(ij))
	      }
	      # If still more than one trajectory, choose based on distance
	      if(dim(ij)[1]>1) {
          d.k <- apply(ij,1,function(x) {
	          if(!is.na(j1.all[x[2]])) {
	            return(d12[i2.all[x[1]],j1.all[x[2]]])
	          } else {
	            return(d23[i2.all[x[1]],i3.all[x[1]]])
	          }
          })
	        k <- which.min(d.k)
	        ij <- ij[k,]
	        if(is.null(dim(ij))) dim(ij) <- c(1,length(ij))
	      }
	      # If still more than one trajectory, choose first one
	      if(!is.null(dim(ij))) ij <- ij[1,]
      }
      k1 <- j1.all[ij[2]]
      k2 <- i2.all[ij[1]]
      k3 <- i3.all[ij[1]]
      if(!is.na(k1)) {
        num.k <- step1$num[k1]
        rank.all[,j1.all==k1 & !is.na(j1.all)] <- NA
      } else {
        num.k <- max(c(n0,step1$num,step2$num,step3$num),na.rm=TRUE)+1
      }
      step2$num[k2] <- num.k
      step2$dx[k2] <- d12[k2,k1]
      rank.all[,j2.all==k2] <- NA
      rank.all[i2.all==k2,] <- NA
      if(!is.na(k3)) {
        step3$num[k3] <- num.k
        step3$dx[k3] <- d23[k2,k3]
        rank.all[i3.all==k3,] <- NA
      }
      if(plot) {
        lon.k <- c(step1$lon[k1],step2$lon[k2],step3$lon[k3])
        lat.k <- c(step1$lat[k1],step2$lat[k2],step3$lat[k3])
        lines(lon.k,lat.k,lwd=3,col="blue",type="l",lty=1)
      }
      n0 <- max(c(n0,step1$num,step2$num,step3$num),na.rm=TRUE)
    }
  }
  
  if(all(is.na(step2$num))) {
    step2$num[is.na(step2$num)] <- seq(sum(is.na(step2$num))) + n0
  } else if (any(is.na(step2$num))) {
    step2$num[is.na(step2$num)] <- seq(sum(is.na(step2$num))) +
        max(c(n0,step2$num),na.rm=TRUE)    
  }
  ## trajectories that should end: tracks that are in step 1 but not in step 3
  nok1 <- (!(step1$num %in% step3$num) & !is.na(step1$num))
  if(any(nok1)) nend <- c(nend,step1$num[nok1])
  nend <- unique(nend[!is.na(nend)])
  
  return(list(step1=step1,step2=step2,step3=step3,nend=nend,
         n0=max(c(step1$num,step2$num,step3$num,n0),na.rm=TRUE)))
}

# NOT EXPORTED
angle <- function(lon1,lat1,lon2,lat2) {
  a <- 360 - (atan2(lat2-lat1,lon2-lon1)*(180/pi) + 360) %% 360
  #a[a>180] <- 360-a[a>180]
  return(a)
  #return(atan2(lat2-lat1,lon2-lon1)*180/pi+90)
}

#' Calculate trajectory statistics
#'
#' The function enumerates the trajectories ("trajectory"),
#' adds time steps and the total number of steps in a trajecotry ("trackcount"),
#' length of trajectory (in km): "tracklength"),
#' and if the distance "dx" beween time steps exists the
#' length of the trajectory from start to end ("tracklength") is also calculated.
#'
#' @param x an \code{events} object
#' @param verbose a boolean; if TRUE print information about progress
#' @return an \code{events} object with statistics describing the trajectories
#'
#' @export
trackstats <- function(x,verbose=FALSE) {
  if(verbose) print("trackstats")
  if (!any("trajectory" %in% names(x))) x <- track(x,verbose=verbose)
  y <- x[order(x$trajectory),]
  y <- attrcp(x,y)
  class(y) <- class(x)
  lons <- y["lon"][[1]]
  lats <- y["lat"][[1]]

  if(verbose) print("Enumerate")
  rnum <- enumerate(y)
  nummax <- max(rnum,na.rm=TRUE)
  rnum[is.na(rnum)] <- nummax+1
  if(verbose) print("trackcount")
  trackcount <- data.frame(table(rnum))

  if(verbose) print("timestep")
  timestep <- rep(1, length(rnum))
  for(k in 1:max(trackcount$Freq)) {
    ik <- duplicated(rnum) & c(1, diff(timestep))==0
    timestep[ik] <- timestep[ik] + 1
  }
  #ts <- unlist(sapply(unique(rnum),function(i) 1:trackcount$Freq[trackcount$rnum==i]))
  #timestep <- rep(NA,length(ts))
  #timestep[order(rnum)] <- ts
  if (any(rnum>nummax)) timestep[rnum==(nummax+1)] <- 1

  if(verbose) print("distance")
  fn <- function(x) {
    if(dim(x)[1]==1) {
      dx <- 0 
    } else if (x[1,1]==x[dim(x)[1],1] & x[1,2]==x[dim(x)[1],2]) {
      dx <- 0
    } else {
      dx <- distAB(x[1,1],x[1,2],x[dim(x)[1],1],x[dim(x)[1],2])
    } 
    return(dx)
  }
  distance <- as.numeric(by(cbind(lons,lats),rnum,fn))*1E-3

  y$trackcount <- trackcount$Freq[rnum]
  y$trackcount[is.na(y$trajectory)] <- 1
  y$timestep <- timestep
  y$distance <- distance[rnum]
  if (length(attr(y,"unit"))<dim(y)[2]) {
    attr(y,"unit") <- c(attr(y,"unit"),c("count","numeration","km"))
  }
  if('dx' %in% colnames(y)) {
    if(verbose) print("dx")
    dx <- Displacement(y)
    if (any(rnum>nummax)) dx[nummax+1] <- 0
    y$dx <- dx
    if(verbose) print("tracklength")
    tracklength <- as.numeric(by(y$dx,rnum,function(x) sum(x,na.rm=TRUE)))
    if (any(rnum>nummax)) tracklength[nummax+1] <- 0
    y$tracklength <- tracklength[rnum]
    if(length(attr(y,"unit"))<dim(y)[2]) {
      attr(y,"unit") <- c(attr(y,"unit"),c("km","km"))
    }
  }
  invisible(y)
}

# NOT EXPORTED
Displacement <- function(x,verbose=FALSE) {
  if(verbose) print("Calculate displacement")
  if(!"trajectory" %in% colnames(x)) {
    if(verbose) print("track.events")
    x <- track.events(x)
  }
  if("dx" %in% colnames(x)) {
    dx <- x$dx
  } else {
    dx <- rep(0,dim(x)[1])
  }
  kvec <- unique(x$trajectory[x$n>1 & x$timestep>1 & dx==0])
  kvec <- kvec[!is.na(kvec)]
  for(k in kvec) {
    ik <- x$trajectory==k & !is.na(x$trajectory)
    if(sum(ik)>1) {
      dk <- mapply(distAB,x$lon[ik][1:(sum(ik)-1)],x$lat[ik][1:(sum(ik)-1)],
                          x$lon[ik][2:sum(ik)],x$lat[ik][2:sum(ik)])
      dk[is.na(dk)] <- 0
      if(length(dk)<sum(ik)) dk <- c(0,dk)
      dkk <- try(dk*1E-3)
      if (!inherits(dkk,"try-error")) {
        dx[ik] <- dk*1E-3
      } else {
        print("oops, something went wrong!")
      }
    } else if (sum(ik)==1) dx[ik] <- 0
  }
  invisible(dx)
}


#' Enumerate trajectories
#'
#' The function enumerates the trajectories ("trajectory"). It is used in the function \code{"trackstats"} and can be applied to re-enumerate cyclone trajectories in an \code{events} object after selecting a subset using the function \code{"subset"}. 
#'
#' @param x an \code{events} object
#' @param param parameter name of the trajectory numbers. Default: trajectory.
#' @param verbose a boolean; if TRUE print information about progress
#' @return an \code{events} object with statistics describing the trajectories
#'
#' @export
enumerate <- function(x,param="trajectory",verbose=FALSE) {
  if(verbose) print("enumerate")
  stopifnot(inherits(x,"data.frame"))
  num <- x[param][[1]]
  if (identical(unique(num),seq_along(unique(num)))) {
    rnum <- num
  } else if (length(num)==1) {
    rnum <- num
  } else {
    if(verbose) print(paste("number of events:",length(num)))
    if(verbose) print(paste("number of trajectories:",length(unique(num))))
    rnum <- diff(num)
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

# NOT EXPORTED
NearestNeighbour <- function(lon1,lat1,lon2,lat2,dmax=1E6,
                             plot=FALSE,verbose=FALSE) {
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
  if (plot) {
    plot(lon1,lat1,type="p",col="black",pch=1:length(lon1),lwd=1,
      xlim=c(-180,180),ylim=c(0,90))
   points(lon2,lat2,col="blue",pch=num,lwd=1)
  }
  invisible(cbind(num,d.num))
}

