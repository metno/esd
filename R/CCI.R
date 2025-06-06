#' Calculus Cyclone identification.
#' 
#' Identifies cyclones (low pressure systems) in a gridded data set using a
#' Calculus Cyclone Identification (CCI) method (EMS04-A-00146, Benestad, R.E.;
#' Sorteberg, A.; Chen, D. 'Storm statistics derived using a calculus-based
#' cyclone identification method',
#' \url{http://www.cosis.net/members/meetings/sessions/oral_programme.php?p_id=110&s_id=1845},
#' European Meteorological Society AC2, Nice, Sept 28, 2004). Storms are
#' identified with longitude, latitude, and date. Also returned are estimates
#' of local minimum pressure, max pressure gradient near storm, max geostrophic
#' and gradient windspeed near storm, and radius of the storm. The storm
#' location is by means of finding where first derivatives of north--south and
#' east--west gradients both are zero. The storm radius is estimated from the
#' points of inflexion along the latitude and longitude lines running trough
#' the centre of the storm.
#' 
#' If a north-south or east-west sea level pressure (p) profile can be
#' approximated as \deqn{p(\theta) = p_0 + \sum_{i=1}^{N_\theta} [ a_\theta(i)
#' \cos(\omega_\theta(i) \theta) + b_\theta(i) \sin(\omega_\theta(i) \theta)
#' ]}{p = p0 + sum [ a cos(w t) + b sin(w t) ]}
#' 
#' Then the pressure gradient can be estimated as: \deqn{\frac{\partial
#' \hat{p}(\theta)}{\partial \theta} = \sum_{i=1}^{N_\theta} \omega_\theta(i)[
#' -\hat{a}_\theta(i) \sin(\omega_\theta(i) \theta) + \hat{b}_\theta(i)
#' \cos(\omega_\theta(i) \theta)]}{dp/dt = sum (w[ -a sin(w t) + b cos(w w)])}
#' 
#' The gradient in x and y directions are found after the transform
#' \deqn{\frac{d\hat{p}(x)}{dx} = \frac{1}{a \cos(\phi)}
#' \frac{d\hat{p}(\theta)}{d\theta}}{dp/dx= 1/[a cos(t)] dp(t)/dt} and
#' \deqn{\frac{d\hat{p}(y)}{dy} = \frac{1}{a}
#' \frac{d\hat{p}(\phi)}{d\phi}}{dp/dy = 1/a dp/dt} (Gill, 1982).
#' 
#' This code is based on the CCI method of the R-package 'cyclones' and has
#' been adapted for 'esd'.
#' 
#' The maximum gradient wind (max.speed) is estimated as described in Fleagle
#' and Businger (1980) p. 163. (eq 4.27).
#' 
#' Reference: Benestad \& Chen (2006) 'The use of a Calculus-based Cyclone
#' Identification method for generating storm statistics', Tellus A, in press.
#' Benestad (2005)
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar data
#' @importFrom graphics points image contour lines
#'
#' @aliases CCI
#'
#' @param Z A field object.
#' @param m Number of harmonics used for fit to profile (Fourier truncation),
#' which decides the degree of smoothing of the field. Lower values of m result
#' in greater smoothing.
#' @param it A list providing time index, e.g. month.
#' @param is A list providing space index, lon and/or lat.
#' @param label Label for ID-purposes.
#' @param cyclones TRUE: identifies cyclones, FALSE: anticyclones.
#' @param mindistance Min distance between cyclones. Unit: m.
#' @param dpmin Min pressure gradient at points of inflection around cyclone.
#' Unit: Pa/m (10hPa/km).
#' @param pmax Max pressure in center of cyclone. Unit: hPa.
#' @param rmin Min average radius of cyclone. Unit: m.
#' @param rmax Max average radius of cyclone. Unit: m.
#' @param nsim Number of simultaneous cyclones identified and saved ordered
#' according to depth/strength (NULL = no limit).
#' @param do.track TRUE: tracks the cyclones with the 'track' function, FALSE:
#' no tracking is applied.
#' @param fname Name of output file.
#' @param plot TRUE: show intermediate results as plots.
#' @param greenwich a boolean; if TRUE longitudes are transformed to a system starting at the Greenwich line (0-360E);
#'        if FALSE longitudes are transformed to a system starting at the date line (180W-180E)
#' @param progress a boolean; if TRUE show progress bar
#' @param allow.open a boolean; if TRUE allow (anti)cyclones that are not closed, i.e.,
#' that have a point of inflexion on only 3 of 4 sides.   
#' @param verbose a boolean; if TRUE print out diagnostics.
#' @param accuracy Not yet finished.
#' @param \dots additional arguments
#'
#' @return The CCI function returns an 'events' object (a data frame with
#' spatio-temporal information) that is organized as follows:
#' \code{as.events(X=data.frame(date,time,lon,lat,pcent,max.dspl,
#' max.speed,radius,qf),label=label,...}, where \code{date} and \code{time} are
#' vectors containing the date and time of each cyclone, \code{lon} and
#' \code{lat} are the longitude and latitude of the cyclone centers (unit:
#' degrees), \code{pcent} is the pressure at the center of each cyclone (unit:
#' hPa), \code{max.dpsl} is the maximum pressure gradient associated with each
#' cyclone (unit: hPa/m), \code{max.speed} is the estimated maximum windspeed
#' (unit: m/s), \code{radius} is the cyclone radius (unit: km), and \code{qf}
#' is a kind of quality flag (1 = OK, 2 = less spatially precise, identified
#' after widening the pressure gradient zero-crossings).
#' @author K.M. Parding & R.E. Benestad
#' @seealso track.events
#' @keywords manip,cyclones,CCI
#' @examples
#' 
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
#' map(Ct)
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
#' @export CCI
CCI <- function(Z,m=12,it=NULL,is=NULL,cyclones=TRUE,greenwich=NULL,
                label=NULL,mindistance=5E5,dpmin=1E-3,
                rmin=1E4,rmax=2E6,nsim=NULL,progress=TRUE,
                fname="cyclones.rda",plot=FALSE,accuracy=NULL,
                allow.open=FALSE,do.track=FALSE,
		anomaly=TRUE,pmax=NULL,verbose=FALSE,...) {
  if(verbose) print("CCI - calculus based cyclone identification")

  stopifnot(inherits(Z,'field'))
  Z <- subset(Z,it=it,is=is,verbose=verbose)
  if(is.null(greenwich) & !is.null(attr(Z,"greenwich"))) {
    greenwich <- attr(Z,"greenwich")
  } else if (is.null(greenwich) & is.null(attr(Z,"greenwich"))) {
    greenwich <- !(any(lon(Z)<0) | !any(lon(Z)>180))
  }  
  Z <- g2dl(Z,greenwich=greenwich)

  if(diff(range(lon(Z)))>=357) {
    ## KMP 2025-02-11: If the field covers a whole spherical cap, split into two overlapping regions
    ##   and do CCI separately, then combine. Otherwise cyclones near the "edge" (0/360 or -180/180) will be missed.
    if(verbose) print("Splitting field into two overlapping regions")
    if(verbose) print("Apply CCI to longitudes [-120, 120]")
    XA <- CCI(subset(Z, is=list(lon=c(-120, 120), lat=range(lat(Z)))), m=m, it=it, is=is, 
                     cyclones=cyclones, mindistance=mindistance, dpmin=dpmin,
                     rmin=rmin, rmax=rmax, nsim=nsim, progress=progress, fname=NULL,
                     accuracy=accuracy, allow.open=allow.open, do.track=FALSE,
                     plot=plot, anomaly=anomaly, pmax=pmax, verbose=verbose, ...)
    if(verbose) print("Apply CCI to longitudes [60, 300]")
    XB <- CCI(subset(Z, is=list(lon=c(60, 300), lat=range(lat(Z)))), m=m, it=it, is=is, 
                     cyclones=cyclones, mindistance=mindistance, dpmin=dpmin,
                     rmin=rmin, rmax=rmax, nsim=nsim, progress=progress, fname=NULL,
                     accuracy=accuracy, allow.open=allow.open, do.track=FALSE,
                     plot=plot, anomaly=anomaly, pmax=pmax, verbose=verbose, ...)
    X <- combine(XA, XB, verbose=verbose)
    if(do.track) X <- track(X,verbose=verbose,...)
    if(!is.null(fname)) save(file=fname, X)
    invisible(X)
  } else {
  if(is.null(pmax)) if(anomaly) pmax <- 0 else pmax <- 1012
  yrmn <- format(index(Z),"%Y")#"%Y%m")
  if (length(unique(yrmn))>2) {
    t1 <- Sys.time()  
    if (progress) pb <- txtProgressBar(style=3)
    X <- NULL
    if (progress) setTxtProgressBar(pb,0/(length(unique(yrmn))))
    for (i in 1:length(unique(yrmn))) {
      if(verbose) print(unique(yrmn)[i])
      Z.y <- subset(Z,it=(yrmn==unique(yrmn)[i]))
      if(dim(Z.y)[1]>0) {
        X.y <- CCI(Z.y,m=m,cyclones=cyclones,
                label=label,mindistance=mindistance,dpmin=dpmin,
                pmax=pmax,rmin=rmin,rmax=rmax,nsim=nsim,progress=FALSE,
                fname=NULL,plot=plot,accuracy=accuracy,verbose=verbose)
        if(progress) setTxtProgressBar(pb,i/(length(unique(yrmn))))
        if(is.null(X)) {
          X <- X.y
        } else {
          X <- merge(X,X.y,all=TRUE)
          X <- attrcp(X.y,X)
          class(X) <- class(X.y)
        }
      }
      if(!is.null(fname)) save(file=fname,X)
    }
    t2 <- Sys.time()
    if (verbose) print(paste('CCI took',round(as.numeric(t2-t1,units="secs")),'s'))
    if(!is.null(fname)) save(file=fname,X)
    invisible(X)
  } else {

  ## Rearrange time index
  t <- as.numeric(format(index(Z),format="%Y%m%d%H%M"))

  ## Anomaly of the pressure field
  if(anomaly) Z <- as.anomaly(Z)

  ## Calculate first and second derivative
  if(verbose) print("Calculate first and second derivative")
  #browser()
  resx <- dX(Z,m=m,accuracy=accuracy,verbose=verbose,progress=progress)
  resy <- dY(Z,m=m,accuracy=accuracy,verbose=verbose,progress=progress)
  
  ## Spatial resolution
  dx <- abs(diff(resx$lon))[1]
  dy <- abs(diff(resx$lat))[1]
  
  ## Search for zero crossings of the first derivative:
  ## Reorganize and reshape
  if(verbose) print("Reshape arrays")
  nx <- length(resx$lon); ny <- length(resx$lat); nt <- dim(Z)[1]
  lonXY <- rep(0.5*(resx$lon[2:nx]+resx$lon[1:(nx-1)]),ny-1)
  dim(lonXY) <- c(nx-1,ny-1)
  latXY <- rep(0.5*(resx$lat[2:ny]+resx$lat[1:(ny-1)]),nx-1)
  dim(latXY) <- c(ny-1,nx-1); latXY <- t(latXY)
  Zx <- as.matrix(resx$Z.fit); dim(Zx) <- c(nt,nx,ny)
  Zy <- as.matrix(resy$Z.fit); dim(Zy) <- c(nt,nx,ny)
  dslpdy <- as.matrix(resy$dZ); dim(dslpdy) <- c(nt,nx,ny)
  dslpdx <- as.matrix(resx$dZ); dim(dslpdx) <- c(nt,nx,ny)
  dslpdx2 <- as.matrix(resx$dZ2); dim(dslpdx2) <- c(nt,nx,ny)
  dslpdy2 <- as.matrix(resy$dZ2); dim(dslpdy2) <- c(nt,nx,ny)
  px <- 0.25*(Zx[,1:(nx-1),2:ny] + Zx[,2:nx,2:ny] +
              Zx[,1:(nx-1),1:(ny-1)] + Zx[,2:nx,1:(ny-1)])
  py <- 0.25*(Zy[,1:(nx-1),2:ny] + Zy[,2:nx,2:ny] +
              Zy[,1:(nx-1),1:(ny-1)] + Zy[,2:nx,1:(ny-1)])
  dim(px) <- c(nt,nx-1,ny-1)
  dim(py) <- c(nt,nx-1,ny-1)
      
  ## Clear temporary objects from working memory
  rm("Zx","Zy","resx","resy"); gc(reset=TRUE)

  ## Check units
  if (attr(Z,"unit")=='Pa') {
    if (verbose) print("transform pressure units: Pa -> hPa")
    px <- px*1E-2
    py <- py*1E-2
  }
   
  ## Search for zero crossing in y-direction
  if(verbose) print("Find zero crossing of first derivative in y-direction")
  P.lowy <- rep(0,nt*(nx-1)*(ny-1)); dim(P.lowy) <- c(nt,nx-1,ny-1) 
  dy11 <- 0.5*(dslpdy[,2:nx,2:ny]+dslpdy[,1:(nx-1),2:ny])
  dy12 <- 0.5*(dslpdy[,2:nx,1:(ny-1)]+dslpdy[,1:(nx-1),1:(ny-1)])
  dy21 <- 0.5*(dslpdy2[,2:nx,2:ny]+dslpdy2[,1:(nx-1),2:ny])
  dy22 <- 0.5*(dslpdy2[,2:nx,1:(ny-1)]+dslpdy2[,1:(nx-1),1:(ny-1)])
  if (cyclones) { i.low <- (dy11*dy12 < 0) & (dy21+dy22 > 0) &
              is.finite(dy11+dy12+dy21+dy22)
  } else { i.low <- (dy11*dy12 < 0) & (dy21+dy22 < 0) &
              is.finite(dy11+dy12+dy21+dy22)
  }
  P.lowy[i.low] <- 1
  DY <- 0.5*(dy11+dy12)
  DY2 <- 0.5*(dy21+dy22)

  ## Zero crossing in x-direction
  if(verbose) print("Find zero crossing of first derivative in x-direction")
  P.lowx <- rep(0,nt*(nx-1)*(ny-1)); dim(P.lowx) <- c(nt,nx-1,ny-1) 
  dx11 <- 0.5*(dslpdx[,2:nx,2:ny] + dslpdx[,2:nx,1:(ny-1)])
  dx12 <- 0.5*(dslpdx[,1:(nx-1),2:ny] + dslpdx[,1:(nx-1),1:(ny-1)])
  dx21 <- 0.5*(dslpdx2[,2:nx,2:ny] + dslpdx2[,2:nx,1:(ny-1)])
  dx22 <- 0.5*(dslpdx2[,1:(nx-1),2:ny] + dslpdx2[,1:(nx-1),1:(ny-1)])
  if (cyclones) { i.low <- (dx11*dx12 < 0) & (dx21+dx22 > 0) &
              is.finite(dx11 + dx12 + dx21 + dx22)
  } else { i.low <- (dx11*dx12 < 0) & (dx21+dx22 < 0) &
              is.finite(dx11 + dx12 + dx21 + dx22)
  }
  P.lowx[i.low] <- 1
  DX <- 0.5*(dx11+dx12)
  DX2 <- 0.5*(dx21+dx22)
  
  ## Clear temporary objects from working memory
  rm("dx11","dx12","dx21","dx22","dy11","dy12","dy21","dy22",
     "dslpdx","dslpdx2","dslpdy","dslpdy2","i.low"); gc(reset=TRUE)
  
  ## Find zero crossings in both directions
  ## Plowx & P.lowy are matrices with 0's and 1's.
  lows1 <- (P.lowy & P.lowx)
  
  ## Find cyclones that were missed in the first search
  ## using widened masks of zero crossings
  if(verbose) print("Find extra cyclones using widened masks")
  ## Mask the cylcones already selected
  P.lowx2 <- P.lowx; P.lowy2 <- P.lowy
  ## Widen mask in x-direction
  P.lowx2[,2:(nx-1),] <- P.lowx2[,1:(nx-2),] | P.lowx2[,2:(nx-1),]
  P.lowx2[,1:(nx-2),] <- P.lowx2[,1:(nx-2),] | P.lowx2[,2:(nx-1),]
  ## Widen mask in y-direction
  P.lowy2[,,2:(ny-1)] <- P.lowy2[,,1:(ny-2)] | P.lowy2[,,2:(ny-1)]
  P.lowy2[,,1:(ny-2)] <- P.lowy2[,,1:(ny-2)] | P.lowy2[,,2:(ny-1)]
  ## Find new zero crossings
  P.lowx2[lows1] <- 0; P.lowy2[lows1] <- 0
  lows2 <- (P.lowy2 & P.lowx2)
  ## Clear temporary objects from working memory
  rm("P.lowx2","P.lowy2"); gc(reset=TRUE)
  
  ## Cyclones should be deeper than the zonal average slp
  if(verbose) print("Compare slp to zonal average")
  dP <- 0.5*(px+py)
  P.zonal <- apply(dP,c(1,3),FUN="mean",na.rm=TRUE)
  for(i in seq(dim(dP)[2])) dP[,i,] <- dP[,i,] - P.zonal
  if(cyclones) {
    lows1[dP>0] <- FALSE
    lows2[dP>0] <- FALSE
  } else {
    lows1[dP<0] <- FALSE
    lows2[dP<0] <- FALSE
  }
  rm("dP","P.zonal"); gc(reset=TRUE)

  # Exclude identified depressions in high altitude regions
  if(verbose) print("Penalty factor for high altitude")
  etopo5 <- NULL
  data('etopo5', envir = environment())
  fn <- function(lon=0,lat=60) {
    i.lon <- which.min(abs(attr(etopo5,"longitude")-lon))
    i.lat <- which.min(abs(attr(etopo5,"latitude")-lat))
    h <- etopo5[i.lon,i.lat]
    h[h<0] <- 0
    nlon <- max(lon(etopo5))
    nlat <- max(lat(etopo5))
    dhdx <- abs(etopo5[min(i.lon+1,nlon),i.lat]-etopo5[max(i.lon-1,1),i.lat])/
            (distAB(lon+dx,lat,lon-dx,lat)*1E-3)
    dhdy <- abs(etopo5[i.lon,min(i.lat+1,nlat)]-etopo5[i.lon,max(i.lat-1,1)])/
            (distAB(lon,lat+dx,lon,lat-dx)*1E-3)
    return(c(mean(h),mean(dhdx),mean(dhdy)))
  }

  ## Penalty factor for topography
  hmax <- 1000
  h <- mapply(fn,lonXY,latXY)
  dh <- apply(h[2:3,],2,max)  
  pf.h <- 1 - 0.25*h[1,]/hmax - 0.25*dh/200
  pf.h[h[1,]<0] <- 1
  pf.h[h[1,]>hmax] <- 0
  pf.h <- rep(pf.h,nt)
  dim(pf.h) <- c(nx-1,ny-1,nt)
  pf.h <- aperm(pf.h,c(3,1,2))
 
  ## Quality flag to keep track of cyclones found with widened mask
  qf <- matrix(rep(0,length(lows1)),dim(lows1))
  qf[lows1] <- 1
  qf[lows2] <- 2

  ## Make sure dimensions are correct
  dim(DX) <- c(nt, nx-1, ny-1)
  dim(DX2) <- dim(DX)
  dim(DY) <- dim(DX)
  dim(DY2) <- dim(DX)
  
  ## Lat, lon, and dates of cyclones
  lon<-rep(lonXY,nt); dim(lon)<-c(nx-1,ny-1,nt); lon<-aperm(lon,c(3,1,2)) 
  lat<-rep(latXY,nt); dim(lat)<-c(nx-1,ny-1,nt); lat<-aperm(lat,c(3,1,2))
  date<-rep(t,(nx-1)*(ny-1)); dim(date)<-c(nt,nx-1,ny-1)
  ## Cyclones found with ordinary method
  lon1 <- lon[lows1]; lat1 <- lat[lows1]; date1 <- date[lows1]
  pcent1 <- 0.5*(px[lows1]+py[lows1])
  strength1 <- rank(pcent1)
  if (!cyclones) strength1 <- rank(-pcent1)
  ## Cyclones found after widening zero crossing masks
  lon2 <- lon[lows2]; lat2<-lat[lows2]; date2 <- date[lows2]
  pcent2 <- 0.5*(px[lows2]+py[lows2])
  strength2 <- rank(pcent2)
  if (!cyclones) strength2 <- rank(-pcent2)

  ## Keep only cyclones that are deeper than pmax
  if(verbose) print("Remove shallow cyclones")
  del1 <- rep(TRUE,length(date1))
  del2 <- rep(TRUE,length(date2))
  if(!is.null(pmax)) {
    if(is.function(pmax)) {
      pmax <- pmax(Z,na.rm=T)
    } else if(is.character(pmax)) {
      pmax <- eval(parse(text=paste(pmax,"(Z,na.rm=T)",sep="")))
    }
    ok1 <- pcent1 < pmax
    if(!cyclones) ok1 <- pcent1 > pmax
    del1[!ok1] <- FALSE
    ok2 <- pcent2 < pmax
    if(!cyclones) ok2 <- pcent2 > pmax
    del2[!ok2] <- FALSE
    rm("ok1","ok2"); gc(reset=TRUE)
  }
  lows1[lows1] <- del1
  lows2[lows2] <- del2
  lon1 <- lon[lows1]; lat1 <- lat[lows1]; date1 <- date[lows1]
  pcent1 <- 0.5*(px[lows1] + py[lows1])
  strength1 <- rank(pcent1)
  if (!cyclones) strength1 <- rank(-pcent1)
  lon2 <- lon[lows2]; lat2 <- lat[lows2]; date2 <- date[lows2]
  pcent2 <- 0.5*(px[lows2] + py[lows2])
  strength2 <- rank(pcent2)
  if (!cyclones) strength2 <- rank(-pcent2)
 
  ## Remove secondary cyclones near a deeper one (same cyclonic system):
  if(verbose) print("Remove secondary cyclones")
  del1 <- rep(TRUE,length(date1))
  del2 <- rep(TRUE,length(date2))
  for (d in t) {
    
    i1 <- which(date1==d)
    i2 <- which(date2==d)

    ## Remove secondary cyclones identified with the widened masks,
    ## requiring 1000 km distance to nearest neighbouring cyclone
    if(any(i2) & any(i1)) {
      distance <- apply(cbind(lon2[i2],lat2[i2]),1,
       function(x) suppressWarnings(distAB(x[1],x[2],lon1[i1],lat1[i1])))
      if(length(i1)>1) {
        del.i2 <- unique(which(distance<mindistance,arr.ind=TRUE)[,2])
      } else {
        del.i2 <- unique(which(distance<mindistance,arr.ind=TRUE))
      }
      del2[i2[del.i2]] <- FALSE
      i2 <- i2[!1:length(i2) %in% del.i2]
    }
    if(length(i2)>1) {
      distance <- apply(cbind(lon2[i2],lat2[i2]),1,
       function(x) suppressWarnings(distAB(x[1],x[2],lon2[i2],lat2[i2])))
      diag(distance) <- NA; distance[lower.tri(distance)] <- NA
      del.i2 <- which(distance<mindistance,arr.ind=TRUE)
      if(any(del.i2)) {
        col.del <- rep(1,length(del.i2)/2)
        if (is.null(dim(del.i2))) {
          s1 <- strength2[i2][del.i2[1]]
          s2 <- strength2[i2][del.i2[2]]
          col.del[s1<=s2] <- 2
          del.i2 <- unique(del.i2[col.del])
        } else {
          s1 <- strength2[i2][del.i2[,1]]
          s2 <- strength2[i2][del.i2[,2]]
          col.del[s1<=s2] <- 2
          del.i2 <- unique(del.i2[cbind(seq(1,dim(del.i2)[1]),col.del)])
        }
        del2[i2[del.i2]] <- FALSE
        i2 <- i2[!1:length(i2) %in% del.i2]
      }
    }

    ## Remove secondary cyclones identified with the standard method,
    ## located close (<1000km) to other stronger cyclones
    if (length(i1)>1) {
      distance <- apply(cbind(lon1[i1],lat1[i1]),1,
       function(x) suppressWarnings(distAB(x[1],x[2],lon1[i1],lat1[i1])))
      diag(distance) <- NA; distance[lower.tri(distance)] <- NA
      del.i1 <- which(distance<mindistance,arr.ind=TRUE)
      if(any(del.i1)) {
        col.del <- rep(1,length(del.i1)/2)
        if (is.null(dim(del.i1))) {
           s1 <- strength1[i1][,del.i1[1]]
           s2 <- strength1[i1][,del.i1[2]]
           col.del[s1<=s2] <- 2
           del.i1 <- unique(del.i1[col.del])
        } else {
           s1 <- strength1[i1][del.i1[,1]]
           s2 <- strength1[i1][del.i1[,2]]
           col.del[s1<=s2] <- 2
           del.i1 <- unique(del.i1[cbind(seq(1,dim(del.i1)[1]),col.del)])
        }
        del1[i1[del.i1]] <- FALSE
      }
      i1 <- i1[!1:length(i1) %in% del.i1]
    }
  }
  lows1[lows1] <- del1
  lows2[lows2] <- del2
 
  ## Add the primary and secondary cyclones together,
  ## keep track of which is the two groups in the quality flag qf
  lows <- lows1 | lows2
  if(sum(lows)==0) {
    print("No cyclones identified!")
    X <- data.frame(date=NA,time=NA,lon=NA,lat=NA,pcent=NA,
         max.gradient=NA,max.speed=NA,radius=NA,closed=NA,accuracy=NA)
  } else {
    ## Clear temporary objects from working memory
    lon <- lon[lows]
    lat <- lat[lows]
    date <- date[lows]
    pcent <- 0.5*(px[lows]+py[lows])
    qf <- qf[lows]
    pf.h <- pf.h[lows]
   
    rm("lows","lows1","lows2","lon1","lon2","lat1","lat2",
     "date1","date2","strength1","strength2",
     "pcent1","pcent2","del1","del2"); gc(reset=TRUE)
  
    ## Put cyclones in order of date
    i <- order(date)
    lon <- lon[i]
    lat <- lat[i]
    date <- date[i]
    pcent <- pcent[i]
    qf <- qf[i]
    pf.h <- pf.h[i]
    strength <- rank(pcent)
    if (!cyclones) strength <- rank(-pcent)

    ## Keep only the nsim strongest cyclones each time step
    if(!is.null(nsim)) {
      s <- sapply(date,function(d) rank(strength[date==d]))
      lon <- lon[s<nsim]
      lat <- lat[s<nsim]
      date <- date[s<nsim]
      pcent <- pcent[s<nsim]
      qf <- qf[s<nsim]
      pf.h <- pf.h[s<nsim]
    }

    ## Pressure gradient
    if(verbose) print("Pressure gradient")
    rho <- 1.2922
    dpsl <- sqrt(DX^2+DY^2)
    if (attr(Z,"unit") %in% c("hPa","mbar")) dpsl <- dpsl*100

    ## Remove temporary variables and release the memory:
    rm('DX','DY'); gc(reset=TRUE)

    if(allow.open) nmin <- 3 else nmin <- 4
        
    # Find points of inflexion (2nd derivative==0) to estimate
    # the storm radius and maximum speed and pressure gradient
    t1 <- Sys.time()
    if(verbose) print("Find points of inflexion")
    NX <- dim(lonXY)[1]; NY <- dim(latXY)[2]
    lonXX <- rep(0.5*(lonXY[2:NX,2:NY]+lonXY[1:(NX-1),2:NY]),nt)
    dim(lonXX) <- c(NX-1,NY-1,nt); lonXX <- aperm(lonXX,c(3,1,2))
    latXX <- rep(0.5*(latXY[2:NX,2:NY]+latXY[2:NX,1:(NY-1)]),nt)
    dim(latXX) <- c(NX-1,NY-1,nt); latXX <- aperm(latXX,c(3,1,2))
    radius <- rep(NA,length(date))
    max.gradient <- rep(NA,length(date))
    max.speed <- rep(NA,length(date))
    max.vg <- rep(NA,length(date))
    closed <- rep(0,length(date))
    ok <- rep(TRUE,length(date))
    for (i in seq(1,length(date))) {
      if(verbose) print(paste(i,date[i]))
      inflx <- DX2[date[i]==t,2:NX,latXY[1,]==lat[i]]*
        DX2[date[i]==t,1:(NX-1),latXY[1,]==lat[i]]
      infly <- DY2[date[i]==t,lonXY[,1]==lon[i],2:NY]*
        DY2[date[i]==t,lonXY[,1]==lon[i],1:(NY-1)]
      if(length(inflx)>nrow(lonXY)) browser()
      lon.infl <- lonXY[inflx<0,1]
      lat.infl <- latXY[1,infly<0]
      dlon <- lon.infl-lon[i]
      dlat <- lat.infl-lat[i]
      ## KMP 2024-08-23: The edges have only positive or negative dlon/dlat values
      ##  and the two lines below results in an annoying warning
      #dlat_min <- c(min(dlat[dlat>0]),max(dlat[dlat<0]))])
      #dlon_min <- c(min(dlon[dlon>0]),max(dlon[dlon<0]))])
      if(any(dlat>0) & any(dlat<0)) dlat_min <- c(min(dlat[dlat>0]), max(dlat[dlat<0])) else 
        if(any(dlat>0)) dlat_min <- min(dlat[dlat>0]) else dlat_min <- max(dlat[dlat<0])
      if(any(dlon>0) & any(dlon<0)) dlon_min <- c(min(dlon[dlon>0]), max(dlon[dlon<0])) else 
        if(any(dlon>0)) dlon_min <- min(dlon[dlon>0]) else dlon_min <- max(dlon[dlon<0])
      ilat <- which(latXY[1,] %in% lat.infl[dlat %in% dlat_min])
      ilon <- which(lonXY[,1] %in% lon.infl[dlon %in% dlon_min])
      nlon <- length(ilon)
      nlat <- length(ilat)
      ilon <- c(rep(which(lonXY[,1]==lon[i]),nlat),ilon)
      ilat <- c(ilat,rep(which(latXY[1,]==lat[i]),nlon))
      oki <- sum(!is.na(ilon))>=nmin
      if(oki) {
        dpi <- mapply(function(i1,i2) dpsl[t==date[i],i1,i2],ilon,ilat)
        oki <- sum(!is.na(dpi) & (pf.h[i]*dpi)>dpmin/2)>=nmin &
               (pf.h[i]*mean(dpi,na.rm=TRUE))>dpmin
      }
      if (oki) {
        ri <- distAB(lon[i],lat[i],lonXY[ilon,1],latXY[1,ilat])
        fi <- 2*7.29212*1E-5*sin(pi*abs(latXY[1,ilat])/180)
        vg <- dpi/(fi*rho)
        v.grad <- -0.5*fi*ri*(1 - sqrt(1 + 4*vg/(fi*ri)))
          
        radius[i] <- mean(ri,na.rm=TRUE)
        max.gradient[i] <- mean(dpi,na.rm=TRUE)
        max.speed[i] <- mean(v.grad,na.rm=TRUE)
        max.vg[i] <- mean(vg,na.rm=TRUE)
        closed[i] <- floor(length(dpi>dpmin & !is.na(dpi))/4)
      } else {
        ok[i] <- FALSE
      }
    }
    t2 <- Sys.time()
    if (verbose) print(paste('finding points of inflexion took',
                             round(as.numeric(t2-t1,units="secs")),'s'))
    if (verbose) print("transform pressure gradient units: Pa/m -> hPa/km")
    max.gradient <- max.gradient*1E-2*1E3
    if(verbose) print("remove cyclones according to rmin, rmax, dpmin")
    
    ok <- ok
    if(!is.null(rmin)) ok <- ok & radius>=rmin
    if(!is.null(rmax)) ok <- ok & radius<=rmax
    lon <- lon[ok]
    lat <- lat[ok]
    date <- date[ok]
    pcent <- pcent[ok]
    qf <- qf[ok]
    closed <- closed[ok]
    #dslp <- dslp[ok]
    max.gradient <- max.gradient[ok]
    max.speed <- max.speed[ok]
    radius <- radius[ok]
    strength <- rank(pcent)
    if (!cyclones) strength <- rank(-pcent)

    if (plot) {
      if(verbose) print("plot example of cyclone identification")
      geoborders <- NULL
      data('geoborders',envir=environment())
      i <- length(date)/2
      inflx <- DX2[date[i]==t,2:NX,latXY[1,]==lat[i]]*
        DX2[date[i]==t,1:(NX-1),latXY[1,]==lat[i]]
      infly <- DY2[date[i]==t,lonXY[,1]==lon[i],2:NY]*
        DY2[date[i]==t,lonXY[,1]==lon[i],1:(NY-1)]
      pxi <- px[date[i]==t,,];  pyi <- py[date[i]==t,,]
      xi <- lonXY[,1]; yi <- latXY[1,]; zi <- pxi
      if (!(all(diff(xi)>0))) {xi <- rev(xi); zi <- apply(zi,2,rev)}
      if (!(all(diff(yi)>0))) {yi <- rev(yi); zi <- t(apply(t(zi),2,rev))}
      dev.new()
      plot(lonXY[,1],pxi[,latXY[1,]==lat[i]],lty=1,type="l",main=date[i],
         xlab="lon",ylab="slp (hPa)")
      points(lon[i],pxi[lonXY[,1]==lon[i],latXY[1,]==lat[i]],col="blue",pch=19)
      points(lonXY[inflx<0,1],pxi[inflx<0,latXY[1,]==lat[i]],col="red",pch=1)
      #dev.copy2eps(file="cyclones.lon.eps", paper="letter")#; dev.off()
      dev.new()
      plot(latXY[1,],pyi[lonXY[,1]==lon[i],],lty=1,type="l",main=date[i],
         xlab="lat",ylab="slp (hPa)")
      points(lat[i],pyi[lonXY[,1]==lon[i],latXY[1,]==lat[i]],col="blue",pch=19)
      points(latXY[1,infly<0],pyi[lonXY[,1]==lon[i],infly<0],col="red",pch=1)
      #dev.copy2eps(file="cyclones.lat.eps", paper="letter")#; dev.off()
      dev.new()
      image(xi,yi,zi,main=date[i],col=colscal(pal="budrd",n=14,rev=FALSE),
          xlab="lon",ylab="lat",breaks=seq(940,1080,10))
      contour(xi,yi,zi,add=TRUE,col='Grey40',lty=1,zlim=c(940,1010),nlevels=6)
      contour(xi,yi,zi,add=TRUE,col='Grey40',lty=2,zlim=c(1020,1080),nlevels=5)
      if(!greenwich) {
        lines(geoborders,col="grey40")
      } else {
        gb.w <- geoborders
        gb.e <- geoborders
        lon.gb <- gb.w[,1]
        lon.gb[!is.na(lon.gb) & lon.gb<0] <-
          lon.gb[!is.na(lon.gb) & lon.gb<0] + 360
        lon.w <- lon.gb
        lon.e <- lon.gb
        lon.w[!is.na(lon.w) & lon.w>=180] <- NA
        lon.e[!is.na(lon.e) & lon.e<180] <- NA
        gb.w[,1] <- lon.w
        gb.e[,1] <- lon.e
        lines(gb.w,col="grey40")
        lines(gb.e,col="grey40")
      }
      a <- which(P.lowx[t==date[i],,]==1,arr.ind=TRUE)
      b <- which(P.lowy[t==date[i],,]==1,arr.ind=TRUE)
      lon.a <- mapply(function(i1,i2) lonXY[i1,i2],a[,1],a[,2])
      lat.a <- mapply(function(i1,i2) latXY[i1,i2],a[,1],a[,2])
      lon.b <- mapply(function(i1,i2) lonXY[i1,i2],b[,1],b[,2])
      lat.b <- mapply(function(i1,i2) latXY[i1,i2],b[,1],b[,2])
      points(lon.a,lat.a,col="black",pch=1,cex=0.5,lwd=0.5)
      points(lon.b,lat.b,col="black",pch="|",cex=0.5,lwd=0.5)
      j <- date==date[i]
      col <- rep("black",sum(j,na.rm=TRUE))
      sz <- rep(2,sum(j,na.rm=TRUE))
      col[closed[j]==0] <- "grey50"
      sz[qf[j]==2] <- 1
      points(lon[j],lat[j],pch=21,lwd=2,bg="white",col=col,cex=sz)
      points(lon[i],lat[i],pch=4,lwd=2,col="black",cex=1)
      #dev.copy2eps(file="cyclones.map.eps", paper="letter")#; dev.off()
    }

    ## Remove temporary variables and release the memory:
    rm('lonXY','latXY','inflx','infly','DX2','DY2','px','py'); gc(reset=TRUE)
  
    ## Arrange results
    dd <- round(date*1E-4)
    hh <- round((date-dd*1E4)*1E-2)
    X <- data.frame(date=dd,time=hh,lon=lon,lat=lat,pcent=pcent,
         #dslp=dslp,
         max.gradient=max.gradient,max.speed=max.speed,
         radius=radius*1E-3,closed=closed,accuracy=qf)
  }
  unit <- c("date","hour CET","degrees","degrees","hPa","hPa",
            "hPa/km","m/s","km","TRUE/FALSE","grid points")
  if (cyclones) {
    longname <- "low pressure systems identified with CCI method"
    param <- "cyclones"
  } else {
    longname <- "high-pressure systems identified with CCI method"
    param <- "anti-cyclones"
  }
  X <- as.events(X,unit=unit,longname=longname,greenwich=greenwich,
         param=param,src=attr(Z,"source"),file=attr(Z,"file"),
         calendar=attr(Z,"calendar"),
         method="calculus based cylone identification, CCI",
         version="CCI in esd v1.0 (after October 6, 2015)",
         reference="Benestad & Chen, 2006, The use of a calculus-based cyclone identification method for generating storm statistics, Tellus A 58(4), 473-486.",
         url="http://onlinelibrary.wiley.com/doi/10.1111/j.1600-0870.2006.00191.x/abstract")
  if(do.track) X <- track(X,verbose=verbose,...)
  if(!is.null(fname)) save(file=fname,X)
  invisible(X)
  }
  }
}

#library(esd)
### MONTHLY SLP DATA:
##slp <- slp.ERAINT()
##cyclones <- CCI(slp,is=list(lon=c(-180,180),lat=c(0,90)),it='djf')
### 6-HOURLY SLP DATA:
## ncfile <- "/home/kajsamp/data/ecmwf/ERAinterim_slp_1979.nc"
## ncid <- open.ncdf(ncfile)
## lon <- get.var.ncdf(ncid,"longitude")
## lat <- get.var.ncdf(ncid,"latitude")
## time <- get.var.ncdf(ncid,"time")
## tunits <- att.get.ncdf(ncid, "time", "units")
## slp <- get.var.ncdf(ncid,"msl")
## close.ncdf(ncid)
## time <- as.POSIXlt(time*60*60,origin="1900-01-01")
## nt <- length(time); nlon <- length(lon); nlat <- length(lat)
## slp <- aperm(slp,c(3,1,2)); dim(slp) <- c(nt,nlon*nlat)
## slp <- zoo(slp,order.by=time)
## slp <- as.field(slp,lon=lon,lat=lat,
##               unit="Pa",longname="mean sea level pressure",
##               param="slp",quality=NULL,src="ERAinterim")
## #slp <- subset(slp,it="december")
## #Z <- slp
## cyclones <- CCI(slp,it="december",fname="cyclones.ERAint.1979.12.rda")

