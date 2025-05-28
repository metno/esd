#' Calculate density of trajectories
#'
#' Internal function used in trajectory2density and events2field
#'
#' @param lons longitudes of trajectory
#' @param lats latitudes of trajectory
#' @param track trajectory number
#' @param dx spatial resolution of output field in east-west direction (unit: degrees east)
#' @param dy spatial resolution of output field in north-south direction (unit: degrees north)
#' @param radius radius within which to look for trajectories for each grid point (unit: m)
#' @param type "track" or "trajectory": calculate density of trajectories;
#' "genesis", "cyclogenesis" or "start": calculate density of cyclogenesis  events;
#' "lysis", "cyclolysis" or "end": calculate density of cyclolysis events
#' @param verbose if TRUE print progress
#'
#' @export
trackdensity <- function(lons,lats,track=NULL,dx=NULL,dy=NULL,
                         radius=5E5,type="track",verbose=FALSE) {
  if (is.null(dx)) dx <- min(diff(sort(unique(lons))))
  if (is.null(dy)) dy <- min(diff(sort(unique(lats))))
  fn <- function(A) {
    A <- unique(A)
    if(dim(A)[1]>1) {
      if(type %in% c("track","trajectory")) {
        B <- approx.lonlat(A$lon,A$lat,n=50)
        lon <- B[,1]
        lat <- B[,2]
      } else if (type %in% c("genesis","cyclogenesis","start")) {
        lon <- A[1,1]
        lat <- A[1,2]
      } else if (type %in% c("lysis","cyclolysis","end")) {
        lon <- A[nrow(A),1]
        lat <- A[nrow(A),2]          
      } else stop(paste("invalid input type =",type))
    } else { 
      lon <- A$lon
      lat <- A$lat
    } 
    xvec <- seq(round(min(lon)/dx)*dx-round(5+20*max(abs(lat))/90)*dx,
                round(max(lon)/dx)*dx+round(5+20*max(abs(lat))/90)*dx,dx)
    yvec <- seq(round(min(lat)/dy)*dy-ceiling(5/dy)*dy,
                round(max(lat)/dy)*dy+ceiling(5/dy)*dy,dy)
    xx <- as.vector(sapply(xvec,function(x) rep(x,length(yvec))))
    yy <- rep(yvec,length(xvec))
    if(any(xx > 360) | any(xx < -180)) {
      xx[xx > 360] <- xx[xx > 360] - 360
      xx[xx < -180] <- xx[xx < -180] + 360
    }
    if(length(lon)>1) {
      i <- lapply(1:length(lon),function(i) distAB(lon[i],lat[i],xx,yy)<radius)
      rx <- unlist(lapply(i,function(j) xx[j]))
      ry <- unlist(lapply(i,function(j) yy[j]))
    } else {
      i <- distAB(lon,lat,xx,yy)<radius
      rx <- xx[i]
      ry <- yy[i]
    }
    xy <- unique(data.frame(x=rx,y=ry))    
    invisible(xy)
  }
  if (is.null(track)) track <- rep(1,length(lons))
  lonlat <- do.call(rbind,by(data.frame(lon=lons,lat=lats),track,fn))
  hits <- as.data.frame(table(lon=lonlat[,1],lat=lonlat[,2]))
  hx <- factor2numeric(hits$lon)
  hy <- factor2numeric(hits$lat)
  d <- hits$Freq
  invisible(data.frame(lon=hx,lat=hy,density=d))
}
