read.best.track <- function(url='http://www.usno.navy.mil/NOOC/nmfc-ph/RSS/jtwc/best_tracks/',domain='io',
                            start=1945,end=2014,n=20,verbose=TRUE) {
  ## domain = c('io','sh','wp')
  ## Format  	
  ## BASIN,CY, YYYYMMDDHH,TECHNUM, TECH,TAU, LatN/S,LonE/W, VMAX,MSLP
  ## TY,RAD, WINDCODE,RAD1, RAD2,RAD3, RAD4,RADP, RRP,MRD, GUSTS,EYE,
  ## SUBREGION,MAXSEAS, INITIALS,DIR, SPEED,STORMNAME, DEPTH,SEAS, SEASCODE,SEAS1, SEAS2,SEAS3, SEAS4 
  
  TC <- list(info='best-track cyclones')
  storms <- c()
  ii <- 1
  for (yr in start:end) {
    for (i in 1:30) {
      if (i < 10) ci <- paste(0,i,sep='') else ci <- as.character(i)
      urldata <- paste(url,yr,'/',yr,'s-b',domain,'/b',domain,ci,yr,'.txt',sep='')
      print(urldata)


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
#' @aliases Cylcone tracking algorithm track.default track.events Track
#' Trackstats
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
#' @param plot TRUE: Show plots of trajectories for selected time steps.
#' @param progress TRUE: Show progress bar.
#' @param verbose TRUE: Print out diagnosics.
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
#' Cstorms <- CCI(slp.ERA5, m=20, label='ERA5', pmax=1000, verbose=TRUE, plot=TRUE)
#' 
#' ## Cyclone tracking
#' Ctracks <- track(Cstorms, plot=TRUE, verbose=TRUE)
#' 
#' ## Map with points and lines showing the cyclone centers and trajectories
#' map(Ctracks, type=c("trajectory","points"), col="blue")
#' ## Map with only the trajectory and start and end points
#' map(Ctracks, type=c("trajectory","start","end"), col="red")
#' ## Map showing the cyclone depth (slp) as a color scale (rd = red scale)
#' map(Ctracks, param="pcent", type=c('trajectory','start'), 
#'     colbar=list(pal="rd", rev=TRUE, breaks=seq(980,1010,5)), alpha=0.9)
#' 
#' ## Select only the long lasting trajectories...
#' Ct <- subset(Ctracks, ic=list(param='trackcount', pmin=12) )
#' map(Ct)
#' ## ...or only the long distance ones...
#' Ct <- subset(Ctracks, ic=list(param='tracklength', pmin=3000) )
#' map(Ct)
#' ## ...or only the deep cyclones
#' Ct <- subset(Ctracks, ic=list(param='pcent', pmax=980) )
#' map(Ct)
#' 
#' ## Map of cyclone trajectories with the slp field in background
#' cb <- list(pal="budrd",breaks=seq(990,1040,5))
#' map(Ctracks, slp.ERA5, it=as.POSIXct("2016-09-30 19:00"), colbar=cb, verbose=TRUE)
#' 
#' ## Transform the cyclones into a 'trajectory' object which takes up less space
#' Ctraj <- as.trajectory(Ctracks)
#' map(Ctraj)
#' print(object.size(Ctracks), units="auto")
#' print(object.size(Ctraj), units="auto")
#' 
#' @export track
      track <- try(read.table(urldata,sep=','),silent=TRUE)
      if (inherits(track,'try-error')) {
        urldata <- paste(url,yr,'/',yr,'s-b',domain,'/b',domain,ci,yr,'.dat',sep='')
        track <- try(read.table(urldata,sep=','),silent=TRUE)
      }
      if (!inherits(track,'try-error')) {
        tc <- paste(domain,yr,ii,sep='.')
        TC[[tc]] <- track
        storms <- c(storms,tc)
      }
    }
  }
  
  ## Create an esd trajectory-object from the tropical cyclone storm tracks
  ns <- length(storms)
  Y <- matrix(rep(NA,(3*n+3)*ns),ns,3*n+3)
  for (i in 1:ns) {
    z <- TC[[storms[i]]]
    south <- regexpr('S',z$V7) > 0
    west <- regexpr('W',z$V8) > 0
    Lon <- as.numeric(sub('W','',sub('E','',as.character(z$V8))))/10
    Lat <- as.numeric(sub('S','',sub('N','',as.character(z$V7))))/10
    if (south) Lat <- -Lat
    if (west) Lon <- -Lon
    np <- length(Lon)
    lons <- approx(1:np,Lon,1:n)$y
    lats <- approx(1:np,Lat,1:n)$y
    slps <- approx(1:np,z$V9,1:n)$y
    Y[i,] <- c(lons,lats,slps,start,end,n)
  }
  Y[Y <= -999] <- NA
  ok <- is.finite(rowMeans(Y))
  Y <- Y[ok,]
  colnames(Y) <- c(rep('lon',n),rep('lat',n),rep('slp',n),'start','end','n')
  
  # add attributes to trajectory matrix X
  attr(Y, "location")= NA
  attr(Y, "variable")= 'storm tracks'
  attr(Y, "longname")= 'Tropical cyclone storm tracks'
  attr(Y, "quality")= NA
  attr(Y, "calendar")= "gregorian"
  attr(Y, "source")= 'US NAVY JTWC'
  attr(Y, "URL")= url
  attr(Y, "unit")= NA
  attr(Y, "type")= "analysis"
  attr(Y, "aspect")= "interpolated"
  attr(Y, "reference")= ''
  attr(Y, "info")= 'http://www.usno.navy.mil/NOOC/nmfc-ph/RSS/jtwc/best_tracks/TC_bt_report.html'
  attr(Y, "method")= NA
  attr(Y,"lon") <- NA
  attr(Y,"lat") <- NA
  attr(Y,"alt") <- NA
  attr(Y,"cntr") <- NA
  attr(Y,"stid") <- NA
  attr(Y,'domain') <- switch(domain,
                             'io'='Northern Indian Ocean',
                             'sh'='Southern Hemisphere','wp'='Northwestern Pacific')
  attr(Y, "history")= history.stamp()
  attr(Y,'storm name') <- storms
  class(Y) <- c('trajectory','matrix')
  invisible(Y)
  return(Y)
}
