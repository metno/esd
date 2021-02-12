#' Calculate and plot the cyclone budget
#' 
#' Calculate and plot the cyclone budget: total - tracks all grid boxes
#' visited, can be visited many times system - tracks all grid boxes visited,
#' but only the first visit genesis - record the position of the first step
#' lysis - record the position of the last step outN,E,S,W - from where did the
#' cyclone come to the present grid box? inN,E,S,W - to where will the cyclone
#' go from the present grid box?
#' 
#' 
#' @aliases calculate.cyclonebudget plot.cyclonebudget
#' @param traj A 'trajectory' or 'event' object of cyclone trajectories
#' @param it A list or data.frame providing time index, e.g. months, season, year range
#' @param is A list providing space index, e.g.,
#' list(lon=c(-50,50),lat=c(45,70)
#' @param resolution.lon Longitudinal resolution
#' @param resolution.lat Latitudinal resolution
#' @param progress Show progress bar. TRUE or FALSE.
#' @param verbose Print out diagnostics. TRUE or FALSE.
#' @return A 'cyclonebudget' object: a list of various aspects of the cyclone
#' budget.
#' @author K. Parding, MET Norway
#' @keywords cyclonebudget
#' @examples
#' 
#' \dontrun{
#' data(storms)
#' storms.deep <- trackfilter(storms,param="pcent",pmax=970,FUN="any")
#' storms.deep <- trackfilter(storms.deep,param="max.gradient",pmin=2.5e-2,FUN="any")
#' bud <- calculate.cyclonebudget(storms.deep)
#' plot(bud,col=colscal(n=9,pal="bu"))
#' }
#'
#' @export
calculate.cyclonebudget <- function(traj,is=NULL,it=NULL,
                                    resolution.lon=12,resolution.lat=6,
                                    progress=TRUE,verbose=FALSE){

  if(verbose) print("calculate.cyclonebudget")
  stopifnot(inherits(traj,c("trajectory","events")))
  if(inherits(traj,"events")) traj <- as.trajectory(traj)
  # helper functions:
  
  box.coordinates <- function(x,type){
    # from geographic coordinates to rows and columns of the matrix  
    
    if(type=="lon"){
      return( floor((x-minlon)/resolution.lon ) )
    }else{
      return( floor((x-minlat)/resolution.lat ) )
    }
  }
  
  inside <- function(lon0,lat0,lon1,lat1,lons,lats){
    # test if two locations are in the same box
    
    x0 <- box.coordinates(lon0,"lon")
    x1 <- box.coordinates(lon1,"lon")
    
    y0 <- box.coordinates(lat0,"lat")
    y1 <- box.coordinates(lat1,"lat")
    
    if(x0==x1 & y0==y1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  # if wrong class, do something!
  if (inherits(traj,"events")) {
    if(verbose) print("object is of type 'events' - transform to 'trajectory' object")
    traj <- as.trajectory(traj)
  } else if (!inherits(traj,"trajectory")) {
    warning("wrong class for budget")
  }

  ## select a subset of trajectories based on the space index is
  ## (a list with 'lon' and 'lat' ranges) and time index it
  ## (can be a range of years, a season or month ('djf', 'jan'),
  ## or logical or numerical indexing
  if(verbose) {
    print("select subset")
    print(it)
    print(is)
  }
  traj <- subset(traj,is=is,it=it)
  
  # which columns have the longitude (and latitude) values?
  lon.cols <- colnames(traj)=="lon"
  lat.cols <- colnames(traj)=="lat"
  
  # make nice minimum and maximum values for the area we are interested in,
  # get spatial range from is if not null
  if(!is.null(is$lat)) {
    maxlat <- max(is$lat)
    minlat <- min(is$lat)
  } else {
    minlat <- floor(min(traj[,lat.cols])/resolution.lat)*resolution.lat
    maxlat <- ceiling(max(traj[,lat.cols])/resolution.lat)*resolution.lat
  }
  if(!is.null(is$lat)) {
    maxlon <- max(is$lon)
    minlon <- min(is$lon)
  } else {
    minlon <- floor(min(traj[,lon.cols])/resolution.lon)*resolution.lon
    maxlon <- ceiling(max(traj[,lon.cols])/resolution.lon)*resolution.lon
  }
  # latitudes and logitudes of the area boxes
  lats <- seq(minlat,maxlat,by=resolution.lat)
  lons <- seq(minlon,maxlon,by=resolution.lon)
  
  # create a matrix of boxes for "genesis" ...
  genesis <- matrix(0, nrow = length(lons), ncol = length(lats))
  lysis <- genesis
  total <- genesis
  system <- genesis
  inN <- genesis
  inS <- genesis
  inW <- genesis
  inE <- genesis
  outN <- genesis
  outS <- genesis
  outW <- genesis
  outE <- genesis

  # each storm is a row in traj
  if(verbose) print(paste("loop over",nrow(traj),"storms"))
  if(progress) pb <- txtProgressBar(style=3)
  if (progress) setTxtProgressBar(pb,0/nrow(traj))
  for(storm in 1:nrow(traj)){
    # some advanced counter here would be nice...
    X <- NULL
    if(progress) setTxtProgressBar(pb,storm/nrow(traj))
    n <- sum(lon.cols)
    lon <- traj[storm,lon.cols]
    lat <- traj[storm,lat.cols]
    xx <- sapply(lon,box.coordinates,"lon")
    yy <- sapply(lat,box.coordinates,"lat")
    ins <- xx[2:n]==xx[1:(n-1)] & yy[2:n]==yy[1:(n-1)]
    dlon <- lon[2:n] - lon[1:(n-1)]
    dlat <- lat[2:n] - lat[1:(n-1)]
    d <- atan2(dlat,dlon)*180/pi+90
    for (i in seq(1,(n-1))[!ins]) {
      if(d[i] >= 315 | d[i] < 45) {
        inN[xx[i+1],yy[i+1]] <- inN[xx[i+1],yy[i+1]]+1
        outS[xx[i],yy[i]] <- outS[xx[i],yy[i]]+1
      } else if(d[i] >= 45 & d[i] < 135) {
        inW[xx[i+1],yy[i+1]] <- inW[xx[i+1],yy[i+1]]+1
        outE[xx[i],yy[i]] <- outE[xx[i],yy[i]]+1
      } else if(d[i] >= 135 & d[i] < 225) {
        inS[xx[i+1],yy[i+1]] <- inS[xx[i+1],yy[i+1]]+1
        outN[xx[i],yy[i]] <- outN[xx[i],yy[i]]+1
      } else if(d[i] >= 225 & d[i] < 315) {
        inE[xx[i+1],yy[i+1]] <- inE[xx[i+1],yy[i+1]]+1
        outW[xx[i],yy[i]] <- outW[xx[i],yy[i]]+1
      }
    }
  }
  
  lon.total <- traj[,which(lon.cols)]; dim(lon.total) <- length(lon.total)
  lat.total <- traj[,which(lat.cols)]; dim(lat.total) <- length(lat.total)
  xx.total <- sapply(lon.total,box.coordinates,"lon")
  yy.total <- sapply(lat.total,box.coordinates,"lat")
  xx.genesis <- sapply(traj[,which(lon.cols)[1]],box.coordinates,"lon")
  yy.genesis <- sapply(traj[,which(lat.cols)[1]],box.coordinates,"lat")
  xx.lysis <- sapply(traj[,which(lon.cols)[n]],box.coordinates,"lon")
  yy.lysis <- sapply(traj[,which(lat.cols)[n]],box.coordinates,"lat")
  dim(xx.total) <- c(nrow(traj),n)
  dim(yy.total) <- c(nrow(traj),n)
  xy.system <- apply(cbind(xx.total,yy.total),1,
     function(x) unique(cbind(x[1:n],x[(n+1):(2*n)])))
  xx.system <- unlist(lapply(xy.system,function(x) x[,1]))
  yy.system <- unlist(lapply(xy.system,function(x) x[,2]))
  for (i in 1:dim(genesis)[1]) {
    for (j in 1:dim(genesis)[2]) {
      total[i,j] <- sum(xx.total==i & yy.total==j)
      system[i,j] <- sum(xx.system==i & yy.system==j)
      genesis[i,j] <- sum(xx.genesis==i & yy.genesis==j)
      lysis[i,j] <- sum(xx.lysis==i & yy.lysis==j)
    }
  }
  bud=list(total=total,system=system,
           genesis=genesis,lysis=lysis,
           inN=inN,inE=inE,inS=inS,inW=inW,
           outN=outN,outE=outE,outS=outS,outW=outW,lats=lats,lons=lons)
  # making this S3 class might useful, or not? What do you think?
  for(v in names(bud)) {
    x <- bud[[v]]
    dim(x) <- c(1,length(x))
    x <- as.field(zoo(x,order.by=1),lon=bud$lons,lat=bud$lats,
                param=v,unit="cyclones")
    bud[[v]] <- x
  }
  class(bud) <- "cyclonebudget"
  return(bud)
}

#' @exportS3Method
#' @export
plot.cyclonebudget = function(bud,budnames=NULL,new=TRUE,
    colbar=list(pal="precip",n=10,show=TRUE),
    xlim=NULL,ylim=NULL,projection="sphere",
    lonR=0,latR=90,verbose=FALSE,...){

  #x <- unlist(bud[1:(length(bud)-2)])
  #colbar$breaks <- pretty(seq(0,q95(x[x>0], n=colbar$n)))
  #colbar <- colbar.ini(bud,colbar=colbar)
  #col <- colscal(n=length(colbar$breaks),pal=colbar$pal,
  #                   rev=colbar$rev,alpha=colbar$alpha,
  #                   verbose=FALSE)
  #colbar$col <- col[1:(length(col)-1)]
  #colmax <- col[length(col)]

  if(is.null(budnames) | !any(budnames %in% names(bud))) {
    budnames <- names(bud)[!(names(bud) %in% c("lons","lats"))]
  } else {
    budnames <- budnames[budnames %in% names(bud) & !(budnames %in% c("lons", "lats"))]
  }
  nr <- ceiling(length(budnames)/3)
  nc <- ceiling(length(budnames)/nr)
  if(new) dev.new(width=nc*2.2,height=nr*2.3)
  par(mar=c(1.8,1.5,4.0,1.5),mgp=c(0.5,0.5,0.5))
  
  for(i in 1:length(budnames)){
    v <- budnames[i]
    cb <- colbar.ini(bud[[v]],colbar=colbar)
    ir <- nr + 1 - ceiling(i/nc)
    ic <- i - nc*(nr - ir)
    par(new=!(i==1),fig=c((ic-1)/nc,ic/nc,0.95*(ir-1)/nr+0.05,0.95*ir/nr+0.05))
    x <- bud[[v]]
    attr(x,"variable") <- NULL
    map(x,colbar=cb,projection=projection,xlim=xlim,ylim=ylim,
        new=FALSE,lonR=lonR,latR=latR,main=v,adj=1)
  }
}

