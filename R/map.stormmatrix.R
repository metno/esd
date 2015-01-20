## Author 	 Kajsa Parding
## Last update   20.01.2015
## Includes	 map.stormtracks()
## Require 	 geoborders.rda

map.stormmatrix <- function(x,FUN=NULL,
      lonR=10,latR=90,alpha=0.2,new=TRUE) {

  # Stuff that need to be fixed:
  # 1) Stormtracks that cross the dateline do not look good, are now excluded
  # 2) The color of the tracks should be set by FUN
  # 3) Use a different color palette?
  # 4) Add choice of alternative projection (not polar)
  
  ## Data to be plotted:
  if (inherits(x,"stormmatrix")) {
    x0 <- x
    lons <- x[,colnames(x)=='lon']
    lats <- x[,colnames(x)=='lat']
    t <- strptime(x[,colnames(x)=="start"],format="%Y%m%d%H")
    yr <- as.numeric(strftime(t,"%Y"))

    # Rotatio:
    if (is.null(lonR)) lonR <- mean(lons)
    if (is.null(latR)) latR <- mean(lats)
     
    # coastline data:
    data("geoborders",envir=environment())
    ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  
    # rotate coastline data
    theta <- pi*geoborders$x[ok]/180; phi <- pi*geoborders$y[ok]/180
    x <- sin(theta)*cos(phi)
    y <- cos(theta)*cos(phi)
    z <- sin(phi)
    a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
    a <- rotM(x=latR,y=0,z=0) %*% a
    x <- a[1,]; y <- a[2,]; z <- a[3,]
  
    # rotate stormtrack data
    fn <- function(X) {
      lon <- X[1:10]; lat <- X[11:20]
      theta <- pi*lon/180; phi <- pi*lat/180
      x <- sin(theta)*cos(phi)
      y <- cos(theta)*cos(phi)
      z <- sin(phi)
      a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
      a <- rotM(x=latR,y=0,z=0) %*% a
      x <- a[1,]; y <- a[2,]; z <- a[3,]
      return(c(x,y,z))
    }
    A <- apply(x0,1,fn)
    X <- A[1:10,]; Y <- A[11:20,]; Z <- A[21:30,]
    
    # Define colour palette:
    # (colmap should be filled using FUN but for now 
    # the colorscale will represent the start year of the storm)
    map <- yr
    n <- length(unique(map))
    col <- colscal(n=n)
    col = adjustcolor(col,alpha.f=alpha)
    colmap <- rep(NA,length(map))
    for (i in 1:n) {
      iy <- (map==unique(map)[i])
      colmap[iy] <- col[i]
    }

    if (new) dev.new()
    par(bty="n",xaxt="n",yaxt="n",new=TRUE)
    plot(x,z,pch=".",col="white",xlab="",ylab="")

    # Exclude stormtracks crossing the dateline
    OK <- apply(lons,1,function(x) !(any(x < -90) & any(x > 90)))
    Xv <- X; Xv[Y < 0] <- NA
    Zv <- Z; Zv[Y < 0] <- NA
    matlines(Xv[,OK],Zv[,OK],lty=1,col=colmap[OK])

    # Problems plotting stormtracks that cross the dateline
    #matlines(Xv[,!OK],Zv[,!OK],col=adjustcolor('black',alpha.f=0.4),lty=1)

    ## Plot the coast lines  
    visible <- y > 0
    points(x[visible],z[visible],pch=".")
    lines(cos(pi/180*1:360),sin(pi/180*1:360),col="black")

  }
}
  
