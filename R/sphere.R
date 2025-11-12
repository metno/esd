enclosing_circle <- function(xlim, ylim, lonR, latR, axiR, 
                             stereographic = FALSE, mask_horizon=FALSE) {
  
  # Create a grid of boundary points (lon, lat) based on xlim and ylim
  lon_bounds <- seq(xlim[1], xlim[2], length.out = 100)  # Longitudes along the edges
  lat_bounds <- seq(ylim[1], ylim[2], length.out = 100)  # Latitudes along the edges
  
  # Combine lon-lat into a grid of boundary points
  lon_grid <- c(lon_bounds, rep(xlim[1], length(lat_bounds)), rep(xlim[2], length(lat_bounds)))
  lat_grid <- c(rep(ylim[1], length(lon_bounds)), lat_bounds, lat_bounds)
  
  # Project the boundary points using your existing function
  projection <- cartesian2sphere(lon_grid, lat_grid, lonR = lonR, latR = latR, axiR=axiR,
                                 stereographic = stereographic, mask_horizon = mask_horizon)
  
  # Extract projected coordinates
  x <- projection$X
  y <- projection$Y
  
  # Calculate the maximum distance from the center to any boundary point
  distances <- sqrt(x^2 + y^2)
  radius <- max(distances, na.rm = TRUE)
  
  # Generate points for the enclosing circle
  angles <- seq(0, 2 * pi, length.out = 500)
  x_circle <- radius * cos(angles)
  y_circle <- radius * sin(angles)
  
  return(data.frame(x = x_circle, y = y_circle))
}


#' Map station data on a spherical projection
#'
#' The function sphere is used inside map.station to plot station data on a spherical projection 
#'
#' @seealso map.station
#' 
#' @param x an input object of class 'station'
#' @param FUN a function, see \code{\link{aggregate.zoo}}. Default: "mean"
#' @param lonR Center longitude of viewing angle
#' @param latR Center latitude of viewing angle
#' @param stereographic if TRUE return a stereographic projection, else skip the last projection step and return rotated Cartesian coordinates
#' @param mask_horizon if TRUE mask areas beyond the horizon so that only one hemisphere can be shown at once. This is only 
#'                     relevant to the stereographic projection which can expand beyond the hemisphere and show the whole globe
#'                     (with a distorted perspective far from the center)  
#' @param stereographic If TRUE use a stereographic projection, otherwise use rotated Cartesian coordinates
#' @param verbose a boolean; if TRUE print information about progress
#' @param xlim range of map in the longitude direction 
#' @param ylim range of map in the latitude direction
#' @param gridlines=TRUE
#' @param cex magnification factor for point symbols
#' @param cex.axis magnification factor for axis text
#' @param cex.lab magnification factor for label text
#' @param cex.main magnification factor for main title
#' @param pch point symbol
#' @param colbar List containing colorbar specifications: 
#'    list(pal='t2m',rev=FALSE,n=10,col=NULL,breaks=NULL,type="p",cex=2,h=0.6,v=1,pos=0.1,show=TRUE).
#'    See \link{col.bar} for details.  
#' @param new If new create new plot window
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional argument
#' 
#' @return A visualization of the data on spherical coordinates and a data frame 
#         with coordinates on the rotated projection
#'
#' @keywords utilities
#' 
#' @export sphere
sphere <- function(x,FUN="mean",lonR=NULL,latR=NULL,axiR=0,
                   stereographic=FALSE, mask_horizon=TRUE,
                   xlim=NULL,ylim=NULL,gridlines=TRUE,cex=0.2,
                   cex.axis=1,cex.lab=1,cex.main=1.5,pch=21,
                   colbar= list(pal='t2m',col=NULL,rev=FALSE,n=10,
                                breaks=NULL,type="p",cex=2,h=0.6, v=1,
                                pos=0.1,show=TRUE),
                   new=TRUE,verbose=FALSE,...) {
  if(verbose) print(paste('sphere:',lonR,latR))
  x0 <- x
  
  ## KMP 2016-12-21: To handle xlim in greenwich format, e.g., 180-360
  if(!is.null(xlim)) {
    greenwich <- (min(xlim)>0 & max(xlim)>180)
    g2dl(x,greenwich=greenwich,verbose=verbose)
  } else greenwich <- NULL
  
  ## Data to be plotted:
  if (inherits(x,"stationmeta")) {
    lon <- x$longitude
    lat <- x$latitude
    param <- param2ele(x$ele)
    unit <- " "  
  } else if (inherits(x,"station")) {
    lon <- attr(x,'longitude')
    lat <- attr(x,'latitude')
    param <- as.character(levels(factor(attr(x,'parameter'))))
  }
  ## KMP 2016-12-21: To handle xlim in greenwich format, e.g., 180-360
  if(is.null(greenwich)) greenwich <- (min(lon)>0 & max(lon)>180)
  
  ## To deal with grid-conventions going from north-to-south or east-to-west:
  ##srtx <- order(attr(x,'longitude')); lon <- lon[srtx]
  ##srty <- order(attr(x,'latitude')); lat <- lat[srty]
  
  if (!is.null(FUN)) {
    map <- apply(as.matrix(x),2,FUN,na.rm=TRUE) ##map <- x[srtx,srty]
  } else {
    #map <- x
    if(length(x)==length(lon(x))) map <- x else map <- rep(1, length(lon(x)))
  }
  
  ## KMP 2025-11-11: the input col and bg are not used at all in sphere
  ## Is this what we want? It's the same in map2sphere
  ## Initialise colbar
  colbar <- colbar.ini(map, colbar=colbar, verbose=verbose)
  breaks <- colbar$breaks
  colb <- colbar$col
  col <- colb[findInterval(map, breaks)]
  bg <- col
  nc <- length(colb)
  
  # Center of projection
  if(is.null(axiR)) axiR <- 0
  if(!is.null(xlim)) {
    lonR <- mean(xlim, na.rm=TRUE)
  } else if (is.null(lonR)) {
    lonR <- mean(lon, na.rm=TRUE)  
  }
  if (is.null(latR)) {
    if(!is.null(ylim)) {
      latR <- mean(ylim, na.rm=TRUE) 
    } else if (is.null(latR)) {
      latR <- mean(lat, na.rm=TRUE)
    }
  }
  
  # Rotate xlim and ylim
  if(is.null(xlim)) xlim <- c(-180, 180)
  if(is.null(ylim)) {
    if((latR - 45) < -90) ylim <- c(-90, 0) else
      if((latR + 45) > 90) ylim <- c(0, 90) else
        ylim <- c(latR - 45, latR + 45)
  }
  if(is.null(axiR)) axiR <- 0
  
  ## KMP 2025-10-09: Create a lon/lat grid from xlim/ylim, then transform 
  ##  the grid to stereographic coordinates and get Xlim/Ylim based on the 
  ##  range of the rotated grid
  circle <- enclosing_circle(xlim, ylim, lonR, latR, axiR, 
               stereographic = stereographic, mask_horizon=mask_horizon)
  Xlim <- range(circle$x, na.rm=TRUE)
  Ylim <- range(circle$y, na.rm=TRUE)
  
  ## KMP 2025-10-09: Adding space under the map to fit the color bar
  dy <- 0.3*diff(Ylim)
  if(colbar$show) Ylim <- Ylim + c(-1,0)*dy else Ylim <- Ylim + c(-0.1,0)*dy
  
  # coastline data:
  geoborders <- NULL # KMP 2019-10-11: create dummy to avoid warning during CHECK
  data("geoborders",envir=environment())
  ## If the input data is in dateline format (longitude 0 - 360)
  ## the geoborders longitude also needs to be transformed to avoid problems
  if(max(lon)>180 & min(lon)>=0) {
    gx <- geoborders$x
    gx[!is.na(gx) & gx < 0] <- gx[!is.na(gx) & gx < 0] + 360
    geoborders$x <- gx
  }
  
  gx <- geoborders$x
  gy <- geoborders$y
  ok <- is.finite(gx) & is.finite(gy)
  ## KMP 2024-08-09: Apply xlim and ylim to the geoborders
  if (!is.null(xlim)) ok <- ok & gx>=min(xlim) & gx<=max(xlim)
  if (!is.null(ylim)) ok <- ok & gy>=min(ylim) & gy<=max(ylim)
  
  ## KMP 2025-02-06: Calculating spherical coordinates with a separate function
  xy_geoborders <- cartesian2sphere(gx[ok], gy[ok], lonR=lonR, latR=latR, axiR=axiR,
                                    stereographic=stereographic, mask_horizon=mask_horizon)
  x_geoborders <- xy_geoborders$X[xy_geoborders$visible]
  y_geoborders <- xy_geoborders$Y[xy_geoborders$visible]
  
  ## KMP 2025-02-11: Rotate coordinates using the function cartesian2sphere
  xy_stations <- cartesian2sphere(lon, lat, lonR=lonR, latR=latR, axiR=axiR,
                                  stereographic=stereographic, mask_horizon=mask_horizon)
  X <- xy_stations$X
  Y <- xy_stations$Y
  Visible <- xy_stations$visible
  
  ## Define colour pal:
  if (!is.null(FUN)) {
    if(!is.null(colbar$n)) n <- colbar$n else n <- 10
    if(!is.null(colbar$pal)) pal <- colbar$pal else varid(x0)[1]
    if (is.null(col)) col <- colscal(n=n, pal=pal) else
      if (length(col)==1) {
        col <- colscal(pal=col, n=n)
      }
    nc <- length(col)
    index <- round( nc*( map - min(map) )/
                      ( max(map) - min(map) ) )
  }
  
  # Plot the results:
  if(new) dev.new()
  par(bty="n",xaxt="n",yaxt="n",new=TRUE)
  plot(Xlim, Ylim, pch=".", col="white", xlab="", ylab="",
       Xlim=Xlim, Ylim=Ylim, cex.axis=cex.axis, cex.lab=cex.lab)
  
  points(X[Visible], Y[Visible], cex=cex, pch=pch, 
         col=col[Visible], bg=bg[Visible])
  ## Add contour lines
  ## Plot the coast lines
  points(x_geoborders, y_geoborders, pch=".")
  ## Add circle around globe. How does this work in stereographic coordinates?
  lines(x = circle$x, y = circle$y, col = "black", lwd = 2)
  
  # Colorbar
  if (!is.null(FUN) & colbar$show) {      
    ## Generate a label (variable name and unit)
    label <- generate_varlabel(x0)
    ## Where to place colorbar
    xlim <- par()$usr[1:2]
    ylim <- par()$usr[3:4]
    dy <- diff(ylim)*0.1
    #below <- c(min(xlim), min(ylim)-dy/2, max(xlim), min(ylim)+dy/2)
    below <- c(min(xlim), min(ylim)+dy, max(xlim), min(ylim)+2*dy)
    dy_below <- below[4]-below[2]
    rect(below[1], below[2], 
         below[3], below[4]-dy_below*0.2, 
         col = "white", border = "white")
    ## Add colorbar and title
    col.bar(below[1],below[2]+dy_below*0.1,below[3],below[4]-dy_below*0.1,
            colbar$breaks,horiz=TRUE,pch=15,v=1,h=1,
            col=colbar$col,cex=2,cex.lab=colbar$cex.lab,
            type=colbar$type,verbose=FALSE,vl=1,border=FALSE)
    title(sub = label, line = -1, cex.sub = cex.lab)
  }
  
  if (inherits(x0,"stationmeta")) {
    result <- data.frame(x=X[Visible], y=Y[Visible])
  } else if (inherits(x0,"station")) {
    result <- data.frame(x=X[Visible], y=Y[Visible], z=map[Visible])
  }
  invisible(result)
}


# OLD SPHERE FUNCTION
# sphere <- function(x,n=30,FUN="mean",lonR=10,latR=45,axiR=0,xlim=NULL,ylim=NULL,
#                    gridlines=TRUE,col="green",bg="darkgreen",cex=0.2,
#                    cex.axis=1,cex.lab=1,cex.main=1.5,pch=".",
#                    colbar= list(pal='t2m',col=NULL,rev=FALSE,n=10,
#                                 breaks=NULL,type="p",cex=2,h=0.6, v=1,
#                                 pos=0.1,show=TRUE),
#                    new=TRUE,verbose=FALSE,...) {
#   if(verbose) print(paste('sphere:',lonR,latR))
#   x0 <- x
#   
#   ## KMP 2016-12-21: To handle xlim in greenwich format, e.g., 180-360
#   if(!is.null(xlim)) {
#     greenwich <- (min(xlim)>0 & max(xlim)>180)
#     g2dl(x,greenwich=greenwich,verbose=verbose)
#   } else greenwich <- NULL
#   
#   ## Data to be plotted:
#   if (inherits(x,"stationmeta")) {
#     lon <- x$longitude
#     lat <- x$latitude
#     param <- param2ele(x$ele)
#     unit <- " "  
#   } else if (inherits(x,"station")) {
#     lon <- attr(x,'longitude')
#     lat <- attr(x,'latitude')
#     param <- as.character(levels(factor(attr(x,'parameter'))))
#   }
#   ## KMP 2016-12-21: To handle xlim in greenwich format, e.g., 180-360
#   if(is.null(greenwich)) greenwich <- (min(lon)>0 & max(lon)>180)
#   
#   ## To deal with grid-conventions going from north-to-south or east-to-west:
#   ##srtx <- order(attr(x,'longitude')); lon <- lon[srtx]
#   ##srty <- order(attr(x,'latitude')); lat <- lat[srty]
#   
#   if (!is.null(FUN)) {
#     map <- apply(as.matrix(x),2,FUN,na.rm=TRUE) ##map <- x[srtx,srty]
#   } else {
#     #map <- x
#     if(length(x)==length(lon(x))) map <- x else map <- rep(1, length(lon(x)))
#   }
#   
#   # Rotatio:
#   # longitudinal rotation
#   if(!is.null(xlim)) {
#     lonR <- mean(xlim, na.rm=TRUE)
#   } else if (is.null(lonR)) {
#     lonR <- mean(lon, na.rm=TRUE)  
#   }
#   # latitudinal rotation
#   if(!is.null(ylim)) {
#     latR <- mean(ylim, na.rm=TRUE) 
#   } else if (is.null(latR)) {
#     latR <- mean(lat, na.rm=TRUE)
#   }
#   # axiR: rotation of Earth's axis
#   
#   # coastline data:
#   data("geoborders",envir=environment())
#   ## KMP 10-11-2015: apply xlim and ylim
#   gx <- geoborders$x
#   gy <- geoborders$y
#   ok <- is.finite(gx) & is.finite(gy)
#   if(greenwich) gx[gx<0 & ok] <- gx[gx<0 & ok] + 360
#   ## KMP 2025-02-11: Rotate coordinates using the function cartesian2sphere
#   xyz_geoborders <- cartesian2sphere(gx[ok], gy[ok], lonR=lonR, latR=latR)
#   x <- xyz_geoborders$X
#   y <- xyz_geoborders$Y
#   z <- xyz_geoborders$Z
#   #theta <- pi*gx[ok]/180
#   #phi <- pi*gy[ok]/180
#   #ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
#   #theta <- pi*geoborders$x[ok]/180; phi <- pi*geoborders$y[ok]/180
#   #x <- sin(theta)*cos(phi)
#   #y <- cos(theta)*cos(phi)
#   #z <- sin(phi)
#   
#   ## KMP 2025-02-11: Rotate coordinates using the function cartesian2sphere
#   xyz_stations <- cartesian2sphere(lon, lat, lonR=lonR, latR=latR)
#   X <- xyz_stations$X
#   Y <- xyz_stations$Y
#   Z <- xyz_stations$Z
#   #Theta <- pi*lon/180; Phi <- pi*lat/180
#   #
#   ## Transform -> (X,Y,Z):
#   #X <- sin(Theta)*cos(Phi)
#   #Y <- cos(Theta)*cos(Phi)
#   #Z <- sin(Phi)
#   #print(c( min(x),max(x)))
#   
#   ## Define colour pal:
#   if (!is.null(FUN)) {
#     if (is.null(col)) col <- colscal(n=n,pal=varid(x)) else
#       if (length(col)==1) {
#         col <- colscal(pal=col,n=n)
#       }
#     nc <- length(col)
#     index <- round( nc*( map - min(map) )/
#                       ( max(map) - min(map) ) )
#   }
#   ## KMP 2025-02-11: The coordinates are transformed and rotated above
#   # Rotate coastlines:
#   #a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
#   #a <- rotM(x=latR,y=0,z=0) %*% a
#   #x <- a[1,]; y <- a[2,]; z <- a[3,]
#   
#   # Grid coordinates:
#   #d <- dim(X)
#   #print(d)
#   
#   ## KMP 2025-02-11: The coordinates are transformed and rotated above
#   # Rotate data grid:  
#   #A <- rotM(x=0,y=0,z=lonR) %*% rbind(c(X),c(Y),c(Z))
#   #A <- rotM(x=latR,y=0,z=0) %*% A
#   #X <- A[1,]; Y <- A[2,]; Z <- A[3,]
#   #dim(X) <- d; dim(Y) <- d; dim(Z) <- d
#   #print(dim(rbind(X,Z)))
#   
#   # Rotate xlim and ylim
#   if(!is.null(xlim) & !is.null(ylim)) {
#     xyz_lim <- cartesian2sphere(xlim, ylim, lonR=lonR, latR=latR)
#     Xlim <- xyz_lim$X
#     Ylim <- xyz_lim$Y
#     Zlim <- xyz_lim$Z
#     #thetalim <- pi*xlim/180
#     #philim <- pi*ylim/180
#     #Xlim <- sin(thetalim)*cos(philim)
#     #Ylim <- cos(thetalim)*cos(philim)
#     #Zlim <- sin(philim)
#     #Alim <- rotM(x=0,y=0,z=lonR) %*% rbind(c(Xlim),c(Ylim),c(Zlim))
#     #Alim <- rotM(x=latR,y=0,z=0) %*% Alim
#     #Xlim <- Alim[1,]; Ylim <- Alim[2,]; Zlim <- Alim[3,]
#   } else {
#     Xlim <- range(x, na.rm=TRUE)
#     Zlim <- range(z, na.rm=TRUE)
#   }
#   
#   # Plot the results:
#   if(new) dev.new()
#   par(bty="n",xaxt="n",yaxt="n",new=TRUE)
#   plot(Xlim,Zlim,pch=".",col="white",xlab="",ylab="",
#        cex.axis=cex.axis,cex.lab=cex.lab)
#   #plot(x,z,pch=".",col="white",xlab="",ylab="",
#   #     cex.axis=cex.axis,cex.lab=cex.lab)
#   #par0 <- par()
#   
#   # plot the grid boxes, but only the gridboxes facing the view point:
#   ##Visible <- Y > 0 ##colMeans(Y) > 0
#   ##X <- X[,Visible]; Y <- Y[,Visible]; Z <- Z[,Visible]
#   ##index <- index[Visible]
#   ##apply(rbind(X,Z,index),2,gridbox,cols)
#   # c(W,E,S,N, colour)
#   # xleft, ybottom, xright, ytop
#   
#   ##if (!is.null(FUN)) {
#   ##    colb <- colscal(n=length(breaks)) 
#   ##    col <- colb[findInterval(map,breaks)]
#   ##    bg <- col
#   ##    nc <- length(colb)
#   ##}
#   
#   ## Initialise colbar
#   if (verbose) {cat('Plotting'); str(colbar)}
#   colbar <- colbar.ini(map, colbar=colbar, verbose=verbose)
#   #colbar <- colbar.ini(map, FUN=FUN)
#   breaks <- colbar$breaks
#   colb <- colbar$col
#   col <- colb[findInterval(map, breaks)]
#   bg <- col
#   nc <- length(colb)
#   visible <- Y > 0
#   points(X[visible],Z[visible],cex=cex,pch=pch,col=col,bg=bg)
#   
#   ## Add contour lines?
#   ## Plot the coast lines  
#   visible <- y > 0
#   points(x[visible],z[visible],pch=".")
#   #plot(x[visible],y[visible],type="l",xlab="",ylab="")
#   if(is.null(xlim) & is.null(ylim)) lines(cos(pi/180*1:360),sin(pi/180*1:360),col="black")
#   ## Add grid ?
#   
#   # Colorbar
#   if (!is.null(FUN) & colbar$show) {      
#     ## Generate a label (variable name and unit)
#     label <- generate_varlabel(x0)
#     ## Where to place colorbar
#     ylim <- par()$usr[1:2]
#     ylim <- par()$usr[3:4]
#     dy <- diff(ylim)*0.1
#     below <- c(min(xlim), min(ylim)-dy/2, max(xlim), min(ylim)+dy/2)
#     dy_below <- below[4]-below[2]
#     rect(below[1], below[2], 
#          below[3], below[4]-dy_below*0.2, 
#          col = "white", border = "white")
#     ## Add colorbar and title
#     col.bar(below[1],below[2]+dy_below*0.1,below[3],below[4]-dy_below*0.1,
#             colbar$breaks,horiz=TRUE,pch=15,v=1,h=1,
#             col=colbar$col,cex=2,cex.lab=colbar$cex.lab,
#             type=colbar$type,verbose=FALSE,vl=1,border=FALSE)
#     title(sub = label, line = 1, cex.sub = cex.lab)
#     #par(fig = c(0.3, 0.7, 0.05, 0.10),mar=rep(0,4),cex=0.8,
#     #    new = TRUE, mar=c(1,0,0,0), xaxt = "s",yaxt = "n",bty = "n")
#     #print("colourbar")
#     ##breaks <- round( nc*(seq(min(map),max(map),length=nc)- min(map) )/                 ( max(map) - min(map) ) )
#     #bar <- cbind(breaks,breaks)
#     #image(seq(breaks[1],breaks[length(breaks)],length=nc),
#     #      c(1,2),bar,col=col,cex.axis=cex.axis)
#     #
#     #par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
#     #    fig=c(0,1,0,1),new=TRUE)
#     #plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
#     #text(0.1,0.95,param,cex=cex.main,pos=4)
#     #text(0.72,0.002,unit,pos=4)  
#   }
#   ##result <- data.frame(x=colMeans(Y),y=colMeans(Z),z=c(map))
#   if (inherits(x0,"stationmeta")) {
#     result <- data.frame(x=Y,y=Z)
#   } else if (inherits(x0,"station")) {
#     result <- data.frame(x=Y,y=Z,z=map)
#   }
#   invisible(result)
# }