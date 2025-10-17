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


sphere <- function(x,n=30,FUN="mean",lonR=NULL,latR=NULL,axiR=0,
                   stereographic=FALSE, mask_horizon=TRUE,
                   xlim=NULL,ylim=NULL,
                   gridlines=TRUE,col="green",bg="darkgreen",cex=0.2,
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
  #xgrid <- seq(min(xlim), max(xlim), diff(range(xlim))/20)
  #ygrid <- seq(min(ylim), max(ylim), diff(range(ylim))/20)
  #nx <- length(xgrid)
  #ny <- length(ygrid)
  #xgrid <- rep(xgrid, ny)
  #ygrid <- as.vector(sapply(ygrid, rep, nx))
  #xy_grid <- cartesian2sphere(xgrid, ygrid, lonR=lonR, latR=latR, axiR=axiR, 
  #                            stereographic=stereographic, mask_horizon=mask_horizon)
  #Xlim <- range(xy_grid$X[xy_grid$visible])
  #Ylim <- range(xy_grid$Y[xy_grid$visible])
  
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
    if (is.null(col)) col <- colscal(n=n,pal=varid(x0)[1]) else
      if (length(col)==1) {
        col <- colscal(pal=col,n=n)
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
