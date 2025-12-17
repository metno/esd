# enclosing_circle <- function(xlim, ylim, lonR, latR, axiR, 
#                              stereographic = FALSE, mask_horizon=FALSE) {
#   
#   # Create a grid of boundary points (lon, lat) based on xlim and ylim
#   lon_bounds <- seq(xlim[1], xlim[2], length.out = 100)  # Longitudes along the edges
#   lat_bounds <- seq(ylim[1], ylim[2], length.out = 100)  # Latitudes along the edges
#   
#   # Combine lon-lat into a grid of boundary points
#   lon_grid <- c(lon_bounds, rep(xlim[1], length(lat_bounds)), rep(xlim[2], length(lat_bounds)))
#   lat_grid <- c(rep(ylim[1], length(lon_bounds)), lat_bounds, lat_bounds)
#   
#   # Project the boundary points using your existing function
#   projection <- cartesian2sphere(lon_grid, lat_grid, lonR = lonR, latR = latR, axiR=axiR,
#                                  stereographic = stereographic, mask_horizon = mask_horizon)
#   
#   # Extract projected coordinates
#   x <- projection$X
#   y <- projection$Y
#   
#   # Calculate the maximum distance from the center to any boundary point
#   distances <- sqrt(x^2 + y^2)
#   radius <- max(distances, na.rm = TRUE)
#   
#   # Generate points for the enclosing circle
#   angles <- seq(0, 2 * pi, length.out = 500)
#   x_circle <- radius * cos(angles)
#   y_circle <- radius * sin(angles)
#   
#   return(data.frame(x = x_circle, y = y_circle))
# }


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
#' @param gridlines Default = TRUE - for adding gridlines on map
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
#' @export station2sphere
station2sphere <- function(x,FUN="mean",lonR=NULL,latR=NULL,axiR=0,
                   stereographic=FALSE, mask_horizon=TRUE,
                   xlim=NULL,ylim=NULL,gridlines=TRUE,cex=0.2,
                   cex.axis=1,cex.lab=1,cex.main=1.5,pch=21,
                   colbar= list(pal='t2m',col=NULL,rev=FALSE,n=10,
                                breaks=NULL,type="p",cex=2,h=0.6, v=1,
                                pos=0.1,show=TRUE),
                   new=TRUE,add=FALSE,verbose=FALSE,...) {
  if(verbose) print(paste('station2sphere:',lonR,latR))
  x0 <- x
  
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
  
  ## KMP 2025-11-12: Transform to Greenwich format (to avoid xlim problems)
  greenwich <- (all(x > 0) & any(x > 180))
  
  # Default xlim/ylim to global if not provided
  if(is.null(xlim)) xlim <- c(-180, 180)
  if(is.null(ylim)) ylim <- c(-90, 90)
  
  if(greenwich) xlim[xlim < 0] <- xlim[xlim < 0] + 360 else 
    xlim[xlim > 180] <- xlim[xlim > 180] - 360
  
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
  
  # Center of projection
  if(is.null(axiR)) axiR <- 0
  if(is.null(lonR)) {
    if(!is.null(xlim)) {
      lonR <- mean(xlim, na.rm=TRUE)
    } else {
      lonR <- mean(lon, na.rm=TRUE)
    }
  }
  if (is.null(latR)) {
    if(!is.null(ylim)) {
      latR <- mean(ylim, na.rm=TRUE) 
    } else if (is.null(latR)) {
      latR <- mean(lat, na.rm=TRUE)
    }
  } else if(latR > 0) {
    latR <- min(latR, 90)
  } else {
    latR <- max(latR, -90)
  }
  
  if(verbose) print(paste("lonR", lonR, "latR", latR))
  ## Avoid singularity at poles
  pole_offset <- 1e-3
  if (ylim[1] <= (-90 + pole_offset)) ylim[1] <- -90 + pole_offset
  if (ylim[2] >= (90 - pole_offset)) ylim[2] <- 90 - pole_offset
  # Restrict ylim to visible latitudes
  if (mask_horizon) {
    phi_S_horizon <- latR - 90 + pole_offset # most southerly visible latitude
    phi_N_horizon <- latR + 90 - pole_offset # most northerly visible latitude
    if (ylim[1] < phi_S_horizon) flag_adjusted_ylim1 <- TRUE
    if (ylim[2] > phi_N_horizon) flag_adjusted_ylim2 <- TRUE
    ylim[1] <- max(-90 + pole_offset, max(phi_S_horizon, ylim[1]))
    ylim[2] <- min(90 - pole_offset, min(phi_N_horizon, ylim[2]))
  }

  if (verbose) print("Calculating plot limits from geographical bounds")
  # Define the geographical boundaries of the plot area
  n_pts_boundary <- 100
  lon_x_min <- rep(xlim[1], n_pts_boundary) # Left side longitude
  lon_x_max <- rep(xlim[2], n_pts_boundary) # Right side longitude
  lat_y_min <- rep(ylim[1], n_pts_boundary) # Bottom side latitude
  lat_y_max <- rep(ylim[2], n_pts_boundary) # Top side latitude
  
  lat_seq <- seq(ylim[1], ylim[2], length.out = n_pts_boundary)
  lon_seq <- seq(xlim[1], xlim[2], length.out = n_pts_boundary)
  
  lat_bounds <- sapply(lat_seq, function(lat) rep(lat, length(lon_seq)))
  lon_bounds <- rep(lon_seq, length(lat_seq))
  #lon_bounds <- c(lon_seq, lon_seq, lon_x_min, lon_x_max)
  #lat_bounds <- c(lat_y_min, lat_y_min, lat_seq, lat_seq)
  
  # Project the boundaries
  proj_bounds <- cartesian2sphere(lon_bounds, lat_bounds, lonR = lonR, latR = latR, axiR=axiR,
                                  stereographic = stereographic, mask_horizon = mask_horizon)
  
  # Calculate the projected Xlim and Ylim (only from visible points)
  Xlim <- range(proj_bounds$X[proj_bounds$visible], na.rm = TRUE)
  Ylim <- range(proj_bounds$Y[proj_bounds$visible], na.rm = TRUE)
  
  # Process coastline data (geoborders)
  geoborders <- NULL
  data("geoborders",envir=environment())
  
  # Transform from dateline to greenwich longitude format
  gx <- geoborders$x
  if(greenwich) gx[!is.na(gx) & gx < 0] <- gx[!is.na(gx) & gx < 0] + 360 else
    gx[!is.na(gx) & gx > 180] <- gx[!is.na(gx) & gx > 180] - 360
  gy <- geoborders$y
  ok <- is.finite(gx) & is.finite(gy) # Only filter NAs initially
  
  # Project ALL finite geoborders data
  xy_geoborders_all <- cartesian2sphere(gx[ok], gy[ok], lonR=lonR, latR=latR, axiR=axiR,
                                        stereographic=stereographic, mask_horizon=mask_horizon)
  
  # Apply clipping masks based on xlim and ylim
  # Projected horizon mask (Z' >= 0)
  mask_horizon_geoborders <- xy_geoborders_all$visible
  # Geographical longitude mask (xlim cut)
  mask_xlim_geoborders <- gx[ok] >= min(xlim) & gx[ok] <= max(xlim)
  # Geographical latitude mask (ylim cut)
  mask_ylim_geoborders <- gy[ok] >= min(ylim) & gy[ok] <= max(ylim)
  # Combine all three masks
  mask_final_geoborders <- mask_horizon_geoborders & mask_xlim_geoborders & mask_ylim_geoborders
  # Apply the final combined mask
  x_geoborders <- xy_geoborders_all$X[mask_final_geoborders]
  y_geoborders <- xy_geoborders_all$Y[mask_final_geoborders]
  
  # Process station data
  xy_stations <- cartesian2sphere(lon, lat, lonR=lonR, latR=latR, axiR=axiR,
                                  stereographic=stereographic, mask_horizon=mask_horizon)
  X_all <- xy_stations$X
  Y_all <- xy_stations$Y
  Visible_horizon <- xy_stations$visible
  
  # Create geographical masks for stations
  mask_xlim_stations <- lon >= min(xlim) & lon <= max(xlim)
  mask_ylim_stations <- lat >= min(ylim) & lat <= max(ylim)
  
  # Combine all three masks for stations
  Visible_final <- Visible_horizon & mask_xlim_stations & mask_ylim_stations
  
  ## KMP 2025-10-09: Adding space under the map to fit the color bar
  dy <- 0.3*diff(Ylim)
  if(colbar$show) Ylim_adjusted <- Ylim + c(-1,0)*dy else 
    Ylim_adjusted <- Ylim + c(-0.1,0)*dy
  
  # Plotting
  if(new) dev.new()
  par(bty="n",xaxt="n",yaxt="n",new=add)
  plot(Xlim, Ylim, pch=".", col="white", xlab="", ylab="",
       xlim=Xlim, ylim=Ylim_adjusted, cex.axis=cex.axis, cex.lab=cex.lab)
  
  # Plot Station Data
  points(X_all[Visible_final], Y_all[Visible_final], cex=cex, pch=pch,
         col=col[Visible_final], bg=bg[Visible_final])
  
  # Plot Coastlines
  points(x_geoborders, y_geoborders, pch=".")
  
  # Plot bounding box
  if (verbose) print("Drawing projected box boundary")
  n_points <- 200
  lat_seq_boundary <- seq(ylim[1], ylim[2], length.out = n_points)
  lon_seq_boundary <- seq(xlim[1], xlim[2], length.out = n_points)
  
  if (verbose) {
    print("Projected boundaries and visibility:")
    print(range(proj_bounds$X))
    print(range(proj_bounds$Y))
    print(range(proj_bounds$visible))
  }
  if(diff(range(xlim)) > 350) {
    # Bounding circle
    distances <- sqrt(proj_bounds$X[proj_bounds$visible]^2 + proj_bounds$Y[proj_bounds$visible]^2)
    radius <- max(distances, na.rm = TRUE)
    # Generate points for the enclosing circle
    angles <- seq(0, 2 * pi, length.out = 500)
    x_circle <- radius * cos(angles)
    y_circle <- radius * sin(angles)
    lines(x_circle, y_circle)
  } else {
    # Left boundary
    b1 <- cartesian2sphere(rep(xlim[1], n_points), lat_seq_boundary, lonR, latR, axiR, stereographic, mask_horizon)
    lines(b1$X[b1$visible], b1$Y[b1$visible], col = "black", lwd = 1)
    # Right boundary
    b2 <- cartesian2sphere(rep(xlim[2], n_points), lat_seq_boundary, lonR, latR, axiR, stereographic, mask_horizon)
    lines(b2$X[b2$visible], b2$Y[b2$visible], col = "black", lwd = 1)

    # Bottom and top boundary
    if(ylim[1]!=latR & ylim[1] > (-90 + pole_offset)) {
      b3 <- cartesian2sphere(lon_seq_boundary, rep(ylim[1], n_points), lonR, latR, axiR, stereographic, mask_horizon)
      lines(b3$X[b3$visible], b3$Y[b3$visible], col = "black", lwd = 1)
    } 
    if(ylim[2]!=latR & ylim[2] < (90 - pole_offset)) {
      b4 <- cartesian2sphere(lon_seq_boundary, rep(ylim[2], n_points), lonR, latR, axiR, stereographic, mask_horizon)
      lines(b4$X[b4$visible], b4$Y[b4$visible], col = "black", lwd = 1)
    }
  }
  
  # Add colorbar
  if (!is.null(FUN) & colbar$show) {
    label <- generate_varlabel(x0)
    xlim_usr <- par()$usr[1:2]
    ylim_usr <- par()$usr[3:4]
    dy_total <- diff(ylim_usr)*0.1
    below <- c(min(xlim_usr), min(ylim_usr)+dy_total, max(xlim_usr), min(ylim_usr)+2*dy_total)
    dy_below <- below[4]-below[2]
    
    # Masking rectangle (if needed for background)
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
  
  # Return results
  if (inherits(x0,"stationmeta")) {
    result <- data.frame(x=X_all[Visible_final], y=Y_all[Visible_final])
  } else if (inherits(x0,"station")) {
    result <- data.frame(x=X_all[Visible_final], y=Y_all[Visible_final], 
                         z=as.vector(map[Visible_final]))
  }
  invisible(result)
}
