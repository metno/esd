# Documentation in map.R - presents a map on a sphere
#' @export
map2sphere <- function(x,it=NULL,is=NULL,new=TRUE,style="plain",
                       colbar= list(pal='t2m.IPCC',rev=FALSE,n=10,
                           breaks=NULL,type="p",cex=2, cex.axis=0.9,
                           cex.lab=0.9, cex.sub=0.8, h=0.6, v=1,pos=0.05, srt=45),
                       lonR=NULL,latR=NULL,axiR=0, 
                       nx=100, ny=100, type="fill", #c("fill","contour"),
                       col_contour="grey70", breaks_contour=NULL,
                       pos="top",gridlines=TRUE,fancy=TRUE,fig=NULL,add=FALSE,
                       main=NULL,xlim=NULL,ylim=NULL,
                       text.sub=NULL,verbose=FALSE,...) {
  if (verbose) print(paste('map2sphere:',lonR,latR,axiR))
  if (verbose) {print(lon(x)); print(lat(x))}
  if (!is.null(it) | !is.null(is)) x <- subset(x,it=it,is=is,verbose=verbose)
  
  x0 <- x
  ## KMP 2024-08-09: Apply xlim and ylim to the data by using subset
  if (!is.null(xlim) | !is.null(ylim)) x <- subset(x, is=list(lon=xlim, lat=ylim))
  
  # Data to be plotted:
  lon <- lon(x)  # attr(x,'longitude')
  lat <- lat(x)  # attr(x,'latitude')
  if (verbose) {print(summary(lon)); print(summary(lat))}
  # To deal with grid-conventions going from north-to-south or east-to-west:
  if (inherits(x,'field')) map <- map(x,plot=FALSE) else map <- x
  if (verbose) print(class(map))
  srtx <- order(lon); lon <- lon[srtx]
  srty <- order(lat); lat <- lat[srty]
  if (length(dim(map))!= 2) dim(map) <- c(length(srtx),length(srty))
  map <- map[srtx,srty]
  param <- attr(x,'variable')
  unit <- attr(x,'unit')[1]
  if (!is.null(unit) & !is.expression(unit)) if (unit =='%') unit <- "'%'"
  
  ## KMP 10-11-2015: prepare unit and parameter labels
  if(!is.null(param) & !inherits(param,'expression'))
    param <- gsub(" ","~",param)
  if(!is.null(unit) & !inherits(param,'expression'))
    unit <- gsub(" ","~",unit)
  if(length(param)>1) param <- param[1]
  if(length(unit)>1) unit <- unit[1]
  # if (is.T(x)) {
  #   unit <- "degrees*C"
  # }
 
  # Rotation:
  if (is.null(lonR)) lonR <- round(mean(lon,na.rm=TRUE),4)  # longitudinal rotation
  if (is.null(latR)) latR <- round(mean(lat,na.rm=TRUE),4)  # Latitudinal rotation
  if (verbose) print(paste('lonR=',lonR,'latR=',latR))
  # axiR: rotation of Earth's axis

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
  #ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  #theta <- pi*geoborders$x[ok]/180; phi <- pi*geoborders$y[ok]/180

  gx <- geoborders$x
  gy <- geoborders$y
  ok <- is.finite(gx) & is.finite(gy)
  ## KMP 2024-08-09: Apply xlim and ylim to the geoborders
  if (!is.null(xlim)) ok <- ok & gx>=min(xlim) & gx<=max(xlim)
  if (!is.null(ylim)) ok <- ok & gy>=min(ylim) & gy<=max(ylim)
  ## REB 2023-03-10
  xrng <- range(lon)
  yrng <- range(lat)
  ok <- ok & gx>=min(xrng) & gx<=max(xrng)
  ok <- ok & gy>=min(yrng) & gy<=max(yrng)
  
  ## KMP 2025-02-06: Calculating spherical coordinates with a separate function
  xyz_geoborders <- cartesian2sphere(gx[ok], gy[ok], lonR=lonR, latR=latR)
  x <- xyz_geoborders$X
  y <- xyz_geoborders$Y
  z <- xyz_geoborders$Z
  #theta <- pi*gx[ok]/180
  #phi <- pi*gy[ok]/180
  #
  #x <- sin(theta)*cos(phi)
  #y <- cos(theta)*cos(phi)
  #z <- sin(phi)

  # Calculate contour lines if requested...  
  if (verbose) print('contour lines...')
  lonxy <- rep(lon,length(lat))
  latxy <- sort(rep(lat,length(lon)))
  map <- c(map)
  
  # Remove grid boxes with missing data:
  ok <- is.finite(map)
  #print(paste(sum(ok)," valid grid point"))
  lonxy <- lonxy[ok]; latxy <- latxy[ok]; map <- map[ok]

  # Define the grid box boundaries:
  dlon <- min(abs(diff(lon)),na.rm=TRUE); dlat <- min(abs(diff(lat)),na.rm=TRUE)
  Lon <- rbind(lonxy - 0.5*dlon,lonxy + 0.5*dlon,
               lonxy + 0.5*dlon,lonxy - 0.5*dlon)
  Lat <- rbind(latxy - 0.5*dlat,latxy - 0.5*dlat,
               latxy + 0.5*dlat,latxy + 0.5*dlat)
  
  ## KMP 2025-02-06: Calculating spherical coordinates with a separate function
  xyz_gridboxes <- cartesian2sphere(Lon, Lat, lonR=lonR, latR=latR)
  X <- xyz_gridboxes$X
  Y <- xyz_gridboxes$Y
  Z <- xyz_gridboxes$Z
  
  #Theta <- pi*Lon/180; Phi <- pi*Lat/180
  #
  # Transform -> (X,Y,Z):
  #X <- sin(Theta)*cos(Phi)
  #Y <- cos(Theta)*cos(Phi)
  #Z <- sin(Phi)
  
  ## AM commented
  ## KMP 2019-10-11: uncommented color palette definition
  ## because otherwise map2sphere doesn't work
  ## AM 2021-06-02 There is still sth wrong with color palette definition but I cannot figure out what is the problem
  # Define colour palette:
  if (verbose) print('colbar...')
  if (is.null(colbar)) {
    colbar$show <- FALSE
  } else if (is.null(colbar$show)) {
    colbar$show <- TRUE
  }
  if (is.null(colbar$rev)) colbar$rev <- FALSE
  if (is.null(colbar$breaks)) {
    colbar$breaks <- pretty(c(map),n=31)
    colbar$n <- length(colbar$breaks)-1
  } else {
    colbar$n <- length(colbar$breaks) -1
  }
  if (is.null(colbar$srt)) colbar$srt <- 0.45
  ## REB 2023-02-09: flexibility to deal with old habit of confusing colbar$col with colbar$pal
  if ( is.null(colbar$pal) & is.character(colbar$col) ) {colbar$pal <- colbar$col; colbar$col <- NULL}
  nc <- length(colbar$col)
  if (is.null(colbar$col)) {
    colbar <- colbar.ini(map,colbar=colbar)
    col <- colscal(n=colbar$n)
  }
  
  ## AM 2021-06-03: Moved this before index
  ## REB 2015-11-25: Set all values outside the colour scales to the colour scale extremes
  if (verbose) print('Clip the value range to extremes of colour scale')
  toohigh <- map>max(colbar$breaks)
  if (sum(toohigh)>0) map[toohigh] <- max(colbar$breaks)
  toolow <- map<min(colbar$breaks)
  if (sum(toolow)>0) map[toolow] <- min(colbar$breaks)
  if (verbose) print(paste(sum(toohigh),'set to highest colour and',sum(toolow),'to lowest'))
    
  index <- findInterval(map,colbar$breaks,all.inside=TRUE)
  ## where all.inside does to the indices what the clipping does to the values.
 
  ## KMP 2025-02-06: Spherical coordinates are calculated with a separate function (see line 98)
  if (verbose) {print('map2sphere: set colours'); print(colbar)}
  # Rotate coastlines:
  #if (verbose) {print('Rotation:');print(dim(rotM(x=0,y=0,z=lonR))); print(dim(rbind(x,y,z)))}
  #a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
  #a <- rotM(x=latR,y=0,z=0) %*% a
  #x <- a[1,]; y <- a[2,]; z <- a[3,]

  ## KMP 2025-02-06: Spherical coordinates are calculated with a separate function (see line 98)
  # Grid coordinates:
  #d <- dim(X)
  #print(d)
  
  # Rotate data grid:  
  #A <- rotM(x=0,y=0,z=lonR) %*% rbind(c(X),c(Y),c(Z))
  #A <- rotM(x=latR,y=0,z=0) %*% A
  #X <- A[1,]; Y <- A[2,]; Z <- A[3,]
  #dim(X) <- d; dim(Y) <- d; dim(Z) <- d
  
  # Plot the results:
  if (new) {
    if(verbose) print("Create new graphic device")
    dev.new()
    par(mgp=c(2,0.5,0), mar=c(4,1,2,1))
    add <- FALSE
  }
  if(!is.null(fig)) par(fig=fig,new=(add & dev.cur()>1))
  par(bty="n")
  
  ## KMP 2024-08-09: Adding space under the map to fit the color bar
  zlim <- range(Z, na.rm=TRUE)
  dz <- 0.3*diff(zlim)
  if(colbar$show) zlim <- zlim + c(-1,0)*dz else zlim <- zlim + c(-0.1,0)*dz

  ## KMP 2025-02-04: trying to fix issues with map being cut off
  xlim <- range(X) + 0.05*diff(range(X))*c(-1, 1)

  ## REB 2023-03-10
  ## REB: Why 'zlim' and not 'ylim'? KMP 2024-08-08: Because the rotated coordinates are called x and z (plot(x, z)). 
  ## KMP 2025-02-04: xlim was previously commented out for some reason but this caused trouble with the map, 
  ##  not showing the whole area, so I'm adding it back in and hopefully this will do the trick
  if(!add) plot(x, z, xaxt="n", yaxt="n", pch=".", col="grey90", xlim=xlim, 
                ylim=zlim, xlab="", ylab="", main=main)
  
  # plot the grid boxes, but only the gridboxes facing the view point:
  if (verbose) print('Visible grid boxes')
  Visible <- colMeans(Y) > 0
  X <- X[,Visible]; Y <- Y[,Visible]; Z <- Z[,Visible]
  index <- index[Visible]
  
  ## REB 2020-01-26
  if (style=='night') {
    if (verbose) print('Add night-day shading')
    ## Add shadow effect to colours
    brightness <- cos(pi*Lon[1,Visible]/180 - pi*lonR/180)
    
  } else brightness <- rep(1,length(index))
  alpha <- rep(1,length(index))
  
  if (verbose) {
    print(c(length(X),length(Z),length(index),
            length(brightness),length(alpha)))
    print(dim(X))
  }
  ## ADD COLOR FILL IF REQUESTED
  if (sum(is.element(tolower(type),'fill'))>0) {
    apply(rbind(X,Z,index,brightness,alpha),2,gridbox,colbar$col)
  } else {
    colbar$show <- FALSE
  }
  # Plot the coast lines  
  if (verbose) print('plot the coast lines')
  visible <- y > 0
  if (verbose) print(paste(sum(visible),'coast-line points'))
  points(x[visible],z[visible],pch=".")
  if (verbose) {print(summary(x)); print(summary(y))}
  
  ## ADD CONTOURS IF REQUESTED
  if (sum(is.element(tolower(type),'contour'))>0) {
    map_contour <- map
    dim(map_contour) <- c(length(lon), length(lat))
    attr(map_contour, "longitude") <- lon
    attr(map_contour, "latitude") <- lat
    if(is.null(breaks_contour)) {
      breaks_contour <- colbar$breaks
      if(length(breaks_contour)>5) breaks_contour <- 
          breaks_contour[seq(1, length(breaks_contour), floor(length(breaks_contour)/5))]
    }
    spherical_contour(map_contour, lonR = lonR, latR = latR, col=col_contour,
          nx = nx, ny = ny, breaks = breaks_contour, add = TRUE, verbose = verbose)
  }
  
  ## KMP 2024-08-19: This line looks bad when the data is subset in space.
  ##  Perhaps it can be adapted to follow the spatial subset, but for now 
  ##  I'm just excluding it when xlim or ylim is specified
  if(is.null(xlim) & is.null(ylim)) lines(cos(pi/180*1:360),sin(pi/180*1:360))
  if (colbar$show) {
    if (verbose) print('plot colourbar')
    #if (is.null(breaks))
    #  breaks <- round( nc*(seq(min(map),max(map),length=nc)- min(map) )/
    #                  ( max(map) - min(map) ) )
    #bar <- cbind(breaks,breaks)
    #image(breaks,c(1,2),bar,col=col)
    #
    #par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
    #    xaxt = "n",fig=par0$fig,mar=par0$mar,new=TRUE)
    #
    # Adopt from map.station
    par(xaxt="s",yaxt="s",cex.lab=cex.lab,cex.axis=cex.axis)
    if (fancy) {
      if (verbose) print("fancy colbar")
      col.bar(min(xlim,na.rm=TRUE), 
              min(zlim,na.rm=TRUE), 
              max(xlim,na.rm=TRUE), 
              min(zlim,na.rm=TRUE) + dz/2,
              colbar$breaks,horiz=TRUE,pch=21,v=colbar$v,h=colbar$h,
              col=colbar$col,cex=2,cex.lab=colbar$cex.lab,
              cex.axis=colbar$cex.axis,srt=colbar$srt,
              type=colbar$type,verbose=FALSE,vl=1,border=FALSE)
    } else {
      if (verbose) print("regular colbar")
      image.plot(col=colbar$col,breaks=colbar$breaks,
                 lab.breaks=colbar$breaks,
                 horizontal = TRUE, legend.only = TRUE,
                 zlim = range(colbar$breaks),
                 pal = colbar$pal, nlevel=length(colbar$breaks)-1, 
                 legend.width = 1,rev = colbar$rev,
                 axis.args = list(cex.axis = colbar$cex.axis),border=FALSE,
                 verbose=verbose, ...)
        ##image.plot(lab.breaks=colbar$breaks,horizontal = TRUE,
        ##             legend.only = T, zlim = range(colbar$breaks),
        ##             col = colbar$col, legend.width = 1,
        ##             axis.args = list(cex.axis = 0.8), border = FALSE)
    }
    if(!is.null(fig)) par(fig=fig)
  }

  ## plot(range(x,na.rm=TRUE),range(z,na.rm=TRUE),type="n",
  ##     xlab="",ylab="",add=FALSE)
  if(is.null(text.sub) & !is.null(param)) {
    param <- as.character(param)
    unit <- as.character(unit)
    text.sub <- param
    if(!is.null(unit) & (unit!='')) text.sub <- paste(param,'~(',unit,')') else
      if(!is.null(unit)) text.sub <- param
  }
  text.sub <- eval(parse(text=paste('expression(',text.sub,')')))
  if(length(text.sub)>0) {
    if(grepl("top", pos)) text(min(X)+diff(range(X))*0.0, max(Z)+diff(range(Z))*0.01,
           text.sub, cex=cex.sub, pos=4) else 
             text(min(X)+diff(range(X))*0.0, min(Z)-diff(range(Z))*0.05,
                  text.sub, cex=cex.sub, pos=4)
  }
  #result <- data.frame(x=colMeans(Y),y=colMeans(Z),z=c(map))
  
  attr(Z,'longitude') <- X
  attr(Z,'latitude') <- Y
  attr(Z,'variable') <- esd::varid(x)
  attr(Z,'unit') <- esd::unit(x)
  attr(Z,'colbar') <- colbar
  attr(Z,'time') <- range(index(x))
  invisible(Z)
}

