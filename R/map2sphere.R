# Documentation in map.R
#' @export
map2sphere <- function(x,it=NULL,is=NULL,new=TRUE,style="plain",
                       colbar= list(col='t2m',rev=FALSE,n=10,
                           breaks=NULL,type="p",cex=2, cex.axis=0.9,
                           cex.lab = 0.9, h=0.6, v=1,pos=0.05),
                       lonR=NULL,latR=NULL,axiR=0,
                       type=c("fill","contour"),                      
                       gridlines=TRUE,fancy=FALSE,
                       main=NULL,xlim=NULL,ylim=NULL,verbose=FALSE,...) {
  if (verbose) print(paste('map2sphere:',lonR,latR,axiR))
  if (verbose) {print(lon(x)); print(lat(x))}
  if (!is.null(it) | !is.null(is)) x <- subset(x,it=it,is=is,verbose=verbose)
  
  ## KMP 10-11-2015: apply xlim and ylim
  is <- NULL
  if (!is.null(xlim)) is$lon <- xlim
  if (!is.null(ylim)) is$lat <- ylim
  x <- subset(x,is=is)
  
  # Data to be plotted:
  lon <- lon(x)  # attr(x,'longitude')
  lat <- lat(x)  # attr(x,'latitude')
  if (verbose) {print(summary(lon)); print(summary(lat))}
  # To deal with grid-conventions going from north-to-south or east-to-west:
  if (inherits(x,'field')) map <- map(x,plot=FALSE) else map <- x
  srtx <- order(lon); lon <- lon[srtx]
  srty <- order(lat); lat <- lat[srty]
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
  if (is.T(x)) {
    unit <- "degrees*C"
  }
 
  # Rotation:
  if (is.null(lonR)) lonR <- mean(lon,na.rm=TRUE)  # logitudinal rotation
  if (is.null(latR)) latR <- mean(lat,na.rm=TRUE)  # Latitudinal rotation
  if (verbose) print(paste('lonR=',lonR,'latR=',latR))
  # axiR: rotation of Earth's axis

  # coastline data:
  geoborders <- NULL # KMP 2019-10-11: create dummy to avoid warning during CHECK
  data("geoborders",envir=environment())
  #ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  #theta <- pi*geoborders$x[ok]/180; phi <- pi*geoborders$y[ok]/180

  ## KMP 10-11-2015: apply xlim and ylim
  gx <- geoborders$x
  gy <- geoborders$y
  ok <- is.finite(gx) & is.finite(gy)
  if (!is.null(xlim)) ok <- ok & gx>=min(xlim) & gx<=max(xlim)
  if (!is.null(ylim)) ok <- ok & gy>=min(ylim) & gy<=max(ylim)
  theta <- pi*gx[ok]/180
  phi <- pi*gy[ok]/180

  x <- sin(theta)*cos(phi)
  y <- cos(theta)*cos(phi)
  z <- sin(phi)

# Calculate contour lines if requested...  
#contourLines
  lonxy <- rep(lon,length(lat))
  latxy <- sort(rep(lat,length(lon)))
  map <- c(map)

# Remove grid boxes with missign data:
  ok <- is.finite(map)
  #print(paste(sum(ok)," valid grid point"))
  lonxy <- lonxy[ok]; latxy <- latxy[ok]; map <- map[ok]

# Define the grid box boundaries:
  dlon <- min(abs(diff(lon)),na.rm=TRUE); dlat <- min(abs(diff(lat)),na.rm=TRUE)
  Lon <- rbind(lonxy - 0.5*dlon,lonxy + 0.5*dlon,
               lonxy + 0.5*dlon,lonxy - 0.5*dlon)
  Lat <- rbind(latxy - 0.5*dlat,latxy - 0.5*dlat,
               latxy + 0.5*dlat,latxy + 0.5*dlat)
  Theta <- pi*Lon/180; Phi <- pi*Lat/180

# Transform -> (X,Y,Z):
  X <- sin(Theta)*cos(Phi)
  Y <- cos(Theta)*cos(Phi)
  Z <- sin(Phi)
  #print(c( min(x),max(x)))

  ## AM commented
  ## KMP 2019-10-11: uncommented colour palette defintion
  ## because otherwise map2sphere doesn't work
  # Define colour palette:
  if (is.null(colbar$rev)) colbar$rev <- FALSE
  if (is.null(colbar$breaks)) {
    colbar$breaks <- pretty(c(map),n=31)
    colbar$n <- length(colbar$breaks)-1
  } else {
    colbar$n <- length(colbar$breaks) -1
  }
  nc <- length(colbar$col)
  if (is.null(colbar$col)) {
    colbar <- colbar.ini(map,colbar=colbar)
    col <- colscal(n=colbar$n-1) 
  } else if (nc==1) {
    col <- colscal(pal=colbar$col,n=colbar$n-1)
  }
  if (colbar$rev) col <- rev(col)
  ## AM commented
  ## OL 2018-01-26: The following line assumes that breaks are regularly spaced
  #index <- round( nc*( map - min(colbar$breaks) )/
  #                  ( max(colbar$breaks) - min(colbar$breaks) ) )
  ## The findInterval implementation can use irregularly spaced breaks.
  ## (If a point has the same value as a break it will be assigned to the bin above it.)
  index <- findInterval(map,colbar$breaks,all.inside=TRUE)
  ## where all.inside does to the indices what the clipping does to the values.
  
  ## REB 2015-11-25: Set all values outside the colour scales to the colour scale extremes
  if (verbose) print('Clip the value range to extremes of colour scale')
  toohigh <- map>max(colbar$breaks)
  if (sum(toohigh)>0) map[toohigh] <- max(colbar$breaks)
  toolow <- map<min(colbar$breaks)
  if (sum(toolow)>0) map[toolow] <- min(colbar$breaks)
  if (verbose) print(paste(sum(toohigh),'set to highest colour and',sum(toolow),'to lowest'))
  
  ## KMP 2015-09-29: extra colors if higher/lower values occur  # REB: this gives strange colour bars
  #crgb <- col2rgb(colbar$col)
  #if(any(map>max(colbar$breaks))) {
  #  cmax <- crgb[,nc] + (crgb[,nc]-crgb[,nc-1])*0.5
  #  crgb <- cbind(crgb,cmax)
  #  index[index>nc] <- nc+1
  #  colbar$breaks <- c(colbar$breaks,max(map))
  #}
  #if(any(map<min(colbar$breaks))) {
  #  cmin <- crgb[,1] + (crgb[,1]-crgb[,2])*0.5
  #  crgb <- cbind(cmin,crgb)
  #  index[index>nc] <- nc+1
  #  colbar$breaks <- c(min(map),colbar$breaks)
  #}
  #crgb[crgb>255] <- 255; crgb[crgb<0] <- 0
  #colbar$col <- rgb(t(crgb),maxColorValue=255)
  #colbar$n <- length(colbar$col)-1
  #if (min(colbar$breaks)<min(map)) index[map<min(colbar$breaks)] <- 1
  #if (max(colbar$breaks)>max(map)) index[map>max(colbar$breaks)] <- nc
  if (verbose) {print('map2sphere: set colours'); print(colbar)}
  
  # Rotate coastlines:
  a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
  a <- rotM(x=latR,y=0,z=0) %*% a
  x <- a[1,]; y <- a[2,]; z <- a[3,]

  # Grid coordinates:
  d <- dim(X)
  #print(d)

  # Rotate data grid:  
  A <- rotM(x=0,y=0,z=lonR) %*% rbind(c(X),c(Y),c(Z))
  A <- rotM(x=latR,y=0,z=0) %*% A
  X <- A[1,]; Y <- A[2,]; Z <- A[3,]
  dim(X) <- d; dim(Y) <- d; dim(Z) <- d
  #print(dim(rbind(X,Z)))
  
  # Plot the results:
  if (new) {
    dev.new()
    par(fig=c(0,1,0.1,1), mgp=c(2,0.5,0), mar=c(4,1,2,1))
  }
  par(bty="n") ## ,xaxt="n",yaxt="n")
  plot(x,z,xaxt="n",yaxt="n",pch=".",col="grey90",xlab="",ylab="",main=main)
  
  # plot the grid boxes, but only the gridboxes facing the view point:
  Visible <- colMeans(Y) > 0
  X <- X[,Visible]; Y <- Y[,Visible]; Z <- Z[,Visible]
  index <- index[Visible]
  
  ## REB 2020-01-26
  if (style=='night') {
    ## Add shadow effect to collours
    brightness <- cos(Theta[1,Visible] + lonR)
    
  } else brightness <- rep(1,length(index))
  alpha <- rep(1,length(index))
  
  if (verbose) {print(c(length(X),length(Z),length(index),length(brightness),length(alpha)))
    print(dim(X))}
  
  apply(rbind(X,Z,index,brightness,alpha),2,gridbox,col)
  # c(W,E,S,N, colour)
  # xleft, ybottom, xright, ytop
  # Plot the coast lines  
  visible <- y > 0
  points(x[visible],z[visible],pch=".")
  #plot(x[visible],y[visible],type="l",xlab="",ylab="")
  lines(cos(pi/180*1:360),sin(pi/180*1:360))
  # Add contour lines?
  # Add grid ?
  # Colourbar:
  if (!is.null(colbar)) {
    if (verbose) print('plot colourbar')
    #par0 <- par()
    #par(fig = c(0.3, 0.7, 0.05, 0.10),cex=0.8,
    #    new = TRUE, mar=c(1,0,0,0), xaxt = "s",yaxt = "n",bty = "n")
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
    if(is.null(colbar$show)) colbar$show <- TRUE
    par(xaxt="s",yaxt="s",cex.lab=0.7,cex.axis=0.9)
    if(colbar$show) {
    if (fancy & !is.null(colbar)) {
      if (verbose) print("fancy colbar")
      col.bar(colbar$breaks,horiz=TRUE,pch=21,v=colbar$v,h=colbar$h,#v=1,h=1,
              col=colbar$col,cex=2,cex.lab=colbar$cex.lab,
              cex.axis=colbar$cex.axis,
              type=type,verbose=FALSE,vl=1,border=FALSE)
    } else if (!is.null(colbar)) {
      if (verbose) print("regular colbar")
      image.plot(breaks=colbar$breaks,
                 lab.breaks=colbar$breaks,
                 horizontal = TRUE, legend.only = TRUE,
                 zlim = range(colbar$breaks),
                 pal = colbar$col, nlevel=length(colbar$breaks)-1, 
                 legend.width = 1,rev = colbar$rev,
                 axis.args = list(cex.axis = colbar$cex.axis),border=FALSE,
                 verbose=verbose, ...)
                 #xaxp=c(range(colbar$breaks),colbar$n)),
                 #border = FALSE,...)
      ##image.plot(lab.breaks=colbar$breaks,horizontal = TRUE,
      ##             legend.only = T, zlim = range(colbar$breaks),
      ##             col = colbar$col, legend.width = 1,
      ##             axis.args = list(cex.axis = 0.8), border = FALSE)
    }
    }  
  }

  ## plot(range(x,na.rm=TRUE),range(z,na.rm=TRUE),type="n",
  ##     xlab="",ylab="",add=FALSE)
  txt <- param
  if (!is.null(param)) {
    param <- as.character(param); unit <- as.character(unit)
    if(!is.null(unit) & (unit!='')) txt <- paste(param,'~(',unit,')') else
      if(!is.null(unit)) txt <- param
    text(min(x),max(z),eval(parse(text=paste('expression(',txt,')'))),
         cex=1.5,pos=4)
  }
  #result <- data.frame(x=colMeans(Y),y=colMeans(Z),z=c(map))
  result <- NULL # For now...
  invisible(result)
}

