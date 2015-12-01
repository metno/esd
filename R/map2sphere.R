# This function maps a longitude-latitude grid onto a sphere. This
# function will be embeeded in map as an option.
#" c("lonlat","sphere","NP","SP")
# First transform longitude (theta) and latitude (phi) to cartesian
# coordinates:
#
# (theta,phi) -> (x,y,z) = r
#
# Use a transform to rotate the sphere through the following matrix
# product:
#
# ( x' )    ( a1  b1  c1 )  ( x )
# ( y' ) =  ( a2  b2  c2 )  ( y )
# ( z' )    ( a3  b3  c3 )  ( z )
#
# There can be three different rotation, about each of the axes
# For 90-degree rotatins, e.g. aroud the z-axis:
# Z = 90 (X = Y = 0): x' = y; y' = -x: z'= z
# a1 = cos(Z) ; b1=sin(Z); c1=cos(Z);
# a2 = -sin(Z); b2=cos(Z); c2=cos(Z)
# a3 = cos(Z) ; b3=cos(Z); c3=sin(Z)
# Y = 90 (Z = X = 0): x' = -z; z' = x; y' = y
# a1= cos(Y); b1 = cos(Y); c1=-sin(Y);
# a2= cos(Y); b2 = sin(Y); c2= cos(Y);
# a3= sin(Y); b3 = cos(Y); c3= cos(Y)
# X = 90 (Y = Z = 0): x' = x; z' = y; y' = - z
# a1 = sin(X); b1 = cos(X); c1 = cos(X)
# a2 = cos(X); b2 = cos(X); c2 = -sin(X);
# a3 = cos(X); b3 = sin(X), c3 = cos(X)
#
# a1 =  cosYcosZ; b1 = cosXsinZ; c1= -cosXsinY
# a2 = -cosXsinZ; b2 = cosXcosZ; c2= -sinXcosY
# a3 =  cosXsinY; b3 = sinXcosZ; c3=  cosXcosY
#
# Rasmus E. Benestad
# 14.04.2013


rotM <- function(x=0,y=0,z=0) {
  X <- -pi*x/180; Y <- -pi*y/180; Z <- -pi*z/180
  cosX <- cos(X); sinX <- sin(X); cosY <- cos(Y)
  sinY <- sin(Y); cosZ <- cos(Z); sinZ <- sin(Z)
  rotM <- rbind( c( cosY*cosZ, cosY*sinZ,-sinY*cosZ ),
                 c(-cosX*sinZ, cosX*cosZ,-sinX*cosZ ),
                 c( sinY*cosX, sinX*cosY, cosX*cosY ) )
  return(rotM)
}

gridbox <- function(x,col,density = NULL, angle = 45) {
#  if (length(x) != 5) {print(length(x)); stop("gridbox")}
# W,E,S,N
# xleft, ybottom, xright, ytop
  i <- round(x[9])
  if(i==0) {border <- rgb(0,0,0,0.1)
  } else border <- col[i]
#  rect(x[1], x[3], x[2], x[4],col=cols[i],border=cols[i])
  polygon(c(x[1:4],x[1]),c(x[5:8],x[5]),
          col=col[i],border=border)
}

#angleofview <- function(r,P) {
#  angle <- acos( (P%*%r)/( sqrt(P%*%P)*sqrt(r%*%r) ) )
#  return(angle)
#}

map2sphere <- function(x,it=NULL,is=NULL,new=TRUE,
                       colbar= list(pal='t2m',rev=FALSE,n=10,
                           breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                       lonR=NULL,latR=NULL,axiR=0,
                       type=c("fill","contour"),                      
                       gridlines=TRUE,fancy=FALSE,
                       xlim=NULL,ylim=NULL,verbose=FALSE,...) {

  
  if (verbose) print(paste('map2sphere:',lonR,latR,axiR))
  if (verbose) {print(lon(x)); print(lat(x))}
  ## If only a few items are provided in colbar - then set the rest to the default
  if (!is.null(colbar)) {
    if (verbose) print('sort out the colours')
    if (is.null(colbar$pal)) colbar$pal <- 't2m'
    if (is.null(colbar$rev)) colbar$rev <- FALSE
    if (is.null(colbar$n)) colbar$n <- 10
    if (is.null(colbar$breaks)) {
      colbar$breaks <- pretty(x,n=colbar$n)
    } else if (length(colbar$breaks)==2)
      colbar$breaks <- seq(colbar$breaks[1],colbar$breaks[2],length=colbar$n)
    colbar$n <- length(colbar$breaks)-1
    if (is.null(colbar$type)) colbar$type <- 'p'
    if (is.null(colbar$cex)) colbar$cex <- 2
    if (is.null(colbar$h)) colbar$h <- 0.6
    if (is.null(colbar$v)) colbar$v <- 1
    if (is.null(colbar$pos)) colbar$pos <- 0.05
    if (verbose) print(colbar)
    colbar$col <- colscal(n=colbar$n,col=colbar$pal,
                          rev=colbar$rev,verbose=verbose)
    if (verbose) print(paste("length(col) =",length(colbar$col)))
  }
  ## if (!is.null(colbar$col)) col <- colbar$col else col <- NULL
  ## if (!is.null(colbar$breaks)) breaks <- colbar$breaks else breaks <- NULL
  if (!is.null(it) | !is.null(is))
    x <- subset(x,it=it,is=is,verbose=verbose)

  ## KMP 10-11-2015: apply xlim and ylim
  is <- NULL
  if (!is.null(xlim)) is$lon <- xlim
  if (!is.null(ylim)) is$lat <- ylim
  x <- subset(x,is=is)

  # Data to be plotted:
  lon <- lon(x)  # attr(x,'longitude')
  lat <- lat(x)  # attr(x,'latitude')
  # To deal with grid-conventions going from north-to-south or east-to-west:
  srtx <- order(lon); lon <- lon[srtx]
  srty <- order(lat); lat <- lat[srty]
  map <- x[srtx,srty]
  param <- attr(x,'variable')
  unit <- attr(x,'unit')
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
 
  # Rotatio:
  if (is.null(lonR)) lonR <- mean(lon)  # logitudinal rotation
  if (is.null(latR)) latR <- mean(lat)  # Latitudinal rotation
  if (verbose) print(paste('lonR=',lonR,'latR=',latR))
  # axiR: rotation of Earth's axis

  # coastline data:
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
  map<- c(map)

# Remove grid boxes with missign data:
  ok <- is.finite(map)
  #print(paste(sum(ok)," valid grid point"))
  lonxy <- lonxy[ok]; latxy <- latxy[ok]; map <- map[ok]

# Define the grid box boundaries:
  dlon <- min(abs(diff(lon))); dlat <- min(abs(diff(lat)))
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

# Define colour palette:
  ##if (is.null(breaks)) {
  ##  n <- 30
  ##  breaks <- pretty(c(map),n=n)
  ##} else n <- length(breaks)
  ## AM commented
  ##if (is.null(col)) col <- colscal(n=n) else
  ##if (length(col)==1) {
  ##    palette <- col
  ##    col <- colscal(col=palette,n=n)
  ## }
  nc <- length(colbar$col)
  ## AM commented
  index <- round( nc*( map - min(colbar$breaks) )/
                    ( max(colbar$breaks) - min(colbar$breaks) ) )

  ## REB 2015-11-25: Set all values outside the colour scales to the colour scale extremes
  print('Clip the value range to extremes of colour scale')
  toohigh <- map>max(colbar$breaks)
  if (sum(toohigh)>0) map[toohigh] <- max(colbar$breaks)
  toolow <- map<min(colbar$breaks)
  if (sum(toolow)>0) map[toolow] <- min(colbar$breaks)
  print(paste(sum(toohigh),'set to highest colour and',sum(toolow),'to lowest'))
  
  ## KMP 2015-09-29: extra colors if higher/lower values occur  # REB: this gives strange colour bars
  crgb <- col2rgb(colbar$col)
  if(any(map>max(colbar$breaks))) {
    cmax <- crgb[,nc] + (crgb[,nc]-crgb[,nc-1])*0.5
    crgb <- cbind(crgb,cmax)
    index[index>nc] <- nc+1
    colbar$breaks <- c(colbar$breaks,max(map))
  }
  if(any(map<min(colbar$breaks))) {
    cmin <- crgb[,1] + (crgb[,1]-crgb[,2])*0.5
    crgb <- cbind(cmin,crgb)
    index[index>nc] <- nc+1
    colbar$breaks <- c(min(map),colbar$breaks)
  }
  crgb[crgb>255] <- 255; crgb[crgb<0] <- 0
  colbar$col <- rgb(t(crgb),maxColorValue=255)
  colbar$n <- length(colbar$col)-1
  #if (min(colbar$breaks)<min(map)) index[map<min(colbar$breaks)] <- 1
  #if (max(colbar$breaks)>max(map)) index[map>max(colbar$breaks)] <- nc
  if (verbose) print('set colours')
  
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
  if (new) dev.new()
  par(bty="n") ## ,xaxt="n",yaxt="n")
  plot(x,z,xaxt="n",yaxt="n",pch=".",col="grey90",xlab="",ylab="")
  
# plot the grid boxes, but only the gridboxes facing the view point:
  Visible <- colMeans(Y) > 0
  X <- X[,Visible]; Y <- Y[,Visible]; Z <- Z[,Visible]
  index <- index[Visible]
  apply(rbind(X,Z,index),2,gridbox,colbar$col)
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
#    #print(breaks)
#    par0 <- par()
#    par(fig = c(0.3, 0.7, 0.05, 0.10),cex=0.8,
#        new = TRUE, mar=c(1,0,0,0), xaxt = "s",yaxt = "n",bty = "n")
#  #print("colourbar")
#    if (is.null(breaks))
#      breaks <- round( nc*(seq(min(map),max(map),length=nc)- min(map) )/
#                      ( max(map) - min(map) ) )
#    bar <- cbind(breaks,breaks)
#    image(breaks,c(1,2),bar,col=col)
#
#    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
#        xaxt = "n",fig=par0$fig,mar=par0$mar,new=TRUE)
#
    # Adopt from map.station
    par(xaxt="s",yaxt="s",cex.lab=0.7,cex.axis=0.7)
    if (fancy & !is.null(colbar)) {
      if (verbose) print("fancy colbar")
      col.bar(colbar$breaks,horiz=TRUE,pch=21,v=1,h=1,
              col=colbar$col, cex=2,cex.lab=colbar$cex.lab,
              type=type,verbose=FALSE,vl=1,border=FALSE)
    } else if (!is.null(colbar)) {
      if (verbose) print("regular colbar")
      image.plot(breaks=colbar$breaks,
                 lab.breaks=colbar$breaks,
                 horizontal = TRUE,legend.only = T,
                 zlim = range(colbar$breaks),
                 col = colbar$col, legend.width = 1,
                 axis.args = list(cex.axis = 0.8),border=FALSE,...)
                 #xaxp=c(range(colbar$breaks),colbar$n)),
                 #border = FALSE,...)
      ##image.plot(lab.breaks=colbar$breaks,horizontal = TRUE,
      ##             legend.only = T, zlim = range(colbar$breaks),
      ##             col = colbar$col, legend.width = 1,
      ##             axis.args = list(cex.axis = 0.8), border = FALSE)
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

#map2sphere(x)

vec <- function(x,y,it=10,a=1,r=1,ix=NULL,iy=NULL,new=TRUE,nx=150,ny=80,
                projection='lonlat',lonR=NULL,latR=NULL,axiR=0,...) {
  x <- subset(x,it=it); y <- subset(y,it=it)
  d <- attr(x,'dimensions')
  #print(d); print(dim(x))
  if (is.null(ix)) ix <- pretty(lon(x),n=nx)
  if (is.null(iy)) iy <- pretty(lat(x),n=ny)
  #print(c(d[2],d[1]))
  X <- coredata(x); Y <- coredata(y)
  dim(X) <- c(d[1],d[2])
  dim(Y) <- c(d[1],d[2])
  #X <- t(X); Y <- t(Y)
  x0 <- rep(ix,length(iy))
  y0 <- sort(rep(iy,length(ix)))
  ij <- is.element(ix,lon(x))
  ji <- is.element(iy,lat(x))
  #print(ix); print(lon(x)); print(sum(ij))
  #print(iy); print(lat(x)); print(sum(ji))
  dim(x0) <- c(length(ij),length(ji)); dim(y0) <- dim(x0)
  x0 <- x0[ij,ji]
  y0 <- y0[ij,ji]
  ii <- is.element(lon(x),ix)
  jj <- is.element(lat(x),iy)
  x1 <- a*X[ii,jj]; y1 <- a*Y[ii,jj]
  #print(dim(x1)); print(c(length(x0),sum(ii),sum(jj)))
  x1 <- x0 + x1; y1 <- y0 + y1
  if (projection=='sphere') {
    # Rotate data grid:
    # The coordinate of the arrow start:
    
    theta <- pi*x0/180; phi <- pi*y0/180
    x <- r*c(sin(theta)*cos(phi))
    y <- r*c(cos(theta)*cos(phi))
    z <- r*c(sin(phi))
    a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
    x0 <- a[1,]; y0 <- a[3,]
    invisible <- a[2,] < 0
    x0[invisible] <- NA; y0[invisible] <- NA
    #The coordinate of the arrow end:
    theta <- pi*x1/180; phi <- pi*y1/180
    x <- r*c(sin(theta)*cos(phi))
    y <- r*c(cos(theta)*cos(phi))
    z <- r*c(sin(phi))
    a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
    x1 <- a[1,]; y1 <- a[3,]
    invisible <- a[2,] < 0
    x1[invisible] <- NA; y1[invisible] <- NA
  }    

  if (new) {
    dev.new()
    plot(range(x0,x1),range(y0,y1),xlab='',ylab='')
    data(geoborders)
    lines(geoborders$x,geoborders$y)
  }
  arrows(x0, y0, x1, y1,...)
}
