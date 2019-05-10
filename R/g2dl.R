g2dl <- function(x,greenwich=TRUE,verbose=FALSE,...)
  UseMethod("g2dl")

g2dl.default <- function(x,greenwich=TRUE,verbose=FALSE,...,lon=NULL,lat=NULL,d=NULL) {
  if(verbose) {print("g2dl.default"); str(x)}
  if (is.null(lon)) lon <- attr(x,'longitude')
  if (is.null(lat)) lat <- attr(x,'latitude')
  if (is.null(d)) d <- attr(x,'dimensions') 
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  y <- x
  if(length(lon)>1) {
    xsrt <- order(lon)
    xsrt <- xsrt[!xsrt %in% which(duplicated(lon))]
    dim(y) <- d
    if (length(dim(y))==3) {
      ## For a field object
      y <- y[xsrt,,] 
      dim(y) <- c(length(lon)*length(lat),d[3])
    } else if (length(dim(y))==2) {
      ## For a matrix
      y <- y[xsrt,] 
    } else stop(paste('Problem in gfdl.default - x has more than 3 or 3 dimensions',dim(x),collapse=' '))
    lon <- lon[xsrt]
  }
  y <- attrcp(x,y)
  attr(y,'longitude') <- lon
  class(y) <- class(x)
  return(y)
}

g2dl.stationmeta <- function(x,greenwich=TRUE,verbose=FALSE,...) {
  if(verbose) print("g2dl.stationmeta")
  lon <- x$lon                          
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  y <- x
  y$lon <- lon
  attr(y,'greenwich') <- as.logical(greenwich)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

g2dl.field <- function(x,greenwich=TRUE,verbose=FALSE,...) {
  if(verbose) print("g2dl.field")
  attr(x,'longitude') -> lon
  attr(x,'latitude') -> lat
  d <- attr(x,'dimensions')
  
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  
  xsrt <- order(lon)
  xsrt <- xsrt[!duplicated(lon)]
  X <- t(coredata(x))
  dim(X) <- d
  X <- X[xsrt,,]
  dim(X) <- c(length(lon)*d[2],d[3])
  y <- zoo(t(X),index(x))
  lon <- sort(lon)
  
  y <- attrcp(x,y,ignore='names')
  #nattr <- softattr(x,ignore=c('greenwich','longitude'))
  #for (i in 1:length(nattr))
  #  attr(y,nattr[i]) <- attr(x,nattr[i])
  attr(y,'dimensions') <- attr(x,'dimensions')
  attr(y,'longitude') <- lon
  attr(y,'greenwich') <- as.logical(greenwich)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

g2dl.eof <- function(x,greenwich=TRUE,verbose=FALSE,...) {
  if(verbose) print("g2dl.eof")
  attr(x,'longitude') -> lon
  attr(x,'latitude') -> lat
  d <- attr(x,'dimensions')
  X <- attr(x,'pattern')
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  xsrt <- order(lon)
  X <- X[xsrt,,]
  lon <- sort(lon)
  X -> attr(x,'pattern')
  attr(x,'greenwich') <- greenwich
  attr(x,'longitude') <- lon
  attr(x,'latitude') <- lat
  return(x)
}

g2dl.corfield <- function(x,greenwich=TRUE,verbose=FALSE,...) {
  if(verbose) print("g2dl.corfield")
  attr(x,'longitude') -> lon
  attr(x,'latitude') -> lat
  d <- attr(x,'dimensions')
  #print(d)
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  
  y <- x
  xsrt <- order(lon)
  dim(y) <- d
  y <- y[xsrt,]
  y <- c(y)
  lon <- sort(lon)
  
  y <- attrcp(x,y,ignore='names')
  #nattr <- softattr(x,ignore=c('greenwich','longitude'))
  #for (i in 1:length(nattr))
  #  attr(y,nattr[i]) <- attr(x,nattr[i])
  attr(y,'dimensions') <- attr(x,'dimensions')
  attr(y,'longitude') <- lon
  attr(y,'greenwich') <- as.logical(greenwich)
  class(y) <- class(x)
  invisible(y)
}

g2dl.events <- function(x,greenwich=TRUE,verbose=FALSE,...) {
  if(verbose) print("g2dl.events")
  lon <- x$lon                          
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  y <- x
  y$lon <- lon
  attr(y,'greenwich') <- as.logical(greenwich)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

g2dl.trajectory <- function(x,greenwich=TRUE,verbose=FALSE,...) {
  if(verbose) print("g2dl.trajectory")
  lon <- x[,colnames(x)=="lon"]                         
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  y <- x
  y[,colnames(y)=="lon"] <- lon
  attr(y,'greenwich') <- as.logical(greenwich)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}