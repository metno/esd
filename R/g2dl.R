#' Transform longitudes between a system of 0-360E and 180W-180E 
#' 
#' Transform longitudes between a system starting at the Greenwich line (\code{greenwich}=TRUE), going from 0 to 360 degrees, 
#' and one starting at the date line (\code{greenwich}=FALSE) going from -180 to 180 degrees.
#' 
#' \code{g2dl} is an S3 method and will redirect to a fitting function depending on the input. 
#' The output of \code{g2dl} is of the same class and format as the input. The attribute 'greenwich' (\code{attr(x,'greenwich')})
#' holds information about the longitude system of an object.
#' 
#' @aliases g2dl.default g2dl.stationmeta g2dl.field g2dl.eof g2dl.corfield g2dl.events g2dl.trajectory
#' 
#' @param x the input object
#' @param greenwich a boolean; if TRUE longitudes are transformed to a system starting at the Greenwich line (0-360E); if FALSE longitudes are transformed to a system starting at the date line (180W-180E)
#' @param verbose a boolean; if TRUE print information about progress
#' @param lon longitudes
#' @param lat latitudes
#' @param d dimensions
#' @param ... other arguments
#' 
#' @return an object of the same type as the input
#' 
#' @export
g2dl <- function(x,greenwich=TRUE,verbose=FALSE,...) UseMethod("g2dl")

#' @exportS3Method
#' @export g2dl.default
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

#' @exportS3Method
#' @export g2dl.station
g2dl.station <- function(x,greenwich=TRUE,verbose=FALSE,...,lon=NULL) {
  if(verbose) {print("g2dl.default"); str(x)}
  if (is.null(lon)) lon <- attr(x,'longitude')
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
    y <- y[,xsrt]
  }
  y <- attrcp(x,y)
  attr(y,'longitude') <- lon
  class(y) <- class(x)
  return(y)
}


#' @exportS3Method
#' @export g2dl.stationmeta
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

#' @exportS3Method
#' @export g2dl.field
g2dl.field <- function(x,greenwich=TRUE,verbose=FALSE,...) {
  if(verbose) print("g2dl.field")
  attr(x,'longitude') -> lon
  attr(x,'latitude') -> lat
  ## REB 2024-08-26 - more robust code using the actual object size
  #d <- attr(x,'dimensions')
  d <- c(length(index(x)),length(lon),length(lat))
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  
  xsrt <- order(lon)
  xsrt <- xsrt[!duplicated(lon)]
  if (verbose) print(paste('g2dl: length(lon)=',length(lon),length(xsrt)))
  X <- t(coredata(x))
  dim(X) <-  c(d[2],d[3],d[1])
  X <- X[xsrt,,]
  if (verbose) {print(dim(X)); print(d)}
  dim(X) <- c(d[2]*d[3],d[1])
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

#' @exportS3Method
#' @export g2dl.eof
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

#' @exportS3Method
#' @export g2dl.corfield
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

#' @exportS3Method
#' @export g2dl.events
g2dl.events <- function(x,greenwich=TRUE,verbose=FALSE,...) {
  if(verbose) print("g2dl.events")
  lon <- as.numeric(x$lon)
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

#' @exportS3Method
#' @export g2dl.trajectory
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