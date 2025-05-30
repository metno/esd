# Uses bi-linear interpolation to regrid one field onto the next
# First the routine computes a set of weight, then performs a matrix multiplication
# to map the original data onto the new grid.
# The weights are based on the distance between points, taking longitude & latitude
# and use distAB() to estimate the geographical distance in km.
#
# X(i,j) is a i-j matrix containing the data on a grid with i logitudes and j latitudes
# X(i,j) -> Y(k,l)
# Y = beta X
# beta(i*j,k*l)
#  ( Y(1,1) )   (beta(1,1), beta(2,1), beta(3,1), ... ) ( X(1,1) )
#  ( Y(1,2) ) = (beta(1,2), beta(2,2), beta(3,2), ... ) ( X(1,2) )
#  ( .....  )   (beta(1,3), beta(2,3), beta(3,3), ... ) ( X(1,3) )
#
# Most of the elements in Beta are zero, and only 4 weights are needed to compute
# the interpolated value within a grid with 4 values at the corners.
#
# The weights are estimated first, and are then applied to each time slice through
# a sparse matrix multiplication (hence keeping track of indices).
#
# The weights are derived according to:
# http://en.wikipedia.org/wiki/Bilinear_interpolation
#
# R.E. Benestad, 
# rasmus.benestad@met.no
#
#------------------------------------------------------------------------

#' @export
sparseMproduct <- function(beta,x) {
  # b contains the weights and i the indexes
  b <- beta[1:4]; i <- beta[5:8]
  ok <- (i > 0)
  y <- sum(b[ok]*x[i][ok])  
  y
}


#' Regrid
#' 
#' Fast transform data from one longitude-latitude grid to another through
#' bi-linear interpolation. The regridding is done by first calculating a set
#' of weights. This is a "QUICK & DIRTY" way of getting approximate results.
#' More sophisticated methods exist (e.g. Kriging - LatticeKrig).
#' 
#' Let X(i,j) be a i-j matrix containing the data on a grid with i logitudes
#' and j latitudes. We want to transform this to a different grid with k
#' longitudes and l latitudes:
#' 
#' X(i,j) -> Y(k,l)
#' 
#' First the routine computes a set of weight, then performs a matrix
#' multiplication to map the original data onto the new grid.  The weights are
#' based on the distance between points, taking longitude & latitude and use
#' distAB() to estimate the geographical distance in km.
#' 
#' The matrix operation is: Y = beta X
#' 
#' beta is a matrix with dimensions (i*j,k*l)
#' 
#' ( Y(1,1) ) (beta(1,1), beta(2,1), beta(3,1), ... ) ( X(1,1) ) ( Y(1,2) ) =
#' (beta(1,2), beta(2,2), beta(3,2), ... ) ( X(1,2) ) ( .....  ) (beta(1,3),
#' beta(2,3), beta(3,3), ... ) ( X(1,3) )
#' 
#' Most of the elements in Beta are zero!
#' 
#' 
#' @aliases regrid regrid.default regrid.field regrid.station
#' regrid.matrix regrid.eof regrid.pca sparseMproduct
#' 
#' @param xo Old x-coordinates (longitudes)
#' @param yo Old y-coordinates (latitudes)
#' @param xn New x-coordinates (longitudes)
#' @param yn New y-coordinates (latitudes)
#' @param beta The matrix of interpolation weights
#' @param x a field object.
#' @param is A list holding the coordinates lon and lat, a field object, an eof
#' object, or a station object - for the latter three, the field x is
#' interpolated to the longitude/latitude held by is.
#' @param approach 'station' or 'pca2station'. If 'pca2station', the stations
#' are turned into PCAs before regridding and then converted back to station
#' objects.
#' @param verbose If TRUE, print out diagnostics
#' @return A field object
#' @author R.E. Benestad and A. Mezghani
#' @keywords utilities
#' @examples
#' 
#' # Use regrid to interpolate to station location:
#' t2m <- t2m.DNMI()
#' data(Oslo)
#' z.oslo <- regrid(t2m,is=Oslo)
#' plot(Oslo)
#' lines(z.oslo)
#' 
#' # Regrid t2m onto the grid of the gcm
#' gcm <- t2m.NorESM.M()
#' Z <- regrid(t2m,is=gcm)
#' map(Z)
#' 
#' # Example using regrid on a matrix object:
#' t2m.mean <- as.pattern(t2m,FUN='mean')
#' z <- regrid(t2m.mean,is=list(lon=seq(min(lon(t2m)),max(lon(t2m)),by=0.5),
#'                              lat=seq(min(lat(t2m)),max(lat(t2m),by=0.5))))
#' image(lon(z),lat(z),z)
#' # Add land borders on top
#' data(geoborders)
#' lines(geoborders)
#' 
#' \dontrun{
#' ## Regrid station data using weights defined by the distance of the 4
#' ## nearest stations: quick and dirty method
#' if (!file.exists("stationsVALUE_exp1a.rda")) {
#'   download.file("http://files.figshare.com/2085591/value_predictands4exp1a.R",
#'                 "value_predictands4exp1a.R")
#'   source("value_predictands4exp1a.R")
#' }
#'    
#' load("stationsVALUE_exp1a.rda")
#' TX <- regrid(Tx,is=list(lon=seq(-8,30,by=1),lat=seq(40,60,by=0.5)))
#' map(TX)
#' }
#' 
#' @seealso aggregate.grid
#' @export regrid
regrid <- function(x,is=NULL,...)
  UseMethod("regrid")

regridweights <- function(xo,yo,xn,yn,verbose=FALSE) {
# Compute weights for regular grids: (xo,yo) - the field class
  
  lindist <- function(x,X) {
    # Linear distance
    # This function works for non-overlapping coordinates:
    I <- 1:length(X)
    lt <- X < x; gt <- X > x
    if ( (sum(lt)>0) & (sum(gt)>0) ) {
      i1 <-  max(I[lt]); i2 <- min(I[gt])
      x1 <- X[i1]; x2 <- X[i2]
      dx <- x2 - x1
      w1 <- x2 - x
      w2 <- x - x1
      output <- c(w1,w2,i1,i2,x,x1,x2,dx)
    } else output <- c(NA,NA,1,1,x,NA,NA,diff(X)[1])
    names(output) <- c("w1","w2","i1","i2","x","x1","x2","dx")
    return(output)
  }
  
  t1 <- Sys.time()
  nx <- length(xn); ny <- length(yn)
  nxo <- length(xo)
  dim(xn) <- c(nx,1)
  nyo <- length(yo)
  dim(yn) <- c(ny,1)
  
  # Split the search for weights in two:
  # 1) Find the overlapping coordinates:
  if (verbose) print(paste("Set up weights: nx=",nx,"ny=",ny))
  Wx <- matrix(rep(0,nx*8),8,nx)         # REB 29.01.2014
  Wy <- matrix(rep(0,ny*8),8,ny)
  rownames(Wx) <- c("w1","w2","i1","i2","x","x1","x2","dx")
  rownames(Wy) <- c("w1","w2","i1","i2","x","x1","x2","dx")
  olx <- is.element(xn,xo); dx <- min(diff(xo))
  oly <- is.element(yn,yo); dy <- min(diff(yo))
#  print(c(sum(olx),sum(oly),length(olx),length(oly),dx,dy))
  if (sum(olx)>0) {
     II <- (1:nxo)[is.element(xo,xn)]
     Wxol <- rbind(rep(1,sum(olx))*dx, rep(0,sum(olx)),  # Weights
                   II,    II+1,                          # indices
                   xn[olx],    xn[olx],   xn[olx]+dx,
                   rep(dx,sum(olx)))
     #str(Wxol); str(Wx)
     # Ensure that the indexing doesn't go beyond matrix size
     ifix <- Wxol[4,] > nxo
     if (sum(ifix)>0) Wxol[4,ifix] <- nxo
     Wx[,olx] <- Wxol[,]
  }
  if (sum(oly)>0) {
     JJ <-  (1:nyo)[is.element(yo,yn)]
     Wyol <- rbind(rep(1,sum(oly))*dy,  rep(0,sum(oly)),   # Weights
                   JJ,     JJ+1,                           # indices
                   yn[oly],     yn[oly],   yn[oly]+dy,
                   rep(dy,sum(oly)))
     #str(Wyol); str(Wy)
      # Ensure that the indexing doesn't go beyond matrix size
     jfix <- Wyol[4,] > nyo
     if (sum(jfix)>0) Wyol[4,jfix] <- nyo
     Wy[,oly] <- Wyol[,]
   }
   #print(dim(Wx)); print(dim(Wy))
  
  # 2) Find the non-overlapping coordinates: no 
#  Wx <- apply(xn,1,lindist,xo) 
#  Wy <- apply(yn,1,lindist,yo)
  if (verbose) print("Non-overlapping coordinates")
  xno <- xn[!olx]; dim(xno) <- c(sum(!olx),1)
  yno <- yn[!oly]; dim(yno) <- c(sum(!oly),1)
  if (sum(!olx)>0) {
    Wxno <- apply(xno,1,lindist,xo) # REB 29.01.2014
    Wx[,!olx] <- Wxno[,]
  }
  if (sum(!oly)>0) {
    Wyno <- apply(yno,1,lindist,yo)  # REB 29.01.2014              
    #str(Wxno); str(Wx)
    Wy[,!oly] <- Wyno[,]
  }
  #print(dim(Wx)); print(dim(Wy))
  if (verbose) print("test: colsums: sum of weights should be 1")
  if (verbose) {
    if (dim(Wx)[2]>1) print(colSums(Wx[1:2,])/Wx[8,]) else
                            print(Wx)
    if (dim(Wx)[2]>1) print(colSums(Wy[1:2,])/Wy[8,]) else
                            print(Wy)
  }
  if (verbose) {print(table(Wy[1,],Wy[2,])); print(table(Wy[3,],Wy[4,]))}  
  #print("beta:")
  
  srty <- order(rep(Wy[5,],nx))
  beta <- cbind(rep(Wx[1,],ny)*rep(Wy[1,],nx)[srty],
                rep(Wx[2,],ny)*rep(Wy[1,],nx)[srty],
                rep(Wx[1,],ny)*rep(Wy[2,],nx)[srty],
                rep(Wx[2,],ny)*rep(Wy[2,],nx)[srty])
  #print("denom")
  denom <- rep(Wx[8,],ny)*rep(Wy[8,],nx)[srty]
  beta <- beta/denom

  #print("indx")
  indx <- cbind(rep(Wx[3,],ny)+(rep(Wy[3,],nx)[srty]-1)*nxo,
                rep(Wx[4,],ny)+(rep(Wy[3,],nx)[srty]-1)*nxo,
                rep(Wx[3,],ny)+(rep(Wy[4,],nx)[srty]-1)*nxo,
                rep(Wx[4,],ny)+(rep(Wy[4,],nx)[srty]-1)*nxo)
  beta[!is.finite(indx)] <- 0
  indx[!is.finite(indx)] <- 1
  attr(beta,'index') <- indx
  attr(beta,'Wx') <- Wx
  attr(beta,'Wy') <- Wy
  if (verbose) {
    print('longitude from indx:')
    print(rbind(xo[Wx[3,]],xo[Wx[4,]]))
    print('latitude from indx:')
    print(rbind(yo[Wy[3,]],yo[Wy[4,]]))
    #par(mfcol=c(2,1))
    #image(xo[Wx[3,]]%o%yo[Wy[3,]])
    #image(xo[Wx[4,]]%o%yo[Wy[4,]])
  }
  t2 <- Sys.time()
  if (verbose) {print("Computation time:"); print(t2 - t1)}
  invisible(beta)
}


regridtemporal <- function(x,it,verbose=FALSE) {
  if (verbose) print('regridtemporal')
  stopifnot(is.field(x))
  if (verbose) {print(index(x)); print(it)}
  if (verbose) print(paste(sum(is.finite(x)),'data-boxes with valid data and',
                           sum(!is.finite(x)),'with no data'))
  ## Only aply interpolation on grid-boxes with (some) valid data.
  ok <- apply(coredata(x),2,FUN='nv') >= 2
  if (verbose) print(paste('Interpolating',sum(ok),'grid-boxes'))
  y <- zoo(coredata(x)[,ok],order.by=index(x))
  if (verbose) print(summary(c(coredata(y))))
  #zc <- apply(y,2,function(x) approx(index(x),coredata(x),it)$y)
  z <- matrix(rep(NA,length(ok)*length(it)),length(it),length(ok))
  for (i in 1:dim(y)[2])
    z[,(1:length(ok))[ok][i]] <- approx(index(y),coredata(y[,i]),it)$y
  if (verbose) print(summary(c(z)))
  if (verbose) {print(dim(z[,ok])); print(dim(x)); print(dim(z))}
  #z[,ok] <- zc
  z <- zoo(z,order.by=it)
  z <- attrcp(x,z)
  class(z) <- class(x)
  if (verbose) print('finished regridtemporal')
  return(z)
}

#' @exportS3Method
#' @export
regrid.field <- function(x,is=NULL,...,it=NULL,verbose=FALSE,approach="field",clever=FALSE) {

  if (verbose) print("regrid.field ")
  stopifnot(inherits(x,'field'))
 
  if (approach=="eof2field") {
    y <- regrid.eof2field(x,is)
    return(y)
  }
  x <- sp2np(x)
  
  ## If it is provided, also regrid in time
  if (!is.null(it)) x <- regridtemporal(x,it,verbose=verbose)
  if (is.null(is)) return(x)

  ## case wether lon or lat is given in is i.e. regrid on these values along the other dimension
  if (length(is)==1) {
      nm <- names(is)
      if (length(nm)==0) print("Please rewrite 'is' as in is=list(lon=...) or is=list(lat=...)")
      if (grepl("lon",nm))
          is <- list(lon=rep(is[[1]],length(lon(x))),lat=lat(x))
      else if (grepl("lat",nm))
          is <- list(lon=lon(x),lat=is[[2]])            
  }
  ## The coordinates which to grid: lon.new & lat.new
  if ( (is.data.frame(is)) | (is.list(is)) ) {lon.new <- is[[1]]; lat.new <- is[[2]]} else
  if ( (inherits(is,'station')) | (inherits(is,'field')) | (inherits(is,'eof'))| (inherits(is,'pca')) ) {
    if (verbose) print(paste('Use the coordinates from a',class(is)[1],'object'))
    lon.new <- sort(lon(is)); lat.new <- sort(lat((is)))
    if (verbose) print(paste(sum(duplicated(lon.new)),'duplicated longitudes'))
    lon.new <- lon.new[!duplicated(lon.new)]
    if (verbose) print(paste(sum(duplicated(lat.new)),'duplicated latitudes'))
    lat.new <- lat.new[!duplicated(lat.new)]
  }
  
  greenwich <- attr(x,'greenwich')
  if (verbose) print(paste('greenwich=',greenwich))
  # REB 13.05.2014
  if ( (min(lon.new) <= 0) & (max(lon.new) <= 180) ) x <- g2dl(x,greenwich=FALSE) else
  if ( (min(lon.new) >= 0) & (max(lon.new) <= 360) ) x <- g2dl(x,greenwich=FALSE) else
     stop(paste('Bad longitude range: ',min(lon.new),'-',max(lon.new))) 
  
  if (verbose) {print("New coordinates"); print(lon.new); print(lat.new)}
  
  #print("regrid.field before subset:");print(lon(x)); print(lat(x));print("---")
  ## AM 20-04-2015 the value "2" has been changed to "10"
  x <- subset(x,is=list(lon=c(floor(min(lon.new))-10,ceiling(max(lon.new))+10),
                        lat=c(floor(min(lat.new))-10,ceiling(max(lat.new))+10)))
  if (verbose) {print("regrid.field after subset:");print(lon(x));print(lat(x));print("---")}
  if (verbose) {print("Longitude range"); print(range(lon.new)); print(range(lon(x)))}
  if (verbose) {print("Latitude range"); print(range(lat.new)); print(range(lat(x)))}
  
  mlon <- MATCH(attr(x,'longitude'), lon.new)
  mlat <- MATCH(attr(x,'latitude'),  lat.new)
  #print(mlon); print(mlat)
  # If those coordinates are the same as the original data, then
  # return with the original data:
  if ( sum(is.na(c(mlon,mlat))) ==0 ) {
    if (verbose) print(summary(mlon,mlat))
    if ( max( c(diff(mlon),diff(mlat)) ) == 1 & 
         identical(lon.new, attr(x, "longitude")) & 
         identical(lat.new, attr(x, "latitude")) ) return(x)
  }
  
  if (verbose) print(paste("regrid.field: from",
                           length(attr(x,'longitude')),"x",
                           length(attr(x,'latitude')),"to",
                           length(lon.new),"x",length(lat.new)))  
  
# Apply the grid-based subset - only use the nearest longitudes/latitudes
#  for the regridding. 

  nx <- length(attr(x,'longitude')); ny <- length(attr(x,'latitude'))

# Clever solution (not default): fix edge problems...  
  if (clever) {
    if (nx < length(attr(x,'longitude'))) {
      lonj2w <- rep(NA,nx); lonj2e <- lonj2w
      latj2n <- rep(NA,ny); latj2s <- latj2n
      for (i in 1:nx) {
        lonj2w[i] <- max(lon.new[lon.new <= attr(x,'longitude')[i]])
        lonj2e[i] <- min(lon.new[lon.new >= attr(x,'longitude')[i]])
      }
      lon.old <- union(lonj2w,lonj2e)
    } else lon.old <- attr(x,'longitude')
    if (ny < length(attr(x,'latitude'))) {
      for (j in 1:ny) {
        latj2s[j] <- max(lat.new[lat.new <= attr(x,'latitude')[j]])
        latj2n[j] <- min(lat.new[lat.new >= attr(x,'latitude')[j]])
      }
      lat.old <- union(latj2s,latj2n)
    } else lat.old <- attr(x,'latitude')
  } else {
    lon.old <- attr(x,'longitude')
    lat.old <- attr(x,'latitude')
  }
  #print("...");print(lonj2w)
  #print(lon.old); print(lat.old)
  if (verbose) print(paste("subsample the original grid and use",length(lon.old),
                      "longitudes and",length(lat.old),"latitudes"))
  #print("union:"); print(union(attr(x,'longitude'),lon.nearest))
  #print("...")

  # What does this do?
  # x <- subset(x,is=list(i=lon.old,j=lat.old))
  
#  beta <- regridweights(lon.old,lat.old,lon.new,lat.new,verbose=verbose)
  beta <- regridweights(lon.old,lat.old,lon.new,lat.new,verbose=verbose)
  if (verbose) {print("Weight matrix");  print(dim(beta))}
                                          
  d <- dim(x)
  X <- x; 
  D <- c(length(lon.new),length(lat.new))
  y <- matrix(rep(NA,D[1]*D[2]*d[1]),D[1]*D[2],d[1])

  #print(dim(cbind(beta,attr(beta,'index'))))
  
  if (verbose) pb <- txtProgressBar(style=3)
  for (i in 1:d[1]) {
    #if (verbose) cat(".")
    if (verbose) setTxtProgressBar(pb,i/d[1])  
    #z <- apply(cbind(beta,attr(beta,'index')),1,sparseMproduct,coredata(x[i,]))
    #if (verbose) print(c(i,d[1],length(z),length(y[,i]),NA,dim(x),dim(y)))
    #y[,i] <- z
    M <- as.matrix(X[i,attr(beta,"index")])
    dim(M) <- c(D[1] * D[2],4)
    y[, i] <- rowSums(as.matrix(beta) * M)      
  }
  if (verbose) print("set attributes:")
  y <- zoo(t(y),order.by=index(x))
  if (inherits(is,'station')) {
    if (verbose) print('select the station coordinates which is *not* a regular grid')
    xlons <- lon(is); dim(xlons) <- c(length(xlons),1)
    ix <- unlist(apply(xlons,1,function(x) (1:D[1])[is.element(lon.new,x)][1]))
    xlats <- lat(is); dim(xlats) <- c(length(xlats),1)
    iy <- unlist(apply(xlats,1,function(x) (1:D[2])[is.element(lat.new,x)][1]))
    y <- y[,ix + (iy-1)*D[1]]
    #plot.zoo(y[,1]); lines(as.monthly(is)[,1],col='red')
    lon.new <- lon(is)
    lat.new <- lat(is)
    attr(y,'altitude') <- alt(is)
    attr(y,'location') <- loc(is)
    if (verbose) {
      print(dim(y)); print(dim(is)); print(lon(is)); print(lat(is))
    }
  }

  if ( (is.station(is)) | (is.field(is)) ) {
    if (verbose) print(paste('Set appropriate class:',class(is)[1]))
    class(y) <- class(is)
    class(y)[2] <- class(x)[2]
  }  else class(y) <- class(x)

  #mostattributes(y) <- attributes(x)
  #nattr <- softattr(x,ignore=c('longitude','latitude','dimensions'))
  #for (i in 1:length(nattr))
  #  attr(y,nattr[i]) <- attr(x,nattr[i])
  y <- attrcp(x,y)
  attr(y,'longitude') <- lon.new
  attr(y,'latitude') <- lat.new
  attr(y,'dimensions') <- c(D,d[1])
  if (verbose) print(paste("New dimensions:",attr(y,'dimensions')[1],"x",
             attr(y,'dimensions')[2],"x",attr(y,'dimensions')[3] ))
  #attr(y,'history') <- c(attr(y,'history'),'regrid.field')
  #attr(y,'date-stamp') <- date()
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  #print("---")

  #if (attr(y,'greenwich') != greenwich) y <- g2dl(y,greenwich)
  attr(y,'greenwich') <- attr(x,'greenwich')
  invisible(y)
}


#' @exportS3Method
#' @export
regrid.matrix <- function(x,is=NULL,...,verbose=FALSE) {
# assumes that dimensions of x are [x,y,(t)] and that the coordinates are
# provided as attributes as in field
  
  if (verbose) print(match.call())
  d <- dim(x)
  if (length(d)==2) d <- c(d,1)  # if only 2-dimensional maps, then one dimension of 1 is added
  lon.old <- lon(x)
  lat.old <- lat(x)

  if ( (is.data.frame(is)) | (is.list(is)) )
    {lon.new <- is[[1]]; lat.new <- is[[2]]} else
  if ( (inherits(is,'station')) | (inherits(is,'field')) |
      (inherits(is,'eof')) ) {
    lon.new <- lon(is); lat.new <- lat(is)
  } else if ( is.matrix(is) & !is.null(lon(is)) &
              !is.null(lat(is)) ) {
    lon.new <- lon(is); lat.new <- lat(is)
  }
  
  if (verbose) {str(lon.old);str(lat.old);str(lon.new);str(lat.new)}
  if (is.null(lon.old) | is.null(lat.old) | is.null(lon.new) | is.null(lat.new)) return(NULL)

  beta <- regridweights(lon.old,lat.old,lon.new,lat.new,verbose=verbose)
  if (verbose) {print("Weight matrix");  print(dim(beta)); print(attr(beta,"index"))}
  X <- x;
  dim(X) <- c(d[1]*d[2],d[3])
  D <- c(length(lon.new),length(lat.new))
  y <- matrix(rep(NA,D[1]*D[2]*d[3]),D[1]*D[2],d[3])
  for (i in 1:d[3]) {
    if (verbose) cat(".")
    M <- as.matrix(X[attr(beta,"index"),i])
    dim(M) <- c(D[1] * D[2],4)
    y[, i] <- rowSums(as.matrix(beta) * M)      
  }
  
  if (d[3]>1) dim(y) <- c(D[1],D[2],d[3]) else
              dim(y) <- c(D[1],D[2])
  y <- attrcp(x,y)
  attr(y,'longitude') <- lon.new
  attr(y,'latitude') <- lat.new
  attr(y,'dimensions') <- dim(y)
  if (verbose) print(paste("New dimensions:",attr(y,'dimensions')[1],"x",
             attr(y,'dimensions')[2],"x",attr(y,'dimensions')[3] ))
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}



#' @exportS3Method
#' @export
regrid.eof <- function(x,is=NULL,...,verbose=FALSE) {
  if(verbose) print("regrid.eof")
  stopifnot(inherits(x,'eof'))
  if (is.list(is)) {lon.new <- is[[1]]; lat.new <- is[[2]]} else
  if ( (inherits(is,'station')) | (inherits(is,'field')) | (inherits(is,'eof')) ) {
    lon.new <- attr(is,'longitude'); lat.new <- attr(is,'latitude')
  }
  lon.old <- attr(x,'longitude'); lat.old <- attr(x,'latitude')

  greenwich <- attr(x,'greenwich')
  if ( greenwich & (min(lon.new < 0)) ) x <- g2dl(x,greenwich=FALSE)
  
  if (sum(is.na(c( MATCH(lon.old, lon.new),
                   MATCH(lat.old, lat.new) ))) == 0) return(x)
  
  if (verbose) print(paste("regrid.field: from",length(lon.old),"x",
              length(lat.old),"to",
              length(lon.new),"x",length(lat.new)))  

  beta <- regridweights(lon.old,lat.old,lon.new,lat.new)
                                          
  #print(dim(beta))
  X <- attr(x,'pattern')
  d <- dim(X)
  dim(X) <- c(d[1]*d[2],d[3])
  D <- c(length(lon.new),length(lat.new))
  Y <- matrix(rep(NA,D[1]*D[2]*d[3]),D[1]*D[2],d[3])
  for (i in 1:d[3]) {
    cat(".")
#    z <- apply(cbind(beta,attr(beta,'index')),1,sparseMproduct,coredata(X[i,]))
#   print(c(i,d[1],length(z),length(y[,i]),NA,dim(X),dim(y)))
#    Y[,i] <- z
 #   old lines
 #z <- apply(b, 1, sparseMproduct, X[, i])
 #print(c(i, d[1], length(z), length(Y[, i]), NA, dim(X), dim(Y)))
 #Y[, i] <- z
 # New lines
    M <- as.matrix(X[attr(beta,"index"),i])
    dim(M) <- c(D[1] * D[2],4)
    Y[, i] <- rowSums(as.matrix(beta) * M)      
    #print(c(i, d[1], length(z), length(Y[, i]), NA, dim(X), dim(Y)))
    #Y[, i] <- z
  }
  #print("set attributes:")

  y <- x
  if (attr(y,'greenwich') != greenwich) y <- g2dl(y,greenwich)
  #mostattributes(y) <- attributes(x)
  #nattr <- softattr(x,ignore=c('longitude','latitude','dimensions'))
  #for (i in 1:length(nattr))
  #  attr(y,nattr[i]) <- attr(x,nattr[i])
  #dim(Y) <- c(D[1],D[2],d[3])
  y <- attrcp(x,y)
  Y ->  attr(y,'pattern')
  attr(y,'longitude') <- lon.new
  attr(y,'latitude') <- lat.new
  attr(y,'dimensions') <- c(D,d[3])
  if (verbose) print(paste("New dimensions:",attr(y,'dimensions')))
  #attr(y,'history') <- c(attr(y,'history'),'regrid.field')
  #attr(y,'date-stamp') <- date()
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

#' @export
regrid.eof2field <- function(x,is=NULL,...,verbose=FALSE) {
  if(verbose) print("regrid.eof2field")
  stopifnot(inherits(x,'field'),inherits(is,'list'))

  if ( greenwich & (min(lon.new < 0)) ) x <- g2dl(x,greenwich=FALSE)
  eof0 <- EOF(x)
  eof1 <- regrid(eof0,is)
  print("Reconstruct field from EOFs")
  y <- eof2field(eof1)
  if (attr(y,'greenwich') != greenwich) y <- g2dl(y,greenwich)
  #attr(y,'history') <- c(attr(X,'history'),'regrid.eof2field')
  #attr(y,'date-stamp') <- date()
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  return(y)
}






