regrid.irregweights <- function(xo,yo,xn,yn,verbose=FALSE) {
# Compute the weights for irreguolar grids (xo,yo) - the station class
  
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
  if (nxo != nyo) stop(paste('regrid.irregweights requires',nx,'=',ny))
  
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
    ## II has different length to the other elements.
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
  
  # The length of lon/lat is the same as the input data itself, 

  beta <- cbind(rep(Wx[1,],ny),
                rep(Wx[2,],ny),
                rep(Wy[1,],nx),
                rep(Wy[2,],nx))
#print("denom")
  denom <- rep(Wx[8,],ny)*rep(Wy[8,],nx)
  beta <- beta/denom
#print("indx")
  indx <- cbind(rep(Wx[3,],ny),
                rep(Wx[4,],ny),
                rep(Wy[3,],nx),
                rep(Wy[4,],nx))
  beta[!is.finite(indx)] <- 0
  indx[!is.finite(indx)] <- 1
  attr(beta,'index') <- indx
  attr(beta,'Wx') <- Wx
  attr(beta,'Wy') <- Wy
  t2 <- Sys.time()
  if (verbose) {print("Computation time:"); print(t2 - t1)}
  invisible(beta)
}




regrid.station <- function(x,is,approach="station",clever=FALSE,verbose=FALSE) {

  stopifnot(inherits(x,'station'))
  print('regrid.station')  
  if (approach=="pca2station") {
    y <- regrid.pca2station(x,is)
    return(y)
  }
  
  #print("regrid.field ")
  #x <- sp2np(x)
  
  # The coordinates which to grid: lon.new & lat.new
  if ( (is.data.frame(is)) | (is.list(is)) ) {lon.new <- is[[1]]; lat.new <- is[[2]]} else
  if ( (inherits(is,'station')) | (inherits(is,'field')) | (inherits(is,'eof')) | (inherits(is,'pca')) ) {
    lon.new <- attr(is,'longitude'); lat.new <- attr(is,'latitude')
  }
  
  #greenwich <- attr(x,'greenwich')
  #if (verbose) print(paste('greenwich=',greenwich))
  # REB 13.05.2014
  #if ( (min(lon.new) < 0) & (max(lon.new) <= 180) ) x <- g2dl(x,greenwich=FALSE) else
  #if ( (min(lon.new) > 0) & (max(lon.new) <= 360) ) x <- g2dl(x,greenwich=FALSE) else
  #   stop(paste('Bad longitude range: ',min(lon.new),'-',max(lon.new))) 
  
  if (verbose) {print("New coordinates"); print(lon.new); print(lat.new)}

  #print("regrid.field before subset:");print(lon(x)); print(lat(x));print("---")
  #x <- subset(x,is=list(lon=c(floor(min(lon.new))-2,ceiling(max(lon.new))+2),
  #                      lat=c(floor(min(lat.new))-2,ceiling(max(lat.new))+2)))
  #if (verbose) {print("regrid.field after subset:");print(lon(x));print(lat(x));print("---")}
  if (verbose) {print("Longitude range"); print(range(lon.new)); print(range(lon(x)))}
  if (verbose) {print("Latitude range"); print(range(lat.new)); print(range(lat(x)))}
  
  mlon <- MATCH(attr(x,'longitude'), lon.new)
  mlat <- MATCH(attr(x,'latitude'),  lat.new)
  #print(mlon); print(mlat)
  # If thhose corordinates are the same as the original data, then
  # return with the original data:
  if ( sum(is.na(c(mlon,mlat))) ==0 ) {
    if (verbose) print(summary(mlon,mlat))
     if ( max( c(diff(mlon),diff(mlat)) ) == 1) return(x)
  }
  
  
  ns <- length(attr(x,'longitude'))
  if (verbose) print(paste("regrid.field: from",
                           ns,"locations to", length(lon.new),"x",length(lat.new),'regular grid'))  
  
# Apply the grid-based subset - only use the nearest longitudes/latitudes
#  for the regridding. 

# Clever solution (not default): fix edge problems...  
 
  lon.old <- attr(x,'longitude')
  lat.old <- attr(x,'latitude')
 
  #print("...");print(lonj2w)
  #print(lon.old); print(lat.old)
  if (verbose) print(paste("subsample the original grid and use",length(lon.old),
                      "longitudes and",length(lat.old),"latitudes"))
  #print("union:"); print(union(attr(x,'longitude'),lon.nearest))
  #print("...")

  # What does this do?
  # x <- subset(x,is=list(i=lon.old,j=lat.old))
  
#  beta <- regrid.weights(lon.old,lat.old,lon.new,lat.new,verbose=verbose)
  beta <- regrid.irregweights(lon.old,lat.old,lon.new,lat.new,verbose=verbose)
  if (verbose) {print("Weight matrix");  print(dim(beta))}
                                          
  d <- dim(x)
  X <- x; 
  D <- c(length(lon.new),length(lat.new))
  y <- matrix(rep(NA,D[1]*D[2]*d[1]),D[1]*D[2],d[1])

  #print(dim(cbind(beta,attr(beta,'index'))))
  
  if (verbose) pb <- txtProgressBar(style=3)
  #str(X); str(beta)
  
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
  print("set attributes:")
  y <- zoo(t(y),order.by=index(x))
  class(y) <- c('field',class(x)[2:3]) 
  #mostattributes(y) <- attributes(x)
  #nattr <- softattr(x,ignore=c('longitude','latitude','dimensions'))
  #for (i in 1:length(nattr))
  #  attr(y,nattr[i]) <- attr(x,nattr[i])
  y <- attrcp(x,y)
  attr(y,'longitude') <- lon.new
  attr(y,'latitude') <- lat.new
  attr(y,'altitude') <- lat.new+NA
  attr(y,'location') <- NULL
  attr(y,'dimensions') <- c(D,d[1])
  if (verbose) print(paste("New dimensions:",attr(y,'dimensions')[1],"x",
             attr(y,'dimensions')[2],"x",attr(y,'dimensions')[3] ))
  #attr(y,'history') <- c(attr(y,'history'),'regrid.field')
  #attr(y,'date-stamp') <- date()
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  #print("---")

  #if (attr(y,'greenwich') != greenwich) y <- g2dl(y,greenwich)
  invisible(y)
}


regrid.pca <- function(x,is,verbose=FALSE) {
  stopifnot(inherits(x,'pca'))
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

  beta <- regrid.irregweights(lon.old,lat.old,lon.new,lat.new)
                                          
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



regrid.pca2station <- function(x,is) {
  stopifnot(inherits(x,'station'),inherits(is,'list'))
  print("regrid.pca.station - estimate EOFs")

 #greenwich <- attr(x,'greenwich')
  if ( greenwich & (min(lon.new < 0)) ) x <- g2dl(x,greenwich=FALSE)
  pca0 <- PCA(x)
  pca1 <- regrid(pca0,is)
  print("Reconstruct field from EOFs")
  y <- pca2station(pca1)
  if (attr(y,'greenwich') != greenwich) y <- g2dl(y,greenwich)
  #attr(y,'history') <- c(attr(X,'history'),'regrid.eof2field')
  #attr(y,'date-stamp') <- date()
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  return(y)
}
