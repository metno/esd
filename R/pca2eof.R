## Function for gridding station data Y.

## Check if you need to get the esd-package:
#install.LK <- ("LatticeKrig" %in% rownames(installed.packages()) == FALSE)
install.LK <- ("LatticeKrig" %in% rownames(installed.packages()) == FALSE)

if (install.LK) {
  print('Need to install the LatticeKrig package')
  ## Need online access.
  install.packages('LatticeKrig')
  print('The latest version of LatticeKrig has been installed from CRAN')
}

gridstation <- function(Y,i=1,verbose=FALSE,xlim=NULL,ylim=NULL) {
  if (verbose) print(paste('gridstation'))
  if (!requireNamespace("LatticeKrig", quietly = TRUE)) {
    stop("Package 'LatticeKrig' needed to use 'gridstation'. Please install it.")
  } else {
  ## Instead of importing LatticeKrig, just call the functions that we need from the external package
    #require(LatticeKrig)

    if (is.null(xlim)) xlim <- range(lon(Y))
    if (is.null(ylim)) ylim <- range(lat(Y))
    
    ## Get data on the topography on the 5-minute resolution
    data(etopo5)
    etopo5 <- subset(etopo5, is=list(lon=range(lon(Y))+c(-1,1),
                             lat=range(lat(Y))+c(-1,1)))
    ## Mask the sea: elevations below 1m below sea level is masked.
    etopo5[etopo5<=-1] <- NA
    
    ## Set the grid to be the same as that of etopo5:
    grid <- LatticeKrig::structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')
    
    ## Flag dubplicated stations:
    ok <- !(duplicated(lon(Y)) & duplicated(lat(Y)))
    
    if (verbose) print(paste('Use',sum(ok),'locations from',length(lon(Y))))
    
    obj <- LatticeKrig::LatticeKrig( x=cbind(lon(Y)[ok],lat(Y)[ok]),
                        y=Y[ok,i],Z=alt(Y)[ok])
    
    w <- LatticeKrig::predictSurface(obj, grid.list = grid,Z=etopo5)
    w$z[is.na(etopo5)] <- NA
    
    ## Get rid of packages that have functions of same name:
    #detach("package:LatticeKrig")
    #detach("package:fields")
    #detach("package:spam")
    #detach("package:grid")
    #detach("package:maps")
    
    ## Convert the results from LatticeKrig to esd:
    W <- w$z
    attr(W,'variable') <- varid(Y)[1]
    attr(W,'unit') <- unit(Y)[1]
    attr(W,'longitude') <- w$x
    attr(W,'latitude') <- w$y
    
    ## Make a projection that zooms in on the Barents region
    
    invisible(W)
  }
}

## Script that uses gridding to transform the station map into a regular
## map making use of elevation.
pca2eof <- function(x,verbose=FALSE,xlim=NULL,ylim=NULL) {
  if (verbose) print('pca2eof')
  stopifnot(inherits(x,'pca'))
  y <- x
  if (!is.null(z <- attr(x,'pattern'))) {
    z <- attr(x,'pattern')
    attr(z,'longitude') <- lon(x)
    attr(z,'latitude') <- lat(x)
    attr(z,'altitude') <- alt(x)
  } else if (!is.null(x$pca)) z <- x$pca else
                              stop('Do not know how to handle this object!')
  d <- dim(z)
  Z <- list()
  if (verbose) print('Grid the modes')
  for (i in 1:d[2]) {
    Z[[i]] <- gridstation(z,i,verbose=verbose)
  }
  if (verbose) print('Grid the mean')
  zc <- attr(x,'mean'); dim(zc) <- c(length(zc),1)
  attr(zc,'longitude') <- lon(x)
  attr(zc,'latitude') <- lat(x)
  attr(zc,'altitude') <- alt(x)  
  clim <- gridstation(zc,1,verbose=verbose)
  z <- unlist(Z)
  dim(z) <- c(dim(Z[[1]]),d[2])
  z -> attr(y,'pattern')
  clim  -> attr(y,'mean')
  attr(y,'longitude') <- lon(Z[[1]])
  attr(y,'latitude') <- lat(Z[[1]])
  attr(y,'old_longitude') <- lon(zc)
  attr(y,'old_latitude') <- lat(zc)
  attr(y,'old_altitude') <- alt(zc)
  attr(y,'dimensions') <- c(dim(Z[[1]]),d[2])
  attr(y,'variable') <- varid(x)[1]
  attr(y,'unit') <- unit(x)[1]
  attr(y,'longname') <- attr(x,'longname')[1]
  attr(y,'greenwich') <- TRUE
  class(y) <- c('eof','field',class(x)[-c(1,2)])
  return(y)
}

## A function that converts PCA-based DSensemble objects to EOF-based results (gridded)
as.eof.dsensemble.pca <- function(X,is=NULL,it=NULL,ip=NULL,verbose=FALSE,...) {
  if (verbose) print('as.eof.dsensemble.pca')
  stopifnot(inherits(X,"dsensemble") & inherits(X,"pca"))
  if (inherits(X,"eof")) {
      invisible(X)
  } else {
    eof <- pca2eof(X$pca)
    eof <- subset(eof,ip=ip)
    if (!is.null(is)) eof <- subset(eof,is=is,it=it,verbose=verbose)
    X$eof <- eof 
    class(X) <- c("dsensemble", "eof", "list")
    invisible(X)
  }
}


## Function for convertin station data to field data vie the computation of PCAs
## grididng to EOFs and then transforming the EOFs to field object.
station2field <- function(x,verbose=FALSE) {
    if (verbose) print('station2field')
    stopifnot(inherits(x,'station'))
    x <- pcafill(x,verbose=verbose)
    pca <- PCA(x)
    eof <- pca2eof(pca,verbose=verbose)
    X <- eof2field(eof,verbose=verbose)
    return(X)
}
