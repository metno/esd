#' Multivariate Empirical Orthogonal Functions (EOFs).
#' 
#' Computes EOFs (a type of principal component analysis) for an \code{mvcomb} object which is
#' the combination of two fields of different variables and output of the \code{mvcombine} function.
#'
#' @seealso EOF
#'
#' @param X An \code{mvcomb} object or \code{field} object
#' @param Y A field object. If two \code{field} objects X and Y are provided, \code{mvcombine} is applied to combine the two. 
#' @param it A list or data.frame providing time index, e.g. a range of years like c(1979,2010), a season ('djf'), or a month ('dec' or 'december').
#' @param is A list or data.frame providing space index, e.g. a list of longitude and latitude range like list(lon=c(0,60), lat=c(35,60)).
#' @param n Number of EOFs.
#' @param lon set longitude range - see \code{\link{t2m.NCEP}}
#' @param lat set latitude range.
#' @param verbose If TRUE print diagnostics
#' @param \dots Additional arguments
#' 
#' @exportS3Method
#' @export
EOF.mvcomb <- function(X, Y=NULL, it=NULL, is=NULL, n=20, lon=NULL, lat=NULL,
                       geoweighting=FALSE, verbose=FALSE, ...) {
  if(verbose) print("EOF.mvcomb")
  attr(X,'dimnames') <- NULL
  stopifnot(!missing(X), is.matrix(X),
            inherits(X,c("field","zoo")))
  if(!inherits(X, "mvcomb") & !is.null(Y)) {
    X <- mvcombine(X, Y, it=it, is=is, verbose=verbose)
  }
  
  SF <- function(x) {sum(is.finite(x))}
  
  x <- X
  dates <- index(x)
  if(verbose) print(dates)
  
  cls <- class(x)
  id <- attr(X, "id")
  d <- attr(X, "dimensions")
  #d <- dim(X)
  
  Y <- t(coredata(x))
  Y[!is.finite(Y)] <- NA
  
  if (verbose) print('Apply geographical weighting')
  if(geoweighting) {
    for(i in seq(length(unique(id)))) {
      Y.i <- Y[id==unique(id)[i],]
      stdv <- sd(c(Y.i),na.rm=TRUE)  # Account for mixed fields with different magnitudes...
      Wght <- matrix(nrow=d[[i]][1], ncol=d[[i]][2])
      for (j in 1:d[[i]][1]) Wght[j,] <- sqrt(abs(cos(pi*attr(X,'latitude')[[i]]/180)))
      dim(Wght) <- c(d[[i]][1]*d[[i]][2])
      if(verbose) {print(length(Wght)); print(dim(Y.i)); print(d[[i]][3])}
      for (k in 1:d[[i]][3]) Y.i[,k] <- (Wght/stdv)*Y.i[,k]
      Y[id==unique(id)[i],] <- Y.i 
    }
  }
  
  ## KMP 2020-02-18: Switching off this part for now because other information (lon,lat,id)
  ## has to be updated if grid points or time slices are exluded.
  ## Exclude the missing values 'NA' and grid points with sd == 0 for all times:
  sd0 <- apply(as.matrix(Y),2,sd,na.rm=TRUE)
  nf <- apply(as.matrix(Y),2,SF)
  if (verbose) print(paste('Exclude the missing values/zero-sd:',
                           sum(sd0>0.0),sum(nf > 0)))
  if(any((sd0<0) | (nf < 0))) browser()
  #y <- Y[,(sd0>0.0) & (nf > 0)]
  y <- Y
  
  # Exclude the time slices with missing values:
  skip <- apply(as.matrix(y),1,SF)
  npts <- dim(y)[2]
  if(any((skip!=npts))) browser()
  #y <- as.matrix(y)[skip==npts,]
  
  npca <- min(dim(y)) 
  ny <- min(c(dim(y),20))
  
  if(verbose) print("Apply SVD to joined field")
  SVD <- try(svd(y, nu=min(c(ny,npca)), nv=min(c(ny,npca))))
  if (inherits(SVD,"try-error")) {
    if (verbose) print("svd(x) failed. try svd(t(x))")
    SVD <- try(svd(t(X1), nu=min(c(ny,npca)),nv=min(c(ny,npca))))
    temp <- SVD$u
    SVD$u <- SVD$v
    SVD$v <- temp
  }
  if (inherits(SVD,"try-error") & verbose) print("both svd(x) and svd(t(x) failed.")
  
  autocor <- 0
  if (verbose) print(paste("Find max autocorr in the grid boxes."))
  for (i in 1:npts) {
    vec <- as.vector(y[,i])
    i.bad <- is.na(vec)
    if (sum(i.bad) == 0) {
      ar1 <- acf(vec[],plot=FALSE)
      autocor <- max(c(autocor,ar1$acf[2,1,1]),na.rm=TRUE)
    }
  }
  n <- min(n,length(SVD$d))
  
  eof <- zoo(SVD$v[,1:n], order.by=dates)
  invert <- apply(SVD$u[,1:n], 2, mean) < 0
  eof[,invert] <- -eof[,invert]
  eof <- attrcp(X, eof)
  names(eof) <- paste("X.",1:n,sep="")
  
  patterns <- list()
  for(i in seq_along(unique(id))) {
    pattern.i <- matrix(rep(NA,d[[i]][1]*d[[i]][2]*n), d[[i]][1]*d[[i]][2], n)
    i.id <- id==unique(id)[i]
    pattern.i[,] <- SVD$u[i.id,1:n]
    
    # Some data points may have been excluded due to missing values.
    # Need to insert the results for valid data onto the original grid.
    #pattern <- matrix(rep(NA,d[1]*d[2]*n),d[1]*d[2],n)
    #pattern[skip == npts,] <- SVD$u[,1:n]
    
    # Make all the EOF vectors havine the same sense rather than
    # being random:
    if (verbose) print(paste("Invert EOF",(1:length(invert))[invert],collapse=' '))
    pattern.i[,invert] <- -pattern.i[,invert]
    dim(pattern.i) <- c(d[[i]][1],d[[i]][2],n)
    patterns[[paste0("pattern.",unique(id)[i])]] <- pattern.i
  }
  
  attr(eof,'pattern') <- patterns
  attr(eof,'dimensions') <- d
  attr(eof,'max.autocor') <- autocor
  attr(eof,'eigenvalues') <- SVD$d[1:n]
  attr(eof,'sum.eigenv') <- sum(SVD$d)
  attr(eof,'tot.var') <- sum(SVD$d^2)
  attr(eof,'history') <- history.stamp(X)
  attr(eof,'aspect') <- 'anomaly'
  attr(eof,'dimnames') <- NULL
  class(eof) <- c("eof",cls)
  return(eof)
}
