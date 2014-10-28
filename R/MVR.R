# Multi-variate regression (MVR) and MVR-based predictions for EOFs, PCA
# and gridded fiels/data matrix (zoo).
#
# R.E. Benestad, met.no, Oslo, Norway 20.08.2013
# rasmus.benestad@met.no
#------------------------------------------------------------------------
# The equation:
# y = x %psi + %xi

MVR <-function(Y,X,...) UseMethod("MVR")

MVR.default <- function(Y,X,...) {
  print("Don't know what to do - the classes are not the ones I know how to handle")
}


MVR.field <- function(Y,X,SVD=TRUE,LINPACK=FALSE) {
  # Synchronise the two time series objects:
  print("MVR.field")
  history <- attr(X,'history')
  cls <- class(Y)
  colnames(Y) <- paste("Y",1:dim(Y)[2],sep=".")
  colnames(X) <- paste("X",1:dim(Y)[2],sep=".")
  yX <- merge(y,X,all=FALSE)
  Z <- attr(y,pattern)

  ys <- vars[grep('Y',vars)]
  Xs <- vars[grep('X',vars)]
  ix <- is.element(vars,Xs)
  iy <- is.element(vars,ys)
  X <- coredata(comb[,ix])
  Y <- coredata(comb[,iy])
  t <- index(comb)

  if (SVD) {
    if (LINPACK) UWV <-svd(X) else 
                 UWV <-La.svd(X)
    V <- UWV$v; U <- UWV$u; D <- UWV$d

    if (dim(V)[2] != dim(V)[1]) {
      print("V is not a square matrix (no. spatial pts < no. temporal pts)")
      print("Need to swap and transpose SVD products to get correct dimensions")
      if (!LINPACK) {
           V <- UWV$v; U <- UWV$u
      } else {
           V <- t(UWV$v); U <- t(UWV$u)
      }
      psi <- chol2inv(t(V) %*% diag(D^2) %*% V) %*% t(X) %*% Y
    }
  } else {
      print("Warning: may not always be well-posed.")
      psi <- solve(t(X) %*% X) %*% t(X) %*% Y # may be close to singular
  }
      
  Yhat <- X %*% psi
  rmse <- colMeans( (Y - Yhat)^2 )
  R2 <- colSums( Yhat^2/Y^2 )
  
  mvr  <- list(model=psi,fitted.values=Yhat,
               residual=Y - Yhat, r.squared=R2,rmse=rmse)
  attr(mvr,'history') <- history.stamp(x)
  #attr(mvr,'call') <- match.call()
  #attr(mu.pca,'history') <- c('MVR.field',history)
  #attr(mu.pca,'date-stamp') <- date()
  class(mvr) <- c("MVR",class(y))
  invisible(mvr,cls)
}


MVR.eof <- function(Y,X, SVD = SVD, LINPACK = LINPACK) {
  print("MVR.eof")
  history <- attr(X,'history')
  Z <- Y
  cls <- class(Y)
  # Synchronise the two time series objects:
  y <- zoo(coredata(Y),order.by=index(Y))
  x <- zoo(coredata(X),order.by=index(X))
  colnames(y) <- paste("Y",1:dim(y)[2],sep=".")
  colnames(x) <- paste("X",1:dim(x)[2],sep=".")
  YX <- merge(y,x,all=FALSE)
  #str(YX); str(Y); str(X)
  vars <- names(YX)
  #print(vars)
  ys <- vars[grep('Y',vars)]
  Xs <- vars[grep('X',vars)]
  ix <- is.element(vars,Xs)
  iy <- is.element(vars,ys)
  x <- coredata(YX[,ix])
  y <- coredata(YX[,iy])
  Y <- YX[,iy]

  print(dim(x)); print(dim(y))
#  psi <- solve(t(x) %*% x) %*% t(x) %*% y
  xtx <- t(x) %*% x
  print(round(xtx,3))
  psi <- solve(xtx) %*% t(x) %*% y
  Yhat <- zoo(x %*% psi,order.by=index(YX))
  res <- Y - Yhat
  
  nattr <- softattr(Z)
  if (length(nattr)>0) {
    for (i in 1:length(nattr)) {
      attr(Yhat,nattr[i]) <- attr(Z,nattr[i])
      attr(res,nattr[i]) <- attr(Z,nattr[i])
    }
}
  #attr(Yhat,'history') <- c('MVR',attr(Z,'history'))
  #attr(Yhat,'date-stamp') <- date()
  #attr(Yhat,'call') <- match.call()
  attr(Yhat,'history') <- history.stamp(x)
  class(Yhat) <- cls
  #attr(res,'history') <- c('MVR',attr(Z,'history'))
  #attr(res,'date-stamp') <- date()
  #attr(res,'call') <- match.call()
  attr(res,'history') <- history.stamp(x)
  class(res) <- cls
  
  rmse <- colMeans( (coredata(Y) - coredata(Yhat))^2 )
  R2 <- colMeans( coredata(Yhat)^2 )/colMeans( coredata(Y)^2 )
  
  class(psi) <- "matrix"
  mvr  <- list(model=psi,fitted.values=Yhat,data=X,
               residual=res, r.squared=R2,rmse=rmse)

  attr(mvr,'type') <- "MVR"
  #attr(mvr,'call') <- match.call()
  #attr(mvr,'history') <- c('MVR.field',history)
  #attr(mvr,'date-stamp') <- date()
  attr(mvr,'history') <- history.stamp(x)
  class(mvr) <- 'mvr'
  invisible(mvr)
}

MVR.pca <- function(Y,X,SVD=TRUE,LINPACK=FALSE) {
  print("MVR.pca")
  mvr <- MVR.eof(Y,X,SVD=SVD,LINPACK=LINPACK)
  invisible(mvr)
}


# To get one predictor pattern, use predict with newdata set to
# a vector where most variables are set to zero apart from one
# variable set to unity for the identification of teleconnection pattern.

predict.mvr <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) newdata <- object$data
  x <- newdata
  
  psi <- object$model
  Z <- object$fitted.values
  if (inherits(newdata,'zoo'))
    Yhat <- zoo(coredata(x) %*% psi,order.by=index(x)) else
    if (is.vector(x)) Yhat <- t(x) %*% psi
  
  nattr <- softattr(Z)
  for (i in 1:length(nattr)) 
    attr(Yhat,nattr[i]) <- attr(Z,nattr[i])
  class(Yhat) <- class(object$fitted.values)

  if (inherits(Yhat,'eof')) Yhat <- eof2field(Yhat)
  #attr(Yhat,'history') <- c('predict.MVR',attr(Z,'history'))
  #attr(Yhat,'date-stamp') <- date()
  #attr(Yhat,'call') <- match.call()
  attr(Yhat,'history') <- history.stamp(x)
  invisible(Yhat)
}

