#' Multi-variate regression
#' 
#' MVR solves the equation \deqn{Y = \Psi X}{Y = Psi X} and estimates
#' \deqn{\Psi}{Psi} by inverting the equation. Predictions give the value of
#' Y, given this matrix and some input. MVR is useful for data where Y contains
#' several time series where the spatial coherence/covariance is important to
#' reproduce. For instance, Y may be a combination of stations, the two wind
#' components from one station, or a set of different elements from a group of
#' stations.
#' 
#' @aliases MVR MVR.default MVR.field MVR.pca MVR.eof predict.MVR plot.MVR
#'
#' @param Y An object with climate data: field, eof, or pca.
#' @param X Same as Y or any zoo object.
#' @param SVD Use a singular value decomposition as a basis for the PCA.
#' @param LINPACK an option for \code{\link{svd}}.
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @return A CCA object: a list containing a.m, b.m, u.k, v.k, and r,
#' describing the Canonical Correlation variates, patterns and correlations.
#' a.m and b.m are the patterns and u.k and v.k the vectors (time evolution).
#' @author R.E. Benestad
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' # Example for using EOF and MVR
#' slp <- slp.NCEP(lat=c(-40,40),anomaly=TRUE)
#' sst <- sst.NCEP(lat=c(-40,40),anomaly=TRUE)
#' eof.1 <- EOF(slp,mon=1)
#' eof.2 <- EOF(sst,mon=1)
#' mvr <- MVR(eof.1,eof.2)
#' plot(mvr)
#' 
#' # Example for using PCA and MVR
#' oslo <- station(src="NACD",loc="Oslo")
#' bergen <- station.nacd("Bergen")
#' stockholm <- station.nacd("Stockholm")
#' copenhagen <- station.nacd("Koebenhavn")
#' helsinki <- station.nacd("Helsinki")
#' reykjavik <- station.nacd("Stykkisholmur")
#' edinburgh <- station.nacd("Edinburgh")
#' debilt <- station.nacd("De_Bilt")
#' uccle <- station.nacd("Uccle")
#' tromso <- station.nacd("Tromsoe")
#' falun <- station.nacd("Falun")
#' stensele <- station.nacd("Stensele")
#' kuopio <- station.nacd("Kuopio")
#' valentia <- station.nacd("Valentia")
#' X <- combine(oslo,bergen,stockholm,copenhagen,helsinki,reykjavik,
#'            edinburgh,debilt,uccle,tromso,falun,stensele,kuopio,valentia)
#' pca <- PCA(X)
#' slp <- slp.NCEP(lon=c(-20,30),lat=c(30,70))
#' eof <- EOF(slp)
#' mvr <- MVR(pca,eof)
#' plot(mvr)
#' 
#' # Find the teleconnection pattern to the NAO 
#' data("NAOI")
#' data("sunspots")
#' data("NINO3.4")
#' X <- merge(NAOI,sunspots,NINO3.4,all=FALSE)
#' 
#' mvr <- MVR(pca,X)
#' 
#' # Find the pattern for NAOI:
#' teleconnection <- predict(mvr,newdata= c(1,0,0))
#' map(teleconnection,cex=2)
#' }
#' 
#' @export
MVR <- function(Y,X,SVD=TRUE,LINPACK=FALSE,verbose=FALSE) UseMethod("MVR")

#' @export
MVR.default <- function(Y,X,SVD=TRUE,LINPACK=FALSE,verbose=FALSE) {
  print("Don't know what to do - the classes are not the ones I know how to handle")
}

#' @export
MVR.field <- function(Y,X,SVD=TRUE,LINPACK=FALSE,verbose=FALSE) {
  # Synchronise the two time series objects:
  if(verbose) print("MVR.field is not finished. Converting to EOF and redirecting to MVR.eof.")
  eof.Y <- EOF(Y)
  eof.X <- EOF(X)
  mvr <- MVR.eof(Y,X,SVD=SVD,LINPACK=LINPACK,verbose=verbose)
  invisible(mvr)
  # history <- attr(X,'history')
  # cls <- class(Y)
  # colnames(Y) <- paste("Y",1:dim(Y)[2],sep=".")
  # colnames(X) <- paste("X",1:dim(Y)[2],sep=".")
  # YX <- merge(Y,X,all=FALSE)
  # Z <- attr(Y,"pattern")
  # 
  # ys <- vars[grep('Y',vars)]
  # Xs <- vars[grep('X',vars)]
  # ix <- is.element(vars,Xs)
  # iy <- is.element(vars,ys)
  # X <- coredata(comb[,ix])
  # Y <- coredata(comb[,iy])
  # t <- index(comb)
  # 
  # if (SVD) {
  #   if (LINPACK) UWV <-svd(X) else 
  #                UWV <-La.svd(X)
  #   V <- UWV$v; U <- UWV$u; D <- UWV$d
  # 
  #   if (dim(V)[2] != dim(V)[1]) {
  #     print("V is not a square matrix (no. spatial pts < no. temporal pts)")
  #     print("Need to swap and transpose SVD products to get correct dimensions")
  #     if (!LINPACK) {
  #          V <- UWV$v; U <- UWV$u
  #     } else {
  #          V <- t(UWV$v); U <- t(UWV$u)
  #     }
  #     psi <- chol2inv(t(V) %*% diag(D^2) %*% V) %*% t(X) %*% Y
  #   }
  # } else {
  #     print("Warning: may not always be well-posed.")
  #     psi <- solve(t(X) %*% X) %*% t(X) %*% Y # may be close to singular
  # }
  #     
  # Yhat <- X %*% psi
  # rmse <- colMeans( (Y - Yhat)^2 )
  # R2 <- colSums( Yhat^2/Y^2 )
  # 
  # mvr  <- list(model=psi,fitted.values=Yhat,
  #              residual=Y - Yhat, r.squared=R2,rmse=rmse)
  # attr(mvr,'history') <- history.stamp(x)
  # #attr(mvr,'call') <- match.call()
  # #attr(mu.pca,'history') <- c('MVR.field',history)
  # #attr(mu.pca,'date-stamp') <- date()
  # class(mvr) <- c("MVR",cls)
  # invisible(mvr)
}

#' @export
MVR.eof <- function(Y, X, SVD=SVD, LINPACK=LINPACK, verbose=FALSE) {
  if(verbose) print("MVR.eof")
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
  #if(verbose) print(vars)
  ys <- vars[grep('Y',vars)]
  Xs <- vars[grep('X',vars)]
  ix <- is.element(vars,Xs)
  iy <- is.element(vars,ys)
  x <- coredata(YX[,ix])
  y <- coredata(YX[,iy])
  Y <- YX[,iy]

  if(verbose) print(dim(x)); print(dim(y))
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
  #attr(res,'history') <- c('MVR',attr(Z,'history'))n
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

#' @export
MVR.pca <- function(Y,X,SVD=TRUE,LINPACK=FALSE,verbose=FALSE) {
  if(verbose) print("MVR.pca")
  mvr <- MVR.eof(Y,X,SVD=SVD,LINPACK=LINPACK,verbose=verbose)
  invisible(mvr)
}

