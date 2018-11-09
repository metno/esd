# Canonical correlation analysis (CCA) and CCA-based predictions
#
# R.E. Benestad, met.no, Oslo, Norway 20.08.2013
# rasmus.benestad@met.no
#------------------------------------------------------------------------
# Y - first data set
# X - second data set

CCA <-function(Y,X,...) UseMethod("CCA")

CCA.default <- function(Y,X,...) {
  print("Don't know what to do - the classes are not the ones I know how to handle")
}

CCA.eof <- function(Y,X,ip=1:8,verbose=FALSE) {

  if (verbose) print("CCA.eof")
  history <- attr(X,'history')
  Z <- Y
  cls <- class(Y)

  if (inherits(index(Y),c('numeric','integer')))
    index(Y) <- as.Date(paste(index(Y),'01-01',sep='-'))
  if (inherits(index(X),c('numeric','integer')))
    index(X) <- as.Date(paste(index(X),'01-01',sep='-'))
  # Synchronise the two time series objects:
  y <- zoo(coredata(Y),order.by=as.Date(format(index(Y),'%Y-%m-01')))
  x <- zoo(coredata(X),order.by=as.Date(format(index(X),'%Y-%m-01')))
  
  oky <- is.finite(rowMeans(Y)); okx <- is.finite(rowMeans(X))
  #print(length(oky))
  y <- y[oky,]; x <- x[okx,]
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
  #plot(YX)
  #print(dim(YX)); str(x); str(y)

  n.eof1 <- dim(y)[2];   n.eof2 <- dim(x)[2]
  ip <- ip[(ip <= n.eof1) & (ip <= n.eof2)]

  yy <- y[,ip]*attr(Y,'eigenvalues')[ip]
  xx <- x[,ip]*attr(X,'eigenvalues')[ip]
  #print(dim(X1)); print(dim(X2))
  if (verbose) print("Barnett-Preisendorfer CCA")
  S.yx <- cov(yy,xx)
  #print(round(S.yx,4))
  S.yy <- cov(yy,yy)
  #print(round(S.yy,4))
  S.xx <- cov(xx,xx)
  #print(round(S.xx,4)); print("---")
  if (inherits(Y,'eof')) {
    U <- attr(Y,'pattern')[,,ip]
    dU <- dim(U)
    dim(U) <- c(dU[1]*dU[2],dU[3])
  } else if (inherits(Y,'pca')) {
    U <- attr(Y,'pattern')[,ip]
  }
  if (inherits(X,'eof')) {
    V <- attr(X,'pattern')[,,ip]
    dV <- dim(V)
    dim(V) <- c(dV[1]*dV[2],dV[3])
  } else if (inherits(X,'pca')) {
    V <- attr(X,'pattern')[,ip]
  }
  LY <- attr(Y,'eigenvalues')[ip]; LX <- attr(X,'eigenvalues')[ip]

# After Wilks, 1995, p. 401
  info <- "(BP CCA - after Wilks (1995))"
  M.y <- solve(S.yy) %*% S.yx %*% solve(S.xx) %*% t(S.yx)
  M.x <- solve(S.xx) %*% t(S.yx) %*% solve(S.yy) %*% S.yx
  
  b.m <- eigen(M.y)
  a.m <- eigen(M.x)

  #str(b.m)
  B.m <- U %*% diag(LY) %*% Re(t(b.m$vectors)) 
  A.m <- V %*% diag(LX) %*% Re(t(a.m$vectors))
  #str(B.m)
  #print(dim(b.m$vectors))
  w.m <- t( t(b.m$vectors) %*% t(yy[,ip])) 
  v.m <- t( t(a.m$vectors) %*% t(xx[,ip]))
  #print(dim(w.m)); print(dim(y)); print(diag(cor(w.m,v.m)))
  R <- sqrt(Re(b.m$values))
  #print(Re(b.m$values)); print(Re(a.m$values))

  r <- diag(cor(w.m,v.m))
  s <- r < 0
  B.m[,s] <- -B.m[,s]; b.m$vectors[,s] <- -b.m$vectors[,s];
  w.m[,s] <- -w.m[,s]  
  cca <- list(A.m=A.m, B.m = B.m, a.m = a.m, b.m =b.m,
              w.m= w.m, v.m = v.m, r=R,index=index(YX),
              Y=Y,X=X,info=info,ip=ip)

  class(cca) <- c("cca", class(Y)[2])
  
  if (round(R[1],2) != round(cor(w.m[,1],v.m[,1]),2)) {
    print("WARNING: The correlations are not internally consistent!")
    print(paste("CCA: leading canonical correlation=", round(R[1],2),
                " actual correlation=",round(cor(w.m[,1],v.m[,1]),2)))
  }
  attr(cca,'variable') <- paste(varid(Y)[1],varid(X)[1],sep='-')
  attr(cca,'history') <- history.stamp(X,Y)
  invisible(cca)
}


CCA.pca <- function(Y,X,ip=1:8,verbose=FALSE) {
  if (verbose) print("CCA.pca")
  cca <- CCA.eof(Y,X,ip)
  invisible(cca)
}

CCA.field <- function(Y,X,ip=1:8,verbose=FALSE) {
  
 if(verbose) print("CCA.field")
 if(verbose) "print performing EOF analysis and redirecting to CCA.eof"
 eof.y <- EOF(Y,verbose=verbose)
 eof.x <- EOF(X,verbose=verbose)
 cca <- CCA.eof(eof.y,eof.x,ip=ip,verbose=verbose)
 invisible(cca)
}

## KMP 2018-11-02: this code (CCA.field) doesn't work
#   history <- attr(X,'history')
#   Z <- Y
#   cls <- class(Y)
#   # Synchronise the two time series objects:
#   y <- zoo(coredata(Y),order.by=as.Date(format(index(Y),'%Y-%m-01')))
#   x <- zoo(coredata(X),order.by=as.Date(format(index(X),'%Y-%m-01')))
#   colnames(y) <- paste("Y",1:dim(y)[2],sep=".")
#   colnames(x) <- paste("X",1:dim(x)[2],sep=".")
#   YX <- merge(y,x,all=FALSE)
#   #str(YX); str(Y); str(X)
#   vars <- names(YX)
#   #print(vars)
#   ys <- vars[grep('Y',vars)]
#   Xs <- vars[grep('X',vars)]
#   ix <- is.element(vars,Xs)
#   iy <- is.element(vars,ys)
#   x <- coredata(YX[,ix])
#   y <- coredata(YX[,iy])
#   Y <- YX[,iy]
#   
#   x.m <- rowMeans(x); y.m <- rowMeans(y)
#   x <- x - x.m; y <- y - y.m
#   S.yx <- cov(x,y)
#   S.yy <- cov(x,y)
#   S.xx <- cov(x,y)
# 
# # After Wilks, 1995, p. 401
#     ## sub <- paste(sub,"(BP CCA - after Wilks (1995))")
#     M.y <- solve(S.yy) %*% S.yx %*% solve(S.xx) %*% S.yx
#     M.x <- solve(S.xx) %*% S.yx %*% solve(S.yy) %*% S.yx
#     a.m <- eigen(M.y)
#     b.m <- eigen(M.x)
#     # Dimensions of X & Y are ordered as [time,space]
#     if (verbose) {print(dim(a.m)); print(dim(b.m))}
#     w.m <- t(a.m) %*% t(coredata(X))
#     v.m <- t(b.m) %*% t(coredata(Y))
#     R <- sqrt(Re(x.m$values))
#  
#   dim(a.m) <- c(d1[1],d.1[2],d.1[3]); dim(b.m) <- c(d2[1],d.2[2],d.2[3])
# 
#   cca <- list(a.m = a.m, b.m =b.m, w.m= w.m, v.m = v.m, r=R,
#               Y=Y,X=X, info=info)
# 
#   class(cca) <- c("cca", class(Y)[2])
# 
#   if (round(R[1],2) != round(cor(w.m[i1,1],v.m[i2,1]),2)) {
#     print("WARNING: The correlations are not internally consistent!")
#     print(paste("CCA: leading canonical correlation=", round(R[1],2),
#                 " actual correlation=",round(cor(w.m[i1,1],v.m[i2,1]),2)))
#   }

Psi <- function(cca,verbose=FALSE) {
  if(verbose) print("Psi")
  G <- cca$A.m; #d1 <- dim(G); dim(G) <- c(d1[1],d1[2]*d1[3])
  H <- cca$B.m; #d2 <- dim(H); dim(H) <- c(d2[1],d2[2]*d2[3])
  M <- diag(cca$r)
  V <- cca$w.m
  U <- cca$v.m
  G <- t(G); H <- t(H)
  
  # print(dim(G)); print(dim(M)); print(dim(H))
  if (class(cca)[1] =="cca") Psi <- G %*% M %*% solve(t(H) %*% H) %*% t(H)
  #if (class(cca)[1] =="svd") Psi <- G %*% M %*% solve(Cxx) %*% t(H)
  class(Psi) <- paste(class(cca)[1],"model",sep=".")
  attr(Psi,"dims") <- dim(Psi)
  attr(Psi,"lon") <- cca$x1$lon
  attr(Psi,"lat") <- cca$x1$lat
  return(Psi)
}

predict.cca <- function(x, newdata=NULL, ...) {
  if (!is.null(newdata)) X <- newdata else X <- x$X
#predict.CCA <- function(Psi,X) {
  
  #if ( (class(X)[1]!="eof") & (class(X)[1]!="field")) stop('Need a field or EOF object!')
  #type <- class(X)
  #if (type[1]=="eof") X <- EOF2field(X)
  #X <- field$dat
  #d <- dim(X); dim(X) <- c(d[1],d[2]*d[3])
  #X <- t(X)
  #print(dim(Psi)); print(dim(X)); print(d)
  Y.hat <-  Psi(x) %*% X
  #field$dat <- t(Y.hat)
  #print(dim(field$dat))
  #d1 <- attr(Psi,"dims")
  #dim(field$dat) <- c(d[1],d1[2],d1[3])
  #field$lon <- attr(Psi,"lon"); nx <- length(field$lon)
  #field$lat <- attr(Psi,"lat"); ny <- length(field$lat)
  #field$id.x <- rep("CCA",nx*ny)
  #field$id.lon <- rep("CCA",nx)
  #field$id.lat <- rep("CCA",ny)
  #field$id.t <- rep("CCA",d[1])
  #print("HERE")
  #if (type[1]=="eof") result <- EOF(field) else result <- field
  #result
  Y.hat <- attrcp(Y.hat,X)
  Y.hat
}


