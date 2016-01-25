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

CCA.eof <- function(Y,X,i.eofs=1:8) {

  print("CCA.eof")
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
  i.eofs <- i.eofs[(i.eofs <= n.eof1) & (i.eofs <= n.eof2)]

  yy <- y[,i.eofs]*attr(Y,'eigenvalues')[i.eofs]
  xx <- x[,i.eofs]*attr(X,'eigenvalues')[i.eofs]
  #print(dim(X1)); print(dim(X2))
  print("Barnett-Preisendorfer CCA")
  S.yx <- cov(yy,xx)
  #print(round(S.yx,4))
  S.yy <- cov(yy,yy)
  #print(round(S.yy,4))
  S.xx <- cov(xx,xx)
  #print(round(S.xx,4)); print("---")
  if (inherits(Y,'eof')) {
    U <- attr(Y,'pattern')[,,i.eofs]
    dU <- dim(U)
    dim(U) <- c(dU[1]*dU[2],dU[3])
  } else if (inherits(Y,'pca')) {
    U <- attr(Y,'pattern')[,i.eofs]
  }
  if (inherits(X,'eof')) {
    V <- attr(X,'pattern')[,,i.eofs]
    dV <- dim(V)
    dim(V) <- c(dV[1]*dV[2],dV[3])
  } else if (inherits(X,'pca')) {
    V <- attr(X,'pattern')[,i.eofs]
  }
  LY <- attr(Y,'eigenvalues')[i.eofs]; LX <- attr(X,'eigenvalues')[i.eofs]

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
  w.m <- t( t(b.m$vectors) %*% t(yy[,i.eofs])) 
  v.m <- t( t(a.m$vectors) %*% t(xx[,i.eofs]))
  #print(dim(w.m)); print(dim(y)); print(diag(cor(w.m,v.m)))
  R <- sqrt(Re(b.m$values))
  #print(Re(b.m$values)); print(Re(a.m$values))

  r <- diag(cor(w.m,v.m))
  s <- r < 0
  B.m[,s] <- -B.m[,s]; b.m$vectors[,s] <- -b.m$vectors[,s];
  w.m[,s] <- -w.m[,s]  
  cca <- list(A.m=A.m, B.m = B.m, a.m = a.m, b.m =b.m,
              w.m= w.m, v.m = v.m, r=R,index=index(YX),
              Y=Y,X=X,info=info,i.eofs=i.eofs)

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


CCA.pca <- function(Y,X,i.eofs=1:8) {
  print("CCA.pca")
  cca <- CCA.eof(Y,X,i.eofs)
  invisible(cca)
}

CCA.field <- function(Y,X,i.eofs=1:8) {
  
  print("CCA.field")
  history <- attr(X,'history')
  Z <- Y
  cls <- class(Y)
  # Synchronise the two time series objects:
  y <- zoo(coredata(Y),order.by=as.Date(format(index(Y),'%Y-%m-01')))
  x <- zoo(coredata(X),order.by=as.Date(format(index(X),'%Y-%m-01')))
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
  
  x.m <- rowMeans(x); y.m <- rowMeans(y)
  x <- x - x.m; y <- y - y.m
  S.yx <- cov(x,y)
  S.yy <- cov(x,y)
  S.xx <- cov(x,y)

# After Wilks, 1995, p. 401
    ## sub <- paste(sub,"(BP CCA - after Wilks (1995))")
    M.y <- solve(S.yy) %*% S.yx %*% solve(S.xx) %*% S.yx
    M.x <- solve(S.xx) %*% S.yx %*% solve(S.yy) %*% S.yx
    a.m <- eigen(M.y)
    b.m <- eigen(M.x)
    # Dimensions of X & Y are ordered as [time,space]
    print(dim(a.m)); print(dim(b.m))
    w.m <- t(a.m) %*% t(coredata(X))
    v.m <- t(b.m) %*% t(coredata(Y))
    R <- sqrt(Re(x.m$values))
 
  dim(a.m) <- c(d1[1],d.1[2],d.1[3]); dim(b.m) <- c(d2[1],d.2[2],d.2[3])

  cca <- list(a.m = a.m, b.m =b.m, w.m= w.m, v.m = v.m, r=R,
              Y=Y,X=X, info=info)

  class(cca) <- c("cca", class(Y)[2])

  if (round(R[1],2) != round(cor(w.m[i1,1],v.m[i2,1]),2)) {
    print("WARNING: The correlations are not internally consistent!")
    print(paste("CCA: leading canonical correlation=", round(R[1],2),
                " actual correlation=",round(cor(w.m[i1,1],v.m[i2,1]),2)))
  }

  invisible(cca)
}






test.cca <- function(method="CCA",reconstr=FALSE,mode=1,test=TRUE,LINPACK=TRUE,
                     SVD=TRUE,n.pc=4,synthetic=TRUE) {
  print("version 0.1:")
  data(eof.slp,envir=environment())
  print(dim(eof.slp$EOF)); print(dim(eof.slp$PC)); print(length(eof.slp$W))
  eof.slp$EOF[!is.finite(eof.slp$EOF)] <- 0
  print(summary(c(eof.slp$EOF)))
  print(summary(c(eof.slp$PC)))
  print(summary(c(eof.slp$W)))
  eof.slp$clim <- eof.slp$EOF[1,]*0
  dim(eof.slp$clim) <- c(73,144)
  if (synthetic) {
    nt <- 200
    eof.slp$tim <- 1:nt; eof.slp$yy <- 2000 + floor((1:nt)/360)
    eof.slp$mm <- mod(floor((1:nt)/30),12)+1; eof.slp$dd <- mod(1:nt,30)+1 
    eof.slp$size[1,1] <- nt
    #print(rbind(eof.slp$tim,eof.slp$yy,eof.slp$mm,eof.slp$dd))
    eof.slp$PC <- matrix(rep(0,n.pc*nt),nt,n.pc)
    eof.slp$id.t <- rep("test",nt)
  } else nt <- length(eof.slp$tim)
  eof1 <- eof.slp; eof2 <- eof.slp; rm(eof.slp); gc(reset=TRUE)
  eof1$PC <- eof1$PC[,1:n.pc]; eof1$EOF <- eof1$EOF[1:n.pc,]; eof1$W <- eof1$W[1:n.pc]
  eof2$PC <- eof2$PC[,1:n.pc]; eof2$EOF <- eof2$EOF[1:n.pc,]; eof2$W <- eof2$W[1:n.pc]
  eof1$dW <- eof1$dW[1:n.pc]; eof2$dW <- eof2$dW[1:n.pc]
  eof1$var.eof <- eof1$var.eof[1:n.pc]; eof2$var.eof <- eof2$var.eof[1:n.pc]

  print(dim(eof1$EOF)); print(dim(eof1$PC)); print(length(eof1$W))
  y.1 <- EOF2field(eof1,anomalies=TRUE)$dat[,1,1]
  print(summary(y.1))
  
  if (synthetic) {
    modes <- 1:n.pc; modes <- modes[-mode]
    print("signal: sin")
    eof1$PC[,mode] <- 50*sin(seq(-12*pi,12*pi,length=nt))
    eof2$PC[,mode] <- 50*sin(seq(-12*pi,12*pi,length=nt))
    print("noise: rnorm")
    for (i in modes) {
      eof1$PC[,i] <- rnorm(nt)
      eof2$PC[,i] <- rnorm(nt)
    }
  }
  if (reconstr) {
    print("Reconstruct fields...")
    x1 <- EOF2field(eof1)
    x2 <- EOF2field(eof2)
    print(class(x1))
    print("Run test...")
    print(paste(method,"(x1,x2,test=",test,
                 ",LINPACK=",LINPACK,",SVD=",SVD,")",sep=""))
    cca.test <- eval(parse(text=paste(method,"(x1,x2,test=",test,
                 ",LINPACK=",LINPACK,",SVD=",SVD,")",sep="")))
  } else {
    print(paste(method,"(eof1,eof2,test=",test,
                 ",LINPACK=",LINPACK,",SVD=",SVD,")",sep=""))
    cca.test <- eval(parse(text=paste(method,"(eof1,eof2,test=",test,
                 ",LINPACK=",LINPACK,",SVD=",SVD,")",sep="")))
  }
  invisible(cca.test)
}

predict.cca <- function(object, newdata=NULL, ...) {


#predict.CCA <- function(Psi,X) {


  Psi <- function(cca) {
    G <- cca$a.m; d1 <- dim(G); dim(G) <- c(d1[1],d1[2]*d1[3])
    H <- cca$b.m; d2 <- dim(H); dim(H) <- c(d2[1],d2[2]*d2[3])
    M <- diag(cca$r)
    V <- cca$w.m
    U <- cca$v.m
    G <- t(G); H <- t(H)
 
  # print(dim(G)); print(dim(M)); print(dim(H))
    if (class(cca)[1] =="CCA") Psi <- G %*% M %*% solve(t(H) %*% H) %*% t(H)
  #if (class(cca)[1] =="SVD") Psi <- G %*% M %*% solve(Cxx) %*% t(H)
    class(Psi) <- paste(class(cca)[1],"model",sep=".")
    attr(Psi,"dims") <- d1
    attr(Psi,"lon") <- cca$x1$lon
    attr(Psi,"lat") <- cca$x1$lat
    Psi
  }
  
  if ( (class(X)[1]!="eof") & (class(X)[1]!="field")) stop('Need a field or EOF object!')
  type <- class(X)
  if (type[1]=="eof") field <- EOF2field(X)
  X <- field$dat
  d <- dim(X); dim(X) <- c(d[1],d[2]*d[3])
  X <- t(X)
  #print(dim(Psi)); print(dim(X)); print(d)
  Y.hat <-  Psi %*% X
  field$dat <- t(Y.hat)
  #print(dim(field$dat))
  d1 <- attr(Psi,"dims")
  dim(field$dat) <- c(d[1],d1[2],d1[3])
  field$lon <- attr(Psi,"lon"); nx <- length(field$lon)
  field$lat <- attr(Psi,"lat"); ny <- length(field$lat)
  field$id.x <- rep("CCA",nx*ny)
  field$id.lon <- rep("CCA",nx)
  field$id.lat <- rep("CCA",ny)
  field$id.t <- rep("CCA",d[1])
  #print("HERE")
  if (type[1]=="eof") result <- EOF(field) else result <- field
  result
}


