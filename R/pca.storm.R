# Author 	Kajsa Parding
# Last update   19.02.2015

# Principle Component Analysis (PCA) of storm tracks

PCA.storm <- function(X,neofs=20,param=c('lon','lat','slp'),
                      anomaly=TRUE,verbose=FALSE) {
  stopifnot(!missing(X), inherits(X,"storm"))

  X <- sort.storm(X)
  if (anomaly) { X <- anomaly.storm(X,param)#anomaly.storm(X,param)
  } else {
    i.lon <- colnames(X)=='lon'
    i.dateline <- apply(X,1,function(x) (max(x[i.lon])-min(x[i.lon]))>180)
    lon.dateline <- X[i.dateline,i.lon]
    lon.dateline[lon.dateline<0] <- lon.dateline[lon.dateline<0]+360
    X[i.dateline,i.lon] <- lon.dateline
  }
  i <- sapply(param,function(p) which(colnames(x) %in% p))
  i <- array(i,length(i))
  #xy <- X[,is.element(colnames(X),param)]
  xy <- X[,i]
  D <- dim(xy)

  xyt <- t(coredata(xy))
  neofs <- min(neofs,D[2])
  ok.time <- is.finite(colMeans(xyt))
  z <- xyt[,ok.time]
  ok.site <- is.finite(rowMeans(z))
  z <- z[ok.site,]
 
  pca <- svd(xyt)
  U <- matrix(rep(NA,D[2]*neofs),D[2],neofs)
  U[ok.site,] <- pca$u[,1:neofs]
  V <- matrix(rep(NA,D[1]*neofs),D[1],neofs)
  V[ok.time,] <- pca$v[,1:neofs]
  y <- zoo(V,order.by=index(X)[ok.time])
  names(y) <- paste("X.",1:neofs,sep="")

  invert <- apply(U,2,mean) < 0
  U[,invert] <- -U[,invert]
  y[,invert] <- -y[,invert]

  y <- attrcp(X,y)
  if (anomaly) attr(y,'mean') <- attr(X,'mean') 
  attr(y,'start') <- X[,colnames(X)=='start']
  attr(y,'end') <- X[,colnames(X)=='end']
  attr(y,'n') <- X[,colnames(X)=='n']
  attr(y,'colnames') <- colnames(xy)
  attr(y,'pattern') <- U
  attr(y,'dimensions') <- D
  attr(y,'eigenvalues') <- pca$d[1:neofs]
  attr(y,'sum.eigenv') <- sum(pca$d)
  attr(y,'tot.var') <- sum(pca$d^2)
  attr(y,'history') <- history.stamp(X)
  class(y) <- c("pca",class(X))
  invisible(y)
}


pca2storm <- function(X) {
  stopifnot(!missing(X), inherits(X,"pca"))
  print('pca2storm')
  
  pca <- X
  cls <- class(pca)
  U <- attr(pca,'pattern')
  d <- dim(U)
  W <- attr(pca,'eigenvalues')
  V <- coredata(pca)
  V[!is.finite(V)] <- 0
  x <-U %*% diag(W) %*% t(V)

  x <- cbind(t(x),attr(pca,'start'),attr(pca,'end'),attr(pca,'n'))
  colnames(x) <- c(attr(pca,"colnames"),'start','end','n')

  if (any("anomaly" %in% aspect(pca))) {
    for (i in 1:length(attr(pca,'mean'))) {
      param.i <- names(attr(pca,'mean'))[i]
      mean.i <- unlist(attr(pca,'mean')[i])
      if (length(mean.i)==1) {
        x[,colnames(x)==param.i] <- x[,colnames(x)==param.i] + mean.i
      } else {
        x[,colnames(x)==param.i] <- x[,colnames(x)==param.i] +
           matrix( rep(array(mean.i),sum(colnames(x)==param.i)),
           length(mean.i), sum(colnames(x)==param.i) )
      }
    }
  }

  lon <- x[,colnames(x)=='lon']
  lon[lon>180] <- lon[lon>180]-360
  x[,colnames(x)=='lon'] <- lon

  x <- attrcp(pca,x)
  attr(x,'aspect') <- attr(pca,'aspect')[attr(pca,'aspect')!="anomaly"]
  attr(x,'history') <- history.stamp(pca)
  class(x) <- cls[-1]
  invisible(x)
}

plot.pca.storm <- function(X,cex=1.5,new=TRUE,m=2,param=c('lon','lat')) {

  stopifnot(!missing(X), inherits(X,"storm"))
  if (inherits(X,'pca')) {
    pca <- X; X <- pca2storm(pca)
  } else pca <- PCA.storm(X,param=param)
  
  colvec <- c('red3','mediumblue','darkolivegreen3',
              'darkturquoise','darkorange')
  U <- attr(pca,'pattern')
  R2 <- round(100*attr(pca,'eigenvalues')^2/attr(pca,'tot.var'),2)

  if (!is.null(m)) m <- min(m,dim(U)[2])
  else m <- min(3,sum(R2>=2))
    
  date <- strptime(attr(pca,'start'),'%Y%m%d%H')
  while (sum(duplicated(date))>0) {
    date[duplicated(date)] <- date[duplicated(date)]+60
  }
  V <- zoo(coredata(pca),order.by=date)
  V.mn <- aggregate(V,FUN="mean",by=as.yearmon(index(V)))
  V.yr <- aggregate(V,FUN="mean",by=strftime(index(V),"%Y"))

  if (new) dev.new()
  par( oma=c(1.5,1,1,1.0), mar=c(4,4,2,1) , bty='n' )
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
   widths=c(1.5,1), heights=c(2.5,2))

  param <- unique(attr(pca,"colnames"))
  
  # Patterns - space
  if (length(param)==1) {
    uy <- U
    ux <- matrix(rep(1:(dim(U)[1]),m),dim(U)[1],m)
    xlab <- ""; ylab <- param
  } else if (length(param)==2) {
    ux <- U[attr(pca,"colnames")==param[1],]
    uy <- U[attr(pca,"colnames")==param[2],]
    xlab <- param[1]; ylab <- param[2]
  } else if (length(param)==3) {
    ux <- U[attr(pca,"colnames")==param[1],]
    uy <- U[attr(pca,"colnames")==param[2],]
    uz <- U[attr(pca,"colnames")==param[3],]
    xlab <- param[1]; ylab <- param[2]; zlab <- param[3]
  }
 
  xlim <- c(min(ux),max(ux))
  ylim <- c(min(uy),max(uy))
  plot(0,0,type='n',xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
       main="PCA components")
  for (i in 1:m) {
    lines(ux[,i],uy[,i],lty=1,col=colvec[i])
    points(ux[,i],uy[,i],pch='o',col=colvec[i])
    points(ux[1,i],uy[1,i],col=colvec[i],pch=19)
  }

  # Explained variance - R2
  plot(0,0,type='n',xlim=c(0.5,10),ylim=c(0,100),xlab='EOF #',
       ylab="(%)",main='Explained variance')
  points(R2,type='b',pch=20,col='black')
  for (i in 1:m) {
    points(i,R2[i],col=colvec[i],pch=19)
  }

  # time 
  plot(V.yr[,1],type='n',ylim=c(min(V.yr[,1:m])-5e-3,max(V.yr[,1:m])+5e-3),
       xlab="Time",ylab="",main="")
  lines(index(V.mn),rep(0,length(index(V.mn))),col='grey80',lwd=1.4)
  for (i in 1:m) {
    lines(V.yr[,i],col=colvec[i],lty=1)
    points(V.mn[,i],col=colvec[i],pch=20)
  }
}

