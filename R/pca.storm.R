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
  xy <- X[,is.element(colnames(X),param)]
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

plot.pca.storm <- function(X,cex=1.5,new=TRUE,m=2) {
  stopifnot(!missing(X), inherits(X,"pca"))

  pca <- X
  colvec <- c('red3','mediumblue','darkolivegreen3',
              'darkturquoise','darkorange')
  U <- attr(pca,'pattern')
  R2 <- round(100*attr(pca,'eigenvalues')^2/attr(pca,'tot.var'),2)

  if (!is.null(m)) m <- min(m,dim(U)[2])
  else m <- sum(R2>=5)
    
  date <- strptime(attr(pca,'start'),'%Y%m%d%H')
  while (sum(duplicated(date))>0) {
    date[duplicated(date)] <- date[duplicated(date)]+60
  }
  V <- zoo(coredata(pca),order.by=date)
  V.avg <- aggregate(V,FUN="mean",by=as.yearmon(index(V)))
                                     #strftime(index(V),"%Y"))

  if (new) dev.new()
  par( oma=c(1.5,1,1,1.0), mar=c(4,4,2,1) , bty='n' )
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
   widths=c(1.5,1), heights=c(2.5,2))
  
  # Patterns - space
  xlim <- c(min(U[1:10,1:m]),max(U[1:10,1:m]))
  ylim <- c(min(U[11:20,1:m]),max(U[11:20,1:m]))
  plot(0,0,type='n',xlab="",ylab="",xlim=xlim,ylim=ylim,
       main="PCA components")
  for (i in 1:m) {
    lines(U[1:10,i],U[11:20,i],lty=1,col=colvec[i])
    points(U[1:10,i],U[11:20,i],pch='o',col=colvec[i])
    points(U[1,i],U[11,i],col=colvec[i],pch=19)
  }

  # Explained variance - R2
  plot(0,0,type='n',xlim=c(0.5,10),ylim=c(0,100),xlab='EOF #',
       ylab="(%)",main='Explained variance')
  points(R2,type='b',pch=20,col='black')
  for (i in 1:m) {
    points(i,R2[i],col=colvec[i],pch=19)
  }

  # time 
  plot(V.avg[,1],type='n',ylim=c(min(V.avg[,1:m])-5e-3,max(V.avg[,1:m])+5e-3),
       xlab="Time",ylab="",main="")
  for (i in 1:m) {
    lines(V.avg[,i],col=colvec[i])
  }
}

