#' Principle component analysis of trajectory objects.
#' 
#' Computes principal component analysis for trajectory data, e.g., storm
#' tracks.  Add some reference and details about the method.  The PCA is based
#' on \code{\link{svd}}.
#' 
#' @aliases PCA.trajectory plot.pca.trajectory
#' 
#' @param X a 'trajectory' object
#' @param verbose TRUE - clutter the screen with messages
#' @param anomaly logical. If TRUE, subtract the first latitude/longitude from
#' each trajectory.
#' @param param parameters to include in principle component analysis.
#' 
#' @keywords spatial ts multivariate
#' 
#' @examples
#' # Simple EOF for annual mean SST:
#' data(imilast.M03)
#' x <- subset(imilast.M03,is=list(lon=c(-20,20),lat=c(50,70)))
#' # PCA of longitude and latitude
#' pca <- PCA.trajectory(x,param=c('lon','lat'))
#' plot.pca.trajectory(pca)
#' map.pca.trajectory(pca,projection='latlon')
#' 
#' # latitude only
#' pca <- PCA.trajectory(x,param=c('lat'))
#' plot.pca.trajectory(pca)
#' 
#' @export
PCA.trajectory <- function(X,...,neofs=20,param=c('lon','lat'),
                      anomaly=TRUE,verbose=FALSE) {
  if(verbose) print("PCA.trajectory")
  stopifnot(!missing(X), inherits(X,"trajectory"))

  X <- sort(X)
  if (anomaly) {
    p <- param[!param %in% names(attr(X,'mean'))]
    if(length(p)>0) {
      if(verbose) print('calculating anomaly')
      if(verbose) print(p)
      X <- anomaly.trajectory(X,param=p)
    }
  } else if (!anomaly & "anomaly" %in% attr(X,'aspect')) {
    X <- anomaly2trajectory(X)
  }
  if ("lon" %in% param & !anomaly) {
    i.lon <- colnames(X)=='lon'
    i.dateline <- apply(X,1,function(x) (max(x[i.lon])-min(x[i.lon]))>180)
    lon.dateline <- X[i.dateline,i.lon]
    lon.dateline[lon.dateline<0] <- lon.dateline[lon.dateline<0]+360
    X[i.dateline,i.lon] <- lon.dateline
  }
  i <- sapply(param,function(p) which(colnames(X) %in% p))
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

  U.anom <- U
  for (p in param) {
    for (i in 1:neofs) {
      U.anom[colnames(xy)==p,i] <- U[colnames(xy)==p,i]-U[colnames(xy)==p,i][1]
    }
  }  
  invert <- apply(U.anom,2,mean) < 0#apply(U,2,mean) < 0
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

pca2trajectory <- function(X,verbose=FALSE) {
  stopifnot(!missing(X), inherits(X,"pca"))
  if(verbose) print('pca2trajectory')
  
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
  x <- attrcp(pca,x)
  attr(x,'aspect') <- attr(pca,'aspect')
  attr(x,'history') <- history.stamp(pca)
  class(x) <- 'trajectory'

  if (any("anomaly" %in% attr(x,'aspect'))) {
    x <- anomaly2trajectory(x)
  }
  invisible(x)
}

#' @export
plot.pca.trajectory <- function(x,...,cex=1.5,new=TRUE,m=2,param=c('lon','lat'),
                           main=NULL,verbose=FALSE) {

  if(verbose) print("plot.pca.trajectory")
  stopifnot(!missing(x), inherits(x,"trajectory"))
  if (inherits(x,'pca')) {
    pca <- x; x <- pca2trajectory(pca)
  } else pca <- PCA.trajectory(x,param=param)
  
  if(verbose) print("Extractpatterns, PCs and eigenvalues")
  colvec <- c('red3','mediumblue','darkolivegreen3',
              'darkturquoise','darkorange')
  U <- attr(pca,'pattern')
  R2 <- round(100*attr(pca,'eigenvalues')^2/attr(pca,'tot.var'),2)

  if (!is.null(m)) m <- min(m,dim(U)[2])
  else m <- min(3,sum(R2>=2))

  if(verbose) print("Aggregate time series")
  date <- strptime(attr(pca,'start'),'%Y%m%d%H')
  while (sum(duplicated(date))>0) {
    date[duplicated(date)] <- date[duplicated(date)]+60
  }
  V <- zoo(coredata(pca),order.by=date)
  V.mn <- aggregate(V,FUN="mean",by=as.yearmon(index(V)))
  V.yr <- aggregate(V,FUN="mean",by=strftime(index(V),"%Y"))

  param <- unique(attr(pca,"colnames"))

  if(verbose) print("Arrange patterns")
  # Patterns - space
  if (length(param)==1) {
    uy <- U
    ux <- matrix(rep(1:(dim(U)[1]),m),dim(U)[1],m)
    xlab <- "step"
    ylab <- paste(param,"component")
  } else if (length(param)==2) {
    ux <- U[attr(pca,"colnames")==param[1],]
    uy <- U[attr(pca,"colnames")==param[2],]
    xlab <- paste(param[1],"component")
    ylab <- paste(param[2],"component")
  } else if (length(param)==3) {
    ux <- U[attr(pca,"colnames")==param[1],]
    uy <- U[attr(pca,"colnames")==param[2],]
    uz <- U[attr(pca,"colnames")==param[3],]
    xlab <- paste(param[1],"component")
    ylab <- paste(param[2],"component")
    zlab <- paste(param[3],"component")
  }

  xlim <- c(min(ux[,1:m]),max(ux[,1:m]))
  ylim <- c(min(uy[,1:m]),max(uy[,1:m]))
 
  if (new) dev.new()
  par( oma=c(1.5,1,1,1.0), mar=c(4,4,2,1) , bty='n' )
  if (length(param)<3) {
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
     widths=c(1.5,1), heights=c(2.5,2))
  } else if (length(param)==3) {
    layout(matrix(c(1,1,1,2,2,2,3,3,4,4,4,4), 2, 6, byrow = TRUE),
     widths=c(1,1,1,1,1,1), heights=c(2,2.5))
  }

  if(verbose) print("Plot map")
  plot(0,0,type='n',xlab=xlab,
       ylab=ylab,xlim=xlim,ylim=ylim,main="loading pattern")
  for (i in 1:m) {
    lines(ux[,i],uy[,i],lty=1,col=colvec[i])
    points(ux[,i],uy[,i],pch='o',col=colvec[i])
    points(ux[1,i],uy[1,i],col=colvec[i],pch=19)
  }

  if (length(param)==3) {
    zlim <- c(min(uz[,1:m]),max(uz[,1:m]))
    plot(0,0,type='n',xlab=zlab,ylab=ylab,
         xlim=zlim,ylim=ylim,main="loading pattern")
    for (i in 1:m) {
      lines(uz[,i],uy[,i],lty=1,col=colvec[i])
      points(uz[,i],uy[,i],pch='o',col=colvec[i])
      points(uz[1,i],uy[1,i],col=colvec[i],pch=19)
    }
  }
    
  # Explained variance - R2
  plot(0,0,type='n',xlim=c(0.5,10),ylim=c(0,100),xlab='EOF #',
     ylab="(%)",main='explained variance')
  points(R2,type='b',pch=20,col='black')
  for (i in 1:m) {
    points(i,R2[i],col=colvec[i],pch=19)
  }

  # time 
  plot(V.yr[,1],type='n',ylim=c(min(V.mn[,1:m])-5e-3,max(V.mn[,1:m])+5e-3),
     xlab="Time",ylab="PC",main="")
  #lines(index(V.mn),rep(0,length(index(V.mn))),col='grey80',lwd=1.4)
  lines(as.Date(c("1900-01-01","2020-01-01")),c(0,0),col='grey80',lwd=1)
  for (i in 1:m) {
    #lines(V.yr[,i],col=colvec[i],lty=1)
    #points(V.mn[,i],col=colvec[i],pch=20)
    lines(V.mn[,i],col=colvec[i],lty=1)
    #lines(trend(V.mn[,i]),col=colvec[i],lty=2,lwd=2)
  }

  title(main=main,outer=T)
  invisible(pca)
}

