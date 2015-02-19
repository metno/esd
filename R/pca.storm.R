# Author 	Kajsa Parding
# Last update   19.02.2015

# Principle Component Analysis (PCA) of storm tracks

PCA.storm <- function(X,neofs=20,param=c('lon','lat','slp'),
                      anomaly=TRUE,verbose=FALSE) {
  stopifnot(!missing(X), inherits(X,"storm"))

  X <- sort.storm(X)
  if (anomaly) X <- anomaly.storm(X,param)
  else {
    i.lon <- colnames(X)=='lon'
    i.dateline <- apply(X,1,function(x) (max(x[i.lon])-min(x[i.lon]))>180)
    lon.dateline <- X[i.dateline,i.lon]
    lon.dateline[lon.dateline<0] <- lon.dateline[lon.dateline<0]+360
    X[i.dateline,i.lon] <- lon.dateline
  }
  xy <- X[,is.element(colnames(X),param)]
  D <- dim(xy)

  xyt <- t(coredata(xy))
  neofs <- min(neofs,D[1])
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
      if (param.i=='slp') {
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



plot.pca.storm <- function(X,cex=1.5,new=TRUE,m=3) {
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


#============= THE END ==============

  # Show patterns in units of lon-lat. 
  #W <- attr(y,'eigenvalues')
  #dev.new()
  #plot(0,0,type='n',xlim=c(-90,90),ylim=c(-90,90))
  #for (i in 1:m) {
  #  PCi <- median(V[,i])*(U[,i]*W[i])
  #  lines(PCi[1:10],PCi[11:20],lty=1,col=colvec[i])
  #  points(PCi[1],PCi[11],col=colvec[i],pch=19)
  #  PCi <- max(V[,i])*(U[,i]*W[i])
  #  lines(PCi[1:10],PCi[11:20],lty=3,col=colvec[i])
  #  points(PCi[1],PCi[11],col=colvec[i],pch=20)
  #  PCi <- min(V[,i])*(U[,i]*W[i])
  #  lines(PCi[1:10],PCi[11:20],lty=4,col=colvec[i])
  #  points(PCi[1],PCi[11],col=colvec[i],pch=20)
  #}
#}

 # Set scale for colour scheme
  #str(y)
  ## a.T <- matrix(rep(NA,4*N),4,N)
  ## ax <- quantile(abs(attr(y,'mean')),0.99,na.rm=TRUE)
  ## if (min(attr(y,'mean'))<0) scale0 <- seq(-ax,ax,length=nc) else
  ##                            scale0 <- seq(0,ax,length=nc)
  ## ax <- quantile(abs(attr(y,'pattern')),0.99,na.rm=TRUE)
  ## scale <- seq(-ax,ax,length=nc)

  ## #print("here")
  ## for (i in 1:N) {
  ##   a.T[1,i] <-  sum(attr(y,'mean')[i] > scale0)
  ##   for (j in 1:m) 
  ##     a.T[j+1,i] <-  sum(attr(y,'pattern')[i,j] > scale)
  ## }
  ## a.T[a.T < 1] <- 1; a.T[a.T > 100] <- 100
  
  #if (new) dev.new(width=7,height=9)
  #par(mfrow=c(3,2),mar=c(3.5,3,3.5,3),bty="n",xaxt="n",yaxt="n")

  ## plot(lon,lat,
  ##      main="Climatology",
  ##      col=col[a.T[1,]],pch=19,xlab="",ylab="",cex=cex)
  ## points(lon,lat,cex=cex)
  ## data(geoborders,envir=environment())
  ## lines(geoborders,col='grey40')
  ## lines(geoborders$x - 360,geoborders$y,col='grey40')
  ## points(lon,lat,cex=cex,col=col[a.T[1,]],pch=19)

  ## plot(lon,lat,
  ##      main=paste("EOF #1:",R2[1],"% of variance"),
  ##      col=col[a.T[2,]],pch=19,xlab="",ylab="",cex=cex)
  ## points(lon,lat,cex=cex)
  ## lines(geoborders)
  ## lines(geoborders$x - 360,geoborders$y)
  ## points(lon,lat,cex=cex,col=col[a.T[2,]],pch=19)

  ## plot(lon,lat,
  ##      main=paste("EOF #2:",R2[2],"% of variance"),
  ##      col=col[a.T[3,]],pch=19,xlab="",ylab="",cex=cex)
  ## points(lon,lat,cex=cex)
  ## lines(geoborders,col='grey40')
  ## lines(geoborders$x - 360,geoborders$y,col='grey40')
  ## points(lon,lat,cex=cex,col=col[a.T[3,]],pch=19)

  ## plot(lon,lat,
  ##      main=paste("EOF #3:",R2[3],"% of variance"),
  ##      col=col[a.T[4,]],pch=19,xlab="",ylab="",cex=cex)
  ## points(lon,lat,cex=cex)
  ## lines(geoborders,col='grey40')
  ## lines(geoborders$x - 360,geoborders$y,col='grey40')
  ## points(lon,lat,cex=cex,col=col[a.T[4,]],pch=19)

  ## par(mar=c(1,0,0,0),fig=c(0.1,0.3,0.665,0.695),new=TRUE,cex.axis=0.6)
  ## image(cbind(1:nc,1:nc),col=col)
  ## nl <- pretty(scale0)
  ## par(xaxt="s")
  ## axis(1,at=seq(0,1,length=length(nl)),label=nl)

  ## par(mar=c(1,0,0,0),fig=c(0.1,0.3,0.32,0.35),new=TRUE,cex.axis=0.6,xaxt="n")
  ## image(cbind(1:nc,1:nc),col=col)
  ## nl <- pretty(scale)
  ## par(xaxt="s")
  ## axis(1,at=seq(0,1,length=length(nl)),label=nl)

  ## par(mar=c(1,0,0,0),fig=c(0.6,0.8,0.665,0.695),new=TRUE,cex.axis=0.6,xaxt="n")
  ## image(cbind(1:nc,1:nc),col=col)
  ## nl <- pretty(scale)
  ## par(xaxt="s")
  ## axis(1,at=seq(0,1,length=length(nl)),label=nl)

  ## par(mar=c(1,0,0,0),fig=c(0.6,0.8,0.32,0.35),new=TRUE,cex.axis=0.6,xaxt="n")
  ## image(cbind(1:nc,1:nc),col=col)
  ## nl <- pretty(scale)
  ## par(xaxt="s")
  ## axis(1,at=seq(0,1,length=length(nl)),label=nl)
  
  ## par(mfcol=c(1,1),fig=c(0,1,0,0.33),new=TRUE,xaxt="s",yaxt="n",bty="n",
  ##     mar=c(2,2,1,1))
  ## ylim <- 2*range(coredata(y[,1:m]),na.rm=TRUE)
  ## plot(y[,1]+0.5*ylim[2],lwd=2,ylim=ylim)
  ## grid()
  ## col <- c("red","blue")
  ## for (j in 1:m) lines(y[,j+1]+(1-j)*0.5*ylim[2],lwd=2,col=col[j])
  ## legend(index(y)[1],ylim[1],c('PC 1','PC 2','PC 3'),
  ##        col=c('black','red','blue'),bty='n',lwd=2)
  ## invisible(a.T)
#}
