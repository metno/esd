test.regrid <- function(xn=seq(-93,60,by=2),yn=seq(27,80,by=2)) {
#  require(akima)
  load("slp.DNMI.rda")
  x <- slp.DNMI
  X <- t(coredata(x[1,]))
  print(dim(X))
  dim(X) <- attr(x,'dimensions')[1:2]
  beta <- regrid.weights(attr(x,'longitude'),attr(x,'latitude'),xn,yn)
  print("Matrix holding the interpolation weights")
  str(beta); print(dim(beta))
  print(summary(c(beta)))
  print(summary(attr(beta,'npts')))
  print(summary(attr(beta,'chksum')))
#  image(beta)
  #hist(c(beta[beta > 0])); dev.new()
  
  D <- c(length(xn),length(yn))
  y <- apply(cbind(beta,attr(beta,'index')),1,sparseMproduct,X)
  str(y)
  
  #dev.new()
  print(length(y)); print(dim(y)); print(c(length(xn),length(yn)))
  #image(xn,yn,y)
  print(paste("The regridded results: length(y)=",length(y),
              "dimensions=",D[1],"x",D[2],"=",D[1]*D[2]))
  class(y) <- class(x)
#  mostattributes(y) <- attributes(x)
  attr(y,'longitude') <- xn
  attr(y,'latitude') <- yn
  attr(y,'history') <- c(attr(y,'history'),'regrid.field')
  attr(y,'dimensions') <- c(D[1],D[2],1)
  
  print("map")
  map(y)
  print(attr(x,'dimensions')); print(dim(x))

  contour(attr(x,'longitude'),attr(x,'latitude'),X)

  x1 <- attr(slp.DNMI,'longitude');  nx1 <- length(x1)
  y1 <- attr(slp.DNMI,'latitude');   ny1 <- length(y1)
  xy1 <- rep(x1,ny1); yx1 <- sort(rep(y1,nx1))
  print(c(length(xy1),length(yx1),length(X)))

#  z <- interp(xy1,yx1,X,lon,lat)$z
#  contour(xn,yn,z,add=TRUE,col="red",lty=2)

  dim(y) <- c(D[1]*D[2])
  invisible(y)

}

test.regrid2station <- function(x=NULL,y=NULL) {
  if (is.null(x)) {
    load("t2m.DNMI.rda")
    x <- t2m.DNMI
  }
  if (is.null(y)) {
    load("Oslo.rda")
    y <- Oslo
  }
  X <- coredata(x)
  beta <- regrid.weights(attr(x,'longitude'),attr(x,'latitude'),
                         attr(y,'longitude'),attr(y,'latitude'))
  print(dim(beta))
  print(summary(c(beta)))
  print(summary(attr(beta,'npts')))
  print(summary(attr(beta,'chksum')))
  d <- dim(x)
  Z <- rep(NA,d[1])
  for (i in 1:d[1]) {
    z <- apply(beta,1,sparseMproduct,coredata(X[i,]))
    #print(c(i,d[1],length(z),length(y[,i]),NA,dim(x),dim(y)))
    Z[i] <- z
  }
  print(summary(Z)); print(length(Z))
  Z <- zoo(Z,order.by=index(x))
  #print(format(index(y),'%Y'))
  plot(Z,lty=3)
  lines(aggregate(Z,by=as.numeric(format(index(Z),'%Y')),mean),lwd=3)
  year <- as.numeric(format(index(y),'%Y'))
  #print(table(year))
  class(y) <- "zoo"
  print(aggregate(y,by=year,mean))
  lines(y,col="red",lty=3)
  lines(aggregate(y,by=year,mean),col="red",lty=2)
}
