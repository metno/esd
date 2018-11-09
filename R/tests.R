test.cca <- function(method="CCA",reconstr=FALSE,mode=1,test=TRUE,LINPACK=TRUE,
                     SVD=TRUE,n.pc=4,synthetic=TRUE,verbose=FALSE) {
  if(verbose) print("test.cca")
  ## KMP 2018-11-09: Don't think this works! The 
  data(eof.slp,envir=environment())
  print(dim(eof.slp$EOF))
  print(dim(eof.slp$PC))
  print(length(eof.slp$W))
  eof.slp$EOF[!is.finite(eof.slp$EOF)] <- 0
  print(summary(c(eof.slp$EOF)))
  print(summary(c(eof.slp$PC)))
  print(summary(c(eof.slp$W)))
  eof.slp$clim <- eof.slp$EOF[1,]*0
  dim(eof.slp$clim) <- c(73,144)
  if (synthetic) {
    nt <- 200
    eof.slp$tim <- 1:nt; eof.slp$yy <- 2000 + floor((1:nt)/360)
    eof.slp$mm <- floor((1:nt)/30)%%12+1; eof.slp$dd <- 1:nt%%30+1 
    eof.slp$size[1,1] <- nt
    #print(rbind(eof.slp$tim,eof.slp$yy,eof.slp$mm,eof.slp$dd))
    eof.slp$PC <- matrix(rep(0,n.pc*nt),nt,n.pc)
    eof.slp$id.t <- rep("test",nt)
  } else {
    nt <- length(eof.slp$tim)
  }
  eof1 <- eof.slp
  eof2 <- eof.slp
  rm(eof.slp)
  gc(reset=TRUE)
  eof1$PC <- eof1$PC[,1:n.pc]
  eof1$EOF <- eof1$EOF[1:n.pc,]
  eof1$W <- eof1$W[1:n.pc]
  eof2$PC <- eof2$PC[,1:n.pc]
  eof2$EOF <- eof2$EOF[1:n.pc,]
  eof2$W <- eof2$W[1:n.pc]
  eof1$dW <- eof1$dW[1:n.pc]
  eof2$dW <- eof2$dW[1:n.pc]
  eof1$var.eof <- eof1$var.eof[1:n.pc]
  eof2$var.eof <- eof2$var.eof[1:n.pc]
  print(dim(eof1$EOF))
  print(dim(eof1$PC))
  print(length(eof1$W))
  y.1 <- as.field(eof1,anomalies=TRUE)
  print(summary(y.1[,1]))
  
  if (synthetic) {
    modes <- 1:n.pc
    modes <- modes[-mode]
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
    x1 <- as.field(eof1)
    x2 <- as.field(eof2)
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

test.ghcnm <- function(verbose=FALSE) {
  stid <- 89001000
  obs <- t2m.GHCNM(stid = stid , verbose=verbose)
  str(obs)
  plot(obs)
}

test.coherence <- function(x=NULL,y=NULL,verbose=FALSE) {
  if(verbose) print("test.coherence")
  default <- FALSE
  if ( (is.null(x)) & (is.null(y)) ) { 
    N <- 1000
    default=TRUE 
  } else if (!is.null(x)) {
    N <- length(x) 
  } else if (!is.null(y)) {
    N <- length(y)
  }
  if ( (!is.null(x)) & (!is.null(y)) & (length(x)!=length(y)) ) {
    stop("testcoh: arguments x and y must have same length")
  }
  if (is.null(x)) x <- cos(pi * seq(0,6,length=N)) +
                       0.5*cos(pi * seq(0,40,length=N)) +
                       0.1*rnorm(N)
  if (is.null(y)) y <- 0.2*sin(pi * seq(0,6,length=N)) +
                       0.3*sin(pi * seq(0,40,length=N)) +
                       0.6*sin(pi * seq(0,20,length=N)) +
                       0.05*rnorm(N)
  plot(x,type="l",main="Test data"); lines(y,col="red")
  if (default) mtext("Common time scales= 167 & 25",1,col="grey")
  dev.new()
  coherence(x,y) -> coh
  if (default) mtext("Common time scales= 167 & 25",1,col="grey")
  invisible(coh)
}

## Would be better to use an example file included with esd
test.regrid <- function(xn=seq(-93,60,by=2),yn=seq(27,80,by=2),verbose=FALSE) {
  if(verbose) print("test.regrid")
  x <- slp.DNMI()
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
  print(paste("The regridded results: length(y)=",length(y),
              "dimensions=",D[1],"x",D[2],"=",D[1]*D[2]))
  class(y) <- class(x)
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

  dim(y) <- c(D[1]*D[2])
  invisible(y)
}

test.regrid2station <- function(x=NULL,y=NULL,verbose=FALSE) {
  if(verbose) print("test.regrid2station")
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

## Name		: test.retrieve.ncdf4
## Description	: This routine contains a series of test functions which compute the global mean 2m-temperature anomalies and other predifined regions based on both CMIP3 and CMIP5 experiments.
## Author 	: Abdelkader Mezghani, METNO
## contact 	: abdelkaderm@met.no
## Created	: 21-03-2013
## Last Updates	: 25-04-2013 ; 22-10-2013 ; 14-05-2014

## Comparing cmip3 to cmip5 over the Arctic
plot.cmip35.global <- function(saveplot=TRUE) { 
  data(global.t2m.cmip5, envir = environment())
  data(global.t2m.cmip3, envir = environment())
  dev.new()
  fig1.zoo(z=list(z1=arctic.t2m.cmip3,z2=arctic.t2m.cmip5),select="noresm")
  if (saveplot) dev.copy2pdf(file="Global_mean_t2m_anomaly_CMIP3-5_1986-2005.pdf")
}

## Comparing cmip3 to cmip5 over Scandinavia
plot.cmip35.global <- function(saveplot=TRUE) { 
  data(scandinavia.t2m.cmip5, envir = environment())
  data(scandinavia.t2m.cmip3, envir = environment())
  dev.new()
  fig1.zoo(z=list(z1=arctic.t2m.cmip3,z2=arctic.t2m.cmip5),select="noresm")
  if (saveplot) dev.copy2pdf(file="Scandinavia_mean_t2m_anomaly_CMIP3-5_1986-2005.pdf")
}

## Comparing cmip3 to cmip5 over the Arctic
plot.cmip35.arctic <- function(saveplot=TRUE) { 
  data(arctic.t2m.cmip5, envir = environment())
  data(arctic.t2m.cmip3, envir = environment())
  dev.new()
  fig1.zoo(z=list(z1=arctic.t2m.cmip3,z2=arctic.t2m.cmip5),select="noresm")
  if (saveplot) dev.copy2pdf(file="Arctic_mean_t2m_anomaly_CMIP3-5_1986-2005.pdf")
}

## test for cmip3 over Scandinavia
test.cmip3.scandinavia <- function(verbose=FALSE) {
  scandinavia.t2m.cmip3 <- test.retrieve.ncdf4(path="CMIP3.monthly/20c3m-sresa1b",param="tas",
                           cntr = "Scandinavia", lon=c(-10,30),lat=c(50,70),
                           experiment="CMIP3",climatology=c(1986,2005),
                           fig=TRUE,saveinfile=NULL,verbose=verbose) 
  save(scandinavia.t2m.cmip3,file="scandinavia.t2m.cmip3.rda")
}

## test for cmip5 over Scandinavia
test.cmip5.scandinavia <- function(verbose=FALSE) {
  scandinavia.t2m.cmip5 <- test.retrieve.ncdf4(path="CMIP5.monthly/rcp45/",param="tas",
                                               cntr = "Scandinavia",lon=c(-10,30),lat=c(50,70),
                                               experiment="CMIP5",climatology=c(1986,2005),fig=TRUE,
                                               saveinfile=NULL,verbose=verbose) 
  save(scandinavia.t2m.cmip5,file="scandinavia.t2m.cmip5.rda")
}

## test cmip3 Global
test.cmip3.global <- function(verbose=FALSE) {
  global.t2m.cmip3 <- test.retrieve.ncdf4(path="CMIP3.monthly/20c3m-sresa1b",param="tas",
                                          cntr = "Global", lon=NULL,lat=NULL,experiment="CMIP3",
                                          climatology=c(1986,2005),fig=TRUE,saveinfile=NULL,
                                          verbose=FALSE) 
  save(global.t2m.cmip3,file="global.t2m.cmip3.rda")
}

## test cmip5 Global
test.cmip5.global <- function(verbose=FALSE) {
  global.t2m.cmip5 <- test.retrieve.ncdf4(path="CMIP5.monthly/rcp45/",param="tas",
                                          cntr = "Global", lon=NULL,lat=NULL,experiment="CMIP5",
                                          climatology=c(1986,2005),fig=TRUE,saveinfile =NULL,
                                          verbose=verbose) 
  save(global.t2m.cmip5,file="global.t2m.cmip5.rda")
}

## test for cmip5 over Arctic
test.cmip5.arctic <- function(...) {
  arctic.t2m.cmip5 <- test.retrieve.ncdf4(path="~/CMIP5.monthly/rcp45/",
                                          param="tas",cntr = "Arctic", lon=NULL,lat=c(60,90),
                                          experiment="CMIP5",climatology=c(1986,2005),fig=TRUE,
                                          saveinfile =TRUE,verbose=FALSE) 
  save(arctic.t2m.cmip5,file="arctic.t2m.cmip5.rda")
}

## test for cmip3 over Arctic
test.cmip3.arctic <- function(...) {
  arctic.t2m.cmip3 <- test.retrieve.ncdf4(path="CMIP3.monthly/20c3m-sresa1b",
                                          param="auto",cntr = "Arctic", lon=NULL,lat=c(60,90),
                                          experiment="CMIP3",climatology=c(1986,2005),
                                          fig=TRUE,saveinfile =TRUE,outfile=NULL) 
  save(arctic.t2m.cmip3,file="arctic.t2m.cmip3.rda")
}

## Function test # 
test.retrieve.ncdf4 <- function(path="CMIP3.monthly/20c3m-sresa1b",param="auto",
                                cntr = "Global", lon=NULL,lat=NULL,experiment="CMIP3",
                                climatology=c(1986,2005),fig=TRUE,saveinfile=TRUE,
                                outfile=NULL,verbose=FALSE) {
  if(verbose) print("test.retrieve.ncdf4")
  # the dash is added to avoid selecting tasmin and tasmax !
  ncfiles <- list.files(path=path,pattern=paste("tas","_",sep=""),full.names=TRUE) 
  if (length(ncfiles)<1) ncfiles <- path
  
  ## Initialize
  vmodel <- rep("-",length(ncfiles))
  for (k in 1:length(ncfiles)) {
    if (substr(ncfiles[k],nchar(ncfiles[k])-1,nchar(ncfiles[k])) != "nc")
      stop("Netcdf file has not been found or the path has not been set correctly !")
    ncid <- nc_open(ncfiles[k])
    print("-------------------------------------------")
    print(ncid$filename)
    gcm <- retrieve.ncdf4(ncfile=ncid,param="auto", lon=lon,lat=lat, lev = NULL,time = NULL,ncdf.check=TRUE,miss2na = TRUE,verbose = verbose)
    project_id <- attr(gcm,"project_id")
    model_id <- attr(gcm,"model_id")
    realization <- attr(gcm,"realization")
    experiment_id <- attr(gcm,"experiment_id")
    title <- attr(gcm,"title")
    ## Compute the spatial average along lon and lat
    gcm2 <- spatial.avg.field(gcm)
    ## Compute the annual mean 
    gcm.am <- aggregate(x=gcm2,by=format(time(gcm2),"%Y"),FUN=mean)
    year <- time(gcm.am)   
    ## Compute the anomalies relative to the climatology
    agcm.ave <- gcm.am-mean(gcm.am[is.element(as.numeric(year),c(climatology[1]:climatology[2]))])
    if (k == 1) {
      X <- matrix(NA,length(year),length(ncfiles)) 
      if (fig)
        plot(agcm.ave,col="steelblue" ,xlim=c(1900,2100),ylim=c(-5,5)) 
    } else if (fig) lines(agcm.ave,col="steelblue")
    if (!is.null(model_id)) { 
      vmodel[k] <- paste(model_id,realization,sep="-")
    }
    X[,k] <- agcm.ave 
  }
  ## output format
  z <- zoo(x=X,order.by=year)
  attr(z,"country") <- cntr
  attr(z,"variable") <- "t2m"
  attr(z,"longname") <- "2m-temperature" 
  attr(z,"timeunit") <- "year"
  attr(z,"climatology") <- paste(climatology,collapse="-")
  attr(z,"model_id") <- vmodel
  attr(z,"project_id") <- project_id 
  attr(z,"experiment") <- experiment
  attr(z,"title") <- title
  attr(z,"author") <- "A. Mezghani, MET Norway , Email : abdelkaderm@met.no"
  attr(z,"Call") <- match.call()
  attr(z,"date") <- date()
  attr(z, "anomaly")        <- TRUE
  attr(z, "time:method")    <- "annual:mean"
  attr(z, "spatial:method") <- "cell:mean"
  attr(z, "title")	    <- "Spatial average of the annual 2m-temperature anomalies"
  attr(z, "URL") 	    <- "http://climexp.knmi.nl/"
  if (!is.null(saveinfile) | (is.logical(saveinfile))) {
    if (is.logical(saveinfile)) {
      outvar <- tolower(paste(cntr,param,project_id,sep="."))
      outfile <-paste(outvar,"rda",sep=".")
    }
    else if (is.character(saveinfile)) outfile <- saveinfile
    attr(z,"file") <- tolower(outfile)
    eval(parse(text=paste(outvar," <- z",sep="")))
    eval(parse(text=paste("save(",outvar,",file='",outfile,"')",sep="")))
  }
  invisible(z)
}


sumisnotna <- function(x) sum(!is.na(x))/length(x)
