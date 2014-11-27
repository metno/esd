# Computes Empirical Orthogonal Functions (EOFs)
#
# R.E. Benestad, 
# rasmus.benestad@met.no
#
#------------------------------------------------------------------------

require(zoo)

EOF<-function(X,it=NULL,is=NULL,n=20,lon=NULL,lat=NULL,verbose=FALSE,...)
  UseMethod("EOF")

EOF.default <- function(X,it=NULL,is=NULL,n=20,lon=NULL,lat=NULL,
                        area.mean.expl=FALSE,verbose=FALSE,...) {
  # Verify Arguments
  if (verbose) print("EOF.default")
  stopifnot(!missing(X), is.matrix(X),inherits(X,"zoo"))

  if ( !zeros(inherits(X,c("comb","zoo"),which=TRUE)) )
         eof <- EOF.comb(X,it=it,is=is,n=n,
                         area.mean.expl=area.mean.expl,
                         verbose=verbose) else
  if ( !zeros(inherits(X,c("field","zoo"),which=TRUE)) )
         eof <- EOF.field(X,it=it,is=is,n=n,
                       area.mean.expl=area.mean.expl,
                          verbose=verbose) 
  return(eof)
}


# Apply EOF analysis to the monthly mean field values:

EOF.field <- function(X,it=NULL,is=NULL,n=20,lon=NULL,lat=NULL,
                      area.mean.expl=FALSE,verbose=FALSE) {

  SF <- function(x) {sum(is.finite(x))}

  if (verbose) print("EOF.field")
  stopifnot(!missing(X), is.matrix(X),
            inherits(X,c("field","zoo")))

  # REB: 29.04.2014
  if (area.mean.expl) {
    if (verbose) print('area.mean.expl')
    A <- aggregate.area(X,FUN='mean')
    x <- X - A
    x <- attrcp(X,x)
    attr(x,'dimensions') <- attr(X,'dimensions')
    class(x) <- class(X)
    X <- x
  }
  
  X <- sp2np(X)
  dates <- index(X)
  d <- attr(X,'dimensions')
  cls <- class(X)
  #print(cls)
  browser()
  x <- subset(X,it=it,is=is)
#  if (!is.null(it)) {
#    if (verbose) print(paste('temporal subset: it=',it))
#    #print(it)
#    #print(table(as.POSIXlt(dates)$mon+1))
#    # Select a subset of the months
#    if ( (min(it) > 0) & (max(it) < 13) & (inherits(X,c("month"))) ) {
#      #keepm <- as.numeric(format(index(X),"%m"))==it
#      #print("Monthly aggregated field")
#      keepm <- is.element(as.POSIXlt(dates)$mon+1,it)
#    } else 
#    if ( (min(it) > 0) & (max(it) < 5) & (inherits(X,c("season"))) ) {
#      #print("Seasonally aggregated field")
#      #print(table(as.POSIXlt(dates)$mon+1))
#      keepm <- is.element(as.POSIXlt(dates)$mon+1,c(1,4,7,10)[it])
#      #print(c(it,sum(keepm)))
#    } else
#    ## Select a range of years (interval)
#    if ( (length(it)>1) & (sum(is.element(it,1500:2500))>0) ) {
#      if (length(it)==2) it <- it[1]:it[2]
#      ##print(it); print(as.POSIXlt(dates)$year+1900); print(dates)
#      keepm <- is.element(as.numeric(format(index(X),"%Y")),it)
#      ## AM replacement 13-11-2013 old line keepm <- is.element(as.POSIXlt(dates)$year+1900,it)
#    } else ## if (inherits(it,"POSIXt"))
#    keepm <- is.element(as.Date(dates),as.Date(it))
#    ## AM replacement 13-11-2013 old line keepm <- is.element(as.POSIXlt(dates),as.POSIXlt(it))
#    dates <- dates[keepm]
#    x <- zoo(X[keepm,],order.by = dates)
#    d[3] <- sum(keepm)
#    #print(d[3])
#  } else x <- X
  Y <- t(coredata(x))
  #print(dim(Y)); print(d)
  
#  # to select geographicla regions, the zoo aspects are no longer needed:
#  # expand into lon lat dimensions in addition to time:
#  dim(Y) <- d
#  #print(d); A <- Y[,,1]; image(t(A)); dev.new()
#  if (!is.null(lon)) {
#    if (length(lon) != 2)
#      warning("EOF: argument 'lon' must be a range") else { 
#        # Select a subset of the longitudes        
#        keepx <- (attr(X,'longitude') >= min(lon)) &
#                 (attr(X,'longitude') <= max(lon))
#        attr(X,'longitude') <- attr(X,'longitude')[keepx]
#        d[1] <- sum(keepx)
#        Y <- Y[keepx,,]
#      }
#  }
#  if (!is.null(lat)) {
#    if (length(lat) != 2)
#      warning("EOF: argument 'lat' must be a range") else { 
#        # Select a subset of the latitudes
#        keepy <- (attr(X,'latitude') >= min(lat)) &
#                 (attr(X,'latitude') <= max(lat))
#        attr(X,'latitude') <- attr(X,'latitude')[keepy]
#        d[2] <- sum(keepy)
#        Y <- Y[,keepy,]
#      }
#  }
#  d -> attr(X,'dimensions')
#  
#  dim(Y) <- c(d[1]*d[2],d[3])
  
  # Apply geographical weighting to account for different grid area at
  # different latitudes:
  if (verbose) print('svd:')
  stdv <- sd(c(Y),na.rm=TRUE)  # Account for mixed fields with different
                               # magnitudes...

  #print(d); print(dim(Y)); str(attr(X,'latitude'))
  Wght <-matrix(nrow=d[1],ncol=d[2])
  for (i in 1:d[1])  Wght[i,]<-sqrt(abs(cos(pi*attr(X,'latitude')/180)))
  #plot(attr(X,'latitude'),colMeans(Wght),type="l"); stop("HERE")
  dim(Wght) <- c(d[1]*d[2])
  #print(length(Wght)); print(dim(Y)); print(d[3])
  #for (it in 1:d[3]) Y[,it] <- (Wght/stdv)*Y[,it]
  
  # Exclude the missing values 'NA' and grid points with sd == 0 for all times:
  sd0 <- apply(Y,2,sd,na.rm=TRUE)
  nf <- apply(Y,2,SF)
  y <- Y[,(sd0>0.0) & (nf > 0)]
  browser()
  # Exclude the time slices with missing values:
  skip <- apply(y,1,SF); npts <- dim(y)[2]
  y <- y[skip == npts,]
  
  # Remove the mean value - center the analysis:
  if (verbose) print('center the data')
  ave <- rowMeans(y)
  #print(c(length(ave),dim(y)))
  y <- y - ave
  npca <- min(dim(y)) 
  
  # Apply the SVD decomposition: see e.g. Strang (1988) or Press et al. (1989)
  SVD <- svd(y,nu=min(c(20,npca)),nv=min(c(20,npca)))

  #print("---"); print(dim(SVD$u)); print(dim(SVD$v)); print("---")

  autocor <- 0
  if (verbose) print(paste("Find max autocorr in the grid boxes."))
  #print(dim(y))
  for (i in 1:npts) {
  vec <- as.vector(y[,i])
  i.bad <- is.na(vec)
    if (sum(i.bad) == 0) {
      ar1 <- acf(vec[],plot=FALSE)
      autocor <- max(c(autocor,ar1$acf[2,1,1]),na.rm=TRUE)
    }
  }

  n <- min(n,length(SVD$d))
  
  #a <- y[,1]; dim(a) <- c(d[1],d[2]); image(a)
  eof <- zoo(SVD$v[,1:n],order.by = dates)
  invert <- apply(SVD$u[,1:n],2,mean) < 0

  # Some data points may have been excluded due to missing values.
  # Need to insert the results for valid data onto the original grid.
  pattern <- matrix(rep(NA,d[1]*d[2]*n),d[1]*d[2],n)
  pattern[skip == npts,] <- SVD$u[,1:n]

  Ave <- rep(NA,d[1]*d[2])
  Ave[skip == npts] <- ave

  if (area.mean.expl) {
    ave <- ave + mean(A,na.rm=TRUE)
    A <- A - mean(A,na.rm=TRUE)
    s <- sd(A,na.rm=TRUE)
    A <- A/s
    n <- n + 1
    #print(dim(pattern))
    pattern <- cbind(rep(1,d[1]*d[2]),pattern)
    SVD$d <- c(d[1]*d[2]*s,SVD$d)
    #print(dim(pattern))
    eof <- merge(zoo(A,order.by=index(eof)),eof)
  }
  #print(dim(pattern)); print(d); print(n)

  # Make all the EOF vectors havine the same sense rather than
  # being random:
  pattern[,invert] <- -pattern[,invert]
  eof[,invert] <- -eof[,invert]
  
  dim(pattern) <- c(d[1],d[2],n)
  dim(Ave) <- c(d[1],d[2])
  eof <- attrcp(X,eof)
  #nattr <- softattr(X)
  #for (i in 1:length(nattr))
  #  attr(eof,nattr[i]) <- attr(X,nattr[i])

  names(eof) <- paste("X.",1:n,sep="")
  attr(eof,'pattern') <- pattern
  attr(eof,'dimensions') <- d
  attr(eof,'mean') <- Ave
  attr(eof,'max.autocor') <- autocor
  attr(eof,'eigenvalues') <- SVD$d[1:n]
  attr(eof,'sum.eigenv') <- sum(SVD$d)
  attr(eof,'tot.var') <- sum(SVD$d^2)
  attr(eof,'history') <- history.stamp(X)
  attr(eof,'aspect') <- 'anomaly'
  if (area.mean.expl) attr(eof,'area.mean.expl') <- TRUE else
                      attr(eof,'area.mean.expl') <- FALSE
  class(eof) <- c("eof",cls)
  #str(eof)
  return(eof)
}


EOF.comb <- function(X,it=NULL,is=NULL,n=20,
                     area.mean.expl=FALSE,verbose=FALSE) {

  iv <- function(x) return(sum(is.finite(x)))
  
  n.app <- attr(X,'n.apps')
  if (verbose) print(paste("EOF.comb: ",n.app,"additional field(s)"))

  # Extract the data into an ordinary matrix, remove the respective
  # mean values from the different data sets, and combine the data into
  # one matrix. Also keep track of dates, but this is monthly data, and
  # the day is used to reflect the different data set, where day==1 for
  # the original (first field), day==2 for the first field appended,
  # day ==3, for the second, and so on.

  ## AM 11-11-2013 added lines begin

  #print('subset')
  if (!is.null(is) | !is.null(it))
    X <- subset(X,it=it,is=is)
  #else if (!is.null(it))
  #  X <- subset(X,it=it,is=is)
  ## AM 11-11-2013 added lines end
  #print(dim(X))

  #print('soutpole-to-northpole')
  X <- sp2np(X)
  YYt <- t(coredata(X)); 
  clim <- rowMeans(YYt,na.rm=TRUE);
  YYt <- YYt - clim
  YY <- t(YYt)
  d <- attr(X,'dimensions')
  
  #print('time house keeping')
  t <- index(X)
  datetype <- class(t)
  if (datetype=="Date") {
    fakedates <- paste(format(t,'%Y-%m'),'-01',sep='')
    realdates <- paste(format(t,'%Y-%m'),'-01',sep='')
  } else
  if (datetype=="numeric") {
    fakedates <- paste(t,'-01-01',sep='')
    realdates <- paste(t,'-01-01',sep='')
  }

  #print(realdates)
  # Keep track of the different fields:
  if (is.null(attr(X,'source'))) attr(X,'source') <- "0"
  id.t <- rep(attr(X,'source'),length(index(X)))
  ID.t <- attr(X,'source') 
  endsofar <- max(as.numeric(format(as.Date(fakedates),'%Y')))

  for (i in 1:n.app) {
    if (verbose) print(paste("Additional field",i,endsofar))
    YYY <- attr(X,paste('appendix.',i,sep=""))

    ttt <- index(YYY)
    ##print(class(ttt)); print(ttt[1:5])
    if (inherits(ttt,'Date')) {
        year <- as.numeric(format(ttt,'%Y')) 
        month <- format(ttt,'%m')
    } else 
    if (inherits(ttt,'numeric')) {
        year <- ttt
        month <- rep(1,length(year))
    }
    yearf <- year - min(year) + endsofar + 10
    #print(year)
    #str(YYYt)
    d <- rbind(d,attr(YYY,'dimensions'))
    Zt <- t(coredata(YYY))
    clim <- rowMeans(Zt,na.rm=TRUE)
    #str(clim)
    eval(parse(text=paste('clim.',i,' <- clim',sep='')))
    Zt <- Zt  - clim
    #print(paste("Temporal-spatial mean values:",
    #            round(mean(YY,na.rm=TRUE),3),
    #            round(mean(clim,na.rm=TRUE),3)))
    #print(dim(YY)); print(dim(t(Zt)))
    YY <- rbind(YY,t(Zt))
    #print(length(index(YYY)))

    # Keep track of the different fields:
    if (is.null(attr(YYY,'source'))) attr(YYY,'source') <- as.character(i)
    src <- paste(attr(YYY,'source'),i,sep="+")
    id.t <- c(id.t,rep(src,length(index(YYY))))
    ID.t <- c(ID.t,src)
    fakedates <- c(fakedates,paste(yearf,"-",month,'-01',sep=''))
    #print(fakedates)
    realdates <- c(realdates,paste(year,"-",month,'-01',sep=''))
    endsofar <- max(as.numeric(format(as.Date(fakedates),'%Y')))
    #print('YY:'); print(dim(YY)); print("d:"); print(d)
  }
  #print(fakedates)

  #print(dates)
  
  # Synthetise a new object with combined data that looks like a
  # field object, and then call the ordinary EOF method:

  Y <- zoo(YY,order.by=as.Date(fakedates))
  #plot(rowMeans(YY,na.rm=TRUE),type="l")

  # Discard time slices with no valid data, e.g. DJF in the beginning of the record
  ngood <- apply(coredata(Y),1,iv)
  realdates <- realdates[ngood>0]
  fakedates <- fakedates[ngood>0]
  id.t <- id.t[ngood>0]
  Y <- Y[ngood>0,]
  Y <- attrcp(X,Y)
  class(Y) <- class(X)[-1]
  #print(class(Y)); print(index(Y)[1:24])
  
  attr(Y,'dimensions') <- c(d[1,1],d[1,2],sum(ngood>0))
  #print(dim(Y)); print(attr(Y,'dimensions'))
  #browser()

  eof <- EOF.field(Y,it=it,is=is,n=n,
                   area.mean.expl=area.mean.expl,verbose=verbose)

#  print("Computed the eofs...")
  # After the EOF, the results must be reorganised to reflect the different
  # data sets.
  ## browser()
  ceof <- eof
  ii <- is.element(id.t,ID.t[1])
  if (verbose) {print("Check:"); print(sum(ii)); print(ID.t); print(table(id.t))
                print(realdates[ii]); print(dim(eof))}
  ceof <- zoo(eof[ii,],order.by=as.Date(realdates[ii])) 
  ceof <- attrcp(eof,ceof)
  dim(clim) <- attr(X,'dimensions')[1:2]
  attr(ceof,'mean') <- clim
  attr(ceof,'dimensions') <- attr(X,'dimensions')
  
  for (i in 1:n.app) {
    jj <- is.element(id.t,ID.t[i+1])
    if (verbose) print(paste(ID.t[i+1],' -> appendix.',i,' data points=',sum(jj),sep=''))
    #cline <- paste("Z <- attr(X,'appendix.",i,"')",sep="")
    #print(cline)
    #eval(parse(text=cline))
    #print("Z:"); print(names(attributes(Z)))
    z <- zoo(eof[jj,],order.by=as.Date(realdates[jj]))
#    eval(parse(text=paste("XXX <- attr(X,'appendix.",i,"')",sep="")))
#    z <- zoo(eof[jj,],order.by=index(XXX))
#    rownames(z) <- as.Date(realdates[jj])
    eval(parse(text=paste("yyy <- attr(X,'appendix.",i,"')",sep="")))
    z <- attrcp(yyy,z)
    # attr(z,'clim') <-  eval(parse(text=paste('clim.',i,sep=""))) # REB attr 'mean'
    attr(ceof,paste('appendix.',i,sep="")) <- z
  }
  attr(ceof,'n.apps') <- n.app
  attr(ceof,'history') <- history.stamp(X)
  attr(ceof,'aspect') <- 'anomaly'
  class(ceof) <- c("eof",class(X))
  invisible(ceof)
}



eof2field <- function(x,it=NULL,is=NULL,anomaly=FALSE) {
  #print("HERE"); print(lon); print(lat)
  greenwich <- attr(x,'greenwich')
#  if (!is.null(lon)) lon.rng <- range(lon) else lon.rng <- NULL
#  if (!is.null(lat)) lat.rng <- range(lat) else lat.rng <- NULL
  if ( !is.null(it) | !is.null(is) ) 
    eof <- subset(x,it=it,is=is) else
#  if (!is.null(is))
#    eof <- subset(x,is=is) else
    eof <- x
  #print(c(greenwich,attr(eof,'greenwich')))
                                        # REB 04.12.13 comment below 
  U <- attr(eof,'pattern')
  d <- dim(U); dim(U) <- c(d[1]*d[2],d[3])
  W <- attr(eof,'eigenvalues')
  V <- coredata(eof)
  y <-U %*% diag(W) %*% t(V)

  if (!anomaly) y <- y + c(attr(eof,'mean'))
  y <- t(y)
  y <- as.field.default(y,index(eof),
                        lon=attr(eof,'longitude'),lat=attr(eof,'latitude'),
                param=attr(eof,'variable'),unit=attr(eof,'unit'),
                longname=attr(eof,'longname'),src=attr(eof,'source'),
                url=attr(eof,'url'),reference=attr(eof,'reference'),
                info=attr(eof,'info'),calendar=attr(eof,'calendar'),
                greenwich=attr(eof,'greenwich'))
  if (!is.null(lon) | !is.null(lat)) {
    if (is.null(lon)) lon <- c(-180,180)
    if (is.null(lat)) lat <- c(-90,90)
  }
  attr(y,'history') <- history.stamp(x)
  if (anomaly) attr(y,'aspect') <- 'anomaly' else
               attr(y,'aspect') <- 'original'
  class(y) <- class(eof)[-1]
  
  invisible(y)
}


PCA <- function(X,...) UseMethod("PCA")

PCA.default <- function(X,...) {
  stop("Don't know how to handle objects other than station")
}

PCA.station <- function(X,neofs=20,na.action='fill',verbose=FALSE) {

  if (na.action=='fill') {
    if (verbose) print('Fill missing data gaps')
    # Use interpolation to till missing data gaps. OK for small glitches.
    d <- dim(X)
    for (i in 1:d[2]) {
      x <- coredata(X[,i]); ok <- is.finite(x)
      X[,i] <- approx(y=x[ok],x=(1:d[1])[ok],xout=1:d[1])$y
    }
  }
  
  # Re-order dimensions: space x time
  x <- t(coredata(X))
  D <- dim(x)
  neofs <- min(neofs,D[1])
  ok.time <- is.finite(colMeans(x))
  z <- x[,ok.time]
  ok.site <- is.finite(rowMeans(z))
  z <- z[ok.site,]
  if (verbose) print(paste('Extract',sum(ok.site),'stations',
                             'and',sum(ok.time),'time steps',
                             'dimension=',dim(z)[1],'x',dim(z)[2]))
  X.clim <- rowMeans(x,na.rm=TRUE)
  pca <- svd(z - rowMeans(z))
  #str(pca); print(neofs)

  autocor <- 0
  #print(paste("Find max autocorr in the grid boxes."))
  #print(dim(y))
  for (i in 1:dim(z)[2]) {
    vec <- as.vector(coredata(z[,i]))
     ar1 <- acf(vec[],plot=FALSE)
     autocor <- max(c(autocor,ar1$acf[2,1,1]),na.rm=TRUE)
   }  

  # Recover the original matrix size, and insert the results
  # where there was valid data
  U <- matrix(rep(NA,D[1]*neofs),D[1],neofs)
  U[ok.site,] <- pca$u[,1:neofs]
  V <- matrix(rep(NA,D[2]*neofs),D[2],neofs)
  #print(D); print(dim(V)); print(dim(pca$v))
  #print(dim(V[ok.time,])); print(dim(pca$v[1:neofs,]))
  V[ok.time,] <- pca$v[,1:neofs]
  y <- zoo(V,order.by=index(X))
  names(y) <- paste("X.",1:neofs,sep="")

  invert <- apply(U,2,mean) < 0
  U[,invert] <- -U[,invert]
  y[,invert] <- -y[,invert]
  
  y <- attrcp(X,y)
  #nattr <- softattr(X)
  #for (i in 1:length(nattr))
  #  attr(y,nattr[i]) <- attr(X,nattr[i])
  attr(y,'pattern') <- U
  attr(y,'dimensions') <- D
  attr(y,'mean') <- X.clim
  attr(y,'max.autocor') <- autocor
  attr(y,'eigenvalues') <- pca$d[1:neofs]
  attr(y,'sum.eigenv') <- sum(pca$d)
  attr(y,'tot.var') <- sum(pca$d^2)
  attr(y,'aspect') <- 'anomaly'
  attr(y,'history') <- history.stamp(X)
  class(y) <- c("pca",class(X))
  invisible(y)
}

# Transfer PCA back to station data
pca2station <- function(X,lon=NULL,lat=NULL,anomaly=FALSE) {
  stopifnot(!missing(X), inherits(X,"pca"))
  if (inherits(X,'ds')) class(X) <- class(X)[-1]
  #print('pca2station')
  
  pca <- X
  cls <- class(pca)
  U <- attr(pca,'pattern')
  d <- dim(U)
  W <- attr(pca,'eigenvalues')
  V <- coredata(pca)
  V[!is.finite(V)] <- 0
  #str(U); str(W); str(V)
  x <-U %*% diag(W) %*% t(V)
  #str(x)
  
  if (!anomaly)
    x <- x + c(attr(pca,'mean'))
  x <- zoo(t(x),order.by=index(pca))
  
  if (anomaly) attr(x,'aspect') <- 'anomaly' else
               attr(x,'aspect') <- 'original'
#  attr(x,'call') <- match.call()
#  class(x) <- class(pca)[-1]

  names(x) <- attr(pca,'location') # AM 30.07.2013 added

  x <- attrcp(pca,x)

  # REB 2014-10-27: if the object is DS-results, then look for
  # cross-validation
  if (!is.null(attr(pca,'evaluation'))) {
    cval <- attr(pca,'evaluation')
    d.cval <- dim(cval)
    V.x <- coredata(cval)
    # The evaluation data are stored as the original calibration
    # predictand followed by the prediction, i.e. station 1,1,2,2,3,3
    ii1 <- seq(1,d[2]-1,by=2); ii2 <- seq(2,d[2],by=2)
    # Recover the staiton data from the original data x and the cross-validation prediction z
    # seperately using the same spatial PCA pattern and eiganvalues:
    x.cvalx <-U %*% diag(W) %*% t(V.x[,ii1])
    x.cvalz <-U %*% diag(W) %*% t(V.x[,ii2])
    # Combine the two together and then sort so that the prediction of the first station follows the observation
    # from the first station:
    ii <- order(seq(1,d[2],by=1),seq(1,d[2],by=1)+0.5)
    browser()
    x.cval <- rbind(x.cvalx,x.cvalz)[ii]
    mpca <- c(attr(pca,'mean'))
    jj <- order(c(1:length(mpca),1:length(mpca)+0.5))
    #browser()
    if (!anomaly) x.cval <- x.cval + rep(mpca,2)[jj]
    if (anomaly) attr(x.cval,'aspect') <- 'anomaly' else
                 attr(x.cval,'aspect') <- 'original'
    attr(x,'evaluation') <- zoo(t(x.cval),order.by=index(cval)) 
  }
  
  #nattr <- softattr(pca)
  #for (i in 1:length(nattr))
  #  attr(x,nattr[i]) <- attr(pca,nattr[i])
  # Add meta data as attributes:
  attr(x,'longitude') <- attr(pca,'longitude')
  attr(x,'latitude') <- attr(pca,'latitude')
  attr(x,'history') <- history.stamp(pca)
  class(x) <- cls[-1]
  #print("HERE")
  invisible(x)
}
