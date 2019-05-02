## Computes Empirical Orthogonal Functions (EOFs)
##
## R.E. Benestad, 
## rasmus.benestad@met.no
## 
## ------------------------------------------------------------------------

EOF<-function(X,it=NULL,is=NULL,n=20,lon=NULL,lat=NULL,verbose=FALSE,anomaly=TRUE,...)
  UseMethod("EOF")

EOF.default <- function(X,it=NULL,is=NULL,n=20,lon=NULL,lat=NULL,verbose=FALSE,anomaly=TRUE,...) {
  # Verify Arguments
  if (verbose) print("EOF.default")
  stopifnot(!missing(X), is.matrix(X),inherits(X,"zoo"))
  
  if ( !zeros(inherits(X,c("comb","zoo"),which=TRUE)) )
    eof <- EOF.comb(X,it=it,is=is,n=n,anomaly=anomaly,
                    verbose=verbose) else
                      if ( !zeros(inherits(X,c("field","zoo"),which=TRUE)) )
                        eof <- EOF.field(X,it=it,is=is,n=n,
                                         anomaly=anomaly,
                                         verbose=verbose)
  attr(eof,'dimnames') <- NULL   # REB 2016-03-04
  return(eof)
}


# Apply EOF analysis to the monthly mean field values:

EOF.field <- function(X,it=NULL,is=NULL,n=20,lon=NULL,lat=NULL,verbose=FALSE,anomaly=TRUE,...) {
  
  SF <- function(x) {sum(is.finite(x))}
  
  if (verbose) print("EOF.field")
  attr(X,'dimnames') <- NULL
  stopifnot(!missing(X), is.matrix(X),
            inherits(X,c("field","zoo")))
  
  # Remove time slices with missing data:
  #  nok <- !is.finite(rowMeans(X))
  # The regridded appendix may contain some NA's if its domain exceeds that of the original field.
  # Get rid of time slices with all NAs.
  nok <- apply(X,1,nv) < 0.5*dim(X)[2]
  if (sum(nok)> 0) {
    it <- (1:length(nok))[!nok]
    if (verbose) print(paste('removing ',sum(nok),'NA time slices'))
  }
  
  x <- subset(X,it=it,is=is,verbose=verbose)
  x <- sp2np(x)
  dates <- index(x)
  if (verbose) print(dates)
  
  d <- attr(x,'dimensions')
  if ((length(d) != 3) | min(d) == 1) {
    stop(paste('EOF.field: too small data dimensions'))
  }    
  cls <- class(x)
  
  Y <- t(coredata(x))
  Y[!is.finite(Y)] <- NA
  
  # Apply geographical weighting to account for different grid area at
  # different latitudes:
  if (verbose) print('svd:')
  stdv <- sd(c(Y),na.rm=TRUE)  # Account for mixed fields with different
  # magnitudes...
  
  #print(d); print(dim(Y)); str(attr(X,'latitude'))
  Wght <- matrix(nrow=d[1],ncol=d[2])
  for (i in 1:d[1])  Wght[i,]<-sqrt(abs(cos(pi*attr(X,'latitude')/180)))
  #plot(attr(X,'latitude'),colMeans(Wght),type="l"); stop("HERE")
  dim(Wght) <- c(d[1]*d[2])
  #print(length(Wght)); print(dim(Y)); print(d[3])
  #for (it in 1:d[3]) Y[,it] <- (Wght/stdv)*Y[,it]

  # Exclude the missing values 'NA' and grid points with sd == 0 for all times:
  sd0 <- apply(as.matrix(Y),2,sd,na.rm=TRUE)
  nf <- apply(as.matrix(Y),2,SF)
  if (verbose) print(paste('Exclude the missing values/zero-sd:',
                           sum(sd0>0.0),sum(nf > 0)))
  y <- Y[,(sd0>0.0) & (nf > 0)]
  
  # Exclude the time slices with missing values:
  skip <- apply(as.matrix(y),1,SF); npts <- dim(y)[2]
  y <- as.matrix(y)[skip == npts,]

  # Remove the mean value - center the analysis:
  if (anomaly) {
    if (verbose) print('center the data')
    ave <- rowMeans(y)
    #print(c(length(ave),dim(y)))
    y <- y - ave
  } else ave <- rowMeans(y)*0
  npca <- min(dim(y)) 
  ny <- min(c(dim(y),20))

  # REB 2015-05-21
  # Apply the SVD decomposition: see e.g. Strang (1988) or Press et al. (1989)
  #SVD <- svd(y,nu=min(c(20,npca)),nv=min(c(20,npca)))
  ## KMP 23-11-2015
  ## When running DSensemble, this convergence error comes up occasionally:
  ## "Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'."
  ## For some reason rounding the matrix before performing svd helps.
  ## Another alternative is to do svd on the matrix transpose
  #y <- round(y,digits=10)
  #SVD <- svd(y,nu=min(c(ny,npca)),nv=min(c(ny,npca)))
  SVD <- try(svd(y,nu=min(c(ny,npca)),nv=min(c(ny,npca))))
  if (inherits(SVD,"try-error")) {
    if (verbose) print("svd(x) failed. try svd(t(x))")
    SVD <- try(svd(t(y),nu=min(c(ny,npca)),nv=min(c(ny,npca))))
    temp <- SVD$u
    SVD$u <- SVD$v
    SVD$v <- temp
  }
  if (inherits(SVD,"try-error") & verbose) print("both svd(x) and svd(t(x) failed.")
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
  # Make all the EOF vectors havine the same sense rather than
  # being random:
  if (verbose) print(paste("Invert EOF",(1:length(invert))[invert],collapse=' '))
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
  attr(eof,'dimnames') <- NULL   # REB 2016-03-04
  class(eof) <- c("eof",cls)
  #str(eof)
  return(eof)
}


EOF.comb <- function(X,it=NULL,is=NULL,n=20,lon=NULL,lat=NULL,verbose=FALSE,anomaly=TRUE,...) {
  
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
    X <- subset(X,it=it,is=is,verbose=verbose)
  #else if (!is.null(it))
  #  X <- subset(X,it=it,is=is)
  ## AM 11-11-2013 added lines end
  #print(dim(X))
  
  #print('soutpole-to-northpole')
  X <- sp2np(X)
  YYt <- t(coredata(X));
  clim <- rowMeans(YYt,na.rm=TRUE)
  clim.0 <- clim
  YYt <- YYt - clim
  YY <- t(YYt)
  d <- attr(X,'dimensions')
  
  ## KMP 2016-12-28: time housekeeping in EOF.comb creates problems
  ## for DS when applied to annual data. The predictand ends up with
  ## years as time index but the EOF has dates (YYYY-01-01).
  ## Do the realdates and fakedates need to be in date format?
  if (verbose) print('time house keeping')
  t <- index(X)
  datetype <- class(t)
  if (datetype=="Date") {
    fakedates <- paste(format(t,'%Y-%m'),'-01',sep='')
    realdates <- paste(format(t,'%Y-%m'),'-01',sep='')
    endsofar <- max(as.numeric(format(as.Date(fakedates),'%Y')))
  } else if (datetype=="numeric") {
    fakedates <- t#paste(t,'-01-01',sep='')#
    realdates <- t#paste(t,'-01-01',sep='')#
    endsofar <- max(t)#max(as.numeric(format(as.Date(fakedates),'%Y')))#
  }
  
  #print(realdates)
  # Keep track of the different fields:
  if (is.null(attr(X,'source'))) attr(X,'source') <- "0"
  id.t <- rep(attr(X,'source'),length(index(X)))
  ID.t <- attr(X,'source')
  
  for (i in 1:n.app) {
    if (verbose) print(paste("Additional field",i,endsofar))
    YYY <- attr(X,paste('appendix.',i,sep=""))
    
    ttt <- index(YYY)
    ##print(class(ttt)); print(ttt[1:5])
    if (inherits(ttt,'Date')) {
      year <- as.numeric(format(ttt,'%Y')) 
      month <- format(ttt,'%m')
    } else if (inherits(ttt,c('numeric','integer'))) {
      year <- ttt
      month <- rep(1,length(year))
    }
    yearf <- year - min(year) + endsofar + 10
    if (inherits(ttt,'Date')) {
      fakedates <- c(fakedates,paste(yearf,"-",month,'-01',sep=''))
      realdates <- c(realdates,paste(year,"-",month,'-01',sep=''))
      endsofar <- max(as.numeric(format(as.Date(fakedates),'%Y')))
      #print(fakedates)
    } else if (inherits(ttt,c('numeric','integer'))) {
      fakedates <- c(fakedates,yearf)
      realdates <- c(realdates,year)
      endsofar <- max(as.numeric(fakedates))
    }
    #print(year)
    #str(YYYt)
    d <- rbind(d,attr(YYY,'dimensions'))
    Zt <- t(coredata(YYY))
    clim <- rowMeans(Zt,na.rm=TRUE)
    #str(clim)
    eval(parse(text=paste('clim.',i,' <- clim',sep='')))
    Zt <- Zt - clim
    #print(paste("Temporal-spatial mean values:",
    #            round(mean(YY,na.rm=TRUE),3),
    #            round(mean(clim,na.rm=TRUE),3)))
    #print(dim(YY)); print(dim(t(Zt)))
    YY <- rbind(YY,t(Zt))
    #print(length(index(YYY)))
    
    # Keep track of the different fields:
    if (is.null(attr(YYY,'source'))) attr(YYY,'source') <- as.character(i) else
      if (is.na(attr(YYY,'source')))  attr(YYY,'source') <- as.character(i)
    src <- paste(attr(YYY,'source'),i,sep="+")
    id.t <- c(id.t,rep(src,length(index(YYY))))
    ID.t <- c(ID.t,src)
    #print('YY:'); print(dim(YY)); print("d:"); print(d)
  }
  
  # Synthetise a new object with combined data that looks like a
  # field object, and then call the ordinary EOF method:
  
  if (verbose) print("combine original and appended fields")
  if(is.character(fakedates)){
    Y <- zoo(YY,order.by=as.Date(fakedates))
  } else {
    Y <- zoo(YY,order.by=fakedates)#as.Date(fakedates))
  }
  #plot(rowMeans(YY,na.rm=TRUE),type="l")
  
  # Discard time slices with no valid data, e.g. DJF in the beginning of the record
  ngood <- apply(coredata(Y),1,nv)
  if (verbose) print(summary(ngood))
  if (verbose) {print("Check:"); print(table(id.t)); print(dim(Y))}
  
  if (verbose) print(paste('Remove missing data gaps: ngood <- apply(coredata(Y),2,nv):',ngood))
  realdates <- realdates[ngood>0]
  fakedates <- fakedates[ngood>0]
  id.t <- id.t[ngood>0]
  Y <- Y[ngood>0,]
  Y <- attrcp(X,Y)
  class(Y) <- class(X)[-1]
  
  attr(Y,'dimensions') <- c(d[1,1],d[1,2],sum(ngood>0))
  if (verbose) {print(dim(Y)); print(attr(Y,'dimensions'))}
  
  if (verbose) print('Ordinary EOF')
  eof <- EOF.field(Y,it=it,is=is,n=n,lon=lon,lat=lat,anomaly=anomaly,verbose=verbose)
  
  if (verbose) print("Computed the eofs:")
  # After the EOF, the results must be reorganised to reflect the different
  # data sets.
  ceof <- eof
  ii <- is.element(id.t,ID.t[1])
  if (verbose) {print("Check:"); print(sum(ii)); print(ID.t); print(table(id.t))
    print(realdates[ii]); print(dim(eof))}
  
  if(is.character(realdates)) {
    ceof <- zoo(eof[ii,],order.by=as.Date(realdates[ii]))
  } else {
    ceof <- zoo(eof[ii,],order.by=realdates[ii])
  }
  
  ## Finalise - set the metadata
  if (verbose) {print("Copy attributes"); print(names(attributes(eof)))}
  ceof <- attrcp(eof,ceof)
  clim <- clim.0
  dim(clim) <- attr(X,'dimensions')[1:2]
  attr(ceof,'mean') <- clim
  attr(ceof,'dimensions') <- attr(X,'dimensions')
  
  ## Set the metadata for the appended data: climatology etc.
  for (i in 1:n.app) {
    jj <- is.element(id.t,ID.t[i+1])
    if (verbose) print(paste(ID.t[i+1],' -> appendix.',i,' data points=',sum(jj),sep=''))
    if(is.character(realdates)) {
      z <- zoo(eof[jj,],order.by=as.Date(realdates[jj]))
    } else {
      z <- zoo(eof[jj,],order.by=realdates[jj])
    }
    cline1 <- paste("yyy <- attr(X,'appendix.",i,"')",sep="")
    if (verbose) print(cline1)
    eval(parse(text=cline1))
    z <- attrcp(yyy,z)
    cline2 <- paste("clim <- clim.",i,sep="")
    if (verbose) print(cline2)
    eval(parse(text=cline2))
    dim(clim) <- attr(X,'dimensions')[1:2]
    if (verbose) print('add information about mean')
    attr(z,"mean") <- clim
    attr(ceof,paste('appendix.',i,sep="")) <- z
  }
  attr(ceof,'n.apps') <- n.app
  attr(ceof,'history') <- history.stamp(X)
  attr(ceof,'aspect') <- 'anomaly'
  attr(ceof,'dimnames') <- NULL   # REB 2016-03-04
  class(ceof) <- c("eof",class(X))
  invisible(ceof)
}



eof2field <- function(x,it=NULL,is=NULL,ip=NULL,anomaly=FALSE,verbose=FALSE) {
  if (verbose) {print("eof2field"); if (!is.null(is)) print(is)}
  greenwich <- attr(x,'greenwich')
  if ( !is.null(it) | !is.null(is) ) {
    eof <- subset(x,it=it,is=is,verbose=verbose)
  } else {
    eof <- x
  }
  U <- attr(eof,'pattern')
  d <- dim(U) 
  if (verbose) {str(U); print(d)}
  dim(U) <- c(d[1]*d[2],d[3])
  W <- attr(eof,'eigenvalues')
  V <- coredata(eof)
  ### ==================================================
  ## KMP 2016-01-15: added selection of patterns (ip)
  if(is.null(ip)) {
    ip <- seq(length(W))
  } else if(any(ip %in% seq(length(W)))) {
    ip <- ip[ip>0 & ip<length(W)]
  } else {
    stop(paste("Error in input ip =",paste(ip,collaps=", ")))
  }
  ### ==================================================
  U <- U[,ip]; W <- W[ip]; V <- V[,ip]
  y <-U %*% diag(W) %*% t(V)
  
  if (!anomaly) {
    if (verbose) print('Anomalies')
    y <- y + c(attr(eof,'mean'))
  }
  y <- t(y)
  y <- as.field.default(y,index=index(eof),
                        lon=attr(eof,'longitude'),lat=attr(eof,'latitude'),
                        param=attr(eof,'variable'),unit=attr(eof,'unit'),
                        longname=attr(eof,'longname'),src=attr(eof,'source'),
                        url=attr(eof,'url'),reference=attr(eof,'reference'),
                        info=attr(eof,'info'),calendar=attr(eof,'calendar'),
                        greenwich=attr(eof,'greenwich'),verbose=verbose)
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