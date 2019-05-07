# R.E .Benestad

# I'm unsure about the use of 'as' functions - perhaps just use 'monthly.station' or 'aggregate'?
# as.Date(,format='%m/%d/%Y'), weekdays()


# Time conversion tools:



aggregate.station <- function(x,by,FUN = 'mean', na.rm=TRUE, ...,
                              regular = NULL, frequency = NULL) {

  args <- list(...)
  ix0 <- grep('threshold',names(args))
  iv0 <- grep('verbose',names(args))
  if (length(ix0)>0) threshold <- args[[ix0]] else threshold <- 1
  if (length(iv0)>0) verbose <- args[[iv0]] else verbose <- FALSE

  if (verbose) {print('aggregate.station'); print(names(args)); print(threshold)}
  class(x) -> cls
  #print(deparse(substitute(by)))
  class(x) <- "zoo"
  
# if (by=='yearmon') {
#    yyyymm <- format(index(Parea.merra),'%Y-%m-%d')
#    by <- yyyymm
#  }
  
  if ( (sum(is.element(names(formals(FUN)),'na.rm')==1)) |
       (sum(is.element(FUN,c('mean','min','max','sum','quantile')))>0 ) )
    y <- aggregate(x, by, FUN, na.rm=na.rm, ...,
                   regular = regular, frequency = frequency) else
    y <- aggregate(x, by, FUN, ..., regular = regular, frequency = frequency)

 # if (inherits(by[1],'character')) index(y) <- as.Date(index(y))
    
  if (class(index(y))=="Date") {
  dy <- day(y); mo <- month(y); yr <- year(y)
    if (dy[2] - dy[1] > 0) cls[length(cls) - 1] <- "day" else
    if (mo[2] - mo[1] == 1) cls[length(cls) - 1] <- "month" else
    if (mo[2] - mo[1] == 3) cls[length(cls) - 1] <- "season" else
    if (yr[2] - yr[1] > 0) cls[length(cls) - 1] <- "annual"
   } else
  if (class(index(y))=="yearmon") cls[length(cls) - 1] <- "month" else
  if (class(index(y))=="yearqtr") cls[length(cls) - 1] <- "qtr" else
  if (class(index(y))=="numeric") cls[length(cls) - 1] <- "annual" else
  if (class(index(y))=="character") cls[length(cls) - 1] <- "annual"
  if ( (length(index(y)) <= 12) & (class(index(y))=="numeric") & 
       (min(index(y)) >= 1) & (max(index(y)) <= 12) ) cls[length(cls) - 1] <- "seasonalcycle"
  class(y) <- cls
  y <- attrcp(x,y)

   #print(FUN)
  if (substr(FUN,1,4)=="count")  {
    #print("Count")
    attr(y,'unit') <- paste("counts | X >",threshold," * ",attr(x,'unit'))
  } else if (substr(FUN,1,4)=="freq") {
    #print("Frequency")
    attr(y,'variable') <- 'f'
    attr(y,'unit') <- paste("frequency | X >",threshold," * ",attr(x,'unit'))
  } else if (FUN=="wetfreq") {
    #print("Wet-day frequency")
    attr(y,'variable')[] <- 'f[w]'
    attr(y,'longname')[] <- 'Wet-day frequency'
    attr(y,'unit')[] <- paste("frequency | X >",threshold," * ",attr(x,'unit'))
  } else if (FUN=="wetmean") {
    #print("Wet-day mean")
    attr(y,'variable')[] <- 'mu'
    attr(y,'longname')[] <- 'Wet-day mean precipitation'
    attr(y,'unit')[] <- 'mm/day'
    n <- aggregate(x,by,FUN='nv', ...,
                   regular = regular, frequency = frequency)
    std.err <- 2*coredata(y)/sqrt(coredata(n)-1)
    attributes(std.err) <- NULL
    dim(std.err) <- dim(y)
    attr(y,'standard.error') <- zoo(std.err,order.by=index(y))
  } else if (FUN=="mean") {
    #print("Mean")
    sigma <- aggregate(x, by, FUN='sd', na.rm=na.rm, ...,
                       regular = regular, frequency = frequency)
    n <- aggregate(x, by, FUN='nv', na.rm=na.rm, ...,
                   regular = regular, frequency = frequency)
    #n <- count(x,threshold=threshold)
    std.err <- 2*coredata(sigma)/sqrt(coredata(n)-1)
    attributes(std.err) <- NULL
    dim(std.err) <- dim(y)
    #Finite population correction factor
    #http://en.wikipedia.org/wiki/Standard_error#Standard_error_of_the_mean
    #N <- length(n)
    #fpc <- sqrt((N-n)/(N-1))
    #std.err <- fpc*std.err
    attr(y,'standard.error') <- zoo(std.err,order.by=index(y))
  } else if (FUN=="HDD") {
    attr(y,'variable') <- 'HDD'
    attr(y,'unit') <- 'degree-days'
  } else if (FUN=="CDD") {
    attr(y,'variable') <- 'CDD'
    attr(y,'unit') <- 'degree-days'
  } else if (FUN=="GDD") {
    attr(y,'variable') <- 'GDD'
    attr(y,'unit') <- 'degree-days'
  } else attr(y,'unit') <- attr(x,'unit')

  attr(y,'history') <- history.stamp(y)
  return(y)
}

aggregate.comb <- function(x,by,FUN = 'mean', ...,
                              regular = NULL, frequency = NULL) {
  # Also need to apply the aggregation to the appended fields
  #if (verbose) print("aggregate.comb")
  #print(class(x))

  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB140314
  if (deparse(substitute(by))=="year")
    by <- as.Date(strptime(paste(year(x),1,1,sep='-'),'%Y-%m-%d'))
  
  x <- aggregate.field(x,by=by,FUN=FUN, ...,
                       regular = regular, frequency = frequency)
  n <- attr(x,'n.apps')

  print("appended fields...")
  for (i in 1:n) {
    eval(parse(text=paste("z <- attr(x,'appendix.",i,"')",sep="")))
    #print(class(z)); print(by); print(index(z))
    z <- aggregate(z,by=by,FUN=FUN, ...,
                   regular = regular, frequency = frequency)
    print("update appended field")
    eval(parse(text=paste("z -> attr(x,'appendix.",i,"')",sep="")))
  }
  return(z)
}

aggregate.field <- function(x,by,FUN = 'mean', ...,
                              regular = NULL, frequency = NULL) {

  args <- list(...)
  ix0 <- grep('threshold',names(args))
  iv0 <- grep('verbose',names(args))
  if (length(ix0)>0) threshold <- args[[ix0]] else threshold <- 0
  if (length(iv0)>0) verbose <- args[[iv0]] else verbose <- FALSE
  
  if (verbose) {print('aggregate.station'); print(names(args)); print(threshold)}
  #verbose <- TRUE; str(...)
  #if (verbose) print("aggregate.field")
  class(x) -> cls
  #print(class(index(x)))
  class(x) <- "zoo"
  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB140314
  #str(X); plot(X)
  #print("here...")
  #print(class(by))
  
  if (!is.list(by)) {
  # Temporal aggregation:
    #print("HERE")
    #print(deparse(substitute(by)))
    #clsy2 <- switch(deparse(substitute(by)),
    #                     "as.yearmon"="month",
    #                     "as.yearqtr"="quarter",
    #                     "as.annual"="annual",
    #                     "year"="annual",
    #                     "by" = "by")
    #if (is.null(clsy2)) clsy2 <- deparse(substitute(by))
    #print(clsy2)
    #if (deparse(substitute(by))[1]=="year") {
    #  ## KMP 2017-05-07: annual mean should have year as index, not date
    #  #by <- as.Date(strptime(paste(year(x),1,1,sep='-'),'%Y-%m-%d'))
    #  by <- year(x)
    #  index(x) <- year(x)
    #}
    ## REB - 'what do the following lines do?'year' changed to 'month' in the if-statement
    #if (deparse(substitute(by))[1]=="month") {
    #  print('fixed bug')
    #  by <- month(x)
    #  x <- as.data.frame(x)
    #  index(x) <- month(x)
    #}
    ## BER
    #print(deparse(substitute(by)))
    #print(class(x))
    #print(class(index(x)))
    #print('aggregate')
    ## y <- aggregate(x, by, match.fun(FUN), ...) ## AM quick fix replaced by
    y <- aggregate(x, by, FUN, ...)
    class(x) <- cls; 
    
    if (class(index(y))=="Date") {
      dy <- day(y); mo <- month(y); yr <- year(y)
      if (dy[2] - dy[1] > 0) cls[length(cls) - 1] <- "day" else
        if (mo[2] - mo[1] == 1) cls[length(cls) - 1] <- "month" else
          if (mo[2] - mo[1] == 3) cls[length(cls) - 1] <- "season" else
            if (yr[2] - yr[1] > 0) cls[length(cls) - 1] <- "annual"
    } else
      if (class(index(y))=="yearmon") cls[length(cls) - 1] <- "month" else
        if (class(index(y))=="yearqtr") cls[length(cls) - 1] <- "qtr" else
          if (class(index(y))=="numeric") cls[length(cls) - 1] <- "annual" else
            if (class(index(y))=="character") cls[length(cls) - 1] <- "annual"
    if ( (length(index(y)) <= 12) & (class(index(y))=="numeric") & 
         (min(index(y)) >= 1) & (max(index(y)) <= 12) ) cls[length(cls) - 1] <- "seasonalcycle"
    class(y) <- cls
    #class(y)[2] <- clsy2
    
    y <- attrcp(x,y)
    #nattr <- softattr(x)
    #for (i in 1:length(nattr))
    #  attr(y,nattr[i]) <- attr(x,nattr[i])
    attr(y,'dimensions') <- c(attr(x,'dimensions')[-3],length(index(y)))
    #attr(y,'history') <- c(attr(x,'history'),'aggregate')
    #attr(y,'date-stamp') <- date()
    #attr(y,'call') <- match.call()
    attr(y,'history') <- history.stamp(x)
    attr(y,'dimnames') <- NULL
    return(y)
  } else {
    #print("there")
  # spatial aggregation
    Lon <- by[[1]]; Lat=by[[2]]
    dx <- diff(Lon)[1]; dy <- diff(Lat)[1]
    t <- index(x)
    d <- attr(x,'dimensions')
    lon <- dx*floor(attr(x,'longitude')/dx)
    lat <- dy*floor(attr(x,'latitude')/dy)
    xy <- paste(rep(lon,d[2]),sort(rep(lat,d[1])),sep="/")
#    xy <- paste(rep(lat,d[1]),sort(rep(lon,d[2])),sep="/")
#    xy <- paste(sort(rep(lon,d[2])),rep(lat,d[1]),sep="/")
    lon <- as.numeric(rownames(table(lon)))
    lat <- as.numeric(rownames(table(lat)))
    D <-  c(length(lon),length(lat),length(t))
    print(paste("spatial aggregation:",d[1],"x",d[2]," ->",D[1],"x",D[2]))
    Z <- t(coredata(x)); attributes(Z) <- NULL
    dim(Z) <- c(d[1]*d[2],d[3]); rownames(Z) <- xy   
    ## z0 <- aggregate(Z,by=list(xy), match.fun(FUN),simplify=TRUE) ## AM 14-04-2015 replaced by
    z0 <- aggregate(Z,by=list(xy), FUN,simplify=TRUE)

    # The aggregate function rearranges the order of lon-lat:
    lonlat <- z0$Group.1
    ll <- strsplit(lonlat,split='/')
    lxy <- rep(NA,length(ll))
    for (i in 1:length(ll))
      lxy[i] <- as.numeric(ll[[i]][1]) + as.numeric(ll[[i]][2])*100
    #print(xy[1:100]); print(lonlat[1:100]); print(lonlat[order(lxy)][1:100])
    z <- as.matrix(z0[,2:(length(t)+1)])
    dim(z) <- c(D[1]*D[2],D[3])
    z <- z[order(lxy),]
    y <- zoo(t(as.matrix(z)),order.by=t)
    y <- attrcp(x,y)
    #nattr <- softattr(x)
    #for (i in 1:length(nattr))
    #  attr(y,nattr[i]) <- attr(x,nattr[i])
    attr(y,'dimensions') <- D
    attr(y,'latitude') <- lat
    attr(y,'longitude') <-lon
    #attr(y,'history') <- c(attr(x,'history'),'aggrgate')
    #attr(y,'date-stamp') <- date()
    #attr(y,'call') <- match.call()
    attr(y,'history') <- history.stamp(x)
    attr(y,'dimnames') <- NULL
    class(y) <- cls
    return(y)
    }
}


aggregate.area <- function(x,is=NULL,it=NULL,FUN='sum',
                           na.rm=TRUE,smallx=FALSE,verbose=FALSE,
                           a= 6378, x0=NULL) {
  # Estimate the area-aggregated values, e.g. the global mean (default)
  if (verbose) print(paste("aggregate.area",FUN))
  if (verbose) {
    if (FUN=='sum') print(rowSums(coredata(x),na.rm=TRUE)) else
                    print(rowMeans(coredata(x),na.rm=TRUE))
  }
  if (inherits(x,'eof')) {
    if (verbose) print('aggregate.area for EOF')
    y <- as.pattern(x)
    ya <- aggregate.area(y,is=is,FUN=FUN,na.rm=na.rm,smallx=smallx,verbose=verbose,a=a,x0=x0)
    if (verbose) {print(length(ya)); print(length(attr(x,'eigenvalues'))); print(t(dim(coredata(x))))}
    z <- apply(diag(ya*attr(x,'eigenvalues')) %*% t(coredata(x)),2,FUN='sum')
    if (is.zoo(x)) z <- zoo(x=z,order.by=index(x))
    attr(z,'history') <- history.stamp(x)
    return(z)
  }
  x <- subset(x,is=is,it=it,verbose=verbose)
  if ( (verbose) & (!is.null(is) | !is.null(it)) ) {
    if (FUN=='sum') print(rowSums(coredata(x),na.rm=TRUE)) else
                    print(rowMeans(coredata(x),na.rm=TRUE))
  }
  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB140314
  if (!is.null(attr(x,'dimensions'))) d <- attr(x,'dimensions') else d <- c(dim(x),1)
  if (verbose) print(paste('dimensions',paste(d,collapse='-')))
  #image(attr(x,'longitude'),attr(x,'latitude'),area)
  #print(c(length(colSums(area)),length(attr(x,'latitude')),sum(colSums(area))))
  #lon <- rep(lon(x),d[2])
  ## KMP 2017-10-18: this doesn't look right. length(d) is 3 for all fields. 
  #if (inherits(x,'pattern') | length(d)==3) {
  if (inherits(x,'pattern')) {
    if (verbose) print('need to make the pattern look like field')
    dim(x) <- c(d[1]*d[2],d[3])
    x <- t(x)
  }
  
  srtlat <- order(rep(lat(x),d[1]))
  dY <- a*diff(pi*lat(x)/180)[1]
  dtheta <- diff(pi*lon(x)/180)[1]
  ## The first assumes a global field and the second is for a limited longitude range
  #if (diff(range(lon(x)))> 350) aweights <- rep(dY * 2*pi/d[1] * a*cos(pi*lat(x)/180),d[1])[srtlat] else
  aweights <- rep(dY * dtheta * a*cos(pi*lat(x)/180),d[1])[srtlat]
  if (verbose) print(sum(aweights))
  if (FUN=='mean') {
    aweights <- aweights/sum(aweights,na.rm=TRUE)
    FUN <- 'sum'
  }
  if (verbose) print(paste('Sum of aweights should be area or 1:',round(sum(aweights))))

  
  ## REB: For sum, we also need to consider the area:
  if (FUN %in% c('sum','area','exceedance','exceedence','lessthan')) {
    if (FUN=='area') {
      ## Estimate the area of the grid boxes
      coredata(x) -> cx
      ## REB: 2018-11-13: minor fix - added 'cx[!is.finite(cx)] <- 0'
      if (is.null(x0)) {cx[is.finite(cx)] <- 1; cx[!is.finite(cx)] <- 0} else {cx[cx < x0] <- 0; cx[cx >= x0] <- 1}
      coredata(x) <- cx; rm('cx'); gc(reset=TRUE)
      FUN <- 'sum'
    } else if ( (FUN %in% c('exceedance','exceedence')) & !is.null(x0) ) {
      # Estimate the sum of grid boxes with higher value than x0
      coredata(x) -> cx
      cx[cx < x0] <- NA
      coredata(x) <- cx; rm('cx'); gc(reset=TRUE)
      FUN <- 'sum'
    } else if ( (FUN == 'lessthan') & !is.null(x0) ) {
      # Estimate the sum of grid boxes with lower value than x0
      coredata(x) -> cx
      cx[cx >= x0] <- NA
      coredata(x) <- cx; rm('cx'); gc(reset=TRUE)
      FUN <- 'sum'
    } 
    attr(x,'unit') <- paste(attr(x,'unit'),' * km^2')
  }
  
  if (smallx) {
    if (sum(!is.finite(x))==0) X <- coredata(x)%*%diag(aweights) else {
      if (verbose) print('Need to account for missing data in the area weighting')
        Aweights <- rep(aweights,length(index(x))); dim(Aweights) <- dim(x)
        print('This is incomplete - needs vchecking!')
        browser()
        Aweights[!is.finite(coredata(x))] <- NA
        Aweights <- Aweights/apply(Aweights,1,FUN='sum',na.rm=TRUE)
        X <- coredata(X)*Aweights
    }
    y <- zoo(apply(X,1,FUN,na.rm=na.rm),order.by=index(x))
  } else {
    X <-coredata(x) 
    if (d[3]==1) dim(X) <- c(1,length(X)) ## If only one map, then set the dimensions right to get a matrix.
    if (verbose) {print(dim(X)); print(length(aweights))}
    for (i in 1:d[3]) {
      ## Temporary weights to account for variable gaps of missing data
      aweights2 <- aweights
      aweights2[!is.finite(coredata(x[i,]))] <- NA
      aweights2 <- aweights2/sum(aweights2,na.rm=TRUE)
      X[i,] <- X[i,]*aweights2
    }
    y <- zoo(apply(X,1,FUN,na.rm=na.rm),order.by=index(x))
  }
  if (verbose) print(y)
  
  Y <- as.station(y,loc=paste('area',FUN,'of',src(x)),
                  param=attr(x,'variable'),
                  unit=attr(x,'unit'),
                  lon=range(lon(x)),lat=range(lat(x)),alt=NA,cntr=NA,
                  longname=paste(FUN,attr(x,'longname')),stid=NA,quality=NA,
                  src=attr(x,'source'),url=attr(x,'url'),
                  reference=attr(x,'reference'),info=attr(x,'info'),
                  method=paste(FUN,attr(x,'method')),type='area aggregate',
                  aspect=attr(x,'aspect'))
  if (verbose) attr(Y,'aweights') <- aweights
  attr(Y,'history') <- history.stamp(x)
  return(Y)
}

