
as.station <- function(x,...) UseMethod("as.station")

as.station.zoo <- function(x,loc=NA,param=NA,unit=NA,lon=NA,lat=NA,alt=NA,
                          cntr=NA,longname=NA,calendar=NA,stid=NA,quality=NA,src=NA,
                          url=NA,reference=NA,info=NA, method= NA,type=NA,
                           aspect=NA,verbose=FALSE) {
  #print(c(length(X),length(index)))
  y <- zoo(x,order.by=index(x))

  if (is.null(loc)) loc <- NA
  if ((is.na(loc[1])) & !is.null(attr(x,'location')))
    loc <- attr(x,'location')
  attr(y,'location') <- loc
  if (is.null(param)) param <- NA
  if ((is.na(param[1])) & !is.null(attr(x,'variable')))
    param <- attr(x,'variable')
  attr(y,'variable') <- param
  if (is.null(unit)) unit <- NA
  if ((is.na(unit[1])) & !is.null(attr(x,'unit')))
    unit <- attr(x,'unit')
  attr(y,'unit') <- unit
  if (is.null(lon)) lon <- NA
  if ((is.na(lon[1])) & !is.null(attr(x,'longitude')))
    lon <- attr(x,'longitude')  
  attr(y,'longitude') <- lon
  if (is.null(lat)) lat <- NA
  if ((is.na(lat[1])) & !is.null(attr(x,'latitude')))
    lat <- attr(x,'latitude')
  attr(y,'latitude') <- lat
  if (is.null(alt)) alt <- NA
  if ((is.na(alt[1])) & !is.null(attr(x,'altitude')))
    alt <- attr(x,'altitude')    
  attr(y,'altitude') <- alt
  if (is.null(cntr)) cntr <- NA
  if ((is.na(cntr[1])) & !is.null(attr(x,'country')))
    cntr <- attr(x,'country')      
  attr(y,'country') <- cntr
  if (is.null(longname)) longname <- NA
  if ((is.na(longname[1])) & !is.null(attr(x,'longname')))
    longname <- attr(x,'longname')    
  attr(y,'longname') <- longname
  if (is.null(stid)) stid <- NA
  if ((is.na(stid[1])) & !is.null(attr(x,'station_id')))
    stid <- attr(x,'station_id')    
  attr(y,'station_id') <- stid
  if (is.null(quality)) quality <- NA
  if ((is.na(quality[1])) & !is.null(attr(x,'quality')))
    quality <- attr(x,'guality')    
  attr(y,'quality') <- quality
  if (is.na(calendar)) {
    if(!is.null(attr(x,'calendar'))) {
      calendar <- attr(y,'calendar')
    } else {
      calendar <- 'gregorian'
    }
  }
  attr(y,'calendar') <- calendar
  if (is.null(src)) src <- NA 
  if ((is.na(src[1])) & !is.null(attr(x,'source')))
    src <- attr(x,'source')    
  attr(y,'source') <- src
  if (is.null(url)) url <- NA 
  if ((is.na(url[1])) & !is.null(attr(x,'URL')))
    url <- attr(x,'URL')    
  attr(y,'URL') <- url
  #attr(y,'history') <- 'as.station.data.frame'
  #attr(y,'date-stamp') <- date()
  if (is.null(type)) type <- NA 
  if ((is.na(type)) & !is.null(attr(x,'type')))
    type <- attr(x,'type')    
  attr(y,'type') <- type
  if (is.null(aspect)) aspect <- NA 
  if ((is.na(aspect[1])) & !is.null(attr(x,'aspect')))
    aspect <- attr(x,'aspect')    
  attr(y,'aspect') <- aspect
  if (is.null(reference)) reference <- NA 
  if ((is.na(reference[1])) & !is.null(attr(x,'reference')))
    reference <- attr(x,'reference')    
  attr(y,'reference') <- reference
  if (is.null(info)) info <- NA 
  if ((is.na(info[1])) & !is.null(attr(x,'info')))
    info <- attr(x,'info')    
  attr(y,'info') <- info
  if (is.null(method)) method <- NA
  if ((is.na(method[1])) & !is.null(attr(x,'method')))
    method <- attr(x,'method')    
  attr(y,'method') <- method
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)

  ## Make sure that the index is a Date object:
  ##print(index(y))
  if (is.numeric(index(y)) | is.integer(index(y)) |
      (is.character(index(y)) & (max(nchar(index(y))) ==4))) {
    tdate <- paste(index(y),'-01-01',sep='')
    ##print(tdate)
    index(y) <- as.Date(tdate)
  }
  dfi <- diff(index(y))
  if (length(dfi)>0) {
      dt <- as.numeric(median(dfi))
      if (dt==1)
          tscale <- 'day'
      else if (((dt>=28) & (dt <=31)) | (dt < 0.1))
          tscale <- 'month'
      else if ( (dt>=89) & (dt <=93) )
          tscale <- 'season' else 
      if ((dt>=360) & (dt <=366))
          tscale <- 'annual'
      else tscale <- 'annual'
      class(y) <- c("station",tscale,"zoo")
  }
  else if (verbose) print("Warning - A single value has been recorded in the time index of the data")
  return(y)
}

as.station.data.frame <-  function (x, loc = NA, param = NA, unit = NA,
                                    lon = NA, lat = NA,
    alt = NA, cntr = NA, longname = NA, stid = NA, quality = NA,
    src = NA, url = NA, reference = NA, info = NA, method= NA,type=NA,aspect=NA)
{
    cnam <- names(x)
    year <- sort(rep(x[[1]], 12))
    month <- rep(1:12, length(x[[1]]))
    index <- as.Date(paste(year, month, 1, sep = "-"))
    X <- c(t(as.matrix(x)[, 2:13]))
    y <- zoo(X, order.by = index)
    attr(y, "location") <- loc
    attr(y, "variable") <- param
    attr(y, "unit") <- unit
    attr(y, "longitude") <- lon
    attr(y, "latitude") <- lat
    attr(y, "altitude") <- alt
    attr(y, "country") <- cntr
    attr(y, "longname") <- longname
    attr(y, "station_id") <- stid
    attr(y, "quality") <- quality
    attr(y, "calendar") <- "gregorian"
    attr(y, "source") <- src
    attr(y, "URL") <- url
    #attr(y, "history") <- "as.station.data.frame"
    #attr(y, "date-stamp") <- date()
    attr(y, "type") <- type
    attr(y, "aspect") <- aspect
    attr(y, "reference") <- reference
    attr(y, "info") <- info
    attr(y, "method") <- method
    #attr(y, "call") <- match.call()
    attr(y,'history') <- history.stamp(x)
    class(y) <- c("station", "month", "zoo")
    return(y)
}

as.station.ds <- function(x) {
  if (inherits(x,'pca')) {
    class(x) <- class(x)[-1]
    y <- as.station.pca(x)
  } else {
    y <- zoo(coredata(x),order.by=index(x))
    if (!is.null(attr(x,"original_data")))
      y <- attrcp(attr(x,"original_data"),y) else
      y <- attrcp(x,y)
    ##if (!is.null(attr(x,"original_data")))
      class(y) <- class(attr(x,"original_data"))## else
    ##  class(y) <-class(x)
  }
  attr(y,'history') <- history.stamp(x)
  attr(y,'method') <- attr(x,'method')
  attr(y,'info') <- attr(x,'info')
  class(y) <- c("station", class(x)[-(1:(length(class(x))-2))])
  return(y)
}

as.station.pca <- function(x,...) {
  if (inherits(x,"dsensemble")) {
    y <- as.station.dsensemble.pca(x,...)
  } else {
    y <- pca2station(x,...)
    if (!is.null(attr(x,"pca"))) {
      fit <- attr(x,'pca')
      coredata(fit) <- coredata(attr(x,'fitted_values'))
      attr(y,'fitted_values') <- fit
      attr(y,'original_data') <- attr(x,'original_data')
    }
  }
  return(y)
}


as.station.list <- function(x,verbose=FALSE) {
  if(verbose) print("as.station.ds")
#  Jan <- x$Jan + attr(x$Jan,'mean')
#  Feb <- x$Feb + attr(x$Feb,'mean')
#  Mar <- x$Mar + attr(x$Mar,'mean')
#  Apr <- x$Apr + attr(x$Apr,'mean')
#  May <- x$May + attr(x$May,'mean')
#  Jun <- x$Jun + attr(x$Jun,'mean')
#  Jul <- x$Jul + attr(x$Jul,'mean')
#  Aug <- x$Aug + attr(x$Aug,'mean')
#  Sep <- x$Sep + attr(x$Sep,'mean')
#  Oct <- x$Oct + attr(x$Oct,'mean')
#  Nov <- x$Nov + attr(x$Nov,'mean')
#  Dec <- x$Dec + attr(x$Dec,'mean')
#  cline <- "merge.zoo("  # AM/REB 2014-12-01
  cline <- "rbind("
  if (is.list(x)) {
    for (i in 1:length(x)) {
        ave <- switch(attr(x[[i]],'aspect'),
                                 'original'=0,
                                 'downscaled'=0, ## AM 2014-12-02 x[[i]] could be downscaled results
                    'anomaly'=attr(x[[i]],'mean'))
      attr(x[[i]],'type') <- 'downscaled results'
      z <- x[[i]] + ave
      eval(parse(text=paste("ds.",i," <- z",sep="")))
      cline <- paste(cline,"ds.",i,",",sep="")
    }
 
    cline <- paste(substr(cline,1,nchar(cline)-1),')-> ALL')
    #print(cline)
    eval(parse(text=cline))
    y <- zoo(coredata(ALL),order.by=index(ALL))
    #print(names(y))
  } else {
    if (inherits(x,'ds')) {
      ave <- switch(attr(x,'aspect'),
                    'original'=0,
                    'anomaly'=attr(x,'mean'))
      attr(x,'type') <- 'downscaled results'
      z <- x + ave
      class(z) <- "zoo"
      y <- zoo(coredata(z),order.by=index(z))
      x <- list(x)
    }
  }
  
                     
#  merge.zoo(Jan,Feb,Mar,Apr,May,Jun,
#            Jul,Aug,Sep,Oct,Nov,Dec)-> ALL
  
  attr(y,'location') <- attr(x[[1]],'location')
  attr(y,'variable') <- attr(x[[1]],'variable')
  attr(y,'unit') <- attr(x[[1]],'unit')
  attr(y,'longitude') <- attr(x[[1]],'longitude')
  attr(y,'latitude') <- attr(x[[1]],'latitude') 
  attr(y,'altitude') <- attr(x[[1]],'altitude')
  attr(y,'country') <- attr(x[[1]],'country')
  attr(y,'longname') <- attr(x[[1]],'longname')
  attr(y,'station_id') <- attr(x[[1]],'station_id')
  attr(y,'quality') <- attr(x[[1]],'quality')
  attr(y,'calendar') <- attr(x[[1]],'calendar')
  attr(y,'source') <- attr(x[[1]],'source')
  attr(y,'URL') <- attr(x[[1]],'URL')
  #attr(y,'history') <- c('as.station.ds',attr(x$Jan,'history'))
  #attr(y,'date-stamp') <- attr(x,'date-stamp')
  attr(y,'type') <-  attr(x,'type')
  attr(y,'aspect') <- attr(x,'aspect')
  attr(y,'reference') <- attr(x,'reference')
  attr(y,'info') <- attr(x,'info')
  #attr(y,'date') <- date()
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  if (attr(x[[i]],'type')== 'downscaled results')
      class(y) <- class(x[[1]])
  else
      class(y) <- c("station",class(x[[1]]))
  return(y)
}

as.station.field <- function(x,is=NULL,verbose=FALSE) {
  index <- index(x)
  stopifnot(!missing(x),
           !zeros(inherits(x,c("field","zoo"),which=TRUE) ))
  if (!is.null(is)) y <- regrid.field(x,is=is,verbose=verbose) else y <- x

  ## Expand the grid ...
  lon <- expand.grid(lon(y),lat(y))$Var1
  lat <- expand.grid(lon(y),lat(y))$Var2
  ns <- length(lon) 

  y <- as.station.zoo(y,loc=NA,param=rep(varid(y),ns),unit=rep(unit(y),ns),
                      lon=lon,lat=lat,
                      alt=rep(NA,ns),cntr=rep(NA,ns),
                      longname=rep(attr(y,'longname'),ns),
                      stid=rep(NA,ns),quality=rep(NA,ns),
                      src=rep(src(y),ns),
                      url=rep(attr(y,'URLs'),ns),
                      reference=rep(attr(y,'reference'),ns),
                      info=rep(attr(y,'info'),ns), method="field",
                      type=NA,aspect=NA)
  index(y) <- index
  attr(y,'history') <- history.stamp(x)
  class(y)[2] <- class(x)[2]
  if (dim(y)[2]==1) y <- subset(y,is=1)
  invisible(y)
}

as.station.spell <- function(x) {
  y <- coredata(x)
  ok1 <- is.finite(y[,1])
  above <- zoo(y[ok1,1],index(x)[ok1])
  ok2 <- is.finite(y[,2])
  below <- zoo(y[ok2,2],index(x)[ok2])
#  if (attr(x,'variable')=='t2m') suffix <- c("warm","cold") else
#                                 suffix <- c("wet","dry")
  suffix <- attr(x,'variable')
  y <- merge.zoo(above,below,suffixes=suffix,fill=0)
  #plot(y)
  y <- attrcp(x,y)
  attr(y,'unit') <- c("days","days")
  attr(y,'variable') <- suffix
  class(y) <- c('station',class(x))
  attr(y,'history') <- history.stamp(x)
  return(y)
}


as.station.eof <- function(x,ip=1:10) {
  stopifnot(!missing(x),inherits(x,'eof'))
  z <- zoo(x[,ip],order.by=index(x))
  y <- as.station.zoo(z,loc=paste('PC[',ip,']',sep=''),
                      param=varid(x),unit=unit(x),
                      lon=lon(x),lat=lat(x),alt=NA,
                      cntr=NA,longname='principal component',stid=NA,quality=NA,
                      src=src(x),url=attr(x,'url'),reference=NA,info=NA, method="field",
                      type='index',aspect='anomaly')
  attr(y,'history') <- history.stamp(x)
  class(y)[2] <- class(x)[2]
  if (dim(y)[2]==1) y <- subset(y,is=1)
  invisible(y)
}


as.station.dsensemble <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.station.dsensemble")
  if (!is.null(x$pca) & !inherits(x,"pca") & inherits(x,"eof")) {
    class(x) <- gsub("eof","pca",class(x)) ## REB 2016-12-13: to also work on gridded versions.
  }
  if (inherits(x,"pca")) {
    y <- as.station.dsensemble.pca(x,verbose=verbose,...)
  } else if (inherits(x,c("station","list"))) {
    y <- as.station.dsensemble.station(x,verbose=verbose,...)
  } else {
    print(paste('unexpected class - do not know how to handle:',class(x)))
    y <- x
  }
  if (verbose) print(class(y))
  return(y)
}

as.station.dsensemble.pca <- function(x,is=NULL,ip=NULL,verbose=FALSE,...) {
  X <- x ## quick fix
  if (verbose) print('as.station.dsensemble.pca')
  ## REB: need to remove the EOF object if it is present:
  if (!is.null(X$eof)) X$eof <- NULL
  if (inherits(X,"station")) return(X)
  stopifnot(inherits(X,"dsensemble") & inherits(X,"pca"))
  if (inherits(X,"station")) {
      invisible(X)
  } else {
    #if (is.null(is)) is <- 1:length(loc(X$pca)) 
    if (verbose) print('Extract the results model-wise')
    ## Find the size of the PC matrices representing model projections
    d <- apply(sapply(X[3:length(X)],dim),1,min)
    ## The PCs from the list are extracted into the matrix V 
    V <- array(unlist(lapply( X[3:length(X)],
      function(x) coredata(x[1:d[1],1:d[2]]))),dim=c(d,length(X)-2))
    if (verbose) print(paste('dim V=',paste(dim(V),collapse='-')))
    ## Select number of patterns

    ## REB 2016-11-03
    ## If there is only one single station, avoid collapse of dimension
    if (is.null(dim(attr(X$pca,'pattern'))))
      dim(attr(X$pca,'pattern')) <- c(1,length(attr(X$pca,'pattern')))
    if (is.null(ip)) {
      U <- attr(X$pca,'pattern')
      W <- attr(X$pca,'eigenvalues')
    } else {
    ## If ip is specified, use a sub set of the PCA modes.
      U <- attr(X$pca,'pattern')[,ip]
      W <- attr(X$pca,'eigenvalues')[ip]
      V <- V[,ip,]
    }    
      
    ## Multi-station case (REB 2016-11-03)
    if (verbose) print('multiple stations')
    d <- dim(U)
    S <- apply(V, 3, function(x) U %*% diag(W) %*% t(x))
    dim(S) <- c(dim(U)[1], dim(V)[1], length(X)-2)
    for (i in seq(1:dim(S)[1])) {
      S[i,,] <- S[i,,] + c(attr(X$pca,'mean'))[i]
    }
    S <- lapply(split(S, arrayInd(seq_along(S),dim(S))[,1]),
                array,dim=dim(S)[-1])
    S <- lapply(S,function(x) zoo(x,order.by=index(X[[3]])))
    if (verbose) print('Set attributes')
    Y <- as.station(X$pca,verbose=verbose)
    locations <- gsub("[[:space:][:punct:]]","_",tolower(attr(Y,"location")))
    locations <- gsub("__","_",locations)
    ##locations <- paste(paste("i",attr(X$pca,"station_id"),sep=""),
    ##                   locations,sep=".")
    S <- setNames(S,locations)
    param <- attr(X$pca,"variable")[1]
    longname <- attr(X$pca,"longname")[1]
    unit <- attr(X$pca,"unit")[1]
    lons <- attr(X$pca,"longitude")
    lats <- attr(X$pca,"latitude")
    alts <- attr(X$pca,"altitude")
    stid <- attr(X$pca,"station_id")
    locs <- attr(X$pca,"location")
    gcms <- sub("^i[0-9]{1,3}_","",names(X)[3:length(X)])
    for (i in 1:length(S)) {
      yi <- Y[,i]
      class(yi) <- class(Y)
      yi <- attrcp(Y,yi)
      attr(yi,"longitude") <- lons[i]
      attr(yi,"latitude") <- lats[i]
      attr(yi,"altitude") <- alts[i]
      attr(yi,"station_id") <- stid[i]
      attr(yi,"location") <- locs[i]
      attr(yi,"unit") <- unit
      attr(yi,"variable") <- param
      attr(yi,"longname") <- longname
      attr(S[[i]],"station") <- yi
      attr(S[[i]],'aspect') <- 'original'
      attr(S[[i]],"longitude") <- lons[i]
      attr(S[[i]],"latitude") <- lats[i]
      attr(S[[i]],"altitude") <- alts[i]
      attr(S[[i]],"station_id") <- stid[i]
      attr(S[[i]],"location") <- locs[i]
      attr(S[[i]],'model_id') <- gcms
      class(S[[i]]) <- c('dsensemble','zoo')
    }
    if (!is.null(is)) S <- subset(S,is=is,verbose=verbose)
    if (length(S)>1) {
      class(S) <- c("dsensemble",class(X$pca)[2:3],"list") 
    } else {
      S <-  S[[1]]
    }
    #REB 2018-03-02: The line below causes big problems. Besides, I don't understand why it's there
    #if ( (is.list(S)) & (length(S)==length(locs)) ) names(S) <- locs
    attr(S,"unit") <- unit
    attr(S,"variable") <- param
    attr(S,"longname") <- longname
    attr(S,"aspect") <- "dsensemble.pca transformed to stations"
    attr(S,"history") <- history.stamp()
    invisible(S)
  }
}

as.station.dsensemble.station <- function(x,is=NULL,it=NULL,FUN='mean',verbose=FALSE,...) {

    if (verbose) print('as.station.dsensemble.station')
    ns <- length(x)
    ## Find the size of the PC matrices representing model projections
    d <- apply(sapply(x,dim),1,min); nt <- d[1]
    ## The PCs from the list are extracted into the matrix V 
    #V <- array(unlist(lapply( X[3:length(X)],
    #  function(x) coredata(x[1:d[1],1:d[2]]))),dim=c(d,length(X)-2))
    V <- unlist(lapply(x,FUN=
             function(x,FUN) apply(coredata(x),1,FUN=FUN),FUN))
    dim(V) <- c(nt,ns)
    if (verbose) str(V)
    loc <- unlist(lapply(x,esd::loc)); lon <- unlist(lapply(x,esd::lon)); lat <- unlist(lapply(x,esd::lat))
    alt <- unlist(lapply(x,esd::alt)); stid <- unlist(lapply(x,esd::stid));
    param <- attr(x,"variable")#unlist(lapply(x,varid))
    unit <- attr(x,"unit")#unlist(lapply(x,unit))
    longname <- attr(x,"longname")#unlist(lapply(x,function(x) attr(x,'longname')))
    y <- as.station(zoo(V,order.by=index(x[[1]])),loc=loc, param=param,unit=unit,
                    lon=lon,lat=lat,alt=alt,stid=stid,longname=longname)
    attr(y,"history") <- history.stamp()
    return(y)
}

as.pca <- function(x) UseMethod("as.pca")

as.pca.ds <- function(x) {
  stopifnot(inherits(x,'pca'))
  class(x) <- class(x)[-1]
  return(x)
}

as.pca.station <- function(x) {
  y <- PCA(x)
  return(y)
}


as.ds <- function(x) {
  y <- zoo(x,order.by=index)
  attr(y,'location') <- attr(x,'location')
  attr(y,'variable') <- attr(x,'variable')
  attr(y,'unit') <- attr(x,'unit')
  attr(y,'longitude') <- attr(x,'longitude')
  attr(y,'latitude') <- attr(x,'latitude')
  attr(y,'altitude') <- attr(x,'altitude')
  attr(y,'country') <- attr(x,'country')
  attr(y,'longname') <- attr(x,'longname')
  attr(y,'station_id') <- attr(x,'station_id')
  attr(y,'quality') <- attr(x,'quality') 
  attr(y,'calendar') <- attr(x,'calendar')
  attr(y,'source') <- attr(x,'source')
  attr(y,'URL') <- attr(x,'URL')
  #attr(y,'history') <- attr(x,'history')
  #attr(y,'date-stamp') <- attr(x,'date-stamp')
  attr(y,'type') <- attr(x,'type')
  attr(y,'aspect') <- attr(x,'aspect')
  attr(y,'reference') <- attr(x,'reference')
  attr(y,'info') <- attr(x,'info')
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  class(y) <- c("station","month","zoo")
  return(y)
}

as.comb <- function(x,...) UseMethod("as.comb")

as.comb.eof <- function(x,...) {
  stopifnot(inherits(x,'eof'))
  #print("as.field.eof")
  y <- eof2field(as.eof(x))
  #print("HERE")
  n <- attr(x,'n.apps')
  if (!is.null(n)) {
    for ( i in 1:n ) {
      z <- as.eof(x,i)
      #plot(z)
      Z <- eof2field(z)
      eval(parse(text=paste("Z -> attr(y,'appendix.",i,"')",sep="")))
    }
  }
  n -> attr(y,'n.apps')
  attr(y,'history') <- history.stamp(x)
  attr(y,'quality') <- c(attr(x,'quality'),'EOF-filtered')
  return(y)
}



as.field <- function(x,...) UseMethod("as.field")

as.field.zoo <- function(x,lon,lat,param,unit,
                         longname=NA,quality=NA,src=NA,url=NA,
                         reference=NA,info=NA,calendar='gregorian',
                         greenwich=TRUE, method= NA,type=NA,aspect=NA,
                         verbose=FALSE) {
  if(verbose) print("as.field.zoo")
  #print("lon"); print(lon); print("lat"); print(lat); print(param); print(unit)
  #print(c(length(lon),length(lat),length(index)))

  t <- index(x)
  #print(length(t))
  #dyr <- as.numeric(format(t[2],'%Y')) - as.numeric(format(t[1],'%Y')) 
  #dmo <- as.numeric(format(t[2],'%m')) - as.numeric(format(t[1],'%m')) 
  #dda <- as.numeric(format(t[2],'%d')) - as.numeric(format(t[1],'%d'))
  if (length(year(x))!=1) {
      dyr <- median(diff(year(x)))#diff(year(x))[1]
      dmo <- median(diff(month(x)))#diff(month(x))[1]
      dda <- median(diff(day(x)))#diff(day(x))[1]
      timescale <- "annual"
      if (dmo>0)  timescale <- "month"
      if (dmo==3)  timescale <- "season"
      if (dda>0)  timescale <- "day"
      if (dyr==0 & dmo==0 & dda==0) timescale <- "sub-daily"
      ##print(timescale)
  } else { 
    timescale <- "day"
  }
  # Add attributes to x
  attr(x,"variable") <- param
  attr(x,"longname") <- longname
  attr(x,"unit") <- unit
  attr(x,"source") <- src
  attr(x,"dimensions") <- c(length(lon),length(lat),length(t))
  attr(x,"longitude") <- lon
  attr(x,"latitude") <- lat
  attr(x,"greenwich") <- greenwich
  attr(x,"calendar") <- calendar 
  attr(x,'type') <- type
  attr(x,'aspect') <- aspect
  attr(x,'reference') <- reference
  attr(x,'info') <- attr(x,'info')
  attr(x,'history') <- history.stamp(x)
  class(x) <- c("field",timescale,"zoo")
  invisible(x)
}



as.field.default <- function(x,index,lon,lat,param,unit,
                         longname=NA,quality=NA,src=NA,url=NA,
                         reference=NA,info=NA,calendar='gregorian',
                         greenwich=TRUE, method= NA,type=NA,aspect=NA,
                         verbose=FALSE) {

if(verbose) print("as.field.default")
#create a zoo object z
  z <- zoo(x=x,order.by=index)
  x <- as.field.zoo(z,lon=lon,lat=lat,param=param,unit=unit,
                    longname=longname,quality=quality,src=src,url=url,
                    reference=reference,info=info,calendar=calendar,
                    greenwich=greenwich, method=method,type=type,
                    aspect=aspect, verbose=verbose)
  invisible(x)
}


as.field.field <- function(x,...) {
  if (inherits(x,'comb')) x <- as.field.comb(x,...)
  return(x)
}

as.field.comb <- function(x,iapp=NULL,verbose=FALSE,...) {
  if(verbose) print("as.field.comb")
  if (is.null(iapp)) {
    # Drop the appendend fields:
    n <- attr(x,'n.apps')
    for ( i in 1:n ) eval(parse(text=paste("attr(x,'appendix.",i,"') <- NULL",sep="")))
    attr(x,'n.apps') <- NULL
    class(x) <- class(x)[-1]
    return(x)
  }
  # Select one of the appended fields
  eval(parse(text=paste("y <- attr(x,'appendix.",iapp,"')",sep="")))
  attr(y,'longitude') <- attr(x,'longitude')
  attr(y,'latitude') <- attr(x,'latitude')
  attr(y,'history') <- history.stamp(x)
  return(y)  
}

as.field.eof <- function(x,iapp=NULL,anomaly=FALSE,verbose=FALSE,...) {
  if(verbose) print("as.field.eof")
  if (inherits(x,'dsensemble')) {
    y <- as.field.dsensemble.eof(x,verbose=verbose,...)
  } else if (!inherits(x,'comb')) {
    y <- eof2field(x,verbose=verbose,...)
  } else {
    y <- as.eof(x,iapp,verbose=verbose)
    y <- eof2field(y,verbose=verbose,anomaly=anomaly,...)
  }
  return(y)
}

as.field.ds <- function(x,iapp=NULL,verbose=FALSE,...) {
  if(verbose) print("as.field.ds")
  if (inherits(x,'eof')) {
    class(x) <- class(x)[-1]
    ## REB a few lines to catch cases where ds has not caught the comb-aspects.
    if (!is.null(iapp)) {
      if (!is.null(attr(x,'n.apps'))) {
        class(x)[length(class(x))+1]<-'comb'
        y <- as.field.eof(x,iapp,...)
        return(y)}
    }else{
    y <- as.field.eof(x,iapp,...)
    ## The residuals
    fit <- attr(x,'fitted_values')
    fit <- attrcp(attr(x,'eof'),fit)
    class(fit) <- class(attr(x,'eof'))
    attr(y,'fitted_values') <- fit
    attr(y,'original_data') <- attr(x,'original_data')
    attr(y,'calibration_data') <- attr(x,'calibration_data')
  }} else y <- NULL
  return(y)
}

as.field.station <- function(x,lon=NULL,lat=NULL,nx=30,ny=30,
                             verbose=FALSE,...) {
  if(verbose) print("as.field.station")
  if (is.null(lon)) lon <- seq(min(lon(x)),max(lon(x)),length=nx)
  if (is.null(lat)) lat <- seq(min(lat(x)),max(lat(x)),length=ny)
  y <- regrid(x,is=list(lon=lon,lat=lat))
  attr(y,'history') <- history.stamp(x)
  return(y)  
}

as.field.dsensemble.eof <- function(X,is=NULL,ip=NULL,im=NULL,anomaly=FALSE,verbose=FALSE,...) {
  if (verbose) print('as.field.dsensemble.eof')
  stopifnot(inherits(X,"dsensemble") & inherits(X,"eof"))
  if (inherits(X,"field")) {
      invisible(X)
  } else {
    #if (is.null(is)) is <- 1:length(loc(X$pca))
    if (verbose) print('Extract the results model-wise')
    ## KMP 2016-01-04: select ensemble members with im
    if(is.null(im)) {
      ix <- 3:length(X)
    } else {
      ix <- im[im>0 & im<(length(X)-2)] + 2
    }
    #ix <- 3:length(X)
    d <- apply(sapply(X[ix],dim),1,min)
    V <- array(unlist(lapply( X[ix],
      function(x) coredata(x[1:d[1],1:d[2]]))),dim=c(d,length(ix)))
    if (is.null(ip)) {
      U <- attr(X$eof,'pattern')
      W <- attr(X$eof,'eigenvalues')
    } else {
    ## If ip is specified, use a subset of the PCA modes.
      U <- attr(X$eof,'pattern')[,,ip]
      W <- attr(X$eof,'eigenvalues')[ip]
      V <- V[,ip,]
    }    
    d <- dim(U)
    dim(U) <- c(d[1]*d[2],d[3])
    S <- apply(V, 3, function(x) U %*% diag(W) %*% t(x))
    dim(S) <- c(dim(U)[1], dim(V)[1], dim(V)[3])
    if(!anomaly) {
      for (i in seq(1:dim(S)[1])) {
        S[i,,] <- S[i,,] + c(attr(X$eof,'mean'))[i]
      }
    }

    S <- aperm(S,c(3,2,1))
    S <- lapply(split(S, arrayInd(seq_along(S),dim(S))[,1]),
                array,dim=dim(S)[2:3])
    
    if (verbose) print('Set attributes')
    Y <- as.field(X$eof)
    gcms <- sub(".*_","",names(X)[ix])
    S <- setNames(S,gcms)
    for (i in seq_along(ix)) {
      S[[i]] <- as.field(S[[i]],index=index(X[[ix[i]]]),
                 lon=attr(Y,"longitude"),lat=attr(Y,"latitude"),
                 param=varid(Y),unit=unit(Y),
                 longname=paste('fitted',attr(Y,'longname')),
                 greenwich=attr(Y,'greenwich'),aspect='fitted')
      attr(S[[i]],'unitarea') <- attr(X,"unitarea")
      class(S[[i]]) <- class(Y)
    }
    if (!is.null(is)) S <- subset(S,is=is,verbose=verbose)
    class(S) <- c("dsensemble","field","list")
    S <- attrcp(X,S)
    attr(S,"unit") <- attr(Y,"unit")
    attr(S,"variable") <- attr(Y,"variable")
    attr(S,"longname") <- attr(Y,"longname")
    attr(S,"aspect") <- "dsensemble.eof transformed to field"
    if(anomaly) {
      attr(S,"aspect") <- c(attr(S,"aspect"), "anomaly")
      attr(S,"mean") <- attr(X$eof,"mean")
    }
    attr(S,"history") <- history.stamp()
    invisible(S)
  }
}


#as.annual <- function(x) {
#  year <- format(index), '%Y') 
#  return(year)
#  y <- annual(x)
#}
as.annual <- function(x, ...) UseMethod("as.annual")
#as.annual.default <- function(x, ...) as.annual(as.numeric(x))
as.annual.default <- function(x, ...) annual(x,...)
as.annual.numeric <- function(x, ...) annual(x)
as.annual.integer <- function(x, ...) structure(x, class = "annual")
as.annual.yearqtr <- function(x, frac = 0, ...) {
    if (frac == 0) annual(as.numeric(x)) else
    as.annual(as.Date(x, frac = frac), ...)
}
as.annual.station <- function(x, ...) annual.station(x,...)

as.annual.spell <- function(x, ...) annual.spell(x,...)

#as.annual <- function(x,FUN=mean,...) {
#  y <- annual(x,FUN=match.fun(FUN), ...)
#  return(y)
#}

as.monthly <- function(x,...) UseMethod("as.monthly")

yyyymm <- function(x) ym <- as.Date(paste(year(x),month(x),'01',sep='-'))

as.monthly.default <- function(x,...) {
  y <- aggregate(x,by=yyyymm,...)
  return(y)
}


as.monthly.field <- function(x,FUN='mean',...) {
if (inherits(x,'month')) return(x)
  y <- aggregate(as.zoo(x), yyyymm, #function(tt) as.Date(as.yearmon(tt)),
                 FUN=FUN,...)
  y <- attrcp(x,y)
  attr(y,"dimensions") <- c(attr(x,"dimensions")[1:2],length(index(y)))
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  class(y)[2] <- "month" 
  return(y)
}

## This is a dublicate of that in as.R
as.monthly.station <- function (x, FUN = "mean", ...) 
{
    y <- aggregate(zoo(x), yyyymm, #function(tt) as.Date(as.yearmon(tt)), 
                   FUN = FUN, ...)
    y <- attrcp(x, y)
    attr(y, "history") <- history.stamp(x)
    class(y) <- class(x)
    class(y)[2] <- "month"
    return(y)
}


# Not to confuse with season
# This function extracts a given seasonal interval and aggregates a given statistic
as.seasons <- function(x,start='01-01',end='12-31',FUN='mean', ...) {
  IV <- function(x) sum(is.finite(x))
  yrs <- year(x); d <- dim(x)
  # ns = number of stations
  if (is.null(d)) ns <- 1 else ns <- d[2]
  years <- as.numeric(rownames(table(yrs))); n <- length(years)
  y <- matrix(rep(NA,n*ns),n,ns); k <- y
  start.1 <- as.numeric(as.Date(paste(years[1],start,sep='-')))
  end.1 <- as.numeric(as.Date(paste(years[1],end,sep='-')))
  if (start.1 > end.1) twoyears <- 1 else twoyears <- 0
  
  for (i in 1:n) {
    z <- coredata(window(x,start=as.Date(paste(years[i],start,sep='-')),
                           end=as.Date(paste(years[i]+twoyears,end,sep='-'))))
    k[i,] <- apply(matrix(z,length(z),ns),2,IV)
    y[i,] <- apply(matrix(z,length(z),ns),2,FUN, ...)
  }
  y <- zoo(y,order.by=as.Date(paste(years,start,sep='-')))
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  if (twoyears==0) attr(y,'season.interval') <- paste(start,'to',end) else
                   attr(y,'season.interval') <- paste(start,'to',end,'the following year')
  attr(y,'n.valid') <- k
  class(y) <- class(x)
  class(y)[2] <- "annual"
  return(y)
}


## Original script : as.R
## Author          : Rasmus
## Major updates   : 2013-12-05 Abdelkader 

## Comments        : ! fourseasons() is no longer used 

as.4seasons <- function(x,...) UseMethod("as.4seasons")

# Aggregate for the seasons DJF (winter), MAM (spring), JJA (summer), and
# SON (autumn). Need to combine december with the following year's Jan-Feb.
# Best solution: fourseasons is a function applied to index(x)?



as.4seasons.default <- function(x,FUN='mean',slow=FALSE,verbose=FALSE,nmin=NULL,...) {
  if(verbose) print('as.4seasons.default')
  if (inherits(x,'season')) return(x)
  attr(x,'names') <- NULL
  d <- dim(coredata(x))
  #print(d)
  if (is.null(d)) d <- c(length(x),1)
  if (!slow) {
    if (inherits(x,"month")) {
      if ( (is.null(d)) | (d[2]==1) ) {
        X <- c(NA,coredata(x)[1:length(x)-1]) # shift the coredata by 1 to start on December. This works only for monthly data !!!  
      } else {
        X <- rbind(rep(NA,d[2],1),coredata(x)[1:d[1]-1,])
      }
      #print(dim(X))
      X <- zoo(X,order.by=index(x))
      ##yrseas <- fourseasons(ix)
      ##print(yrseas)
      #print('aggregate')
      #print(names(list(...)))
      yq <- function(t) as.yearqtr(year(t) + 0.25*floor((month(t)-1)/3))
      y <- aggregate(x=as.zoo(X),by=yq,#as.yearqtr,
                     FUN=match.fun(FUN),...)
      # convert yearqtr to yearmon
      y <- zoo(x=y,order.by=as.Date(as.yearmon(index(y))))
    } else y <- as.4seasons.day(x,FUN=FUN,nmin=nmin,...)
    #y <- as.4seasons.day(x,FUN=match.fun(FUN),...)
    
    #print(dim(y))
    ok <- length(index(y))
    #print(summary(c(coredata(y))))
  } else {
    yr <- sort(rep(as.integer(rownames(table(year(x)))),4))
    n <- length(yr)

    q <- rbind(c(12,1,2),3:5,6:8,9:11)
    X <- matrix(rep(NA,n*d[2]),n,d[2]) 
    t <-rep(NA,n)

    #print("start loop")
    for (i in 1:n) {
      iq <- (i-1) %% 4 + 1
      if (iq == 1) 
        ii <- (is.element(year(x),yr[i]) &
               is.element(month(x),c(1,2))) |
              (is.element(year(x),yr[i]-1) &
               is.element(month(x),12))  else
        ii <-  is.element(year(x),yr[i]) &
               is.element(month(x),q[iq,])
     if (d[2]==1) cline <- paste(FUN,"(coredata(x[ii,]),...)",sep="") else
                  cline <- paste("apply(coredata(x[ii,]),2,",FUN,", ...)",sep="")
      #print(cline)
      if ( (inherits(x,'day')) & (sum(ii)>=85) |
           (inherits(x,'month')) & (sum(ii)>=3) ) {
        X[i,] <- eval(parse(text=cline))
        t[i] <- yr[i] + (iq-1)/4
        print(c(i,yr[i],round(mean(X[i,]),2),iq,sum(ii),round(X[i],2),t[i]))
      }
    } 
    #print("end loop")
    ok <- is.finite(rowMeans(X,na.rm=TRUE)) & is.finite(t)
    #print(table(as.yearqtr(t[ok])))
    #print(summary(c(X)))
    y <- zoo(X[ok,],order.by=as.Date(as.yearqtr(t[ok])))
    #names(y) <- FUN
  }
  if (!is.null(dim(y))) {
    ok <- is.finite(rowMeans(y,na.rm=TRUE))
    y <- y[ok,]
  } 
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  if (inherits(x,'field'))
    attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],sum(ok))
  class(y) <- class(x)
  class(y)[2] <- "season"
  #plot(y)
  return(y) 
}

as.4seasons.day <- function(x,FUN='mean',na.rm=TRUE,dateindex=TRUE,nmin=85,...) {
    ##browser()
  IV <- function(x) sum(is.finite(x))
    if (inherits(x,'month')) nmin <- 3 # AM 06-07-2015
  #print('as.4seasons.day')
  attr(x,'names') <- NULL  
  t <- index(x)
  year <- year(t) #as.numeric(format(t,'%Y'))
  month <- month(t) #as.numeric(format(t,'%m'))
  day <- day(t) # as.numeric(format(t,'%d'))
  #shift the time stamps by one month, sneaking December into the subsequent year
  month <- month + 1
  dec <- is.element(month,13)
  year[dec] <- year[dec] + 1
  month[dec] <- 1
  # Change the day to avoid warning that the calendar is wrong (e.g. due to
  # too many days in February). Since the data is aggregated, the exact day
  # in the month doesn't matter here.
  hour <- 12*(day - 2*trunc(day/2))
  day <- trunc(day/2) + 1
  #print(table(year)); print(table(month));
  #print(table(day)); print(table(hour));   
  #tshifted <- as.Date(paste(year,month,day,sep="-"))
  tshifted <-  ISOdate(year=year,month=month,day=day,hour=hour)
  #print(summary(tshifted))
  X <- zoo(coredata(x),order.by=tshifted)
 
  # Test for the presens of 'na.rm' in argument list - this is a crude fix and not a
  # very satisfactory one. Fails for FUN==primitive function.
  if (is.function(FUN)) test.na.rm <- FALSE else
      test.na.rm <- (sum(is.element(FUN,c('mean','min','max','sum','quantile')))>0)
  if ( (sum(is.element(names(formals(FUN)),'na.rm')==1)) | (test.na.rm) )
     y <- aggregate(X,as.yearqtr,FUN=match.fun(FUN),...,na.rm=na.rm) else
     y <- aggregate(X,as.yearqtr,FUN=match.fun(FUN),...)

  # Set to missing for seasons with small data samples:
  nd <- aggregate(X,as.yearqtr,FUN=IV)
  ok <- nd >= nmin  
  coredata(y)[!ok] <- NA
  # dateindex: convert "1775 Q1" to "1775-01-01"
  if (dateindex)
    y <- zoo(coredata(y),order.by=as.Date(index(y)))
  unit <- attr(y,'unit')
  y <- attrcp(x,y,ignore=c("unit","names"))
  unit -> attr(y,'unit')
  #str(y); print(unit)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  class(y)[2] <- "season"
  invisible(y)
}

as.4seasons.station <- function(x,FUN='mean',...) {
  #print('as.4seasons.station')
  y <- as.4seasons.default(x,FUN,...)
#  y <- attrcp(x,y)
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- class(x)
#  class(y) <- gsub("month","season",class(x))
  return(y) 
}

as.4seasons.spell <- function(x,FUN='mean',...) {
  y <- as.4seasons.default(as.station(x),FUN,...)
#  y <- attrcp(x,y)
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- class(x)
#  class(y) <- gsub("month","season",class(x))
  return(y) 
}


as.4seasons.field <- function(x,FUN='mean',verbose=FALSE,...) {
  if(verbose) print("as.4seasons.field")
  d <- attr(x,"dimensions")
  y <- as.4seasons.default(x,FUN,verbose=verbose,...)
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  ## Update dimensions
  attr(y, "dimensions") <- c(d[1],d[2],dim(y)[1])
  class(y) <- class(x)
  class(y) <- gsub("month","season",class(x))
  return(y)
}

as.4seasons.dsensemble <- function(x,FUN='mean',...) {
    ##browser()
    cls <- class(x)
    class(x) <- c("station",cls[2],"zoo") ## AM 06-07-2015 Quick fix here, time step added into the class of x
    attrx <- attributes(x)
    y <- as.4seasons.station(x,FUN,...)
    ##attributes(y) <- attrx
    
    y <- attrcp(x,y)
    attr(y,"station") <- as.4seasons.station(attr(x,"station"))
   
    attr(y,'history') <- history.stamp(x)
    
    class(y) <- c("dsensemble","season","zoo")
    return(y)
}

as.anomaly <- function(x,...) UseMethod("as.anomaly")

# REB 2015-03-23 Tidy up - use anomaly

as.anomaly.default <- function(x,ref=NULL,na.rm=TRUE) anomaly.default(x)

#
#as.anomaly.default <- function(x,ref=NULL,na.rm=TRUE) {
## The argument monthly can be used to force the method to be
## julian-day regression-based or based on monthly mean
 
##  print('as.anomaly.default')
##  yr <- as.integer(format(index(x),'%Y'))
##  mon <- as.integer(format(index(x),'%m'))
##  dy <- as.integer(format(index(x),'%d'))
##  str(x)
##  browser()
#  if ((is.numeric(x)) & (!is.null(attr(x, "names"))))
#  if ((!is.zoo(x)) & (!is.null(attr(x, "names"))))
#    x <- zoo(x,order.by=as.Date(attr(x, "names")))
#  yr <- year(x);  mon <- month(x);  dy <- day(x); seas <- season(x,format='numeric')
#  if (is.null(ref))
#    ref <- seq(min(yr,na.rm=TRUE),max(yr,na.rm=TRUE),by=1)
#  # Check whether the object is a time series of monthly data.
#  ndd <- length(table(dy)); nmm <- max(diff(mon))
#  #print(c(ndd,nmm))
#  if ((ndd==0) & (nmm==0)) {
#    # Only one month/season
#    return(x - mean(x,na.rm=TRUE))
#  }
#  
#  #print(ndd); print(class(x))
#  #if ( (ndd==1) & is.null(monthly) ) monthly <- TRUE else
#  #                                   monthly <- FALSE
#
#  if ( (inherits(x,'month')) | ((ndd==1) & (nmm==1)) ) {
#    # Monthly data
#    print('monthly')
#    monthly <- TRUE
#    nm <- 12
#  } else monthly <- FALSE
#  if ( (inherits(x,'seasonal')) | ((ndd==1) & (nmm==3)) ) {
#    # seasonal data
#    #print('seasonal')
#    seasonal <- TRUE
#    mon <- seas
#    nm <- 4
#  } else seasonal <- FALSE
#  
#  # Check if x contains one or more stations
#  d <- dim(x)
#
#  if (is.null(d)) {
#    # single records
#    #cat('.')
#    y <- coredata(x)
#    
#    if ( (monthly | seasonal) ) {
#    # If monthly, subtract the monthly mean
##    print(mon)
##    if (is.null(dim(x))) {
#      print("1D")
#      clim <- rep(0,nm)
#      for (i in 1:nm) {
#        im <- is.element(mon,i)
#        z <- mean(y[im & is.element(yr,ref)],na.rm=na.rm)
#        clim[i] <- z
#        y[im] <- y[im] - clim[i]
#      }
#      
#      y <- zoo(y,order.by=index(x))
##    } else {
##      #print("2D")
##      y <- t(y)
##      clim <- rep(0,length(y[,1])*12); dim(clim) <- c(length(y[,1]),12)
##      for (i in 1:12) {
##        im <- is.element(mon,i)
##        z <- rowMeans(y[,im & is.element(yr,ref)],na.rm=na.rm) 
##        #print(c(length(z),length(clim[,1]))); plot(z)
##        clim[,i] <- z
##        y[,im] <- y[,im] - clim[,i]
##        #print(c(i,sum(im),mean(clim[,i])))
##      }
##      y <- zoo(t(y),order.by=index(x))
#
#    } else {
##      print("daily"); print(class(x))
#      t0 <- julian(index(x)[is.element(yr,ref)]) -
#            julian(as.Date(paste(yr[is.element(yr,ref)],"-01-01",sep="")))
#      t <- julian(index(x)) -
#         julian(as.Date(paste(yr,"-01-01",sep="")))
#      c1 <- cos(pi*t0/365.25); s1 <- sin(pi*t0/365.25)
#      c2 <- cos(2*pi*t0/365.25); s2 <- sin(2*pi*t0/365.25)
#      c3 <- cos(3*pi*t0/365.25); s3 <- sin(3*pi*t0/365.25)
#      c4 <- cos(4*pi*t0/365.25); s4 <- sin(4*pi*t0/365.25)
#      C1 <- cos(pi*t/365.25);   S1 <- sin(pi*t/365.25)
#      C2 <- cos(2*pi*t/365.25); S2 <- sin(2*pi*t/365.25)
#      C3 <- cos(3*pi*t/365.25); S3 <- sin(3*pi*t/365.25)
#      C4 <- cos(4*pi*t/365.25); S4 <- sin(4*pi*t/365.25)
#      cal <- data.frame(y=coredata(x),c1=c1,c2=c2,c3=c3,c4=c4,
#                        s1=s1,s2=s2,s3=s3,s4=s4)
#      pre <- data.frame(c1=C1,c2=C2,c3=C3,c4=C4,
#                        s1=S1,s2=S2,s3=S3,s4=S4)
#      i1 <- is.element(year(x),year(x)[1])
#      pre1 <- data.frame(c1=C1[i1],c2=C2[i1],c3=C3[i1],c4=C4[i1],
#                         s1=S1[i1],s2=S2[i1],s3=S3[i1],s4=S4[i1])
#      acfit <- lm(y ~ c1 + s1 + c2 + s2 + c3 + s3 + c4 + s4,data=cal)
#      clim <- predict(acfit,newdata=pre)
#      y <- zoo(coredata(x) - clim,order.by=index(x))
#      clim <-  predict(acfit,newdata=pre1)
#    }
#  } else {
#    #print("many stations")
#    rownames(x) <- as.character(index(x))
#    y <- apply(x,2,FUN='as.anomaly.default',ref=ref)
#    y <- zoo(y,order.by=index(x))
#    clim <- x - y
#  }
#  #print("attributes")
#  y <- attrcp(x,y)
#  attr(y,"aspect") <- 'anomaly'
#  attr(y,"anomaly_method") <- monthly
#  attr(y,"climatology") <- clim
#  #attr(y,"date") <- date()
#  #attr(y,"call") <- match.call()
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- class(x)
#  invisible(y)
#}

as.anomaly.zoo <- function(x,ref=NULL,na.rm=TRUE,...) {
  y <- as.anomaly.station(x,ref=ref,na.rm=na.rm,...)
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.anomaly.list <- function(x,ref=NULL,na.rm=TRUE,...) {
  y <- lapply(x,anomaly(x))
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.anomaly.station <- function(x,ref=NULL,na.rm=TRUE,...) {
  y <- as.anomaly.default(x,ref=ref,na.rm=na.rm,...)
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.anomaly.field<- function(x,ref=NULL,na.rm=TRUE,...) {
   y <- anomaly.default(x,ref=ref,na.rm=na.rm,...)
   attr(y,'history') <- history.stamp(x)
   attr(y,'dimensions') <- attr(x,'dimensions')
   invisible(y)
}

# Handy conversion algorithms:
as.climatology <- function(x,...) {
    ya <- as.anomaly(x,...)
    clim <- coredata(attr(ya,'climatology'))
    if (!is.null(dim(clim)))
        len.clim <- dim(clim)[1]
    else
        len.clim <- length(clim)
    y <- zoo(clim,order.by=1:len.clim)      
    y <- attrcp(x,y)
    attr(y,'aspect') <- 'climatology'
    attr(y,'history') <- history.stamp(x)
    class(y) <- class(x)
    invisible(y)
}

as.residual <- function(x,...) UseMethod("as.residual")

as.residual.ds <- function(x,verbose=FALSE){
  if (verbose) print('as.residual.ds')
  if (is.ds(x)) {
    ## If the predictand was originally an EOF or PCA product, then
    ## the residual needs to inherits their attributes
      if (verbose) print('Re-construct and re-compute')
      ## Need to reconstruct the data matrix and re-calculate the EOFs/PCAs 
      if (is.eof(x)) {
        if (verbose) print('eof/field')
        z0 <- as.field(attr(x,'original_data'))
        z1 <- as.field(x)
        y <- z1 - z0
        y <- attrcp(z0,y); class(y) <- class(z0)
      } else
      if (is.pca(x)) {
        if (verbose) print('pca/station')
        z0 <- as.station(attr(x,'original_data'))
        z1 <- as.station(x)
        y <- z1 - z0
        y <- attrcp(z0,y); class(y) <- class(z0)
      } else
      if (is.station(x)) {
        if (verbose) print('station')
        z0 <- attr(x,'original_data')
        z1 <- x
        y <- z1 - z0
        y <- attrcp(z0,y); class(y) <- class(z0)
      }      
   } else
  ## If the results are a field object, then the residuals are stored as EOFs.
  if (is.field(x)) {
    if (verbose) {print('x is a field object'); print(class(x))}
    y <- as.field(attr(x,'original_data')) -
         as.field(attr(x,'fitted_values'))
    y <- attrcp(attr(x,'original_data'),y)
    class(y) <- class(attr(x,'original_data'))
    y <- as.field(y)
    attr(y,'aspect') <- 'residual'
  }
  attr(y,'aspect') <- 'residual'
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.residual.station <- function(x){
  if (!is.null(attr(x,'calibration_data')))
    y <- as.residual.ds(x) else y <- NULL
  invisible(y)
}

as.residual.field <- function(x){
  if (!is.null(attr(x,'calibration_data')))
    y <- as.residual.ds(x) else y <- NULL
  invisible(y)
}


as.calibrationdata <- function(x) UseMethod("as.calibrationdata")

as.calibrationdata.ds <- function(x) {
  y <- attr(x,'calibration_data')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'calibration_data'))
  invisible(y)
}

as.calibrationdata.station <- function(x) {
   if (!is.null(attr(x,'calibration_data')))
    y <- as.calibrationdata.ds(x) else y <- NULL
  invisible(y)
}

as.fitted.values <- function(x) UseMethod("as.fitted.values")

as.fitted.values.ds <- function(x) {
  y <- attr(x,'fitted_values')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'fitted_values'))
  invisible(y)
}

as.fitted.values.station <- function(x) {
   if (!is.null(attr(x,'fitted.values')))
    y <- as.fitted.values.ds(x) else y <- NULL
  invisible(y)
}



as.original.data <- function(x) UseMethod("as.original.data")

as.original.data.ds <- function(x) {
  y <- attr(x,'original_data')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'original_data'))
  invisible(y)
}

as.original.data.station <- function(x) {
  y <- as.original.data.ds(x)
  invisible(y)
}


as.pattern <- function(x,...) UseMethod("as.pattern")

as.pattern.ds <- function(x) {
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)
}

as.pattern.eof <- function(x) {
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)  
}

as.pattern.mvr <- function(x) {
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)  
}

as.pattern.cca <- function(x) {
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)  
}

as.pattern.trend <- function(x) {
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)  
}

as.pattern.matrix <- function(x) x 

as.pattern.array <- function(x) x 

as.pattern.field <- function(x,FUN=NULL,...) {
  if (!is.null(FUN)) {
    y <- apply(x,2,FUN,...)
    dim(y) <- attr(x,'dimension')[1:2]
  } else {
    y <- t(coredata(x))
    dim(y) <- attr(x,'dimension')
  }
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- varid(x)
  attr(y,'unit') <- unit(x)
  class(y) <- c('pattern',class(y))
  attr(y,'time') <- index(x) 
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.pattern.corfield <- function(x) {
  y <- coredata(x)
  dim(y) <- attr(x,'dimension')[1:2]
  y <- attrcp(x,y)
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- varid(x)
  attr(y,'unit') <- unit(x)
  class(y) <- c('pattern',class(y))
  attr(y,'time') <- index(x)
  invisible(y)
}

as.eof <- function(x,...) UseMethod("as.eof")

as.eof.zoo <- function(x,...) {
  class(x) <- c('eof','zoo')
  return(x)
}

as.eof.ds <- function(x,iapp=NULL) {
  y <- as.eof(attr(x,'eof'),iapp) 
  return(y)
}

as.eof.eof <-function(x,iapp=NULL,...) {
  if (inherits(x,'comb')) x <- as.eof.comb(x,iapp) else
  return(x)
}
  
as.eof.comb <- function(x,iapp=NULL) {
  #print("as.eof.comb")
  stopifnot(inherits(x,'comb'))

  # if x is a 'field'
  if (!inherits(x,'eof')) x <- EOF(x)

  # assume x from now on is an 'eof'
  if (!is.null(iapp)) {
    y <- as.eof.appendix(x,iapp)
    return(y)
  }
  class(x) <- class(x)[-grep('comb',class(x))]
  napps <- attr(x,'n.apps')
  for (i in seq(napps)) {
    eval(parse(text=paste("attr(x,'appendix.",i,"') <- NULL",sep="")))
  }
  attr(x,'n.apps') <- NULL
  attr(x,'history') <- history.stamp(x)
  return(x)
}

as.eof.field <- function(x,iapp=NULL,...) {
  y <- EOF(x,...)
  if (!is.null(iapp)) y <- as.eof.appendix(y,iapp)
  return(y)
}

as.eof.appendix <- function(x,iapp=1,verbose=FALSE) {
  if (verbose) print("as.eof.appendix")
  clim <- eval(parse(text=paste("attr(attr(x,'appendix.",iapp,"'),'climatology')",sep="")))
  aveg <- eval(parse(text=paste("attr(attr(x,'appendix.",iapp,"'),'mean')",sep="")))
  stopifnot(inherits(x,'comb'))
  y <- eval(parse(text=paste("attr(x,'appendix.",iapp,"')",sep="")))
  x <- as.eof.comb(x)
  y <- attrcp(x,y)
  if (!is.null(clim)) attr(y,'climatology') <- clim 
  if (!is.null(aveg)) attr(y,'mean') <- aveg
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  return(y)
}

as.eof.list <- function(x,verbose=FALSE) {
  stopifnot(inherits(x,'list'),inherits(x[[1]],'eof'))
  if (verbose) print('as.eof.list')
  
  wPC <- function(z,iapp=NULL) {
    eigv <- attr(z,'eigenvalues')
    w <- eigv/sum(eigv)
    if (is.null(iapp)) Z <- z %*% diag(w) else
                       Z <- attr(z,paste('appendix.',iapp,sep='')) %*% diag(w)
    Z <- zoo(Z,order.by=index(z))
    return(Z)
  }

  if (verbose) try(print(summary(x)))
  if (inherits(x[[1]],'character')) x[[1]] <- NULL
  if (inherits(x[[1]],'eof')) {eof <- x[[1]]; x[[1]] <- NULL}
  X.list <- lapply(x,wPC)
  X <- do.call("merge", X.list)
  if (verbose) print(summary(X))
  t <- index(X)
  udv <- svd(coredata(X))
  eof <- zoo(udv$u[,1:20],order.by=t)
  attr(eof,'eigenvalues') <- udv$d
  pattern <- rep(1,dim(udv$v)[1])
  names(pattern) <- names(X)
  attr(eof,'pattern') <- pattern
  if (inherits(x[[1]],'comb')) {
    if (verbose) print('Combined field: appendix.1')
    for (i in 1:attr( attr(x[[1]],'n.apps'))) {
      z.list <- lapply(x,wPC,iapp=i)
      udv1 <- svd(coredata(do.call("merge", z.list)))
      attr(eof,paste('appendix.',i,sep='')) <- zoo(udv1$u[,1:20],
               order.by=index(attr(x,paste('appendix.',i,sep=''))))
      names(attr(eof,paste('appendix.',i,sep=''))) <- paste("X.",1:20,sep="")
    }
  }
  attr(eof,'original.list.of.eofs') <- x
  attr(eof,'udv') <- udv
  id <- c()
  for (i in 1:length(x)) id <- c(id,rep(i,length(attr(x[[i]],'eigenvalues'))))
  attr(eof,'id') <- id
  names(eof) <- paste("X.",1:20,sep="")
  class(eof) <- class(x[[1]])
  return(eof)
}

as.eof.dsensemble <- function(x,FUN='mean',verbose=FALSE) {
  ## R.E. Benestad, 2017-05-19
  ## Convert the dsensemble object to an EOF of the multi-model mean
  stopifnot(inherits(x,'dsensemble'),inherits(x[[2]],'eof')|inherits(x[[2]],'pca'))
  if (verbose) print('as.eof.dsensemble')
  eof0 <- x[[2]]; x[[2]] <- NULL
  x[[1]] -> info; x[[1]] <- NULL
  d <- c(dim(x[[1]]),length(x))
  y <- unlist(x)
  dim(y) <- c(d[1]*d[2],d[3])
  Y <- apply(y,1,FUN)
  dim(Y) <- c(d[1],d[2])
  eof <- zoo(Y,order.by=index(x[[1]]))
  eof <- attrcp(eof0,eof)
  class(eof) <- class(eof0)
  attr(eof,'info') <- info
  attr(eof,'history') <- history.stamp()
  return(eof)
}



as.appended <- function(x,...) UseMethod("as.appended")

as.appended.ds.comb <- function(x,iapp=1,verbose=FALSE) {
  if(verbose) print("as.appended.ds.comb")
  eval(parse(text=paste("X <- attr(x,'appendix.",iapp,"')",sep="")))
  X <- attrcp(x,X,ignore='appendix')
  attr(X,'history') <- history.stamp(x)
  invisible(X)
}

as.appended.eof.comb <- function(x,iapp=1) {
  X <- as.appended.ds.comb(x,iapp=iapp)
  invisible(X)
}

as.appended.field.comb <- function(x,iapp=1) {
  X <- as.appended.ds.comb(x,iapp=iapp)
  invisible(X)
}

as.stand <- function(x,...) UseMethod("as.stand")

as.stand.station <- function(x,verbose=FALSE,na.rm=TRUE) {
  if(verbose) print("as.stand.station")
  if (is.precip(x)) {
    mu <- apply(x,2,mean,na.rm=na.rm)
    X <- 100*x/mu
    attr(X,'clim') <- mu
    attr(X,'aspect') <- 'proportional'
    attr(X,'unit') <- '%'
    attr(X,'oldunit') <- attr(x,'unit')
  } else if (is.T(x)) {
    mu <- apply(x,2,mean,na.rm=na.rm)
    sigma <- apply(x,2,sd,na.rm=na.rm)
    X <- (x - mu)/sigma
    attr(X,'mean') <- mu
    attr(X,'sigma') <- sigma
    attr(X,'aspect') <- 'standardised'    
  }
  attr(X,'history') <- history.stamp(x)
  return(X)
}



as.original <- function(x) UseMethod("as.original")

as.original.station <- function(x) {
  if (attr(x,'aspect')=='proportional') {
    X <- attr(x,'clim')*x/100
    attr(X,'clim') <- NULL
    attr(X,'unit') <- attr(x,'oldunit')
    attr(X,'oldunit') <- NULL
    attr(X,'aspect') <- 'original'
  } else if (attr(x,'aspect')=='standardised') {
    X <- x * attr(x,'sigma') + attr(x,'mean')
    attr(X,'mean') <- NULL
    attr(X,'sigma') <- NULL
    attr(X,'aspect') <- 'original'     
  } else X <- x
  attr(X,'history') <- history.stamp(x)
  return(X)
}

as.events <- function(x,...) UseMethod("as.events")

as.events.default <- function(x,label=NULL,dx=NULL,dy=NULL,
                      units=NULL,longname=NULL,variable=NULL,calendar=NULL,
                      qflabel=NULL,method=NULL,src=NULL,reference=NULL,
                      file=NULL,version=NULL,url=NULL,verbose=FALSE) {
  if (verbose) print("as.events")
  X <- data.frame(x)
  n <- names(X)
  if (!all(c("date","time","lon","lat") %in% names(X))) {
    print(paste("Missing input:",
     names(X)[!c("date","time","lon","lat")%in%names(X)]))
  }
  attr(X,"label") <- label
  attr(X,"dx") <- dx
  attr(X,"dy") <- dy
  attr(X,"longname") <- longname
  attr(X,"calendar") <- calendar
  attr(X,"variable") <- variable
  attr(X,"quality") <- qflabel
  attr(X,"units") <- units
  attr(X,"source") <- src
  attr(X,"file") <- file
  attr(X,"version") <- version
  attr(X,"method") <- method
  attr(X,"URL") <- url
  attr(X,"reference") <- reference
  class(X) <- c("events",class(X))
  attr(X,"history") <- history.stamp(X)
  invisible(X)
}

as.events.trajectory <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.events.trajectory")
  stopifnot(inherits(x,"trajectory"))
  invisible(x)
}

as.field.events <- function(x,...) {
  y <- events2field(x,...)
  return(y)
}

as.field.trajectory <- function(x,...) {
  y <- trajectory2field(x,...)
  return(y)
}

as.station.trajectory <- function(x,...) {
  y <- trajectory2station(x,...)
  return(y)
}
