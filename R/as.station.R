#' Coerce input to a \code{station} object
#' 
#' Transform an input object into the esd class \code{station}. 
#' \code{as.station} is an S3 method and will redirect to a fitting function depending on the type of input data.
#'
#' \code{as.station.zoo} and \code{as.station.data.frame} adds attributes and changes the class to transform the input to a 'station' object.
#'
#' \code{as.station.field} returns an object where every grid box is represented as one station.
#'
#' \code{as.station.pca} transforms a 'pca' object (see \code{\link{PCA}}) to a 'station' object using the method \code{\link{pca2station}}.
#'
#' \code{as.station.eof} represents the principle components of an 'eof' object as different stations.
#'
#' \code{as.station.dsensemble} transform a \code{dsensemble} object to a \code{station} object in one of two ways:
#' i) If the input is of class \code{dsensemble} \code{pca}, you are redirected to \code{as.station.dsensemble.pca}
#' which calculates the downscaled ensemble for each station based on the downscaled ensemble of principle components,
#' returning a \code{dsensemble} \code{station} or \code{dsensemble} \code{list} object.
#' ii) If the input is a \code{dsensemble} \code{station} or \code{dsensemble} \code{list} object, you are redirected to \code{as.station.dsensemble.station}
#' which returns a \code{station} object holding the ensemble mean (or another statistical characteristic of the ensemble, see input argument \code{FUN})
#' of the downscaled results for each station.
#'
#' \code{as.station.events} and \code{as.station.trajectory} aggregate an 'event' or 'trajectory' object to a 'station' object by
#' aggregating some aspect of the cyclones/anti-cyclones (or other type of event). By default, the total number of events/trajectories per month is calculated 
#' but the method can also estimate some other characteristic, e.g., the monthly mean sea level pressure at the center of the cyclones ().
#' 
#' @aliases as.station as.station.zoo as.station.data.frame as.station.ds as.station.pca as.station.list
#' as.station.spell as.station.events as.station.dsensemble.station as.station.eof as.station.field
#' @seealso as.station.dsensemble as.station.dsensemble.pca as.station.dsensemble.station 
#' 
#' @param x the input object
#' @param loc location(s), e.g, "Manchester" or c("Oslo","Bergen")
#' @param param short name of variable
#' @param unit unit of variable, e.g., 't2m'
#' @param lon longitude(s), a numerical or numerical vector
#' @param lat latitudes(s), a numerical or numerical vector
#' @param alt altitude(s), a numerical or numerical vector
#' @param cntr country, a character string or vector of character strings
#' @param longname long name of variable, e.g, 'temperature at 2m'
#' @param calendar calendar type
#' @param stid station id, a numerical or numerical vector
#' @param quality quality flag
#' @param src source of data
#' @param url url to website where data can be downloaded
#' @param reference reference describing data set
#' @param info additional information
#' @param method method applied to data
#' @param type type of data
#' @param aspect aspect describing data, e.g., 'original', 'anomaly', 'climatology'
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{station} object
#' 
#' @examples
#' # How to generate a new 'station' object
#' data <- round(matrix(rnorm(20*12),20,12),2)
#' colnames(data) <- month.abb
#' x <- data.frame(year=1981:2000,data)
#' X <- as.station(x,loc="",param="noise",unit="none")
#'
#' # Transform a field object into a station object
#' slp.field <- slp.DNMI(lon=c(-20,20), lat=c(50,70)) # get example SLP data
#' slp.station <- as.station(slp.field) # coerce SLP field to a station object
#' cb <- list(pal="burd", breaks=seq(1000,1020,2)) # specify color bar for maps
#' map(slp.field, FUN="mean", colbar=cb) # show map of SLP field
#' map(slp.station, FUN="mean", colbar=cb) # show map of SLP as station data
#'
#' @export as.station
as.station <- function(x,...) UseMethod("as.station")

#' @export as.station.zoo
as.station.zoo <- function(x,...,loc=NA,param=NA,unit=NA,lon=NA,lat=NA,alt=NA,
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

#' @export as.station.data.frame
as.station.data.frame <-  function (x,...,loc=NA,param=NA,unit=NA,lon=NA,lat=NA,alt=NA,
                                    cntr=NA,longname=NA,stid=NA,quality=NA,src=NA,url=NA,
                                    reference=NA,info=NA,method=NA,type=NA,aspect=NA,verbose=FALSE) {
  if(verbose) print("as.station.data.frame")
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

#' @export as.station.ds
as.station.ds <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.station.ds")
  if (inherits(x,'pca')) {
    class(x) <- class(x)[-1]
    y <- as.station.pca(x,...,verbose=verbose)
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

#' @export as.station.pca
as.station.pca <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.station.pca")
  if (verbose) print(names(attributes(attr(x,'original_data'))))
  if (inherits(x,"dsensemble")) {
    if (verbose) print('data type: dsensemble')
    y <- as.station.dsensemble.pca(x,...)
  } else {
    if (verbose) print('data type: pca')
    y <- pca2station(x,...)
    if (!is.null(attr(x,"pca"))) {
      fit <- attr(x,'pca')
      coredata(fit) <- coredata(attr(x,'fitted_values'))
      attr(y,'fitted_values') <- fit
      attr(y,'original_data') <- attr(x,'original_data')
    }
    y <- attrcp(attr(x,'original_data'),y)
    if (verbose) print(names(attributes(y)))
  }
  return(y)
}

#' @export as.station.list
as.station.list <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.station.list")
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

    ALL <- NULL
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
  if (attr(x[[i]],'type')== 'downscaled results') {
    class(y) <- class(x[[1]])
  } else {
    class(y) <- c("station",class(x[[1]]))
  }
  return(y)
}

#' @export as.station.field
as.station.field <- function(x,...,is=NULL,verbose=FALSE) {
  if(verbose) print("as.station.field")
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

#' @export as.station.spell
as.station.spell <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.station.spell")
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
  y <- attrcp(attr(x,'original_data'),y) # REB 2019-08-05
  attr(y,'history') <- history.stamp(x)
  return(y)
}

#' @export as.station.eof
as.station.eof <- function(x,...,ip=1:10,verbose=FALSE) {
  if(verbose) print("as.station.eof")
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
  y <- attrcp(attr(x,'original_data'),y) # REB 2019-08-05
  invisible(y)
}

#' Coerce input to a \code{station} object
#' 
#'
#' @param x the input object of class \code{dsensemble} 
#' @param ... other arguments
#' @param verbose a boolean; if TRUE print information about progress
#' 
#' @return a \code{station} object or a \code{dsensemble} \code{list} or \code{dsensemble} \code{station} object
#' 
#' @seealso as.station as.station.dsensemble.pca as.station.dsensemble.station DSensemble PCA
#' 
#' @export as.station.dsensemble
as.station.dsensemble <- function(x,...,verbose=FALSE) {
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
  y <- attrcp(attr(x,'original_data'),y) # REB 2019-08-05
  return(y)
}

#' Coerce input into a station object
#'
#' Transform a \code{dsensemble pca} object into a \code{dsensemble station} object
#' by extracting the results model wise and using the downscaled principle components
#' to reconstruct station time series.
#'
#' @param x a dsensemble pca object
#' @param is A list or data.frame providing space index, e.g. a list of longitude and
#' latitude range like list(lon=c(0,60), lat=c(35,60)).
#' @param ip selection of patterns in PCA or EOF (used for e.g. filtering the data)
#' @param verbose if TRUE print progress
#' @param \dots additional arguments
#'
#' @seealso as.station as.station.dsensemble DSensemble
#' 
#' @export as.station.dsensemble.pca
as.station.dsensemble.pca <- function(x,...,is=NULL,ip=NULL,verbose=FALSE) {
  if(verbose) print("as.station.dsensemble.pca")
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
    if (is.null(dim(attr(X$pca,'pattern')))) {
      dim(attr(X$pca,'pattern')) <- c(1,length(attr(X$pca,'pattern')))
    }
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
    if ( (is.list(S)) & (length(S)==length(locs)) ) names(S) <- locs
    attr(S,"unit") <- unit
    attr(S,"variable") <- param
    attr(S,"longname") <- longname
    attr(S,"aspect") <- "dsensemble.pca transformed to stations"
    attr(S,"history") <- history.stamp()
    invisible(S)
  }
}

#' @export as.station.dsensemble.station
as.station.dsensemble.station <- function(x,...,is=NULL,it=NULL,FUN='mean',verbose=FALSE) {

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

'# @export as.station.trajectory
as.station.trajectory <- function(x,...) {
  y <- trajectory2station(x,...)
  return(y)
}

'# @export as.station.events
as.station.events <- function(x,...) {
  y <- events2station(x,...)
  invisible(y)
}
