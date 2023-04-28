#' Coerce input to a \code{field} object
#' 
#' Transform an input object into the esd class \code{field}. 
#' \code{as.field} is an S3 method and will redirect to a fitting function depending on the output. 
#' The way in which the transformation is performed depends on the type of input data.
#' 
#' \code{as.field.events} redirects to \code{\link{events2field}}.
#' \code{as.field.trajectory} redirects to \code{\link{trajectory2field}}.
#' 
#' @seealso as.field.default as.field.zoo as.field.eof as.field.comb as.field.field as.field.ds as.field.station as.field.events as.field.trajectory as.field.dsensemble.eof as.field.matrix
#'
#' @importFrom stats median setNames
#'
#' @param x the input object
#' @param ... other arguments
#' 
#' @return a \code{field} object
#'
#' @examples
#' # how to generate a new field object.
#' year <- sort(rep(1991:2000,12))
#' month <- rep(1:12,length(1991:2000))
#' n <-length(year)
#' lon <- seq(-30,40,by=5); nx <- length(lon)
#' lat <- seq(40,70,by=5); ny <- length(lat)
#' # Time dimension should come first, space second.
#' y <- matrix(rnorm(nx*ny*n),n,nx*ny)
#' index <- as.Date(paste(year,month,1,sep="-"))
#' Y <- as.field(y,index=index,lon=lon,lat=lat,param="noise",unit="none")
#' map(Y)
#'
#' @export as.field
as.field <- function(x,...) UseMethod("as.field")

#' Coerce input to a \code{field} object
#' 
#' Transform a \code{zoo} object into a \code{field} object
#' 
#' @param x the input object of class \code{zoo} typically containing data from one or several measurement stations
#' @param lon longitude(s), a numerical or numerical vector
#' @param lat latitudes(s), a numerical or numerical vector
#' @param param short name of variable
#' @param unit unit of variable, e.g., 't2m'
#' @param longname long name of variable, e.g, 'temperature at 2m'
#' @param quality quality flag
#' @param src source of data
#' @param url url to website where data can be downloaded
#' @param reference reference describing data set
#' @param info additional information
#' @param calendar calendar type
#' @param greenwich a boolean; if TRUE center map on the Greenwich line (0E)
#' @param method method applied to data
#' @param type type of data
#' @param aspect aspect describing data, e.g., 'original', 'anomaly', 'climatology'
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{field} object
#' 
#' @seealso as.field
#' 
#' @exportS3Method
#' @export
as.field.zoo <- function(x,...,lon,lat,param,unit,
                         longname=NA,quality=NA,src=NA,url=NA,
                         reference=NA,info=NA,calendar='gregorian',
                         greenwich=TRUE, method= NA,type=NA,aspect=NA,
                         verbose=FALSE) {
  if(verbose) print("as.field.zoo")
  t <- index(x)
  if (length(year(x))!=1) {
    dyr <- median(diff(year(x)))#diff(year(x))[1]
    dmo <- median(diff(month(x)))#diff(month(x))[1]
    dda <- median(diff(day(x)))#diff(day(x))[1]
    timescale <- "annual"
    if (dmo>0)  timescale <- "month"
    if (dmo==3)  timescale <- "season"
    if (dda>0)  timescale <- "day"
    if (dyr==0 & dmo==0 & dda==0) timescale <- "sub-daily"
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

#' Coerce input to a \code{field} object
#' 
#' Transform an input object into the esd class \code{field}. The function first transforms the input object \code{x} into a \code{zoo} object (\code{zoo(x, order.by=index)}) and then applies \code{as.field.zoo} to obtain a \code{field} object.)  
#' 
#' @param x the input object
#' @param time index
#' @param lon longitude(s), a numerical or numerical vector
#' @param lat latitudes(s), a numerical or numerical vector
#' @param param short name of variable
#' @param unit unit of variable, e.g., 't2m'
#' @param longname long name of variable, e.g, 'temperature at 2m'
#' @param quality quality flag
#' @param src source of data
#' @param url url to website where data can be downloaded
#' @param reference reference describing data set
#' @param info additional information
#' @param calendar calendar type
#' @param greenwich a boolean; if TRUE center map on the Greenwich line (0E)
#' @param method method applied to data
#' @param type type of data
#' @param aspect aspect describing data, e.g., 'original', 'anomaly', 'climatology'
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other input arguments
#' 
#' @return a \code{field} object
#' 
#' @seealso as.field as.field.zoo zoo
#'
#' @examples
#' # how to generate a new field object.
#' year <- sort(rep(1991:2000,12))
#' month <- rep(1:12,length(1991:2000))
#' n <-length(year)
#' lon <- seq(-30,40,by=5); nx <- length(lon)
#' lat <- seq(40,70,by=5); ny <- length(lat)
#' # Time dimension should come first, space second.
#' y <- matrix(rnorm(nx*ny*n),n,nx*ny)
#' index <- as.Date(paste(year,month,1,sep="-"))
#' Y <- as.field(y,index=index,lon=lon,lat=lat,param="noise",unit="none")
#' map(Y)
#'
#' @exportS3Method
#' @export
as.field.default <- function(x,...,index,lon,lat,param,unit,
                             longname=NA,quality=NA,src=NA,url=NA,
                             reference=NA,info=NA,calendar='gregorian',
                             greenwich=TRUE, method=NA,type=NA,aspect=NA,
                             verbose=FALSE) {
  if(verbose) print("as.field.default")
  z <- zoo(x=x,order.by=index)
  x <- as.field.zoo(z,lon=lon,lat=lat,param=param,unit=unit,
                    longname=longname,quality=quality,src=src,url=url,
                    reference=reference,info=info,calendar=calendar,
                    greenwich=greenwich, method=method,type=type,
                    aspect=aspect, verbose=verbose)
  invisible(x)
}

#' @exportS3Method
#' @export
as.field.matrix <- function(x,...,index,lon,lat,param,unit,
                            longname=NA,quality=NA,src=NA,url=NA,
                            reference=NA,info=NA,calendar='gregorian',
                            greenwich=TRUE, method=NA,type=NA,aspect=NA,
                            verbose=FALSE) {
  x <- as.field.default(x=x,index=index,lon=lon,lat=lat,param=param,unit=unit,
                        longname=longname,quality=quality,src=src,url=url,
                        reference=reference,info=info,calendar=calendar,
                        greenwich=greenwich, method=method,type=type,
                        aspect=aspect, verbose=verbose)
  invisible(x)
}

#' Coerce input to a \code{field} object
#' 
#' Transform an input object into the esd class \code{field}. If the input is a combined field, redirect to \code{as.field.comb}, else return input. 
#' 
#' @param x the input object of class \code{field}
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{field} object
#' 
#' @seealso as.field as.field.comb
#' 
#' @exportS3Method
#' @export
as.field.field <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.field.field")
  if (inherits(x,'comb')) x <- as.field.comb(x,...)
  return(x)
}

#' Coerce input to a \code{field} object
#' 
#' Transform a combined field object into a \code{field} object, either dropping the appendices (if ip=NULL) or selecting one of the appended fields.
#' 
#' @param x the input object of class \code{field} \code{comb}
#' @param iapp a numerical; an index representing the appendix to extract
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{field} object
#' 
#' @seealso as.field
#' 
#' @exportS3Method
#' @export
as.field.comb <- function(x,...,iapp=NULL,verbose=FALSE) {
  if(verbose) print("as.field.comb")
  if (is.null(iapp)) {
    # Drop the appendend fields:
    n <- attr(x,'n.apps')
    for ( i in 1:n ) eval(parse(text=paste("attr(x,'appendix.",i,"') <- NULL",sep="")))
    attr(x,'n.apps') <- NULL
    class(x) <- class(x)[-1]
    return(x)
  } else {
    # Select one of the appended fields
    eval(parse(text=paste("y <- attr(x,'appendix.",iapp,"')",sep="")))
    attr(y,'longitude') <- attr(x,'longitude')
    attr(y,'latitude') <- attr(x,'latitude')
    attr(y,'history') <- history.stamp(x)
    return(y)
  }
}

#' Coerce input to a \code{field} object
#' 
#' Transform an \code{eof} object into the esd class \code{field}. If the input is a \code{dsensemble} object, it will be redirected to \code{as.field.dsensemble}. Otherwise the object is transformed using the function \code{eof2field}.
#' 
#' @param x the input object of class \code{eof}
#' @param iapp a numerical; an index representing the appendix to extract
#' @param anomaly a boolean; if TRUE return anomalies (climatology is in attribute 'mean') 
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{field} object
#' 
#' @seealso as.field eof2field EOF as.field.dsensemble
#' 
#' @exportS3Method
#' @export
as.field.eof <- function(x,...,iapp=NULL,anomaly=FALSE,verbose=FALSE) {
  if(verbose) print("as.field.eof")
  if (inherits(x,'dsensemble')) {
    y <- as.field.dsensemble(x,verbose=verbose,...)
  } else if (!inherits(x,'comb')) {
    y <- eof2field(x,verbose=verbose,...)
  } else {
    y <- as.eof(x,iapp=iapp,verbose=verbose)
    y <- eof2field(y,verbose=verbose,anomaly=anomaly,...)
  }
  return(y)
}

#' Coerce input to a \code{field} object
#' 
#' Transform \code{ds} object (output of the function \code{DS}) into the esd class \code{field}. 
#' 
#' @param x the input object of class \code{ds}
#' @param iapp a numerical; an index representing the appendix to extract
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{field} object
#' 
#' @seealso as.field as.field.eof DS
#' 
#' @exportS3Method
#' @export
as.field.ds <- function(x,...,iapp=NULL,verbose=FALSE) {
  if(verbose) print("as.field.ds")
  if (inherits(x,'eof')) {
    class(x) <- class(x)[-1]
    ## REB a few lines to catch cases where ds has not caught the comb-aspects.
    if (!is.null(iapp)) {
      if (!is.null(attr(x,'n.apps'))) {
        class(x)[length(class(x))+1]<-'comb'
        y <- as.field.eof(x,iapp=iapp,...)
        return(y)}
    } else {
      y <- as.field.eof(x,iapp=iapp,...)
      ## The residuals
      fit <- attr(x,'fitted_values')
      fit <- attrcp(attr(x,'eof'),fit)
      class(fit) <- class(attr(x,'eof'))
      attr(y,'fitted_values') <- fit
      attr(y,'original_data') <- attr(x,'original_data')
      attr(y,'calibration_data') <- attr(x,'calibration_data')
    }
  } else y <- NULL
  return(y)
}

#' Coerce input to a \code{field} object
#' 
#' Transform a \code{station} object (output of the function \code{DS}) into the esd class \code{field} using the function \code{regrid}. 
#' 
#' @param x the input object of class \code{ds}
#' @param lon a numerical vector of longitudes defining the grid of the output field
#' @param lat a numerical vector of latitudes defining the grid of the output field
#' @param nx number of grid points in the longitude direction - only used if lon=NULL
#' @param ny number of grid points in the latitude direction - only used if lat=NULL
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{field} object
#' 
#' @seealso as.field regrid
#' 
#' @exportS3Method
#' @export
as.field.station <- function(x,...,lon=NULL,lat=NULL,nx=30,ny=30,
                             verbose=FALSE) {
  if(verbose) print("as.field.station")
  if (is.null(lon)) lon <- seq(min(lon(x)),max(lon(x)),length=nx)
  if (is.null(lat)) lat <- seq(min(lat(x)),max(lat(x)),length=ny)
  y <- regrid(x,is=list(lon=lon,lat=lat),verbose=verbose)
  attr(y,'history') <- history.stamp(x)
  return(y)  
}

#' Coerce input to a \code{field} object
#' 
#' Transform a \code{dsensemble} \code{eof} object (output of the function \code{DS}) into the esd class \code{dsensemble} \code{field} \code{list}.
#' The function uses the downscaled principle components with the corresponding spatial EOF patterns and eigenvalues to calculate downscaled fields.
#' This is done for all selected ensemble members (see input \code{im}).
#' 
#' @param x the input object of class \code{ds}
#' @param is a list or data.frame providing a space index, e.g., station record or a lon(gitude) and lat(itude) range. If NULL include all.
#' @param ip a numerical or numerical vector with indices of the principle components to be included. If NULL include all.
#' @param im a numerical or numerical vector with indices of the ensemble members to be included. If NULL include all.
#' @param anomaly a boolean; if TRUE return anomalies (climatology is in attribute 'mean') 
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{dsensemble} \code{field} \code{list} object, i.e., a list where each element is a downscaled field corresponding to an ensemble member
#' 
#' @seealso as.field DSensemble EOF
#' 
#' @exportS3Method
#' @export
as.field.dsensemble <- function(x,...,is=NULL,ip=NULL,im=NULL,
                                    anomaly=FALSE,verbose=FALSE) {
  if (verbose) print('as.field.dsensemble')
  stopifnot(inherits(x,"dsensemble") & (inherits(x,"eof")|inherits(x,"pca")))
  if(inherits(x,"pca")) x <- as.eof(x,is=is,ip=ip,aggregate=FALSE,verbose=verbose,...)
  if (inherits(x,"field")) {
    invisible(x)
  } else {
    if (verbose) print('Extract the results model-wise')
    if(is.null(im)) {
      ix <- 3:length(x)
    } else {
      ix <- im[im>0 & im<(length(x)-2)] + 2
    }
    d <- apply(sapply(x[ix],dim),1,min)
    V <- array(unlist(lapply( x[ix],
                              function(x) coredata(x[1:d[1],1:d[2]]))),dim=c(d,length(ix)))
    if (is.null(ip)) {
      U <- attr(x$eof,'pattern')
      W <- attr(x$eof,'eigenvalues')
    } else {
      U <- attr(x$eof,'pattern')[,,ip]
      W <- attr(x$eof,'eigenvalues')[ip]
      V <- V[,ip,]
    }    
    d <- dim(U)
    dim(U) <- c(d[1]*d[2],d[3])
    S <- apply(V, 3, function(x) U %*% diag(W) %*% t(x))
    dim(S) <- c(dim(U)[1], dim(V)[1], dim(V)[3])
    if(!anomaly) {
      for (i in seq(1:dim(S)[1])) {
        S[i,,] <- S[i,,] + c(attr(x$eof,'mean'))[i]
      }
    }
    
    S <- aperm(S,c(3,2,1))
    S <- lapply(split(S, arrayInd(seq_along(S),dim(S))[,1]),
                array,dim=dim(S)[2:3])
    
    if (verbose) print('Set attributes')
    Y <- as.field(x$eof, verbose=verbose)
    gcms <- names(x)[ix]
    S <- setNames(S,gcms)
    for (i in seq_along(ix)) { 
      S[[i]] <- as.field(S[[i]],index=index(x[[ix[i]]]),
                         lon=attr(Y,"longitude"),lat=attr(Y,"latitude"),
                         param=varid(Y),unit=attr(Y,'unit'),
                         longname=paste('fitted',attr(Y,'longname')),
                         greenwich=attr(Y,'greenwich'),aspect='fitted')
      attr(S[[i]],'unitarea') <- attr(x,"unitarea")
      class(S[[i]]) <- class(Y)
    }
    if (!is.null(is)) S <- subset(S,is=is,verbose=verbose)
    class(S) <- c("dsensemble","field","list")
    S <- attrcp(x,S)
    if(!is.null(attr(x,"model_id"))) attr(S,"model_id") <- attr(x,"model_id")[ix]
    attr(S,"unit") <- attr(Y,"unit")
    attr(S,"variable") <- attr(Y,"variable")
    attr(S,"longname") <- attr(Y,"longname")
    attr(S,"aspect") <- "dsensemble.eof transformed to field"
    if(anomaly) {
      attr(S,"aspect") <- c(attr(S,"aspect"), "anomaly")
      attr(S,"mean") <- attr(x$eof,"mean")
    }
    attr(S,"history") <- history.stamp()
    invisible(S)
  }
}

#' @exportS3Method
#' @export
as.field.events <- function(x,...) {
  y <- events2field(x,...)
  return(y)
}

#' @exportS3Method
#' @export
as.field.trajectory <- function(x,...) {
  y <- trajectory2field(x,...)
  return(y)
}
