# Combines a group of station objects to form a set of time series
# Combines field objects; either in time (concatinates) or in space (mixing fields)
# Can be used to combine EOFS

#' Extension of rbind for field objects
#'
#' @param \dots input arguments
#' @param verbose if TRUE print information on progress
#' 
#' @exportS3Method
#' @export rbind.field
rbind.field <- function(...,verbose=FALSE) {
  if(verbose) print("rbind.field")
  if(verbose) print('note: the results inherit the metadata from the first argument')
  x <- list(...)
  y <- rbind.zoo(...)
  y <- attrcp(x[[1]],y)
  attr(y,'dimensions') <- c(attr(x[[1]],'dimensions')[1:2],length(index(y)))
  attr(y,'history') <- history.stamp(x[[1]])
  class(y) <- class(x[[1]])
  return(y)
}

#' Extension of rbind for field objects
#'
#' @param \dots input arguments
#' @param verbose if TRUE print information on progress
#'
#' @exportS3Method
#' @export cbind.field
cbind.field <- function(...,verbose=FALSE) {
  if(verbose) print("cbind.field")
  if(verbose) print('note: the results inherit the metadata from the first argument')
  x <- list(...)
  y <- cbind.zoo(...)
  y <- attrcp(x[[1]],y)
  attr(y,'dimensions') <- c(attr(x[[1]],'dimensions')[1:2],length(index(y)))
  attr(y,'history') <- history.stamp(x[[1]])
  class(y) <- class(x[[1]])
  return(y)
}


#' Combine
#' 
#' \code{combine} is a S3 method for combining esd objects, e.g. into groups of
#' stations, stations and eof object, or fields. The function is based on
#' \code{\link[zoo]{merge.zoo}}, and is also used to synchronise the esd objects.
#' 
#' For fields, \code{combine.field} is used to append different data sets, e.g.
#' for the purpose of computing common EOFs (seeo \code{\link{EOF}} or for
#' mixing fields (coupled EOFs).
#' 
#' For stations, \code{combine.station} can work tow ways: (1) to combine a set
#' of stations and group them into one data object; (2) combine series with
#' different monthly values for one specific site into one record for the
#' monthly data. E.g. January, February, ..., December months can be combined
#' into one complete series of monthly data.
#' 
#' For DS-results, \code{combine.ds} is based on \code{combine.station}, but
#' also takes care of the additional meta data (the original series and
#' predictor patterns). For instance, this method can combine seperate
#' downscaled results for each calendar months at a single location into one
#' complete time series.
#' 
#' \code{g2dl} transform objects between grid starting at the grenwich
#' (\code{greenwich=TRUE}) and the data line (\code{greenwich=FALSE}).
#' 
#' \code{sp2np} re-arranges field objects accroding to a grid going from 90S
#' (South Pole) to 90N (Noth Pole) for \code{SP2NP=TRUE}. Otherwise, the object
#' is arranged from 90N to 90S.
#' 
#' Other operations, such as \code{c(...)}, \code{rbind(...)} (combine along
#' the time dimension), and \code{cbind(...)} (combine along the space
#' dimension) also work.
#' 
#' @aliases combine combine.default combine.station combine.stations
#' combine.zoo combine.ds combine.ds.comb combine.ds.station
#' combine.ds.station.eof combine.ds.station.field combine.ds.pca
#' combine.station.eof combine.field combine.filed.station
#' combine.station.field sp2np combine.trajectory combine.events
#' combine.field.station combine.list
#'
#' @param x station, eof, or field object
#' @param all See \code{link{merge.zoo}}
#' @param orig.format TRUE: the result will the formatted the same way as the
#' input.
#' @param dimension Which dimension to combine - in time or in space
#' @param approach How to combine
#' @param greenwich TRUE: center map on the Greenwich line (0E)
#' @param SP2NP TRUE: order from south pole (bottom of plot) to north pole (top
#' of plot)
#' @param ignore List of attributes to ignore.
#'
#' @return A field object
#'
#' @keywords utilities
#'
#' @examples 
#' T2m_DNMI <- t2m.DNMI(lon=c(-40,40),lat=c(30,70))
#' T2m_NorESM <- t2m.NorESM.M(lon=c(-40,40),lat=c(30,70))
#' 
#' # Combine in time to compute common EOFs:
#' X <- combine(T2m_DNMI,T2m_NorESM)
#' ceof <- EOF(X,it="Jan")
#' plot(ceof)
#' 
#' # Use combine to synchronise field and station data:
#' data("Oslo")
#' y <- combine(Oslo,T2m_DNMI)
#' plot(y$y)
#' 
#' @export combine
combine <- function(...) UseMethod("combine")

#' @exportS3Method
#' @export combine.default
combine.default <- function(x=NULL,y=NULL,...,all=FALSE,orig.format=TRUE,verbose=FALSE) {
  if(verbose) print("combine.default")
  stopifnot(!missing(x))
  if (is.null(x)) return(y)
  if (is.null(y)) return(x)
  if ( !zeros(c(inherits(x,c("station","zoo"), which=TRUE),
                inherits(y,c("station","zoo"), which=TRUE))) ) { 
    X <- combine.station(x=x,y=y,verbose=verbose)
  } else if ( !zeros( c(inherits(x,c("station","zoo"), which=TRUE),
                        inherits(y,c("eof","zoo"), which=TRUE)) ) |
              !zeros( c(inherits(y,c("station","zoo"), which=TRUE),
	                inherits(x,c("eof","zoo"), which=TRUE)) ) ) {
    X <- combine.station.eof(x=x,y=y,all=all,orig.format=orig.format,verbose=verbose)
  } else if ( !zeros( c(inherits(x,c("station","zoo"), which=TRUE),
                        inherits(y,c("field","zoo"), which=TRUE)) ) |
              !zeros( c(inherits(y,c("station","zoo"), which=TRUE),
	                inherits(x,c("field","zoo"), which=TRUE)) ) ) {
    X <- combine.station.field(x=x,y=y,all=all,orig.format=orig.format,verbose=verbose)
  } else { 
    print("combine.default - don't know what to do :-(")
    X <- NULL
  }
  attr(X,'history') <- history.stamp(x)
  invisible(X)
}

# combine.station can be used to either combine a group of stations into
# one data object or time series from one stations with different monthly
# values into one time series with all months
#' @exportS3Method
#' @export combine.station
combine.station <- function(...,all=TRUE,verbose=FALSE) {
  if(verbose) print("combine.station")
  cl <- as.list(match.call())
  args <- list(...)
  n <- length(args)
  stid <- NULL
  #stid <- rep(NA,n)
  allstations <- TRUE
  for (i in 1:n) {
    Z <- args[[i]]
    if (inherits(Z,'station')) {
      stid <- c(stid,attr(Z,'station_id'))
      #stid[i] <- attr(Z,'station_id')#[i]
    } else {
      allstations <- FALSE
      stid <- c(stid,NA)
    }
  }
  #print(allstations)
  if (allstations) {
    ns <- length(table(stid))
    # If only one site, then combine the months into one series, otherwise
    # group the stations into one multivariate object
    if (ns==1) {
      X <- combine.station.month(...,verbose=verbose)
    } else {
      X <- combine.stations(...,all=all,verbose=verbose)
    }
    class(X) <- class(Z)
  } else {
    X <- combine.default(...,all=all,verbose=verbose)
  }
  attr(X,'history') <- history.stamp(X)
  invisible(X)
}

#' @exportS3Method
combine.station.month <- function(...,verbose=FALSE) {
  if(verbose) print("combine.station.month")
  cl <- as.list(match.call())
  args <- list(...)
  #print(summary(args)); print(length(args))
  z <- args[[1]]
  #plot.zoo(z); dev.new()
  Z <- merge.zoo(...,all=TRUE)
  #plot(Z);str(Z)
  X <- zoo(rowMeans(coredata(Z),na.rm=TRUE),order.by=as.Date(index(Z)))
  #print(summary(X))
  #print(names(attributes(z)))
  X <- attrcp(z,X,ignore=c('mean','calibration_data','fitted_values',
                           'original_data','aspect'))
  class(X) <- class(z)
  attr(X,'history') <- history.stamp(X)
  invisible(X)
}

#' @exportS3Method
#' @export combine.zoo
combine.zoo <- function(...,verbose=FALSE) {
  if(verbose) print("combine.zoo")
  Z <- merge.zoo(...,all=TRUE)
  #plot(Z);str(Z)
  invisible(Z)
}

# combine.stations is used to combine a group of stations into one object
#' @exportS3Method
#' @export combine.stations
combine.stations <- function(...,all=TRUE,verbose=FALSE) {
  if(verbose) print("combine.stations")
  # If either of the arguments is NULL, then return the x -
  # useful for looping
  
  # REB: rewritten to work with merge and long lists of stations
  
  cl <- as.list(match.call())
  #str(cl)
  args <- list(...)
  #if (verbose) print(args)
  #str(args)
  X <- merge.zoo(...,all=all)
  #plot(X)
  #str(X)
  #print(length(args))
  n <- length(args)
  cls <- class(args[[1]])
  
  loc <- NULL; cn <- loc; ID <- NULL; unit <- loc
  lon <- ID; lat <- ID; alt <- ID; param <- loc; lname <- loc
  src <- loc; qlty <- loc; url <- loc; ref <- loc
  info <- loc; ele <- asp <- ID; th <- NULL; thu <- NULL
  if (verbose) print('Organise the attributes')
  for (i in 1:n) {
    Z <- args[[i]]
    #attr <- softattr(Z)
    loc <- c(loc,attr(Z,'location'))
    cn  <- c(cn,attr(Z,'country'))
    ID <- c(ID,attr(Z,'station_id'))
    unit <- c(unit,attr(Z,'unit'))
    lon <- c(lon,attr(Z,'longitude'))
    lat <- c(lat,attr(Z,'latitude'))
    alt <- c(alt,attr(Z,'altitude'))
    param <- c(param,attr(Z,'variable'))
    lname <- c(lname,attr(Z,'longname'))
    src <- c(src,attr(Z,'source'))
    qlty <- c(qlty,attr(Z,'quality'))
    url <- c(url,attr(Z,'URL'))
    ref <- c(ref,attr(Z,'reference'))
    info <- c(info,attr(Z,'info'))
    ele <- c(ele,attr(Z,'element'))
    asp <- c(asp,attr(Z,'aspect'))
    th <- c(th,attr(Z,'threshold'))
    thu <- c(thu,attr(Z,'threshold.unit'))
  }
  if (dim(X)[2]==length(loc)) colnames(X) <- loc
  if(any(is.na(colnames(X)))) colnames(X) <- paste("X",seq(ncol(X)),sep=".")
  attr(X,'location') <- loc
  attr(X,'country') <- cn
  attr(X,'station_id') <- ID
  attr(X,'longitude') <- lon
  attr(X,'latitude') <- lat
  attr(X,'altitude') <- alt
  attr(X,'variable') <- param
  attr(X,'longname') <- lname
  attr(X,'unit') <- unit
  attr(X,'aspect') <- asp      ## AM
  attr(X,'source') <- src
  attr(X,'element') <- ele
  attr(X,'quality') <- qlty
  attr(X,'threshold') <- th
  attr(X,'threshold.unit') <- thu
  attr(X,'URL') <- url
  attr(X,'history') <- history.stamp(Z)
  #attr(X,'date-stamp') <- date()
  attr(X,'reference') <- ref
  attr(X,'info') <- info
  class(X) <- cls
  invisible(X)
}

#' @exportS3Method
#' @export combine.ds
combine.ds <- function(...,all=TRUE,verbose=FALSE) {
  if(verbose) print("combine.ds")
  cl <- as.list(match.call())
  ##str(cl)
  args <- list(...)
  #str(args)
  z <- args[[1]]
  n <- length(args)
  #print(n); print(summary(z)); print(names(attributes(z[[1]])))
  #print(class(z))
  if (n > 1) {
    #  Arguments consisting of several objects to be combined
    #print("Several stations")
    if (inherits(z,'pca')) {
      X <- combine.ds.pca(...,all=all) 
    } else {
      X <- combine.station(...,all=TRUE)
    }
  } else if (n==1) {
    #print("One station")
    #print(class(z)); print(attr(z[[1]],'type')); print(summary(z))
    if ( (attr(z[[1]],'type')=="downscaled results") & is.list(z) ) {
      # The output of DS.t2m.month.field: one list with several months
      # The main results looks like a station object:
      #    X <- as.station.ds(z)
      X <- as.station.list(z)
      #X <- rbind(z[[1]],z[[2]],z[[3]],z[[4]]) # AM 27.11.2014 quick fix must be included in as.station.list ..
      #print("HERE"); print(names(X))
      #plot(X,col="red")
      
      # The original data & fitted values:
      #print("combine the training data ...")
      m <- length(z)
      for (i in 1:m) {
        #model <- attr(z[[i]],'model')
        y <- attr(z[[i]],'fitted_values')
        x <- attr(z[[i]],'calibration_data')
        x0 <- attr(z[[i]],'original_data')
        xval <- crossval(z[[i]])
        #print(summary(y - y0))
        attr(y,'station_id') <- attr(X,'station_id')
        attr(x,'station_id') <- attr(X,'station_id')
        attr(x0,'station_id') <- attr(X,'station_id')
        eval(parse(text=paste("y.",i," <- y",sep="")))
        eval(parse(text=paste("x.",i," <- x",sep="")))
        eval(parse(text=paste("x0.",i," <- x0",sep="")))
        eval(parse(text=paste("xv.",i," <-  zoo(xval)",sep="")))
        eof <- attr(z[[i]],'eof')  # REB 12.02.2014
        if (i==1) {
          argsx <- 'x.1'
          argsy <- 'y.1'
          argsx0 <- 'x0.1'
          argsxY <- 'xv.1'
          pattern <- attr(z[[1]],'pattern')
          d <- dim(pattern)
          dim(pattern) <- c(1,d[1]*d[2])
          # REB - also capture the diagnostics - from the commone EOF
          if (!is.null(attr(eof,'diagnose')))
            diag <- list(s.1=attr(eof,'diagnose'))
        } else {
          argsx <- paste(argsx,",x.",i,sep="")
          argsy <- paste(argsy,",y.",i,sep="")
          argsx0 <- paste(argsx0,",x0.",i,sep="") ## AM 27.11.2014 argsx replaced by argsx0
          argsxY <- paste(argsxY,",xv.",i,sep="") ## AM 27.11.2014 argsx replaced by argsx0
          patt <- attr(z[[i]],'pattern')
          dim(patt) <- c(1,d[1]*d[2])
          pattern <- rbind(pattern,patt)
          # REB - also capture the diagnostics - from the commone EOF
          if (!is.null(attr(eof,'diagnose')))
            eval(parse(text=paste("diag$s.",i," <- attr(eof,'diagnose')",sep="")))
        }
      }
      #print(argsx)
      #str(x.1); print(class(x.1))
      #print("---")
      X0 <- NULL
      cline1 <- parse(text=paste("X0 <- rbind(",argsx,")")) ## cline1 <- parse(text=paste("X0 <- merge(",argsx,",all=TRUE)"))
      #print(cline1)
      eval(cline1)
      #lines(X0)    
      
      #print("here"); print(argsy)
      Y <- NULL
      cline2 <- parse(text=paste("Y <- rbind(",argsy,")")) ## AM cline2 <- parse(text=paste("Y <- merge(",argsy,",all=TRUE)"))
      #print(cline2)
      eval(cline2)
      #lines(Y,col="pink")
      
      #print("HERE")
      X0.0 <- NULL
      cline3 <- parse(text=paste("X0.0 <- rbind(",argsx0,")")) ## cline3 <- parse(text=paste("X0.0 <- merge(",argsx0,",all=TRUE)"))
      #print(cline3)
      eval(cline3)
      #lines(Y,col="pink")
      
      # Cross-validation
      xY <- NULL
      #str(x.1); print(class(xv.1))
      cline4 <- parse(text=paste("xY <- rbind(",argsxY,")")) ## cline4 <- parse(text=paste("xY <- merge(",argsxY,",all=TRUE)"))
      eval(cline4)
      
      attr(X,'calibration_data') <- X0
      attr(X,'fitted_values') <- Y
      attr(X,'evaluation') <- xY
      attr(X,'original_data') <- X0.0
      attr(X,'diagnose') <- diag
      #print(dim(pattern))
      dim(pattern) <- c(m,d[1],d[2])
      attr(pattern,'dimensions') <- c('month','longitude','latitude')
      attr(pattern,'month') <- month.abb
      #            attr(pattern,'longitude') <- attr(z[[1]],'longitude')
      #            attr(pattern,'latitude') <- attr(z[[1]],'latitude')
      attr(pattern,'longitude') <- lon(attr(z[[1]],'pattern'))
      attr(pattern,'latitude') <- lat(attr(z[[1]],'pattern'))
      
      attr(X,'pattern') <- pattern
      #print(class(X))
    } else X <- NULL 
  }
  
  #print(names(attributes(z)))
  #str(X)
  # Copy the attributes from z, but not those already used for station
  # objects, as those are taken care of in combine.station
  #print("Copy attributes")
  ## X <- attrcp(z[[1]],X,ignore=c('location','country','station_id','longitude',
  ##                  'latitude','altitude','variable','longname','unit','evaluation',
  ##                  'source','element','quality','URL','reference','info','names'))
  #print("...")
  attr(X,'type') <- 'downscaled results'
  attr(X,'history') <- history.stamp(X)
  attr(X,'old_class') <- class(z[[1]])
  #print("exit combine.ds")
  #print("HERE"); print(names(X))
  invisible(X)
}

#' @exportS3Method
#' @export 
combine.list <- function(...,all=TRUE,verbose=FALSE) {
  if(verbose) print("combine.list")
  args <- list(...)
  z <- args[[1]]
  #print(class(z[[1]]))
  if (inherits(z[[1]],'comb'))
    y <- combine.ds.comb(...,all=all) else
      if (inherits(z[[1]],'ds'))
        y <- combine.ds(...,all=all) else
          y <- NULL
  return(y)
}

#' @exportS3Method
#' @export combine.ds.comb
combine.ds.comb <- function(...,all=TRUE,verbose=FALSE) {
  if(verbose) print("combine.ds.comb")
  cl <- as.list(match.call())
  #str(cl)
  args <- list(...)
  #print(args)  # [[1]][[1:4]]
  #print(summary(args))
  z <- args[[1]]
  #print(summary(z)); print(class(z)); print(length(z))
  k <- length(args)
  m <- length(...) 
  #print(k); print(m); print(summary(z)); print(names(attributes(z[[1]])))
  #print(class(z))
  
  if (k > 1) {
    #  Arguments consisting of several objects to be combined
    #print("Several stations")
    if (inherits(z,'pca'))
      X <- combine.ds(...,all=all)
  } else {
    #print("One station")
    #print(class(z))
    X <- combine.ds(...,all=all)
    n <- attr(z[[1]],'n.apps')
    
    # Combine the different appendixes for the different list objects:
    # appendix.1 in each list object is combined with corresponding
    # appendix.1 in the others and so on.
    #print("appendix.x")
    for (i in 1:n) {
      for (j in 1:m) {
        Z <- z[[j]]
        #print(summary(Z))
        #print(mean(attr(Z,paste('appendix.',i,sep=''))))
        x <- attr(Z,paste('appendix.',i,sep=''))
        #attr(x,'station_id') <- attr(Z,'station_id')
        #attr(x,'aspect') <- 'original'
        #class(x) <- c('station',class(Z))
        #print(summary(x))
        class(x) <- 'zoo'
        eval(parse(text=paste("x.",j," <- x",sep="")))
        if (j==1) {
          argsx <- 'x.1'
        } else {
          argsx <- paste(argsx,",x.",j,sep="")
        }
      }
      
      #print(argsx)
      #print("combine.station -> appendix.x")
      #plot(c(x.1,x.2,x.3,x.4)); dev.new()  
      #      cline1 <- parse(text=paste("attr(X,'appendix.",i,
      #                                 "') <- combine.zoo(",argsx,")",sep=""))
      cline1 <- parse(text=paste("attr(X,'appendix.",i,"') <- c(",argsx,")",sep=""))  
      #print(cline1)
      eval(cline1)
      #lines(X0)    
      #print("exit combine.ds.comb") 
    }
  }
  invisible(X)
}

# not exported
combine.ds.station <- function(...,all=TRUE,verbose=FALSE) {
  if(verbose) print("combine.ds.station") 
  # Combine downscaled station records. Use combine.station for
  # combining the station values: either as a group of stations
  # or combine different months for one station into one series
  
  #print("combine.ds.station")
  cl <- as.list(match.call())
  #str(cl)
  args <- list(...)
  #str(args)
  X <- combine.station(...,all=all)
  attr(X,'type') <- 'ESD_result'
  attr(X,'history') <- history.stamp(X)
  
  invisible(X)
}

# not exported
combine.ds.pca <- function(...,all=TRUE,verbose=FALSE) {
  if(verbose) print("combine.ds.pca")  
  # Combine downscaled PCA: i.e. the different principal components
  # Assmue that the pattern is the same for all these
  #print("combine.ds.pca")
  cl <- as.list(match.call())
  #str(cl)
  args <- list(...)
  #str(args)
  X <- combine.ds.station(...,all=all)
  attr(X,'pattern') <- attr(args[[1]],'pattern')
  attr(X,'eigenvalues') <-  attr(args[[1]],'eigenvalues')
  attr(X,'mean') <-  attr(args[[1]],'mean')
  attr(X,'sum.eigenv') <- attr(args[[1]],'sum.eigenv')
  attr(X,'tot.var') <- attr(args[[1]],'tot.var')
  X <- attrcp(args[[1]],X)
  attr(X,'history') <- history.stamp(X)
  class(X) <- class(args[[1]])
  invisible(X)
}

#' @exportS3Method
#' @export combine.station.eof
combine.station.eof <- function(x=NULL,y=NULL,...,all=FALSE,orig.format=TRUE,verbose=FALSE) {
  if(verbose) print("combine.station.eof")
  if (is.null(x)) return(y)
  if (is.null(y)) return(x)
  
  #print("combine.station.eof")
  #print(dim(y)); print(dim(x))
  
  # Keep track of which is an eof object and which is a station record:
  if ( inherits(x,c("station")) & inherits(y,c("eof"))) {
    yy <- x
    x <- y
    y <- yy
  } 
  
  #str(y); str(x)
  clsx <- class(x)
  clsy <- class(y)
  index(y) <- as.Date(index(y))
  #print(class(x))
  index(x) <- as.Date(index(x))
  #print("combine: dates")
  #print(index(x)[1]); print(index(y)[1]); 
  # REB 20.08.2013: added colnames to x & y before merge to keep track of
  # y & x.
  
  dy <- dim(y)
  if (is.null(dy)) dy <- c(length(y),1)
  if (dy[2]>1) colnames(y) <- paste("y",1:dy[2],sep=".")
  dx <- dim(x)
  if (is.null(dx)) dx <- c(length(x),1)
  if (dx[2]>1) colnames(x) <- paste("x",1:dx[2],sep=".") 
  
  if (orig.format) {
    comb <- merge(x,y,all=all)
    vars <- tolower(names(comb))
    ys <- vars[grep('y',vars)]
    Xs <- vars[grep('x',vars)]
    ix <- is.element(vars,Xs)
    iy <- is.element(vars,ys)
    
    XX <- zoo(coredata(comb[,ix]),order.by=index(comb))
    yy <- zoo(comb[,iy],order.by=index(comb))
    
    XX <- attrcp(x,XX,ignore='names')
    yy <- attrcp(y,yy,ignore='names')
    
    
    # REB 29.04.2014
    
    #   nattr1 <- softattr(y)
    #   for (i in 1:length(nattr1))
    #     attr(yy,nattr1[i]) <- attr(y,nattr1[i])
    #   nattr2 <- softattr(x)
    #   for (i in 1:length(nattr2))
    #     attr(XX,nattr2[i]) <- attr(x,nattr2[i])
    
    clsy -> class(yy)
    clsx -> class(XX)
    
    if (!is.null(attr(yy,'standard.error'))) {
      sterr <- attr(yy,'standard.error')
      sterr <- matchdate(sterr,yy)
      attr(yy,'standard.error') <- coredata(sterr)
    }
    X <- list(y=yy,X=XX)
  } else {
    # REB 29.04.2014
    # This option introdices an additional index in the PCs and
    # and additional uniform pattern.
    comb <- merge(y,x,all=all)
    
    d <- dim(y)
    if (is.null(d)) {
      s <- sd(y,na.rm=TRUE)
      z <- ( coredata(y) - mean(y,na.rm=TRUE) )/s
    } else {
      s <- apply(y,2,sd,na.rm=TRUE)
      m <- apply(y,2,mean,na.rm=TRUE)
      z <- ( coredata(y) -  m )/s
    }
    coredata(y) <- z
    comb <- merge(y,x,all=all)
    
    X <- comb
    #    nattr1 <- softattr(x)
    #    for (i in 1:length(nattr1))
    #      attr(X,nattr1[i]) <- attr(x,nattr1[i])
    X <- attrcp(x,X,ignore='names')
    attr(X,'pattern') <- rbind(matrix(rep(s,d[1]*d[2]),d[1],d[2]),
                               attr(X,'pattern'))
    attr(X,'X.attributes') <- attributes(y)
    class(X) <- class(x)
  }
  #print(summary(combined))
  attr(X,'history') <- history.stamp(X)
  invisible(X)
}

#' @exportS3Method
#' @export combine.station.field
combine.station.field <- function(x=NULL,y=NULL,...,all=FALSE,orig.format=TRUE,verbose=FALSE) {
  if(verbose) print("combine.station.field")
  X <- combine.field.station(x=y,y=x,all=all,orig.format=orig.format,verbose=verbose)
  invisible(X) 
}

#' @exportS3Method
#' @export combine.field.station
combine.field.station <- function(x=NULL,y=NULL,...,all=FALSE,
                                  orig.format=TRUE,verbose=FALSE) {
  if (verbose) print("combine.field.station")
  if (is.null(x)) return(y)
  if (is.null(y)) return(x)
  attr(x,'dimnames') <- NULL
  # Keep track of which is an eof object and which is a station record:
  swapped <- FALSE
  if ( inherits(x,c("station")) & inherits(y,c("field"))) {
    yy <- x
    x <- y
    y <- yy
    swapped <- TRUE
  } 
  if(verbose) print(swapped)
  clsx <- class(x)
  clsy <- class(y)
  index(y) <- as.Date(index(y))
  if (verbose) print(class(X)) 
  index(x) <- as.Date(index(x))
  if (length(dim(y))==2) {
    colnames(y) <- paste("y",1:dim(y)[2],sep=".") 
    if (length(dim(x))==2) {
      colnames(x) <- paste("x",1:dim(x)[2],sep=".")
    }
  }
  comb <- merge(x,y,all=all)
  if (verbose) print(summary(comb))
  
  if (orig.format) {
    vars <- names(comb)
    ys <- vars[grep('y',vars)]
    Xs <- vars[grep('x',vars)]
    ix <- is.element(vars,Xs)
    iy <- is.element(vars,ys)
    if (verbose) print(c(sum(ix),sum(iy)))
    
    XX <- zoo(coredata(comb[,ix]),order.by=index(comb))
    yy <- zoo(coredata(comb[,iy]),order.by=index(comb))
    names(XX) <- Xs
    names(yy) <- ys
    
    XX <- attrcp(x,XX,ignore='names')
    yy <- attrcp(y,yy,ignore='names')
    
    attr(XX,'dimensions') <- attr(x,'dimensions')
    
    clsx -> class(XX)
    clsy -> class(yy)
    if (swapped) {
      X <- list(y=XX,X=yy)
    } else {
      X <- list(y=yy,X=XX)
    }
  } else {   
    X <- comb
    X <- attrcp(x,X,ignore='names')
    attr(X,'X.attributes') <- attributes(y)
    class(X) <- c("comb",clsx)
  }
  attr(X,'history') <- history.stamp(X)
  invisible(X)
}
#' @exportS3Method
#' @export combine.field
combine.field <- function(x=NULL,y=NULL,...,all=FALSE,dimension="time",
                          approach="field",orig.format=TRUE,verbose=FALSE) {
  if (verbose) print(paste("combine.field",approach))
  attr(x,'dimnames') <- NULL
  attr(y,'dimnames') <- NULL
  if (inherits(y,'station')) {
    xy <- combine.field.station(x=x,y=y,orig.format=orig.format)
    return(xy)
  }
  stopifnot(!missing(x),inherits(x,'field'),inherits(y,'field'))
  clsx <- class(x); clsy <- class(y)
  hst <- c(attr(x,'history'),attr(y,'history'))
  if ( (unit(x) %in% c('hPa','mbar')) & (unit(y)=='Pa')) {
    if (verbose) print('Resetting unit of x: hPa -> Pa')
    coredata(x) <- 100*coredata(x)
    attr(x,'unit') <- 'Pa'
  } 
  if ( (unit(y) %in% c("hPa","mbar")) & (unit(x)=="Pa")) {
    if (verbose) print('Resetting unit of y: hPa -> Pa')
    coredata(y) <- 100*coredata(y)
    attr(y,'unit') <- 'Pa'
  }
  if ( (unit(x)=='m') & (unit(y)=='mm/day')) {
    if (verbose) print('Resetting unit of x: m -> mm/day')
    coredata(x) <- 1000*coredata(x)
    attr(x,'unit') <- 'mm'
    # if (clsy[2] %in% c('month','season','annual')) {
    #   if (clsy[2]=='month') coredata(y) <- 30*coredata(y)
    #   if (clsy[2]=='season') coredata(y) <- 90*coredata(y)
    #   if (clsy[2]=='annual') coredata(y) <- 365.25*coredata(y)
    #   attr(y,'unit') <- 'mm'
    # }
  } 
  if ( (unit(x)=='m') & (unit(y)=='kg m-2 s-1')) {
    if (verbose) print('Resetting unit of x: m -> mm')
    coredata(x) <- 1000*coredata(x)
    attr(x,'unit') <- 'mm'
    if (verbose) print('Resetting unit of y: m -> mm')
    # if (clsy[2] %in% c('day','month','season','annual')) {
    #   if (clsy[2]=='day') coredata(y) <- 3600*24*coredata(y)
    #   if (clsy[2]=='month') coredata(y) <- 30*3600*24*coredata(y)
    #   if (clsy[2]=='season') coredata(y) <- 90*3600*24*coredata(y)
    #   if (clsy[2]=='annual') coredata(y) <- 365.25*3600*24*coredata(y)
    #   attr(y,'unit') <- 'mm'
    # }
    attr(x,'unit') <- 'mm/day'
  } 
  if (unit(x)=='mbar' & unit(y)=='hPa') attr(x,'unit') <- 'hPa'
  if (unit(y)=='mbar' & unit(x)=='hPa') attr(y,'unit') <- 'hPa'
  if (unit(x)=='deg C') attr(x,'unit') <- 'degC'
  if (unit(y)=='deg C') attr(y,'unit') <- 'degC'
  if (unit(x) != unit(y)) print(paste('Warning - different units:',unit(x),unit(y)))
  
  if (missing(y)) return(x)
  
  ## Check the scales/units
  sx <- mean(coredata(x[1,]),na.rm=TRUE)
  sy <- mean(coredata(y[1,]),na.rm=TRUE)
  test.ratio <- try(abs(log(sx/sy)/log(10)))  ## Needed because some CMIP6 data files are not well conformed...
  if ( (inherits(test.ratio,'try-error')) | (!is.finite(test.ratio)) ) test.ratio <- 99
  if (test.ratio > 2) {
    print(paste('combine.field detected scale issues - sx=',round(sx),'sy=',
                round(sy),esd::unit(x)[1],esd::unit(y),src(y)))
    warning(paste('combine.field detected scale issues - sx=',round(sx),'sy=',
                  round(sy),esd::unit(x)[1],esd::unit(y),src(y)))
  }
    
  dimension <- tolower(dimension)
  approach <- tolower(approach)
  ## Make sure that both fields are presented on the same gid
  if (verbose) print (paste('Greenwich dateline - X:',attr(x,'greenwich'),'Y:',attr(y,'greenwich')))
  # x <- g2dl(x,greenwich=attr(y,'greenwich'))
  # if (verbose) print(rbind(range(lon(x)),range(lon(y))))
  x <- g2dl(x, greenwich = FALSE)
  y <- g2dl(y, greenwich = FALSE)
  x <- sp2np(x)
  y <- sp2np(y)
  
  if (sum(is.element(dimension,"time"))) {
    # Combine the two gridded data sets along the time axis:
    # The y data set is regridded onto the grid of the former.
    
    # Organise the grid coordinates:
    #x1 <- attr(x,'longitude');  nx1 <- length(x1)
    #y1 <- attr(x,'latitude');   ny1 <- length(y1)
    #x2 <- attr(y,'longitude'); nx2 <- length(x2)
    #y2 <- attr(y,'latitude');  ny2 <- length(y2)
    #d1 <- attr(x,'dimensions')
    #d2 <- attr(y,'dimensions')
    #xy2 <- rep(x2,ny2); yx2 <- sort(rep(y2,nx2))
    
    # Work with the raw matrices rather than zoo class to avoid surrises;-)
    #    X <- coredata(x); Y <- coredata(y)
    #    Z <- matrix(rep(NA,d2[3]*d1[1]*d1[2]),d2[3],d1[1]*d1[2])
    
    # Only use the same region as the previous:
    # y <- subset(y,is=list(range(x1),range(y1)))
    if (verbose) print("combine.field after subset:")
    #Z <- regrid(y,is=list(x1,y1))
    Z <- regrid(y,is=x)
    
    iv <- function(x) return(sum(is.finite(x))) # KMP 2016-03-16
    ngood <- apply(coredata(x),2,iv) # KMP 2016-03-16
    maskna <- ngood==0 # KMP 2016-03-16
    #maskna <- !is.finite(colMeans(coredata(x))) # REB 2016-02-18 mask out same missing
    if (sum(maskna)>0) {                                  # data as in the first field
      z <- coredata(Z)
      z[,maskna] <- NA
      coredata(Z) <- z
    }                     
    #attr(Z,'dimensions') <- c(d1[1],d1[2],d2[3])
    #mostattributes(Z) <- attributes(y)
    #class(Z) <- class(y)
    
    # Add the regridded data as an attribute, but change class
    if (is.null(attr(x,'n.apps'))) {
      n.app <- 1 
    } else {
      n.app <- attr(x,'n.apps') + 1
    }
    
    attr(x,paste('appendix.',n.app,sep='')) <- Z
    attr(x,'n.apps') <- n.app
    X <- x
    if (sum(is.element(clsx,"comb"))==0) {
      class(X) <- c("comb",clsx) 
    } else {
      class(X) <- clsx
    }
    
  }
  if (sum(is.element(dimension,"space"))) {
    # This option synchronises the data and expand the spatial grid
    # to include both data sets.
    #print(length(names(x[1,]))); print(length(x[1,])); print(names(x)[1:10])
    #xnm <- paste("x",1:length(x[1,]),sep=".") 
    #ynm <- paste("y",1:length(y[1,]),sep=".") 
    #names(x) <- xnm
    #names(y) <- ynm
    # REB 20.08.2013: added colnames to x & y before merge to keep track of
    # y & x.    
    colnames(y) <- paste("y",1:dim(y)[2],sep=".")
    colnames(x) <- paste("x",1:dim(x)[2],sep=".")
    comb <- merge(x,y,all=all)
    nt <- length(index(comb))
    if (verbose) str(comb)
    if (orig.format) {
      vars <- names(comb)
      if (verbose) print("After merge - original format")
      
      ys <- vars[grep('y',vars)]
      Xs <- vars[grep('x',vars)]
      ix <- is.element(vars,Xs)
      iy <- is.element(vars,ys)
      if (verbose) {print(paste(sum(ix),"xs and",sum(iy),"ys"))
        print(dim(comb))}
      XX <- zoo(coredata(comb[,ix]),order.by=index(comb))
      yy <- zoo(coredata(comb[,iy]),order.by=index(comb))
      names(yy) <- ys
      names(XX) <- Xs
      
      if (verbose) print("set class")
      #print(clsx); print(clsy)
      clsx -> class(XX)
      clsy -> class(yy)
      if (verbose) {print("add attributes")
        print(names(attributes(y)))}
      #nattr1 <- softattr(y)
      attr(y,'dimnames') <- NULL
      attr(x,'dimnames') <- NULL
      attr(yy,'dimnames') <- NULL
      attr(XX,'dimnames') <- NULL
      #for (i in 1:length(nattr1))
      #  attr(yy,nattr1[i]) <- attr(y,nattr1[i])
      yy <- attrcp(y,yy,ignore='names')
      attr(yy,'dimensions') <- c(attr(y,'dimensions')[1:2],nt)
      if (verbose) print(names(attributes(x)))
      #nattr2 <- softattr(x)
      #for (i in 1:length(nattr2))
      #  attr(XX,nattr2[i]) <- attr(x,nattr2[i])
      XX <- attrcp(x,XX,ignore='names')
      attr(XX,'dimensions') <- c(attr(x,'dimensions')[1:2],nt)
      
      #print("into list")
      X <- list(y=yy,X=XX)
      if (verbose) {print(names(attributes(X$y))); print(names(attributes(X$X)))}
    } else {
      X <- comb
    }
    class(X) <- c("comb","list")
  }
  #print("HERE")
  #print(names(attributes(combined$y)))
  #print(names(attributes(combined$X)))
  #print(summary(combined))
  #print("HERE")
  attr(X,'history') <- history.stamp(x)
  invisible(X)
}

#' @exportS3Method
#' @export combine.events
combine.events <- function(x=NULL,y=NULL,...,remove.close=TRUE,mindistance=5E5,FUN=NULL,verbose=FALSE) {
  if(verbose) print("combine.events")
  stopifnot(inherits(x,"events") & inherits(y,"events"))
  
  if(is.null(attr(x,"greenwich"))) attr(x,"greenwich") <- !(any(x$lon<0) & !any(x$lon>180))
  x <- g2dl.events(x,greenwich=attr(x,"greenwich"))
  y <- g2dl.events(y,greenwich=attr(x,"greenwich"))
  
  ## Combine events in x and y
  cn <- colnames(x)[colnames(x) %in% colnames(y)]
  if(!"trajectory" %in% cn) {
    z <- rbind(x[colnames(x) %in% cn],y[colnames(y) %in% cn])
    dt <- as.numeric(z$date)*1E2 + as.numeric(z$time)
    z <- z[order(dt, decreasing=FALSE), ]
  } else {
    if(require("PCICt")) {
      if(!is.null(attr(x,"calendar"))) cal <- attr(x,"calendar") else cal <- "gregorian"
      dh <- difftime(min(as.PCICt(paste(y$date,y$time), format="%Y%m%d %H", cal=cal)),
                     max(as.PCICt(paste(x$date,x$time), format="%Y%m%d %H", cal=cal)),
		     units="hours")
    } else {
      dh <- difftime(min(as.Date(paste(y$date,y$time), format="%Y%m%d %H")),
                     max(as.Date(paste(x$date,x$time), format="%Y%m%d %H")),
		     units="hours")
    }
    if(dh>6) {
      dt <- max(x$trajectory)-min(y$trajectory)+1
      y$trajectory <- y$trajectory+dt
      z <- rbind(x[colnames(x) %in% cn],y[colnames(y) %in% cn])
      remove.close <- FALSE
    } else {
      # If there is 6 hours or less between the end of x and beginning of y
      # check if the tracks of x continue in y. If not, track y again.
      trajectories_x <- x$trajectory[x$date==max(x$date)]
      trajectories_y <- y$trajectory[y$date==min(y$date)]
      if(any(c(trajectories_x, max(trajectories_x)+1) %in% trajectories_y)) {
        z <- rbind(x[colnames(x) %in% cn],y[colnames(y) %in% cn])
      } else {
        y2 <- y[!colnames(y) %in% c("trajectory","dx","trackcount","timestep","distance","tracklength")]
        y2 <- track(y2, x0=x)
        z <- rbind(x[colnames(x) %in% cn],y2[colnames(y2) %in% cn])
        remove.close <- FALSE
      }
    }
  }
  
  if(!any(x$date %in% y$date)) remove.close <- FALSE
  
  ## Remove events located close to other stronger events
  if(remove.close) {
    t <- z$date*1E2 + z$time
    lon <- z$lon
    lat <- z$lat
    ## Define measure of strength.
    ## Decides the which one is thrown out if two events are close together. 
    if(is.null(FUN)) {
      if(attr(x,"variable")=="cyclones") {
        strength <- 1/z$pcent
      } else if (attr(x,"variable")=="anti-cyclones") {
        strength <- z$pcent
      } else {
        strength <- rep(1,nrow(z))
        print(paste("Warning: No measure of strength was provided.",
                    "If two closely located events are detected one will be randomly removed.",
                    "You can turn off the removing of neighbouring events (remove.close=FALSE)",
                    "or provide a function (FUN) that describes the rank/strength of the events."))
      }
    } else {
      strength <- FUN(z)
    }
    if(verbose) print("Remove duplicate cyclones")
    del.z <- rep(TRUE,nrow(z))
    for (d in unique(t)) {
      i <- which(t==d)
      if (length(i)>1) {
        distance <- apply(cbind(lon[i],lat[i]),1,
                          function(x) suppressWarnings(distAB(x[1],x[2],lon[i],lat[i])))
        diag(distance) <- NA; distance[lower.tri(distance)] <- NA
        del.i <- which(distance<mindistance,arr.ind=TRUE)
        if(any(del.i)) {
          col.del <- rep(1,length(del.i)/2)
          if (is.null(dim(del.i))) {
            s1 <- strength[i][,del.i[1]]
            s2 <- strength[i][,del.i[2]]
            col.del[s1<=s2] <- 2
            del.i <- unique(del.i[col.del])
          } else {
            s1 <- strength[i][del.i[,1]]
            s2 <- strength[i][del.i[,2]]
            col.del[s1<=s2] <- 2
            del.i <- unique(del.i[cbind(seq(1,dim(del.i)[1]),col.del)])
          }
          del.z[i[del.i]] <- FALSE
        }
        i <- i[!1:length(i) %in% del.i]
      }
    }
    z <- z[del.z,]
  }
  z <- attrcp(x,z)
  class(z) <- class(x)
  return(z)
}

#' @exportS3Method
#' @export combine.trajectory
combine.trajectory <- function(x=NULL,y=NULL,...,verbose=FALSE) {
  if(verbose) print("combine.trajectory")
  z <- rbind(x,y,all=TRUE)
  z <- z[!duplicated(z),]
  z <- attrcp(x,z)
  class(z) <- class(x)
  return(z)
}

#' @export
sp2np <- function(x,SP2NP=TRUE,verbose=FALSE) {
  if(verbose) print("sp2np")
  # Take care of different latitude orientations: N-> S & S -> N
  ysrt <- order(attr(x,'latitude'))
  d <- attr(x,'dimensions')
  if (SP2NP) {
    if (ysrt[1]==1) return(x)
    y <- t(coredata(x))
    dim(y) <- attr(x,'dimensions')
    y <- y[,ysrt,]
    dim(y) <- c(d[1]*d[2],d[3])
    y <- zoo(t(y),order.by=index(x))
    class(y) <- class(x)
    y <- attrcp(x,y,ignore='names')
    #nattr <- softattr(x,ignore=c('longitude','latitude','dimensions'))
    #for (i in 1:length(nattr))
    #  attr(y,nattr[i]) <- attr(x,nattr[i])
    attr(y,'longitude') <- attr(x,'longitude')
    attr(y,'latitude') <- attr(x,'latitude')[ysrt]
    attr(y,'dimensions') <- d
  } else {
    ysrt <- order(attr(x,'latitude'), decreasing = TRUE)
    if (ysrt[1]==1) return(x)
    y <- t(coredata(x))
    dim(y) <- d
    y <- y[,ysrt,]
    dim(y) <- c(d[1]*d[2],d[3])
    dim(y) <- c(d[1]*d[2],d[3])
    class(y) <- class(x)
    y <- attrcp(x,y,ignore='names')
    #nattr <- softattr(x,ignore=c('longitude','latitude','dimensions'))
    #for (i in 1:length(nattr))
    #  attr(y,nattr[i]) <- attr(x,nattr[i])
    attr(y,'longitude') <- attr(x,'longitude')
    attr(y,'latitude') <- attr(x,'latitude')[ysrt]
    attr(y,'dimensions') <- attr(x,'dimensions')
  }
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

