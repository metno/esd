#' @export
trackfilter <- function(x,...) UseMethod("trackfilter")

#' @export
trackfilter.events <- function(x,param=NULL,pmin=NULL,pmax=NULL,FUN="any",verbose=FALSE) {
  if(verbose) print("trackfilter")
  stopifnot(inherits(x,"events"))
  if(!is.null(param) & (!is.null(pmin) | !is.null(pmax))) {
    #y <- subset(x,ic=list(param=param,pmin=pmin,pmax=pmax,FUN=FUN),verbose=verbose)
    if(!param %in% names(x)) {
      if(verbose) print(paste("Unkown input param =",param))
      y <- x
    } else {
      if(param %in% c("trajectory","trackcount","distance","tracklength")) FUN <- NULL
      if(is.null(pmin)) pmin <- min(x[param],na.rm=TRUE)
      if(is.null(pmax)) pmax <- max(x[param],na.rm=TRUE)
      if(verbose) print(paste(param,"in range",pmin,"-",pmax))
      if(verbose) print(paste("FUN =",FUN))
      if(!"trackcount" %in% names(x)) x <- trackstats(x)
      if(is.null(FUN)) {
        ok <- as.vector(x[param]>=pmin & x[param]<=pmax)
      } else if (FUN=="any") {
        if(!"trajectory" %in% names(x)) x <- track(x)
        ok.ev <- as.vector(x[param]>=pmin & x[param]<=pmax)
        ok <- x$trajectory %in% unique(x$trajectory[ok.ev])
      } else if (FUN=="all") {
        nok.ev <- as.vector(x[param]<pmin | x[param]>pmax)
        ok <- !x$trajectory %in% unique(x$trajectory[nok.ev])
      }
      y <- subset(x,it=ok)
    }
  } else {
    if(verbose) "No filter applied. Input 'param' or range ('pmin', 'pmax') not provided."
    y <- x
  }
  invisible(y)
}

#' @export
trackfilter.trajectory <- function(x,param=NULL,pmin=NULL,pmax=NULL,FUN="any",verbose=FALSE) {
  if(verbose) print("trackfilter")
  stopifnot(inherits(x,"trajectory"))
  if(!is.null(param) & (!is.null(pmin) | !is.null(pmax))) {
    #y <- subset(x,ic=list(param=param,pmin=pmin,pmax=pmax,FUN=FUN),verbose=verbose)
    if(!param %in% colnames(x)) {
      if(verbose) print(paste("Unkown input param =",param))
      y <- x
    } else {
      ip <- colnames(x)==param  
      if(is.null(pmin)) pmin <- min(x[,ip],na.rm=TRUE)
      if(is.null(pmax)) pmax <- max(x[,ip],na.rm=TRUE)
      if(verbose) print(paste(param,"in range",pmin,"-",pmax))
      if(verbose) print(paste("FUN =",FUN))
      if (FUN %in% c("any","all")) {
        fn <- function(x) do.call(FUN,list(x[ip]>=pmin & x[ip]<=pmax))
      } else {
        fn <- function(x) do.call(FUN,list(x[ip]))>=pmin &  do.call(FUN,list(x[ip]))<=pmax
      }
      ok <- apply(x,1,fn)
      y <- subset(x,it=ok)
    }
  } else {
    if(verbose) "No filter applied. Input 'param' or range ('pmin', 'pmax') not provided."
    y <- x
  }
  invisible(y)
}
