#' Extract attributes
#'
#' \code{as.appended} extracts an appendix (e.g, attr(x,'appendix.1')) and copies the other attributes from the main object to the appendix.
#' 
#' \code{as.calibrationdata} extracts the calibration data from a DS object.
#' Note: if DS was performed with rm.trend=TRUE, the calibration data are detrended.
#'
#' \code{as.fitted.values} extracts fitted values from a DS object.
#'
#' \code{as.original.data} extracts the original data (the predictand) from a DS object.
#'
#' \code{as.pattern} extracts a spatial pattern, e.g., from an 'eof' or 'cca' object.
#'
#' @aliases as.appended as.appended.ds.comb as.appended.eof.comb as.appended.field.comb
#' as.calibrationdata as.calibrationdata.ds as.calibrationdata.station
#' as.fitted.values as.fitted.values.ds as.fitted.values.station
#' as.original.data as.original.data.ds as.original.data.station
#' as.pattern as.pattern.ds as.pattern.corfield as.pattern.eof as.pattern.cca as.pattern.mvr as.pattern.field
#'
#' @param x input object with appendix (of class \code{comb})
#' @param iapp index of appendix to extract
#' @param \dots additional arguments
#'
#' @export
as.appended <- function(x,...) UseMethod("as.appended")

#' @export
as.appended.ds.comb <- function(x,...,iapp=1,verbose=FALSE) {
  if(verbose) print("as.appended.ds.comb")
  eval(parse(text=paste("X <- attr(x,'appendix.",iapp,"')",sep="")))
  X <- attrcp(x,X,ignore='appendix')
  attr(X,'history') <- history.stamp(x)
  invisible(X)
}

#' @export
as.appended.eof.comb <- function(x,...,iapp=1) {
  X <- as.appended.ds.comb(x,iapp=iapp)
  invisible(X)
}

#' @export
as.appended.field.comb <- function(x,...,iapp=1) {
  X <- as.appended.ds.comb(x,iapp=iapp)
  invisible(X)
}

#' @export
as.calibrationdata <- function(x) UseMethod("as.calibrationdata")

#' @export
as.calibrationdata.ds <- function(x) {
  y <- attr(x,'calibration_data')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'calibration_data'))
  invisible(y)
}

#' @export
as.calibrationdata.station <- function(x) {
  if (!is.null(attr(x,'calibration_data'))) {
    y <- as.calibrationdata.ds(x)
  } else {
    y <- NULL
  }
  invisible(y)
}

#' @export
as.fitted.values <- function(x) UseMethod("as.fitted.values")

#' @export
as.fitted.values.ds <- function(x) {
  y <- attr(x,'fitted_values')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'fitted_values'))
  invisible(y)
}

#' @export
as.fitted.values.station <- function(x) {
   if (!is.null(attr(x,'fitted.values')))
    y <- as.fitted.values.ds(x) else y <- NULL
  invisible(y)
}

#' @export as.original.data
as.original.data <- function(x) UseMethod("as.original.data")

#' @export
as.original.data.ds <- function(x) {
  y <- attr(x,'original_data')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'original_data'))
  invisible(y)
}

#' @export
as.original.data.station <- function(x) {
  y <- as.original.data.ds(x)
  invisible(y)
}

#' @export
as.pattern <- function(x,...) UseMethod("as.pattern")

#' @export
as.pattern.ds <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.ds")
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

#' @export
as.pattern.eof <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.eof")
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

#' @export
as.pattern.mvr <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.mvr")
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

#' @export
as.pattern.cca <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.cca")
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

#' @export
as.pattern.trend <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.trend")
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

as.pattern.matrix <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.matrix")
  if(verbose) print("Unfinished function - returning input object")
  return(x)
}

as.pattern.array <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.array")
  if(verbose) print("Unfinished function - returning input object")
  return(x)
}

#' @export
as.pattern.field <- function(x,...,FUN=NULL,verbose=FALSE) {
  if(verbose) print("as.pattern.field")
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

#' @export
as.pattern.corfield <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.corfield")
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