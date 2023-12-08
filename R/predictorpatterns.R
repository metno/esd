#' Extract predictor patterns from a dse-object and apply a common EOF evaluation
#'
#' @param x A DSensemble object.
#' @param plot Plot intermediate results if TRUE.
#' @param verbose TRUE for checking and debugging the functions.
#'
#' @return A common EOF object of predictor patterns for different models
#' 
#' @keywords manip
#' @examples
#' 
#' \dontrun{

#' }
#' 
#' @export 
predictorpatterns <- function(x,ip=1,plot=FALSE, verbose=FALSE,...) {
  if (verbose) print('predictorpatterns')
  x[["info"]] <- NULL
  x[["pca"]] <- NULL
  x[["eof"]] <- NULL
  cn <- names(x)
  n <- length(cn)
  if (verbose) print(cn)
  Y <- subset.pattern(attr(x[[1]],'predictor.pattern'),ip=ip,verbose=verbose)
  nxy <- length(Y)
  X <- matrix(rep(NA,n*nxy),n,nxy)
  X[1,] <- c(Y)
  for (i in 1:n) {
    x1 <- x[[cn[i]]]
    z <- subset.pattern(attr(x1,'predictor.pattern'),ip=ip,verbose=verbose)
    z <- regrid(z,is=Y,verbose=TRUE)
    X[i,] <- c(z)
    if (verbose) print(paste(i,': ',cn[i],'lenth=',length(z),'mean=',round(mean(c(z),na.rm=TRUE),2)))
  }
  
  X <- as.field(zoo(X,order.by=1:n),param='weights',unit="none",lon=lon(Y),lat=lat(Y),
                info='record maps: predictorpattern',longname='Multi-model ESD predictor patterns')
  
  ceof <- EOF(X)
  if (plot) {
    map(X)
    map(X,FUN='sd',main='Multi-modal spread (CMIP6)')
    plot(ceof,ip=ip)
  }
  par(new=FALSE)
  attr(ceof,'ID') <- cn 
  srt <- order(ceof[,ip])
  if (verbose) print(cbind(cn[srt],ceof[srt,ip]))
  invisible(ceof)
}