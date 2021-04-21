#' Multivariate combination of fileds
#'
#' \code{mvcombine} combines two fields with different variables by attaching them along the spatial axis.
#' Before combining the two fields, the same time slices are selected from both fields and
#' the data are normalized by subtracting the mean and dividing by the standard deviation.
#' The output of mvcombine is an 'mvcomb' object to which you can apply a multivariate
#' Empirical Orthogonal Function (EOF) using the function \code{EOF}.
#'
#' @seealso EOF.mvcomb EOF 
#'
#' @param X A field object.
#' @param Y A field object.
#' @param is A list or data.frame providing space index, e.g. a list of longitude and latitude range like list(lon=c(0,60), lat=c(35,60)).
#' @param verbose If TRUE, print out diagnosics.
#'
#' @export
mvcombine <- function(X, Y, is=NULL, verbose=FALSE) {
  if(verbose) print("mvcombine")
  if(verbose) print("Subset same time range of the two fields")
  it.12 <- index(X)[index(X) %in% index(Y)]
  X <- subset(X, it=it.12, is=is)
  Y <- subset(Y, it=it.12, is=is)
  
  if(!is.null(attr(X,'greenwich')) & !is.null(attr(Y,'greenwich'))) { 
    if(attr(X, "greenwich")!=attr(Y, "greenwich")) {
      Y <- g2dl(Y, attr(X,'greenwich'))
    }
  }
  
  if(verbose) print("Normalize the data")
  X <- sp2np(X)
  Y <- sp2np(Y)
  mean.x <- mean(apply(X, 2, mean, na.rm=TRUE), na.rm=TRUE)
  mean.y <- mean(apply(Y, 2, mean, na.rm=TRUE), na.rm=TRUE)
  sd.x <- sd(as.vector(X), na.rm=TRUE)
  sd.y <- sd(as.vector(Y), na.rm=TRUE)
  X.norm <- (X-mean.x)/sd.x
  Y.norm <- (Y-mean.y)/sd.y
  if(verbose) print("Join the two fields")
  XY <- cbind(X.norm, Y.norm)
  if(verbose) print("Update attributes of joined field")
  attr(XY, "id") <- c(rep(1, ncol(X)), rep(2, ncol(Y)))
  attr(XY, "mean") <- list(mean.x, mean.y)
  attr(XY, "sd") <- list(sd.x, sd.y)
  attr(XY, "dimnames") <- list(index(XY), paste0("x.",seq(ncol(XY))))
  attrnames <- names(attributes(X))
  attrnames <- attrnames[!attrnames %in% c("dim","dimnames","class","index")]
  for(nm in attrnames) {
    eval(parse(text=paste0('attr(XY, nm) <- list(attr(X, nm), attr(Y, nm))')))
  }
  class(XY) <- c("mvcomb", class(XY))
  return(XY)
}
