#' @exportS3Method
#' @export
map.mvcomb <- function(x, ..., it=NULL, is=NULL, new=FALSE, projection="lonlat",
                       what="eof", xlim=NULL, ylim=NULL, zlim=NULL, FUN="mean", na.rm=TRUE,
                       colbar=list(pal=NULL, rev=FALSE, n=10, breaks=NULL,
                                   pos=0.05, show=TRUE, type="p", cex=2, h=0.6, v=1),
                       type=c("fill", "contour"), gridlines=FALSE,
                       lonR=NULL, latR=NULL, axiR=NULL, verbose=FALSE,
                       iv=1, ip=1, cex=1, plot=TRUE) {
  if(verbose) print("map.mvcomb")
  stopifnot(inherits(x,'mvcomb'))
  projection <- tolower(projection)
  
  if(inherits(x,'ds')) {
    print("Map mvcomb ds")
    patterns <- attr(x, "pattern")
    X <- patterns[[iv]]
    attr(X,"longitude") <- attr(patterns,"longitude")[[iv]]
    attr(X,"latitude") <- attr(patterns,"latitude")[[iv]]
    attr(X,"variable") <- attr(patterns,"variable")[[iv]]
    attr(X,"unit") <- attr(patterns,"unit")[[iv]]
    attr(X,"greenwich") <- attr(patterns,"greenwich")[[iv]]
    attr(X,"time") <- range(index(x))
  } else if(inherits(x,'eof')) {
    print("Map mvcomb eof")
    X <- attr(x, "pattern")[[iv]][, , ip]
    attr(X, "longitude") <- attr(x, "longitude")[[iv]]
    attr(X, "latitude") <- attr(x, "latitude")[[iv]]
    attr(X, "variable") <- attr(x, "variable")[[iv]]
    attr(X, "unit") <- attr(x, "unit")[[iv]]
    attr(X, "greenwich") <- attr(x, "greenwich")[[iv]]
    attr(X,"time") <- range(index(x))
  } else {
    print("Map mvcomb")
    X <- x[,attr(x,"id")==unique(attr(x,"id"))[iv]]
    X <- X*attr(x,"sd")[[iv]] + attr(x,"mean")[[iv]]
    attrnames <- names(attributes(x))
    attrnames <- attrnames[!attrnames %in% c("dimnames","mean","sd")]
    for(nm in attrnames) {
      if(is.list(attr(x, nm)) & length(attr(x, nm))==length(unique(attr(x,"id")))) {
        attr(X, nm) <- attr(x, nm)[[iv]]
      }
    }
    attr(X, "dimnames") <- list(index(x), paste0("x.",seq(ncol(X))))
    class(X) <- class(x)[!class(x)=="mvcomb"]
    map(X, FUN=FUN, it=it, is=is, new=new,
        projection=projection, xlim=xlim, ylim=ylim, zlim=zlim, 
        n=n, colbar=colbar, type=type, gridlines=gridlines,
        lonR=lonR, latR=latR, axiR=axiR, verbose=verbose,
        na.rm=na.rm, plot=plot, add=add, ...)
    return()
  }
  
  if (!is.null(zlim)) {
    d <- dim(X)
    mask <- (X < min(zlim)) | (X > max(zlim))
    X[mask] <- NA
    dim(X) <- d
    if (verbose) {print(zlim); print(dim(X)); print(sum(mask))}
  }
  
  if (plot) {
    if (projection=="lonlat") {
      lonlatprojection(x=X,it=it,xlim=xlim,ylim=ylim,
                       colbar=colbar,new=new,type=type,
                       gridlines=gridlines,verbose=verbose,...)
    } else if (projection=="sphere") {
      map2sphere(x=X,it=it,lonR=lonR,latR=latR,axiR=axiR,
                 xlim=xlim,ylim=ylim,type=type,gridlines=gridlines,
                 colbar=colbar,new=new,verbose=verbose,...)
    } else if (projection=="np") {
      map2sphere(X,it=it,lonR=lonR,latR=90,axiR=axiR,
                 xlim=xlim,ylim=ylim,type=type,gridlines=gridlines,
                 colbar=colbar,new=new,verbose=verbose,...)
    } else if (projection=="sp") {
      map2sphere(X,it=it,lonR=lonR,latR=-90,axiR=axiR,
                 xlim=xlim,ylim=ylim,type=type,gridlines=gridlines,
                 colbar=colbar,new=new,verbose=verbose,...)
    }
  }
  z <- X
  invisible(z)
}
