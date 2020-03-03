#' Validate
#' 
#' The method \code{validate}
#' 
#' @aliases validate.dsensemble validate.eof.comb validate.default validate.cca
#' validate.ds validate.eof validate.eof.field validate.pca
#'
#' @importFrom stats wilcox.test
#'
#' @param x esd object to be validated
#' @param conf.int confidence interval
#' @param colbar for plotting. See \code{colbar}
#' @param plot if TRUE produce plot
#' @param verbose if TRUE print progress
#' @param \dots other arguments
#' @examples
#' 
#'   slp1 <- slp.DNMI(lon=c(-50,50),lat=c(30,70))
#'   slp2 <- slp.DNMI(lon=c(-50,50),lat=c(30,70))
#'   slpcomb <- combine(slp1,slp2)
#'   eofcomb <- EOF(slpcomb)
#'   validate(eofcomb, new=FALSE)
#' 
#' @export
validate <- function(x, ...)  UseMethod("validate")

#' @export validate.eof.comb
validate.eof.comb <- function(x,...,new=TRUE,verbose=FALSE) {
  if(verbose) print("validate.eof.comb")
  zz <- attr(x,'appendix.1')
  if (new) dev.new()
  plot(attr(x,'mean'),type="l",main="Mean values")
  lines(attr(zz,'mean'),col="red")
}

#' @export validate.dsensemble
validate.dsensemble <- function(x, conf.int=c(0.05,0.95),text=FALSE,
                                colbar=list(breaks=seq(0,1,by=0.1),cex=1.5,
                                col=colscal(11,pal="t2m",alpha=0.5)),
				plot=TRUE,verbose=FALSE,...) {
  if (verbose) print('validate.dsensemble')
  ranktest <- function(x) {
    ## Convert to zoo objects and extract data for common times
    obs <- zoo(attr(x,'station'))
    res <- zoo(x)
    ## Make sure both station and results have same date stamp type
    if (is.numeric(index(res))) index(obs) <- year(obs)
    if (is.numeric(index(obs))) index(res) <- year(res)
    comb <- coredata(merge(obs,res,all=FALSE))
    ## Find the rank order of the observation amongst the model results for each time step
    rs <- apply(comb,1,function(x) order(x)[1])
    medH0 <- median(1:dim(comb)[2])
    p.val <- wilcox.test(rs, alternative="two.sided", mu=medH0, conf.int=TRUE)$p.value
    return(p.val)
  }
  
  if (verbose) print('as.station')
  if (inherits(x,'pca')) x <- as.station(x)
  if (verbose) print('ranktest')
  if(is.list(x)) {
    ro <- unlist(lapply(x,ranktest))
    attr(ro,'longitude') <- unlist(lapply(x,function(x) lon(attr(x,'station'))))
    attr(ro,'latitude') <- unlist(lapply(x,function(x) lat(attr(x,'station'))))
  } else {
    ro <- ranktest(x)
    attr(ro,'longitude') <- lon(attr(x,'station'))
    attr(ro,'latitude') <- lat(attr(x,'station'))
  }
  if (verbose) print('attributes')
  attr(ro,'variable') <- 'Wilcox-test score'
  attr(ro,'unit') <- 'probability'
  attr(ro,'longname') <- 'Wilcox-test of how observation ranks amongst model results'
  attr(ro,'history') <- history.stamp(x)
  if(plot) {
  if (verbose) print('visualisation')
  if (is.null(colbar)) colbar <- colbar.ini(ro,verbose=verbose)
  ic <- round(as.numeric(ro)*length(colbar$col))
  ic[ic < 1] <- 1; ic[ic > length(colbar$col)] <- length(colbar$col)
  cols <- colbar$col[ic]
  ## Plot the results
    par0 <- par()
    par(bty='n',fig=c(0,1,0,0.82))
    plot(lon(ro),lat(ro),pch=19,col=cols,cex=colbar$cex,
         xlab='',ylab='',main='Validation of downscaled ensemble',
         sub='p-values from a Wilcox-test: how observation ranks amongst model results')
    data(geoborders, envir = environment())
    lines(geoborders,col='grey')
    points(lon(ro),lat(ro),cex=colbar$cex,col='grey')
    good <- (ro > conf.int[1]) & (ro < conf.int[2])
    points(lon(ro)[good],lat(ro)[good],cex=colbar$cex,lwd=2)
    if (text) text(lon(ro),lat(ro),as.numeric(round(ro*100)),col='grey',cex=0.75,pos=1)
    colbar(colbar$breaks,colbar$col,fig=c(0.12,0.15,0.75,0.90))
    par(bty='n',fig=c(0.7,1,0.8,0.95),new=TRUE,xaxt='n')
    hist(ro,breaks=colbar$breaks,col=colbar$col,main='')
    par(bty='n',fig=c(0.25,0.7,0.8,0.95),new=TRUE,xaxt='n',yaxt='n')
    plot(c(0,1),c(0,1),type='n',xlab='n',ylab='n')
    legend(0,1,c(paste('ci=[',paste(conf.int,collapse=','),']'),
                 expression(p %in% ci),
                 expression(p %notin% ci)),pch=21,lty=0,
           lwd=c(1,2,1),col=c('white','black','grey'),bty='n')
    par0$cin<-NULL; par0$cra<-NULL; par0$csi<-NULL; par0$cxy<-NULL; par0$din<-NULL; par0$page<-NULL
    par(par0)
  }
  ## return the results
  invisible(ro)
}

#' @export validate.default
validate.default <- function(x,...,verbose=FALSE) {
  if(verbose) print("validate.default")
  if(verbose) print("unfinished function - not returning anything")
}

#' @export validate.eof
validate.eof <- function(x,...,verbose=FALSE) {
  if(verbose) print("validate.eof")
  if(inherits(x,"comb")) {
    validate.eof.comb(x,...,verbose=verbose)
  } else {
    if(verbose) print("unfinished function - not returning anything")
  }
}

#' @export validate.pca
validate.pca <- function(x,...,verbose=FALSE) {
  if(verbose) print("validate.eof.field")
  if(verbose) print("unfinished function - not returning anything")
}

#' @export validate.eof.field
validate.eof.field <- function(x,...,verbose=FALSE) {
  if(verbose) print("validate.eof.field")
  if(verbose) print("unfinished function - not returning anything")
}

#' @export validate.ds
validate.ds <- function(x,...,verbose=FALSE) {
  if(verbose) print("validate.ds")
  if(verbose) print("unfinished function - not returning anything")
}

#' @export validate.cca
validate.cca <- function(x,...,verbose=FALSE) {
  if(verbose) print("validate.cca")
  if(verbose) print("unfinished function - not returning anything") 
}
