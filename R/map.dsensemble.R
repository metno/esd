#' expandpca
#' 
#' \code{map.dsensemble} is a wrapper for \code{map} that distills ensemble statistics from PCA-based dsensemble objects. 
#' The argument 'FUNX' is the function that is applied to the ensemble. If 'FUNX' is 'mean', the function \code{expandpca} is 
#' used to aggregate the data, otherwise the function \code{aggregate.dsensemble} is employed. The difference is that 
#' \code{expandpca} applies 'FUNX' to the PCs before expanding the data, while \code{aggregate.dsensemble} expands the data
#' for all ensemble members before aggregating them which is more suitable when calculating e.g. the .
#' 
#' @param x an object of type 'pca'
#' @param it time index (see \code{\link{subset}})
#' @param FUN function applied to aggregate in time
#' @param FUNX function applied aggregate ensemble members
#' @param verbose if TRUE print progress
#' @param anomaly if FALSE add the mean value stored as attribute in x
#' @param test if TRUE perform test on one GCM simulation
#' 
#' @seealso map aggregate.dsensemble expandpca
#' @aliases map.dsensemble

#' @exportS3Method
#' @export
map.dsensemble <- function(x,it=c(2000,2099),is=NULL,im=NULL,ip=NULL,
                           colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                                       show=TRUE,type="p",cex=2,h=0.6,v=1),
                           FUN='mean',FUNX='mean',verbose=FALSE,anomaly=FALSE,test=FALSE,plot=TRUE,...) {
  ## PCA/EOF objects
  
  if (verbose) print('map.dsensemble')
  
  if (inherits(x,c('pca','eof'))) {
    ## Extract a subset of the data
    if (verbose) print(names(x)[2])
    x <- subset(x,is=is,im=im,ip=ip,verbose=verbose)
    ## REB 2016-12-01: Do all the analysis on the PC weights to speed up. Linearity.  
    #    Y <- expandpca(x,it=it,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
    if (FUNX=='mean') 
      Y <- expandpca(x,it=it,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test) else
        Y <- aggregate.dsensemble(x,it=it,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
      ## KMP 2022-04-06: aggregate.dsensemble doesn't have an argument FUN and in expandpca it gives unexpected results   
      #  Y <- expandpca(x,it=it,FUN=FUN,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test) else
      #  Y <- aggregate.dsensemble(x,it=it,FUN=FUN,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
    if (verbose) {str(x[[2]]); str(Y)}
    #    if (plot) map(Y,FUN=FUN,colbar=colbar,verbose=verbose,...)
    if (verbose) {print('Visualise...'); print(dim(Y))}
    if (plot | !is.null(FUN)) {
      if (min(dim(Y))==0) stop(paste('dim(Y)=',paste(dim(Y),collapse='x')))
      z <- map(Y,FUN=FUN,plot=plot,colbar=colbar,verbose=verbose,...)
      attr(z,'ensemble.data') <- Y
    } else z <- Y
    invisible(z)
  } else return(NULL)
}






