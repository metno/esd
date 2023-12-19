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
    z <- regrid(z,is=Y)
    X[i,] <- c(z)
    if (verbose) print(paste(i,': ',cn[i],'lenth=',length(z),'mean=',round(mean(c(z),na.rm=TRUE),2)))
  }
  
  X <- as.field(zoo(X,order.by=1:n),param='weights',unit="none",lon=lon(Y),lat=lat(Y),
                info='record maps: predictorpattern',longname='Multi-model ESD predictor patterns')
  
  ceof <- EOF(X)
  if (plot) {
    map(X)
    map(X,FUN='sd',main='Multi-modal spread (CMIP6)')
    mid <-  xgcmname(cn)
    tmid <- sort(table(mid),decreasing=TRUE)
    print(tmid)
    modelid <- rownames(tmid)
    col <- rep('grey',length(index(mid)))
    cols <- c('blue',"brown","purple","darkgreen","orange","cyan","black","red","darkblue")
    for (i in 1:min(c(length(cols),length(modelid)))) col[grep(modelid[i],mid)] <- cols[i]
    plot(coredata(ceof)[,ip],main=paste('PC#',ip),col=col,xlab='Simulation',ylab='weight',pch=19)
    legend(0,max(coredata(ceof)[,ip]),modelid[1:9],col=cols,pch=19,bty='n',cex=0.7)
  }
  par(new=FALSE)
  attr(ceof,'ID') <- cn 
  srt <- order(ceof[,ip])
  if (verbose) print(cbind(cn[srt],ceof[srt,ip]))
  invisible(ceof)
}

xgcmname <- function(x) {
  if (length(x) > 1) {
    for (i in 1:length(x)) x[i] <- xgcmname(x[i])
    return(x)
  }
  x <- gsub('.','_',x,fix=TRUE)
  x <- strsplit(x,split='_')
  #print(x)
  ns <- length(x[[1]])
  x <- x[[1]]
  x <- x[-1]
  x <- x[-length(x)]
  return(paste(x,collapse='_'))
}