## Script for filling in missing values using PCA. This makes sense
## for situations where the data to a large extent is accounted for
## by a few leading PCs/modes.

## Function for one time series based on multiple regression
## Allow EOFs with different number of PCs.
#' @export
fitpc <- function(y,x,ip=1:4) {
  caldat <- data.frame(y=y,x)
  #print(summary(caldat))
  ip <- ip[ip <=dim(x)[2]]
  fitstr <- paste('y ~ ',paste('X.',ip,sep='',collapse=' + '))
  fit <- eval(parse(text=paste('lm(',fitstr,', data=caldat)')))
  z <- predict(fit,newdata=caldat)
  invisible(z)
}

eoffit <- function(X,U,ip) {
  ## Use regression to project the pattern of observations onto the PCA pattern and estimate
  ## the PCs.
  caldat <- data.frame(X=X,U=U)
  names(caldat) <- c('X',paste('U.',ip,sep=''))
  #print(names(caldat))
  calexpr <- paste('X ~ ',paste('U.',ip,sep='',collapse=' + '))
  eval(parse(text=paste('projection <- lm(',calexpr,',data=caldat)')))
  V <- projection$coefficients
  invisible(V[-1])
}

#' PCA-based missing-value filling
#' 
#' Fills missing (station) values by predicting their values using multiple
#' regression. The regression uses as input pcincipal components from PCA from
#' the same (group of station) data, but where series with missing data have
#' been excluded. This makes sense for (station) data where most of the
#' variability is accounted for by a few leading modes. This method is not
#' expected to be useful when there are many large data gaps.
#' 
#' This function is handy for the downscaling of PCAs. See Benestad, R.E., D.
#' Chen, A. Mezghani, L. Fan, K. Parding, On using principal components to
#' represent stations in empirical-statistical downscaling, Tellus A 28326,
#' accepted.
#' 
#' @aliases pcafill pcafill.test fitpc
#' @param X station data (group of stations)
#' @param insertmiss Used for testing and evaluating. Missing data are
#' introduced to test the predictive capability
#' @param ip Number of EOFs/PCAs to include in filling in. In many cases, it
#' may be useful to keep this to a small set of values.
#' @param mnv Minimum number of valid data points for any given time. Can be
#' used to get around the problem with too many missing data
#' @param complete Use pattern projection between PCA pattern and original data
#' to get a complete record - otherwise a subset of times with sufficient data.
#' @param test Extra test - debugging
#' @param verbose Print diagnostics - debugging
#' @param N Number of runs in Monte-Carlo simulation
#' @param max.miss Maximum NAs to insert (insertmiss) in Monte-Carlo
#' simulations
#' @param x time series for calibrating regression analysis
#' @param y PC input for regression analysis
#' 
#' @return The same as the input - station object with filled-in values
#' 
#' @seealso \code{\link{PCA}}, \code{\link{allgood}}
#' 
#' @keywords PCA missing data
#' 
#' @examples
#' download.file('http://files.figshare.com/2073466/Norway.Tx.rda',
#'               'Norway.Tx.rda')
#' load('Norway.Tx.rda')
#' X <- annual(Tx,FUN='mean',nmin=200)
#' ok<- apply(X,1,nv)
#' X <- subset(X,it=ok > 0)
#' Y <- pcafill(X)
#' 
#' plot(PCA(Y))
#' plot(c(coredata(Y)),c(coredata(X)))
#' 
#' ## Monte-Carlo test with random selection of data points set to NA:
#' Y.test <- pcafill.test(X,max.miss=10,ip=1:3)
#' cor(Y.test)
#' 
#' 
#' @export pcafill
pcafill <- function(X,insertmiss=0,ip=1:4,mnv=0,complete=FALSE,test=FALSE,verbose=FALSE) {
  if (verbose) print('pcafill')
  X0 <- X ## For debugging
  if (verbose) print(dim(X0))
  if (insertmiss>0) {
    ## Test by inserting false missing values in the data
    if (verbose) print(paste('Test: insert',insertmiss,'NAs'))
    x <- c(coredata(X))
    ok <- is.finite(x)
    ## Only replace valid data
    xok <- x[ok]
    ## Pick random samples
    isx <- order(rnorm(sum(ok)))[1:insertmiss]
    if (verbose) print(c(isx,sum(ok)))
    ## Save the original data set to NA
    x0 <- xok[isx]
    ## Set theselected data points ot NA
    xok[isx] <- NA
    ## Replace the data with that with introduced NAs
    x[ok] <- xok
    dim(x) <- dim(X)
    coredata(X) <- x
    if (verbose) print('---')
  }
  nok <- apply(X,1,nv)
  ## Remove the years with no (of little) data
  X <- subset(X,it=nok > mnv)
  if (verbose) print(dim(X))
  mok <- apply(X,2,nv)

  if (verbose) print(paste(sum(nok>mnv),'stations with',
                           sum(mok>0),'data points'))
  if (sum(mok==length(index(X))) <= 1)
    stop('pcafill: Too many missing data or too small set of stations')
  
  ## PCA for stations with complete data
  pca <- PCA(subset(X,is=mok==length(index(X))),n=max(ip))

  ## Y contains the more complete data set - copy all the attributes from X
  Y <- X
  ## Replace the data in Y with predicted based on the PCA
  #  coredata(Y) <- fillmiss(X,pca,neofs=neofs)
  coredata(Y) <- apply(coredata(X),2,fitpc,coredata(pca),ip=ip)

  if (complete) {
    print("UNFINISHED")
    if (verbose) print('Projection to provide complete dataset')
    if (verbose) print(dim(Y))
    ## Projection of patterns onto the original data to get a more complete data set
    ## X0 = U W t(V) -> estimate new t(V) based on U & W from Y
    ## 1/W t(U) X0 = t(V) -> V = t(X0) U /W
    pcax <- PCA(Y); pcax <- subset(pcax,ip=ip)
    U <- attr(pcax,'pattern'); W <- attr(pcax,'eigenvalues')
    U <-  U %*% diag(W)
    XX <- X0
    XX <- t(t(XX) - rowMeans(XX,na.rm=TRUE))
    ## Estimate the PCs for the entire data set through projection:
    Vx <- zoo(t(apply(t(coredata(XX)),2,eoffit,U=U,ip=ip)),order.by=index(X0))
    ## assign the attributes from the PCA
    Vx <- attrcp(pcax,Vx); class(Vx) <- class(pcax) 
    plot.zoo(pca[,1])
    lines(Vx[,1],col='red')
    ## Reconstruct the data:
    Y <- pca2station(Vx,verbose=verbose)
  }
  
  if (insertmiss>0) {
    if (verbose) print('Assess the test results')
    ## For testing, extract the data points replaced by NAs
    if (test) y <- c(coredata(X0)) else ## For debugging
              y <- c(coredata(Y))
    ## Use the same points with valid data
    yok <- y[ok]
    if (verbose) print(sum(!is.finite(xok[isx])))
    y0 <- yok[isx]
    if (sum(insertmiss)>1) {
      xy0 <- cbind(x0,y0); colnames(xy0) <- c('original','predicted')
    } else {
      xy0 <- c(x0,y0); names(xy0) <- c('original','predicted')
    }
    
    if ( (test) & (sum(is.element(x0,y0)) != length(x0)) )
      warning('The test results should be identical!')
    if (verbose) print(xy0)
    if (!is.null(attr(Y,'na.test')))
        attr(Y,'na.test') <- rbind(attr(Y,'na.test'),xy0) else
        attr(Y,'na.test') <- xy0
  }
  invisible(Y)
}

pcafill.test <- function(X,N=100,max.miss=100,ip=1:4,verbose=FALSE) {
  if(verbose) print('pcafill.test')
  insertmiss <- round(runif(N)*max.miss)
  insertmiss[insertmiss<1] <- 2
  par(bty='n')
  plot(range(c(coredata(X)),na.rm=TRUE),
       range(c(coredata(X)),na.rm=TRUE),
       type='l',lwd=3,col='grey30',
       xlab='original',ylab='predicted',
       main='Test of prediction of missing values',
       sub="NA's randomly inserted")
  grid()
  for (i in 1:N) {
    if (verbose) print(insertmiss[i])
    Y.test <- pcafill(X,insertmiss=insertmiss[i],ip=ip,verbose=verbose)
    points(attr(Y.test,'na.test'),pch=19,col=rgb(0,0,0.6,0.1))
    if (i==1) test.res <- attr(Y.test,'na.test') else
              test.res <- rbind(test.res,attr(Y.test,'na.test'))
  }
  invisible(test.res)
}
