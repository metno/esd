#' Prediction based on DS or CCA model
#' 
#' Apply an empirical-statistical downscaling model to new data
#' 
#' \code{predict} is similar to the predict function in R
#' 
#' \code{project] returns projection of climate
#' 
#' @aliases predict.ds predict.ds.eof predict.ds.comb predict.mvr predict.cca project.ds
#'
#' @param x A ds object
#' @param newdata An eof object containing the new data sets on which the
#' prediction is made. 
#' @param addnoise If TRUE, will add an attribute called "noise" to the ouput
#' based on WG
#' @param n Number of runs to be generated, used only if addnoise is set to
#' TRUE 
#' 
#' @return Predicted ds values.
#' @seealso \code{\link{DS}} 
#' 
#' @examples
#' 
#' # Get predictor
#' ## Get reanalysis
#' X <- t2m.DNMI(lon=c(-40,50),lat=c(40,75))
#' ## Get Gcm output
#' Y <- t2m.NorESM.M(lon=c(-40,50),lat=c(40,75))
#' ## Combine
#' XY <- combine(X,Y)
#' # Compute common eof for January
#' ceof <- EOF(XY,it='jan')
#' # Get predictand
#' data(Oslo)
#' # Do the downscaling
#' ds <- DS(Oslo,ceof)
#' # Plot ds results
#' plot(ds)
#' # Do the prediction based on the calibration (or the fitted values)
#' ds.pre <- predict(ds)
#' # Plot predicted results based on ds object
#' plot(ds.pre)
#' # Display the attribute "aspect"
#' attr(ds.pre, "aspect")
#' ## Extract the projected results
#' plot(project.ds(ds))
#' 
#' @export
predict.ds <- function(x,newdata=NULL,...,addnoise=FALSE,n=100,verbose=FALSE) {
  if (verbose) print(paste("predict.ds",paste(class(x),collapse='-')))
  stopifnot(!missing(x),inherits(x,"ds"))

  if ( inherits(x,c('eof','comb')) & (is.null(newdata) | is.logical(newdata)) ) {
    if(verbose) print("no new predictor data is provided")
    if(verbose) print("predictand is an EOF")
    if (inherits(x,'comb')) {
      if(verbose) print("predictand is a combined EOF")
      y <- predict.ds.comb(x,newdata=newdata,addnoise=addnoise,
                           n=n,verbose=verbose)
    } else {
      if(verbose) print("predictand is an EOF, but use PCA call")
      y <- predict.ds.pca(x,newdata=newdata,addnoise=addnoise,
                          n=n,verbose=verbose)
    }
  } else if (inherits(x,'field')) {
    if (verbose) print("predictand is a field object")
    
    if (inherits(x,'station')) {
      y <- predict.ds.eof(x,newdata=newdata,addnoise=addnoise,
                          n=n,verbose=verbose) 
    } else {
      y <- predict.ds.pca(x,newdata=newdata,addnoise=addnoise,
                          n=n,verbose=verbose)
    }
  } else if (inherits(x,'eof')) {
    if (verbose) print("predictand is an EOF") 
    if(verbose) print("new predictor data is provided")
    if (inherits(newdata,'comb')) {
      if (verbose) print("new predictor data is a combined EOF") 
      y <- predict.ds.comb(x,newdata=newdata,addnoise=addnoise,
                           n=n,verbose=verbose)
    } else if (inherits(newdata,'eof')) {
      if (verbose) print("new predictor data is an EOF") 
      #      y <- predict.ds.eof(x,newdata=newdata,addnoise=addnoise,
      #                          n=n,verbose=verbose)
      y <- predict.ds.pca(x,newdata=newdata,addnoise=addnoise,
                          n=n,verbose=verbose)
    }
  } else if (inherits(x,'pca')) {
    if(verbose) print("predictand is PCA")
    y <- predict.ds.pca(x,newdata=newdata,addnoise=addnoise,n=n,verbose=verbose)
  } else if (inherits(x,'station')) {
    y <- predict.ds.station(x,newdata=newdata,addnoise=addnoise,n=n,verbose=verbose)
  }
  y <- attrcp(attr(x,'original_data'),y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

#' @export
predict.ds.station <- function(x,newdata=NULL,...,addnoise=FALSE,n=100,verbose=FALSE) {
  if (verbose) print(paste("predict.ds.pca",paste(class(x),collapse='-')))
  if (is.null(names(newdata))) names(newdata) <- paste('X',1:dim(newdata)[2],sep='.')
  t <- index(as.eof(newdata))
  newdata <- data.frame(coredata(as.eof(newdata)))
  ## Fudge fix for when the names of the data are x.1, x.2,,, instead of X.1, X.2,.. 
  if (nchar(names(newdata))[1]!=3) names(newdata) <- paste('X',1:dim(newdata)[2],sep='.')
  names(newdata) <- toupper(names(newdata))
  d <- dim(newdata)
  if (is.null(d)) d <- c(length(x),1)
  model <- attr(x,'model')
  y <- predict(model,newdata)
  y <- zoo(y, order.by=t)
  # KMP 2018-12-30: Changed to include model in attrcp. Is there any reason not to copy model to y? 
  y <- attrcp(x, y, ignore=c("name", "n.apps", "appendix", "dimnames","aspect"))
  class(y) <- class(x)[-1]
  invisible(y)
}

#' @export
predict.ds.eof <- function(x,newdata=NULL,...,addnoise=FALSE,n=100,verbose=FALSE) {
  stopifnot(!missing(x),inherits(x,"ds"))
  if (verbose) print(paste("predict.ds.eof",paste(class(x),collapse='-')))
  X <- as.eof(x)
  if (verbose) print(paste(class(X),collapse='-'))
  #print(dim(X))
  if (is.null(newdata)) neofs <- length(attr(X,'eigenvalues')) else {
    if (inherits(newdata,'eof')) 
      neofs <- length(attr(newdata,'eigenvalues')) else {
        if (!is.null(dim(newdata))) neofs <- dim(newdata)[2] else
          neofs <- 1
      }
  }
  
  # For some reason, the column names of newdata is muddled here,
  # and hence Xnames is used to enforce names 'X.1', 'X.2', 'X.3', ...
  Xnames <- paste("X.",1:neofs,sep="")
  if (verbose) print(Xnames)
  if (is.null(newdata)) {
    if (verbose) print('Use calibration data')
    newdata <- data.frame(X=coredata(X))
    src <- attr(X,'source')
    idx <- index(X)
  } else {
    if (verbose) print('Use new data')
    idx <- index(newdata)
    src <- attr(newdata,'source')
    newdata <- as.data.frame(newdata)
    if ((length(attr(X,'eigenvalues'))) != neofs)
      warning(paste('Warning: newdata and X have different number of EOFs:',length(attr(X,'eigenvalues')),neofs))
  }
  #print(summary(newdata))
  names(newdata) <- Xnames 
  model <- attr(x,'model')
  if (verbose) {print('names(newdata):'); print(Xnames)}
  ## AM 04-04-2015 model is always a list object - Quick fix here ...
  ## KMP 19-11-2015 the if (!is.list(model)) solution
  ##                does not work when there is only one model
  ##                which is a list of coefficients, residuals, ...
  if (is.model(model,verbose=verbose)) {
    # if (names(model)[1]=="coefficients") {  ## REB 2016-01-12: changed to the line above.
    y <- predict(model,newdata=newdata) + attr(x,'mean')
  } else {
    #    if (!is.null(newdata)) {
    y <- lapply(model,predict,newdata)
    #    } else {
    #        y <- lapply(model,predict)
    #    }
    y <- matrix(unlist(y),nrow=length(idx),ncol=length(model))
  }
  ##  predict - phase scramble of residual
  ## There is a bug in the following lines, works only if model is not a list object -- need fixes here ...  
  residual <- model$residuals
  if (addnoise) {
    if (verbose) print('add noise')
    l <- length(index(x))
    noise <- matrix(rep(NA,n*l),n,l)
    for (i in 1:n)
      noise[i,] <- FTscramble(noise)
    noise <- zoo(t(noise),order.by(index(x)))
    attr(y,'noise') <- noise
  }
  
  y <- zoo(y,order.by=idx)
  attr(y,'source') <- src
  if (!is.null(residual)) attr(y,'residual.mean') <- mean(residual,na.rm=TRUE)
  if (!is.null(residual)) attr(y,'residual.sd') <- sd(residual,na.rm=TRUE)
  class(y) <- class(x)[-1] ## AM remove 'ds' from output class
  # KMP 2018-12-30: Changed to include model in attrcp. Is there any reason not to copy model to y? 
  y <- attrcp(x, y, ignore=c("name", "n.apps", "appendix", "dimnames","aspect"))
  if (verbose) print('predict.ds.eof complete')
  invisible(y)
}

#' @export
predict.ds.pca <- function(x,newdata=NULL,...,addnoise=FALSE,n=100,verbose=FALSE) {
  if (verbose) print(paste("predict.ds.pca",paste(class(x),collapse='-')))
  if (is.null(newdata)) {
    newdata <- data.frame(coredata(as.eof(x)))
    t <- index(as.eof(x))
  } else {
    if (is.null(names(newdata))) names(newdata) <- paste('X',1:dim(newdata)[2],sep='.')
    t <- index(as.eof(newdata))
    newdata <- data.frame(coredata(as.eof(newdata)))
  }
  ## Fudge fix for when the names of the data are x.1, x.2,,, instead of X.1, X.2,.. 
  if (nchar(names(newdata))[1]!=3) names(newdata) <- paste('X',1:dim(newdata)[2],sep='.')
  names(newdata) <- toupper(names(newdata))
  
  d <- dim(newdata)
  if (is.null(d)) d <- c(length(x),1)
  model <- attr(x,'model')
  y <- lapply(model,predict,newdata)
  y <- matrix(unlist(y),nrow=d[1],ncol=length(model))

  ## Replace 
  y <- zoo(y, order.by=t)
  y <- attrcp(x, y, ignore=c("name", "n.apps", "appendix", "dimnames","aspect"))
  class(y) <- class(x)[-1]
  invisible(y)
}

#' @export
predict.ds.comb <- function(x,newdata=NULL,...,addnoise=FALSE,n=100,verbose=FALSE) {
  ## based on predict.ds.eof function
  stopifnot(!missing(x),inherits(x,"ds"))
  if (verbose) print("predict.ds.comb")
  ## If newdata is set as NULL, reassign it to FALSE
  if (is.null(newdata)) newdata <- FALSE
  
  if (is.logical(newdata)) {
    if (!newdata) {
      X <- attr(x,'eof')
    } else {
      X <- attr(x,'appendix.1') 
    }
  } else {
    X <- newdata ## newdata must be an eof object
  }
  
  if (is.null(X)) neofs <- 1 else neofs <- dim(X)[2]
  model <- attr(x,'model')
  
  ## If newdata is provided as a zoo object or an EOF:
  if (inherits(newdata,'zoo')) Y <- predict.ds.eof(x=x,newdata=newdata,verbose=verbose)
  
  ## In newdata is set as TRUE or FALSE -------------------------------------------------
  if (is.logical(newdata)) {
    # For some reason, the column names of newdata are muddled here,
    # and hence Xnames is used to enforce names 'X.1', 'X.2', 'X.3', ...
    #Xnames <- paste("X.",1:neofs,sep="")
    if (newdata) {
      ## newdata <- data.frame(X=coredata(X))
      n.app <- attr(x,'n.apps') 
    } else {
      ## X <- newdata
      n.app <- 1
    }
    rownm <- rep("",n.app+1)
    rownm[1] <- attr(x,'source')
    
    if (!newdata) {
      Y <- predict.ds.eof(x=x,verbose=verbose)
      attr(Y,'aspect') <- 'fitted'
    } else {
      for (i in 1:n.app) {
        X <- attr(x,paste('appendix.',i,sep=""))
        if (inherits(x,c('pca','eof'))) 
          y <- predict.ds.eof(x=x,newdata=X,verbose=verbose) else
          y <- X
        if (i==1) 
          Y <- y
        else
          Y <- merge(y,Y,all=TRUE)
      }
      attr(Y,'aspect') <- 'predicted'
    }
    
  } ## End of the section dealing with predict if newdata is TRUE/FALSE ----------------
  
  residual <- model$residuals
  if (addnoise) {
    l <- length(residual)
    noise <- matrix(rep(NA,n*l),n,l)
    for (i in 1:n)
      noise[i,] <- FTscramble(noise)
    attr(Y,'noise') <- noise
    names(Y) <- rownm
  }
  Y <- zoo(Y,order.by=index(Y))
  attr(Y,'source') <- attr(X,'source')
  attr(Y,'residual.mean') <- mean(residual,na.rm=TRUE)
  attr(Y,'residual.sd') <- sd(residual,na.rm=TRUE)
  ## KMP 2018-12-30: Changed to include model in attrcp. Is there any reason not to copy model to y? 
  Y <- attrcp(x,Y,ignore = c("name", "n.apps", "appendix", "dimnames","aspect"))
  attr(Y,'aspect') <- 'predicted'
  invisible(Y)
}

#' @export
project.ds <- function(x,newdata=NULL,...,addnoise=FALSE,n=100,verbose=FALSE) {
  ## based on predict.ds.eof function
  stopifnot(!missing(x),inherits(x,"ds")) ## ,inherits(x,'comb')
  if (verbose) print("project.ds")
  
  if (is.null(newdata)) {
    X <- attr(x,'eof')
  } else {
    X <- newdata ## newdata must be an eof object
  }
  
  neofs <- length(attr(X,'eigenvalues'))
  
  # For some reason, the column names of newdata are muddled here,
  # and hence Xnames is used to enforce names 'X.1', 'X.2', 'X.3', ...
  Xnames <- paste("X.",1:neofs,sep="")
  if (verbose) print(Xnames)
  if (is.null(newdata)) {
    n.app <- attr(X,'n.apps')
  } else {
    ## X <- newdata
    n.app <- attr(newdata,'n.apps')
  }
  
  model <- attr(x,'model')
  if (verbose) {print(summary(model)); print(n.app)}
  
  rownm <- rep("",n.app+1)
  rownm[1] <- attr(x,'source')
  for (i in 1:n.app) {
    if(is.null(newdata)) newdata <- attr(X,paste('appendix.',i,sep=""))
    if(verbose) {
      print("Data for obs:")
      print(summary(newdata))
    }
    names(newdata) <- Xnames
    newdata <- attrcp(X,newdata)
    y <- predict.ds.eof(x=x,newdata=newdata,addnoise=FALSE,n=100,verbose=verbose)
    
    names(newdata) <- Xnames
    if (verbose) {
      print("Data for GCM:")
      print(summary(newdata))
    }
    if (verbose) {
      print("DS values:")
      print(summary(predict(model,newdata=newdata)))
    }
    
    if (i==1) {
      Y <- y
    } else {
      Y <- merge(y,Y,all=TRUE)
    }
  }
  
  residual <- model$residuals
  if (addnoise) {
    l <- length(residual)
    noise <- matrix(rep(NA,n*l),n,l)
    for (i in 1:n)
      noise[i,] <- FTscramble(noise)
    attr(Y,'noise') <- noise
  }
  Y <- zoo(Y,order.by=index(Y))
  #print(dim(Y)); print(rownm)
  #print(names(Y))
  names(Y) <- rownm
  attr(Y,'source') <- attr(X,'source')
  attr(Y,'residual.mean') <- mean(residual,na.rm=TRUE)
  attr(Y,'residual.sd') <- sd(residual,na.rm=TRUE)
  # KMP 2018-12-30: Changed to include model in attrcp. Is there any reason not to copy model to y? 
  y <- attrcp(x, y, ignore=c("name", "n.apps", "appendix", "dimnames","aspect"))
  attr(Y,'aspect') <- 'projected'
  Y <- as.station(Y)
  invisible(Y)
}

# To get one predictor pattern, use predict with newdata set to
# a vector where most variables are set to zero apart from one
# variable set to unity for the identification of teleconnection pattern.
#' @export
predict.mvr <- function(x, newdata=NULL, ..., verbose=FALSE) {
  if(verbose) print("predict.mvr")
  object <- x
  if (is.null(newdata)) newdata <- object$data
  x <- newdata
  
  psi <- object$model
  Z <- object$fitted.values
  if (inherits(newdata,'zoo')) {
    Yhat <- zoo(coredata(x) %*% psi,order.by=index(x)) 
  } else if (is.vector(x)) {
    Yhat <- t(x) %*% psi
  }
  
  nattr <- softattr(Z)
  for (i in 1:length(nattr)) {
    attr(Yhat,nattr[i]) <- attr(Z,nattr[i])
  }
  class(Yhat) <- class(object$fitted.values)
  
  if (inherits(Yhat,'eof')) Yhat <- eof2field(Yhat)
  #attr(Yhat,'history') <- c('predict.MVR',attr(Z,'history'))
  #attr(Yhat,'date-stamp') <- date()
  #attr(Yhat,'call') <- match.call()
  attr(Yhat,'history') <- history.stamp(x)
  invisible(Yhat)
}

#' @export
predict.cca <- function(x, newdata=NULL, ..., verbose=FALSE) {
  if(verbose) print("predict.cca")
  if (!is.null(newdata)) X <- newdata else X <- x$X
  #predict.CCA <- function(Psi,X) {
  
  #if ( (class(X)[1]!="eof") & (class(X)[1]!="field")) stop('Need a field or EOF object!')
  #type <- class(X)
  #if (type[1]=="eof") X <- EOF2field(X)
  #X <- field$dat
  #d <- dim(X); dim(X) <- c(d[1],d[2]*d[3])
  #X <- t(X)
  #print(dim(Psi)); print(dim(X)); print(d)
  Y.hat <-  Psi(x) %*% X
  #field$dat <- t(Y.hat)
  #print(dim(field$dat))
  #d1 <- attr(Psi,"dims")
  #dim(field$dat) <- c(d[1],d1[2],d1[3])
  #field$lon <- attr(Psi,"lon"); nx <- length(field$lon)
  #field$lat <- attr(Psi,"lat"); ny <- length(field$lat)
  #field$id.x <- rep("CCA",nx*ny)
  #field$id.lon <- rep("CCA",nx)
  #field$id.lat <- rep("CCA",ny)
  #field$id.t <- rep("CCA",d[1])
  #print("HERE")
  #if (type[1]=="eof") result <- EOF(field) else result <- field
  #result
  Y.hat <- attrcp(Y.hat,X)
  Y.hat
}
