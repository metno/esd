# Predict can take an eof or field, projected onto the EOFs,R then
# apply the DS model.
# Rasmus Benestad

predict.ds <- function(x,newdata=NULL,addnoise=FALSE,n=100,verbose=FALSE) {
  if (verbose) print(paste("predict.ds",paste(class(x),collapse='-')))
  stopifnot(!missing(x),inherits(x,"ds"))
  if ( (inherits(x,'eof')) & (is.null(newdata)) ) {
    if(verbose) print("no new predictor data is provided")
    if(verbose) print("predictand is an EOF")
    if (inherits(x,'comb')) {
      if(verbose) print("predictand is a combined EOF")
      y <- predict.ds.comb(x,newdata=newdata,addnoise=addnoise,
                           n=n,verbose=verbose)
    }
  } else if (inherits(x,'field')) {
    if (verbose) print("predictand is a field object") 
    y <- predict.ds.eof(x,newdata=newdata,addnoise=addnoise,
                        n=n,verbose=verbose)
  } else if (inherits(x,'eof')) {
    if (verbose) print("predictand is an EOF") 
    if(verbose) print("new predictor data is provided")
    if (inherits(newdata,'comb')) {
      if (verbose) print("new predictor data is a combined EOF") 
      y <- predict.ds.comb(x,newdata=newdata,addnoise=addnoise,
                           n=n,verbose=verbose)
    } else if (inherits(newdata,'eof')) {
      if (verbose) print("new predictor data is an EOF") 
      y <- predict.ds.eof(x,newdata=newdata,addnoise=addnoise,
                          n=n,verbose=verbose)
    }
  } else if (inherits(x,'pca')) {
    if(verbose) print("predictand is PCA")
    y <- predict.ds.pca(x,newdata=newdata,addnoise=addnoise,n=n,verbose=verbose)
  }
  
  y <- attrcp(attr(x,'original_data'),y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

predict.ds.eof <- function(x,newdata=NULL,addnoise=FALSE,n=100,verbose=FALSE) {
    stopifnot(!missing(x),inherits(x,"ds"))
  if (verbose) print(paste("predict.ds.eof",paste(class(x),collapse='-')))
  X <- as.eof(x)
  if (verbose) print(paste(class(X),collapse='-'))
  #print(dim(X))
  if (is.null(newdata)) neofs <- length(attr(X,'eigenvalues')) else
                        neofs <- length(attr(newdata,'eigenvalues'))
  
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
      if ((length(attr(X,'eigenvalues'))) != (length(attr(newdata,'eigenvalues'))))
        warning('newdata and X have different number of EOFs')
  }
  #print(summary(newdata))
  names(newdata) <- Xnames 
  model <- attr(x,'model')
  if (verbose) print(summary(model))
  ## AM 04-04-2015 model is always a list object - Quick fix here ...
  ## KMP 19-11-2015 the if (!is.list(model)) solution
  ##                does not work when there is only one model
  ##                which is a list of coefficients, residuals, ...
  #if (!is.list(model)) {
  if (names(model)[1]=="coefficients") {
    y <- predict(model,newdata) + attr(x,'mean')
  } else {
    if (!is.null(newdata)) {
        y <- lapply(model,predict,newdata)
    } else {
        y <- lapply(model,predict)
    }
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
  y <- attrcp(x,y)
  if (verbose) print('predict.ds.eof complete')
  invisible(y)
}


predict.ds.pca <- function(x,newdata=NULL,addnoise=FALSE,n=100,verbose=FALSE) {
  ## REB: modified the code 2015-04-09
  if (verbose) print(paste("predict.ds.pca",paste(class(x),collapse='-')))
  if (is.null(newdata)) {
    newdata <- data.frame(coredata(as.eof(x)))
    t <- index(as.eof(x))
  } else {
    t <- index(as.eof(newdata))
    newdata <- data.frame(coredata(as.eof(newdata)))
  }
  
  d <- dim(newdata)
  if (is.null(d)) d <- c(length(x),1)
  model <- attr(x,'model')
#  browser()
  y <- lapply(model,predict,newdata)
  y <- matrix(unlist(y),nrow=d[1],ncol=length(model))
#  Z <- list()
#  for (i in 1:npca) {
#    y <- zoo(x[,i])
#    attr(y,'model') <- attr(x,'model')[[i]]
#    attr(y,'eof') <- attr(x,'eof')[[i]]
#    attr(y,'eof') <- as.eof(x)
#    attr(y,'mean') <- 0
#    class(y) <- c('ds','eof','zoo')
#    Z[[i]] <- predict.ds.eof(y,newdata=newdata,addnoise=addnoise,n=n,verbose=verbose)
#  }
  ## Copy the original object and only change the predicted values
 
  ## Replace 
  #browser()
  y <- zoo(y, order.by=t)
  y <- attrcp(x,y)
  class(y) <- class(x)[-1]
  invisible(y)
}

predict.ds.comb <- function(x,newdata=NULL,addnoise=FALSE,n=100,verbose=FALSE) {
    ## based on predict.ds.eof function
    stopifnot(!missing(x),inherits(x,"ds"))
  if (verbose) print("predict.ds.comb")

  
  if (is.null(newdata))
      X <- attr(x,'eof')
  else
      X <- newdata ## newdata must be an eof object

  neofs <- length(attr(X,'eigenvalues'))
  
  # For some reason, the column names of newdata is muddled here,
  # and hence Xnames is used to enforce names 'X.1', 'X.2', 'X.3', ...
  Xnames <- paste("X.",1:neofs,sep="")
  if (is.null(newdata)) {
      ## newdata <- data.frame(X=coredata(X))
      n.app <- attr(x,'n.apps') 
  } else {
      ## X <- newdata
      n.app <- attr(newdata,'n.apps')
  }
  
  #print(Xnames)
  
  #print(names(attributes(x)))
  model <- attr(x,'model')
  
  #print(summary(model))
  #print("Data for obs:"); print(summary(newdata))
  ## Y <- zoo(predict(model,newdata=newdata)+attr(x,'mean'),order.by=index(X))
  #print("Y:"); print(summary(coredata(Y)))

  rownm <- rep("",n.app+1)
  rownm[1] <- attr(x,'source')
  
  for (i in 1:n.app) {
      if (is.null(newdata))
          newdata <- attr(X,paste('appendix.',i,sep="")) 
      names(newdata) <- Xnames 
      y <- predict.ds.eof(x=x,newdata=newdata,addnoise=FALSE,n=100)
      
    #print(dim(X))
    #print(names(attributes(X)))
    ## rownm[i+1] <- attr(Z,'source')
    ## newdata <- data.frame(X=coredata(Z))
    #print(summary(newdata))
    names(newdata) <- Xnames
    #print("Data for GCM:"); print(summary(newdata))
    #print("DS values:"); print(summary(predict(model,newdata=newdata)))
    ## y <- zoo(predict(model,newdata=newdata)+ attr(x,'mean'),order.by=index(Z))
      if (i==1) 
          Y <- y
      else
          Y <- merge(y,Y,all=TRUE)
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
  #print("HERE")
    Y <- attrcp(x,Y)
    attr(Y,'aspect') <- 'predicted'
    invisible(Y)
}

project.ds <- function(x,newdata=NULL,addnoise=FALSE,n=100,verbose=FALSE) {
    ## based on predict.ds.eof function
    stopifnot(!missing(x),inherits(x,"ds")) ## ,inherits(x,'comb')
  if (verbose) print("predict.ds.comb")

  
  if (is.null(newdata))
      X <- attr(x,'eof')
  else
      X <- newdata ## newdata must be an eof object

  neofs <- length(attr(X,'eigenvalues'))
  
  # For some reason, the column names of newdata is muddled here,
  # and hence Xnames is used to enforce names 'X.1', 'X.2', 'X.3', ...
  Xnames <- paste("X.",1:neofs,sep="")
  if (is.null(newdata)) {
      ## newdata <- data.frame(X=coredata(X))
      n.app <- attr(x,'n.apps') 
  } else {
      ## X <- newdata
      n.app <- attr(newdata,'n.apps')
  }
  
  #print(Xnames)
  
  #print(names(attributes(x)))
  model <- attr(x,'model')
  
  #print(summary(model))
  #print("Data for obs:"); print(summary(newdata))
  ## Y <- zoo(predict(model,newdata=newdata)+attr(x,'mean'),order.by=index(X))
  #print("Y:"); print(summary(coredata(Y)))

  rownm <- rep("",n.app+1)
  rownm[1] <- attr(x,'source')
  
  for (i in 1:n.app) {
      if (is.null(newdata))
          newdata <- attr(X,paste('appendix.',i,sep="")) 
      names(newdata) <- Xnames 
      y <- predict.ds.eof(x=x,newdata=newdata,addnoise=FALSE,n=100)
      
    #print(dim(X))
    #print(names(attributes(X)))
    ## rownm[i+1] <- attr(Z,'source')
    ## newdata <- data.frame(X=coredata(Z))
    #print(summary(newdata))
    names(newdata) <- Xnames
    #print("Data for GCM:"); print(summary(newdata))
    #print("DS values:"); print(summary(predict(model,newdata=newdata)))
    ## y <- zoo(predict(model,newdata=newdata)+ attr(x,'mean'),order.by=index(Z))
      if (i==1) 
          Y <- y
      else
          Y <- merge(y,Y,all=TRUE)
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
  #print("HERE")
    Y <- attrcp(x,Y)
    attr(Y,'aspect') <- 'projected'
    Y <- as.station(Y)
    invisible(Y)
}

