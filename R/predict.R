# Predict can take an eof or field, projected onto the EOFs,R then
# apply the DS model.
# Rasmus Benestad

predict.ds <- function(x,newdata=NULL,addnoise=FALSE,n=100) {
  stopifnot(!missing(x),inherits(x,"ds"))
  if ( (inherits(x,'eof')) & (is.null(newdata)) ) {
  if (inherits(x,'comb'))
      y <- predict.ds.comb(x,newdata=newdata,addnoise=addnoise,n=n) else
  if (inherits(x,'field'))
      y <- predict.ds.eof(x,newdata=newdata,addnoise=addnoise,n=n)
  } else if (inherits(x,'eof')) {
     if (inherits(newdata,'comb'))
       y <- predict.ds.comb(x,newdata=newdata,addnoise=addnoise,n=n) else
     if (inherits(newdata,'eof'))
       y <- predict.ds.eof(x,newdata=newdata,addnoise=addnoise,n=n)
  }

  y <- attrcp(attr(x,'original_data'),y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

predict.ds.eof <- function(x,newdata=NULL,addnoise=FALSE,n=100) {
  stopifnot(!missing(x),inherits(x,"ds"))
  #print("predict.ds.eof")
  X <- attr(x,'eof')
  W <- attr(x,'eigenvalues')
  #print(dim(X))
  neofs <- length(attr(X,'eigenvalues'))
  
  # For some reason, the column names of newdata is muddled here,
  # and hence Xnames is used to enforce names 'X.1', 'X.2', 'X.3', ...
  Xnames <- paste("X.",1:neofs,sep="")
  if (is.null(newdata)) {
    newdata=data.frame(X=coredata(X))
  } else X <- newdata
  #print(summary(newdata))
  names(newdata) <- Xnames 
  
  model <- attr(x,'model')
  y <- predict(model,newdata) + attr(x,'mean')
  
#  predict - phase scramble of residual
  residual <- model$residuals
  if (addnoise) {
    l <- length(index(x))
    noise <- matrix(rep(NA,n*l),n,l)
    for (i in 1:n)
      noise[i,] <- FTscramble(noise)
    noise <- zoo(t(noise),order.by(index(x)))
    attr(y,'noise') <- noise
  }

  y <- zoo(y,order.by=index(X))
  attr(y,'source') <- attr(X,'source')
  attr(y,'residual.mean') <- mean(residual,na.rm=TRUE)
  attr(y,'residual.sd') <- sd(residual,na.rm=TRUE)
  invisible(y)
}

predict.ds.comb <- function(x,newdata=NULL,addnoise=FALSE,n=100) {
  stopifnot(!missing(x),inherits(x,"ds"))
  #print("predict.ds.comb")
  X <- attr(x,'eof')
  neofs <- length(attr(X,'eigenvalues'))
  
  # For some reason, the column names of newdata is muddled here,
  # and hence Xnames is used to enforce names 'X.1', 'X.2', 'X.3', ...
  Xnames <- paste("X.",1:neofs,sep="")
  if (is.null(newdata)) {
    newdata=data.frame(X=coredata(X))
  } else X <- newdata
  names(newdata) <- Xnames 
  #print(Xnames)
  
  #print(names(attributes(x)))
  model <- attr(x,'model')
  n.app <- attr(x,'n.apps')

  #print(summary(model))
  #print("Data for obs:"); print(summary(newdata))
  Y <- zoo(predict(model,newdata=newdata)+attr(x,'mean'),order.by=index(X))
  #print("Y:"); print(summary(coredata(Y)))

  rownm <- rep("",n.app+1)
  rownm[1] <- attr(x,'source')
    
  for (i in 1:n.app) {
    Z <- attr(X,paste('appendix.',i,sep=""))
    #print(dim(X))
    #print(names(attributes(X)))
    rownm[i+1] <- attr(Z,'source')
    newdata <- data.frame(X=coredata(Z))
    #print(summary(newdata))
    names(newdata) <- Xnames
    #print("Data for GCM:"); print(summary(newdata))
    #print("DS values:"); print(summary(predict(model,newdata=newdata)))
    y <- zoo(predict(model,newdata=newdata)+ attr(x,'mean'),
             order.by=index(Z))
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
  invisible(Y)
}
