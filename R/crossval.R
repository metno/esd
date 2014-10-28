# Cross-validation:
# Estimate both the correlation, RMSE, R^2, and the coeffieients for
# each iteration - leave-out sample. Collect the statistics on the
# coeffieients - estimate and error - and assess the consistency over
# the iterations.

crossval <- function(x, m=5, ...) UseMethod("crossval")

crossval.ds <- function(x, m=5, ...) {
  # Repeat the regression from DS, but through several iterations with
  # leave-m-out. These are masked by setting them to NA before the
  # regression.
  print("crossval.ds")

  CALDAT <- attr(x,'calibration_data')
  calstr <- attr(CALDAT,'calibration_expression')
  #print(calstr)
  swsm <- attr(CALDAT,'stepwise_screening')
  if (is.null(attr(x,'standard.error'))) weighted <- FALSE
  
  nt <- length(CALDAT$y); nv <- dim(CALDAT)[2]

  if (is.character(m)) {
    print(paste('predefined set-up:',m))
    yr <- switch(tolower(m),
                 'cordex-esd-exp1'=c(1979,1984,1989,1994,1999))
    ii <- (1:length(year(x)))[is.element(year(x),yr)]
    print(year(x)[ii]); print(ii)
    m <- switch(tolower(m),
                 'cordex-esd-exp1'=c(5,5,5,5,5))
  } else if (is.numeric(m) | is.integer(m)) {
     ii <- seq(1,nt,by=m)
     m <- rep(m,length(ii))
  }
  beta <- matrix(rep(NA,nv*length(ii)),nv,length(ii))
  
  Y <- rep(NA,nt); Z <- Y; 
  k <- 0
  #print(ii); print(m)
  for (i in ii) {
    k <- k + 1
    caldat <- as.data.frame(coredata(CALDAT))
    #str(caldat)
    j <- i:(i+m[is.element(ii,i)]-1)
    j <- j[j <= nt] # do not exceed the index
    #print(j)
    Y[j] <- caldat$y[j]
    caldat$y[j] <- NA
    MODEL <- eval(parse(text=calstr))
    # Stepwise regression
    if (!is.null(swsm)) {
      cline <- paste("model <- ",swsm,"(MODEL,trace=0)",sep="")
    #print(paste("HERE: stepping",cline))
      eval(parse(text=cline))
    } else
      model <- MODEL
    # Need to also include the intercept:
    terms <- c(1,as.integer(gsub('X.','',attr(model$terms,'term.labels')))+1)
    beta[terms,k] <- model$coefficients
    #print(terms); print(model$coefficients)
    #print(summary(model))
    z <- predict(model,newdata=caldat)
    Z[j] <- z[j]
  }
  #str(Z);str(Y)
  Y <- Y + attr(x,'mean')
  Z <- Z + attr(x,'mean')
  X <- zoo(cbind(Y,Z),order.by=index(CALDAT))
  attr(X,'location') <- attr(x,'location')
  attr(X,'longitude') <- attr(x,'longitude')
  attr(X,'latitude') <- attr(x,'latitude')
  attr(X,'altitude') <- attr(x,'altitude')
  attr(X,'variable') <- attr(x,'variable')
  attr(X,'original_model') <- attr(x,'model')
  attr(X,'m') <- attr(x,'m')
  attr(X,'unit') <- attr(x,'unit')
  attr(X,'beta') <- beta
  attr(X,'correlation') <- round(cor(Z,Y),2)
  attr(X,'rmse') <- round(sum( (Z - Y)^2 )/nt,2)
  attr(X,'fitted_values_all') <- attr(x,'fitted_values')
  attr(X,'call') <- match.call()
  attr(X,'date') <- date()
  class(X) <- c("xval","zoo")
  invisible(X)
}

crossval.list <- function(x, m=5, ...) {
  elements <- names(x)
  for (i in 1:length(elements)) {
    ds <- x[[i]]
    xval <- crossval.ds(ds)
    attr(x[[i]],'evaluation') <- xval
  }
  invisible(x)
}

