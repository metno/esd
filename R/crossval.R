#' Cross-validation
#' 
#' Applies a cross-validation of DS results, using the same strategy as in the
#' DS exercise. Any step-wise screening is applied for each iteration
#' independently of that used to identify the subset of skillful predictors in
#' the original analysis. The model coeffiecients (beta) is saved for each
#' iteration, and both correlation and root-mean-squared-error are returned as
#' scores.
#' 
#' \code{crossval.dsensemble} will make use of the \code{evaluation} attribute
#' with cross-validation results and returns the correlation.
#' 
#' 
#' @aliases crossval crossval.ds crossval.list crossval.dsensemble
#'
#' @param x The results from \code{\link{DS}}.
#' @param m window with - leave m-out for each iteration. There are also some
#' pre-set options: 'cordex-esd-exp1', 'value-exp1', and 'loo' for experiments
#' defined at CORDEX-ESD, COST-VALUE, and leave-one-out ('loo')
#' cross-validation.
#' @param verbose if TRUE print progress
#' @param \dots additional arguments
#'
#' @return Cross-validation object.
#'
#' @keywords manip
#'
#' @examples 
#' data(Oslo)
#' t2m <- t2m.DNMI(lon=c(-20,40),lat=c(45,65))
#' eof <- EOF(t2m)
#' 
#' ds <- DS(Oslo,eof)
#' xv <- crossval(ds)
#' plot(xv)
#' 
#' @export
crossval <- function(x, m=5, verbose=FALSE, ...) UseMethod("crossval")

#' @exportS3Method
#' @export crossval.ds
crossval.ds <- function(x, m=5, verbose=FALSE, ...) {
  # Repeat the regression from DS, but through several iterations with
  # leave-m-out. These are masked by setting them to NA before the
  # regression.
  if(verbose) print("crossval.ds")
  
  CALDAT <- attr(x,'calibration_data')
  calstr <- attr(CALDAT,'calibration_expression')
  #print(calstr)
  swsm <- attr(CALDAT,'stepwise_screening')
  if (is.null(attr(x,'standard.error'))) weighted <- FALSE
  
  nt <- length(CALDAT$y); nv <- dim(CALDAT)[2]
  
  if (is.character(m)) {
    ## The settings for CORDEX-ESD experiment 1 Tier 2:
    ## 5-fold cross-validation: [1979,1983], [1984,1988],[1989,1993],[1994,1998],[1999,2003]
    ## The experiment protocol:
    ## http://wcrp-cordex.ipsl.jussieu.fr/images/pdf/guidelines/CORDEX_ESD_Experiment1.pdf
    if (verbose) print(paste('predefined set-up:',m))
    segments <- switch(tolower(m),
                       'cordex-esd-exp1'= c(1979,1984,1989,1994,1999),
                       'value-exp1'=c(1979,1985,1991,1997,2003),
                       'loo'=year(x))
    ## For winter season Dec-Feb, use the year corresponding to Jan-Dec.
    if (season(x)[1]=='djf') segments <- segments + 1
    yr <-  segments
    ii <- (1:length(year(x)))[is.element(year(x),yr)]
    if (verbose) {
      print(paste('x spans over',min(year(x)),'-',max(year(x))))
      print(year(x)[ii]); print(ii)
    }
    m <- switch(tolower(m),
                 'cordex-esd-exp1'=c(5,5,5,5,5),
                 'value-exp1'=c(6,6,6,6,6),
                 'loo'=rep(1,length(x)))
  } else if (is.numeric(m) | is.integer(m)) {
     ii <- seq(1,nt,by=m)
     m <- rep(m,length(ii))
  }
  beta <- matrix(rep(NA,nv*length(ii)),nv,length(ii))
  
  Y <- rep(NA,nt); Z <- Y; 
  k <- 0
  for (i in ii) {
    k <- k + 1
    caldat <- as.data.frame(coredata(CALDAT))
    j <- i:(i+m[is.element(ii,i)]-1)
    j <- j[j <= nt] # do not exceed the index
    Y[j] <- caldat$y[j]
    caldat$y[j] <- NA
    MODEL <- eval(parse(text=calstr))
    if (!is.null(swsm)) {
      cline <- paste("model <- ",swsm,"(MODEL,trace=0)",sep="")
      eval(parse(text=cline))
    } else {
      model <- MODEL
    }
    terms <- c(1,as.integer(gsub('X.','',attr(model$terms,'term.labels')))+1)
    beta[terms,k] <- model$coefficients
    #print(terms); print(model$coefficients)
    #print(summary(model))
    z <- predict(model,newdata=caldat)
    Z[j] <- z[j]
  }
  
  #str(Z);str(Y)
  if (!inherits(x,'pca')) { 
    Y <- Y + attr(x,'mean')
    Z <- Z + attr(x,'mean')
  }
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

#' @exportS3Method
#' @export crossval.list
crossval.list <- function(x, m=5, verbose=FALSE, ...) {
  if(verbose) print("crossval.list")
  elements <- names(x)
  for (i in 1:length(elements)) {
    ds <- x[[i]]
    xval <- crossval.ds(ds,m=m,verbose=verbose,...)
    attr(x[[i]],'evaluation') <- xval
  }
  invisible(x)
}

#' @exportS3Method
#' @export crossval.dsensemble
crossval.dsensemble <- function(x,m=NULL,verbose=FALSE,...,
                                mar=c(3,4,1.5,0.5),plot=TRUE,xlim=c(-0,1)) {
  if(verbose) print("crossval.dsensemble")
  X <- x
  n <- m
  if (is.null(n)) n <- dim(x$pca)[2]
  #m <- (n-1)%%3+1; k <- (n-1)%/%m+1
  k <- ceiling(n/3)
  m <- ceiling(n/k)
  mfrow=c(m,k)
  X$info <- NULL; X$eof <- NULL; X$pca <- NULL
  xval <- lapply(X,function(x) diag(cor(attr(x,'evaluation'))[seq(2,2*n,by=2),seq(1,2*n-1,by=2)]))
  if (plot) {
    par(mfrow=mfrow,mar=mar)
    for (i in 1:n) { 
      hist(unlist(lapply(xval,function(x) x[i])),col='grey',lwd=2,xlim=xlim,
         main=paste('X-validation correlation for PCA',i),xlab='correlation')
    }
  }
  invisible(xval)
}

