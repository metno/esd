#' Test function for DS.field
#'
#' @param x a \code{ds} \code{eof} object
#' @param verbose a boolean; if TRUE print information on progress
#'
#' @return a \code{field} object with the difference between the original field and the field
#' reconstructed from the independent downscaled principle components from cross-validation
#'
#' @export
test.ds.field <- function(x,verbose=FALSE) {
  if (verbose) print('fieldtest')
  stopifnot (inherits(x,'eof') & inherits(x,'ds'))
  if (verbose) print(colnames(attr(x,'evaluation')))
  isel <- is.element(colnames(attr(x,'evaluation')),paste('X.PCA',1:dim(x)[2],sep='.'))
  eof1 <- attr(x,'evaluation')[,isel]
  isel <- is.element(colnames(attr(x,'evaluation')),paste('Z.PCA',1:dim(x)[2],sep='.'))
  eof2 <- attr(x,'evaluation')[,isel]
  if (verbose) {str(eof1); print(names(attributes(x)))}
  attr(eof1,'variable') <- varid(x)
  attr(eof1,'unit') <- unit(x)
  attr(eof1,'longname') <- attr(x,'longname')
  attr(eof1,'pattern') <- attr(x,'pattern')
  attr(eof1,'eigenvalues') <- attr(x,'eigenvalues')
  attr(eof1,'dimensions') <- c(dim(attr(x,'evaluate'))[1],dim(attr(x,'pattern'))[3])
  attr(eof1,'mean') <- attr(x,'mean')
  attr(eof1,'longitude') <- lon(x)
  attr(eof1,'latitude') <- lat(x)
  attr(eof1,'max.autocor') <- attr(x,'max.autocor')
  attr(eof1,'eigenvalues') <- attr(x,'eigenvalues')
  attr(eof1,'sum.eigenv') <- attr(x,'sum.eigenv')
  attr(eof1,'tot.var') <- attr(x,'tot.var')
  attr(eof1,'aspect') <- 'anomaly'
  attr(eof1,'dimnames') <- NULL   # REB 2016-03-04
  class(eof1) <- c("eof",class(x))[3:6]
  eof2 <- attrcp(eof1,eof2)
  class(eof2) <- c("eof",class(x))[3:6]
  if (verbose) print(c(dim(eof1),dim(eof2)))
  if (verbose) print('estimated EOFs - now get the fields..')
  x1 <- as.field(eof1)
  x2 <- as.field(eof2)
  if (verbose) print(c(dim(x1),dim(x2)))
  coredata(x1) <- coredata(x1 - x2)
  attr(x,'history') <- history.stamp()
  attr(x,'info') <- 'original - downscaled'
  invisible(x1)
}


#' Downscale
#' 
#' Identifies statistical relationships between large-scale spatial climate
#' patterns and local climate variations for monthly and daily data series.
#' 
#' The function calibrates a linear regression model using step-wise screening
#' and common EOFs (\code{\link{EOF}}) as basis functions. It then valuates the
#' statistical relationship and predicts the local climate parameter from
#' predictor fields.
#' 
#' The function is a S3 method that Works with ordinary EOFs, common EOFs
#' (\code{\link{combine}}) and mixed-common EOFs.  DS can downscale results for
#' a single station record as well as a set of stations. There are two ways to
#' apply the downscaling to several stations; either by looping through each
#' station and caryying out the DS individually or by using \code{\link{PCA}}
#' to describe the characteristics of the whole set. Using PCA will preserve
#' the spatial covariance seen in the past. It is also possible to compute the
#' PCA prior to carrying out the DS, and use the method \code{DS.pca}.
#' \code{DS.pca} differs from the more generic \code{DS} by (default) invoking
#' different regression modules (\code{link{MVR}} or \code{\link{CCA}}).
#' 
#' The rationale for using mixed-common EOFs is that the coupled structures
#' described by the mixed-field EOFs may have a more physical meaning than EOFs
#' of single fields [Benestad et al. (2002), "Empirically downscaled
#' temperature scenarios for Svalbard", \emph{Atm. Sci. Lett.},
#' doi.10.1006/asle.2002.0051].
#' 
#' The function \code{DS()} is a generic routine which in principle works for
#' when there is any real statistical relationship between the predictor and
#' predictand. The predictand is therefore not limited to a climate variable,
#' but may also be any quantity affected by the regional climate. \emph{It is
#' important to stress that the downscaling model must reflect a
#' well-understood (physical) relationship.}
#' 
#' The routine uses a step-wise regression (step) using the leading EOFs. The
#' calibration is by default carried out on de-trended data [ref: Benestad
#' (2001), "The cause of warming over Norway in the ECHAM4/OPYC3 GHG
#' integration", \emph{Int. J. Clim.}, 15 March, \bold{vol 21}, p.371-387.].
#' 
#' \code{DS.list} can take a list of predictors and perform a \code{DS} on each
#' of them, seperately, at once. First, \code{DS} is used on the first
#' predictor, then, it is repeated by applying \code{DS} on the residuals from
#' the first step. The DS is repeated for all predictors. The final DS output
#' is list containing as many \code{DS} object as the number of predictors. To
#' get the final DS object, a summation of the different values in the list
#' data object must be done.
#' 
#' \code{DS.seasonalcycle} is an experimental set-up where the calibration is
#' carried out based on the similarity of the seasonal variation to make most
#' use of available information on a 'worst-case' basis, taking the upper limit
#' view that at most, all the seasonal cycle is connected to the corresponding
#' seasonal cycle in the predictor. See Benestad (2009) 'On Tropical Cyclone
#' Frequency and the Warm Pool Area' Nat. Hazards Earth Syst. Sci., 9, 635-645,
#' 2009
#' \url{http://www.nat-hazards-earth-syst-sci.net/9/635/2009/nhess-9-635-2009.html}.
#' 
#' The function \code{biasfix} provides a type of 'bias correction' based on
#' the method \code{\link{diagnose}} which estimates the difference in the mean
#' for the PCs of the calibration data and GCMs over a common period in
#' addition to the ratio of standard deviations and lag-one autocorrelation.
#' This 'bias correction' is described in Imbert and Benestad (2005),
#' \emph{Theor. Appl. Clim.} \url{http://dx.doi.org/10.1007/s00704-005-0133-4}.
#' 
#' 
#' @aliases DS DS.default DS.station DS.list DS.eof DS.comb DS.field
#' DS.t2m.month.field DS.t2m.season.field DS.precip.season.field DS.freq
#' DS.spell DS.pca DS.seasonalcycle biasfix sametimescale
#' @param y The predictand - the station series representing local climate
#' parameter
#' @param X The predictor - an \code{\link{EOF}} object or a list of
#' \code{\link{EOF}} objects representing the large-scale situation.
#' @param method Model type, e.g. \code{\link{lm}} og \code{\link{glm}}
#' @param swsm Stepwise screening, e.g. \code{\link{step}}. NULL skips stepwise
#' screening
#' @param rmtrend TRUE for detrending the predicant and predictors (in the PCs)
#' before calibrating the model
#' @param ip Which EOF modes to include in the model training.
#' @param plot TRUE: plot the results
#' @param verbose TRUE: suppress output to the terminal.
#' @param station TRUE: convert monthly data to station class using
#' \code{\link{combine.ds}}, else return a list of different monthly
#' DS-results.
#' @param area.mean.expl When TRUE, subtract the area mean for the domain and
#' use as a the first co-variate before the PCs from the EOF analysis.
#' @param m passed on to \code{\link{crossval}}. A NULL value suppresses the
#' cross-validation, e.g. for short data series.
#' @param weighted TRUE: use the attribute '\code{error.estimate}' as weight
#' for the regresion analysis.
#' @return The downscaling analysis returns a time series representing the
#' local climate, patterns of large-scale anomalies associated with this,
#' ANOVA, and analysis of residuals. Care must be taken when using this routine
#' to infer local scenarios: check the R2 and p-values to check wether the
#' calibration yielded an appropriate model. It is also important to examine
#' the spatial structures of the large-scale anomalies assocaiated with the
#' variations in the local climate: do these patterns make physical sense?
#' 
#' It is a good idea to check whether there are any structure in the residuals:
#' if so, then a linear model for the relationship between the large and
#' small-scale structures may not be appropriate. It is furthermore important
#' to experiment with predictors covering different regions [ref: Benestad
#' (2001), "A comparison between two empirical downscaling strategies",
#' \emph{Int. J. Climatology}, \bold{vol 21}, Issue 13, pp.1645--1668. DOI
#' 10.1002/joc.703].
#' 
#' There is a cautionary tale for how the results can be misleading if the
#' predictor domain in not appropriate: domain for northern Europe used for
#' sites in Greenland [ref: Benestad (2002), "Empirically downscaled
#' temperature scenarios for northern Europe based on a multi-model ensemble",
#' \emph{Climate Research}, \bold{vol 21 (2)}, pp.105--125.
#' \url{http://www.int-res.com/abstracts/cr/v21/n2/index.html}]
#' @author R.E. Benestad
#' @keywords models multivariate ts spatial
#' @examples
#' 
#' # One exampe doing a simple ESD analysis:
#' X <- t2m.DNMI(lon=c(-40,50),lat=c(40,75))
#' data(Oslo)
#' #X <- OptimalDomain(X,Oslo)
#' eof <- EOF(X,it='jan')
#' Y <- DS(Oslo,eof)
#' plot(Y, new=FALSE)
#' str(Y)
#' 
#' # Look at the residual of the ESD analysis
#' y <- as.residual(Y)
#' plot.zoo(y,new=FALSE)
#' 
#' # Check the residual: dependency to the global mean temperature?
#' T2m <- t2m.DNMI()
#' yT2m <- merge.zoo(y,T2m)
#' plot(coredata(yT2m[,1]),coredata(yT2m[,2]))
#' 
#' # Example: downscale annual wet-day mean precipitation -calibrate over
#' # part of the record and use the other part for evaluation.
#' T2M <- as.annual(t2m.NCEP(lon=c(-10,30),lat=c(50,70)))
#' cal <- subset(T2M,it=c(1948,1980))
#' pre <- subset(T2M,it=c(1981,2013))
#' comb <- combine(cal,pre)
#' X <- EOF(comb)
#' data(bjornholt)
#' y <- as.annual(bjornholt,FUN="exceedance")
#' z <- DS(y,X)
#' plot(z, new=FALSE)
#' 
#' ## Example on using common EOFs as a framework for the downscaling:
#' lon <- c(-12,37)
#' lat <- c(52,72)
#' ylim <- c(-6,6)
#' t2m <- t2m.NCEP(lon=lon,lat=lat)
#' T2m <- t2m.NorESM.M(lon=lon,lat=lat)
#' data(Oslo)
#' X <- combine(t2m,T2m)
#' eof <- EOF(X,it='Jul')
#' ds <- DS(Oslo,eof)
#' plot(ds)
#' 
#' ## Example downscaling statistical parameters: mean and standard deviation
#' ## using different predictors
#' data(ferder)
#' t2m <- t2m.NCEP(lon=c(-30,50),lat=c(40,70))
#' slp <- slp.NCEP(lon=c(-30,50),lat=c(40,70))
#' T2m <- as.4seasons(t2m)
#' SLP <- as.4seasons(slp)
#' X <- EOF(T2m,it='Jan')
#' Z <- EOF(SLP,it='Jan')
#' y <- ferder
#' sametimescale(y,X) -> z
#' ym <- as.4seasons(y,FUN="mean")
#' ys <- as.4seasons(y,FUN="sd")
#' dsm <- DS(ym,X)
#' plot(dsm)
#' dss <- DS(ys,Z)
#' plot(dss)
#' 
#' ## Example for downscaling with missing data
#' data(Oslo)
#' dnmi <- t2m.DNMI(lon=c(-10,20),lat=c(55,65))
#' y <- subset(Oslo,it='jan')
#' X <- EOF(subset(dnmi,it='jan'))
#' ds <- DS(y,X)
#' plot(ds) # Looks OK
#' # Now we replace some values of y with missing data:
#' y2 <- y
#' set2na <- order(rnorm(length(y)))[1:50]
#' y2[set2na] <- NA
#' ds2 <- DS(y2,X)
#' plot(ds2) 
#' 
#' ## Use downscale results to fill in missing data:
#' y3 <- predict(ds2,newdata=X)
#' ## Plot a subset of y based on dates in predicted y3
#' plot(subset(y,it=range(index(y3))),col='grey80',lwd=4,map.show=FALSE)
#' points(as.station(predict(ds2)))
#' # The downscaled 
#' lines(y3,lty=2)
#' 
#' 
#' @export
DS <- function(y,X,verbose=FALSE,plot=FALSE,...) UseMethod("DS")

                                        # The basic DS-function, used by other methods
                                        # REB HERE!

#' @export
DS.default <- function(y,X,verbose=FALSE,plot=FALSE,...,it=NULL,
                       method="lm",swsm="step",m=5,rmtrend=TRUE,ip=1:7,weighted=TRUE) {
    if (verbose) { print('--- DS.default ---'); print(summary(coredata(y)))}
    #print('err(y)'); print(err(y))
    if (verbose) {print('index(y)'); print(index(y))}
    if (verbose) {print(class(y)); print(class(X))}
    swapped <- FALSE
    if ( inherits(y,c("eof")) & inherits(X,c("station"))) {
      if (verbose) print('SWAP y & X')
        yy <- X
        X <- y
        y <- yy
        swapped <- TRUE
    }
    stopifnot(!missing(y),!missing(X), is.matrix(X),
              inherits(X,c("eof","field")),inherits(y,"station"))
    if(!inherits(X,"eof")) X <- EOF(X)
    if (class(index(y)) != (class(index(X)))) {
      warning(paste('DS.default: different indices:', class(index(y)),class(index(X))))
      if (is.numeric(index(y))) index(X) <- year(X)
      if (is.numeric(index(X))) index(y) <- year(y)
    }
    
    y0 <- y
    X0 <- X
    ip <- ip[ip <= length(attr(X,'eigenvalues'))]
    
    if (verbose) {print(paste(sum(!is.finite(coredata(y))),'missing values in y'))}
    if (verbose)  {print('index and y before removing missing values:'); print(zoo(y))}
    y <- subset(y,it=is.finite(coredata(y)))
    W <- attr(X,'eigenvalues')
    cls <- c(class(y)[1],class(X))
 
    if (verbose) print('Ensure matching time scale')
    if (verbose) {print('index(y) before sametimescale:'); print(index(y))}
    y <- sametimescale(y,X,verbose=verbose)
    ## REB is needed to ensure that y is annual if X is annual
    #if (verbose) {print('index(y) after sametimescale:'); print(index(y))}
    #if (verbose) print('Match dates')
    
    y <- matchdate(y,X,verbose=verbose) ##
    #if (verbose) {print('index(y) after matchdate'); print(index(y))}
    X <- matchdate(X,y,verbose=verbose) # REB 2015-01-14
    if (verbose) {print("index(y) & index(X) after synch:");
                  print(index(y)); print(index(X))}
    
    if (!is.null(it)) y <- subset(y,it=it)
    
    ## synchronise the series: use the 'zoo' merge function through combine:
    ##print(index(y)[1:24]); print(index(X)[1:24]);

    ##
    
    #yX <- combine.station.eof(y,X)
    #y <- yX$y; X <- yX$X
    X0 <- X; y0 <- y
    #year <- as.numeric( format(index(y), '%Y') ) 
    #month <- as.numeric( format(index(y), '%m') )
    year <- year(y)
    month <- month(y)
    ##print(length(y)); print(table(year)); print(table(month))
    
    ## De-trend the data used for model calibration:
    if (rmtrend) {
        if (verbose) print('detrend')
        offset <- mean(y,na.rm=TRUE)
        y <- trend(y,result="residual")
        offset <- offset - mean(y,na.rm=TRUE)
        X <- trend(X,result="residual")
    } else offset <- 0

    ##str(y); print(class(y))

    ## REB: 2014-10-03: add weights if available
    weights <- rep(1,length(y))
    if (!is.null(attr(y,'standard.error'))) {
      if (sum(is.finite(attr(y,'standard.error')))>0) weights <- 1/coredata(attr(y,'standard.error'))
      weights[!is.finite(weights)] <- 0
      if (is.null(attr(y,'standard.error'))) weighted <- FALSE
      if (verbose) {print(paste('weights',weighted)); print(weights)}
    }
    #
    ##if (length(index(X)) == length(index(y)))
    caldat <- data.frame(y=coredata(y),X=as.matrix(coredata(X)),
                           weights=weights) 
    
    predat <- data.frame(X=as.matrix(coredata(X0)))
    colnames(predat) <- paste("X",1:ncol(predat),sep=".")#length(colnames(predat)),sep=".")
    if (is.null(names(X))) names(X) <- 1:dim(X)[2]
    Xnames <- paste("X.",1:length(names(X)),sep="")
    colnames(caldat) <- c("y",Xnames,'weights')
    Xnames <- Xnames[ip]
    ## REB 2014-10-03:
    if (weighted) {
      calstr <- paste(method,"(y ~ ",paste(Xnames,collapse=" + "),
                      ", weights=weights, data=caldat, ...)",sep="") 
    } else {
      calstr <- paste(method,"(y ~ ",paste(Xnames,collapse=" + "),
                      ", data=caldat, ...)",sep="")
    }
    
    MODEL <- eval(parse(text=calstr))
    FSUM <- summary(MODEL)
    if (verbose) print(FSUM)

    ## Stepwise regression
    if (!is.null(swsm)) {
      cline <- paste("model <- ",swsm,"(MODEL,trace=0)",sep="")
      eval(parse(text=cline))
    } else {
      model <- MODEL
    }
    terms1 <- attr(model$terms,'term.labels')

    if (verbose) print(summary(model))
    fsum <- summary(model)
    COEFS=FSUM$coefficients
    COEFS[,1] <- 0; 
    COEFS[,2:4] <- NA; 
    coefs=fsum$coefficients
    TERMS <- attr(FSUM$terms,'term.labels')
    terms <- attr(fsum$terms,'term.labels')
    ii <- is.element(attr(COEFS,"dimnames")[[1]],attr(coefs,"dimnames")[[1]])
    COEFS[ii,1] <- coefs[,1]
    COEFS[ii,2] <- coefs[,2]
    COEFS[ii,3] <- coefs[,3]
    COEFS[ii,4] <- coefs[,4]
    dc <- dim(COEFS)

    U <- attr(X,'pattern'); du <- dim(U)
### REB 2015-01-19: Also allow for patterns consisting of vectors - weights of mixed predictors
### See DS.list.
    if (verbose) {print('pattern dimension'); print(du); str(U)}
    if (length(du)==3) dim(U) <- c(du[1]*du[2],du[3])
    if (!is.null(du)) {
      pattern <- t(COEFS[2:dc[1],1]) %*%
          diag(attr(X,'eigenvalues')[ip]) %*% t(U[,ip])
      dim(pattern) <- c(du[1],du[2]) 
    } else pattern <- c(COEFS[2:dc[1],1]) * attr(X,'eigenvalues')[ip]
                                                 
    ##  ds <- zoo(predict(model),order.by=index(X)) + offset
    ##  ds <- zoo(predict(model,newdata=caldat),order.by=index(X)) + offset
    if (verbose) print('predict')
    ds <- zoo(predict(model,newdata=predat),order.by=index(X)) + offset
                                        #plot(y,lwd=4); lines(ds,col="red",lwd=2)
                                        #lines(zoo(model$fitted.values,order.by=index(ds)),col="blue",lwd=2)
                                        #lines(yX$y,col="darkgreen",lwd=2)
    
                                        #nattry <- softattr(y0)
                                        #nattrX <- softattr(X0,ignore=c('longitude','latitude')) 
                                        #for (i in 1:length(nattry))
                                        #  attr(ds,nattry[i]) <- attr(y0,nattry[i])
                                        #for (i in 1:length(nattrX))
                                        #  attr(pattern,nattrX[i]) <- attr(X0,nattrX[i])
    if (verbose) {
      print('set attibutes')
      print(names(attributes(y0)))
      print(names(attributes(X0)))
    }
    ## Make sure that predictions have the same index class (time units) as the original data 
    if (class(index(ds)) != class(index(y0))) {
      ## For annual data:
      if ( (class(index(ds))=='Date') & (class(index(y0))=='numeric') & inherits(y0,'annual') ) 
        index(ds) <- year(index(ds))
      if ( (class(index(ds))=='numeric') & (class(index(y0))=='Date') & inherits(y0,'annual') ) 
        index(ds) <- as.Date(paste(index(ds),'01-01',sep='-'))
    }
    ds <- attrcp(y0,ds,ignore='names')
    pattern <- attrcp(X0,pattern,ignore=c('longitude','latitude','names','dimnames'))
    caldat <- zoo(as.matrix(caldat),order.by=index(X))
    if (verbose) print('Set attributes')
    attr(caldat,'calibration_expression') <- calstr
    attr(caldat,'stepwise_screening') <- swsm
    attr(ds,'calibration_data') <- caldat
    attr(ds,'fitted_values') <- zoo(model$fitted.values +
                                    offset,order.by=index(X))
    class(attr(ds,'fitted_values')) <- class(y0)
    attr(ds,'original_data') <- y0
    r2 <- var(coredata(model$fitted.values))/var(y,na.rm=TRUE)
    attr(r2,'description') <- ' var(fitted.values))/var(y)'
    attr(ds,'quality') <- r2
    attr(ds,'variable') <- attr(y0,'variable')
    attr(ds,'model') <- model
    attr(ds,'mean') <- offset
    attr(ds,'method') <- method
    attr(ds,'eof') <- X0
    attr(pattern,'longitude') <- attr(X0,'longitude')
    attr(pattern,'latitude') <- attr(X0,'latitude')
    attr(ds,'pattern') <- pattern
    attr(ds,'dimensions') <- c(du[1],du[2])
    attr(ds,'longitude') <- attr(y0,'longitude')
    attr(ds,'latitude') <- attr(y0,'latitude')
    ##attr(ds,'source') <- paste(attr(y0,'source'),attr(X0,'source'),sep="-")
    attr(ds,'aspect') <- 'downscaled'
    attr(ds,'source') <- attr(X0,'source')
    attr(ds,'type') <- "downscaled results"
    attr(ds,'history.predictand') <- attr(y0,'history')
    attr(ds,'history') <- history.stamp(X0)
                                        #print("HERE"); print(cls)
    class(ds) <- c("ds",cls[-2])
    ## KMP 2019-04-29: Added crossval in DS.default. Any reason why it shoudn't be here?
    if (!is.null(m))  {
      if (verbose) print("Cross-validation")
      xval <- crossval(ds,m=m)
      attr(ds,'evaluation') <- zoo(xval)
    } else attr(ds,'evaluation') <- NULL
    rm("y0","X0")
                                        #print("Completed")
                                        #lines(ds,col="darkred",lwd=2,lty=2)
                                        #lines(attr(ds,'original_data'),col="green",lwd=2,lty=2)
    if (verbose) print('--- exit DS.default ---')
    if (plot) plot(ds)
    invisible(ds)
}




#' @export
DS.station <- function(y,X,verbose=FALSE,plot=FALSE,...,it=NULL,biascorrect=FALSE,
                       method="lm",swsm="step",m=5,rmtrend=TRUE,ip=1:7,weighted=TRUE,pca=FALSE,npca=20) {
    
    stopifnot(!missing(y),!missing(X),inherits(y,"station"))
    if (verbose) { print('--- DS.station ---'); print(summary(coredata(y)))}
    #print('err(y)'); print(err(y))
    #print('index(y)'); print(index(y))
    y <- matchdate(y,X)
    X <- matchdate(X,y)
    
    ## Used for extracting a subset of calendar months
    if (!is.null(it)) {
      if (verbose) print(paste('it=',it))
      if (inherits(y,'station')) y <- subset(y,it=it)
    }
    
     if ( (!inherits(y,'seasonalcycle')) & (inherits(X,'seasonalcycle')) ) {
                                        #print("HERE")
        ds <- DS.seasonalcycle(y=y,X=X,ip=ip,verbose=verbose,...) 
        return(ds)
    }
    
    if ( (!inherits(X,'eof')) & (inherits(X,'field')) ) {
                                        #print("HERE")
      ds <- DS.field(y=y,X=X,biascorrect=biascorrect,
                     method=method,swsm=swsm,m=m,
                     rmtrend=rmtrend,ip=ip,verbose=verbose,
                     weighted=weighted,pca=pca,npca=npca,...) 
      return(ds)
    } else if (is.list(X)) {
      if (verbose) print("The predictor is a list")
      ds <- DS.list(y=y,X=X,biascorrect=biascorrect,
                    method=method,swsm=swsm,m=m,
                    rmtrend=rmtrend,ip=ip,verbose=verbose,
                    weighted=weighted,pca=pca,npca=npca,...) 
      return(ds)
    } 

    ## REB inserted lines to accomodate for multiple stations
    d <- dim(y)
    if (!is.null(d)) ns <- d[2] else ns <- 1
    if (!is.null(d)) {
      if (verbose) print(paste('predictand contains',ns,'stations'))
      ## More than one station
      Y <- y
      if (pca) {
        if (verbose) print("PCA")
        ## PCA is used when y represents a group of variables and when
        ## their covariance is a concern: it preserves the covariance
        ## as seen in the past, as the PCs and eigenvectors are orthogonal.
        ## Useful for spatial response, wind vectors, temperature/precipitation
        Y <- PCA(y,npca)
      }
      ns <- dim(Y)[2]
    } else {
      Y <- y
    }

    ## Loop over the different PCs or stations
    for (i in 1:ns) {
        if (verbose) print(paste("i=",i))
        z <- subset(Y,is=i)
        if (verbose) {print(class(z)); print(names(attributes(z)))}
        if (inherits(X,'eof')) {
          if (verbose) print("The predictor is some kind of EOF-object")
          ## Call different functions, depending on the class of X:
          #print("the predictor is an EOF-object")
          if (inherits(X,'comb')) {
            if (verbose) print("*** Comb ***")
            ## X is combined EOFs
            ds <- DS.comb(y=z,X=X,biascorrect=biascorrect,
                          method=method,swsm=swsm,it=it,
                          rmtrend=rmtrend,ip=ip,verbose=verbose,...)
            if (verbose) print("---")
          } else {
            if (verbose) print("*** EOF ***")
            ## X is ordinary EOF
            ds <- DS.default(y=z,X=X,
                             method=method,swsm=swsm,
                             rmtrend=rmtrend,ip=ip,verbose=verbose,...)
            if (verbose) print("+++")
          }
        } else if (inherits(X,'field')) {
            if (verbose) print("the predictor is a field-object")
            ## X is a field
            ds <- DS.field(y=z,X=X,biascorrect=biascorrect,
                           method=method,swsm=swsm,
                           rmtrend=rmtrend,ip=ip,verbose=verbose,...)
        }
        ## May need an option for coombined field: x is 'field' + 'comb'
        #if (is.null(ds)) browser()
        
        ## Unless told not to - carry out a cross-validation
        if (!is.null(m))  {
          if (verbose) print("Cross-validation")
          xval <- crossval(ds,m=m)
          attr(ds,'evaluation') <- zoo(xval)
        } else attr(ds,'evaluation') <- NULL

        if (verbose) print(names(attributes(ds)))
        if (ns==1) dsall <- ds else {
            if (i==1) dsall <- list(ds.1=ds) else
            eval(parse(text=paste('dsall$ds.',i,' <- ds',sep='')))
        }

    }
    
    ## If PCA was used to transform the predictands to preserve the
    ## spatial covariance, then do the inverse to recover the results
    ## in a structure comparable to the original stations.
    if ( (ns>1) & pca ) {
      attr(dsall,'pattern') <- attr(Y,'pattern')
      attr(dsall,'eigenvalues') <- attr(Y,'eigenvalues')    
      attr(dsall,'mean') <- attr(Y,'mean')    
      ds.results <- pca2station(dsall)
    } else {
      ds.results <- dsall
    }
    
    if (verbose) print("--- exit DS.station ---")
    if (plot) plot(ds.results)
    invisible(ds.results)  
}


## DS for combined fields - to make predictions not based on the
## calibration data
#' @export
DS.comb <- function(y,X,verbose=FALSE,plot=FALSE,...,biascorrect=FALSE,
                    method="lm",swsm="step",m=5,rmtrend=TRUE,it=NULL,ip=1:7, 
                    weighted=TRUE,pca=FALSE,npca=20) {
    if (verbose) { print('--- DS.comb ---'); print(summary(coredata(y)))}
    ##print('index(y)'); print(index(y))
    ##print('err(y)'); print(err(y))
    if ( inherits(y,c("eof")) & inherits(X,c("station"))) {
        yy <- X
        X <- y
        y <- yy
        swapped <- TRUE
        rm("yy")
    }
    X0 <- X
    
    stopifnot(!missing(y),!missing(X),is.matrix(X),
              inherits(X,"comb"),inherits(y,"station"))
    
    ## For combined fields/common EOFs, do the DS-fitting once, and then
    ## use the model n times to predict the values associated with the
    ## appended fields:
    
    if (class(index(y)) != (class(index(X)))) {
      warning(paste('DS.comb: different indices:', class(index(y)),class(index(X))))
      if (is.numeric(index(y))) index(X) <- year(X)
      if (is.numeric(index(X))) index(y) <- year(y)
    }
    
    if (!inherits(X,"eof")) X <- EOF(X,it=it)
    
    if (biascorrect) {
        if (verbose) print("Bias correction - bias-fix common EOF")
        X <- biasfix(X)
    }
    
    ds <- DS.default(y,X,method=method,swsm=swsm,m=m,
                     rmtrend=rmtrend,ip=ip,weighted=weighted,verbose=verbose,...)

    ## For combined fields, make sure to add the appended PCs to
    ## the results.
                                        #print(names(attributes(X)))
    model <- attr(ds,'model')
    n.app <- attr(X,'n.apps')
    attr(ds,'n.apps') <- n.app
    if (verbose) {print(paste('For combined fields, n.apps=',n.app));
                  print(names(attributes(X))[grep('app',names(attributes(X)))]) }
    for (i in 1:n.app) {
                                        #print("HERE")
        Z <- eval(parse(text=paste("attr(X0,'appendix.",i,"')",sep="")))
                                        #
        newdata <- data.frame(X=Z)
        colnames(newdata) <- paste("X",1:ncol(X),sep=".") 
                                        #print(summary(newdata))
        z <- predict(model,newdata=newdata) + attr(ds,'mean')
        Y <- zoo(z,order.by=index(Z))
        Y <- attrcp(Z,Y)

        # Store in a similar structure as in the common EOF:
        eval(parse(text=paste("attr(ds,'appendix.",i,"') <- Y",sep="")))
    }
    rm("X0"); gc(reset=TRUE)
    if (plot) plot(ds)
    invisible(ds)
}


## This function takes care of downscaling based on a field-object X.
## X can be a combined field. This function calls more primitive DS methods,
## depending on the time scale represented in X (monthly or seasonal).
#' @export
DS.field <- function(y,X,verbose=FALSE,plot=FALSE,...,biascorrect=FALSE,
                     method="lm",swsm="step",m=5,rmtrend=TRUE,ip=1:7,weighted=TRUE) {
    if (verbose) { print('--- DS.field ---'); print(summary(coredata(y)))}
    ## Keep track of which is an eof object and which is a station record:
    swapped <- FALSE
    if ( inherits(y,c("field")) & inherits(X,c("station"))) {
        yy <- X
        X <- y
        y <- yy
        swapped <- TRUE
    }
    if (class(index(y)) != (class(index(X)))) {
      warning(paste('DS.field: different indices:', class(index(y)),class(index(X))))
      if (is.numeric(index(y))) index(X) <- year(X)
      if (is.numeric(index(X))) index(y) <- year(y)
    }
                                        #print(class(X)); print(attr(y,'variable'))
    if (verbose) {print(class(X)); print(varid(y))}
    if (sum(is.element(tolower(attr(y,'variable')),c('t2m','tmax','tmin'))) >0) {
      if (inherits(X,'month')) {
        ds <- DS.default(y=y,X=X,method=method,swsm=swsm,m=m,
                         rmtrend=rmtrend,ip=ip,
                         verbose=verbose)
        #ds <- DS.t2m.month.field(y=y,X=X,biascorrect=biascorrect,
        #                         method=method,swsm=swsm,m=m,
        #                         rmtrend=rmtrend,ip=ip,
        #                         verbose=verbose)
      } else if (inherits(X,'season')) {
        # the DS.t2m.season.field is not in working order
        #ds <- DS.t2m.season.field(y=y,X=X,biascorrect=biascorrect,
        #                          method=method,swsm=swsm,m=m,
        #                          rmtrend=rmtrend,ip=ip,
        #                          verbose=verbose)
        ds <- DS.default(y=y,X=X,method=method,swsm=swsm,m=m,
                         rmtrend=rmtrend,ip=ip,verbose=verbose)
      } else if (inherits(X,'annual')) {
        ds <- DS.default(y=y,X=X,method=method,swsm=swsm,m=m,
                         rmtrend=rmtrend,ip=ip,verbose=verbose)
        #ds <- DS.t2m.annual.field(y=y,X=X,biascorrect=biascorrect,
        #                          method=method,swsm=swsm,m=m,
        #                          rmtrend=rmtrend,ip=ip,
        #                          verbose=verbose)
      }
    } else if (tolower(attr(y,'variable'))=='precip') {
      ds <- DS.default(y=y,X=X,method=method,swsm=swsm,m=m,
                       rmtrend=rmtrend,ip=ip,verbose=verbose)
      #ds <- DS.precip.season.field(y=y,X=X,biascorrect=biascorrect,
      #                             method=method,swsm=swsm,m=m,
      #                             rmtrend=rmtrend,ip=ip,
      #                             verbose=verbose)
    } else {
      ds <- DS.default(y=y,X=X,method=method,swsm=swsm,m=m,
                       rmtrend=rmtrend,ip=ip,verbose=verbose)
    }
    if (verbose) print('return downscaled results')
    if (plot) plot(ds)
    invisible(ds)
}


## Downscales all 12-calendar months for a field by stepping through the months
## and compute the EOFs before applying the default DS method.
# DS.t2m.month.field <- function(y,X,biascorrect=FALSE,it=NULL,
#                                method="lm",swsm="step",m=m,
#                                rmtrend=TRUE,ip=1:7,
#                                verbose=FALSE,weighted=TRUE,station=TRUE) {
#   if (verbose) { 
#     print('--- DS.t2m.month.field ---')
#     print(summary(coredata(y)))
#   }
#   if (inherits(X,'comb')) {
#     type <- 'eof.comb' 
#   } else {
#     type <- "eof.field"
#   }
#   cls <- class(y)
#   ds <- list()
#   if (is.null(it)) X <- subset(X,it=it)
#   mon <-  (1:12)[is.element(month.abb,rownames(table(month(X))))]
#   
#   for (i in mon) {
#     eof <- EOF(X,it=month.abb[i])
#     if (biascorrect) eof <- biasfix(eof)
#     cline <- paste("ds$",month.abb[i],
#                    "<- DS.station(y,eof,method=method,swsm=swsm,m=m,",
#                    "rmtrend=rmtrend,ip=ip,verbose=verbose)",
#                    sep="")
#     if (verbose) print(cline)
#     eval(parse(text=cline))
#   }
# 
#   if (station) {
#     ds <- combine.ds(ds) 
#   } else {
#     cls <- c("list","dsfield",cls)
#   }
# 
#   attr(ds,'predictand.class') <- cls
#   invisible(ds)
# }
# 
# 
# DS.t2m.season.field <- function(y,X,biascorrect=FALSE,
#                                 method="lm",swsm="step",m=5,
#                                 rmtrend=TRUE,ip=1:7,
#                                 verbose=FALSE,weighted=TRUE,station=TRUE) {
#   ## Downscale seasonal mean and standard deviation
#     if (verbose) { print('--- DS.t2m.season.field ---'); print(summary(coredata(y)))}
# 
#     Z1 <- EOF(subset(X,it='djf'))
#     if (verbose) print("downscale DJF")
#     ds1 <- DS(y,Z1,biascorrect=biascorrect,ip=ip)
#     Z2 <- EOF(subset(X,it='mam'))
#     if (verbose) print("downscale MAM")
#     ds2 <- DS(y,Z2,biascorrect=biascorrect,ip=ip)
#     Z3 <- EOF(subset(X,it='jja'))
#     if (verbose) print("downscale JJA")
#     ds3 <- DS(y,Z3,biascorrect=biascorrect,ip=ip)
#     
#     
#     Z4 <- EOF(subset(X,it='son'))
#     if (verbose) print("downscale SON")
#     ds4 <- DS(y,Z4,biascorrect=biascorrect,ip=ip)
#     if (verbose) print("Combine the 4 seasons")
#     ds <- combine(list(ds1,ds2,ds3,ds4))
#     z <- c(crossval(ds1,m=m),crossval(ds2,m=m),
#            crossval(ds3,m=m),crossval(ds4,m=m))
#     attr(ds,'evaluation') <- z
#     save(file='inside.ds.seas.f.1.rda',y,Z1,ds1,Z2,ds2,Z3,ds3,Z4,ds4,X)
#                                         #
#     invisible(ds)
# }
# 
# DS.t2m.annual.field <- function(y,X,biascorrect=FALSE,
#                                 method="lm",swsm="step",m=5,
#                                 rmtrend=TRUE,ip=1:7,
#                                 verbose=FALSE,weighted=TRUE,station=TRUE) {
#   ## Downscale seasonal mean and standard deviation
#     if (verbose) { print('--- DS.t2m.annual.field ---'); print(summary(coredata(y)))}
# 
#     
#     Z <- EOF(annual(X))
#     ds <- DS(annual(y),Z,biascorrect=biascorrect)
#     invisible(ds)
# }
# 
# 
# 
# DS.precip.season.field <- function(y,X,biascorrect=FALSE,threshold=1,
#                                    method="lm",swsm="step",m=5,
#                                    rmtrend=TRUE,ip=1:7,
#                                    verbose=FALSE,weighted=TRUE,...) {
# 
#   ## Computes the annual mean values for wet-day mean mu, wet-day frequency, and spell.
#   ## Also computes seasonal variations from PCA X[year,calendar months].
#   ## One PC for each year.
# 
#     if (verbose) { print('--- DS.precip.season.field ---'); print(summary(coredata(y)))}
#     mu <- as.4seasons(y,FUN="exceedance",threshold=threshold)
#     fw <- as.4seasons(y,FUN="exceedance",fun="freq")
#     wL <- as.4seasons(spell(y))
# 
#     if (!inherits(X,'season')) X <- as.4seasons(X)
# 
#     for (i in 1:4) {
#         x <- EOF(X,it=c('djf','mam','jja','son')[i])
#         if (biascorrect) x <- biasfix(x)
#         ds.mu <- DS.default(mu,x,method=method,swsm=swsm,m=m,
#                             rmtrend=rmtrend,ip=ip,
#                             verbose=verbose,...)
#         ds.fw <- DS.freq(fw,x,rmtrend=rmtrend,ip=ip,m=m,
#                          verbose=verbose,...)
#         ds.wL <- DS.spell(wL,x,rmtrend=rmtrend,ip=ip,m=m,
#                           verbose=verbose,...)
#     }
#     
#     ## DS mu, fw, spell, total
#     ds <- NULL
#     invisible(ds)
# }
# 
# ## Use family='gaussian' for sample sizes gt 30 - > central limit theorem
# DS.freq <- function(y,X,threshold=1,biascorrect=FALSE,method="glm",
#                     family="gaussian",swsm="step",m=5,
#                     rmtrend=TRUE,ip=1:7,verbose=FALSE,weighted=TRUE,...) {
#     if (verbose) { print('--- DS.freq ---'); print(summary(coredata(y)))}
#     if (inherits(X,'month'))
#         Z <- aggregate(y,as.yearmon,FUN="wetfreq",threshold=threshold) else
#     if (inherits(X,'season'))
#         Z <- as.4seasons(y,FUN=wetfreq,threshold=threshold) else
#     if (inherits(X,'annual'))
#         Z <- annual(y,FUN=wetfreq,threshold=threshold)
#     
#     ds <- DS.default(Z,X,method=method,swsm=swsm,m=m,rmtrend=rmtrend,ip=ip,verbose=verbose,...)
#     return(ds)
# }
# 
# 
# DS.spell <- function(y,X,threshold=1,biascorrect=FALSE,
#                      method="glm",family="gaussian",swsm="step",m=5,
#                      rmtrend=TRUE,ip=1:7,verbose=FALSE,weighted=TRUE,...) {
#   ## Downscale the spell length using a GLM with poisson family.
#   ##  the mean spell length over a given interval:
#     if (verbose) { print('--- DS.spell ---'); print(summary(coredata(y)))}
#     if (inherits(y,'spell')) z <- as.station(y) else
#     if (inherits(y,'sstation')) {
#         z <- as.station(spell(y))
#     }
#     
#     if (inherits(X,'month')) Z <- aggregate(z,as.yearmon,FUN=mean) else
#     if (inherits(X,'season')) Z <- as.4seasons(z,FUN=mean) else
#     if (inherits(X,'annual')) Z <- annual(z,FUN=mean)
# 
#     ds <- DS(Z,X,biascorrect=biascorrect,method=method,swsm=swsm,m=m,
#              rmtrend=rmtrend,ip=ip,
#              verbose=verbose,...)
#     invisible(ds)
# }


## DS.pca
## This function applies to PCA results for a group of stations.
## Typically some annual statistics (e.g. mu), stored in compressed form
## as a PCA-object (similar to EOF, but not gridded and without the area
## weighting.
## The data may be pre-filtered using CCA.
## Rasmus Benestad, 19.08.2013
#' @export
DS.pca <- function(y,X,verbose=FALSE,plot=FALSE,biascorrect=FALSE,method="lm",swsm=NULL,m=5,ip=1:10,
                   rmtrend=TRUE,weighted=TRUE,pca=TRUE, npca=20,...) {

    if (verbose) {
      print('--- DS.pca ---')
      print(summary(coredata(y)))
      print(class(y))
      print(class(X))
    }
    
    if (class(index(y)) != (class(index(X)))) {
      if (verbose) {print('different class'); summary(coredata(y))}
      warning(paste('DS.pca: different indices:', class(index(y)),class(index(X))))
      if (is.numeric(index(y)) | is.numeric(index(X))) {index(y) <- year(y); index(X) <- year(X)}
      if (verbose) {print('Summary of predictand - intermediate inspection1'); print(zoo(y))}
    }
  
    ## If the predictor is a list, then use DS.list
    if (is.list(X)) {
      if (verbose) print('Predictors represented by a list object')
      z <- DS.list(y,X,biascorrect=biascorrect,
                   method=method,swsm=swsm,m=m,
                   rmtrend=rmtrend,ip=ip,
                   verbose=verbose,weighted=weighted,pca=pca,npca=npca,...)
      return(z)
    }

    ## Check the predictor
    if (!inherits(X,"eof") & inherits(X,"field")) {
      X <- EOF(X)
    } else if (!inherits(X,"eof") & !inherits(X,"field") & inherits(X,"zoo")) {
      ## If the predictor is an index, then use the same code to
      ## estimate teleconnection patterns
      if (verbose) print('Predictor is a zoo object')
      attr(X,'history') <- history.stamp()
      attr(X,'pattern') <- attr(y,'pattern')
      attr(X,'eigenvalues') <- attr(y,'eigenvalues')
      attr(X,'longitude') <- lon(y)
      attr(X,'latitude') <- lat(y)
      attr(X,'variable') <- varid(y)
      attr(X,'unit') <- unit(y)
  
      class(X) <- c('eof',class(X))
      z <- DS.pca(y,X,method=method,swsm=swsm,m=m,
                  ip=ip,rmtrend=rmtrend,verbose=verbose,
                  weighted=weighted,pca=pca,npca=npca,...)
      return(z)
    } else if (verbose) print('Predictor is OK - an EOF object')
   
    ## Check the predictand
    if (inherits(y,"eof") & inherits(y,"field")) {
      if (verbose) print('Make the predictand EOF look like PCAs before downscaling')
      cls0 <- class(y)
      class(y)[1:2] <- c('pca','station')
      z <- DS.pca(y,X,method=method,swsm=swsm,m=m,
                  ip=ip,rmtrend=rmtrend,verbose=verbose,
                  weighted=weighted,...)
      class(z)[2:3] <- c('eof','field')
      attr(z,'pattern') <- attr(y,'pattern')
      attr(z,'eigenvalues') <- attr(y,'eigenvalues')
      attr(z,'longitude') <- lon(y)
      attr(z,'latitude') <- lat(y)
      attr(z,'variable') <- varid(y)
      attr(z,'unit') <- unit(y)      
      return(z)
    }
    
    stopifnot(!missing(y),!missing(X),
              inherits(X,"eof"),inherits(y,"station"))

    cls <- class(y)
    y0 <- y; X0 <- X
                                        #nattr <- softattr(y)

    # synchronise the two zoo objects through 'merge' (zoo)
    if (verbose) { print('Summary of predictand before matchdate'); print(summary(coredata(y)))
      print(index(y)); print(index(X))
    }
    if (verbose) print('predictand y: match date with predictor X')
    y <- matchdate(y,it=X,verbose=verbose) # REB: 2014-12-16
    if (verbose) {print('summary of predictand y after matchdate'); print(summary(coredata(y)))}
    
    if (verbose) print('predictor: match date with predictand')
    X <- matchdate(X,it=y,verbose=verbose) # REB: 2014-12-16
    dy <- dim(y); if (is.null(dy)) dy <- c(length(y),1)
    dx <- dim(X); if (is.null(dx)) dx <- c(length(X),1)
    
    # Use method for downscaling
    #str(y); str(X)
    if (verbose) print(method)
    if (toupper(method)=='mvr') {
        if (verbose) print('MVR')
        if (is.null(dy)) dy <- c(length(y),1)
        if (dy[2]>1) colnames(y) <- paste("y",1:dy[2],sep=".")
        if (is.null(dx)) dx <- c(length(X),1)
        if (dx[2]>1) colnames(X) <- paste("X",1:dx[2],sep=".") 
                                        #colnames(y) <- paste("y",1:dim(y)[2],sep=".")
                                        #colnames(X) <- paste("X",1:dim(y)[2],sep=".")
                                        #    yX <- merge(y,X,all=FALSE)
                                        #
                                        #    vars <- names(yX)
                                        #    ys <- vars[grep('y',vars)]
                                        #    Xs <- vars[grep('X',vars)]
                                        #    ix <- is.element(vars,Xs)
                                        #    iy <- is.element(vars,ys)
                                        #    x <- zoo(coredata(yX[,ix]),order.by=index(yX))
                                        #    y <- zoo(coredata(yX[,iy]),order.by=index(yX))
                                        #    class(y) <- cls
        
        print('WARNING: does not operate on combined objects')
        model <- eval(parse(text=paste(method,"(y,x)",sep="")))
                                        # Reconstruct the entire data, and then apply the PCA to this
        V <- predict(model); W <- attr(y0,'eigenvalues')
        U <- attr(y0,'pattern')
        dU <- dim(U)
        if (length(dU)==3) {
          if (verbose) print('3D -> 2D')
          dim(U) <- c(dU[1]*dU[2],dU[3])
        }
        UWV <- U %*% diag(W) %*% t(V)
        pca <- svd(UWV)
                                        #str(pca)
        ds <- zoo(pca$v,order.by=index(X))
        model$fitted.values <- zoo(model$fitted.values,order.by=index(X)) # +
          ##    attr(y0,'mean') + offset  # REB 04.12.13: included attr(y0,'mean') in
          ##                              # addition to offset
            model$calibration.data <- X
        class(model$fitted.values) <- class(y0)
        attr(ds,'model') <- model
        attr(ds,'quality') <- var(coredata(model$fitted.values))/var(y,na.rm=TRUE)
        attr(ds,'eigenvalues') <- pca$d
        attr(ds,'sum.eigenv') <- sum(pca$d)
        attr(ds,'pattern') <- pca$u
    } else {
        if (verbose) print('Default')
        if (!verbose) pb <- txtProgressBar(style=3)
        ## treat each PC as a station
        if (verbose) print('Prepare output data')
        ## If common EOFs, then accomodate for the additional predictions
        if (!is.null(attr(X0,'n.apps'))) {
          if (verbose) print(attr(X0,'n.apps'))
          Xp <- attr(X0,'appendix.1')
        } else Xp <- X

        
        y.out <- matrix(rep(NA,dx[1]*dy[2]),dx[1],dy[2])
        fit.val <- y.out
        dxp <- dim(Xp); if (is.null(dxp)) dxp <- c(length(Xp),1)
        yp.out <- matrix(rep(NA,dxp[1]*dy[2]),dxp[1],dy[2])
        if (verbose) print(paste('PC',ip,collapse=' '))
                                        # Loop over the PCs...
        ## REB 2015-03-23
        ## The predictor pattern associated with PCA-predictands: one
        ## pattern for each PC. Combine into one matrix. The predictor pattern
        ## for each station can be recovered by multiplying with the PCA pattern
        if (verbose) print('Predictor pattern')
        x0p <- attr(X0,'pattern')
        dp <- dim(x0p)
        if (is.null(dp)) dp <- c(length(x0p),1,1)  # list combining EOFs
        if (length(dp)==2) dp <- c(dp,1)[c(1,3,2)] # if PCA rather than EOF
        if (verbose) {str(x0p); print(dp); print(dy)}
        predpatt <- rep(NA,dp[1]*dp[2]*dy[2])
        dim(predpatt) <- c(dp[1]*dp[2],dy[2])
        dim(x0p) <- c(dp[1]*dp[2],dp[3])
        ## This works for single predictors, but is more complicated for
        ## multiple predictors.
        if (dp[3] == length(attr(X0,'eigenvalues'))) x0p <- x0p %*% diag(attr(X0,'eigenvalues'))
        model <- list(); eof <- list()
        if (verbose) {print('Summary of predictand'); print(summary(coredata(y)))}
        for (i in 1:dy[2]) {
            if (!verbose) setTxtProgressBar(pb,i/dy[2]) 
            ys <- as.station(zoo(y[,i]),loc=loc(y)[i],param=varid(y)[i],
                             unit=unit(y)[i],lon=lon(y)[i],lat=lat(y)[i],
                             alt=alt(y)[i],cntr=cntr(y)[i],stid=stid(y)[i])
            class(ys) <- c('station',class(y)[-c(1:2)])
            
            if (verbose) {print(class(ys)); print(class(X))}
            z <- DS(ys,X,biascorrect=biascorrect,m=m,swsm=swsm,
                    ip=ip,rmtrend=rmtrend,verbose=verbose,...)
            if (verbose) print('--- return to DS.pca ---')

            model[[i]] <- attr(z,'model')
            eof[[i]] <- X
            attr(z,'mean') <- 0 # can't remember why... REB

            ## Check:
            if (verbose) print(paste(i,'y.out[,i]:',
               length(y.out[is.finite(ys),i]),'=',length(z),'?'))
            
            ## Collect the projections in a matrix:
            y.out[is.finite(ys),i] <- coredata(z)
            fit.val[is.finite(ys),i] <- attr(z,'fitted_values')
            if (!is.null(attr(X0,'n.apps')))
                yp.out[,i] <- attr(z,'appendix.1')
            
            ## Also keep the cross-validation
            if (!is.null(attr(z,'evaluation'))) { ## REB 2015-03-27
              if (verbose) print('extract cross-validation')
              if (i==1) cval <- attr(z,'evaluation') else
                        cval <- merge(cval,attr(z,'evaluation'))
            } else cval <- NULL
            ## REB 2015-03-23
            if (verbose) print('Calculate predictor pattern:')
            ## Only if one type of predictor - case with mixed predictors a
            ## bit more complicated -> return NAs.
            ## KMP 2016-01-13 else if only two coefficients
            if ( (dp[3] >= length(attr(z,'model')$coefficients)-1) &
                 (length(attr(z,'model')$coefficients) > 2) ) {
              predpatt[,i] <- x0p[,1:length(attr(z,'model')$coefficients)-1] %*% attr(z,'model')$coefficients[-1]
            } else if (length(attr(z,'model')$coefficients)==2) {
              predpatt[,i] <- x0p[,1:length(attr(z,'model')$coefficients)-1] * attr(z,'model')$coefficients[-1]
            }
        }

        if (verbose) print(paste('Transform back into a PCA-object of dim',
                                 dim(y.out),collapse=' '))
        ds <- zoo(y.out,order.by=index(X))
        ds <- attrcp(y,ds)
        attr(ds,'eigenvalues') <- attr(y0,'eigenvalues')
        attr(ds,'sum.eigenv') <- attr(y0,'sum.eigenv')
        ## REB 2015-03-23
        dim(predpatt) <- c(dp[1],dp[2],dy[2])
        attr(predpatt,'longitude') <- lon(X0)
        attr(predpatt,'latitude') <- lat(X0)
        attr(predpatt,'variable') <- varid(X0)
        attr(predpatt,'unit') <- unit(X0)
        attr(ds,'predictor.pattern') <- predpatt
        ## REB 2015-03-23
        attr(ds,'pattern') <- attr(y0,'pattern')
        if (!is.null(cval))
          names(cval) <- paste(c('X','Z'),'PCA',sort(rep(1:dy[2],2)),sep='.')
        attr(ds,'evaluation') <- cval
        if (!is.null(attr(X0,'n.apps'))) {
            attr(ds,'n.apps') <- 1
            yp.out <- zoo(yp.out,order.by=index(Xp))
            yp.out -> attr(ds,'appendix.1')
        }
    }
    ## Check the 'eof' attribute
    if ( (!is.list(X0)) & is.list(eof) ) {
      eof <- eof[[1]]
      if (verbose) print('Check suggests that eof is stored as list -> eof')
    }
    
    ## Make sure that predictions have the same index class (time units) as the original data 
    if (class(index(ds)) != class(index(y))) {
      ## For annual data:
      if ( (class(index(ds))=='Date') & (class(index(y0))=='numeric') & inherits(y0,'annual') ) 
        index(ds) <- year(index(ds))
      if ( (class(index(ds))=='numeric') & (class(index(y0))=='Date') & inherits(y0,'annual') ) 
        index(ds) <- as.Date(paste(index(ds),'01-01',sep='-'))
    }
    ##print(class(model)); str(model)
    attr(ds,'calibration_data') <- attr(z,'calibration_data')
    attr(ds,'fitted_values') <- zoo(fit.val,order.by=index(attr(z,'fitted_values')))
    class(attr(ds,'fitted_values')) <- class(y0)
    attr(ds,'model') <- model
    attr(ds,'eof') <- eof
    attr(ds,'original_data') <- y0
    attr(ds,'variable') <- varid(y0)
    attr(ds,'mean') <- attr(y0,'mean') # + offset
    attr(ds,'max.autocor') <- attr(y0,'max.autocor')
    attr(ds,'tot.var') <- attr(y0,'tot.var')
                                        #attr(ds,'dimensions') <- c(du[1],du[2])
    attr(ds,'longitude') <- lon(y0)
    attr(ds,'latitude') <- lat(y0)
                                        #  attr(ds,'source') <- paste(attr(y0,'source'),attr(X,'source'),sep="-")
    attr(ds,'source') <- attr(X,'source')
    attr(ds,'history') <- history.stamp(y)
    #attr(ds,'call') <- match.call()
    class(ds) <- c("ds",cls)

                                        #plot(zoo(y[,1],order.by=year(y)),lwd=3)
                                        #lines(zoo(y.out[,1],order.by=year(X)),col='blue',lwd=2)
                                        #lines(zoo(ds[,1],order.by=year(ds)),col='red',lty=2)
    if (verbose) print('--- return ---')
    if (plot) plot(ds)
    invisible(ds)
    
}

' @export
DS.eof <- function(y,X,verbose=FALSE,plot=FALSE,...,biascorrect=FALSE,
                   method="lm",swsm=NULL,m=5,ip=1:10,rmtrend=TRUE,weighted=TRUE,
                   pca=TRUE,npca=20) {
    if (verbose) { print('--- DS.eof ---'); print(summary(coredata(y)))}
    ds <- DS.pca(y,X,biascorrect=biascorrect,
                 method=method,swsm=swsm,m=m,
                 rmtrend=rmtrend,ip=ip,
                 verbose=verbose,weighted=weighted,pca=pca,npca=npca,...)
    if(verbose) print("---return to DS.eof---")
    ## Make sure that predictions have the same index class (time units) as the original data 
    if (class(index(ds)) != class(index(y))) {
      ## For annual data:
      if ( (class(index(ds))=='Date') & (class(index(y))=='numeric') & inherits(y,'annual') ) 
        index(ds) <- year(index(ds))
      if ( (class(index(ds))=='numeric') & (class(index(y))=='Date') & inherits(y,'annual') ) 
        index(ds) <- as.Date(paste(index(ds),'01-01',sep='-'))
    }
    attr(ds,'original_data') <- y
    class(attr(ds,'original_data')) <- class(y)
    class(attr(ds,'fitted_values')) <- class(y)
    if(verbose) print("---return---")
    if (plot) plot(ds)
    invisible(ds)
}

' @export
DS.list <- function(y,X,verbose=FALSE,plot=FALSE,...,biascorrect=TRUE,
                    method="lm",swsm="step",m=5,rmtrend=TRUE,ip=1:7,weighted=TRUE,pca=FALSE,npca=20) {
  ### This method combines different EOFs into one predictor by making a new
  ### data matrix consisting of the PCs, then weight (w) these according to their
  ### eigenvalues (normalised so that each predictor/EOF type carry similar
  ### weight). Then a SVD is applied to this new set of combined PCs to make
  ### an object that looks like on EOF.
  
  if (verbose) { print('--- DS.list ---'); print(summary(coredata(y)))}
  z <- list()
  for (ieof in 1:length(X)) {
    if (verbose) print(names(X)[ieof])
    z[[ieof]] <- DS(y,X[[ieof]],biascorrect=biascorrect,
            method=method,swsm=swsm,m=m,rmtrend=rmtrend,ip=ip,
            verbose=verbose,weighted=weighted,pca=pca,npca=npca,...)
    
    y <- as.residual(z[[ieof]])
  }
  names(z) <-names(X)
  if (plot) plot(z[[1]])
  invisible(z)
}

# NOT EXPORTED
DS.mixedeof <- function(y,X,plot=FALSE,...,it=NULL,biascorrect=TRUE,
                    method="lm",swsm="step",m=5,
                    rmtrend=TRUE,ip=1:7,
                    weighted=TRUE,pca=FALSE,npca=20,verbose=FALSE) {
              ### This method combines different EOFs into one predictor by making a new
              ### data matrix consisting of the PCs, then weight (w) these according to their
              ### eigenvalues (normalised so that each predictor/EOF type carry similar
              ### weight). Then a SVD is applied to this new set of combined PCs to make
              ### an object that looks like on EOF.
    if (verbose) { print('--- DS.mixedeof ---'); print(summary(coredata(y)))}
    preds <- names(X)
    if (verbose) print(preds)
    np <- length(preds)

## Test: if there is only one predictor, use the method for one predictor
    if (np==1) {
        if (verbose) print('Single predictor')
        predictands <- names(y)
        ds <- list(description='ESD')
        for ( i in 1:length(predictands)) {
          ds1 <- DS(y[[i]],X,biascorrect=biascorrect,
                    method=method,swsm=swsm,
                    rmtrend=rmtrend,ip=ip,
                    weighted=TRUE,pca=FALSE,npca=20,...)
          eval(parse(text=paste('ds$',predictands[i],' <- ds1',sep='')))
        }
        return(ds)
        
    } else if (verbose) print('Several predictors')

## REB 2015-04-09: replace the lines below with
      eof <- as.eof.list(X,verbose=verbose)

    if (verbose) print('DS(y,eof,...)')
    ds <- DS(y,eof,biascorrect=biascorrect,
             method=method,swsm=swsm,m=m,
             rmtrend=rmtrend,ip=ip,
             weighted=TRUE,pca=FALSE,verbose=verbose,...)

    ## Now, we need to reconstruct the spatial maps/patterns. There will be
    ## one pattern for each EOF
    if (verbose) print('synthesise patterns and store in a list')
    zpattern <- list()
    if (class(y)[1] != 'pca') {
      ## When y is a pca-object, then pattern is defined for the predictands
      ## and not the predictor.
      if (verbose) str(attr(ds,'pattern'))
      dp <- length(attr(ds,'pattern'))

      ## udv holds the SVD results applied to the EOFs in the list.
      udv <- attr(eof,'udv')
      id <- attr(eof,'id')
      for (i in 1:np) {
            xp <- attr(ds,'pattern')
            if (is.null(dim(xp))) {
              ## PCA-based predictands: the pattern stores the weights of
              ## the regression coefficients
              if (verbose) print('DS pattern: 1D')
              ##diag(c(xp)) %*% udv$v[1:dp,is.element(id,i)] -> eofweights
              ## Y = beta U^T
              xp %*% t(udv$v[is.element(id,i),1:dp]) -> eofweights
            } else {
              ## EOF-based predictand:
              #
              if (verbose) print('DS pattern: 3D')
              xd <- dim(xp)
              dim(xp) <- c(xd[1]*xd[2],xd[3])
              xd %*% t(udv$v[is.element(id,i),1:dp]) -> eofweights
            }
              
            U <- attr(X[[i]],'pattern'); du <- dim(U)
            dim(U) <- c(du[1]*du[2],du[3])
            if (is.null(dim(eofweights))) z <- U %*% diag(eofweights) else
                                          z <- U %*% t(eofweights)
            dim(z) <- c(du[1],du[2])
            attr(z,'longitude') <- lon(X[[i]])
            attr(z,'latitude') <- lat(X[[i]])
                                        #attr(z,'dimensions') <-  attr(X[[i]],'dimensions')[1,2]
            eval(parse(text=paste('zpattern$',preds[i],' <- z',sep='')))
        }
        if (verbose)print(summary(zpattern))
        attr(ds,'pattern') <- zpattern
    }
    if (plot) plot(ds)
    invisible(ds)
}

#' @export
DS.station.pca <- function(y,X,verbose=FALSE,plot=FALSE,...,it=NULL,method="lm",swsm="step",m=5,
                           rmtrend=TRUE,ip=1:7,weighted=TRUE) {
  ## This function does the same as DS.eof
    if (verbose) { print('--- DS.station.pca ---'); print(summary(coredata(y)))}
    z <- DS.default(y=y,X=X,it=it,method=method,swsm=swsm,m=m,
                    rmtrend=trend,ip=ip,verbose=verbose,weighted=weighted)
    if (plot) plot(z)
    return(z)
}

#' @export
DS.trajectory <- function(y,X,verbose=FALSE,plot=FALSE,...,it=NULL,is=NULL,FUN='count',param=NULL,
                       unit=NULL,longname=NULL,loc=NULL,
                       biascorrect=FALSE, method="lm",swsm="step",m=5,
                       rmtrend=TRUE,ip=1:7,
                       weighted=TRUE,pca=FALSE,npca=20) {
   
  if (verbose) { print('--- DS.trajectory ---'); print(summary(coredata(y)))}
  stopifnot(!missing(y),!missing(X),inherits(y,"trajectory"))

  y <- subset(y,it=it,is=is)
  ys <- trajectory2station(y,...)
  
  cls <- class(X)
  if(is.list(X)) cls <- class(X[[1]])
  if(any("season" %in% cls)) {
    ys <- as.4seasons(ys)
  } else if (any("month" %in% cls)) {
    ys <- as.monthly(ys)
  }
  
  ds <- DS(ys,X,biascorrect=biascorrect,method=method,swsm=swsm,m=m,
           rmtrend=rmtrend,ip=ip,verbose=verbose,weighted=weighted,
           pca=pca,npca=npca)
  if (plot) plot(ds)
  invisible(ds)
}
