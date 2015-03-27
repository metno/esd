### Empirical downscaling using EOFs of monthly values from eof.R
## Predictand is a time series of monthly values from NACD or climate station.
###
### Reference: R.E. Benestad et al. (2002),
###            Empirically downscaled temperature scenarios for Svalbard,
###            doi.10.1006/asle.2002.005, September 18.
###
###            R.E. Benestad (2001),
###            A comparison between two empirical downscaling strategies,
###            Int. J. Climatology, 1645-1668, vol. 21, DOI 10.1002/joc.703
###
### R.E. Benestad, met.no, Oslo, Norway 11.04.2013
### rasmus.benestad@met.no
###------------------------------------------------------------------------

## dat -> y; preds -> X
## Needs to work also for daily, annual and seasonal data
##

sametimescale <- function(y,X,FUN='mean',verbose=FALSE) {
 ### Function to ensure that station y has the same time scales as X
  
    if (verbose) print('sametimescale')
    tsx <- class(X)[length(class(X))-1]
    tsy <- class(y)[length(class(y))-1]
    if (verbose) print(c(tsx,tsy))
    if (tsx==tsy) return(y)

    if (verbose) print('Need to aggregate')
    if (tsx=="day") agrscly <- as.Date(index(y)) else
    if (tsx=="month") agrscly <- as.yearmon(index(y)) else
    if (tsx=="annual") agrscly <- year(y) else
    if (tsx=="year") agrscly <- year(y)
    if (verbose) str(agrscly)
    if (tsx !="season")
        y <- aggregate(y, agrscly, match.fun(FUN)) else 
        y <- as.4seasons(y, FUN=match.fun(FUN),dateindex=TRUE)
    return(y)
}



DS<-function(y,X,verbose=TRUE,...) UseMethod("DS")

                                        # The basic DS-function, used by other methods
                                        # REB HERE!

DS.default <- function(y,X,mon=NULL,
                       method="lm",swsm="step",m=5,
                       rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                       verbose=FALSE,weighted=TRUE,...) {
    ##
    if (verbose) print('DS.default')
    #print('err(y)'); print(err(y))
    if (verbose) {print('index(y)'); print(index(y))}
    if (verbose) {print(class(y)); print(class(X))}
    
    swapped <- FALSE
    if ( inherits(y,c("eof")) & inherits(X,c("station"))) {
        yy <- X
        X <- y
        y <- yy
        swapped <- TRUE
    }
    stopifnot(!missing(y),!missing(X), is.matrix(X),
              inherits(X,"eof"),inherits(y,"station"))
    y0 <- y
    X0 <- X

    if (verbose) {print(paste(sum(!is.finite(coredata(y))),'missing values in y'))}
    if (verbose)  {print('index(y) before removing missing values:'); print(index(y))}
    y <- subset(y,it=is.finite(coredata(y)))
    W <- attr(X,'eigenvalues')
    cls <- class(X)
 
    if (verbose) print('Ensure matching time scale')
    if (verbose) {print('index(y) before sametimescale:'); print(index(y))}
    y <- sametimescale(y,X,verbose=verbose)
    ## REB is needed to ensure that maye y annual if X is annual
    #if (verbose) {print('index(y) after sametimescale:'); print(index(y))}
    #if (verbose) print('Match dates')
    y <- matchdate(y,X,verbose=verbose) ##
    #if (verbose) {print('index(y) after matchdate'); print(index(y))}
    X <- matchdate(X,y,verbose=verbose) # REB 2015-01-14
    if (verbose) {print("index(y) & index(X) after synch:");
                  print(index(y)); print(index(X))}
    
    if (!is.null(mon)) y <- subset(y,it=mon)
    
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

                                        #str(y); print(class(y))

    ## REB: 2014-10-03: add weights if available
    if (!is.null(attr(y,'standard.error')))
        weights <- 1/coredata(attr(y,'standard.error')) else
    weights <- rep(1,length(y))
    weights[!is.finite(weights)] <- 0
    if (is.null(attr(y,'standard.error'))) weighted <- FALSE
    if (verbose) {print(paste('weights',weighted)); print(weights)}

    caldat <- data.frame(y=coredata(y),X=as.matrix(coredata(X)),
                         weights=weights)
    predat <- data.frame(X=as.matrix(coredata(X0)))
    colnames(predat) <- paste("X",1:length(colnames(predat)),sep=".")

    if (is.null(names(X))) names(X) <- 1:dim(X)[2]
    Xnames <- paste("X.",1:length(names(X)),sep="")
    colnames(caldat) <- c("y",Xnames,'weights')
    Xnames <- Xnames[eofs]
    #browser()
                                        # REB 2014-10-03:
    if (weighted)
        calstr <- paste(method,"(y ~ ",paste(Xnames,collapse=" + "),
                        ", weights=weights, data=caldat, ...)",sep="") else
        calstr <- paste(method,"(y ~ ",paste(Xnames,collapse=" + "),
                        ", data=caldat, ...)",sep="")
    MODEL <- eval(parse(text=calstr))
    FSUM <- summary(MODEL)
    if (verbose) print(FSUM)

    ## Stepwise regression
    if (!is.null(swsm)) {
        cline <- paste("model <- ",swsm,"(MODEL,trace=0)",sep="")
        eval(parse(text=cline))
    } else
        model <- MODEL
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
          diag(attr(X,'eigenvalues')[eofs]) %*% t(U[,eofs])
      dim(pattern) <- c(du[1],du[2])
    } else pattern <- c(COEFS[2:dc[1],1]) * attr(X,'eigenvalues')[eofs]
                                                 
    
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
    ds <- attrcp(y0,ds,ignore='names')
    pattern <- attrcp(X0,pattern,ignore=c('longitude','latitude','names'))
    
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
    class(ds) <- c("ds",cls)
    rm("y0","X0")
                                        #print("Completed")
                                        #lines(ds,col="darkred",lwd=2,lty=2)
                                        #lines(attr(ds,'original_data'),col="green",lwd=2,lty=2)
    if (verbose) print('exit DS.default')
    invisible(ds)
}





DS.station <- function(y,X,biascorrect=FALSE,mon=NULL,
                       method="lm",swsm="step",m=5,
                       rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                       verbose=FALSE,weighted=TRUE,pca=FALSE,npca=20,...) {
    ## 
    stopifnot(!missing(y),!missing(X),inherits(y,"station"))
    if (verbose) print("DS.station")
    #print('err(y)'); print(err(y))
    #print('index(y)'); print(index(y))
    
    if ( (!inherits(X,'eof')) & (inherits(X,'field')) ) {
                                        #print("HERE")
        ds <- DS.field(y=y,X=X,biascorrect=biascorrect,mon=mon,
                       method=method,swsm=swsm,m=m,
                       rmtrend=rmtrend,eofs=eofs,
                       area.mean.expl=area.mean.expl,verbose=verbose,
                       weighted=TRUE,pca=FALSE,npca=20,...) 
        return(ds)
    } else if (is.list(X)) {
                                        # REB 2014-10-08
        print("The predictor is a list")
        ds <- DS.list(y=y,X=X,biascorrect=biascorrect,mon=mon,
                      method=method,swsm=swsm,m=m,
                      rmtrend=rmtrend,eofs=eofs,
                      area.mean.expl=area.mean.expl,verbose=verbose,
                      weighted=TRUE,pca=FALSE,npca=20,...) 
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
    } else Y <- y
    

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
                ds <- DS.comb(y=z,X=X,biascorrect=biascorrect,mon=mon,
                              method=method,swsm=swsm,
                              rmtrend=rmtrend,eofs=eofs,
                              area.mean.expl=area.mean.expl,verbose=verbose,...)
                if (verbose) print("---")
            } else if (inherits(X,'eof')) {
                if (verbose) print("*** EOF ***")
                ## X is ordinary EOF
                ds <- DS.default(y=z,X=X,mon=mon,
                                 method=method,swsm=swsm,
                                 rmtrend=rmtrend,eofs=eofs,
                                 area.mean.expl=area.mean.expl,
                                 verbose=verbose,...)
                                        #if (verbose) print("+++")
            }
        } else if (inherits(X,'field')) {
            if (verbose) print("the predictor is a field-object")
            ## X is a field
            ds <- DS.field(y=z,X=X,biascorrect=biascorrect,mon=mon,
                           method=method,swsm=swsm,
                           rmtrend=rmtrend,eofs=eofs,
                           area.mean.expl=area.mean.expl,verbose=verbose,...)
        }
        ## May need an option for coombined field: x is 'field' + 'comb'
        if (verbose) print("Cross-validation")

        ## Unless told not to - carry out a cross-validation
        if (!is.null(m))  {
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
    } else ds.results <- dsall
    
    if (verbose) print("--- exit DS.station ---")
    invisible(ds.results)  
}


## DS for combined fields - to make predictions not based on the
## calibration data

DS.comb <- function(y,X,biascorrect=FALSE,mon=NULL,
                    method="lm",swsm="step",m=5,
                    rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                    verbose=FALSE,weighted=TRUE,...) {

    if (verbose) print("DS.comb")
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

    
    if (!inherits(X,"eof")) X <- EOF(X,mon=mon,area.mean.expl=area.mean.expl)
    
    
    if (biascorrect) {
        if (verbose) print("Bias correcion - bias-fix common EOF")
        X <- biasfix(X)
    }
    
    ds <- DS.default(y,X,mon=mon,method=method,swsm=swsm,m=m,
                     rmtrend=rmtrend,eofs=eofs,
                     area.mean.expl=area.mean.expl,verbose=verbose,...)

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
        colnames(newdata) <- paste("X",1:length(colnames(X)),sep=".")
                                        #print(summary(newdata))
        z <- predict(model,newdata=newdata) + attr(ds,'mean')
        Y <- zoo(z,order.by=index(Z))
        Y <- attrcp(Z,Y)

                                        # Store in a similar structure as in the common EOF:
        eval(parse(text=paste("attr(ds,'appendix.",i,"') <- Y",sep="")))
    }
    rm("X0"); gc(reset=TRUE)
    invisible(ds)
}


## This function takes care of downscaling based on a field-object X.
## X can be a combined field. This function calls more primitive DS methods,
## depending on the time scale represented in X (monthly or seasonal).

DS.field <- function(X,y,biascorrect=FALSE,mon=NULL,
                     method="lm",swsm="step",m=5,
                     rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                     verbose=FALSE,weighted=TRUE,...) {
    if (verbose) print("DS.field")
    ## Keep track of which is an eof object and which is a station record:
    swapped <- FALSE
    if ( inherits(y,c("field")) & inherits(X,c("station"))) {
        yy <- X
        X <- y
        y <- yy
        swapped <- TRUE
    }
    
                                        #print(class(X)); print(attr(y,'variable'))
    if (verbose) {print(class(X)); print(varid(y))}
    if (sum(is.element(tolower(attr(y,'variable')),c('t2m','tmax','tmin'))) >0) {
        if (inherits(X,'month')) 
            ds <- DS.t2m.month.field(y=y,X=X,biascorrect=biascorrect,
                                     mon=mon,method=method,swsm=swsm,m=m,
                                     rmtrend=rmtrend,eofs=eofs,
                                     area.mean.expl=area.mean.expl,
                                     verbose=verbose) else
        if (inherits(X,'season')) 
            ds <- DS.t2m.season.field(y=y,X=X,biascorrect=biascorrect,
                                      method=method,swsm=swsm,m=m,
                                      rmtrend=rmtrend,eofs=eofs,
                                      area.mean.expl=area.mean.expl,
                                      verbose=verbose) else
        if (inherits(X,'annual')) 
            ds <- DS.t2m.annual.field(y=y,X=X,biascorrect=biascorrect,
                                      method=method,swsm=swsm,m=m,
                                      rmtrend=rmtrend,eofs=eofs,
                                      area.mean.expl=area.mean.expl,
                                      verbose=verbose)
    } else if (tolower(attr(y,'variable'))=='precip')
        ds <- DS.precip.season.field(y=y,X=X,biascorrect=biascorrect,
                                     method=method,swsm=swsm,m=m,
                                     rmtrend=rmtrend,eofs=eofs,
                                     area.mean.expl=area.mean.expl,verbose=verbose)
    else ds <- DS.default(y=y,X=X,biascorrect=biascorrect,
                          method=method,swsm=swsm,m=m,
                          rmtrend=rmtrend,eofs=eofs,
                          area.mean.expl=area.mean.expl,verbose=verbose)
    if (verbose) print('return downscaled results')
    invisible(ds)
}



## Downscales all 12-calendar months for a field by stepping through the months
## and compute the EOFs before applying the default DS method.
DS.t2m.month.field <- function(y,X,biascorrect=FALSE,mon=NULL,
                               method="lm",swsm="step",m=m,
                               rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                               verbose=FALSE,weighted=TRUE,station=TRUE) {
    if (verbose) print("DS.t2m.month.field")
    if (inherits(X,'comb')) type <- 'eof.comb' else type <- "eof.field"
    cls <- class(y)

    ds <- list()
    if (is.null(mon)) mon <-  1:12
    
    for (i in mon) {
                                        #print(month.name[i])
        eof <- EOF(X,it=i,area.mean.expl=area.mean.expl)
        if (biascorrect) eof <- biasfix(eof)
        cline <- paste("ds$",month.abb[i],
                       "<- DS.station(y,eof,method=method,swsm=swsm,m=m,",
                       "rmtrend=rmtrend,eofs=eofs,",
                       "area.mean.expl=area.mean.expl,verbose=verbose)",
                       sep="")
        if (verbose) print(cline)
        eval(parse(text=cline))
    }
                                        #print(summary(ds))
    
    if (station) ds <- combine.ds(ds) else
    cls <- c("list","dsfield",cls)

                                        #ds <- attrcp(y,ds)
                                        #nattr <- softattr(y)
                                        #for (i in 1:length(nattr))
                                        #  attr(ds,nattr[i]) <- attr(y,nattr[i])

    attr(ds,'predictand.class') <- cls
                                        #class(ds) <- c("list","dsfield",cls)
    invisible(ds)
}


DS.t2m.season.field <- function(y,X,biascorrect=FALSE,
                                method="lm",swsm="step",m=5,
                                rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                                verbose=FALSE,weighted=TRUE,station=TRUE) {
  ## Downscale seasonal mean and standard deviation
    if (verbose) print("DS.t2m.season.field")

    Z1 <- EOF(subset(X,it='djf'),area.mean.expl=area.mean.expl)
    if (verbose) print("downscale DJF")
    ds1 <- DS(y,Z1,biascorrect=biascorrect,eofs=eofs)
    Z2 <- EOF(subset(X,it='mam'),area.mean.expl=area.mean.expl)
    if (verbose) print("downscale MAM")
    ds2 <- DS(y,Z2,biascorrect=biascorrect,eofs=eofs)
    Z3 <- EOF(subset(X,it='jja'),area.mean.expl=area.mean.expl)
    if (verbose) print("downscale JJA")
    ds3 <- DS(y,Z3,biascorrect=biascorrect,eofs=eofs)
    
    
    Z4 <- EOF(subset(X,it='son'),area.mean.expl=area.mean.expl)
    if (verbose) print("downscale SON")
    ds4 <- DS(y,Z4,biascorrect=biascorrect,eofs=eofs)
    if (verbose) print("Combine the 4 seasons")
    ds <- combine(list(ds1,ds2,ds3,ds4))
    z <- c(crossval(ds1,m=m),crossval(ds2,m=m),
           crossval(ds3,m=m),crossval(ds4,m=m))
    attr(ds,'evaluation') <- z
    save(file='inside.ds.seas.f.1.rda',y,Z1,ds1,Z2,ds2,Z3,ds3,Z4,ds4,X)
                                        #
    invisible(ds)
}

DS.t2m.annual.field <- function(y,X,biascorrect=FALSE,
                                method="lm",swsm="step",m=5,
                                rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                                verbose=FALSE,weighted=TRUE,station=TRUE) {
  ## Downscale seasonal mean and standard deviation
    if (verbose) print("DS.t2m.annual.field")

    
    Z <- EOF(annual(X),area.mean.expl=area.mean.expl)
    ds <- DS(annual(y),Z,biascorrect=biascorrect)
    invisible(ds)
}



DS.precip.season.field <- function(y,X,biascorrect=FALSE,threshold=1,
                                   method="lm",swsm="step",m=5,
                                   rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                                   verbose=FALSE,weighted=TRUE,...) {

  ## Computes the annual mean values for wet-day mean mu, wet-day frequency, and spell.
  ## Also computes seasonal variations from PCA X[year,calendar months].
  ## One PC for each year.

    mu <- as.4seasons(y,FUN="exceedance",threshold=threshold)
    fw <- as.4seasons(y,FUN="exceedance",fun="freq")
    wL <- as.4seasons(spell(y))

    if (!inhetits(X,'season')) X <- as.4seasons(X)

    for (i in 1:4) {
        x <- EOF(X,it=i,area.mean.expl=area.mean.expl)
        if (biascorrect) x <- biasfix(x)
        ds.mu <- DS.default(mu,x,method=method,swsm=swsm,m=m,
                            rmtrend=rmtrend,eofs=eofs,
                            verbose=verbose,...)
        ds.fw <- DS.freq(fw,x,rmtrend=rmtrend,eofs=eofs,m=m,
                         verbose=verbose,...)
        ds.wL <- DS.spell(wL,x,rmtrend=rmtrend,eofs=eofs,m=m,
                          verbose=verbose,...)
    }
    
    ## DS mu, fw, spell, total
    ds <- NULL
    invisible(ds)
}

## Use family='gaussian' for sample sizes gt 30 - > central limit theorem
DS.freq <- function(y,X,threshold=1,biascorrect=FALSE,method="glm",
                    family="gaussian",swsm="step",m=5,
                    rmtrend=TRUE,eofs=1:7,verbose=FALSE,weighted=TRUE,...) {
    if (inherits(X,'month'))
        Z <- aggregate(y,as.yearmon,FUN="wetfreq",threshold=threshold) else
    if (inherits(X,'season'))
        Z <- as.4seasons(y,FUN=wetfreq,threshold=threshold) else
    if (inherits(X,'annual'))
        Z <- annual(y,FUN=wetfreq,threshold=threshold)
    
    ds <- DS.default(Z,X,biascorrect=biascorrect,method=method,
                     swsm=swsm,m=m,rmtrend=rmtrend,eofs=eofs,verbose=verbose,...)
    return(ds)
}


DS.spell <- function(y,X,threshold=1,biascorrect=FALSE,
                     method="glm",family="gaussian",swsm="step",m=5,
                     rmtrend=TRUE,eofs=1:7,verbose=FALSE,weighted=TRUE,...) {
  ## Downscale the spell length using a GLM with poisson family.
  ##  the mean spell length over a given interval:
    if (inherits(y,'spell')) z <- as.station(y) else
    if (inherits(y,'sstation')) {
        z <- as.station(spell(y))
    }
    
    if (inherits(X,'month')) Z <- aggregate(z,as.yearmon,FUN=mean) else
    if (inherits(X,'season')) Z <- as.4seasons(z,FUN=mean) else
    if (inherits(X,'annual')) Z <- annual(z,FUN=mean)

    ds <- DS(Z,X,biascorrect=biascorrect,method=method,swsm=swsm,m=m,
             rmtrend=rmtrend,eofs=eofs,
             verbose=verbose,...)
    return(ds)
}


## DS.pca
## This function applies to PCA results for a group of stations.
## Typically some annual statistics (e.g. mu), stored in compressed form
## as a PCA-object (similar to EOF, but not gridded and without the area
## weighting.
## The data may be pre-filtered using CCA.
## Rasmus Benestad, 19.08.2013
DS.pca <- function(y,X,biascorrect=FALSE,mon=NULL,
                   method="lm",swsm=NULL,m=5,eofs=1:10,
                   rmtrend=TRUE,verbose=FALSE,weighted=TRUE,...) {
    
    if (verbose) {print('DS.pca'); print(class(X))}
    if (is.list(X)) {
      if (verbose) print('Predictors represented by a list object')
      z <- DS.list(y,X,biascorrect=biascorrect,mon=mon,
                   method=method,swsm=swsm,m=m,
                   rmtrend=rmtrend,eofs=eofs,area.mean.expl=area.mean.expl,
                   verbose=verbose,weighted=weighted,pca=pca,npca=npca,...)
      return(z)
    }

    cls <- class(y)
    y0 <- y; X0 <- X
                                        #nattr <- softattr(y)

                                        # synchronise the two zoo objects through 'merge' (zoo)
    y <- matchdate(y,it=X,verbose=verbose) # REB: 2014-12-16
    X <- matchdate(X,it=y,verbose=verbose) # REB: 2014-12-16
    dy <- dim(y)
    dx <- dim(X)

                                        # Use method for downscaling
                                        #str(y); str(X)
    if (verbose) print(method)
    if (toupper(method)=='mvr') {
        if (verbose) print('MVR')
        if (is.null(dy)) dy <- c(length(y),1)
        if (dy[2]>1) colnames(y) <- paste("y",1:dy[2],sep=".")
        if (is.null(dx)) dx <- c(length(x),1)
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
        UWV <- U %*% diag(W) %*% t(V)
        pca <- svd(UWV)
                                        #str(pca)
        ds <- zoo(pca$v,order.by=index(X))
        model$fitted.values <- zoo(model$fitted.values,order.by=index(X)) +
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
        yp.out <- matrix(rep(NA,dim(Xp)[1]*dy[2]),dim(Xp)[1],dy[2])

        if (verbose) print(paste('PC',eofs,collapse=' '))
                                        # Loop over the PCs...
        ## REB 2015-03-23
        ## The predictor pattern associated with PCA-predictands: one
        ## pattern for each PC. Combine into one matrix. The predictor pattern
        ## for each station can be recovered by multiplying with the PCA pattern
        if (verbose) print('Predictor pattern')
        x0p <- attr(X0,'pattern') %*% attr(X0,'eigenvalues')
        dp <- dim(x0p)
        if (is.null(dp)) dp <- c(length(x0p),1,1)  # list combining EOFs
        #str(x0p); print(dp); print(dy)
        predpatt <- rep(NA,dp[1]*dp[2]*dy[2])
        dim(predpatt) <- c(dp[1]*dp[2],dy[2])
        dim(x0p) <- c(dp[1]*dp[2],dp[3])
        model <- list(); eof <- list()
        for (i in 1:dy[2]) {
            if (!verbose) setTxtProgressBar(pb,i/dy[2]) 
            ys <- as.station(zoo(y[,i]),loc=loc(y)[i],param=varid(y)[i],
                             unit=unit(y)[i],lon=lon(y)[i],lat=lat(y)[i],
                             alt=alt(y)[i],cntr=cntr(y)[i],stid=stid(y)[i])
            class(ys) <- c('station',class(y)[-c(1:2)])
            
            if (verbose) {print(class(ys)); print(class(X))}
            z <- DS(ys,X,biascorrect=biascorrect,
                    eofs=eofs,rmtrend=rmtrend,verbose=verbose,...)
            model[[i]] <- attr(z,'model')
            eof[[i]] <- X
            if (verbose) print('--- return to DS.pca ---')
            attr(z,'mean') <- 0 # can't remember why... REB

            ## Check:
            if (verbose) print(paste(i,'y.out[,i]:',
                                     length(y.out[is.finite(ys),i]),'=',length(z),'?'))
            
            ## Collect the projections in a matrix:
            y.out[is.finite(ys),i] <- coredata(z)
            fit.val[is.finite(ys),i] <- attr(z,'fitted_values')
            if (!is.null(attr(X0,'n.apps')))
                yp.out[is.finite(ys),i] <- attr(z,'appendix.1')
            
            ## Also keep the cross-validation
            if (!is.null(attr(z,'evaluation'))) { ## REB 2015-03-27
              if (i==1) cval <- attr(z,'evaluation') else
                        cval <- merge(cval,attr(z,'evaluation'))
            } else cval <- NULL
            ## REB 2015-03-23
            if (verbose) print('Calculate predictor pattern:')
            ## Only if one type of predictor - case with mixed predictors a
            ## bit more complicated -> return NAs.
            if ( (dp[3] >= length(attr(z,'model')$coefficients)-1) &
                 (length(attr(z,'model')$coefficients) > 2) )
              predpatt[,i] <- x0p[,1:length(attr(z,'model')$coefficients)-1] %*%
                attr(z,'model')$coefficients[-1]
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
                                        #print(class(model)); str(model)
    attr(ds,'calibration_data') <- attr(z,'calibration_data')
    attr(ds,'fitted_values') <- zoo(fit.val,
                                    order.by=index(attr(z,'fitted_values')))
    class(attr(ds,'fitted_values')) <- class(y0)
    attr(ds,'model') <- model
    attr(ds,'eof') <- eof
    attr(ds,'original_data') <- y
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
    attr(ds,'call') <- match.call()
    class(ds) <- c("ds",cls)

                                        #plot(zoo(y[,1],order.by=year(y)),lwd=3)
                                        #lines(zoo(y.out[,1],order.by=year(X)),col='blue',lwd=2)
                                        #lines(zoo(ds[,1],order.by=year(ds)),col='red',lty=2)

    invisible(ds)
}

DS.eof <- function(y,X,mon=NULL,
                   method="lm",swsm="step",m=5,
                   rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                   verbose=FALSE,weighted=TRUE,pca=TRUE,...) {
    if (verbose) print("DS.eof")
    ds <- DS.pca(y,X,mon=mon,
                 method=method,swsm=swsm,m=m,
                 rmtrend=rmtrend,eofs=eofs,
                 area.mean.expl=area.mean.expl,
                 verbose=verbose,...)
    return(ds)
}


DS.list <- function(y,X,biascorrect=TRUE,mon=NULL,
                    method="lm",swsm="step",m=5,
                    rmtrend=TRUE,eofs=1:7,
                    verbose=FALSE,weighted=TRUE,pca=FALSE,npca=20,...) {
              ### This method combines different EOFs into one predictor by making a new
              ### data matrix consisting of the PCs, then weight (w) these according to their
              ### eigenvalues (normalised so that each predictor/EOF type carry similar
              ### weight). Then a SVD is applied to this new set of combined PCs to make
              ### an object that looks like on EOF.
    if (verbose) print('DS.list')
    preds <- names(X)
    if (verbose) print(preds)
    np <- length(preds)

## Test: if there is only one predictor, use the method for one predictor
    if (np==1) {
        if (verbose) print('Single predictor')
        predictands <- names(y)
        ds <- list(description='ESD')
        for ( i in 1:length(predictands)) {
          ds1 <- DS(y[[i]],X,biascorrect=biascorrect,mon=mon,
                    method=method,swsm=swsm,
                    rmtrend=rmtrend,eofs=eofs,
                    weighted=TRUE,pca=FALSE,npca=20,...)
          eval(parse(text=paste('ds$',predictands[i],' <- ds1',sep='')))
        }
        return(ds)
        
    } else if (verbose) print('Several predictors')

## Combine the different predictors into one matrix: also for comb...
    x <- zoo(X[[1]])
    w <- attr(X[[1]],'eigenvalues')/sum(attr(X[[1]],'eigenvalues'))
    id <- rep(1,length(attr(X[[1]],'eigenvalues')))

## If combined EOF - need to get the appended fields too
    if (inherits(X[[1]],'comb')) {
        n.app <- attr(X[[1]],'n.apps')
        if (n.app > 1) print('This only works with n.app==1')
        print("combined field")
        z <- attr(X[[1]],'appendix.1')
    }
    for (i in 2:np) {
        x <- merge(x,zoo(X[[i]]),all=TRUE)
        w <- c(w,attr(X[[i]],'eigenvalues')/sum(attr(X[[i]],'eigenvalues')))
        id <- c(id,rep(i,length(attr(X[[i]],'eigenvalues'))))
        if (inherits(X[[1]],'comb')) {
            if ( (!inherits(X[[i]],'comb')) | (attr(X[[i]],'n.apps') != n.app) )
                stop('DS.list: the predictors in the list must match')
            z <- merge(z,attr(X[[1]],'appendix.1'))
        }
    }

    if (verbose) print(c(dim(x),length(w)))
    t <- index(x)
    ## apply the weights
    x <- x %*% diag(w)
    xm <- rowMeans(x)
    x <- x[is.finite(xm),]; t <- t[is.finite(xm)]

    ## Apply an SVD to the combined PCs to extract the common signal in the
    ## different predictors - these are more likely to contain real physics
    ## and be related to the predictand.
    if (verbose) print('svd')
    udv <- svd(coredata(x))
    if (verbose) print(summary(udv))

    ## If the predictor is a common EOF, then also combine the appended fields
    ## the same way as the original flield.
    if (inherits(X[[1]],'comb')) {
        z <- z %*% diag(w)
        udvz <- svd(coredata(z))
    }

    eof <- zoo(udv$u[,1:20],order.by=t)

    ## Let the pattern contain the weights for the EOFs in the combined
    ## PC matrix, rather than spatial patterns. The spatial patterns are
    ## then reconstructed from these.
    pattern <- matrix(rep(1,length(udv$v[,1:20])),dim(udv$v[,1:20]))

    ## Do a little 'fake': here the pattern is not a geographical map but weight
    ## for the EOFs.
    dim(pattern) <- c(1,dim(pattern))
    if (verbose) str(pattern)
    attr(eof,'eigenvalues') <- udv$d
    attr(eof,'pattern') <- rep(1,dim(udv$v)[1])
    names(eof) <- paste("X.",1:20,sep="")
    
    class(eof) <- class(X[[1]])
    if (inherits(X[[1]],'comb'))
        attr(eof,'appendix.1') <- z

    #browser()
    if (verbose) print('DS(y,eof,...)')
    ds <- DS(y,eof,biascorrect=biascorrect,
             method=method,swsm=swsm,m=m,
             rmtrend=rmtrend,eofs=eofs,
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
              #browser()
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
    return(ds)
}

DS.station.pca <- function(y,X,mon=NULL,
                           method="lm",swsm="step",m=5,
                           rmtrend=TRUE,eofs=1:7,area.mean.expl=FALSE,
                           verbose=FALSE,weighted=TRUE,...) {
  ## This function does the same as DS.eof
    z <- DS.default(y=y,X=X,mon=mon,method=method,swsm=swsm,m=m,
                    rmtrend=trend,eofs=eofs,area.mean.expl=area.mean.expl,
                    verbose=verbose,weighted=weighted,..)
    return(z)
}

biasfix <- function(x) {
    stopifnot(!missing(x), inherits(x,"eof"),inherits(x,"comb"))
    ## Check if the results already have been bias-corrected
    if (!is.null(attr(x,'diagnose'))) return(x)
    diag <- diagnose(x)
    n <- attr(x,'n.apps')
    for ( i in 1:n ) {
        eval(parse(text=paste("z <- attr(x,'appendix.",i,"')",sep="")))
        Z <- coredata(z)

        ## diagnose: (1 + sd(z))/(1 + sd(x))
        ## x is reanalysis; z is gcm:
        for (j in 1:length(Z[1,]))
            Z[,j] <- Z[,j]/diag$sd.ratio[i,j] + diag$mean.diff[i,j]
        y <- zoo(Z,order.by=index(z))
        y <- attrcp(z,y)
        eval(parse(text=paste("y -> attr(x,'appendix.",i,"')",sep="")))
    }
    attr(x,'history') <- history.stamp(x)
    attr(x,'quality') <- "'bias' corrected -  ref (Imbert & Benestad (2005); Theor. Appl. Clim.; DOI: 10.1007/s00704-005-0133-4"
    attr(x,'diagnose') <- diag
    return(x)
}
