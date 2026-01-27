## DS function for downscaling by calibrating the model on the
## similarity of the seasonal cycle in the predictor and predictand.
## The inputs y & X are expected to be station and field objects

DS.seasonalcycle <- function(y,X,...) {
#DS.seasonalcycle <- function(y,X,mon=NULL,FUN='wetmean',fun='mean',
#                             method="lm",swsm="step",m=1,
#                              ip=1:7,area.mean.expl=FALSE,
#                              verbose=FALSE,weighted=TRUE,...) {
    
  
  ## REB 2026-01-22
  args <- list(...)
  if (!is.null(args$method)) method <- args$method else method <-"lm"
  if (!is.null(args$swsm)) swsm <- args$swsm else swsm <-"step"
  if (!is.null(args$m)) m <- args$m else m <-5
  if (!is.null(args$ip)) ip <- args$ip else ip <-1:7
  if (!is.null(args$rmtrend)) rmtrend <- args$rmtrend else rmtrend <-TRUE
  if (!is.null(args$weightd)) weighted <- args$weighted else weighted <-TRUE
  if (!is.null(args$verbose)) verbose <- args$verbose else verbose <-FALSE
  if (!is.null(args$pca)) pca <- args$pca else pca <-FALSE
  if (!is.null(args$biascorrect)) biascorrect <- args$biascorrect else biascorrect <-FALSE
  if (!is.null(args$npca)) npca <- args$npca else pca <-20
  if (!is.null(args$plot)) plot <- args$plot else plot <-FALSE
  if (!is.null(args$FUN)) FUN <- args$FUN else FUN <- 'wetmean'
  if (!is.null(args$fun)) fun <- args$fun else FUN <- 'mean'
  if (!is.null(args$area.mean.expl)) area.mean.expl <- args$area.mean.expl else area.mean.expl <-FALSE
  mon <- args$mon
  
  if (verbose) print('DS.seasonalcycle')
  stopifnot(inherits(y,'station'),!inherits(y,'annual'),
            inherits(X,'field'))

  ## Keep a copy of the original data
  y0 <- y; X0 <- X
  ## Find the seaonal cycle in predictand:
  if (verbose) print('seasonal cycle in predictand')
  y <- aggregate(y,by=month,FUN=FUN)
  yam <- annual(y0,FUN=FUN)
  ## Estimate weights
  if (weighted) weights <- aggregate(y,by=month,FUN='nv') else
                weights <- rep(1,12)
  ## Find the seaonal cycle in predictor:
  if (verbose) print('seasonal cycle in predictor')

  if (!inherits(X,'eof')) {
    X <- aggregate(X0,by=month,FUN=fun)
    ## Also estimate the annual mean values
    #browser()
    Xam <- annual(X0,FUN=fun)
    if (!area.mean.expl) {
  ## Estimate the EOFs for the seasonal variations
      if (verbose) print('Estimate EOF of annual cycle')
      if (verbose) print('combine seasonal cycle and annual')
      XX <- combine(X,Xam)
      eof <- EOF(XX,n=12,verbose=verbose)
      caldat <- data.frame(y=coredata(y),X=coredata(eof),
                           weights=coredata(weights))
      predat <- data.frame(X=attr(eof,'appendix.1'))
    } else {
  ## Use the area mean as predictor
      if (verbose) print('Use area average')
      caldat <- data.frame(y=coredata(y),X=coredata(aggregate.area(X)),
                           weights=coredata(weights))
      predat <- data.frame(X=coredata(aggregate.area(Xam)))
    }
  } else {
      ## If the predictor is pre-computed
    if (verbose) print('Use pre-computed predictor')
      caldat <- data.frame(y=coredata(y),X=coredata(X),
                           weights=coredata(weights))
      Xam <- attr(X,'appendix.1')
      predat <- data.frame(X=coredata(Xam))
  }

  ## Organising data by the names of the variable.
  if (verbose) print('Organise the data')
  colnames(predat) <- paste("X",1:length(colnames(predat)),sep=".")
  Xnames <- paste("X.",1:dim(X)[1],sep="")
  if (verbose) print(Xnames)
  colnames(caldat) <- c("y",Xnames,'weights')
  colnames(predat) <- Xnames
  Xnames <- Xnames[ip]
  if (weighted)
    calstr <- paste(method,"(y ~ ",paste(Xnames,collapse=" + "),
                    ", weights=weights, data=caldat, ...)",sep="") else
    calstr <- paste(method,"(y ~ ",paste(Xnames,collapse=" + "),
                    ", data=caldat, ...)",sep="")
  ## Calibrate the model
  MODEL <- eval(parse(text=calstr))
  FSUM <- summary(MODEL)
  if (verbose) print(FSUM)

   ## Stepwise regression
  if (!is.null(swsm)) {
    cline <- paste("model <- ",swsm,"(MODEL,trace=0)",sep="")
    eval(parse(text=cline))
  } else model <- MODEL
  ## House keeping
  terms1 <- attr(model$terms,'term.labels')
  if (verbose) print(summary(model))
  fsum <- summary(model)
  ## Set coefficients for variables excluded by stepwise screening to 0,
  ## but keep the original set of coefficients as before the screening
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
  ## Find the spatial patterns implied by the regression coefficients
  ## and EOF patterns
  U <- attr(eof,'pattern'); du <- dim(U)
  if (verbose) {print('pattern dimension'); print(du); str(U)}
  if (length(du)==3) dim(U) <- c(du[1]*du[2],du[3])
  pattern <- NULL
  #browser()
  if ( (!is.null(du)) & (!area.mean.expl) ) {
    pattern <- t(COEFS[2:dc[1],1]) %*%
      diag(attr(eof,'eigenvalues')[ip]) %*% t(U[,ip])
    dim(pattern) <- c(du[1],du[2])
  } else if (!area.mean.expl)
             pattern <- c(COEFS[2:dc[1],1]) * attr(eof,'eigenvalues')[ip]
  if (verbose) print('predict')
  #browser()
  z <- zoo(predict(model,newdata=predat),order.by=index(Xam))
  ds <- z; ds <- attrcp(y0,ds,ignore='names')
#  if (!area.mean.expl) {
#    pattern <- attrcp(eof,pattern,ignore=c('longitude','latitude','names'))
#  } 

  caldat <- zoo(as.matrix(caldat),order.by=index(X))
  if (verbose) print('Set attributes')
  attr(caldat,'calibration_expression') <- calstr
  attr(caldat,'stepwise_screening') <- swsm
  attr(ds,'calibration_data') <- caldat
  attr(ds,'fitted_values') <- zoo(model$fitted.values,
                                  order.by=index(X))
  class(attr(ds,'fitted_values')) <- class(y)
  attr(ds,'original_data') <- yam
  attr(ds,'evaluation') <- merge(zoo(yam),z)
  r2 <- var(coredata(model$fitted.values))/var(y,na.rm=TRUE)
  attr(r2,'description') <- ' var(fitted.values))/var(y)'
  attr(ds,'quality') <- r2
  attr(ds,'variable') <- attr(y,'variable')
  attr(ds,'model') <- model
  attr(ds,'mean') <- mean(yam)
  attr(ds,'method') <- method
  attr(ds,'eof') <- eof
  attr(pattern,'longitude') <- attr(X0,'longitude')
  attr(pattern,'latitude') <- attr(X0,'latitude')
  attr(ds,'pattern') <- pattern
  attr(ds,'dimensions') <- c(du[1],du[2])
  attr(ds,'longitude') <- attr(y,'longitude')
  attr(ds,'latitude') <- attr(y,'latitude')
  attr(ds,'aspect') <- 'downscaled'
  attr(ds,'source') <- attr(X0,'source')
  attr(ds,'type') <- "downscaled results"
  attr(ds,'history.predictand') <- attr(y0,'history')
  attr(ds,'history') <- history.stamp(X0)
  class(ds) <- c("ds",class(y))
  rm("y0","X0")
  if (verbose) print('--- exit DS.default ---')
  invisible(ds)
}
