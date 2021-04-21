#' @exportS3Method
#' @export
DS.mvcomb <- function(y, X, X2=NULL, plot=FALSE, it=NULL, method="lm",
                      swsm="step", m=5, rmtrend=TRUE, ip=1:10, weighted=TRUE, 
                      pca=TRUE, npca=20, biascorrect=FALSE, verbose=FALSE, ...) {
  if(verbose) print("DS.mvcomb")
  if(!inherits(X,"mvcomb")) {
    if(!is.null(X2)) {
      X <- mvcombine(X, X2)
    } else {
    }
  }
  if(!inherits(X,"eof")) X <- EOF(X, verbose=verbose)
  dy <- dim(y); if (is.null(dy)) dy <- c(length(y),1)
  dx <- dim(X); if (is.null(dx)) dx <- c(length(X),1)
  cls <- c(class(y)[!class(y) %in% class(X)], class(X))
  y0 <- y; X0 <- X
  
  if(!inherits(y,"pca")) {
    z.out <- DS(y, X, plot=plot, it=it, method=method,
                swsm=swsm, m=m, rmtrend=rmtrend, ip=ip, weighted=weighted, 
                pca=pca, npca=npca, biascorrect=biascorrect, verbose=verbose, ...)
    if (verbose) print('Predictor pattern')
    if(method=="bas.lm") {
      coefs <- coef(attr(z.out, "model"))$postmean[-1]
    } else {
      coefs <- coef(attr(z.out,'model'))[-1]
    }
    
    predpatt <- list()
    for(i in seq(length(attr(X0,'pattern')))) {
      U <- attr(X,'pattern')[[i]]
      du <- dim(U)
      if (verbose) {print(paste('pattern',i,'dimension')); print(du); str(U)}
      if (length(du)==3) dim(U) <- c(du[1]*du[2],du[3])
      if ( du[3]>=length(coefs) & length(coefs)>1 ) {
        pattern <- U[,1:length(coefs)] %*% coefs
      } else if (length(coefs)==1) {
        pattern <- U[,1:length(coefs)] * coefs
      } else if (length(coefs)==0) {
        pattern <- matrix(0, nrow=nrow(U), ncol=1)
      }
      dim(pattern) <- c(du[1],du[2])
      predpatt[[i]] <- pattern
    }
    attr(predpatt,'longitude') <- attr(X0,'longitude')
    attr(predpatt,'latitude') <- attr(X0,'latitude')
    attr(predpatt,'variable') <- attr(X0,'variable')
    attr(predpatt,'unit') <- attr(X0,'unit')
    attr(predpatt,'greenwich') <- attr(X0,'greenwich')
    attr(z.out, "predictor.pattern") <- predpatt
  } else {
    if (verbose) print('Prepare output data')
    if(is.null(npca)) npca <- dy[2] else npca <- min(npca, dy[2])
    if (!is.null(attr(X0,'n.apps'))) {
      if (verbose) print(attr(X0,'n.apps'))
      Xp <- attr(X0,'appendix.1')
    } else Xp <- X0
    z.out <- matrix(rep(NA,dx[1]*npca),dx[1],npca)
    fit.val <- z.out
    dxp <- dim(Xp); if (is.null(dxp)) dxp <- c(length(Xp),1)
    zp.out <- matrix(rep(NA,dxp[1]*npca),dxp[1],npca)
    model <- list(); predpatt <- list(); cval <- list()
    if (!verbose) pb <- txtProgressBar(style=3)
    if(verbose) print("Perform downscaling")
    for (i in 1:npca) {
      if (!verbose) setTxtProgressBar(pb,i/npca) 
      ys <- as.station(zoo(y[,i]), loc=loc(y)[i], param=varid(y)[i],
                       unit=unit(y)[i], lon=lon(y)[i], lat=lat(y)[i],
                       alt=alt(y)[i], cntr=cntr(y)[i], stid=stid(y)[i])
      class(ys) <- c('station', class(y)[-c(1:2)])
      if (verbose) {print(class(ys)); print(class(X))}
      z.i <- DS(ys, X, plot=plot, it=it, method=method,
                swsm=swsm, m=m, rmtrend=rmtrend, ip=ip, weighted=weighted, 
                pca=pca, npca=npca, biascorrect=biascorrect, verbose=verbose, ...)
      
      ## Collect the projections in a matrix:
      z.out[is.finite(z.i),i] <- coredata(z.i)
      fit.val[is.finite(z.i),i] <- attr(z.i,'fitted_values')
      if (!is.null(attr(X,'n.apps')))
        zp.out[,i] <- attr(z.i,'appendix.1')
      
      model[[i]] <- attr(z.i,'model')
      #eof[[i]] <- X
      
      ## Calculate predictor pattern
      if(method=="bas.lm") {
        coefs <- coef(attr(z.i, "model"))$postmean[-1]
      } else {
        coefs <- coef(attr(z.i,'model'))[-1]
      }
      
      predpatt.i <- list()
      for(j in seq(length(attr(X,'pattern')))) {
        U <- attr(X,'pattern')[[j]]
        du <- dim(U)
        if (verbose) {print(paste('pattern',j,'dimension')); print(du); str(U)}
        if (length(du)==3) dim(U) <- c(du[1]*du[2],du[3])
        if ( du[3]>=length(coefs) & length(coefs)>1 ) {
          pattern.j <- U[,1:length(coefs)] %*% coefs
        } else if (length(coefs)==1) {
          pattern.j <- U[,1:length(coefs)] * coefs
        } else if (length(coefs)==0) {
          pattern.j <- matrix(0, nrow=nrow(U), ncol=1)
        }
        dim(pattern.j) <- c(du[1],du[2])
        predpatt.i[[j]] <- pattern.j
      }
      attr(predpatt.i,'longitude') <- attr(X,'longitude')
      attr(predpatt.i,'latitude') <- attr(X,'latitude')
      attr(predpatt.i,'variable') <- attr(X,'variable')
      attr(predpatt.i,'unit') <- attr(X,'unit')
      attr(predpatt.i,'greenwich') <- attr(X,'greenwich')
      predpatt[[i]] <- predpatt.i
      
      ## Keep cross-validation
      if (!is.null(attr(z.i,'evaluation'))) { ## REB 2015-03-27
        if (verbose) print('extract cross-validation')
        if (i==1) cval <- attr(z.i,'evaluation') else
          cval <- merge(cval,attr(z.i,'evaluation'))
      } else cval <- NULL
      
    }
    if (!is.null(cval)) names(cval) <- paste(c('X','Z'),"PCA",sort(rep(1:npca,2)),sep='.')
    z.out <- zoo(z.out, order.by=index(z.i))
    z.out <- attrcp(y0, z.out)
    attr(z.out, "predictor.pattern") <- predpatt
    attr(z.out, "evaluation") <- cval
    attr(z.out,'calibration_data') <- attr(z.i,'calibration_data')
    attr(z.out,'fitted_values') <- zoo(fit.val,order.by=index(attr(z.i,'fitted_values')))
    class(attr(z.out,'fitted_values')) <- class(y0)
    attr(z.out,'model') <- model
    attr(z.out,'eof') <- X0
    attr(z.out,'original_data') <- y0
    attr(z.out,'dimensions') <- c(du[1],du[2])
    #browser()
    attr(z.out,'pattern') <- attr(y0,'pattern')[,1:npca]
    attr(z.out,'eigenvalues') <- attr(y0,'eigenvalues')[1:npca]
    #attr(z.out,'variable') <- varid(y0)
    #attr(z.out,'unit') <- attr(y0,'unit')
    #attr(z.out,'mean') <- attr(y0,'mean') # + offset
    #attr(z.out,'max.autocor') <- attr(y0,'max.autocor')
    #attr(z.out,'tot.var') <- attr(y0,'tot.var')
    #attr(z.out,'longitude') <- lon(y0)
    #attr(z.out,'latitude') <- lat(y0)
    attr(z.out,'source') <- attr(X,'source')
    attr(z.out,'history') <- history.stamp(y)
  }
  class(z.out) <- c("ds",cls)
  invisible(z.out)
}
