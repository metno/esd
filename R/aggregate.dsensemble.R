## Aggregate ensemble performs similar function to expandpca but applies the functions to the 
## results after the transform from PCA to field and hence is a bit slower.

## Function for extracting the subset from PCs stored as zoo
# not exported
subset.pc <- function(x,ip=NULL,it=NULL,verbose=FALSE) {
  if (verbose) print('subset.pc')
  d <- dim(x)
  if (!is.null(it)) {
    if (verbose) print('subset it')
    if ((is.numeric(it) | is.integer(it)) & is.dates(index(x))) {
      it <- c(as.Date(paste(it,'01-01',sep='-')),
              as.Date(paste(it,'12-31',sep='-')))
    }
    x <- window(x,start=min(it),end=max(it))
  }
  if (!is.null(ip)) {
    if (verbose) print('subset pattern')
    x <- x[,ip]
    d <- dim(x)
  }
  dim(x) <- c(length(index(x)),d[2])
  if(verbose) print(dim(x))
  return(x)
}

#' aggregate.dsensemble
#' 
#' The function \code{aggregate.dsensemble} is used to extract and aggregate information from PCA-based dsensemble-objects, 
#' to expand PCA into station data or EOFs into field data.
#' It applies operations through the argument 'FUNX' to estimate quantities such as the ensemble mean and 
#' takes care of extra house keeping, such as attributes with meta data. \code{aggregate.dsensemble}
#' is similar to \code{expandpca} but applies 'FUNX' for each time step after it has been expanded by matrix multiplication 
#' of the singular vectors (X = U %*% diag(W) %*% t(V)) and is therefore slower. 
#' \code{map.dsensemble} is a wrapper for \code{map} that uses one of the two former routines to distill selected
#' information.
#' @seealso aggregate expandpca map.dsensemble
#' @aliases aggregate.dsensemble
#'
#' @param x an object of type 'dsensemble'
#' @param it time index (see \code{\link{subset}})
#' @param FUN function applied to aggregate in time
#' @param FUNX function applied aggregate ensemble members
#' @param verbose if TRUE print progress
#' @param anomaly if FALSE add the mean value stored as attribute in x
#' @param qc if TRUE perform quality control and exclude ensemble members with incomplete or suspect data
#' @param test if TRUE perform test on one GCM simulation
#' @param eof if TRUE and if the dsensemble object contains a spatially aggregated pattern, x$eof, 
#' expand with x$eof into an ensemble of field objects to which FUNX is applied. 
#' If FALSE, expand using the original pattern x$pca into an ensemble of station objects to be aggregated. 
#'
#' @exportS3Method
#' @export aggregate.dsensemble
aggregate.dsensemble <- function(x,...,it=NULL,im=NULL,FUN=NULL,FUNX='mean',verbose=FALSE,anomaly=FALSE,
                                 qc=TRUE,test=FALSE,eof=TRUE) {
  ## Get the spatial weights
  if (verbose) print('aggregate.ensemble')
  if ((eof) & (!is.null(x$eof))) x$pca <- x$eof
  ##x$eof <- gridmap(x$pca)
  if (test) print('--TEST ON ONE GCM simulation--')
  if (inherits(x,'pca')) UWD <- x$pca else UWD <- x$eof
  if (verbose) print(names(attributes(UWD)))
  ## Eigenvalues
  D <- attr(UWD,'eigenvalues')
  ## Create a matrix with only the GCM time series
  if (verbose) print('PCA/EOF-based ensemble')
  X <- x
  X$info <- NULL; X$pca <- NULL; X$eof <- NULL
  if (verbose) for (ii in 1:length(X)) print(dim(X[[ii]]))
  ## Dimension of each downscaled GCM results
  if (verbose) print(paste('subset.pc, it=',it))
  if(!is.null(it)) {
    V <- lapply(X,FUN='subset.pc',it=it)
  } else {
    V <- coredata(X)
  }
  if (verbose) print(paste('Interval',paste(range(index(V[[1]])),collapse=' - ')))
  gcmnames <- attr(X, "model_id")
  if(!qc) {
    nok <- c()
  } else {
    ## Check if the ensemble members have the same size - if not, only keep the ones with most common sizes
    if (verbose) print('Check ensemble member size')
    n <- length(names(X))
    if (verbose) print(paste('Original length of X is',n))
    memsiz <- rep("?",n)
    for (i in 1:n) memsiz[i] <- paste(dim(X[[i]]),collapse='x')
    memsiztab <- table(memsiz)
    if (verbose) print(memsiztab)
    memkeep <- rownames( memsiztab)[as.numeric(memsiztab)==max(as.numeric(memsiztab))]
    if (verbose) print(memkeep)
    nok <- sort((1:n)[-grep(memkeep,memsiz)],decreasing = TRUE)
    if (verbose) print(nok)
    ## Check for strange values
    memsiz <- rep("?",n)
    for (i in 1:n) memsiz[i] <- paste(dim(V[[i]]),collapse='x')
    memsiztab <- table(memsiz)
    if (verbose) print(memsiztab)
    d <- dim(V[[1]])
    for (i in seq(n,1,by=-1)) {
      if (max(abs(V[[i]]),na.rm=TRUE) > 10)  {
        nok <- c(nok, i)
      }
    }
    ## Exclude ensemble members that failed quality control
    if (verbose) print(paste('Before quality control: original number of members=',n))
    for (i in sort(unique(nok), decreasing=TRUE)) {
      print(range(V[[i]],na.rm=TRUE)); print(dim(V[[i]]))
      print(paste(i,'Remove suspect or incomplete results: ',gcmnames[i]))
      V[[i]] <- NULL
      gcmnames <- gcmnames[-i]
    }
    n <- length(V)
    if (verbose) print(paste('After quality control: number of members=',n))
  }
  #=========================================
  ## KMP 2023-05-23: Changing how 'im' subsets from the ensemble so that the same member can be sampled multiple times
  if(is.null(im)) im <- 1:length(V) else {
    if(is.logical(im)) im <- which(im)
    if(length(nok)>0) {
      if (verbose) print(paste('Before quality control: selected number of members=',length(im)))
      im <- unlist(sapply(im, function(i) which(gcmnames %in% attr(X,"model_id")[i])))
      if (verbose) print(paste('After quality control: selected number of members=',length(im)))
    }
    gcmnames <- gcmnames[im]
  }
  n <- length(im)
  #=========================================
  if (verbose) print(paste('Number of members in selected ensemble=',length(im)))
  #if (verbose) {print(names(V)); print(c(d,n,length(unlist(V)))); print(paste('FUNX=',FUNX))}
  lengths <- rep(NA,length(V))
  for (i in 1:length(V)) lengths[i] <- length(index(V[[i]]))
  #if (length(table(lengths[im]))>1) browser()
  
  ## Aggregate statistics over ensemble members
  if (verbose) print('Aggregate ensemble statistics')
  U <- attr(UWD,'pattern')
  if (!is.null(dim(U))) {
    dU <- dim(U) 
  } else {
    dU <- c(1,1,length(U)) # If there is only one single station
    dim(U) <- dU
  }
  
  #if (verbose) {print(d); print(dU)}
  if (inherits(UWD,'eof')) {
    #if (verbose) {print('eof'); print(dU)}
    dim(U) <- c(dU[1]*dU[2],dU[3])
  }
  if (verbose) {
    #print('Matrix multiplication')
    #str(U); str(D); str(V)
  }
  ## Loop through each time step - aggregate ensemble statistics for each time step
  if(length(FUNX)>1) {
    ty <- seq(min(sapply(V, function(x) min(index(x)))), 
                max(sapply(V, function(x) max(index(x)))))
    eval(parse(text=paste0(paste0("Y.",FUNX,collapse=" <- "),
                           " <- matrix(rep(NA,dU[1]*dU[2]*length(ty)),length(ty),dU[1]*dU[2])")))
    nvalid <- rep(0,length(ty))
    if(verbose) print(paste("Calculating ensemble statistics for",length(im),"ensemble members"))
    for (i.t in seq_along(ty)) {
      if(verbose) print(paste("Calculating ensemble statistics for time step",i.t,"of",length(ty)))
      ## loop through each ensemble member
      z <- matrix(rep(NA,dU[1]*dU[2]*n),dU[1]*dU[2],n)
      #=========================================
      ## KMP 2023-05-23: Changing how 'im' subsets from the ensemble so that the same member can be sampled multiple times
      for (i.m in seq_along(im)) {
        v <- V[[im[i.m]]]
        if(ty[i.t] %in% index(v)) {
          z[,i.m] <- v[index(v)==ty[i.t],] %*% diag(D) %*% t(U)
        }
      }
      #=========================================
      for(f in FUNX) {
        suppressWarnings(eval(parse(text=paste0("yf <- apply(z,1,f,na.rm=TRUE)"))))
        yf[!is.finite(yf)] <- NA
        eval(parse(text=paste0("Y.",f,"[i.t,] <- yf")))
      }
      nvalid[i.t] <- sum(!is.na(yf))
      rm(v,z); gc(reset=TRUE)
    }
    for(f in FUNX) eval(parse(text=paste0("attr(Y.",f,",'FUN') <- '",f,"'")))
    eval(parse(text=paste0("Y <- list(",paste0("'",FUNX,"'=Y.",FUNX, collapse=", "),")")))
    
    for(f in names(Y)) {
      ## Add mean and insert into zoo frame
      if (f != 'mean') anomaly <- TRUE
      y <- Y[[f]]
      if (!anomaly) {
        if (verbose) print('add mean field')
        y <- t(t(y) + c(attr(UWD,'mean')))
      }
      y <- as.zoo(y, order.by=ty)
      Y[[f]] <- y
    }
    Y <- attrcp(UWD,Y)
    attr(Y,'timeindex') <- ty
    attr(Y,'time') <- range(index(subset(X[[1]],it=it)))
    class(Y) <- c("list", "dsensemblestatistics", class(UWD)[length(class(UWD))-1])
    if (inherits(UWD,'eof')) {
      if (verbose) print('Use dimensions and lon/lat from EOFs')
      attr(Y,'dimensions') <- c(attr(x$eof,'dimensions')[1:2],length(index(V)))
      attr(Y,'longitude') <- lon(UWD)
      attr(Y,'latidude') <- lat(UWD)
    }
    attr(Y, 'model_id') <- gcmnames
    attr(Y, 'nvalid') <- nvalid
    attr(Y , 'mean') <- NULL
    
    if(!is.null(x$pca)) {
      if(inherits(x$pca, "season")) {
        attr(Y,'season') <- season(x$pca)[1]
      } else {
        attr(Y,'season') <- class(x$pca)[length(class(x$pca))-1]
      }
    } else {
      if(inherits(x, "season")) {
        attr(Y,'season') <- season(x[[2]])[1]
      } else {
        attr(Y,'season') <- class(x[[2]])[length(class(x[[2]]))-1]
      }
    }
  } else {
    Y <- matrix(rep(NA,dU[1]*dU[2]*d[1]),d[1],dU[1]*dU[2])
    for (i.t in 1:d[1]) { 
      ## loop through each ensemble member
      z <- matrix(rep(NA,dU[1]*dU[2]*n),dU[1]*dU[2],n)
      for (im in 1:n) { 
        v <- V[[im]]
        z[,im] <- v[i.t,] %*% diag(D) %*% t(U)
      }
      Y[i.t,] <- apply(z,1,FUNX)
      rm(v,z); gc(reset=TRUE)
    }
  
    if (!is.null(FUNX)) {
      if (FUNX != 'mean') anomaly <- TRUE; 
      if (verbose) print(c(FUNX,anomaly))
    }
    
    ## Add mean and insert into zoo frame
    if (!anomaly) {
      if (verbose) print('add mean field')
      Y <- t(t(Y) + c(attr(UWD,'mean')))
    }
    # Not right if FUNX is defined and time mean has been applied:
    if(nrow(V[[1]])==length(index(V[[1]]))) {
      Y <- zoo(Y,order.by=index(V[[1]]))
    } else {
      Y <- zoo(Y,order.by=seq(nrow(V[[1]])))
    }
    Y <- attrcp(UWD,Y)
    attr(Y,'time') <- range(index(subset(X[[1]],it=it)))
    class(Y) <- class(UWD)[-1]
    if (inherits(UWD,'eof')) {
      if (verbose) print('Use dimensions and lon/lat from EOFs')
      attr(Y,'dimensions') <- c(attr(x$eof,'dimensions')[1:2],length(index(V)))
      attr(Y,'longitude') <- lon(UWD)
      attr(Y,'latidude') <- lat(UWD)
      class(Y)[1] <- 'field'
    }
    attr(Y, 'model_id') <- gcmnames
    attr(Y,'mean') <- NULL
    if(!is.null(FUN)) Y <- map(Y, FUN=FUN, plot=FALSE)
  }
  if (verbose) {print('exit aggregate.dsensemble'); print(dim(Y))}
  gc(reset=TRUE)
  return(Y)
}
