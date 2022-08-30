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

#' trend.dsensemble
#' 
#' The function \code{trend.dsensemble} is used to calculate trends of dsensemble objects and summarise the results in terms of ensemble statistics: the minimum (min), fifth percentile (q5), median, mean, 95th percentile (q95) and maximum trend (max), as well as the number and fraction of ensemble members with a positive (n.pos, f.pos) and negative (n.neg, f.neg) trend. The results of the trend analysis can be displayed visually with the function \code{trendplot.dsensemble}, which is called from \code{trend.dsensemble} if the argument plot=TRUE.
#'
#' @seealso trendplot.dsensemble map.dsensemble aggregate.dsensemble 
#'
#' @param x an object of type 'dsensemble'
#' @param it time index (see \code{\link{subset}})
#' @param verbose if TRUE print progress
#' @param plot if TRUE show results visually on a map
#' @param eof if TRUE and if the dsensemble object contains a spatially aggregated pattern, x$eof, 
#' expand with x$eof into an ensemble of field objects to which the trend analysis is applied. 
#' If FALSE, expand using the original pattern x$pca into an ensemble of station objects to be aggregated. 
#'
#' @export trend.dsensemble
trend.dsensemble <- function(x,...,it=NULL,verbose=FALSE,plot=FALSE,eof=FALSE) {
  ## Get the spatial weights
  if (verbose) print('trend.dsensemble')
  if (verbose) print('Calculate ensemble statistics of dsensemble trends')
  stopifnot(inherits(x, "dsensemble"))
  if(inherits(x,"station")) {
    if (verbose) print('station ensemble')
    trends <- matrix(NA, ncol=length(x), nrow=ncol(x[[1]]))
    colnames(trends) <- names(x)
    rownames(trends) <- attr(x[[1]], "model_id")
    for(i in seq_along(x)) trends[,i] <- apply(x[[i]], 2, trend.coef)
    trends.stats <- cbind(apply(trends, 2, min), apply(trends, 2, q5), 
                          apply(trends, 2, median), apply(trends, 2, mean),
                          apply(trends, 2, q95), apply(trends, 2, max),
                          apply(trends, 2, function(x) sum(x>0)), 
                          apply(trends, 2, function(x) sum(x<0)),
                          apply(trends, 2, function(x) sum(x>0)/length(x)), 
                          apply(trends, 2, function(x) sum(x<0)/length(x)))
    colnames(trends.stats) <- c("min","q5","median","mean","q95","max",
                                "n.pos","n.neg","f.pos","f.neg")
    rownames(trends.stats) <- attr(x[[1]], "model_id")
    trends.stats <- data.frame(trends.stats)
  } else {
    if (eof) {
      if (is.null(x$eof)) x$eof <- gridmap(x$pca)
      x$pca <- x$eof
    }
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
    ## Check if the ensemble members have the same size - if not, 
    ## only keep the ones with most common sizes
    if (verbose) print('Check ensemble member size')
    n <- length(names(X))
    if (verbose) print(paste('Original length of X is',n))
    gcms <- attr(X, "model_id")
    memsiz <- rep("?",n)
    for (i in 1:n) memsiz[i] <- paste(dim(X[[i]]),collapse='x')
    memsiztab <- table(memsiz)
    if (verbose) print(memsiztab)
    memkeep <- rownames( memsiztab)[as.numeric(memsiztab)==max(as.numeric(memsiztab))]
    if (verbose) print(memkeep)
    im <- sort((1:n)[-grep(memkeep,memsiz)],decreasing = TRUE)
    if (verbose) print(im)
    for (ix in im) {
      X[[ix]] <- NULL
      gcms <- gcms[-ix]
    }
    n <- length(names(X))
    if (verbose) print(paste('New length of X is',n))
    ## Only select the selected time interval - saves time
    if (is.null(it)) it <- range(index(X[[1]]))
    V <- lapply(X,FUN='subset.pc',it=it)
    if (verbose) print(paste('Interval',paste(range(index(V[[1]])),collapse=' - ')))
    memsiz <- rep("?",n)
    for (i in 1:n) memsiz[i] <- paste(dim(V[[i]]),collapse='x')
    memsiztab <- table(memsiz)
    if (verbose) print(memsiztab)
    d <- dim(V[[1]])
    ## Quality control
    if (verbose) print(paste('Before quality control: original number of members=',n))
    for (i in seq(n,1,by=-1)) {
      if(verbose) {print(range(V[[i]],na.rm=TRUE)); print(dim(V[[i]]))}
      if (max(abs(V[[i]]),na.rm=TRUE) > 10)  {
        if(verbose) print(paste(i,'Remove suspect results'))
        V[[i]] <- NULL
        gcms <- gcms[-i]
      }
    }
    n <- length(V)
    if (verbose) print(paste('After quality control: new number of members=',n))
    if (verbose) {print(names(V)); print(c(d,n,length(unlist(V))))}
    lengths <- rep(NA,length(V))
    for (i in 1:length(V)) lengths[i] <- length(index(V[[i]]))
    if (length(table(lengths))>1) browser()
    
    ## Aggregate statistics over ensemble members
    if (verbose) print('Aggregate ensemble statistics')
    
    U <- attr(UWD,'pattern')
    if (!is.null(dim(U))) {
      dU <- dim(U) 
    } else {
      dU <- c(1,1,length(U)) # If there is only one single station
      dim(U) <- dU
    }
    
    if (verbose) {print(d); print(dU)}
    if (inherits(UWD,'eof')) {
      if (verbose) {print('eof'); print(dU)}
      dim(U) <- c(dU[1]*dU[2],dU[3])
    }
    if (verbose) {
      print('Matrix multiplication')
      str(U); str(D); str(V)
    }
    
    ## Loop through each ensemble member
    trends <- matrix(NA, ncol=dim(U)[1], nrow=n)
    rownames(trends) <- names(V)
    for (im in 1:n) {
      if(verbose) print(im)
      z <- V[[im]] %*% diag(D) %*% t(U)
      trends[im,] <- apply(z, 2, trend.coef)
    }
    # trends.stats <- cbind(apply(trends, 2, min), apply(trends, 2, q5), 
    #                       apply(trends, 2, median), apply(trends, 2, mean),
    #                       apply(trends, 2, q95), apply(trends, 2, max),
    #                       apply(trends, 2, function(x) sum(x>0)), 
    #                       apply(trends, 2, function(x) sum(x<0)),
    #                       apply(trends, 2, function(x) sum(x>0)/length(x)), 
    #                       apply(trends, 2, function(x) sum(x<0)/length(x)))
    # colnames(trends.stats) <- c("min","q5","median","mean","q95","max",
    #                             "n.pos","n.neg","f.pos","f.neg")
    # trends.stats <- data.frame(trends.stats)
    # trends.stats <- attrcp(UWD,trends.stats)
    # attr(trends.stats, "model_id") <- attr(x, "model_id")
    # if(eof) {
    #   attr(trends.stats,'dimension') <- attr(UWD,'dimension')[1:2]
    # } else {
    #   attr(trends.stats,'dimension') <- attr(UWD,'dimension')[1]
    # }
    attr(trends,'time') <- range(index(subset(X[[1]],it=it)))
    trends <- attrcp(UWD,trends)
    attr(trends, "model_id") <- gcms#attr(x, "model_id")
    if(eof) {
     attr(trends,'dimension') <- attr(UWD,'dimension')[1:2]
    } else {
     attr(trends,'dimension') <- attr(UWD,'dimension')[1]
    }
    attr(trends,'time') <- range(index(subset(X[[1]],it=it)))
  }
  
  #if(plot) trendplot.dsensemble(trends.stats, ...)
  if(plot) trendplot.dsensemble(trends, ...)
  if (verbose) {print('exit trend.dsensemble')}
  #return(trends.stats)
  return(trends)
}

#' trendplot.dsensemble
#' 
#' The function \code{trendplot.dsensemble} is used to visualise the trends of downscaled ensembles, calculated with \code{trend.dsensemble}. 
#'
#' @seealso trend.dsensemble map.dsensemble aggregate.dsensemble 
#'
#' @param x A \code{data.frame} containing ensemble statistics of linear trends. Output from the function \code{trend.dsensemble}.
#' @param statistic Ensemble statistic to show on the map. By default, statistic=\code{"mean"}, the ensemble mean of the trends. Other alternatives: \code{"median"}, \code{"min"} (minimum), \code{"max"} (maximum), \code{"q5"} (5th percentile), \code{"q95"} (95th percentile), \code{"n.pos"} and \code{"n.neg"} (number of ensemble members with positive or negative trends), and \code{"f.pos"} and \code{"f.neg"} (fraction of ensemble members with positive or negative trends).
#' @param robustness Ensemble statistic to use as a measure of trend robustness. Default: \code{"f"}, the fraction of ensemble members with the same sign of trends. Other options: \code{"n"}, the number of ensemble members with the same sign of trends, and any ensemble statistic in the input \code{x}. If NULL, the robustness is not estimated and displayed.
#' @param threshold Threshold for estimating trend robustness. 
#' @param threshold.lower if TRUE, \code{threshold} is used as a lower threshold for \code{robustnes}, otherwise it is used as an upper threshold. 
#' @param verbose if TRUE print progress
#' @param pch A vector of plotting characters or symbols: see points.
#' @param cex A numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default. This works as a multiple of par("cex"). NULL and NA are equivalent to 1.0. Note that this does not affect annotation: see below. 
#' @param lwd a vector of line widths, see \code{link{par}}. 
#' @param colbar The colour scales defined through colscal. Users can specify the colour ‘pal’*ette (‘pal’), the number of breaks (‘n’), values of ‘breaks’, and the position of the color bar from the main plot (‘pos’). The ‘rev’ argument, will produce a reversed color bar if set to TRUE. The other arguments (‘type’,‘h’ and ‘v’) are more specific to col.bar and are used only if argument ‘fancy’ is set to TRUE (not yet finished). colbar=NULL is used if the colourbar is not to be shown. Also use colbar=NULL to present several maps in one figure (e.g. with par(mfcol=c(2,2))).
#' @param pch.robustness Plotting character to show robustness
#' @param cex.robustness Scaling of plotting characters to show robustness
#' @param lwd.robustness Line widths for plotting characters to show robustness
#' @param col.robustness Color of plotting characters to show robustness. If NULL, the same color is used as for the trend. \code{col} can be one color, e.g. "black", or a list of two colors, one for positive and one for negative trends, e.g., list("pos"="red", "neg"="blue")).
#' @param projection Projections: c("lonlat","sphere","np","sp") - the latter gives stereographic views from the North and south poles.
#'
#' @export trendmap.dsensemble
trendmap.dsensemble <- function(trends,#.stats, 
                                statistic="mean", new=TRUE,
                                robustness="f", threshold=0.9, threshold.lower=TRUE,
                                pch=19, cex=0.9, lwd=1, colbar=list(show=TRUE),
                                pch.robustness=1, cex.robustness=0.9,
                                lwd.robustness=1.2, main=NULL,
                                bg="grey55", col.robustness="black",
                                projection="lonlat", ..., verbose=FALSE) {
  if(verbose) print("trendplot.dsensemble")
  #X <- trends.stats[[statistic]]
  #trends.stats <- cbind(apply(trends, 2, min), apply(trends, 2, q5), 
  #                      apply(trends, 2, median), apply(trends, 2, mean),
  #                      apply(trends, 2, q95), apply(trends, 2, max),
  #                      apply(trends, 2, function(x) sum(x>0)), 
  #                      apply(trends, 2, function(x) sum(x<0)),
  #                      apply(trends, 2, function(x) sum(x>0)/length(x)), 
  #                      apply(trends, 2, function(x) sum(x<0)/length(x)))
  #colnames(trends.stats) <- c("min","q5","median","mean","q95","max",
  #                            "n.pos","n.neg","f.pos","f.neg")
  #trends.stats <- data.frame(trends.stats)
  #trends.stats <- attrcp(UWD,trends.stats)
  #attr(trends.stats, "model_id") <- attr(x, "model_id")
  #if(eof) {
  #  attr(trends.stats,'dimension') <- attr(UWD,'dimension')[1:2]
  #} else {
  #  attr(trends.stats,'dimension') <- attr(UWD,'dimension')[1]
  #}
  #attr(trends.stats,'time') <- range(index(subset(X[[1]],it=it)))
  eval(parse(text=paste0("X <- apply(trends, 2,", statistic, ")")))

  attr(X,'longitude') <- attr(trends,'longitude')
  attr(X,'latitude') <- attr(trends,'latitude')
  attr(X,'variable') <- paste(attr(trends,'variable')[[1]],statistic,'trend')
  attr(X,'unit') <- paste(attr(trends,'unit')[[1]],'/decade')
  #attr(X,'longitude') <- attr(trends.stats,'longitude')
  #attr(X,'latitude') <- attr(trends.stats,'latitude')
  #attr(X,'variable') <- paste(attr(trends.stats,'variable')[[1]],statistic,'trend')
  #attr(X,'unit') <- paste(attr(trends.stats,'unit')[[1]],'/decade')
  if(is.null(dim(X))) {
    if(verbose) print(paste('Plot',statistic,'of the dsensemble trends'))
    dim(X) <- c(length(X), 1)
    class(X) <- c("station", "trend", "dsensemble")
    if(is.null(colbar$pal)) {
      if(all(X<0)) colbar$pal <- "bu" else if (all(X>0)) colbar$pal <- "rd" 
    } 
    if(is.null(colbar$breaks)) {
      if(all(X<0) | all(X>0)) {
        colbar$breaks <- pretty(range(X), n=10)
      } else {
        colbar$breaks <- pretty(c(-1,1)*max(abs(X)), n=10)
      }
    }
    #par(bg="black") ## TRY BLACK BACKGROUND
    # if(!is.null(par)) CONTINUE HERE!!
    
    if(is.null(main)) main <- paste0(attr(X,"variable"),"\n(", attr(X,'unit'),")")
      
    if(new) dev.new()
    par(bg=bg)
    map.X <- map(X, pch=pch, cex=cex, lwd=lwd, FUN="mean", colbar=colbar,
                 main=paste0(attr(X,"variable"),"\n(", attr(X,'unit'),")"),
                 new=FALSE, ...)
    colbar <- colbar.ini(X,colbar=colbar)
    if (verbose) print('Set colour scheme')
    wr <- round(strtoi(paste('0x',substr(colbar$col,2,3),sep=''))/255,2)
    wg <- round(strtoi(paste('0x',substr(colbar$col,4,5),sep=''))/255,2)
    wb <- round(strtoi(paste('0x',substr(colbar$col,6,7),sep=''))/255,2)
    col <- rep(colbar$col[1],length(X))
    for (i in 1:length(X)) {
      ii <- round(approx(0.5*(colbar$breaks[-1]+colbar$breaks[-length(colbar$breaks)]),
                         1:length(colbar$col),
                         xout=X[i],rule=2)$y)
      if (is.finite(ii)) {
        if (ii < 1) ii <- 1
        if (ii > length(colbar$col)) ii <- length(colbar$col)
        col[i] <- rgb(wr[ii],wg[ii],wb[ii],0.7)
      } else col[i] <- rgb(0.5,0.5,0.5,0.2)
    }
    points(lon(X), lat(X), col=col, pch=pch, cex=cex, lwd=lwd)
    if(!is.null(robustness)) {
      col.s <- col
      if(!is.null(col.robustness)) {
        if(is.list(col.robustness)) {
          if(!is.null(col.robustness$pos) &
             !is.null(col.robustness$neg)) {
            col.s <- rep(col.robustness$pos[[1]], length(X))
            col.s[X<0] <- col.robustness$neg[[1]]
          }
        } else {
          col.s <- rep(col.robustness[[1]], length(X))
        }
      }
      s <- NULL
      if(!is.null(robustness)) {
        if(verbose) print('Show statistical robustness of trend ensemble')
        if(is.null(threshold.lower)) threshold.lower <- TRUE
        if(robustness=="f") {
          if(is.null(threshold)) threshold <- 0.9
          #sig <- apply(trends.stats[,c("f.pos","f.neg")], 1, max)
          sig <- apply(trends, 2, function(x) max(sum(x>0), sum(x<0))/length(x))
        } else if(robustness=="n") {
          if(is.null(threshold)) threshold <- floor(length(X)*0.9)
          #sig <- apply(trends.stats[,c("n.pos","n.neg")], 1, max)
          sig <- apply(trends, 2, function(x) max(sum(x>0), sum(x<0)))
        } else {
          eval(parse(text=paste0("sig <- apply(trends, 2, ",robustness,")")))
        }
        if(is.null(threshold)) {
          warning('threshold for trend robustness not defined')
        } else {
          if(threshold.lower) {
            if(verbose) print(paste('robustness measure:',robustness,'>',threshold))
            s <- sig>=threshold
          } else {
            if(verbose) print(paste('robustness measure:',robustness,'<',threshold))
            s <- sig<=threshold
          }
        }
      }
      if(any(s)) points(lon(X)[s], lat(X)[s], col=col.s[s], 
                        pch=pch.robustness, cex=cex.robustness,
                        lwd=lwd.robustness)
    }
  } else {
    ## NOT FINISHED!
    browser()
    x <- X
    if (projection=="lonlat") {
      lonlatprojection(x=x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                       colbar=colbar,type=type,new=new,
                       verbose=verbose,gridlines=gridlines,...)
    } else if (projection=="sphere") {
      map2sphere(x=x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                 lonR=lonR,latR=latR,axiR=axiR,
                 type=type,gridlines=gridlines,
                 colbar=colbar,new=new,verbose=verbose,...)
    } else if (projection=="np") {
      map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                 lonR=lonR,latR=90,axiR=axiR,
                 type=type,gridlines=gridlines,
                 colbar=colbar,new=new,verbose=verbose,...)
    } else if (projection=="sp") {
      map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                 lonR=lonR,latR=-90,axiR=axiR,
                 type=type,gridlines=gridlines,
                 colbar=colbar,new=new,verbose=verbose,...)
    }
  }
}
