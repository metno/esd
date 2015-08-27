## R.E. Benestad
## Plot a map of the station locations, fields, EOFs, CCA results, correlation, composites, ...

##require(zoo)

map <- function(x,it=NULL,is=NULL,new=FALSE,...) UseMethod("map")

map.default<-function(x,FUN='mean',it=NULL,is=NULL,new=FALSE,
                      projection="lonlat",xlim=NULL,ylim=NULL,zlim=NULL,##n=15
                      colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                          show=TRUE,type="p",cex=2,h=0.6,v=1),
                      type=c("fill","contour"),gridlines=FALSE,
                      lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,...) {
    
    ## default with no arguments will produce a map showing available station
    ## data in the esd package.

    if (verbose) print('map.default')
    ## browser()
    if (is.logical(colbar)) colbar <- NULL
    ## If only a few items are provided in colbar - then set the rest to the default
    if (!is.null(colbar)) {
        colbar <- colbar.ini(x,FUN=FUN,colbar=colbar,verbose=verbose)
    } else col <- 'black'
    
    x <- subset(x,it=it,is=is)
    X <- attr(x,'pattern')

    ## if zlim is specified, then mask data outside this range
    if (!is.null(zlim)) {
        d <- dim(X)
        mask <- (X < min(zlim)) | (X > max(zlim))
        X[mask] <- NA
        dim(X) <- d
        if (verbose) {print(zlim); print(dim(X)); print(sum(mask))}
    }
    attr(X,'longitude') <- lon(x)
    attr(X,'latitude') <- lat(x)
    attr(X,'variable') <- attr(x,'variable')
    attr(X,'unit') <- attr(x,'unit')[1]
    attr(X,'source') <- attr(x,'source')
    attr(X,'variable') <- varid(x)
    if (inherits(X,'zoo')) attr(X,'time') <- range(index(x)) else
    if (!is.null(attr(x,'time'))) attr(X,'time') <- attr(x,'time')
    if (projection=="lonlat") lonlatprojection(x=X,xlim=xlim,ylim=ylim,n=n,
            colbar=colbar,verbose=verbose,
            type=type,new=new,
            gridlines=gridlines,...) else
    if (projection=="sphere") map2sphere(x=X,lonR=lonR,latR=latR,axiR=axiR,
            type=type,gridlines=gridlines,
            colbar=colbar,new=new,...) else
    if (projection=="np") map2sphere(X,lonR=lonR,latR=latR,axiR=axiR,
            type=type,gridlines=gridlines,
            colbar=colbar,new=new,...) else
    if (projection=="sp") map2sphere(X,lonR=lonR,latR=latR,axiR=axiR,new=new,
            type=type,gridlines=gridlines,
            colbar=colbar,...)
    
    invisible(X)
}

map.matrix <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                       xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                       colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                           pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                       type=c("fill","contour"),gridlines=FALSE,
                       lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
                       pattern=1,...) {
    
                                        # If x is provided, map only x...

                                        # default with no arguments will produce a map showing the station data in the esd package.

                                        #  image(lon(x),lat(x),x)
    if (verbose) print('map.matrix')
    if (!is.null(is)) x <- subset(x,is=is)  # if is is set, then call subset
    if (inherits(x,'zoo')) attr(x,'time') <- range(index(x)) 
    if (projection=="lonlat") lonlatprojection(x=x,new=new,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,
                                               type=type,gridlines=gridlines,verbose=verbose,...)  else
    if (projection=="sphere") map2sphere(x=x,new=new,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,verbose=verbose,...) else
    if (projection=="np") map2sphere(x,new=new,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,verbose=verbose,...) else
    if (projection=="sp") map2sphere(x,new=new,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,verbose=verbose,...)
    
                                        #map.station(NULL,...)
}

map.array <- function(x,FUN='mean',it=NULL,is=NULL,new=FALSE,
                      projection="lonlat",na.rm=TRUE,
                      xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
                      colbar=list(col=NULL,rev=FALSE,breaks=NULL,pos=0.05,
                          show=TRUE,type="r",cex=2,h=0.6,v=1),
                      type=c("fill","contour"),gridlines=FALSE,
                      lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,...) {
    if (verbose) print('map.array')
    if (!is.null(is)) x <- subset(x,is=is)  # if is is set, then call subset
    if (is.null(it)) {
        ## If it is NULL, then aggregate all of 3rd dimension
        D <- dim(x)
        x2d <- x
        dim(x2d) <- c(D[1]*D[2],D[3])
        z <- apply(x2d,1,FUN,na.rm=na.rm)
        z <- as.matrix(z)
        dim(z) <- c(D[1],D[2])
        str(z)
    } else  z <- x[,,it]
    d <- dim(z)

    ## if it is a vector of indices aggregate the selected indices
    if (length(d)==3) {
        dim(z) <- c(d[1]*d[2],d[3])
        z <- apply(z,2,FUN)
        dim(z) <- c(d[1],d[2])
    }
    attr(z,'longitude') <- lon(x)
    attr(z,'latitude') <- lat(x)
    attr(z,'variable') <- varid(x)
    attr(z,'unit') <- unit(x)[1]

    map(z,new=new,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,
        type=type,gridlines=gridlines,projection=projection,verbose=verbose,...)
}


map.comb <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                     xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
                     colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                         pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                     type=c("fill","contour"),gridlines=FALSE,
                     lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
                     pattern=1,...) {
    if (verbose) print('map.comb')
    stopifnot(inherits(x,'eof'))
    x <- subset(x,it=it,is=is)
    projection <- tolower(projection)
    ## if (is.null(col)) col <- colscal(col=colbar$pal,n=n-1,rev=colbar$rev) else
    ## if (length(col)==1) {
    ##   pal <- col
    ##   col <- colscal(col=pal,n=n-1,rev=colbar$rev)
    ## }
    if (is.null(varid(x))) attr(x,'variable') <- 'NA'
    ## if (tolower(varid(x))=='precip') col <- rev(col) 
    
    map.eof(x=x,xlim=xlim,ylim=ylim,zlim=zlim,pattern=pattern,
            n=n,projection=projection,colbar=colbar,new=new,
            lonR=lonR,latR=latR,axiR=axiR,type=type,
            gridlines=gridlines,verbose=verbose,...) -> result
    invisible(result)
}

map.eof <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                    xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
                    colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                        pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                    type=c("fill","contour"),gridlines=FALSE,
                    lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
                    pattern=1,...) {
    ## browser()
    if (verbose) print('map.eof')
    stopifnot(inherits(x,'eof'))
    ##x <- subset(x,it=it,is=is)
    projection <- tolower(projection)
    tot.var <- attr(x,'tot.var')
    D <- attr(x,'eigenvalues')
    var.eof <- 100* D^2/tot.var
    X <- attr(x,'pattern')[,,pattern]

    ## if zlim is specified, then mask data outside this range
    if (!is.null(zlim)) {
        d <- dim(X)
        mask <- (X < min(zlim)) | (X > max(zlim))
        X[mask] <- NA
        dim(X) <- d
        if (verbose) {print(zlim); print(dim(X)); print(sum(mask))}
    }
    ##str(x)
    attr(X,'longitude') <- attr(x,'longitude')
    attr(X,'latitude') <- attr(x,'latitude')
    attr(X,'variable') <- attr(x,'variable')
    attr(X,'unit') <- attr(x,'unit')[1]
    attr(X,'source') <- attr(x,'source')
    attr(X,'time') <- range(index(x))
    if ( (pattern==1) & !is.null(attr(x, "area.mean.expl")) )
        if (attr(x, "area.mean.expl"))
            type <- "fill"
    if (projection=="lonlat")
        lonlatprojection(x=X,it=it,xlim=xlim,ylim=ylim,
                         n=n,colbar=colbar,new=new,type=type,
                         gridlines=gridlines,verbose=verbose,...)
    else if (projection=="sphere")
        map2sphere(x=X,it=it,lonR=lonR,latR=latR,axiR=axiR,
                   type=type,gridlines=gridlines,
                   col=col,new=new,verbose=verbose,...)
    else if (projection=="np")
        map2sphere(X,it=it,lonR=lonR,latR=90,axiR=axiR,
                   type=type,gridlines=gridlines,
                   col=col,new=new,verbose=verbose,...)
    else if (projection=="sp")
        map2sphere(X,it=it,lonR=lonR,latR=-90,axiR=axiR,
                   type=type,gridlines=gridlines,
                   col=col,new=new,verbose=verbose,...)
    invisible(X)
}

map.ds <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
                   colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                       pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                   type=c("fill","contour"),gridlines=FALSE,
                   lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,...) {
    if (verbose) print('map.ds')
    stopifnot(inherits(x,'ds'))
    x <- subset(x,is=is)

    ## REB 2015-03-26
    if (inherits(x,'pca')) {
        map.pca(x,it=it,verbose=verbose,new=new,
                xlim=xlim,ylim=ylim,projection=projection,
                lonR=lonR,latR=latR,axiR=axiR,gridlines=gridlines,
                colbar=colbar) ##col=col,breaks=breaks)
        return()
    }
    projection <- tolower(projection)
    X <- attr(x,'pattern')
    if (is.list(X)) {
        X <- X[[1]]
    }
    
                                        # Check if there are several patterns: one for each month/seasons
    d <- dim(X)
    if (length(d)>2) {
        dim(X) <- c(d[1],d[2]*d[3])
        X <- colMeans(X)
        dim(X) <- c(d[2],d[3])
        attr(X,'longitude') <- lon(attr(x,'pattern'))
        attr(X,'latitude') <- lat(attr(x,'pattern'))
    }
    attr(X,'variable') <- varid(x)
    attr(X,'unit') <- unit(x)[1]
    
    unit <- attr(x,'unit')
    if ( (is.na(unit) | is.null(unit)) ) unit <- " "
    for (i in 1:length(unit)) {
        if ((unit[i]=='degree Celsius') | (unit[i]=='deg C') | (unit[i]=='degC'))
            unit[i] <- 'degree*C'
    }
    
    attr(X,'unit') <- unit
    attr(X,'source') <- attr(x,'source')
                                        #dim(X) <- attr(x,'dimensions')[1:2]
                                        #print(c(dim(X),length(attr(X,'longitude')),length(attr(X,'longitude'))))
    ##browser()
    if (projection=="lonlat") {
        lonlatprojection(x=X,n=n,colbar=colbar,verbose=verbose,
                         type='fill',gridlines=gridlines,new=new,...)
        if (is.list(attr(x,'pattern'))) {
            Xa <- attr(x,'pattern')
            nms <- names(Xa)
            col <- c('black','darkgreen','grey','yellow','magenta','cyan',
                     'brown','white','green')
                                        #browser()
            for (i in (2:length(nms))) 
                contour(lon(Xa[[i]]),lat(Xa[[i]]),Xa[[i]],add=TRUE,col=col[i])
        } else if (sum(is.element(type,'contour'))>0)
            contour(lon(X),lat(X),X,add=TRUE,col="grey50")
    } else if (projection=="sphere")
        map2sphere(x=X,lonR=lonR,latR=latR,axiR=axiR,
                   type=type,gridlines=gridlines,verbose=verbose,
                   colbar=colbar,new=new,...)
    else if (projection=="np")
        map2sphere(X,lonR=lonR,latR=90,axiR=axiR,
                   type=type,gridlines=gridlines,
                   colbar=colbar,new=new,verbose=verbose,...)
    else if (projection=="sp")
        map2sphere(X,lonR=lonR,latR=-90,axiR=axiR,
                   type=type,gridlines=gridlines,
                   colbar=colbar,new=new,verbose=verbose,...)
    invisible(X)
}


map.field <- function(x,FUN='mean',it=NULL,is=NULL,new=FALSE,
                      projection="lonlat",
                      xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                      colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                          pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                      type=c("fill","contour"),gridlines=FALSE,
                      lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
                      na.rm=TRUE,...) {
    
    stopifnot(inherits(x,'field'))
    if (verbose) print('map.field')

    x <- subset(x,it=it,is=is)
                                        #print(length(x)); print(attr(x,'dimensions')[1:2])
    projection <- tolower(projection)
    if (FUN=='trend') FUN <- 'trend.coef'

    if (!is.null(xlim)) {
        if (xlim[1] < 0) x <- g2dl(x,greenwich=FALSE)
    }
                                        #str(X)
    X <- coredata(x)
    

    ## If one time slice, then map this time slice
    if (dim(X)[1]==1) X <- coredata(x[1,]) else
    if (is.null(X)) X <- coredata(X) else if (inherits(X,"matrix")) {
        ## If several time slices, map the required statistics
        good <- apply(coredata(x),2,nv) > 1
        X <- rep(NA,length(good))
        xx <- coredata(x[,good])
        X[good] <- apply(xx,2,FUN=FUN,na.rm=na.rm)
    }

    ## if zlim is specified, then mask data outside this range
    if (!is.null(zlim)) {
        d <- dim(X)
        mask <- (X < min(zlim)) | (X > max(zlim))
        rng <- range(X,na.rm=TRUE)
        X[mask] <- NA
        dim(X) <- d
        if (verbose) print(paste('zlim=',zlim[1],zlim[2],
                                 '  sum(mask)=',print(sum(mask)),
                                 '  range(X)=',rng[1],rng[2]))
    }
    
                                        #print(length(X))
    attr(X,'longitude') <- attr(x,'longitude')
    attr(X,'latitude') <- attr(x,'latitude')
    attr(X,'variable') <- attr(x,'variable')[1]
                                        #  if (attr(x,'unit')=="deg C") attr(X,'unit') <- expression(degree*C) else
    unit <- attr(x,'unit')[1]
    if ( (is.na(unit) | is.null(unit)) ) unit <- " "
    if ((unit=='degree Celsius') | (unit=='deg C') | (unit=='degC'))
        unit <- 'degree*C'

    attr(X,'unit') <- unit
    attr(X,'source') <- attr(x,'source')
    attr(X,'time') <- range(index(x))
    attr(X,'method') <- FUN
    attr(X,'timescale') <- class(x)[2]
                                        #print(length(X)); print(attr(x,'dimensions'))
    dim(X) <- attr(x,'dimensions')[1:2]
                                        #class(X) <- class(x)
                                        #str(X)
    
    if (projection=="lonlat") lonlatprojection(x=X,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            colbar=colbar,type=type,new=new,
            gridlines=gridlines,verbose=verbose,...) else
    if (projection=="sphere") map2sphere(x=X,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            lonR=lonR,latR=latR,axiR=axiR,
            type=type,gridlines=gridlines,
            colbar=colbar,new=new,verbose=verbose,...) else
    if (projection=="np") map2sphere(X,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            lonR=lonR,latR=90,axiR=axiR,
            type=type,gridlines=gridlines,
            colbar=colbar,new=new,verbose=verbose,...) else
    if (projection=="sp") map2sphere(X,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            lonR=lonR,latR=-90,axiR=axiR,
            type=type,gridlines=gridlines,
            colbar=colbar,new=new,verbose=verbose,...)
    invisible(X)
}


map.corfield <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                         xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                         colbar= list(pal=NULL,rev=FALSE,n=NULL,
                             breaks=seq(-1,1,by=0.05),pos=0.05,show=TRUE,
                             type="p",cex=2,h=0.6,v=1),
                         type=c("fill","contour"),gridlines=FALSE,
                         lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,...) {
    
    if (verbose) print("map.corfield")
    stopifnot(inherits(x,'corfield'))
    x <- subset(x,it=it,is=is,verbose=verbose)
    projection <- tolower(projection)
    dim(x) <- attr(x,'dimensions')[1:2]
    ## if (!is.null(colbar)) colbar$pal <- varid(x)[2] ## AM 08-07-2015 comment 
    attr(x,'variable') <- paste(varid(x),collapse='/')
    
    #if (length(unit(x))>1) attr(x,'unit') <- paste(unit(x),collapse='/')
    attr(x,'unit') <- unit(x)[1]
    
    ## if zlim is specified, then mask data outside this range
    if (!is.null(zlim)) {
        if (verbose) print(zlim)
        d <- dim(x)
        mask <- (x < min(zlim,na.rm=TRUE)) | (x > max(zlim,na.rm=TRUE))
        x[mask] <- NA
        dim(x) <- d
        if (verbose) {print(zlim); print(dim(x)); print(sum(mask))}
    }

    if (verbose) {print(projection); print(dim(x))}
    
    if (projection=="lonlat") lonlatprojection(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            colbar=colbar,type=type,new=new,verbose=verbose,gridlines=gridlines,...) else
    if (projection=="sphere") map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            lonR=lonR,latR=latR,axiR=axiR,type=type,gridlines=gridlines,
            colbar=colbar,new=new,verbose=verbose,...) else
    if (projection=="np") map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            lonR=lonR,latR=90,axiR=axiR,type=type,gridlines=gridlines,
            colbar=colbar,new=new,verbose=verbose,...) else
    if (projection=="sp") map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            lonR=lonR,latR=-90,axiR=axiR,type=type,gridlines=gridlines,
            colbar=colbar,new=new,verbose=verbose,...)

    if (!is.null(attr(x,'x.longitude')) & !is.null(attr(x,'x.latitude')))
        points(attr(x,'x.longitude'),attr(x,'x.latitude'),lwd=2,cex=1.2)
    invisible(x)
}


map.trend <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                      xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                      colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                          pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                      type=c("fill","contour"),gridlines=FALSE,
                      lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,...) {
    if (verbose) print('map.trend')
    stopifnot(inherits(x,'field'),inherits(x,'trend'))
    x <- subset(x,it=it,is=is)
    projection <- tolower(projection)
    X <- attr(x,'pattern')

    ## if zlim is specified, then mask data outside this range
    if (!is.null(zlim)) {
        d <- dim(x)
        mask <- (x < min(zlim)) | (x > max(zlim))
        x[mask] <- NA
        dim(x) <- d
        if (verbose) {print(zlim); print(dim(x)); print(sum(mask))}
    } 
    attr(X,'longitude') <- attr(x,'longitude')
    attr(X,'latitude') <- attr(x,'latitude')
    attr(X,'variable') <- paste(attr(x,'variable'),'trend')
    attr(X,'time') <- range(index(x))
    attr(X,'unit') <- paste('d',attr(x,'unit'),'/decade')
    attr(X,'source') <- attr(x,'source')
    dim(X) <- attr(x,'dimension')[1:2]
                                        #str(X)

    if (projection=="lonlat") lonlatprojection(x=x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            colbar=colbar,type=type,new=new,
            verbose=verbose,
            gridlines=gridlines,...) else
    if (projection=="sphere") map2sphere(x=x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            lonR=lonR,latR=latR,axiR=axiR,
            type=type,gridlines=gridlines,
            colbar=colbar,new=new,verbose=verbose,...) else
    if (projection=="np") map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            lonR=lonR,latR=90,axiR=axiR,
            type=type,gridlines=gridlines,
            colbar=colbar,new=new,verbose=verbose,...) else
    if (projection=="sp") map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
            lonR=lonR,latR=-90,axiR=axiR,
            type=type,gridlines=gridlines,
            colbar=colbar,new=new,verbose=verbose,...)
    invisible(X)
}





map.pca <- function(x,it=NULL,is=NULL,pattern=1,new=FALSE,projection="lonlat",
                    xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
                    colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                        pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                    type=c("fill","contour"),gridlines=FALSE,
                    verbose=FALSE,...) {
    ##
                                        #args <- list(...)
                                        #print(args)
    X <- rbind(attr(x,'pattern')[,pattern],attr(x,'pattern')[,pattern])
                                        #print(dim(X))
                                        #str(x)
    X <- attrcp(x,X)

    ## if zlim is specified, then mask data outside this range
    if (!is.null(zlim)) {
        d <- dim(X)
        mask <- (X < min(zlim)) | (X > max(zlim))
        X[mask] <- NA
        dim(X) <- d
        if (verbose) {print(zlim); print(dim(X)); print(sum(mask))}
    }    
    attr(X,'longitude') <- lon(x)
    attr(X,'latitude') <- lat(x)
    class(X) <- 'station'
    ##if (is.null(colbar$col) | is.null(colbar)) {
    ##  colbar$col <- colscal(30,col=varid(x))
    ##}
    
    map.station(X,new=new,FUN="mean",
                colbar=colbar,
                xlim=xlim,ylim=ylim,zlim=zlim,verbose=verbose,...)
}

map.mvr <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                    xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
                    colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                        pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                    type=c("fill","contour"),gridlines=FALSE,
                    verbose=FALSE,...) {
    x <- subset(x,it=it,is=is)
    map.field(x,new=new,FUN="mean",
              colbar=colbar,
              cex=cex,xlim=xlim,ylim=ylim,verbose=verbose,...)
    
}

map.cca <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                    xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
                    colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                        pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                    type=c("fill","contour"),gridlines=FALSE,
                    verbose=FALSE,...) {
    ##print('map.cca')
    ##x <- subset(x,it=it,is=is)
    ## browser()
    ## For plotting, keep the same kind of object, but replace the patterns in
    ## the eof/pca with the CCA patterns
    Y <- x$Y
    ##print(dim(attr(Y,'pattern'))); print(dim(U))
    ##attr(Y,'pattern') <- U
    U <- x$B.m
    dim(U) <- c(dim(attr(Y,'pattern'))[-length(dim(attr(Y,'pattern')))],
                length(x$i.eofs))
    attr(Y,'pattern') <- U
    attr(Y,'eigenvalues') <- rep(1,length(x$i.eofs))
    attr(Y,'time') <- range(index(x))
    X <- x$X
    ##print(dim(attr(X,'pattern'))); print(dim(V))
    ##attr(X,'pattern') <- V
    V <- x$A.m
    dim(V) <- c(dim(attr(X,'pattern'))[-length(dim(attr(X,'pattern')))],
                length(x$i.eofs))
    attr(X,'pattern') <- V
    attr(X,'eigenvalues') <- rep(1,length(x$i.eofs))
    attr(X,'time') <- range(index(x))

    ## REB removed '...' in the two following map calls.
    ##  map(Y,icca,xlim=xlim,ylim=ylim,type=type,
    ##      projection=projection,lonR=lonR,latR=latR,axiR=axiR,
    ##      gridlines=gridlines,col=col,breaks=breaks,FUN='mean')
    ##  dev.new()
    ##  map(X,icca,xlim=xlim,ylim=ylim,type=type,
    ##      projection=projection,lonR=lonR,latR=latR,axiR=axiR,
    ##      gridlines=gridlines,col=col,breaks=breaks,FUN='mean')
    ##  print('Need to fix breaks and map.station')

    ##  col <- rgb( c(rep(0,15),1-sqrt(seq(0,1,length=15))),
    ##              abs(sin(seq(0,pi,length=30))),
    ##              c(sqrt(seq(0,1,length=15)),rep(1,15)) )
    ##  col <- colscal(30,col=varid(x))
    ##  if (is.precip(X)) col.x <- rev(col) else
    ##                    col.x <- col
    ##  if (is.precip(Y)) col.y <- rev(col) else
    ##                    col.y <- col
    ## REB: removed col=col.y,bg=col.y

    if (sum(is.element(type,'map'))>0)
        par(fig=c(0,0.5,0.5,1)) ## mar=c(0.05,.05,0.05,0.05),
    else 
        par(fig=c(0,0.5,0.5,1),mar=c(0.2,.2,0.2,0.2))

    ##colbar <- list(col=NULL, breaks=NULL, type="r",cex=2, h=0.6, v=1)
    
    map(Y,pattern=icca,xlim=xlim,ylim=ylim,type=type,cex=cex,
        projection=projection,lonR=lonR,latR=latR,axiR=axiR,
        gridlines=gridlines,FUN='mean',verbose=verbose,
        colbar=colbar,showall=FALSE,new=FALSE)
    ## browser()
    if (sum(is.element(type,'ts'))>0)
        par(fig=c(0,1,0.5,1),new=TRUE) else
    par(fig=c(0.5,1,0.5,1),new=TRUE) ## mar=c(0,0,0,0),
    map(X,pattern=icca,xlim=xlim,ylim=ylim,type=type,cex=cex,
        projection=projection,lonR=lonR,latR=latR,axiR=axiR,
        gridlines=gridlines,FUN='mean',verbose=verbose,
        colbar=colbar,showall=FALSE,new=FALSE)
    
    invisible(list(U=U,V=V))
}


                                        # Produce a KMZ-file to show the data in GoogleEarth.
map.googleearth <- function(x) {
}


lonlatprojection <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                             xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                             colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                 pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                             type=c("fill","contour"),gridlines=FALSE,
                             verbose=FALSE,geography=TRUE,fancy=FALSE,
                             main=NA,...) {

    if (verbose) print('lonlatprojection')
    colid <- 't2m'; if (is.precip(x)) colid <- 'precip'
    colorbar <- !is.null(colbar)
    #print(formals(...))
    ## If only a few items are provided in colbar - hen set the rest to the default
    ## browser()
    if (!is.null(colbar)) {
        colbar <- colbar.ini(x,FUN=NULL,colbar=colbar,verbose=verbose)
    } else {
        if (verbose) print('colbar=NULL - set col etc')
        colbar$n <- 25
        colbar$breaks <- pretty(c(x),n=colbar$n)
        if (verbose) print(colbar$breaks)
        if (verbose) print(varid(x))
        colbar$col <- colscal(n=length(colbar$breaks)-1,col=colid)
        ##  if ( (tolower(variable)=='precip') | (tolower(variable)=='tp') )
        if (colid=='precip') col <- rev(col)
        colbar$show <- TRUE
        colbar$pos <- 0.05
    }
 
    ##    par0 <- par()                             # REB 2015-06-25 these lines open an
    ##    fig0 <- par()$fig                         # unused window.
    fig0 <- c(0,1,0,1)                        # REB 2015-06-25
    
    ##    if (!is.null(colbar$pal) & (!is.null(colbar$n) | !is.null(colbar$breaks))) {
    ##        ##colbar$breaks <- pretty(y,n=length(colbar$col))
    ##        ##colbar$n <- length(colbar$breaks) + 1
    ##        if (is.null(colbar$breaks) & !is.null(colbar$n)) {
    ##            colbar$breaks <- pretty(x,n=colbar$n)
    ##            colbar$n <- length(colbar$breaks)-1
    ##        } else if (!is.null(colbar$breaks) & is.null(colbar$n))
    ##            colbar$n <- length(colbar$breaks)-1
    ##        else if (n != (length(colbar$breaks) -1))
    ##            stop('The length of breaks must equal (n-1)')# default
    ##        ##
    ##        if (verbose) print(paste("n=",colbar$n))
    ##        if (verbose) print(paste("breaks",colbar$breaks))
    ##        if (verbose) print(paste("length(breaks) =",length(colbar$breaks)))
    ##        colbar$col <- colscal(n=colbar$n,col=colbar$pal,rev=colbar$rev)
    ##        if (verbose) print(paste("length(col) =",length(colbar$col)))
    ##    }
    ##
    ##    if (!is.null(colbar$col)) col <- colbar$col else col <- NULL
    ##    if (!is.null(colbar$breaks)) breaks <- colbar$breaks else breaks <- NULL
    ##   
    #browser()
    if (colbar$show) { ## AM 14-07-2015
        ##        fig0[3] <- par0$fig[3] + (par0$fig[4]-par0$fig[3])/200##0.05
        fig0[3] <- fig0[3] + colbar$pos ## (fig0[4]-fig0[3])/200##0.05   # REB 2015-06-25
    } else 
        fig0 <- fig0                                       # REB 2015-06-25
    ##        fig0 <- par0$fig
    ##     par(fig=fig0)                                        ## REB 2015-06-25 opens extra window
    data("geoborders",envir=environment())
    if(sum(is.finite(x))==0) stop('No valid data')
    ## To deal with grid-conventions going from north-to-south or east-to-west:
    srtx <- order(lon(x)); lon <- lon(x)[srtx]
    srty <- order(lat(x)); lat <- lat(x)[srty]
    if (verbose) print('meta-stuff')
    unit <- unit(x); variable <- varid(x); varid <- varid(x); isprecip <- is.precip(x)
    if ( (unit=="degC") | (unit=="deg C") | (unit=="degree C") )
        unit <- "degree*C"
    if (unit=="%") unit <- "'%'"
    if ( (tolower(variable)=="t(2m)") | (tolower(variable)=="t2m") |
        (tolower(variable)=="2t") )
        variable <- "T[2*m]"
    varlabel=eval(parse(text=paste('expression(',
                variable," *(",unit,"))",sep="")))
    sub <- attr(x,'source')
    if (sum(is.element(type,'fill'))==0) colbar <- NULL
    
    if (verbose) print('time')
    if (!is.null(attr(x,'timescale'))) {
        if (verbose) print(attr(x,'timescale'))
        timescale <- attr(x,'timescale')
        if (timescale == 'annual') {
            t1 <- year(attr(x,'time'))[1]
            t2 <- year(attr(x,'time'))[2]
        } else
            if (sum(is.element(c('month','season'),timescale))>0) {
                t1 <- paste(year(attr(x,'time'))[1],month(attr(x,'time'))[1])
                t2 <- paste(year(attr(x,'time'))[2],month(attr(x,'time'))[2])
            } else {
                t1 <- attr(x,'time')[1]  
                t2 <- attr(x,'time')[2]
            }
        period <- paste('[',t1,', ',t2,']',sep='')
    } else period <- NULL
    if (verbose) print(paste('period:',period))
    method <- attr(x,'method')
    if (verbose) print(c(dim(x),length(srtx),length(srty)))
                                        #browser()
    x <- x[srtx,srty]
    if (verbose) {print(xlim); str(x)}
    if (!is.null(xlim)) {
        outside <- (lon < xlim[1]) | (lon > xlim[2])
        x[outside,] <- NA
    } else xlim <- range(lon)
    
    if (!is.null(ylim)) {
        outside <- (lat < ylim[1]) | (lat > ylim[2])
        x[,outside] <- NA
    } else ylim=range(lat)
    
    ##print(c(length(breaks),length(col)))
    ##if (is.Date(type))
    
    ##if ( (par()$mfcol[1]> 1) | (par()$mfcol[2]> 1) ) new <- FALSE
    ## browser()
    
    if (new) {
        par(fig=fig0) 
        dev.new()
        par(bty="n",xaxt="n",yaxt="n",xpd=FALSE)
        ## fig=fig0,mar=c(2,1,1,1)) # c(0.05,0.95,0.13,0.95),mar=rep(1,4)
        ##    par(fig=fig0,mar=c(2.5,2,2,2),bty="n") # c(0.05,0.95,0.13,0.95),mar=rep(1,4)
        ##    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
        ##        fig=c(0.05,0.95,0.12,0.95))
    } else {
        par(bty="n",xaxt="n",yaxt="n",xpd=FALSE)
        ## par(bty="n",xaxt="n",yaxt="n",xpd=FALSE)
    }
    
    plot(range(lon),range(lat),type="n",xlab="",ylab="", # REB 10.03
         xlim=xlim,ylim=ylim,main=main, # to sumerimpose.
         xaxt="n",yaxt="n") # AM 17.06.2015
    ##par0 <- par()

    if (sum(is.element(tolower(type),'fill'))>0)   
        image(lon,lat,x,xlab="",ylab="",add=TRUE,
              col=colbar$col,breaks=colbar$breaks,xlim=xlim,ylim=ylim,...)
    
    if (geography) {
        lines(geoborders$x,geoborders$y,col="darkblue")
        lines(attr(geoborders,'borders')$x,attr(geoborders,'borders')$y,col="pink")
        lines(geoborders$x+360,geoborders$y,col="darkblue")
    }
    if (sum(is.element(tolower(type),'contour'))>0)
        contour(lon,lat,x,lwd=1,col="grey70",add=TRUE)
    if (gridlines) grid()
    par(xpd=FALSE)
    dlat <- diff(range(lat))/60
                                        #print(dlat)
    text(lon[1],lat[length(lat)] + dlat,varlabel,pos=4,font=2)
    text(lon[1],lat[1] - dlat,sub,col="grey30",pos=4,cex=0.7)

    if (!is.null(period))
        text(lon[length(lon)],lat[length(lat)] + dlat,period,pos=2,cex=0.7,col="grey30")
    if (!is.null(method))
        text(lon[length(lon)],lat[1] - 0.5*dlat,method,col="grey30",pos=2,cex=0.7)
    
    if (!is.null(colbar)) {
        if (verbose) print('Add colourbar')

        ## =======
        ##  if (!is.null(colbar)) {
        ##    if (verbose) print('Add colourbar')
        ##    par(xaxt="s",yaxt="s")
        ## 1072ef5b4e555d6484178b0115e5d62be3dbd386
        
        ## Old    
        ##    par(xaxt="s",fig=c(0.05,0.95,0.01,1))
        ##    breaks <- round(seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(col)),1)
        ##    image.plot(horizontal=TRUE,legend.only=TRUE,zlim=range(x,na.rm=TRUE),
        ##               lab.breaks=breaks,col=col,axis.args=list(cex.axis=0.8),
        ##               border=FALSE)
        ##
        ##    par(fig=par0$fig,mar=par0$mar,new=TRUE,xaxt="n")
        ##    plot(range(lon),range(lat),type="n",xlab="",ylab="", # REB 10.03
        ##         xlim=xlim,ylim=ylim)                # to sumerimpose.
        ##
        ## Adopt from map.station
        ##
        ##if (is.null(colbar$col)) colbar$col <- colscal(n=n,varid)
        ##if (is.null(colbar$breaks)) colbar$breaks <- pretty(x,n=length(colbar$col))
        
        ##if (isprecip) colbar$col <- rev(colbar$col)
        ##browser()
        par(xaxt="s",yaxt="s",las=1,col.axis='grey',col.lab='grey',
            cex.lab=0.7,cex.axis=0.7)
        axis(2,at=pretty(lat(x)),col='grey')
        axis(3,at=pretty(lon(x)),col='grey')
        grid()

        par(col.axis='black',col.lab='black',
            cex.lab=0.5,cex.axis=0.5)
        
        if (!is.null(colbar))
            if (colbar$show)
                if (fancy)
                    col.bar(colbar$breaks,horiz=TRUE,pch=21,v=1,h=1,
                            col=colbar$col, cex=2,cex.lab=colbar$cex.lab,
                            type=type,verbose=FALSE,vl=1,border=FALSE)
                else {
                                        #par(fig=par0$fig)
                    ##image.plot(lab.breaks=colbar$breaks,horizontal = TRUE,
                    ##           legend.only = T, zlim = range(colbar$breaks),
                    ##           col = colbar$col, legend.width = 1,
                    ##           axis.args = list(cex.axis = 0.8), border = FALSE)
                    image.plot(lab.breaks=colbar$breaks,horizontal = TRUE,
                               legend.only = T, zlim = range(colbar$breaks),
                               col = colbar$col, legend.width = 1,
                               axis.args = list(cex.axis = 0.8,
                                   xaxp=c(range(colbar$breaks),n=colbar$n)),
                               border = FALSE,...)
                }
    }

                                        #par(fig=fig0)

    par(col.axis='black',col.lab='black',cex.lab=1,cex.axis=1,
        xaxt="s",yaxt="s")
    result <- list(x=lon,y=lat,z=x,breaks=colbar$breaks)
                                        #par(fig=par0$fig)
    invisible(result)
}
