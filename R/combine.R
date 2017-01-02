                                        # R.E. Benestad
                                        # Combines a group of station objects to form a set of time series
                                        # Combines field objects; either in time (concatinates) or in space (mixing fields)
                                        # Can be used to combine EOFS

zeros <- function(x) (sum(is.infinite(1/x)) > 0)

rbind.field <- function(...) {
    print('note: the results inherit the metadata from the first argument')
    x <- list(...)
    y <- rbind.zoo(...)
    y <- attrcp(x[[1]],y)
    attr(y,'dimensions') <- c(attr(x[[1]],'dimensions')[1:2],length(index(y)))
    attr(y,'history') <- history.stamp(x[[1]])
    class(y) <- class(x[[1]])
    return(y)
}

cbind.field <- function(...) {
    x <- list(...)
    y <- cbind.zoo(...)
    print('unfinished - need to combine the attributes from all inputs')
    y <- attrcp(x[[1]],y)
    attr(y,'history') <- history.stamp(x[[1]])
    class(y) <- class(x[[1]])
    return(y)
}


combine <- function(x,y,...)
    UseMethod("combine")

combine.default <- function(x=NULL,y=NULL,all=FALSE,orig.format=TRUE) {
                                        # If either of the arguments is NULL, then return the x - useful for looping
    stopifnot(!missing(x))
                                        #print("combine.default")
    if (is.null(x)) return(y)
    if (is.null(y)) return(x)
    
                                        #print(class(x)); print(class(y)); print(summary(x))
                                        #print("dim(y)="); print(dim(y))
    
    if ( !zeros(c(inherits(x,c("station","zoo"),which=TRUE),
                  inherits(y,c("station","zoo"),which=TRUE)) ) )
        X <- combine.station(x,y) else
    if ( !zeros( c(inherits(x,c("station","zoo"),which=TRUE)
                   ,inherits(y,c("eof","zoo"),which=TRUE)) ) |
        !zeros( c(inherits(y,c("station","zoo"),which=TRUE),
                  inherits(x,c("eof","zoo"),which=TRUE)) ) )
        X <- combine.station.eof(x,y,all=all,
                                 orig.format=orig.format) else
    if ( !zeros( c(inherits(x,c("station","zoo"),which=TRUE)
                   ,inherits(y,c("field","zoo"),which=TRUE)) ) |
        !zeros( c(inherits(y,c("station","zoo"),which=TRUE),
                  inherits(x,c("field","zoo"),which=TRUE)) ) )
        X <- combine.station.field(x,y,all=all,
                                   orig.format=orig.format) else { 
                                       print("combine.default - don't know what to do :-(")
                                       Z <- NULL
                                   }
    attr(X,'history') <- history.stamp(x)
    invisible(X)
}

                                        # combine.station can be used to either combine a group of stations into
                                        # one data object or timerseries from one stations with different monthly
                                        # values into one time series with all months
combine.station <- function(...,all=TRUE) {
                                        #print("combine.station")
    cl <- as.list(match.call())
                                        #str(cl)
    args <- list(...)
    n <- length(args)
    stid <- rep(NA,n)
    allstations <- TRUE
    for (i in 1:n) {
        Z <- args[[i]]
                                        #print(class(Z))
        if (inherits(Z,'station'))
            stid[i] <- attr(Z,'station_id')[i]
        else allstations <- FALSE
    }
                                        #print(allstations)

    if (allstations) {
        ns <- length(table(stid))
                                        # If only one site, then combine the months into one series, otherwise
                                        # group the stations into one multivariate object
        if (ns ==1) X <- combine.station.month(...) else
        X <- combine.stations(...,all=all)
    } else {
        X <- combine.default(...,all=all)
    }
    attr(X,'history') <- history.stamp(X)
    class(X) <- class(Z)
    invisible(X)
}

                                        # combine.station.month is used to
                                        # This is causing the problem in combine.ds REB
combine.station.month <- function(...) {
                                        #print("combine.station.month")
    cl <- as.list(match.call())
                                        #str(cl)
    args <- list(...)
                                        #print(summary(args)); print(length(args))
    z <- args[[1]]
                                        #plot.zoo(z); dev.new()
    Z <- merge.zoo(...,all=TRUE)
                                        #plot(Z);str(Z)
    X <- zoo(rowMeans(coredata(Z),na.rm=TRUE),order.by=as.Date(index(Z)))
                                        #print(summary(X))
                                        #print(names(attributes(z)))
    X <- attrcp(z,X,ignore=c('mean','calibration_data','fitted_values',
                        'original_data','aspect'))
    class(X) <- class(z)
    attr(X,'history') <- history.stamp(X)
    invisible(X)
}

combine.zoo <- function(...) {
                                        #print("combine.zoo")
    Z <- merge.zoo(...,all=TRUE)
                                        #plot(Z);str(Z)
    invisible(Z)
}

                                        # combine.stations is used to combine a group of stations into one object

combine.stations <- function(...,all=TRUE) {
                                        #print("combine.stations")
                                        # If either of the arguments is NULL, then return the x -
                                        # useful for looping

                                        # REB: rewritten to work with merge and long lists of stations
    
    cl <- as.list(match.call())
                                        #str(cl)
    args <- list(...)
                                        #str(args)
                                        
    X <- merge.zoo(...,all=all)
                                        #plot(X)
                                        #str(X)
                                        #print(length(args))
    
    n <- length(args)
    cls <- class(args[[1]])
    
    loc <- NULL; cn <- loc; ID <- NULL; unit <- loc
    lon <- ID; lat <- ID; alt <- ID; param <- loc; lname <- loc
    src <- loc; qlty <- loc; url <- loc; ref <- loc
    info <- loc; ele <- asp <- ID
    for (i in 1:n) {
        Z <- args[[i]]
                                        #attr <- softattr(Z)
        loc <- c(loc,attr(Z,'location'))
        cn  <- c(cn,attr(Z,'country'))
        ID <- c(ID,attr(Z,'station_id'))
        unit <- c(unit,attr(Z,'unit'))
        lon <- c(lon,attr(Z,'longitude'))
        lat <- c(lat,attr(Z,'latitude'))
        alt <- c(alt,attr(Z,'altitude'))
        param <- c(param,attr(Z,'variable'))
        lname <- c(lname,attr(Z,'longname'))
        src <- c(src,attr(Z,'source'))
        qlty <- c(qlty,attr(Z,'quality'))
        url <- c(url,attr(Z,'URL'))
        ref <- c(ref,attr(Z,'reference'))
        info <- c(info,attr(Z,'info'))
        ele <- c(ele,attr(Z,'element'))
        asp <- c(asp,attr(Z,'aspect'))
    }

    if (dim(X)[2]==length(loc)) colnames(X) <- loc
    attr(X,'location') <- loc
    attr(X,'country') <- cn
    attr(X,'station_id') <- ID
    attr(X,'longitude') <- lon
    attr(X,'latitude') <- lat
    attr(X,'altitude') <- alt
    attr(X,'variable') <- param
    attr(X,'longname') <- lname
    attr(X,'unit') <- unit
    attr(X,'aspect') <- asp      ## AM
    attr(X,'source') <- src
    attr(X,'element') <- ele
    attr(X,'quality') <- qlty
    attr(X,'URL') <- url
    attr(X,'history') <- history.stamp(Z)
                                        #attr(X,'date-stamp') <- date()
    attr(X,'reference') <- ref
    attr(X,'info') <- info
    class(X) <- cls
    invisible(X)
}

combine.ds <- function(...,all=TRUE) {
    ##print("combine.ds")
    cl <- as.list(match.call())
    ##str(cl)
    args <- list(...)
                                        #str(args)
    z <- args[[1]]
    n <- length(args)
                                        #print(n); print(summary(z)); print(names(attributes(z[[1]])))
                                        #print(class(z))
    if (n > 1) {
                                        #  Arguments consisting of several objects to be combined
                                        #print("Several stations")
        if (inherits(z,'pca'))
            X <- combine.ds.pca(...,all=all) else
        X <- combine.station(...,all=TRUE)
    } else if (n==1) {
                                        #print("One station")
                                        #print(class(z)); print(attr(z[[1]],'type')); print(summary(z))
        if ( (attr(z[[1]],'type')=="downscaled results") & is.list(z) ) {
                                        # The output of DS.t2m.month.field: one list with several months
                                        # The main results looks like a station object:
                                        #    X <- as.station.ds(z)
            X <- as.station.list(z)
            #X <- rbind(z[[1]],z[[2]],z[[3]],z[[4]]) # AM 27.11.2014 quick fix must be included in as.station.list ..
                                        #print("HERE"); print(names(X))
                                        #plot(X,col="red")
           
                                        # The original data & fitted values:
                                        #print("combine the training data ...")
            m <- length(z)
            for (i in 1:m) {
                                        #model <- attr(z[[i]],'model')
                y <- attr(z[[i]],'fitted_values')
                x <- attr(z[[i]],'calibration_data')
                x0 <- attr(z[[i]],'original_data')
                xval <- crossval(z[[i]])
                                        #print(summary(y - y0))
                attr(y,'station_id') <- attr(X,'station_id')
                attr(x,'station_id') <- attr(X,'station_id')
                attr(x0,'station_id') <- attr(X,'station_id')
                eval(parse(text=paste("y.",i," <- y",sep="")))
                eval(parse(text=paste("x.",i," <- x",sep="")))
                eval(parse(text=paste("x0.",i," <- x0",sep="")))
                eval(parse(text=paste("xv.",i," <-  zoo(xval)",sep="")))
                eof <- attr(z[[i]],'eof')  # REB 12.02.2014
                if (i==1) {
                    argsx <- 'x.1'
                    argsy <- 'y.1'
                    argsx0 <- 'x0.1'
                    argsxY <- 'xv.1'
                    pattern <- attr(z[[1]],'pattern')
                    d <- dim(pattern)
                    dim(pattern) <- c(1,d[1]*d[2])
                                        # REB - also capture the diagnostics - from the commone EOF
                    if (!is.null(attr(eof,'diagnose')))
                        diag <- list(s.1=attr(eof,'diagnose'))
                } else {
                    argsx <- paste(argsx,",x.",i,sep="")
                    argsy <- paste(argsy,",y.",i,sep="")
                    argsx0 <- paste(argsx0,",x0.",i,sep="") ## AM 27.11.2014 argsx replaced by argsx0
                    argsxY <- paste(argsxY,",xv.",i,sep="") ## AM 27.11.2014 argsx replaced by argsx0
                    patt <- attr(z[[i]],'pattern')
                    dim(patt) <- c(1,d[1]*d[2])
                    pattern <- rbind(pattern,patt)
                                        # REB - also capture the diagnostics - from the commone EOF
                    if (!is.null(attr(eof,'diagnose')))
                        eval(parse(text=paste("diag$s.",i," <- attr(eof,'diagnose')",sep="")))
                }
            }
                                        #print(argsx)
                                        #str(x.1); print(class(x.1))
                                        #print("---")
            cline1 <- parse(text=paste("X0 <- rbind(",argsx,")")) ## cline1 <- parse(text=paste("X0 <- merge(",argsx,",all=TRUE)"))
                                        #print(cline1)
            eval(cline1)
                                        #lines(X0)    

                                        #print("here"); print(argsy)
            cline2 <- parse(text=paste("Y <- rbind(",argsy,")")) ## AM cline2 <- parse(text=paste("Y <- merge(",argsy,",all=TRUE)"))
                                        #print(cline2)
            eval(cline2)
                                        #lines(Y,col="pink")

                                        #print("HERE")
            cline3 <- parse(text=paste("X0.0 <- rbind(",argsx0,")")) ## cline3 <- parse(text=paste("X0.0 <- merge(",argsx0,",all=TRUE)"))
                                        #print(cline3)
            eval(cline3)
                                        #lines(Y,col="pink")

                                        # Cross-validation
                                        #str(x.1); print(class(xv.1))
            cline4 <- parse(text=paste("xY <- rbind(",argsxY,")")) ## cline4 <- parse(text=paste("xY <- merge(",argsxY,",all=TRUE)"))
            eval(cline4)
                                                    
            attr(X,'calibration_data') <- X0
            attr(X,'fitted_values') <- Y
            attr(X,'evaluation') <- xY
            attr(X,'original_data') <- X0.0
            attr(X,'diagnose') <- diag
                                        #print(dim(pattern))
            dim(pattern) <- c(m,d[1],d[2])
            attr(pattern,'dimensions') <- c('month','longitude','latitude')
            attr(pattern,'month') <- month.abb
#            attr(pattern,'longitude') <- attr(z[[1]],'longitude')
#            attr(pattern,'latitude') <- attr(z[[1]],'latidude')
            attr(pattern,'longitude') <- lon(attr(z[[1]],'pattern'))
            attr(pattern,'latitude') <- lat(attr(z[[1]],'pattern'))
            
            attr(X,'pattern') <- pattern
                                        #print(class(X))
        } else X <- NULL 
    }

                                        #print(names(attributes(z)))
                                        #str(X)
                                        # Copy the attributes from z, but not those already used for station
                                        # objects, as those are taken care of in combine.station
                                        #print("Copy attributes")
    ## X <- attrcp(z[[1]],X,ignore=c('location','country','station_id','longitude',
    ##                  'latitude','altitude','variable','longname','unit','evaluation',
    ##                  'source','element','quality','URL','reference','info','names'))
                                        #print("...")
    attr(X,'type') <- 'downscaled results'
    attr(X,'history') <- history.stamp(X)
    attr(X,'old_class') <- class(z[[1]])
                                        #print("exit combine.ds")
                                        #print("HERE"); print(names(X))
    invisible(X)
}

combine.list <- function(...,all=TRUE) {
                                        #print("combine.list")
    args <- list(...)
    z <- args[[1]]
                                        #print(class(z[[1]]))
    if (inherits(z[[1]],'comb'))
        y <- combine.ds.comb(...,all=all) else
    if (inherits(z[[1]],'ds'))
        y <- combine.ds(...,all=all) else
    y <- NULL
    return(y)
}

combine.ds.comb <- function(...,all=TRUE) {
                                        #print("combine.ds.comb")
    cl <- as.list(match.call())
                                        #str(cl)
    args <- list(...)
                                        #print(args)  # [[1]][[1:4]]
                                        #print(summary(args))
    z <- args[[1]]
                                        #print(summary(z)); print(class(z)); print(length(z))
    k <- length(args)
    m <- length(...) 
                                        #print(k); print(m); print(summary(z)); print(names(attributes(z[[1]])))
                                        #print(class(z))

    if (k > 1) {
                                        #  Arguments consisting of several objects to be combined
                                        #print("Several stations")
        if (inherits(z,'pca'))
            X <- combine.ds(...,all=all)
    } else {
                                        #print("One station")
                                        #print(class(z))
        X <- combine.ds(...,all=all)
        n <- attr(z[[1]],'n.apps')
        
                                        # Combine the different appendixes for the different list objects:
                                        # appendix.1 in each list object is combined with corresponding
                                        # appendix.1 in the others and so on.
                                        #print("appendix.x")
        for (i in 1:n) {
            for (j in 1:m) {
                Z <- z[[j]]
                                        #print(summary(Z))
                                        #print(mean(attr(Z,paste('appendix.',i,sep=''))))
                x <- attr(Z,paste('appendix.',i,sep=''))
                                        #attr(x,'station_id') <- attr(Z,'station_id')
                                        #attr(x,'aspect') <- 'original'
                                        #class(x) <- c('station',class(Z))
                                        #print(summary(x))
                class(x) <- 'zoo'
                eval(parse(text=paste("x.",j," <- x",sep="")))
                if (j==1) {
                    argsx <- 'x.1'
                } else {
                    argsx <- paste(argsx,",x.",j,sep="")
                }
            }

                                        #print(argsx)
                                        #print("combine.station -> appendix.x")
                                        #plot(c(x.1,x.2,x.3,x.4)); dev.new()  
                                        #      cline1 <- parse(text=paste("attr(X,'appendix.",i,
                                        #                                 "') <- combine.zoo(",argsx,")",sep=""))
            cline1 <- parse(text=paste("attr(X,'appendix.",i,"') <- c(",argsx,")",sep=""))  
                                        #print(cline1)
            eval(cline1)
                                        #lines(X0)    
                                        #print("exit combine.ds.comb") 
        }
    }
    invisible(X)
}

combine.ds.station <- function(...,all=TRUE) {
                                        # Combine downscaled station records. Use combine.station for
                                        # combining the station values: either as a group of stations
                                        # or combine different months for one station into one series
    
                                        #print("combine.ds.station")
    cl <- as.list(match.call())
                                        #str(cl)
    args <- list(...)
                                        #str(args)
    X <- combine.station(...,all=all)
    attr(X,'type') <- 'ESD_result'
    attr(X,'history') <- history.stamp(X)
    
    invisible(X)
}

combine.ds.pca <- function(...,all=TRUE) {
                                        # Combine downscaled PCA: i.e. the different principal components
                                        # Assmue that the pattern is the same for all these
                                        #print("combine.ds.pca")
    cl <- as.list(match.call())
                                        #str(cl)
    args <- list(...)
                                        #str(args)
    X <- combine.ds.station(...,all=all)
    attr(X,'pattern') <- attr(args[[1]],'pattern')
    attr(X,'eigenvalues') <-  attr(args[[1]],'eigenvalues')
    attr(X,'mean') <-  attr(args[[1]],'mean')
    attr(X,'sum.eigenv') <- attr(args[[1]],'sum.eigenv')
    attr(X,'tot.var') <- attr(args[[1]],'tot.var')
    X <- attrcp(args[[1]],X)
    attr(X,'history') <- history.stamp(X)
    class(X) <- class(args[[1]])
    invisible(X)
}


combine.station.eof <- function(x,y,all=FALSE,orig.format=TRUE) {
                                        #print("combine.station.eof")
                                        # If either of the arguments is NULL, then return the x -
                                        # useful for looping
    
    if (is.null(x)) return(y)
    if (is.null(y)) return(x)

                                        #print("combine.station.eof")
                                        #print(dim(y)); print(dim(x))
    
                                        # Keep track of which is an eof object and which is a station record:
    if ( inherits(x,c("station")) & inherits(y,c("eof"))) {
        yy <- x
        x <- y
        y <- yy
    } 

                                        #str(y); str(x)
    clsx <- class(x)
    clsy <- class(y)
    index(y) <- as.Date(index(y))
                                        #print(class(x))
    index(x) <- as.Date(index(x))
                                        #print("combine: dates")
                                        #print(index(x)[1]); print(index(y)[1]); 
                                        # REB 20.08.2013: added colnames to x & y before merge to keep track of
                                        # y & x.

    dy <- dim(y)
    if (is.null(dy)) dy <- c(length(y),1)
    if (dy[2]>1) colnames(y) <- paste("y",1:dy[2],sep=".")
    dx <- dim(x)
    if (is.null(dx)) dx <- c(length(x),1)
    if (dx[2]>1) colnames(x) <- paste("x",1:dx[2],sep=".") 
    
    if (orig.format) {
        comb <- merge(x,y,all=all)
        vars <- tolower(names(comb))
        ys <- vars[grep('y',vars)]
        Xs <- vars[grep('x',vars)]
        ix <- is.element(vars,Xs)
        iy <- is.element(vars,ys)

        XX <- zoo(coredata(comb[,ix]),order.by=index(comb))
        yy <- zoo(comb[,iy],order.by=index(comb))

        XX <- attrcp(x,XX,ignore='names')
        yy <- attrcp(y,yy,ignore='names')

        
                                        # REB 29.04.2014
        
                                        #   nattr1 <- softattr(y)
                                        #   for (i in 1:length(nattr1))
                                        #     attr(yy,nattr1[i]) <- attr(y,nattr1[i])
                                        #   nattr2 <- softattr(x)
                                        #   for (i in 1:length(nattr2))
                                        #     attr(XX,nattr2[i]) <- attr(x,nattr2[i])
        
        clsy -> class(yy)
        clsx -> class(XX)
                                        
        if (!is.null(attr(yy,'standard.error'))) {
            sterr <- attr(yy,'standard.error')
            sterr <- matchdate(sterr,yy)
            attr(yy,'standard.error') <- coredata(sterr)
        }
        X <- list(y=yy,X=XX)
    } else {
                                        # REB 29.04.2014
                                        # This option introdices an additional index in the PCs and
                                        # and additional uniform pattern.
        comb <- merge(y,x,all=all)

        d <- dim(y)
        if (is.null(d)) {
            s <- sd(y,na.rm=TRUE)
            z <- ( coredata(y) - mean(y,na.rm=TRUE) )/s
        } else {
            s <- apply(y,2,sd,na.rm=TRUE)
            m <- apply(y,2,mean,na.rm=TRUE)
            z <- ( coredata(y) -  m )/s
        }
        coredata(y) <- z
        comb <- merge(y,x,all=all)
        
        X <- comb
                                        #    nattr1 <- softattr(x)
                                        #    for (i in 1:length(nattr1))
                                        #      attr(X,nattr1[i]) <- attr(x,nattr1[i])
        X <- attrcp(x,X,ignore='names')
        attr(X,'pattern') <- rbind(matrix(rep(s,d[1]*d[2]),d[1],d[2]),
                                   attr(X,'pattern'))
        attr(X,'X.attributes') <- attributes(y)
        class(X) <- class(x)
    }
                                        #print(summary(combined))
    attr(X,'history') <- history.stamp(X)
    invisible(X)
}


combine.station.field <- function(x,y,all=FALSE,orig.format=TRUE) {
                                        #print("combine.station.field")
                                        # If either of the arguments is NULL, then return the x -
                                        # useful for looping
    
    X <- combine.field.station(x=y,y=x,all=all,orig.format=orig.format)
    invisible(X) 
}


combine.field.station <- function(x,y,all=FALSE,
                                  orig.format=TRUE,verbose=FALSE) {
  if (verbose) print("combine.field.station")
    
                                        # If either of the arguments is NULL, then return the x -
                                        # useful for looping
    
    if (is.null(x)) return(y)
    if (is.null(y)) return(x)
    attr(x,'dimnames') <- NULL
                                        # Keep track of which is an eof object and which is a station record:
    swapped <- FALSE
    if ( inherits(x,c("station")) & inherits(y,c("field"))) {
        yy <- x
        x <- y
        y <- yy
        swapped <- TRUE
    } 
    print(swapped)
    clsx <- class(x)
    clsy <- class(y)
    index(y) <- as.Date(index(y))
  if (verbose) print(class(X))
    index(x) <- as.Date(index(x))
                                        #print("HERE...")
                                        # REB 20.08.2013: added colnames to x & y before merge to keep track of
                                        # y & x.
    if (length(dim(y))==2) 
        colnames(y) <- paste("y",1:dim(y)[2],sep=".") 
    if (length(dim(x))==2)
        colnames(x) <- paste("x",1:dim(x)[2],sep=".")
    comb <- merge(x,y,all=all) # AM 2014-12-02 An extra column is added here and the attribute dimensions could not be used anymore !!! so map will not work here ...
                                        if (verbose) print(summary(comb))

    if (orig.format) {
        vars <- names(comb)
        ys <- vars[grep('y',vars)]
        Xs <- vars[grep('x',vars)]
        ix <- is.element(vars,Xs)
        iy <- is.element(vars,ys)
                                        if (verbose) print(c(sum(ix),sum(iy)))
        
        XX <- zoo(coredata(comb[,ix]),order.by=index(comb))
        yy <- zoo(coredata(comb[,iy]),order.by=index(comb))
        names(XX) <- Xs
        names(yy) <- ys
        
                                        #    nattr1 <- softattr(x)
                                        #    for (i in 1:length(nattr1))
                                        #      attr(XX,nattr1[i]) <- attr(x,nattr1[i])
        XX <- attrcp(x,XX,ignore='names')
                                        #    nattr2 <- softattr(y)
                                        #    for (i in 1:length(nattr2))
                                        #      attr(yy,nattr2[i]) <- attr(y,nattr2[i])
        yy <- attrcp(y,yy,ignore='names')
       
        attr(XX,'dimensions') <- attr(x,'dimensions')
        ##attr(yy,'dimensions') <- attr(x,'dimensions')
        
                                        #    mostattributes(yy) <- attributes(x)
                                        #    mostattributes(XX) <- attributes(y)
        clsx -> class(XX)
        clsy -> class(yy)
        
                                        #print(clsx); print(clsy); print(dim(XX)); print(length(yy))
        if (swapped)
            X <- list(y=XX,X=yy)
        else
            X <- list(y=yy,X=XX)
    } else {   
        X <- comb
        X <- attrcp(x,X,ignore='names')
                                        #     nattr2 <- softattr(x)
                                        #     for (i in 1:length(nattr1))
                                        #       attr(X,nattr1[i]) <- attr(x,nattr1[i])
                                        #mostattributes(comb) <- attributes(x)
        attr(X,'X.attributes') <- attributes(y)
        class(X) <- c("comb",clsx)
    }
                                        #print(summary(combined))
    attr(X,'history') <- history.stamp(X)
    
    invisible(X)
}



combine.field <- function(x,y,all=FALSE,dimension="time",
                          approach="field",orig.format=TRUE,verbose=FALSE) {
 
  if (verbose) print(paste("combine.field",approach))
  attr(x,'dimnames') <- NULL
  attr(y,'dimnames') <- NULL
    if (inherits(y,'station')) {
        xy <- combine.field.station(x=x,y=y,orig.format=orig.format)
        return(xy)
    }
    stopifnot(!missing(x),
              inherits(x,'field'),inherits(y,'field'))
    clsx <- class(x); clsy <- class(y)
    hst <- c(attr(x,'history'),attr(y,'history'))
    if ( (unit(x)=="hPa") & (unit(y)=="Pa")) {
        if (verbose) print('Resetting unit of x: hPa -> Pa')
        coredata(x) <- 100*coredata(x)
        attr(x,'unit') <- 'Pa'
    }
    if (unit(x)=='deg C') attr(x,'unit') <- 'degC'
    if (unit(y)=='deg C') attr(y,'unit') <- 'degC'
    if ( (unit(y)=="hPa") & (unit(x)=="Pa")) {
        if (verbose) print('Resetting unit of y: hPa -> Pa')
        coredata(y) <- 100*coredata(y)
        attr(y,'unit') <- 'Pa'
    }
    if (unit(x) != unit(y))
        print(paste('Warning - different units:',unit(x),unit(y)))

    
    if (missing(y)) return(x)
    
    dimension <- tolower(dimension)
    approach <- tolower(approach)
    x <- sp2np(x)
    y <- sp2np(y)
                                        # Make sure that the longitude conventions are the same:
    if ( ( as.logical(attr(x,'greenwich')) &
          !as.logical(attr(y,'greenwich')) ) |
        ( !as.logical(attr(x,'greenwich')) &
         as.logical(attr(y,'greenwich')) ) )
        y <- g2dl(y,attr(x,'greenwich'))
    
    if (sum(is.element(dimension,"time"))) {
                                        # Combine the two gridded data sets along the time axis:
                                        # The y data set is regridded onto the grid of the former.

                                        # Organise the grid coordinates:
                                        #x1 <- attr(x,'longitude');  nx1 <- length(x1)
                                        #y1 <- attr(x,'latitude');   ny1 <- length(y1)
                                        #x2 <- attr(y,'longitude'); nx2 <- length(x2)
                                        #y2 <- attr(y,'latitude');  ny2 <- length(y2)
                                        #d1 <- attr(x,'dimensions')
                                        #d2 <- attr(y,'dimensions')
                                        #xy2 <- rep(x2,ny2); yx2 <- sort(rep(y2,nx2))

                                        # Work with the raw matrices rather than zoo class to avoid surrises;-)
                                        #    X <- coredata(x); Y <- coredata(y)
                                        #    Z <- matrix(rep(NA,d2[3]*d1[1]*d1[2]),d2[3],d1[1]*d1[2])

                                        # Only use the same region as the previous:
                                        # y <- subset(y,is=list(range(x1),range(y1)))
      if (verbose) print("combine.field after subset:")
                                        #Z <- regrid(y,is=list(x1,y1))
        Z <- regrid(y,is=x)
        
        ## KMP 2016-03-16 Issues with DSensemble because appendix.1 is empty.
        ## The problem originates here. Let's try changing maskna.
        iv <- function(x) return(sum(is.finite(x))) # KMP 2016-03-16
        ngood <- apply(coredata(x),2,iv) # KMP 2016-03-16
        maskna <- ngood==0 # KMP 2016-03-16
        #maskna <- !is.finite(colMeans(coredata(x))) # REB 2016-02-18 mask out same missing
        if (sum(maskna)>0) {                                  # data as in the first field
          z <- coredata(Z)
          z[,maskna] <- NA
          coredata(Z) <- z
        }                     
                                        #attr(Z,'dimensions') <- c(d1[1],d1[2],d2[3])
                                        #mostattributes(Z) <- attributes(y)
                                        #class(Z) <- class(y)

                                        # Add the regridded data as an attribute, but change class
        if (is.null(attr(x,'n.apps'))) n.app <- 1 else
        n.app <- attr(x,'n.apps') + 1
        attr(x,paste('appendix.',n.app,sep='')) <- Z
        attr(x,'n.apps') <- n.app

        X <- x
        if (sum(is.element(clsx,"comb"))==0)
            class(X) <- c("comb",clsx) else
        class(X) <- clsx

    }
    if (sum(is.element(dimension,"space"))) {
                                        # This option synchronises the data and expand the spatial grid
                                        # to include both data sets.
                                        #print(length(names(x[1,]))); print(length(x[1,])); print(names(x)[1:10])
                                        #xnm <- paste("x",1:length(x[1,]),sep=".") 
                                        #ynm <- paste("y",1:length(y[1,]),sep=".") 
                                        #names(x) <- xnm
                                        #names(y) <- ynm
                                        # REB 20.08.2013: added colnames to x & y before merge to keep track of
                                        # y & x.    
        colnames(y) <- paste("y",1:dim(y)[2],sep=".")
        colnames(x) <- paste("x",1:dim(x)[2],sep=".")
        comb <- merge(x,y,all=all)
        nt <- length(index(comb))
        if (verbose) str(comb)
        if (orig.format) {
            vars <- names(comb)
            if (verbose) print("After merge - original format")
            
            ys <- vars[grep('y',vars)]
            Xs <- vars[grep('x',vars)]
            ix <- is.element(vars,Xs)
            iy <- is.element(vars,ys)
            if (verbose) {print(paste(sum(ix),"xs and",sum(iy),"ys"))
                          print(dim(comb))}
            XX <- zoo(coredata(comb[,ix]),order.by=index(comb))
            yy <- zoo(coredata(comb[,iy]),order.by=index(comb))
            names(yy) <- ys
            names(XX) <- Xs

            if (verbose) print("set class")
                                        #print(clsx); print(clsy)
            clsx -> class(XX)
            clsy -> class(yy)

            if (verbose) {print("add attributes")
                          print(names(attributes(y)))}
                                        #nattr1 <- softattr(y)
            attr(y,'dimnames') <- NULL
            attr(x,'dimnames') <- NULL
            attr(yy,'dimnames') <- NULL
            attr(XX,'dimnames') <- NULL
                                        #for (i in 1:length(nattr1))
                                        #  attr(yy,nattr1[i]) <- attr(y,nattr1[i])
            yy <- attrcp(y,yy,ignore='names')
            attr(yy,'dimensions') <- c(attr(y,'dimensions')[1:2],nt)
            if (verbose) print(names(attributes(x)))
                                        #nattr2 <- softattr(x)
                                        #for (i in 1:length(nattr2))
                                        #  attr(XX,nattr2[i]) <- attr(x,nattr2[i])
            XX <- attrcp(x,XX,ignore='names')
            attr(XX,'dimensions') <- c(attr(x,'dimensions')[1:2],nt)

                                        #print("into list")
            X <- list(y=yy,X=XX)
            if (verbose) {print(names(attributes(X$y))); print(names(attributes(X$X)))}
        } else X <- comb
        class(X) <- c("comb","list")

    }
                                        #print("HERE")
                                        #print(names(attributes(combined$y)))
                                        #print(names(attributes(combined$X)))
                                        #print(summary(combined))
                                        #print("HERE")
    attr(X,'history') <- history.stamp(x)
    invisible(X)
}

combine.events <- function(x,y,remove.close=TRUE,mindistance=5E5,FUN=NULL,verbose=FALSE) {
  if(verbose) print("combine.events")
  stopifnot(inherits(x,"events") & inherits(y,"events"))

  if(is.null(attr(x,"greenwich"))) attr(x,"greenwich") <- !(any(x$lon<0) & !any(x$lon>180))
  x <- g2dl.events(x,greenwich=attr(x,"greenwich"))
  y <- g2dl.events(y,greenwich=attr(x,"greenwich"))

  ## Combine events in x and y
  cn <- colnames(x)[colnames(x) %in% colnames(y)]
  z <- rbind(x[colnames(x) %in% cn],y[colnames(y) %in% cn])
  
  if(!any(x$date %in% y$date)) remove.close <- FALSE
  
  ## Remove events located close to other stronger events
  if(remove.close) {
    t <- z$date*1E2 + z$time
    lon <- z$lon
    lat <- z$lat
    ## Define measure of strength.
    ## Decides the which one is thrown out if two events are close together. 
    if(is.null(FUN)) {
      if(attr(x,"variable")=="cyclones") {
        strength <- 1/z$pcent
      } else if (attr(x,"variable")=="anti-cyclones") {
        strength <- z$pcent
      } else {
        strength <- rep(1,nrow(z))
        print(paste("Warning: No measure of strength was provided.",
         "If two closely located events are detected one will be randomly removed.",
         "You can turn off the removing of neighbouring events (remove.close=FALSE)",
         "or provide a function (FUN) that describes the rank/strength of the events."))
      }
    } else {
      strength <- FUN(z)
    }
    if(verbose) print("Remove duplicate cyclones")
    del.z <- rep(TRUE,nrow(z))
    for (d in unique(t)) {
      i <- which(t==d)
      if (length(i)>1) {
        distance <- apply(cbind(lon[i],lat[i]),1,
         function(x) suppressWarnings(distAB(x[1],x[2],lon[i],lat[i])))
        diag(distance) <- NA; distance[lower.tri(distance)] <- NA
        del.i <- which(distance<mindistance,arr.ind=TRUE)
        if(any(del.i)) {
          col.del <- rep(1,length(del.i)/2)
          if (is.null(dim(del.i))) {
            s1 <- strength[i][,del.i[1]]
            s2 <- strength[i][,del.i[2]]
            col.del[s1<=s2] <- 2
            del.i <- unique(del.i[col.del])
          } else {
            s1 <- strength[i][del.i[,1]]
            s2 <- strength[i][del.i[,2]]
            col.del[s1<=s2] <- 2
            del.i <- unique(del.i[cbind(seq(1,dim(del.i)[1]),col.del)])
          }
          del.z[i[del.i]] <- FALSE
        }
        i <- i[!1:length(i) %in% del.i]
      }
    }
    z <- z[del.z,]
  }
  z <- attrcp(x,z)
  class(z) <- class(x)
  return(z)
}

g2dl <- function(x,greenwich=TRUE,...)
    UseMethod("g2dl")

g2dl.default <- function(x,greenwich=TRUE,lon=NULL,lat=NULL,d=NULL,verbose=FALSE) {
    if(verbose) print("g2dl.default")
    if (is.null(lon)) lon <- attr(x,'longitude')
    if (is.null(lat)) lat <- attr(x,'latitude')
    if (is.null(d)) d <- attr(x,'dimensions')
    if (greenwich) {
        wh <- lon < 0
        lon[wh] <- lon[wh] + 360
    } else {
        wh <- lon > 180
        lon[wh] <- lon[wh] - 360
    }
    y <- x
    if(length(lon)>1) {
      xsrt <- order(lon)
      xsrt <- xsrt[!xsrt %in% which(duplicated(lon))]
      dim(y) <- d
      y <- y[xsrt,,]
      lon <- lon[xsrt]
      dim(y) <- c(length(lon)*length(lat),d[3])
    }
    y <- attrcp(x,y)
    attr(y,'longitude') <- lon
    class(y) <- class(x)
    return(y)
}

g2dl.stationmeta <- function(x,greenwich=TRUE,verbose=FALSE) {
  if(verbose) print("g2dl.stationmeta")
  lon <- x$lon                          
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  y <- x
  y$lon <- lon
  attr(y,'greenwich') <- as.logical(greenwich)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

g2dl.field <- function(x,greenwich=TRUE,verbose=FALSE) {
    if(verbose) print("g2dl.field")
    attr(x,'longitude') -> lon
    attr(x,'latitude') -> lat
    d <- attr(x,'dimensions')
    
    if (greenwich) {
        wh <- lon < 0
        lon[wh] <- lon[wh] + 360
    } else {
        wh <- lon > 180
        lon[wh] <- lon[wh] - 360
    }
    
    xsrt <- order(lon)
    xsrt <- xsrt[!duplicated(lon)]
    X <- t(coredata(x))
    dim(X) <- d
    X <- X[xsrt,,]
    dim(X) <- c(length(lon)*d[2],d[3])
    y <- zoo(t(X),index(x))
    lon <- sort(lon)
    
    y <- attrcp(x,y,ignore='names')
                                        #nattr <- softattr(x,ignore=c('greenwich','longitude'))
                                        #for (i in 1:length(nattr))
                                        #  attr(y,nattr[i]) <- attr(x,nattr[i])
    attr(y,'dimensions') <- attr(x,'dimensions')
    attr(y,'longitude') <- lon
    attr(y,'greenwich') <- as.logical(greenwich)
    attr(y,'history') <- history.stamp(x)
    class(y) <- class(x)
    invisible(y)
}

g2dl.eof <- function(x,greenwich=TRUE,verbose=FALSE) {
    if(verbose) print("g2dl.eof")
    attr(x,'longitude') -> lon
    attr(x,'latitude') -> lat
    d <- attr(x,'dimensions')
    X <- attr(x,'pattern')
    if (greenwich) {
        wh <- lon < 0
        lon[wh] <- lon[wh] + 360
    } else {
        wh <- lon > 180
        lon[wh] <- lon[wh] - 360
    }
    xsrt <- order(lon)
    X <- X[xsrt,,]
    lon <- sort(lon)
    X -> attr(x,'pattern')
    attr(x,'greenwich') <- greenwich
    return(x)
}

g2dl.corfield <- function(x,greenwich=TRUE,verbose=FALSE) {
    if(verbose) print("g2dl.corfield")
    attr(x,'longitude') -> lon
    attr(x,'latitude') -> lat
    d <- attr(x,'dimensions')
                                        #print(d)
    if (greenwich) {
        wh <- lon < 0
        lon[wh] <- lon[wh] + 360
    } else {
        wh <- lon > 180
        lon[wh] <- lon[wh] - 360
    }

    y <- x
    xsrt <- order(lon)
    dim(y) <- d
    y <- y[xsrt,]
    y <- c(y)
    lon <- sort(lon)

    y <- attrcp(x,y,ignore='names')
                                        #nattr <- softattr(x,ignore=c('greenwich','longitude'))
                                        #for (i in 1:length(nattr))
                                        #  attr(y,nattr[i]) <- attr(x,nattr[i])
    attr(y,'dimensions') <- attr(x,'dimensions')
    attr(y,'longitude') <- lon
    attr(y,'greenwich') <- as.logical(greenwich)
    class(y) <- class(x)
    invisible(y)
}

g2dl.events <- function(x,greenwich=TRUE,verbose=FALSE) {
  if(verbose) print("g2dl.events")
  lon <- x$lon                          
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  y <- x
  y$lon <- lon
  attr(y,'greenwich') <- as.logical(greenwich)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

g2dl.trajectory <- function(x,greenwich=TRUE,verbose=FALSE) {
  if(verbose) print("g2dl.trajectory")
  lon <- x[,colnames(x)=="lon"]                         
  if (greenwich) {
    wh <- lon < 0
    lon[wh] <- lon[wh] + 360
  } else {
    wh <- lon > 180
    lon[wh] <- lon[wh] - 360
  }
  y <- x
  y[,colnames(y)=="lon"] <- lon
  attr(y,'greenwich') <- as.logical(greenwich)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

sp2np <- function(x,SP2NP=TRUE) {
                                        # Take care of different latitude orientations: N-> S & S -> N
    ysrt <- order(attr(x,'latitude'))
    d <- attr(x,'dimensions')
    if (SP2NP) {
        if (ysrt[1]==1) return(x)
        y <- t(coredata(x))
        dim(y) <- attr(x,'dimensions')
        y <- y[,ysrt,]
        dim(y) <- c(d[1]*d[2],d[3])
        y <- zoo(t(y),order.by=index(x))
        class(y) <- class(x)
        y <- attrcp(x,y,ignore='names')
                                        #nattr <- softattr(x,ignore=c('longitude','latitude','dimensions'))
                                        #for (i in 1:length(nattr))
                                        #  attr(y,nattr[i]) <- attr(x,nattr[i])
        attr(y,'longitude') <- attr(x,'longitude')
        attr(y,'latitude') <- attr(x,'latitude')[ysrt]
        attr(y,'dimensions') <- d
    } else {
        ysrt <- order(attr(x,'latitude'), decreasing = TRUE)
        if (ysrt[1]==1) return(x)
        y <- t(coredata(x))
        dim(y) <- d
        y <- y[,ysrt,]
        dim(y) <- c(d[1]*d[2],d[3])
        dim(y) <- c(d[1]*d[2],d[3])
        class(y) <- class(x)
        y <- attrcp(x,y,ignore='names')
                                        #nattr <- softattr(x,ignore=c('longitude','latitude','dimensions'))
                                        #for (i in 1:length(nattr))
                                        #  attr(y,nattr[i]) <- attr(x,nattr[i])
        attr(y,'longitude') <- attr(x,'longitude')
        attr(y,'latitude') <- attr(x,'latitude')[ysrt]
        attr(y,'dimensions') <- attr(x,'dimensions')
    }
    attr(y,'history') <- history.stamp(x)
    invisible(y)
}


combine.trajectory <- function(x,y) {
  z <- rbind(x,y,all=TRUE)
  z <- z[!duplicated(z),]
  z <- attrcp(x,z)
  class(z) <- class(x)
  return(z)
}
