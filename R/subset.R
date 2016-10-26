


subset.field <- function(x,it=NULL,is=NULL,verbose=FALSE) {
  if (is.null(it) & is.null(is)) return(x)
  if (verbose) print("subset.field")
  
  y <- default.subset(x,is=is,it=it,verbose=verbose)
  attr(y,'history') <- history.stamp(x)
  return(y)
}

subset.comb <- function(x,it=NULL,is=NULL,verbose=FALSE) {
    if (verbose) print("subset.comb")
    y <- subset.field(x,it=it,is=is)
    y <- attrcp(x,y)
    n.app <- attr(x,'n.apps')
                                        #print(n.app)
    for (i in 1:n.app) {
        eval(parse(text=paste("z <- attr(x,'appendix.",i,"')",sep="")))
        attr(z,'longitude') <- attr(x,'longitude')
        attr(z,'latitude') <- attr(x,'latitude')
        yz <- subset(z,it=it,is=is)
        yz <- zoo(coredata(yz),order.by=index(yz))
        yz <- attrcp(z,yz)
        attr(yz,'dimensions') <- c(length(attr(z,'longitude')),
                                   length(attr(z,'latitude')),
                                   length(index(yz)))
        eval(parse(text=paste("yz -> attr(y,'appendix.",i,"')",sep="")))
    }
    n.app -> attr(y,'n.apps')
    attr(y,'history') <- history.stamp(x)
    invisible(y)
}

subset.eof <- function(x,ip=NULL,it=NULL,is=NULL,verbose=FALSE) {
    if (verbose) print("subset.eof")
    ## browser()
    if (is.null(is) & is.null(it) & is.null(ip)) return(x)                                    
    if (is.null(it) & is.null(is[1]) & is.null(is[2]) & is.null(ip)) return(x) 
    d <- dim(x); greenwich <- TRUE
    clim <- attr(x,'mean')
    
    # Pattern extracts certain modes/patterns
    if (!is.null(ip)) {
      if (verbose) print(paste('Chose pattern',ip))
      y <- x[,ip]
      y <- attrcp(x,y)
      class(y) <- class(x)
      attr(y,'eigenvalues') <- attr(y,'eigenvalues')[ip]
      attr(y,'pattern') <- attr(y,'pattern')[,,ip]
      if (!is.null(attr(x,'n.apps'))) {
        attr(y,'n.apps') <- attr(x,'n.apps')
        attr(y,'appendix.1') <- attr(x,'appendix.1')
      }
      x <- y
    }
    
    if (is.null(is)) is <- 1:d[length(d)] else

    if (is.list(is)) {
      if (verbose) {print('Select space'); print(is)}
        if (length(is)==1) is[[2]] <- NULL
        if ( (is.null(is[[1]])) | (sum(is.finite(is[[1]])) < 2) ) is[[1]] <- c(-180,360)
        if ( (is.null(is[[2]])) | (sum(is.finite(is[[2]])) < 2) ) is[[2]] <- c(-90,90)
        
                                        # Select a subregion from the EOFs:
        lons <- attr(x,'longitude'); lon.rng <- range(lons)
        lats <- attr(x,'latitude'); lat.rng <- range(lats)
        X <- attr(x, "pattern")
        if ( (length(is[[1]])==2) & (length(is[[2]])==2) ) {
            lon.rng <- range(is[[1]]); lat.rng <- range(is[[2]])
            if ( (min(lon.rng) < 0) & (attr(x,'greenwich')) ) {
                if (verbose) print("convert to non-greenwich")
                lons[lons > 180] <- lons[lons > 180] - 360
                srt <- order(lons)
                X <- X[srt,,]
                clim <- clim[srt,]
                lons <- lons[srt]
                greenwich <- FALSE
            } 
            keepx <- (lons >= lon.rng[1]) & (lons <= lon.rng[2])
            keepy <- (lats >= lat.rng[1]) & (lats <= lat.rng[2])
            if ( (sum(keepx)==0) | (sum(keepy)==0) ) {
                print(is); print(lons); print(lats)
                stop('Check the coordinates')
            }
            X <- X[keepx,keepy,]
            if (verbose) print(paste('New dimensions',sum(keepx),sum(keepy)))
            attr(x, "pattern") <- X
            lons[keepx] -> attr(x,'longitude')
            lats[keepy] -> attr(x,'latitude')
            clim <- clim[keepx,keepy]
            attr(x,'dimensions') <- c(sum(keepx),sum(keepy),d[2])
            is <- 1:d[2]
        }
    } 
                                        #print(it)
    if (is.null(it)) {
      if (verbose) {print('Select time'); print(it)}
        it <- 1:d[1] 
        dates <- index(x)
        keep <- 1 : d[1]
    } else {
        if (is.numeric(it)) {
                                        #print("numeric it")
            if (max(nchar(as.character(it)))<=2)
                keep <- is.element(as.numeric(format(index(x),"%m")),it)
            else if (sum(nchar(as.character(it))==4) == length(it)) {
                if (length(it)==2)
                    if (diff(it)>1)
                        it <- seq(it[1],it[2],1)
                if (is.character(index(x)) | inherits(index(x),'Date'))
                    keep <- is.element(as.numeric(format(index(x),"%Y")),it)
                else
                    keep <- is.element(index(x),it)
            }
        }
        else if (is.character(it))
            keep <- is.element(index(x),as.Date(it))
                                        #print(c(sum(keep),length(keep)))
        dates <- index(x)[keep]
    }

    ## grep for appendices
    nm <- names(attributes(x))
    id <- grep("appendix",nm)
    if (length(id)>0) {
        nm <- nm[id]
        for (i in 1:length(nm)) {
            eval(parse(text=paste("a <- attr(x,'",nm,"')",sep="")))
            cls <- class(a)
            ais <- zoo(coredata(a)[keep,is],order.by=dates)
            ais <- attrcp(a,ais)
            eval(parse(text=paste("attr(x,'",nm,"') <- ais",sep="")))
            rm(a,ais,cls)
        }
    }
    
    class(x) -> cls
    ##keep <- is.element(index(x),it)
    y <- x[keep,]
    
    class(x) <- cls; class(y) <- cls
    y <- attrcp(x,y,ignore=c('greenwich','mean'))
    attr(y,'greenwich') <- greenwich
    clim -> attr(y,'mean')
    ## browser()
    attr(y,'pattern') <- attr(x,"pattern")
    attr(y,'eigenvalues') <-attr(x,"eigenvalues")
                                        #attr(y,'date-stamp') <- date()
                                        #attr(y,'call') <- match.call()
    attr(y, "dimensions") <- dim(attr(x,"pattern"))
    attr(y,'history') <- history.stamp(x)
    return(y)
}

subset.cca <- function(x,it=NULL,is=NULL) {
    if (!is.null(is))  {
        x <- subset.pattern(x)
    }
    x
}

subset.mvr <- function(x,it=NULL,is=NULL) {
    x
}

subset.pattern <- function(x,is,verbose=FALSE) {
  ## Takes a subset of the pattern attribute, e.g. a smaller region.
  if (verbose) print('subset.pattern')
    if (is.list(is)) {
        y <- attr(x,'pattern')
        lons <- attr(x,'longitude')
        lats <- attr(x,'latitude')
        nms <- substr(tolower(names(is)),1,3)
        IS <- 1:length(nms)
        if (verbose) print(nms)
        if (sum(is.element(nms,'lon'))>0) {
          inm <- IS[is.element(nms,'lon')]
            ix <- (lons >= min(is[[inm]])) &
                  (lons <= max(is[[inm]]))
          } else ix <- is.finite(lons)
        if (sum(is.element(nms,'lat'))>0) {
          inm <- IS[is.element(nms,'lat')]
            iy <- (lats >= min(is[[inm]])) &
                  (lats <= max(is[[inm]]))
          } else iy <- is.finite(lats)

        if (!is.null(attr(x,'pattern'))) {
          if (verbose) print('replace the pattern argument')
          y[ix,iy] -> attr(x,'pattern')
          lons[ix] -> attr(x,'longitude')
          lats[iy] -> attr(x,'latitude')
        } else {
          if (verbose) print(paste('subset the matrix:',sum(ix),sum(iy)))
          if (verbose) print(dim(x))
          y <- x[ix,iy]
          y <- attrcp(x,y)
          attr(y,'variable') <- varid(x)
          attr(y,'unit') <- unit(x)
          lons[ix] -> attr(y,'longitude')
          lats[iy] -> attr(y,'latitude')
          x <- y
        }
    } 
    return(x)
}

subset.matrix <- function(x,is,verbose=FALSE)
  subset.pattern(x,is,verbose=verbose)
  

subset.pca <- function(x,ip=NULL,it=NULL,is=NULL,verbose=FALSE) {
  if (verbose) print('subset.pca')
  if (!is.null(ip)) {
    y <- x[,ip]
    y <- attrcp(x,y)
    class(y) <- class(x)
    attr(y,'eigenvalues') <- attr(y,'eigenvalues')[ip]
    attr(y,'pattern') <- attr(y,'pattern')[,ip]
    if (!is.null(attr(x,'n.apps'))) {
      attr(y,'n.apps') <- attr(x,'n.apps')
      attr(y,'appendix.1') <- attr(x,'appendix.1')
    }
    x <- y
  }
  if (!is.null(it)) {
    y <- subset.station(x,it=it,verbose=verbose)
    y <- attrcp(x,y)
    class(y) <- class(x)
    x <- y
  }

  #browser()
  
  return(x)
}

subset.corfield <- function(x,it=NULL,is=NULL,verbose=FALSE) {
    if (verbose) print('subset.corfield')
    stopifnot(inherits(x,"corfield"))
    y <- x
    dim(y) <- c(1,length(x))
    y <- zoo(y,order.by=1)
    attr(y,"dimensions") <- c(length(attr(x,"longitude")),
                              length(attr(x,"latitude")),1)
    y <- attrcp(x,y)
    y <- as.field(y,param=attr(x,"variable"),unit=attr(x,"unit"),
                    lon=attr(x,"longitude"),lat=attr(x,"latitude"),
                    longname=attr(x,"longname"),src=attr(x,"source"),
                    url=attr(x,"url"))
    y <- subset(y,is=is,verbose=verbose)
    dim(y) <- length(y)
    attr(y,"dimensions") <- c(length(attr(y,"longitude")),
                              length(attr(y,"latitude")))
    class(y) <- "corfield"
    return(y)
}

subset.ds <- function(x,ip=NULL,it=NULL,is=NULL,verbose=FALSE) {
    if (verbose) print('subset.ds')
    y <- x
    if (!is.null(it)) {
      if (verbose) print(paste('it=',it))
    }
    if (!is.null(ip)) {
      if (verbose) print(paste('pattern=',ip))
      if (inherits(x,'pca')) {
        y <- subset.pca(x,ip=ip,verbose=verbose)
        attr(y,'eof') <- subset.eof(attr(x,'eof'),ip=ip,verbose=verbose)
        attr(y,'evaluation') <- attr(x,'evaluation')[,c(2*(ip-1)+1,2*(ip-1)+2)]
        if (!is.null(attr(x,'n.apps'))) {
          natt <- attr(x,'n.apps')
          for (i in 1:natt)
            attr(y,paste('appendix.',i,sep='')) <-
              attr(y,paste('appendix.',i,sep=''))[,ip]
        }
        x <- y
      }
    }
    if (!is.null(is))  {
      if (verbose) print(paste('is=',is))
        x <- subset.pattern(x,is,verbose=verbose)
    }
    x
}

subset.trend <- function(x,it=NULL,is=NULL) {
    y <- subset.field(x,it=it,is=is)
    
    pattern <- attr(x, "pattern")
    lons <- attr(x,'longitude'); lon.rng <- range(lons)
    lats <- attr(x,'latitude'); lat.rng <- range(lats)
    if ( (length(is[[1]])==2) & (length(is[[2]])==2) ) {
        lon.rng <- range(is[[1]]); lat.rng <- range(is[[2]])
        if ( (min(lon.rng) < 0) & (attr(x,'greenwich')) ) {
                                        #print("convert to non-greenwich")
            lons[lons > 180] <- lons[lons > 180] - 360
            srt <- order(lons)
            pattern <- pattern[srt,]
            lons <- lons[srt]
            greenwich <- FALSE
        } 
        keepx <- (lons >= lon.rng[1]) & (lons <= lon.rng[2])
        keepy <- (lats >= lat.rng[1]) & (lats <= lat.rng[2])
        if ( (sum(keepx)==0) | (sum(keepy)==0) ) {
            print(is); print(lons); print(lats)
            stop('Check the coordinates')
        }
        attr(y, "pattern") <- pattern
    }
    return(y)
}

subset.dsensemble <- function(x,it=NULL,is=NULL,
                              ensemble.aggregate=TRUE,verbose=FALSE,...) {
    ## browser()

    if (verbose) print('subset.dsensemble')

    if (inherits(x,'list') & inherits(x,c('pca','eof')) &
       (is.null(is)) & ensemble.aggregate) {
      #x <- as.station(x)
      ## Subset the PCA/EOF
      x <- subset.dsensemble.multi(x,it=it,is=is,verbose=verbose,...)
      return(x)
    }

    if (!is.null(is)) x <- as.station(x)
    
    if (inherits(x,'list') & !inherits(x,'zoo')) {
      if (verbose) print('list of elements')
      ## If x is a list of objects search through its elements
      Locs <- unlist(lapply(x,function(x) loc(attr(x,'station'))))
      Locs <- gsub(' ','',Locs)
      Locs <- gsub('-','.',Locs)
      if (is.character(is)) {
        
        if (verbose) print('search on location names')
      ## search on location name
        Locs <- tolower(Locs)
        locs <- substr(Locs,1,min(nchar(is)))
        if (verbose) {print(is); print(locs)}
        is <- substr(is,1,min(nchar(is)))
        illoc <- (1:length(x))[is.element(locs,tolower(is))]
        if (length(illoc)==1) {
          x2 <- x[[illoc]]
          x2 <- subset(x2,it=it,verbose=verbose)
        } else if (length(illoc)>1) {
          x2 <- list()
          for (i in 1:length(illoc)) {
            xx2 <- x[[illoc[i]]]
            xx2 <- subset(xx2,it=it,verbose=verbose)
            eval(parse(text=paste('x2$',Locs[illoc[i]],' <- xx2',sep='')))
            rm('xx2'); gc(reset=TRUE)
          } 
        } else if (length(illoc)==0) return(NULL)
        if (verbose) {print(is); print(loc(x2))}
        return(x2)
      }
    }
    class(x) <- c(class(x)[1],class(attr(x,'station'))[2],"zoo")

    if (is.null(it) & is.null(is) & length(table(month(x)))==1) return(x)
    if (verbose) {print("subset.dsensemble"); print(it)}
    x0 <- x
    d <- dim(x)
    if (verbose) print(d)
    if (is.null(is)) is <- 1:d[2]
    if (!is.null(it)) {
        if (is.character(it)) it <- tolower(it)
        if (verbose) print(table(month(x)))
        if ( (length(rownames(table(month(x))))==1) & (it==0) ) {
          if (verbose) print('Only one season is available')
          return(x)
          }
        ## Different ways of selecting along the time dimension
        if ( inherits(it[1],"logical") & (length(it)==length(x)) ) {
            if (verbose) print('subindexing with boolean index: y <- x[it,is]')
            y <- x[it,is]
        } else if (it[1]==0) {
            if (verbose) print("Annual means")
            if (inherits(x,'season')) {
                djf <- subset(x,it='djf',is=is)
                mam <- subset(x,it='mam',is=is)
                jja <- subset(x,it='jja',is=is)
                son <- subset(x,it='son',is=is)
                yr1 <- year(djf)
                yr2 <- year(mam)
                yr3 <- year(jja)
                yr4 <- year(son)
                yr <- yr1[is.element(yr1,yr2)]
                yr <- yr[is.element(yr,yr3)]
                yr <- yr[is.element(yr,yr4)]
                if (verbose) print(yr)
                i1 <- is.element(yr1,yr)
                i2 <- is.element(yr2,yr)
                i3 <- is.element(yr3,yr)
                i4 <- is.element(yr4,yr)
                if (verbose) print(c(sum(i1),sum(i2),sum(i3),sum(i4)))
                y <- zoo(0.25*(coredata(djf[i1,]) +
                               coredata(mam[i2,]) +
                               coredata(jja[i3,]) +
                               coredata(son[i4,])),order.by=yr)
#                         order.by=as.Date(paste(yr,'01-01',sep='-')))
                y <- attrcp(x0,y)
                class(y) <- class(x0)
            } else if (inherits(x,'month')) {     
                ## browser()
                jan <- subset(x,it='jan',is=is)
                feb <- subset(x,it='feb',is=is)
                mar <- subset(x,it='mar',is=is)
                apr <- subset(x,it='apr',is=is)
                may <- subset(x,it='may',is=is)
                jun <- subset(x,it='jun',is=is)
                jul <- subset(x,it='jul',is=is)
                aug <- subset(x,it='aug',is=is)         
                sep <- subset(x,it='sep',is=is)
                oct <- subset(x,it='oct',is=is)
                nov <- subset(x,it='nov',is=is)
                dec <- subset(x,it='dec',is=is)

                yr1 <- year(jan)
                yr2 <- year(feb)
                yr3 <- year(mar)
                yr4 <- year(apr)
                yr5 <- year(may)
                yr6 <- year(jun)
                yr7 <- year(jul)
                yr8 <- year(aug)
                yr9 <- year(sep)
                yr10 <- year(oct)
                yr11 <- year(nov)
                yr12 <- year(dec)

                yr <- yr1[is.element(yr1,yr2)]
                yr <- yr[is.element(yr,yr3)]
                yr <- yr[is.element(yr,yr4)]
                yr <- yr1[is.element(yr,yr5)]
                yr <- yr[is.element(yr,yr6)]
                yr <- yr[is.element(yr,yr7)]
                yr <- yr1[is.element(yr1,yr8)]
                yr <- yr[is.element(yr,yr9)]
                yr <- yr[is.element(yr,yr10)]
                yr <- yr[is.element(yr,yr11)]
                yr <- yr[is.element(yr,yr12)]
                
                if (verbose) print(yr)
                i1 <- is.element(yr1,yr)
                i2 <- is.element(yr2,yr)
                i3 <- is.element(yr3,yr)
                i4 <- is.element(yr4,yr)
                i5 <- is.element(yr5,yr)
                i6 <- is.element(yr6,yr)
                i7 <- is.element(yr7,yr)
                i8 <- is.element(yr8,yr)
                i9 <- is.element(yr9,yr)
                i10 <- is.element(yr10,yr)
                i11 <- is.element(yr11,yr)
                i12 <- is.element(yr12,yr) 
                
                if (verbose) print(c(sum(i1),sum(i2),sum(i3),sum(i4),
                                     sum(i5),sum(i6),sum(i7),sum(i8),
                                     sum(i9),sum(i10),sum(i11),sum(i12)))
                y <- zoo(1/12*(coredata(jan[i1,])+
                               coredata(feb[i2,])+
                               coredata(mar[i3,])+
                               coredata(apr[i4,])+
                               coredata(may[i5,])+
                               coredata(jun[i6,])+
                               coredata(jul[i7,])+
                               coredata(aug[i8,])+
                               coredata(sep[i9,])+
                               coredata(oct[i10,])+
                               coredata(nov[i11,])+
                               coredata(dec[i12,])),
                         order.by=as.Date(paste(yr,'01-01',sep='-')))
                y <- attrcp(x0,y)
                class(y) <- class(x0)  
            }
        } else if (is.character(it)) {
            ## browser()
            if (verbose) print("it is character - select a season")
            months <- month(x)
            if (sum(is.element(tolower(substr(it,1,3)),tolower(month.abb)))>0) {
                ii <- is.element(months,(1:12)[is.element(tolower(month.abb),
                                                          tolower(substr(it,1,3)))])
                if (verbose) print(ii)
                y <- x[ii,is]
            } else if (sum(is.element(tolower(it),names(season.abb())))>0) {
                if (verbose) print("season")
                sea <- eval(parse(text=paste('season.abb()$',it,sep='')))
                if (verbose) print(sea)
                ii <- is.element(months,sea)
                if (verbose) print(ii)
                y <- x[ii,is]
                ##      }
                ##    if ( (min(it) > 0) & (max(it) < 5) ) {
                ##      # dsensemble is not for monthly data
                ##      #print("Seasonal selection")
                ##      ii <- is.element(month(x),c(1,4,7,10)[it])
                ##      y <- x[ii,is]
            } else if (sum(is.element(tolower(it),tolower(month.abb)))>0) {
                ## browser()
                if (verbose) print("month")
                mon <- which(is.element(tolower(it),tolower(month.abb))>0)
                if (verbose) print(mon)
                ii <- is.element(months,mon)
                if (verbose) print(ii)
                y <- x[ii,is]
            }
        } else {
            if (sum(is.element(it,1600:2200)) > 0) {
                if (verbose) print("it contains year(s)")
                ii <- is.element(year(x),it)
                if (verbose) print(paste("Number of matches=",sum(ii)))
                y <- x[ii,is]
            } else if (is.character(it)) {
                if (verbose) print("Dates")
                x <- matchdate(x,it)
            } else if (inherits(is,c('field','station'))) {
                ## Match the times of another esd-data object
                if (verbose) print("Match date with another object")
                x <- matchdate(x,it)
            }
        }
        if (verbose) print("housekeeping")
        d[3] <- length(index(y))
        class(y) <- class(x0)
        d -> attr(y,'dimensions')
        ##str(y)
        y <- attrcp(x,y,ignore='station')
        ##print("station"); str(attr(x,'station')); print(class(attr(x,'station')))
        ## browser()
        ##plot(subset(attr(x,'station'),it=it))
        ##browser()
        ##        browser()
### BUG       if ( (it!=0) & (!inherits(attr(x,'station'),'annual')) & (it <= max(year(attr(x,'station')))) ) {
### Why '(it <= max(year(attr(x,'station'))))' - causes problems and it's not clear what the intention was     
        if ( (it!=0) & (!inherits(attr(x,'station'),'annual')) ) {
            if (verbose) print('Also extract the same data for the station')
            attr(y,'station') <- subset(attr(x,'station'),
                                        it=it,verbose=verbose) }
        else
            attr(y,'station') <- annual(attr(x,'station'))
                                        #print("HERE")
                                        #nattr <- softattr(x)
                                        #for (i in 1:length(nattr))
                                        #  attr(y,nattr[i]) <- attr(x,nattr[i])
    } else {
        ## Bug-fix AM: 2014-09-22
        y <- x[,is]
        y <- attrcp(x,y)
        attr(y, "model_id") <- attr(x, "model_id")[is]
        attr(y, "scorestats") <- attr(x, "scorestats")[is]
                                        #browser()
        attr(y,"history") <- history.stamp(x)
        if (length(is)==1) class(y) <- c("ds","zoo") else class(y) <- class(x)
    }
    if (verbose) print("exit subset.dsensemble")
    invisible(y)  
}

subset.spell <- function(x,is=NULL,it=NULL) {
    y <- subset.station(x,is=is,it=it)
    good <- is.finite(y)
    y <- zoo(y[good],order.by=index(y)[good])
    attr(y,'location') <- loc(x)[is]
    attr(y,'variable') <- varid(x)[is]
    attr(y,'unit') <- unit(x)[is]
    attr(y,'station_id') <- stid(x)[is]
    attr(y,'longitude') <- lon(x)[is]
    attr(y,'latitude') <- lat(x)[is]
    attr(y,'altitude') <- alt(x)[is]
    attr(y,'longname') <- attr(x,'longname')[is]
    attr(y,'aspect') <- attr(x,'aspect')[is]
    attr(y,'source') <- attr(x,'source')[is]
    attr(y,'URL') <- attr(x,'URL')[is]
    attr(y,'quality') <- attr(x,'quality')[is]
    attr(y,'country') <- attr(x,'country')[is]
    attr(y,'threshold') <- attr(x,'threshold')[is]
    attr(y,'threshold.unit') <- attr(x,'threshold.unit')[is]
    attr(y,'history') <- history.stamp(x)
    class(y) <- class(x)[-1]
    invisible(y)
}

subset.zoo <- function(x,it=NULL,is=NULL) subset.station(x,it=it,is=is)



## Author Rasmus E. Benestad - was initially part of subset.R file
## Modified by A. Mezghani
## Last update 06.01.2014 ; 24-02-2014

## subset.station <- function(x,it = NULL,is=NULL,loc=NULL , param = NULL,
##                            stid = NULL ,lon = NULL, lat = NULL, 
##                            alt = NULL, cntr = NULL, src = NULL , nmin = NULL,
##                            verbose=FALSE) {
    
##     ##
##     if (inherits(it,c('field','station','zoo'))) {
##         ## Match the times of another esd-data object
##         if (verbose) print('it: field/station')
##         x2 <- matchdate(x,it)
##         return(x2)
##     }

##     if (inherits(is,c('field','station','zoo'))) {
##         ## Match the times of another esd-data object
##         if (verbose) print('is: field/station')
##         browser()
##         x2 <- subset(x,loc=loc(is))
##         return(x2)
##     }
    
##     ##print("subset.station")
##     if (is.null(dim(x))) {
##         x2 <- default.subset(x,it=it,is=1,verbose=verbose)
##     } else {
##         ##print("here")
##         x2 <- default.subset(x,it=it,is=is,verbose=verbose)
##         ## 
##         ## extra selection based on meta data
##         ## ss <- select.station(x=x2,loc = loc , param = param,  stid = stid ,lon = lon, lat = lat, alt = alt, cntr = cntr, src = src , nmin = nmin)
##         ## 
##         ## if (!is.null(ss)) {
##         ##    id <- is.element(attr(x2,'station_id'),ss$station_id)
##         ## Keep selected stations only
##         ##    x2 <- station.subset(x2,it=it,is=which(id),verbose=verbose)
##         ##}
##         ##if (!is.null(is)) x2 <- station.subset(x2,it=it,is=is,verbose=verbose)
##     }
##     return(x2)
## }


default.subregion <- function(x,is=NULL,verbose=FALSE) {
  if (verbose) {print("Sub-region"); print(is)}
  
  if ( (is.list(is)) | (is.data.frame(is)) ) {
    if ( (is.null(is[[1]])) | (sum(is.finite(is[[1]])) < 2) ) is[[1]] <- c(-180,360)
    if ( (is.null(is[[2]])) | (sum(is.finite(is[[2]])) < 2) ) is[[2]] <- c(-90,90)
    nam <- names(is)
    if (is.null(nam)) nam <- c("lon","lat")

    # If ranges are provided
    if ( (length(is[[1]])==2) & (length(is[[2]])==2) &
        (nam[1]=="lon") & (nam[2]=="lat") ) {
      if (verbose) print('Select according to longitude-latitude ranges')
                # Select according to longitude-latitude ranges
      lon.rng <- range(is[[1]]); lat.rng <- range(is[[2]])
                #if ( (is.null(lon.rng)) | (sum(is.finite(lon.rng)) < 2) )
                #    lon.rng <- c(-180,180)
                #if ( (is.null(lat.rng))  | (sum(is.finite(lat.rng)) < 2) )
                #    lat.rng <- c(-90,90)
                #print("subset.field: lon.rng/lat.rng"); print(lon.rng); print(lat.rng)
                #print(class(x))
                #print("g2dl"); print(attr(x,'dimensions')); print(dim(x)); print(class(x))
      if ( (min(lon.rng) < 0) & (max(lon.rng) <= 180) )
        x <- g2dl.field(x,greenwich=FALSE) else
      if ( (min(lon.rng) >= 0) & (max(lon.rng) > 180) )
        x <- g2dl.field(x,greenwich=TRUE)
                                        #print("sp2np")
      x <- sp2np(x)
                
                #print("subset.field: lon.rng/lat.rng"); print(lon.rng); print(lat.rng)
      d <- attr(x,'dimensions')
                                        #print(d)
      xy <- rep(attr(x,'longitude'),d[2])
      yx <- sort(rep(attr(x,'latitude'),d[1]))
      inside <- (xy >= lon.rng[1]) &
                (xy <= lon.rng[2]) &
                (yx >= lat.rng[1]) &
                (yx <= lat.rng[2])
                                        #print(c(length(inside),sum(inside)))
      y <- x[,inside]
      ix <- (attr(x,'longitude') >= lon.rng[1]) &
            (attr(x,'longitude') <= lon.rng[2])
      iy <- (attr(x,'latitude') >= lat.rng[1]) &
            (attr(x,'latitude') <= lat.rng[2])
      d <- c(sum(ix),sum(iy),d[3])
      d -> attr(y,'dimensions')
      attr(y,'longitude') <- attr(x,'longitude')[ix]
      attr(y,'latitude') <- attr(x,'latitude')[iy]
      attr(y,'ix') <- ix
      attr(y,'iy') <- iy
      attr(y,'ixy') <- inside
                                        #print(d)
    } else {
      if (verbose) print('Select according to given coordinates')
                                        # Select the specific grid point values from the longitude
                                        # and latitude coordinates
      lon.pick <- is[[1]]; lat.pick <- is[[2]]
      d <- attr(x,'dimensions')
      if (verbose) print(paste("Select ",length(lon.pick),"x",length(lat.pick),
                               "grid points from the",d[1],"x",d[2],"grid"))
                                        #print(d)
      xy <- rep(attr(x,'longitude'),d[2])
      yx <- sort(rep(attr(x,'latitude'),d[1]))
      ixy <- is.element(xy, lon.pick) &
      is.element(yx, lat.pick)
      ix <- is.element(attr(x,'longitude'),lon.pick)
      iy <- is.element(attr(x,'latitude'),lat.pick)
      y <- x[,ixy]
      d <-  c(sum(ix),sum(iy),d[3])
      d -> attr(y,'dimensions')
                         #print(paste("d=",d[1],d[2],d[3]," sum(ixy)=",sum(ixy),"=",d[1]*d[2]))
      attr(y,'longitude') <- attr(x,'longitude')[ix]
      attr(y,'latitude') <- attr(x,'latitude')[iy]
      attr(y,'ix') <- ix
      attr(y,'iy') <- iy
      attr(y,'ixy') <- ixy
    }
  }
  return(y)
} 
    
  

default.subset <- function(x,it=NULL,is=NULL,verbose=FALSE) {

    ## REB: Use select.station to condition the selection index is...
    ## loc - selection by names
    ## lon/lat selection be geography or closest if one coordinate lon/lat
    ##         if two-element vectors, define a region
    ## alt - positive values: any above; negative any below height
    ## cntr - selection by country
  
    ## REB 2015-02-02: renamed to default.subset because this will be used to subset
    ## both field and station objects.
    
    nval <- function(x) sum(is.finite(x))
    
    ## Sometimes 'it' = 'integer(0)' - reset to NULL!
    if (length(it)==0) it <- NULL
    if (length(is)==0) is <- NULL
    ## Return the original value if 'it' and 'is' are not specified
    if (is.null(it) & is.null(is)) return(x)
    
    if (verbose) {print("default.subset"); print(it); print(is); print('---')}
    x0 <- x
    ## 
    d <- dim(x)
    if (is.null(d)) {
        if (verbose)
            print("default.subset: Warning - One dimensional vector has been found in the coredata")
        x <- zoo(as.matrix(coredata(x)),order.by=index(x))
        x <- attrcp(x0,x)
        class(x) <- class(x0)
    } 
    d <- dim(x)
    if (is.null(is)) is <- 1:d[2]
#    if (is.null(it)) it <- 1:d[1] This lines causes a bug if is is given but not it...
    
    ## 
    ##print("HERE")
    ## get time in t
    t <- index(x)
    
    ## KMP 2016-02-03: to solve problem with subset.events 
    if(inherits(t,"Date")) t <- as.Date(round(as.numeric(t)))
    if(!inherits(t,"POSIXt")) ii <- is.finite(t) else ii <- rep(TRUE,length(t))
    if (verbose) {print('default.subset: time index it'); print(it)}

    ##  if (datetype=="Date") {
    if (inherits(t,c("Date","yearmon"))) {
       if (verbose) print('x is a Date or yearmon object')
        ## REB: replaced by lines below:
        ##    year <- as.numeric( format(t, '%Y') ) 
        ##    month <- as.numeric( format(t, '%m') )
        yr <- year(x)
        mo <- month(x)
        dy <- day(x)
    } else if (inherits(t,c("numeric","integer"))) {
        if (verbose) print('X has a numeric index - select by years')
        yr <- t
        mo <- dy <- rep(1,length(t))
    } else if (inherits(t,"POSIXt")) {
        if (verbose) print('X has a POSIXt index')
        yr <- year(t)
        mo <- month(t)
        dy <- day(t)
        hr <- as.numeric(format(t,"%H"))
        mn <- as.numeric(format(t,"%M"))
        if (!inherits(it,"POSIXt")) t <-  as.Date(format(t,"%Y-%m-%d"))
    } else print("Index of x should be a Date, yearmon, or numeric object")
    
    ## Generate sequence of days, months or years if range of it value is given
    if (is.character(it)) {
        if (verbose) print('it is character')
        if ((levels(factor(nchar(it)))==10)) ##  convert as Date
            it <- as.Date(it)
        if ( length(it) == 2 ) {
            if (verbose) print('Between two dates')
            # REB 2015-01-20: Need to convert it to dates -
            # otherwise the code crashes! 
           if (nchar(it[1])==4) it <- as.Date(c(paste(it[1],'-01-01',sep=''), 
                                                paste(it[2],'-12-31',sep='')))
            if (verbose) {print(it); print(class(x))}
            #browser()
            if (inherits(x,"month")) ## it is a month
                it <- seq(it[1],it[2],by='month') else
            if (inherits(x,"season")) ## it is a season
                it <- seq(it[1],it[2],by='season') else
            if (inherits(x,"annual")) ## it is a year
                it <- seq(it[1],it[2],by='year') else 
            if (inherits(x,"day")) ## assume it is a day
                it <- seq(it[1],it[2],by='day')
            ii <- is.element(t,it)
            #print(it); print(t)
            if (verbose) print(paste('sum(ii)=',sum(ii)))
        } else { ## added AM 10-06-2014
            if (verbose) print('it is a string')
            if (sum(is.element(tolower(substr(it,1,3)),tolower(month.abb)))>0) {
                if (verbose) print('Monthly selected')
                ii <- is.element(month(x),(1:12)[is.element(tolower(month.abb),tolower(substr(it,1,3)))])
                                        #y <- x[ii,is] #  REB Not here
            } else if (sum(is.element(tolower(it),names(season.abb())))>0) {
                if (verbose) print("Seasonally selected")
                if (verbose) print(table(month(x)))
                if (verbose) print(eval(parse(text=paste('season.abb()$',it,sep=''))))
                ii <- is.element(month(x),eval(parse(text=paste('season.abb()$',it,sep=''))))
                                        #y <- x[ii,is] # REB Not here
            }
            else if (inherits(it,"Date")) {
                if (verbose) print('it is a Date object')
                ii <- is.element(t,it)
            } else {
                str(it); print(class(it))
                ii <- rep(FALSE,length(t))
                warning("default.subset: did not recognise the selection citerion for 'it'")
            }
        }
        ##    } else if ((class(it)=="numeric") | (class(it)=="integer")) {
##        if (min(it) > 1500) ## it is a year
##            it <- seq(it[1],it[2],by=1)
##        ##print("HERE"); print(it)
##    }
##    ii <- (t >= it[1]) & (t <= it[length(it)])
##
    ## get the subset indices in ii
    } else if ((class(it)=="numeric") | (class(it)=="integer")) {
        if (verbose) print('it is numeric or integer')
# REB bug        nlev <- as.numeric(levels(factor(nchar(it))))
# nchar returns the string length, but these lines need to find the number of different levels/categories
        nlev <- as.numeric(levels(factor(as.character(it)))) # REB 2015-01-15
        if (verbose) {print(nlev); print(it)}
#         if ((length(nlev)==1)) { REB 2015-01-20: the lines below will never happen with this line:
        if (length(it)==2) {
          if ( (min(it) >= 1800) & (max(it) <= 2500) ) {
            if (verbose) print("it most probably contains a years")
            ii <- is.element(yr,year(it[1]):year(it[2]))
          } else if ( (min(it) >= 1) & (max(it) <= length(yr)) ) {
            if (verbose) print("it most probably contains a indices")
            ii <- is.element(1:length(yr),it[1]:it[2])
          } else  if (min(it) >= min(yr)) {
            if (verbose) print("it most probably contains years")
            ii <- is.element(yr,it[1]:max(yr))
          } else  if (max(it) <= max(yr)) {
            if (verbose) print("it most probably contains years")
            ii <- is.element(yr,min(yr):it[2])
          }
        } else if ((length(it)>2) | length(it==1)) {
          # if it is years:
                if (min(it) > length(yr)) {
                  if (verbose) print("match years")
                  ii <- is.element(yr,it)
                } else if (max(it) <= length(yr)) {
                  if (verbose) print("pick by indices")
                  ii <- is.element(1:length(t),it)
                } else {
                  ii <- rep(FALSE,length(t))
                  warning("default.subset: did not reckognise the selection citerion for 'it'")
                }
        }
    } else if (inherits(it,c("Date","yearmon"))) {       
        ##        ii <- is.element(t,it)
        if (verbose) print('it is a date object')
        ii <- (t >= min(it)) & (t <= max(it))
    } else if (inherits(it,"logical") & length(it)==length(yr)) {
        ii <- it
    } else if (inherits(it,"POSIXt")) {
        if (verbose) print('it is a POSIXt date & time object')
        if (!inherits(t,"POSIXt")) it <- as.Date(it)
        ii <- is.element(t,it)
    } else if (!is.null(it)) {
        ii <- rep(FALSE,length(t))
        warning("default.subset: did not reckognise the selection citerion for 'it'")
    } 
    
    ## it <- (1:length(t))[ii]
    ## 

    class(x) -> cls
    ##print(cls)
    ## update the class of x
    #class(x) <- "zoo" 

    n <- dim(x)[2]
    selx <- is.finite(lon(x)); sely <- is.finite(lat(x))
      
    selz <- rep(TRUE,n)
    selc <- selz; seli <- selz; selm <- selz; salt <- selz
    selp <- selz; selF <- selz ; sell <- selz

    # REB 11.04.2014: is can be a list to select region or according to other criterion
    if ( inherits(is,'list') & inherits(x,'station') ) {
        if (verbose) {print('spatial selection'); print(is)}
        nms <- names(is)
        il <- grep('loc',tolower(nms))
        ix <- grep('lon',tolower(nms))
        iy <- grep('lat',tolower(nms))
                                        #print(nms); print(c(ix,iy))
        iz <- grep('alt',tolower(nms))
        ic <- grep('cntr',tolower(nms))
        im <- grep('nmin',tolower(nms))
        ip <- grep('param',tolower(nms))
        id <- grep('stid',tolower(nms))
        iF <- grep('FUN',nms)
        if (length(il)>0) sloc <- is[[il]] else sloc <- NULL
        if (length(ix)>0) slon <- is[[ix]] else slon <- NULL
        if (length(iy)>0) slat <- is[[iy]] else slat <- NULL
                                        #print(slon); print(range(lon(x)))
        if (length(iz)>0) salt <- is[[iz]] else salt <- NULL
        if (length(ic)>0) scntr <- is[[ic]] else scntr <- NULL
        if (length(im)>0) snmin <- is[[im]] else snmin <- NULL
        if (length(ip)>0) sparam <- is[[ip]] else sparam <- NULL        
        if (length(id)>0) sstid <- is[[id]] else sstid <- NULL
        if (length(iF)>0) sFUN <- is[[iF]] else sFUN <- NULL
                                        #print(slat); print(range(lat(x)))
        if (length(sloc)>0) sell <- is.element(tolower(sloc(x)),sloc)

        # The lines for selx/sely differ for station and field objects:
        if (length(slon)==2) selx <- (lon(x) >= min(slon)) & (lon(x) <= max(slon))
        if (length(slat)==2) sely <- (lat(x) >= min(slat)) & (lat(x) <= max(slat))

        
        if (length(salt)==2) selz <- (alt(x) >= min(salt)) & (alt(x) <= max(salt))
        if (length(salt)==1) {
            if (salt < 0) selz <- alt(x) <= abs(salt) else
            selz <- alt(x) >= salt
        }
        if (length(scntr)>0) selc <- is.element(tolower(cntr(x)),scntr)
        if (length(snmin)>0) selm <- apply(coredata(x),2,nval) > snmin
        if (length(sparam)>0) selp <- is.element(tolower(param(x)),sparam)
        if (length(sstid)==2) seli <- (stid(x) >= min(sstid)) & (stid(x) <= max(sstid)) else
        if (length(sstid)>0) seli <- is.element(stid(x),sstid)
        if (length(sFUN)>0) selm <- apply(coredata(x),2,sFUN) # Not quite finished...
        ##
        is <- sell & selx & sely & selz & selc & seli & selm & selp & selF
        if (verbose) print(paste(sum(is),'spatial points'))
        ##
        ## Need to make sure both it and is are same type: here integers for index rather than logical
        ## otherwise the subindexing results in an empty object
    } else if ( inherits(is,'list') & inherits(x,'field') ) {
      ## KMP 2016-10-20 Can we subset across the dateline and greenwich now?
      y <- default.subregion(x,is=is,verbose=verbose)
      if(!any(attr(y,"longitude")<0) & any(attr(y,"longitude")>180)) {
        x <- g2dl.field(x,greenwich=TRUE) }
      is <- attr(y,'ixy'); selx <- attr(y,'ix'); sely <- attr(y,'iy')
    } else
    if ( is.null(is) ) is <- rep(TRUE,d[2]) else
    if ( is.numeric(is)) {
      iss <- rep(FALSE,d[2]); iss[is] <- TRUE
      is <- iss
    }
    
    if (verbose) print(paste('number of points:',sum(ii),sum(is)))
    y <- x[which(ii),which(is)] ## AM 2015-02-17 do not work with logical indexes
    #if (is.logical(is))
    #    is <- (1:length(is))[is]
    ##else 
    ##    is <- is.element(1:d[2],is)
    
    class(x) <- cls; class(y) <- cls
    y <- attrcp(x,y,ignore=c("names"))
    if (inherits(x,'station')) {
      if (verbose) print('station attributes')
      attr(y,'longitude') <- attr(x,'longitude')[is]
      attr(y,'latitude') <- attr(x,'latitude')[is]
      if (!is.null(attr(y,'altitude')))
        attr(y,'altitude') <- attr(x,'altitude')[is]
      if (!is.null(attr(y,'country')))
        attr(y,'country') <- attr(x,'country')[is]
      if (!is.null(attr(y,'source')))
        attr(y,'source') <- attr(x,'source')[is]
      if (!is.null(attr(y,'station_id')))
        attr(y,'station_id') <- attr(x,'station_id')[is]
      if (!is.null(attr(y,'location')))
        attr(y,'location') <- attr(x,'location')[is]
      if (!is.null(attr(y,'quality')))
        attr(y,'quality') <- attr(x,'quality')[is]
    ## attr(y,'history') <- attr(x,'history')[is]
      if (!is.null(attr(y,'variable')))
        attr(y,'variable') <- attr(x,'variable')[is]
    ## attr(y,'element') <- attr(x,'element')[is]
      if (!is.null(attr(y,'aspect')))
        attr(y,'aspect') <- attr(x,'aspect')[is]
      if (!is.null(attr(y,'unit')))
        if (length(attr(x,'unit'))==length(is)) attr(y,'unit') <- attr(x,'unit')[is] else
                                                attr(y,'unit') <- attr(x,'unit')
      if (!is.null(attr(y,'longname')))
        attr(y,'longname') <- attr(x,'longname')[is]
      if (!is.null(attr(y,'reference')))
        attr(y,'reference') <- attr(x,'reference')[is]
      if (!is.null(attr(y,'info')))
        attr(y,'info') <- attr(x,'info')[is]
      if (!is.null(attr(y,'method')))
        attr(y,'method') <- attr(x,'method')[is]
      if (!is.null(attr(y,'type')))
        attr(y,'type') <- attr(x,'type')[is]
      if (!is.null(attr(y,'URL')))
        attr(y,'URL') <- attr(x,'URL')[is]
      if (!is.null(attr(y,'na')))
        attr(y,'na') <- attr(x,'na')[is]
    
      if (!is.null(err(y)))
        attr(y,'standard.error') <- err(x)[ii,is]
    } else {
      attr(y,'longitude') <- attr(x,'longitude')[selx]
      attr(y,'latitude') <- attr(x,'latitude')[sely]
      c(sum(selx),sum(sely),sum(ii,na.rm=TRUE)) -> attr(y,'dimensions')
    }

    if(!any(attr(y,"longitude")<0) & any(attr(y,"longitude")>180)) {
      attr(y,"greenwich") <- TRUE
    } else {
      attr(y,"greenwich") <- FALSE
    }
    
    ##attr(y,'date-stamp') <- date()
    ##attr(y,'call') <- match.call()
    attr(y,'history') <- history.stamp(x)   
    if (verbose) print('exit default.subset')
    if (inherits(y,"annual")) index(y) <- as.numeric(year(index(y)))
    return(y)
}
    

subset.events <- function(x,it=NULL,is=NULL,verbose=FALSE,...) {
  if(verbose) print("subset.events")
  cls <- class(x)
  
  if (length(it)==0) it <- NULL
  if (length(is)==0) is <- NULL
  ii <- rep(TRUE,dim(x)[1])
  
  if(!is.null(it)) {
    dt <- x[,"date"]*1E2 + x[,"time"]
    t <- as.Date(strptime(x[,"date"],format="%Y%m%d"))
    if(verbose) print(paste('length of t',length(t)))
    
    is.datetime <- function(x) all(!is.months(x) &
                            (is.character(x) &
                             !grepl("-",x) &
                             all(levels(factor(nchar(x)))==10)) |
                            (is.numeric(x) &
                             all(levels(factor(nchar(x)))==10)) |
                            inherits(x,c("POSIXt")))
    is.dates <- function(x) all(!is.months(x) & !is.datetime(x) & 
                            (is.character(x) &
                            all(levels(factor(nchar(x)))==10) |
                            all(levels(factor(nchar(x)))==8)) |
                            (is.numeric(x) & all(levels(factor(nchar(x)))==8)) |
                            inherits(x,"Date"))
    is.years <- function(x) all(!is.months(x) & 
                            is.numeric(x) & levels(factor(nchar(x)))==4)
    is.months <- function(x) all(sum(is.element(tolower(substr(x,1,3)),
                                                tolower(month.abb)))>0)
    is.seasons <- function(x) all(sum(is.element(tolower(substr(x,1,3)),
                                                 names(season.abb())))>0)
                            
    if (is.months(it)) {
      if (verbose) print('Monthly selected')
      mo <- month(t)
      ii <- is.element(mo,(1:12)[is.element(tolower(month.abb),
                               tolower(substr(it,1,3)))])
    } else if (is.seasons(it)) {
      if (verbose) print("Seasonally selected")
      mo <- month(t)
      if (verbose) print(table(mo))
      if (verbose) print(eval(parse(text=paste('season.abb()$',it,sep=''))))
      ii <- is.element(mo,eval(parse(text=paste('season.abb()$',it,sep=''))))
    } else if (is.datetime(it)) {
      if (verbose) print("Date and time")
      if (inherits(it,c("POSIXt"))) it <- as.numeric(strftime(it,"%Y%m%d%H"))
      if (is.character(it)) it <- as.numeric(it)
      ii <- is.element(dt,it)
    } else if (is.dates(it)) {
      if (is.character(it) & all(grepl("-",it))) {
        it <- as.Date(it)
      } else if (!inherits(it,"Date")) {
        it <- as.Date(strptime(it,format="%Y%m%d"))
      }
      if ( length(it) == 2 ) {
        if (verbose) print('Between two dates')
        if (verbose) print(it)
        it <- strftime(seq(it[1],it[2],by='day'),format="%Y%m%d")
        t <- strftime(t,format="%Y%m%d")
        ii <- is.element(t,it)
      } else { 
        if (verbose) print('it is a string of dates')
        ii <- is.element(t,it)
      }
   } else if (is.years(it)) {
      yr <- year(t)
      it <- as.integer(it)
      if ( length(it) == 2 ) {
        if (verbose) print('Between two years')
        if (verbose) print(it)          
        ii <- is.element(yr,seq(it[1],it[2],1))
      } else { 
        if (verbose) print('it is a string of years')
        ii <- is.element(yr,it)
      }
    } else if (is.logical(it) & length(it)==length(t)) {
      if (verbose) print('it is a logical array')
      ii <- it
    } else if (is.integer(it) & max(it)<=length(t)) {
      if (verbose) print('it is an index array')
      ii <- rep(FALSE,length(t))
      ii[it] <- TRUE
    } else {
      ii <- rep(FALSE,length(t))
      warning("subset.station: did not recognise the selection citerion for 'it'")
    }
  }  

  jj <- rep(TRUE,dim(x)[1])
  if (inherits(is,'list')) {
    if (verbose) print('is is a list:')
    nm.x <- names(x)
    nm.is <- names(is)
    ok <- sapply(nm.is,function(n) any(grep(n,nm.x)))
    if (verbose) print(nm.is[ok])
    for (n in nm.is[ok]) {
      jj <- jj & x[n][[1]]>=min(is[n][[1]]) &
                 x[n][[1]]<=max(is[n][[1]])
    }
  } else if (is.numeric(is)) {
    jj <- 1:dim(x)[1] %in% it
  } else if (is.logical(is) & length(is)==dim(x)[1]) {
    jj <- is  
  } else if (!is.null(is)){
    jj <- rep(FALSE,dim(x)[1])
    warning("default.subset: did not reckognise the selection citerion for 'is'")
  }
  
  ij <- ii & jj
  y <- x[ij,]
  attr(y,"aspect") <- "subset"
  if (!is.null(is$lat)) attr(y,"lat") <- is$lat
  if (!is.null(is$lon)) attr(y,"lon") <- is$lon
  class(y) <- cls
  invisible(y)
}

## Routine for sorting the order of station series.
sort.station <- function(x,is=NULL,decreasing=TRUE) {
  if (is.null(is)) is <- order(stid(x),decreasing=decreasing)
  y <- zoo(x)[,is]
  y <- attrcp(x,y)
  attr(y,'station_id') <- stid(x)[is]
  attr(y,'location') <- loc(x)[is]
  attr(y,'variable') <- varid(x)[is]
  attr(y,'unit') <- unit(x)[is]
  attr(y,'longitude') <- lon(x)[is]
  attr(y,'latitude') <- lat(x)[is]
  attr(y,'altitude') <- alt(x)[is]
  attr(y,'cntr') <- cntr(x)[is]
  attr(y,'longname') <- attr(x,'longname')[is]
  attr(y,'history') <- history.stamp(x)
  return(y)
}
