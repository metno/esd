                                        #subset <- function(x,it=NULL,is=NULL,...)
                                        #                   UseMethod("subset")

                                        #subset.station <- function(x,it=NULL,is=NULL,
                                        #                           loc=NULL,lon=NULL,lat=NULL,alt=NULL,
                                        #                           cntr=NULL) {
                                        #  # REB: Use select.station to condition the selection index is...
                                        #  # loc - selection by names
                                        #  # lon/lat selection be geography or closest if one coordinate lon/lat
                                        #  #         if two-element vectors, define a region
                                        #  # alt - positive values: any above; negative any below height
                                        #  # cntr - selection by country
                                        #
                                        #  x0 <- x
                                        #  if (is.null(it) & is.null(is) & is.null(loc) & is.null(lon)
                                        #      & is.null(lat) & is.null(alt) & is.null(cntr)) return(x)
                                        #  d <- dim(x)
                                        #  if (is.null(d)) d <- c(length(x),1)
                                        #  if (is.null(it)) it <- 1:d[1]
                                        #  if (is.null(is)) is <- 1:d[2]
                                        #  
                                        #  t <- index(x)
                                        #  datetype <- class(t)
                                        #  if (datetype=="Date") {
                                        #    year <- as.numeric( format(t, '%Y') ) 
                                        #    month <- as.numeric( format(t, '%m') )
                                        #  } else
                                        #  if (datetype=="numeric") year <- t
                                        #  class(x) -> cls
                                        #  #print(cls)
                                        #  class(x) <- "zoo"
                                        #  if ( (class(it)=="logical") & (length(it)==length(x)) )
                                        #      y <- x[it,is] else
                                        #  if ( (min(it) > 0) & (max(it) < 13) & (inherits(x0,c("month"))) ) {
                                        #      #keepm <- as.numeric(format(index(X),"%m"))==it
                                        #      #print("Monthly aggregated field")
                                        #      #keepm <- is.element(as.POSIXlt(dates)$mon+1,it)
                                        #      ii <- is.element(month,it)
                                        #      y <- x[ii,is]
                                        #    } else 
                                        #    if ( (min(it) > 0) & (max(it) < 5) & (inherits(x0,c("season"))) ) {
                                        #      print("Seasonally aggregated field")
                                        #      #print(table(as.POSIXlt(dates)$mon+1))
                                        #      #keepm <- is.element(as.POSIXlt(dates)$mon+1,c(1,4,7,10)[it])
                                        #      #print(c(it,sum(keepm)))
                                        #      ii <- is.element(month,c(1,4,7,10)[it])
                                        #      y <- x[ii,is]
                                        #    } else
                                        #  if ( (min(it) > 0) & (max(it) < 13) & (sum(is.element(it,1:12)) > 0) ) {
                                        #
                                        #      } else
                                        #  if (sum(is.element(it,1600:2200)) > 0) {
                                        #        ii <- is.element(year,it)
                                        #        y <- x[ii,is]
                                        #      } else
                                        #  if (is.character(it)) {
                                        #        dates <- as.Date(it)
                                        #        ii <- is.element(index(x),dates)
                                        #        y <- x[ii,is]
                                        #      } else y <- x[it,is]
                                        #  
                                        #  class(x) <- cls; class(y) <- cls
                                        #  y <- attrcp(x,y)
                                        #  #nattr <- softattr(x)
                                        #  #for (i in 1:length(nattr))
                                        #  #  attr(y,nattr[i]) <- attr(x,nattr[i])
                                        #  mostattributes(y) <- attributes(x)
                                        #  attr(y,'longitude') <- attr(x,'longitude')[is]
                                        #  attr(y,'latitude') <- attr(x,'latitude')[is]
                                        #  attr(y,'altitude') <- attr(x,'altitude')[is]
                                        #  attr(y,'country') <- attr(x,'country')[is]
                                        #  attr(y,'source') <- attr(x,'source')[is]
                                        #  attr(y,'station_id') <- attr(x,'station_id')[is]
                                        #  attr(y,'location') <- attr(x,'location')[is]
                                        #  attr(y,'quality') <- attr(x,'quality')[is]
                                        #  attr(y,'URL') <- attr(x,'URL')[is]
                                        #  attr(y,'history') <- attr(x,'history')[is]
                                        #  attr(y,'variable') <- attr(x,'variable')[is]
                                        #  attr(y,'aspect') <- attr(x,'aspect')[is]
                                        #  attr(y,'unit') <- attr(x,'unit')[is]
                                        #  attr(y,'longname') <- attr(x,'longname')[is]
                                        #  attr(y,'reference') <- attr(x,'reference')[is]
                                        #  attr(y,'info') <- attr(x,'info')[is]
                                        #  attr(y,'date-stamp') <- date()
                                        #  attr(y,'call') <- match.call()
                                        #  return(y)
                                        #}


subset.field <- function(x,it=NULL,is=NULL,verbose=FALSE) {
                                        #print("subset.field")
        
    x0 <- x
    if (is.null(it) & is.null(is)) return(x)
    ## if (is.null(it) & is.null(is[[1]]) & is.null(is[[2]])) return(x) 
    t <- index(x)
    datetype <- class(t)
    years <- year(x)
    months=month(x)
                                        #  if (datetype=="Date") {
                                        #    years <- as.numeric( format(t, '%Y') ) 
                                        #    months <- as.numeric( format(t, '%m') )
                                        #  } else
                                        #  if (datetype=="numeric") years <- t

    class(x) -> cls
                                        #print(cls)
                                        #print(years); print(it)
    
    if (is.null(is)) is <- 1:dim(x)[2]
    
    if (!is.null(it)) {
                                        #print("select time"); print(it)
                                        #  if (sum(is.element(dimension,"time"))) {
        
        class(x) <- "zoo"
        d <- attr(x,'dimensions')
        if ( inherits(it[1],"logical") & (length(it)==length(x)) )
            y <- x[it,]
        else if (is.character(it)) {
            if (levels(factor(nchar(it)))==10)
                it <- as.Date(it)
            if (sum(is.element(tolower(substr(it,1,3)),tolower(month.abb)))>0) {
                ii <- is.element(months,(1:12)[is.element(tolower(month.abb),
                                                          tolower(substr(it,1,3)))])
                y <- x[ii,is]
            } else if (sum(is.element(tolower(it),names(season.abb())))>0) {
                ii <- is.element(months,eval(parse(text=paste('season.abb()$',it,sep=''))))
                y <- x[ii,is]
            } else if (is.element("annual",cls)) {
                ii <- is.element(years,seq(year(it)[1],year(it)[2],1))
                y <- x[ii,is]
            } else if (inherits(it,"Date")) {
                if ( length(it) == 2 ) {
                    if (verbose) print('Between two dates')
                    if (verbose) print(it)          
                    if (is.element("month",cls)) ## it is a month or season
                        it <- seq(it[1],it[2],by='month')
                    else if (is.element("day",cls)) ## it is a day
                        it <- seq(it[1],it[2],by='day')
                    else if (is.element("annual",cls))  ## it is a year
                        it <- seq(it[1],it[2],by='year')

                    ii <- is.element(t,it)                   
                }
                y <- x[ii,is]
            }
        }
        else if (sum(is.element(it,1600:2200)) > 0) {
            if (length(it)==2) ii <- is.element(years,min(it):max(it)) else
            ii <- is.element(years,it)
                                        #print('years')
                                        #print(paste("Number of matches=",sum(ii)))
                                        #print(years[ii]); print(it)
            y <- x[ii,is]
        } 
        else if (is.element('month',cls) & (max(it) <= 12) ) y <- x[months==it,]
        else if (is.element('season',cls)) {
            it <- switch(it,'1'=1,'2'=4,'3'=7,'4'=10)
            y <- x[months==it,]
        }
        else if (is.element('annual',cls)) y <- x[years==it,]
        else if ( (min(it) > 0) & (max(it) <= length(index(x))) ) y <- x[it,] 
        else if (is.character(it)) {
                                        #print("Dates")
            y <- matchdate(x,it)
        } else
            if (inherits(it,c('field','station','zoo'))) {
                                        # Match the times of another esd-data object
                                        # print('field/station')
                y <- matchdate(x,it)
            }
        d[3] <- length(index(y))
        class(y) <- cls
        d -> attr(y,'dimensions')
        y <- attrcp(x,y)
                                        #nattr <- softattr(x)
                                        #for (i in 1:length(nattr))
                                        #  attr(y,nattr[i]) <- attr(x,nattr[i])
                                        # Need to assign the changes to x to make an effect on subsequent spatial subsetting. 
        x <- y
    }

    if (!is.null(is)) {
                                        #print("Sub-region"); print(is)        
                                        #  if (sum(is.element(dimension,"space"))) {
        if ( (is.list(is)) | (is.data.frame(is)) ) {
            if ( (is.null(is[[1]])) | (sum(is.finite(is[[1]])) < 2) ) is[[1]] <- c(-180,360)
            if ( (is.null(is[[2]])) | (sum(is.finite(is[[2]])) < 2) ) is[[2]] <- c(-90,90)
            nam <- names(is)
            if (is.null(nam)) nam <- c("lon","lat")
            if ( (length(is[[1]])==2) & (length(is[[2]])==2) &
                (nam[1]=="lon") & (nam[2]=="lat") ) {
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
                                        #print(d)
            } else {
                                        # Select the specific grid point values from the longitude
                                        # and latitude coordinates
                lon.pick <- is[[1]]; lat.pick <- is[[2]]
                d <- attr(x,'dimensions')
                print(paste("Select ",length(lon.pick),"x",length(lat.pick),
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
            }
        }
    } 
    
    class(x) <- cls; class(y) <- cls
    y <- attrcp(x,y,ignore=c('longitude','latitude'))
                                        #nattr <- softattr(x,ignore=c('longitude','latitude'))
                                        #for (i in 1:length(nattr))
                                        #  attr(y,nattr[i]) <- attr(x,nattr[i])
                                        #mostattributes(y) <- attributes(x)
                                        #attr(y,'date-stamp') <- date()
                                        #attr(y,'call') <- match.call()
    attr(y,'history') <- history.stamp(x)
    return(y)
}

subset.comb <- function(x,it=NULL,is=NULL) {
                                        #print("subset.comb")
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

subset.eof <- function(x,pattern=NULL,it=NULL,is=NULL,verbose=FALSE) {
                                        #print("subset.eof")
    if (is.null(is) & is.null(it) & is.null(pattern)) return(x)                                    
    if (is.null(it) & is.null(is[1]) & is.null(is[2]) & is.null(pattern)) return(x) 
    d <- dim(x); greenwich <- TRUE
    clim <- attr(x,'mean')
    
    # Pattern extracts certain modes/patterns
    if (!is.null(pattern)) {
      y <- x[,pattern]
      y <- attrcp(x,y)
      class(y) <- class(x)
      attr(y,'eigenvalues') <- attr(y,'eigenvalues')[pattern]
      attr(y,'pattern') <- attr(y,'pattern')[,,pattern]
      if (!is.null(attr(x,'n.apps'))) {
        attr(y,'n.apps') <- attr(x,'n.apps')
        attr(y,'appendix.1') <- attr(x,'appendix.1')
      }
      x <- y
    }
    
    if (is.null(is)) is <- 1:d[length(d)] else

    if (is.list(is)) {
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
        it <- 1:d[1] 
        dates <- index(x)
        keep <- 1 : d[1]
    } else {
        if (is.numeric(it)) {
                                        #print("numeric it")
            if (max(nchar(as.character(it)))<=2)
                keep <- is.element(as.numeric(format(index(x),"%m")),it)
            else {
                it <- seq(it[1],it[2],1)
                keep <- is.element(as.numeric(format(index(x),"%Y")),it)      
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
    y <- x[keep,is]
    
    class(x) <- cls; class(y) <- cls
    y <- attrcp(x,y,ignore=c('greenwich','mean'))
    attr(y,'greenwich') <- greenwich
    clim -> attr(y,'mean')
    ## browser()
    attr(y,'pattern') <- attr(x,"pattern")[,,is]
    attr(y,'eigenvalues') <-attr(x,"eigenvalues")[is]
                                        #attr(y,'date-stamp') <- date()
                                        #attr(y,'call') <- match.call()
    attr(y, "dimensions") <- dim(attr(x,"pattern")[,,is])
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

subset.pattern <- function(x,is) {
    if (is.list(is)) {
        y <- attr(x,'pattern')
        lons <- attr(x,'longitude')
        lats <- attr(x,'latitude')
        nms <- substr(tolower(names(is)),1,3)
        if (sum(is.element(nms,'lon'))>0)
            ix <- (lons >= min(is[[is.element(nms,'lon')]])) &
                (lons >= max(is[[is.element(nms,'lon')]])) else
        ix <- is.finite(lons)
        if (sum(is.element(nms,'lat'))>0)
            iy <- (lons >= min(is[[is.element(nms,'lat')]])) &
                (lons >= max(is[[is.element(nms,'lat')]])) else
        iy <- is.finite(lats)
        y[ix,iy] -> attr(x,'pattern')
        lons -> attr(x,'longitude')
        lats -> attr(x,'latitude')
    }
    return(x)
}

subset.pca <- function(x,pattern=NULL,it=NULL,is=NULL) {
  print('subset.pca')
  if (!is.null(pattern)) {
    y <- x[,pattern]
    y <- attrcp(x,y)
    class(y) <- class(x)
    attr(y,'eigenvalues') <- attr(y,'eigenvalues')[pattern]
    attr(y,'pattern') <- attr(y,'pattern')[,pattern]
    if (!is.null(attr(x,'n.apps'))) {
      attr(y,'n.apps') <- attr(x,'n.apps')
      attr(y,'appendix.1') <- attr(x,'appendix.1')
    }
  } else y <- x
  #browser()
  return(y)
}


subset.corfield <- function(x,it=NULL,is=NULL) {
    x
}

subset.ds <- function(x,it=NULL,is=NULL) {
    y <- x
    if (!is.null(it)) {
    }
    if (!is.null(is))  {
        x <- subset.pattern(x)
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

subset.dsensemble <- function(x,it=NULL,is=NULL) {
    if (is.null(it) & is.null(is)) return(x)
                                        #print("subset.dsensemble"); print(it)
    x0 <- x
    d <- dim(x)
    if (is.null(is)) is <- 1:d[2]
    if (!is.null(it)) {
        if (is.character(it)) it <- tolower(it)
        if ( (length(rownames(table(month(x))))==1) & (it==0) )
            return(x)
                                        # Different ways of selecting along the time dimension
        if ( inherits(it[1],"logical") & (length(it)==length(x)) )
            y <- x[it,is] else
        if (it[1]==0) {
                                        #print("Annual means")
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
                                        #print(yr)
            i1 <- is.element(yr1,yr)
            i2 <- is.element(yr2,yr)
            i3 <- is.element(yr3,yr)
            i4 <- is.element(yr4,yr)
                                        #print(c(sum(i1),sum(i2),sum(i3),sum(i4)))
            y <- zoo(0.25*(coredata(djf[i1,]) +
                           coredata(mam[i2,]) +
                           coredata(jja[i3,]) +
                           coredata(son[i4,])),
                     order.by=as.Date(paste(yr,'01-01',sep='-')))
            y <- attrcp(x0,y)
            class(y) <- class(x0)
        } else if (is.character(it)) {
                                        #print("it is character")
            months <- month(x)
            if (sum(is.element(tolower(substr(it,1,3)),tolower(month.abb)))>0) {
                ii <- is.element(months,(1:12)[is.element(tolower(month.abb),
                                                          tolower(substr(it,1,3)))])
                y <- x[ii,is]
            } else if (sum(is.element(tolower(it),names(season.abb())))>0) {
                                        #print("season")
                mon <- eval(parse(text=paste('season.abb()$',it,sep='')))
                                        #print(mon)
                ii <- is.element(months,mon)
                y <- x[ii,is]
                                        #      }
                                        #    if ( (min(it) > 0) & (max(it) < 5) ) {
                                        #      # dsensemble is not for monthly data
                                        #      #print("Seasonal selection")
                                        #      ii <- is.element(month(x),c(1,4,7,10)[it])
                                        #      y <- x[ii,is]
            } } else
                if (sum(is.element(it,1600:2200)) > 0) {
                                        #print("it contains year(s)")
                    ii <- is.element(year(x),it)
                                        #print(paste("Number of matches=",sum(ii)))
                    y <- x[ii,is]
                } else if (is.character(it)) {
                                        #print("Dates")
                    x <- matchdate(x,it)
                } else if (inherits(is,c('field','station'))) {
                                        # Match the times of another esd-data object
                    x <- matchdate(x,it)
                }
                                        #print("housekeeping")
        d[3] <- length(index(y))
        class(y) <- class(x0)
        d -> attr(y,'dimensions')
                                        #str(y)
        y <- attrcp(x,y,ignore='station')
                                        #print("station"); str(attr(x,'station')); print(class(attr(x,'station')))
                                        #plot(subset(attr(x,'station'),it=it))
                                        #browser()
        if ( (it!=0) & (!inherits(attr(x,'station'),'annual')) )
            attr(y,'station') <- subset(attr(x,'station'),it=it) else
        attr(y,'station') <- annual(attr(x,'station'))
                                        #print("HERE")
                                        #nattr <- softattr(x)
                                        #for (i in 1:length(nattr))
                                        #  attr(y,nattr[i]) <- attr(x,nattr[i])
    } else {
                                        # Bug-fix AM: 2014-09-22
        y <- x[,is]
        y <- attrcp(x,y)
        attr(y, "model_id") <- attr(x, "model_id")[is]
        attr(y, "scorestats") <- attr(x, "scorestats")[is]
                                        #browser()
        attr(y,"history") <- history.stamp(x)
        if (length(is)==1) class(y) <- c("ds","zoo") else class(y) <- class(x)
    }
                                        #print("exit")
    return(y)  
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

