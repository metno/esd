
## Author Rasmus E. Benestad - was initially part of subset.R file
## Modified by A. Mezghani
## Last update 06.01.2014 ; 24-02-2014

subset.station <- function(x,it = NULL,is=NULL,loc=NULL , param = NULL,
                           stid = NULL ,lon = NULL, lat = NULL, 
                           alt = NULL, cntr = NULL, src = NULL , nmin = NULL,
                           verbose=FALSE) {
    
    ## 
    if (inherits(it,c('field','station','zoo'))) {
        ## Match the times of another esd-data object
        print('field/station')
        x2 <- matchdate(x,it)
        return(x2)
    }
    ##print("subset.station")
    if (is.null(dim(x))) {
        x2 <- station.subset(x,it=it,is=1,verbose=verbose)
    } else {
        ##print("here")
        x2 <- station.subset(x,it=it,is=is,verbose=verbose)
        ## 
        ## extra selection based on meta data
        ## ss <- select.station(x=x2,loc = loc , param = param,  stid = stid ,lon = lon, lat = lat, alt = alt, cntr = cntr, src = src , nmin = nmin)
        ## 
        ## if (!is.null(ss)) {
        ##    id <- is.element(attr(x2,'station_id'),ss$station_id)
        ## Keep selected stations only
        ##    x2 <- station.subset(x2,it=it,is=which(id),verbose=verbose)
        ##}
        ##if (!is.null(is)) x2 <- station.subset(x2,it=it,is=is,verbose=verbose)
    }
    return(x2)
}

station.subset <- function(x,it=NULL,is=NULL,verbose=FALSE) {

    ## REB: Use select.station to condition the selection index is...
    ## loc - selection by names
    ## lon/lat selection be geography or closest if one coordinate lon/lat
    ##         if two-element vectors, define a region
    ## alt - positive values: any above; negative any below height
    ## cntr - selection by country
    ## 
    
    nval <- function(x) sum(is.finite(x))
    x0 <- x
    if (is.null(it) & is.null(is)) return(x)
    ## 
    d <- dim(x)
    if (is.null(d)) {
        if (verbose)
            print("Warning : One dimensional vector has been found in the coredata")
        x <- zoo(as.matrix(coredata(x)),order.by=index(x))
        x <- attrcp(x0,x)
        class(x) <- class(x0)
    } 
    d <- dim(x)
    if (is.null(is)) is <- 1:d[2]
    if (is.null(it)) it <- 1:d[1]
    
    ## 
    ##print("HERE")
    ## get time in t
    t <- index(x)
    ii <- is.finite(t)
    if (verbose) print(it)                            

    ##  if (datetype=="Date") {
    if (inherits(t,c("Date","yearmon"))) {
       if (verbose) print('years ++')
        ## REB: replaced by lines below:
        ##    year <- as.numeric( format(t, '%Y') ) 
        ##    month <- as.numeric( format(t, '%m') )
        yr <- year(x)
        mo <- month(x)
        dy <- day(x)
    } else if (inherits(t,c("numeric","integer"))) {
        if (verbose) print('years')
        yr <- t
        mo <- dy <- rep(1,length(t))
    } else print("Index of x should be a Date, yearmon, or numeric object")
    
    ## Generate sequence of days, months or years if range of it value is given
    if (is.character(it)) {
        if (verbose) print('it is character')
        if ((levels(factor(nchar(it)))==10)) ##  convert as Date
            it <- as.Date(it)
        if ( length(it) == 2 ) {
            if (verbose) print('Between two dates')
            ##if (nchar(it[1])==4) it[1] <- paste(it[1],'-01-01',sep='')
            ##if (nchar(it[2])==4) it[2] <- paste(it[1],'-12-31',sep='')
            if (verbose) print(it)          
            if (inherits(x,"month")) ## it is a month or season
                it <- seq(it[1],it[2],by='month')
            else if (inherits(x,"day")) ## it is a day
                it <- seq(it[1],it[2],by='day')
            else if (inherits(x,"annual")) ## it is a year
                it <- seq(it[1],it[2],by='year')
            ii <- is.element(t,it)
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
                print('it is a Date object')
                ii <- is.element(t,it)
            } else {
                str(it); print(class(it))
                ii <- rep(FALSE,length(t))
                warning("subset.station: did not recognise the selection citerion for 'it'")
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
         if ((length(nlev)==1)) {
            if (nlev==4) {
                if (verbose) print("it are most probably years")
                if (length(it)==2)
                    ii <- is.element(yr,it[1]:it[2])
                # if it is years:
                else if (min(it)> length(it)) {
                    if (verbose) print("match years")
                    ii <- is.element(yr,it)
                  } 
            } else if (nlev<=4) {
                if (verbose) print("it are most probably seasons")
                if (inherits(x,'season') & (length(it)==1)) {
                    if (verbose)  print(paste("The 'it' value must be a season index between 1 and 4.",
                                              "If not please use character strings instead. e.g. it='djf'"))
                    it <- switch(it,'1'=1,'2'=4,'3'=7,'4'=10,'djf'=1,'mam'=4,'jja'=7,'son'=10)
                    ii <- is.element(mo,it)
                 } else if (inherits(x,'month') | (inherits(x,'day'))) {
                     if (verbose)
                         print("The 'it' value must be a month index. If not please use character strings instead")
                     ii <- is.element(mo,it)
                 }  else {
                    if (verbose) print("it represents indices")
                    ii <- it
                }
            } else if (nlev<=12) {
                if (verbose)
                         print("The 'it' value are most probably a month index. If not please use character strings instead")
                ii <- is.element(mo,it)
            }        
        } else {
            #  length(nlev) > 1
            if (verbose)  print("it most probably holds indices")
            ii <- it
        }
    } else if (inherits(it,c("Date","yearmon"))) {       
        ##        ii <- is.element(t,it)
        if (verbose) print('it is a date object')
        ii <- (t >= min(it)) & (t <= max(it))
    } else {
        ii <- rep(FALSE,length(t))
        warning("subset.station: did not reckognise the selection citerion for 'it'")
    } 
    
    ## it <- (1:length(t))[ii]
    ## 

    class(x) -> cls
    ##print(cls)
    ## update the class of x
    class(x) <- "zoo" 
   
                                        # REB 11.04.2014: is can be a list to select region or according to other criterion
    if (inherits(is,'list')) {
        n <- dim(x)[2]
        selx <- rep(TRUE,n); sely <- selx; selz <- selx
        selc <- selx; seli <- selx; selm <- selx; salt <- selx
        selp <- selx; selF <- selx ; sell <- selx
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
        ##
        ## Need to make sure both it and is are same type: here integers for index rather than logical
        ## otherwise the subindexing results in an empty object
    }

    y <- x[ii,is]
    #if (is.logical(is))
    #    is <- (1:length(is))[is]
    ##else 
    ##    is <- is.element(1:d[2],is)
    
    class(x) <- cls; class(y) <- cls
    y <- attrcp(x,y,ignore=c("names"))
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
        attr(y,'unit') <- attr(x,'unit')[is]
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
    ##attr(y,'date-stamp') <- date()
    ##attr(y,'call') <- match.call()
    attr(y,'history') <- history.stamp(x)   
    if (inherits(y,"annual")) index(y) <- as.numeric(year(index(y)))
    return(y)
}
    
