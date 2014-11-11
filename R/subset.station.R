
## Author Rasmus E. Benestad - was initially part of subset.R file
## Modified by A. Mezghani
## Last update 06.01.2014 ; 24-02-2014
 
subset.station <- function(x,it = NULL,is=NULL,loc=NULL , param = NULL,
                           stid = NULL ,lon = NULL, lat = NULL, 
                           alt = NULL, cntr = NULL, src = NULL , nmin = NULL,
                           verbose=FALSE) {
  ## browser()
  if (inherits(it,c('field','station','zoo'))) {
    ## Match the times of another esd-data object
    print('field/station')
    x2 <- matchdate(x,it)
    return(x2)
  }
  ##print("subset.station")
  if (is.null(dim(x))) {
      x2 <- station.subset(x,it=it,is=1)
  } else {
      ##print("here")
      x2 <- station.subset(x,it=it,is=is)
      ## browser()
      ## extra selection based on meta data
      ## ss <- select.station(x=x2,loc = loc , param = param,  stid = stid ,lon = lon, lat = lat, alt = alt, cntr = cntr, src = src , nmin = nmin)
      ## browser()
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

    nval <- function(x) sum(is.finite(x))
    
    x0 <- x
    if (is.null(it) & is.null(is)) return(x)
    d <- dim(x)
    if (is.null(d)) d <- c(length(x),1)
    if (is.null(it)) it <- index(x)
    if (is.null(is)) is <- 1:d[2]

    ## browser()
    ##print("HERE")
    ## get time in t
    t <- index(x)
    #if (class(it)!=class(t)) print("Index and it class do not match !")
    
    ##  if (datetype=="Date") {
    if (inherits(t,c("Date","yearmon"))) {
        ## REB: replaced by lines below:
        ##    year <- as.numeric( format(t, '%Y') ) 
        ##    month <- as.numeric( format(t, '%m') )
        yr <- year(x)
        mo <- month(x)
        dy <- day(x)
    } else if (inherits(t,c("numeric","integer"))) {
        yr <- t
        mo <- rep(1,length(t)); dy <- mo
    } else print("Index of x should be a Date, yearmon, or numeric object")
    
    ## Generate sequence of days, months or years if range of it value is given
    if ((length(it)>2) & (is.character(it)))
        it <- as.Date(it)
    else if ( length(it) == 2 ) {
        if (is.character(it)) {
            if (inherits(x,"month")) ## it is a month or season
                it <- seq(as.Date(it[1]),as.Date(it[2]),by='month')
            else if (inherits(x,"day")) ## it is a day
                it <- seq(as.Date(it[1]),as.Date(it[2]),by='day')
            else if (inherits(x,"annual")) ## it is a year
                it <- seq(as.Date(it[1]),as.Date(it[2]),by='year')
          } else if ((class(it)=="numeric") | (class(it)=="integer")) {
                if (min(it) > 1500) ## it is a year
                    it <- seq(it[1],it[2],by=1)
                #print("HERE"); print(it)
            }
      }
  
    ## browser()
    ## get the subset indices in ii
    if ((class(it)=="numeric") | (class(it)=="integer")) {    
        if ( ((min(it,na.rm=TRUE) > 0) & (max(it,na.rm=TRUE) < 13)) &
             (inherits(x,"month") | inherits(x,"season")) ) {## it is a month or season
          # REB 23.04.14: need to handle monthly and seasonal object differently
            if (inherits(x,"month")) it.mo <- it else
            if (inherits(x,"season")) it.mo <- c(1,4,7,10)[it]
            ii <- is.element(mo,it.mo)
        }  #
        else if ((min(it,na.rm=TRUE) > 0) & (max(it,na.rm=TRUE) < 31)) {## it is a day
            it.dy <- it
            ii <- is.element(dy,it.dy)
        }
        else if (min(it) > 1500) {## it is a year
            it.yr <- it
            ii <- is.element(yr,it.yr)
        }
    }
    else if (inherits(it,c("Date","yearmon"))) {       
#        ii <- is.element(t,it)
        ii <- (t >= min(it)) & (t <= max(it))
    }
    else if (is.character(it)) { ## added AM 10-06-2014
        if (sum(is.element(tolower(substr(it,1,3)),tolower(month.abb)))>0) {
            ii <- is.element(month(x),(1:12)[is.element(tolower(month.abb),tolower(substr(it,1,3)))])
            y <- x[ii,is]
        } else
            if (sum(is.element(tolower(it),names(season.abb())))>0) {
                if (verbose) print("Seasonally selected")
                if (verbose) print(table(month(x)))
                if (verbose) print(eval(parse(text=paste('season.abb()$',it,sep=''))))
                ii <- is.element(month(x),eval(parse(text=paste('season.abb()$',it,sep=''))))
                y <- x[ii,is]
            }
    }
    else ## keep all values
        ii <- 1:length(t)
    
    class(x) -> cls
    ##print(cls)
    ## update the class of x
    class(x) <- "zoo" 
    
    # REB 11.04.2014: is can be a list to select region or according to other criterioe
    if (inherits(is,'list')) {
      n <- dim(x)[2]
      selx <- rep(TRUE,n); sely <- selx; selz <- selx
      selc <- selx; seli <- selx; selm <- selx; salt <- selx
      selp <- selx; selF <- selx
      nms <- names(is)
      ix <- grep('lon',tolower(substr(nms,1,3)))
      iy <- grep('lat',tolower(substr(nms,1,3)))
      #print(nms); print(c(ix,iy))
      iz <- grep('alt',tolower(substr(nms,1,3)))
      ic <- grep('cntr',tolower(substr(nms,1,3)))
      im <- grep('nmin',tolower(substr(nms,1,3)))
      ip <- grep('param',tolower(substr(nms,1,3)))
      id <- grep('stid',tolower(substr(nms,1,3)))
      iF <- grep('FUN',substr(nms,1,3))
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
      if (length(slon)==2) selx <- (lon(x) >= min(slon)) & (lon(x) <= max(slon))
      if (length(slat)==2) sely <- (lat(x) >= min(slat)) & (lat(x) <= max(slat))
      #browser()
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
      ## browser()
      is <- selx & sely & selz & selc & seli & selm & selp & selF
      #print(c(length(is),sum(is),sum(selx),sum(sely)))
    }
    
    #else if (inherits(x0,c("month"))) {
    #    ii <- is.element(mo,it.mo)
    #    ##print("here")
    #}
    #else if ((inherits(x0,c("day")))) {
        ##browser()
    #    if (verbose) print("Daily Object")
    #    ii <- is.element(dy,it.dy)
    #}
    #else if ((min(it) > 0) & (max(it) < 5) & (inherits(x0,"season"))) {
    #    if (verbose) print("Seasonal Object")
    #    ii <- is.element(mo,c(1,4,7,10)[it.mo]) 
    #}
    #else if ( (min(it) > 0) & (inherits(x0,"annual")) ) {
    #    if (verbose) print("Annual Object")
    #    ## browser()
    #    ii <- is.element(it,it.yr)
    #    ##if (sum(ii)>0) y <- x[ii,is] else stop(paste("it value(s) must be in the range of",paste(range(yr),collapse="-")))
    # }
    ## else if ((min(it) > 0) & (max(it) < 13) & (sum(is.element(it,1:12)) > 0)) {
        
    ## }
    #else if (sum(is.element(it,1600:2200)) > 0) {
    #    ##print(it)
    #    if (length(it)==2) it <- it[1]:it[2]
    #    ii <- is.element(yr,it.yr)
    #   
    #}
    #else if (is.character(it)) {
    #    dates <- as.Date(it)
    #    ii <- is.element(index(x),dates)
    #   
    #}
    #else
    #    ii <- 1:length(t)
    
    #browser()
    y <- x[ii,is]

    class(x) <- cls; class(y) <- cls
    y <- attrcp(x,y,ignore=c("names"))
    ## nattr <- softattr(x)
    ## for (i in 1:length(nattr))
    ##  attr(y,nattr[i]) <- attr(x,nattr[i])
    ## mostattributes(y) <- attributes(x)
    attr(y,'longitude') <- attr(x,'longitude')[is]
    attr(y,'latitude') <- attr(x,'latitude')[is]
    if (length(alt(x)==length(is)) attr(y,'altitude') <- attr(x,'altitude')[is]
    if (length(cntr(x)==length(is)) attr(y,'country') <- attr(x,'country')[is]
    if (length(src(x)==length(is)) attr(y,'source') <- attr(x,'source')[is]
    if (length(stid(x)==length(is)) attr(y,'station_id') <- attr(x,'station_id')[is]
    if (length(loc(x)==length(is)) attr(y,'location') <- attr(x,'location')[is]
    if (length(attr(x,'quality'))==length(is)) attr(y,'quality') <- attr(x,'quality')[is]
    #attr(y,'URL') <- attr(x,'URL')[is]
    #attr(y,'history') <- attr(x,'history')[is]
    if (length(varid(x)==length(is)) attr(y,'variable') <- attr(x,'variable')[is]
    #attr(y,'element') <- attr(x,'element')[is]
    if (length(attr(x,'aspect'))==length(is)) attr(y,'aspect') <- attr(x,'aspect')[is]
    if (length(unit(x)==length(is)) attr(y,'unit') <- attr(x,'unit')[is]
    if (length(attr(x,'longname'))==length(is)) attr(y,'longname') <- attr(x,'longname')[is]
    if (length(attr(x,'reference'))==length(is)) attr(y,'reference') <- attr(x,'reference')[is]
    if (length(attr(x,'info'))==length(is)) attr(y,'info') <- attr(x,'info')[is]
    if (length(attr(x,'method'))==length(is)) attr(y,'method') <- attr(x,'method')[is]
    if (length(attr(x,'type'))==length(is)) attr(y,'type') <- attr(x,'type')[is]
    ##attr(y,'date-stamp') <- date()
    ##attr(y,'call') <- match.call()
    attr(y,'history') <- history.stamp(y)
    return(y)
}
 
