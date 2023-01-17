# Documentation in subset.R
#' @exportS3Method
#' @export subset.station
subset.station <- function(x, it=NULL, is=NULL, loc=NULL, param=NULL,
                           stid=NULL, lon=NULL, lat=NULL, alt=NULL, cntr=NULL,
			   src=NULL, nmin=NULL, verbose=FALSE, ...) {
    
    ##
    if (verbose) print(match.call())
    if (is.null(attr(x,'unit'))) attr(x,'unit') <- NA
    if (verbose) print(paste("subset.station:",c(varid(x),esd::unit(x))))
    d <- dim(x)
    if (is.null(d)) d <- c(length(x),1)
    ## REB 2022-03-31
    ## Make sure that this algorithm keeps track of the type of data (precip, temp, etc)
    if ( (length(esd::unit(x))!=d[2]) | (length(esd::varid(x))!=d[2]) ) {  
      attr(x,'unit') <- rep(esd::unit(x)[1],d[2])
      attr(x,'variable') <- rep(esd::varid(x)[1],d[2])
    }
    ## subset times:
    if (inherits(it,c('field','station','zoo'))) {
        ## Match the times of another esd-data object
        if (verbose) print('subset.station: field/station')
        x2 <- matchdate(x,it)
        x2 <- attrcp(x,x2)
        attr(x2,'history') <- history.stamp()
        return(x2)
    }
    ## subset space:
    if (inherits(is,c('field','station','zoo'))) {
        ## Match the space index of another esd-data object
        if (verbose) print('subset.station: is: field/station')
        x2 <- subset.station(x,is=loc(is))
        return(x2)
    }
    if (is.character(is) | is.character(loc)) {
      if (verbose) print('subset.station: search on location names')
      ## search on location name: foce lower case for all
      locs <- tolower(loc(x))
      if (is.null(loc)) is <- tolower(is) else
                        is <- tolower(loc)
      #mcl <- min(c(nchar(is),nchar(locs)))
      ## Set same lengths for all; minimum character
      #locs <- substr(locs,1,mcl)
      #is <- substr(is,1,mcl)
      #illoc <- is.element(locs,is)
      ns <- length(is)
      for (ii in 1:ns) illoc <- grep(is[ii],locs, ignore.case = TRUE)
      x2 <- subset(x,it=it,is=illoc,verbose=verbose)
      if (verbose) {
        print(paste("subset.station:",is))
        print(paste("subset.station:",loc(x2)))
      }
      return(x2)
    }
    if (is.null(dim(x))) {
        x2 <- station.subset(x,it=it,verbose=verbose)
    } else {
        ##print("subset.station: here")
        x2 <- station.subset(x,it=it,is=is,verbose=verbose)
    }
    ## Check if there is only one series but if the dimension 
    if ( (!is.null(d)) & is.null(dim(x2)) ) 
      if (d[2]==1) dim(x2) <- c(length(x2),1)
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
  if (verbose) {
    print(match.call())
    print(paste("station.subset:",c(varid(x),esd::unit(x))))
  }
  nval <- function(x) sum(is.finite(x))
  x0 <- x
  if (is.null(it) & is.null(is)) return(x)
  
  ## Check whether matrix of vector
  d <- dim(x)
  if (is.null(d)) {
    if (verbose)
      print("station.subset: Warning : One dimensional vector has been found in the coredata")
    x <- zoo(as.matrix(coredata(x)),order.by=index(x))
    x <- attrcp(x0,x)
    class(x) <- class(x0)
  } 
  d <- dim(x)
  
  if (is.null(is)) is <- 1:d[2]
  if (is.null(it)) it <- 1:d[1]
  #if (is.logical(it)) it <- (1:d[1])[it]
  if (is.logical(is)) is <- (1:d[2])[is]
  
  ## get time in t
  t <- index(x)
  ii <- is.finite(t)
  
  if (verbose) print('station.subset: it - temporal indexing')
  if (verbose) print(it)
  
  if (inherits(t,c("Date","yearmon"))) {
    if (verbose) print('station.subset: years ++')
    yr <- year(x)
    mo <- month(x)
    dy <- day(x)
  } else if (inherits(t,c("numeric","integer"))) {
    if (verbose) print('station.subset: years')
    yr <- t
    mo <- dy <- rep(1,length(t))
  } else print("station.subset: Index of x should be a Date, yearmon, or numeric object")
  
  if(is.character(it)) {
    if ((levels(factor(nchar(it)))==10)) it <- as.Date(it)
  }
  
  if(inherits(it,c("POSIXt"))) it <- as.Date(it)
  if(inherits(it,c("Date"))) {
    if (inherits(t,"yearmon")) t <- as.Date(t)
    if ( length(it) == 2 ) {
      if (verbose) print('station.subset: Between two dates')
      if (verbose) print(it)
      ii <- (t >= min(it)) & (t <= max(it))
    } else {
      ii <- is.element(t,it)
    }
  } else if(inherits(it,"yearmon")) {
    ii <- is.element(as.yearmon(t),it)
  } else if (is.character(it)) {
    if (verbose) print('station.subset: it is character')
    if (sum(is.element(tolower(substr(it,1,3)),tolower(month.abb)))>0) {
      if (verbose) print('station.subset: Monthly selected')
      if (is.seasonal(x)) {
        it <- gsub('Dec', 'Jan', it, ignore.case=TRUE)
        it <- gsub('Feb', 'Jan', it, ignore.case=TRUE)
        it <- gsub('Mar', 'Apr', it, ignore.case=TRUE)
        it <- gsub('May', 'Apr', it, ignore.case=TRUE)
        it <- gsub('Jun', 'Jul', it, ignore.case=TRUE)
        it <- gsub('Aug', 'Jul', it, ignore.case=TRUE)
        it <- gsub('Sep', 'Oct', it, ignore.case=TRUE)
        it <- gsub('Nov', 'Oct', it, ignore.case=TRUE)
      }
      ii <- is.element(month(x),(1:12)[is.element(tolower(month.abb),tolower(substr(it,1,3)))])
    } else if (sum(is.element(tolower(it),names(season.abb())))>0) {
      if (verbose) print("station.subset: Seasonally selected")
      if (verbose) print(table(month(x)))
      if (verbose) print(eval(parse(text=paste('season.abb()$',it,sep=''))))
      ii <- is.element(month(x),eval(parse(text=paste('season.abb()$',it,sep=''))))
    }
  } else if (inherits(it,"Date")) {
    if (verbose) print('station.subset: it is a Date object')
    ii <- is.element(t,it)
  } else if (is.logical(it)) {
    ii <- (1:d[1])[it]
  } else if ((class(it)=="numeric") | (class(it)=="integer")) {
    if (verbose) print('station.subset: it is numeric or integer')
    nlev <- as.numeric(levels(factor(nchar(it)))) # REB bug        
    # nchar returns the string length, but these lines need to find the number of different levels/categories
    # AM 2015-02-16 DO not agree    nlev <- as.numeric(levels(factor(as.character(it)))) # REB 2015-01-15
    if (verbose) {print(nlev); print(it)}
    if ((length(nlev)==1)) {
      if (nlev==4) {
        if (verbose) print("station.subset: it are most probably years")
        if (length(it)==2) {
          ii <- is.element(yr,it[1]:it[2])
          if (verbose) print(paste('Subset of',sum(ii),'data points between',
                                   min(yr),'-',max(yr),'total:',length(yr)))
          # if it is years:
        } else if (min(it)> length(it)) {
          if (verbose) print("station.subset: match years")
          ii <- is.element(yr,it)
        } 
      } else if (nlev<=4) {
        if (verbose) print("station.subset: it are most probably seasons")
        if (inherits(x,'season') & (length(it)==1)) {
          if (verbose) print(paste("station.subset: The 'it' value must be a season index between 1 and 4.",
                                   "If not please use character strings instead. e.g. it='djf'"))
          it <- switch(tolower(it),'1'=1,'2'=4,'3'=7,'4'=10,'djf'=1,'mam'=4,'jja'=7,'son'=10)
          ii <- is.element(mo,it)
        } else if ( (inherits(x,'month') | (inherits(x,'day'))) &
                    ( (max(it) <= 12) & (min(it) >= 1) ) ) {
          if (verbose) {
            print(paste("station.subset: The 'it' value must be a month index.",
                        "If not please use character strings instead"))
            print(range(it))
          }
          ii <- is.element(mo,it)
        } else {
          if (verbose) print("station.subset: it represents indices")
          ii <- it
        }
      } else if (nlev<=12  & ( (max(it) <= 12) & (min(it) >= 1) )) {
        if (verbose) {
          print(paste("station.subset: The 'it' value is most probably a month index.",
                      "If not please use character strings instead"))
          print(range(it))
        }
        ii <- is.element(mo,it)
      } else {
        if (verbose) {
          print("station.subset: The 'it' value is most probably an index.")
          print(range(it))
        }
        ii <- it
      }
    } else {
      #  length(nlev) > 1
      if (verbose) print("station.subset: it most probably holds indices")
      ii <- it
    }
  } else {
    ii <- rep(FALSE,length(t))
    warning("station.subset: did not recognise the selection citerion for 'it'")
  }
  
  class(x) -> cls
  ##print(cls)
  ## update the class of x
  class(x) <- "zoo"
  
  if (verbose) print('station.subset: is - spatial indexing')
  ## REB 11.04.2014: is can be a list to select region or according to other criterion
  if (inherits(is,'list')) {
    if (verbose) print("station.subset: 'is' is a list object")
    selx <- rep(TRUE,dim(x)[2]); sely <- selx; selz <- selx
    selc <- selx; seli <- selx; selm <- selx; salt <- selx
    selp <- selx; selF <- selx ; sell <- selx; selj <- selx
    nms <- names(is)
    if (verbose) print(nms)
    ## See which options are provided in the list and get their indices
    il <- grep('loc',tolower(nms), ignore.case = TRUE)
    ix <- grep('lon',tolower(nms))
    iy <- grep('lat',tolower(nms))
    #print(nms); print(c(ix,iy))
    iz <- grep('alt',tolower(nms))
    ic <- grep('cntr',tolower(nms))
    im <- grep('nmin',tolower(nms))
    ip <- grep('param',tolower(nms))
    id <- grep('stid',tolower(nms))
    iF <- grep('FUN',nms)
    ij <- grep('index',tolower(nms))
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
    if (length(ij)>0) sj <- is[[ij]] else sj <- NULL
    #print(slat); print(range(lat(x)))
    if (verbose) {
      print('station.subset:')
      print(sloc); print(slon); print(slat); print(salt); print(scntr); print(smin); print(sparam)
    }
    ## REB 2021-08-19
    #if (length(sloc)>0) sell <- is.element(tolower(loc(x)),tolower(sloc))
    if (length(sloc)>0) sell <- is.element(1:length(loc(x)),grep(sloc,loc(x),ignore.case = TRUE))
    if (length(slon)==2) selx <- as.logical((lon(x) >= min(slon)) & (lon(x) <= max(slon)))
    if (length(slat)==2) sely <- as.logical((lat(x) >= min(slat)) & (lat(x) <= max(slat)))
    if (length(salt)==2) selz <- as.logical((alt(x) >= min(salt)) & (alt(x) <= max(salt)))
    if (length(salt)==1) {
      if (salt < 0) selz <- (alt(x) <= abs(salt)) else
                    selz <- (alt(x) >= salt)
    }
    if (length(scntr)>0) selc <- is.element(tolower(cntr(x)),tolower(scntr))
    if (length(snmin)>0) selm <- apply(coredata(x),2,nval) > snmin
    if (length(sparam)>0) selp <- is.element(tolower(attr(x,"variable")),tolower(sparam))
    if (length(sstid)==2) seli <- (stid(x) >= min(sstid)) & (stid(x) <= max(sstid)) else
      if (length(sstid)>0) seli <- is.element(stid(x),sstid)
    if (length(sFUN)>0) selm <- apply(coredata(x),2,sFUN) # Not quite finished...
    if (length(sparam)>0) selp <- is.element(tolower(attr(x,"variable")),tolower(sparam))
    if(length(sj)>0) selj <- is.element(seq(1,dim(x)[2]),sj)
    ##
    is <- sell & selx & sely & selz & selc & seli & selm & selp & selF & selj
    if (sum(is)==0) {warning('Returning empty station set'); return(NULL)}
    ##
    ## Need to make sure both it and is are same type: here integers for index rather than logical
    ## otherwise the subindexing results in an empty object
    if (verbose) print(paste('station.subset:',sum(is),'locations'))
  } else if(is.numeric(is) | is.integer(is)) {
    if(verbose) print("station.subset: 'is' is numeric or integer.")
    if(all(is<=dim(x)[2])) {
      if(is.null(stid(x))) {
        if(verbose) print("station.subset: is most probably holds indices")
        is <- is.element(seq(1,dim(x)[2]),is)
      } else if (all(is.element(is,stid(x)))) {
        if(verbose) print("station.subset: is most probably holds station id:s")
        is <- is.element(stid(x),is)
      } else {
        if(verbose) print("station.subset: is most probably holds indices")
        is <- is.element(seq(1,dim(x)[2]),is)
      }
    } else {
      if(verbose) print("station.subset: is most probably holds station id:s")
      is <- is.element(stid(x),is)
    }
  } else {
    is <- rep(FALSE,dim(x)[2])
    warning("station.subset: did not recognise the selection citerion for 'is'")
  }
  
  if (verbose) print(paste('station.subset: Subset of',sum(ii),'data points between',
                           min(yr),'-',max(yr),'total:',length(yr),
                           'from',length(is),'locations'))
  
  is[is.na(is)] <- FALSE
  if (is.logical(ii)) ii <- which(ii)
  if (is.logical(is)) is <- which(is)
  
  d <- dim(x)
  if (is.null(d)) {
    warning('station.subset: unsuccessfull subsetting')
    return(x)
  }
  
  if ( (max(ii) <= d[1]) & (max(is) <= d[2]) ) {
    y <- x[ii,is]
  } else { 
    if(verbose) print(is)
    if(verbose) print(dim(x))
    warning('station.subset: unsuccessful subsetting')
    return(x)
  }
  
  if (verbose) print(summary(coredata(y)))
  class(x) <- cls
  class(y) <- cls
  y <- attrcp(x,y,ignore=c("names"))
  for(a in c('location','longitude','latitude','altitude','station_id')) {
    if(a %in% names(attributes(x))) attr(y,a) <- attr(x,a)[is]
  }
  
  if (length(esd::unit(x))== length(x[1,])) attr(y,'unit') <- esd::unit(x)[is] else
                                         attr(y,'unit') <- esd::unit(x)[1]
  if (length(varid(x))== length(x[1,])) attr(y,'variable') <- varid(x)[is] else
                                     attr(y,'variable') <- varid(x)[1]
  if (is.null(attr(y,'unit'))) attr(y,'unit') <- NA
  if (verbose) print(c(varid(x),esd::unit(x)))
  
  if (verbose) print(paste('subset.station: Before subsetting',loc(x)[is],varid(x)[is],
                           esd::unit(x)[is],lon(x)[is],lat(x)[is]))
  if (verbose) print(is)

  attr(y,'longitude') <- attr(x,'longitude')[is]
  attr(y,'latitude') <- attr(x,'latitude')[is]
  
  if (!is.null(attr(y,'altitude')))
    attr(y,'altitude') <- attr(x,'altitude')[is]
  if (!is.null(attr(y,'country')))
    if (length(cntr(x))>1) attr(y,'country') <- cntr(x)[is] else
      attr(y,'country') <- cntr(x)
  if (!is.null(attr(y,'source')))
    if (length(src(x))>1) attr(y,'source') <- src(x)[is] else
      attr(y,'source') <- src(x)
  if (!is.null(attr(y,'station_id')))
    attr(y,'station_id') <- attr(x,'station_id')[is]
  if (!is.null(attr(y,'location')))
    attr(y,'location') <- attr(x,'location')[is]
  if (!is.null(attr(y,'quality')))
    attr(y,'quality') <- attr(x,'quality')[is]
  ## attr(y,'history') <- attr(x,'history')[is]
  if (!is.null(attr(y,'variable')))
    if (length(varid(x))>1) attr(y,'variable') <- varid(x)[is] else
      attr(y,'variable') <- varid(x)
  # attr(y,'element') <- attr(x,'element')[is]
  if (!is.null(attr(y,'aspect')))
    if (length(attr(y,'aspect'))>1) attr(y,'aspect') <- attr(x,'aspect')[is] else
      attr(y,'aspect') <- attr(x,'aspect')
  if (!is.null(attr(y,'unit')))
    if (length(esd::unit(x))>1) attr(y,'unit') <- esd::unit(x)[is] else
      attr(y,'unit') <- esd::unit(x)
  if (!is.null(attr(y,'longname')))
    if (length(attr(y,'longname'))>1) attr(y,'longname') <- attr(x,'longname')[is] else
      attr(y,'longname') <- attr(x,'longname')
  if (!is.null(attr(y,'reference')))
    if (length(attr(y,'reference'))>1) attr(y,'reference') <- attr(x,'reference')[is] else
      attr(y,'reference') <- attr(x,'reference')
  if (!is.null(attr(y,'info')))
    if (length(attr(y,'info'))>1) attr(y,'info') <- attr(x,'info')[is] else
      attr(y,'info') <- attr(x,'info')
  if (!is.null(attr(y,'method')))
    if (length(attr(y,'method'))>1) attr(y,'method') <- attr(x,'method')[is] else
      attr(y,'method') <- attr(x,'method')
  if (!is.null(attr(y,'type')))
    if (length(attr(y,'type'))>1) attr(y,'type') <- attr(x,'type')[is] else
      attr(y,'type') <- attr(x,'type')
  if (!is.null(attr(y,'URL')))
    if (length(attr(y,'URL'))>1) attr(y,'URL') <- attr(x,'URL')[is] else
      attr(y,'URL') <- attr(x,'URL')
  if (!is.null(attr(y,'na')))
    if (length(attr(y,'na'))>1) attr(y,'na') <- attr(x,'na')[is] else
      attr(y,'na') <- attr(x,'na')
  
  if (verbose) print(paste('station.subset: Final:',loc(y),varid(y),esd::unit(y),lon(y),lat(y)))
  if (length(loc(y))==0) warning('station.subset: no location information - loc(y) == 0')
  if (!is.null(err(y))) attr(y,'standard.error') <- err(x)[ii,is]
  ##attr(y,'date-stamp') <- date()
  ##attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  if (inherits(y,"annual")) index(y) <- as.numeric(year(index(y)))
  return(y)
}


# subset.stationmeta <- function(x, it=NULL, is=NULL, ..., verbose=FALSE) {
#   if(verbose) print(paste('subset.stationmeta: ',paste(dim(x),collapse=' - ')))
#   if (is.null(is)) is <- rep(TRUE,dim(x)[1])
#   if (is.list(is)) {
#     i <- rep(TRUE,length(x[[1]]))
#     listnames <- names(is)
#     if (verbose) {print(sum(i)); print(listnames)}
#     if ('lon' %in% listnames) 
#       i <- (lon(x) >= min(is$lon)) & (lon(x) <= max(is$lon)) 
#     if ('lat' %in% listnames) 
#       i <- i & (lat(x) >= min(is$lat)) & (lat(x) <= max(is$lat)) 
#     if ('alt' %in% listnames) 
#       i <- i & (alt(x) >= min(is$alt)) & (alt(x) <= max(is$alt))
#     if ('cntr' %in% listnames) 
#       #i <- i & (is.element(tolower(substr(x$country,1,nchar(is$cntr))),tolower(is$cntr)))
#       i <- i & (1:(length(i)))[grep(is$cntr,x$country,ignore.case=TRUE, ...)]
#     if (verbose) print(table(tolower(substr(x$country,1,nchar(is$cntr)))))
#     if ('src' %in% listnames) 
#       i <- i & (is.element(tolower(substr(x$source,1,nchar(is$src))),tolower(is$src)))
#     if ('param' %in% listnames) 
#       i <- i & (is.element(tolower(substr(x$element,1,nchar(is$param))),esd2ele(is$param)))
#     if ('loc' %in% listnames) 
#       i <- i & (is.element(tolower(substr(x$location,1,nchar(is$loc))),tolower(is$loc)))
#     is <- i
#   } else if (is.numeric(is) | is.integer(is)) 
#     is <- is.element(1:dim(x)[1],is) else if (is.character(is))
#       is <- is.element(tolower(substr(x$location,1,nchar(is))),tolower(is))
#     if (!is.null(it)) {
#       is <- is & (x$start >= min(it)) & (x$end <= max(it))
#     }
#     if (verbose) print(paste('subset of',sum(is),'elements'))
#     is <- (1:(length(i)))[is]
#     if (verbose) print(str(x))
#     y <- as.data.frame(as.matrix(x[is,]))  
#     if (verbose) print(dim(y))
#     class(y) <- class(x)
#     attr(y,'history') <- history.stamp(x)  
#     if (verbose) str(y)
#     return(y)
# }

#' @exportS3Method    
#' @export subset.stationmeta
subset.stationmeta <- function(x, is=NULL, it=NULL, loc=NULL, param=NULL,
                           stid=NULL, lon=NULL, lat=NULL, alt=NULL, cntr=NULL,
                           src=NULL, nmin=NULL, ..., verbose=FALSE) {
  if (verbose) print('subset.stationmeta')
  # > names(x)
  # [1] "station_id" "location"   "country"    "longitude"  "latitude"
  # [6] "altitude"   "element"    "start"      "end"        "source"
  # [11] "variable"   "wmo"        "quality"

  cls <- class(x)
  x <- as.data.frame(x)
  d <- dim(x)
  ii <- rep(TRUE,d[1])
  if (verbose) print(paste('Original number of stations are',sum(ii)))
  if (!is.null(is)) {
    if ( (is.integer(is)) | (is.numeric(is)) ) ii[-is] <- FALSE
    else if (is.list(is)) {
      listnames <- names(is)
      for (nm in listnames) eval(parse(text=paste0(nm,' <- is$',nm)))
    }
  } 

  if (!is.null(loc)) ii[!is.element(tolower(substr(x$location,1,min(nchar(loc)))),tolower(substr(loc,1,min(nchar(loc)))))] <- FALSE
  if (!is.null(param)) {
    ele <- switch(tolower(param),'tmax'=111,'tmin'=121,'precip'=601,'')
    ii[!is.element(tolower(x$element),tolower(ele))] <- FALSE
  }
  if (!is.null(stid)) ii[!is.element(x$station_id,stid)] <- FALSE
  if (!is.null(lon)) ii[(x$longitude < lon[1]) | (x$longitude > lon[2])] <- FALSE
  if (!is.null(lat)) ii[(x$latitude < lat[1]) | (x$latitude > lat[2])] <- FALSE
  if (!is.null(alt)) ii[(x$altitude < alt[1]) | (x$altitude > alt[2])] <- FALSE
  if (!is.null(cntr)) {
    if (verbose) print(table(substr(x$country,1,6)))
    if (length(cntr)>1) cntr <- paste(cntr,collapse='|')
    iselect <- grep(cntr,x$country,ignore.case=TRUE,...) 
    if (sum(iselect) > 0) ii[-iselect] <- FALSE else ii[] <- FALSE
  }
  #if (!is.null(cntr)) ii[!is.element(tolower(x$country),tolower(cntr))] <- FALSE
  if (!is.null(src)) {
    if (length(src)>1) src <- paste(src,collapse='|')
    ii[-grep(src,x$source,ignore.case=TRUE,...)] <- FALSE
  }
  #if (!is.null(src)) ii[!is.element(tolower(x$source),tolower(src))] <- FALSE
  if (!is.null(nmin)) ii[!(as.numeric(x$end) - as.numeric(x$start)) < nmin - 1] <- FALSE
  if (!is.null(it)) {
    ii[(as.numeric(x$start) > min(it)) < (as.numeric(x$end) < max(it))] <- FALSE
  }
  if (verbose) print(paste('Final number of stations are',sum(ii)))
  x <- x[ii,]
  class(x) <- cls
  return(x)
}
