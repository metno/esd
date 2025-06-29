#' Match date of one object with another object or date string
#'
#' @param x input object (e.g., \code{station}, \code{field} or \code{zoo}) with a date index
#' @param it a character string with dates or an object (e.g., \code{station}, \code{field} or \code{zoo}) with a date index
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @return the part of \code{x} that matches the dates provided in \code{it} 
#' 
#' @aliases matchdate matchdate.list matchdate.default
#' @seealso subset
#' 
#' @export
matchdate <-function(x,it,verbose=FALSE) UseMethod("matchdate")

#' @exportS3Method
#' @export matchdate.list
matchdate.list <- function(x,it,verbose=FALSE) {
  if (verbose) print('matchdate.list')
  y <- lapply(x,matchdate.default,it=it,verbose=verbose)
  invisible(y)
}

#' @exportS3Method
#' @export matchdate.default
matchdate.default <- function(x,it,verbose=FALSE) {
  if(verbose) {print("matchdate.default"); print(range(index(x)))}
  ## If it is the list, then use the first element because otherwise will not find the index
  if (is.list(it)) it <- it[[1]]
  cls0 <- class(index(it))
  
  cls <- class(x)
  ## Check the index type of it and change the time scale of x to match it
  if (inherits(it,c('station','field','eof','ds'))) {
    if (verbose) print('> station/field/eof/ds:')
    if (inherits(it,'annual')) x <- annual(x)
    if (inherits(it,'month')) x <- as.monthly(x)
    if (inherits(it,'seasonal')) x <- as.4seasons(x)
    ## if (inherits(it,'day')) x <- aggregate(x,list(as.Date(index(x))),FUN='mean') ## REB 2018-11-20: what does this line actually do? It causes errors.
    if (verbose) print(index(x))
    cls0 <- class(index(it))
  } else cls0 <- class(it)

  t <- index(x)
  t0 <- t
  ## KMP 2019-09-17: the time scale is next to last in the class, not always the second element
  if (inherits(it,c('annual','month','seasonal','day'))) cls[length(cls)-1] <- class(it)[length(class(it))-1]
  
  if (inherits(it,'character')) {
    if (verbose) print('Convert years and incomplete dates to %YYYY-%MM-%DD date format')
    ## If given years but y has dates as index, convert to dates.
    nc <- nchar(it)
    # Simple fix for short date strings
    if ((nc[1]==4) & (length(it)>1)) {
      it <- paste(it,'-01-01',sep='') 
    } else if ((nc[1]==4) & (length(it)==1)) {
      it <- c(paste(it,'-01-01',sep=''),paste(it,'-12-31',sep='')) 
    } else if (nc[1]==7) {
      it <- paste(it,'-01',sep='')
    }
    it <- as.Date(it)
  }

  if (inherits(it,c('field','station','zoo'))) it <- index(it)
  if (is.logical(it)) it <- seq(1,length(it))[it]
  #print(c(t[1],it[1]));   print(c(class(t),class(it)))
  # Convert indeces all to 'Date':
  # The time index of x:
  #
  if ( (is.numeric(t)) | (is.integer(t)) | (is.numeric(t)) |
       ( (is.character(t)) & (nchar(t[1])==4) ) ) t <- as.Date(paste(t,'-01-01',sep='')) 
  if ( (is.character(t)) & (nchar(t[1])==10) )  t <- as.Date(t)
  
  # The time index to match:
  
  if ( (is.numeric(it)) | (is.integer(it)) | (is.numeric(it)) |
       ( (is.character(it)) & (nchar(it[1])==4) ) ) it <- as.Date(paste(it,'-01-01',sep='')) 
  if ( (is.character(it)) & (nchar(it[1])==10) )  it <- as.Date(it)
  
  if (verbose) { 
    print(paste('matchdate: t = [',min(t),'-',max(t),'], it= [',min(it),'-',max(it),']'))
    print(c(class(t),class(it)))
  }
  
  if (length(it)>2) {
    ii <- is.element(as.character(t),as.character(it))
    if (verbose) {
      print(paste('select',sum(ii),'dates'))
      print(t[ii])
    }
    y <- x[ii,]
    
    if (verbose) print(paste('matchdate found',sum(ii),'matching dates'))
    # KMP 2021-04-08: Changed index assignment method to zoo(y, order.by=t0[ii])) 
    # to solve a mysterious problem related to an unexpected change in the class of y. 
    #if (sum(ii)==length(index(y))) index(y) <- t0[ii] # REB 2015-01-14: to ensure same index class as it.
    if (sum(ii)==length(index(y))) y <- zoo(coredata(y), order.by=t0[ii])
    if (verbose) {print('matchdate: index(y)'); print(index(y))}
    #
  } else if (length(it)==2) {
    # Pick an interval
    if (verbose) print(paste('select an interval',it,collapse=' '))
    dates <- range(it)
    y <- window(x,start=dates[1],end=dates[2])
  } else {
    # Pick one date:
    if (verbose) print('one date')
    if (inherits(x,c('day','zoo'))) {
      ii <- is.element(t,it) 
    } else if (inherits(x,c('month','seasonal'))) {
      ii <- is.element(year(t)*100+month(t),year(it)*100+month(it))
    } else if (inherits(x,'annual')) {
      ii <- is.element(year(t),year(it))
    }
    if (sum(ii) > 0) y <- x[ii,] else {
    # Weight the two nearest in time
      if (verbose) print(paste('matchdate: Weight the two nearest in time because no overlaps:',
      	 	   	       ' sum(ii)=',sum(ii),'t = [',max(t),'-',min(t),'], it= [',
                               max(it),'-',min(it),']'))
      i1 <- t <= it; t1 <- max(t[i1])
      i2 <- t > it; t2 <- min(t[i2])
      dt <- t2 - t1
      if (verbose) print(c(max(t[i1]),min(t[i2]),sum(is.element(t,t1))))
      y1 <- (it-t1)/dt*x[is.element(t,t1),]
      y2 <- (t2 - it)/dt*x[is.element(t,t2),]
      y <- zoo(coredata(y1)+coredata(y2),order.by=it)
    }
  }
  
  y <- attrcp(x,y)
  ## REB 2025-06-12: these lines seem to cause a problem for some data
  ## What is their purpose?
  if (!is.null(err(x))) {
    if (verbose) print('match date for error')
    attr(y,'standard.error') <- try(matchdate(err(x),y))
    #str(err(y))
  }
  if (verbose) print('...')
  if (!is.null(attr(x,'n.apps'))) {
    if (verbose) print('Add appendices')
    attr(y,'n.apps') <- attr(x,'n.apps')
    attr(y,'appendix.1') <- attr(x,'appendix.1')
  }
  good <- as.logical(is.finite(index(y)))
  if (verbose) print(paste('Final check: finite values for index of y:',sum(good)))
  class(y) <- cls
  y <- subset(y,it=good)
  nt <- index(y)
  
  attr(y,'history') <- history.stamp(x)
  #print(index(y)); print(class(index(y)))

  ## REB 2021-02-15: fix to ensure similar  representation of date as original
  if ( (cls0=='Date') & (is.numeric(index(y))) ) index(y) <- as.Date(paste(index(y),'-01-01',sep='')) 

  if (inherits(y,'field')) attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],length(nt))
  if (!is.null(attr(y,'count'))) attr(y,'count') <- c(attr(y,'count')[1:2],length(nt))
  if(verbose) print('dates of x have been matched with it.')
  invisible(y) 
}
