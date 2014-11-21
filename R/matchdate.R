matchdate <- function(x,it) {
  #print('matchdate')
  t <- index(x)

  if (inherits(it,'character')) {
    nc <- nchar(it)
    # Simple fix for short date strings
    if ((nc[1]==4) & (length(it)>1)) it <- paste(it,'-01-01',sep='') else
    if ((nc[1]==4) & (length(it)==1))
      it <- c(paste(it,'-01-01',sep=''),paste(it,'-12-31',sep='')) else
    if (nc[1]==7) it <- paste(it,'-01',sep='')
    it <- as.Date(it)
  } 
  if (inherits(it,c('field','station','zoo'))) it <- index(it)
  #print(c(t[1],it[1]));   print(c(class(t),class(it)))
  
  # Convert indeces all to 'Date':
  # The time index of x:
  #browser()
  if ( (is.numeric(t)) | (is.integer(t)) |
     ( (is.character(t)) & (nchar(t[1])==4) ) ) t <- as.Date(paste(t,'-01-01',sep='')) 
  if ( (is.character(t)) & (nchar(t[1])==10) )  t <- as.Date(t)

  # The time index to match:
  if ( (is.numeric(it)) | (is.integer(it)) |
     ( (is.character(it)) & (nchar(it[1])==4) ) ) it <- as.Date(paste(it,'-01-01',sep='')) 
  if ( (is.character(it)) & (nchar(it[1])==10) )  it <- as.Date(it)
  
  #print(c(t[1],it[1]))
  
  if (length(it)>2) {
    y <- x[is.element(t,it),]
  } else if (length(it)==2) {
  # Pick an interval
    #print('an interval')
    dates <- range(it)
    y <- window(x,start=dates[1],end=dates[2])
  } else {
  # Pick one date:
    #print('one date')
    if (inherits(x,c('day','zoo'))) ii <- is.element(t,it) else
    if (inherits(x,c('month','seasonal')))
      ii <- is.element(year(t)*100+month(t),year(it)*100+month(it)) else
    if (inherits(x,'annual')) ii <- is.element(year(t),year(it)*100)
    #print(class(x)); print(sum(ii))
    if (sum(ii) > 0) y <- x[ii,] else {
  # Weight the two nearest in time      
       i1 <- t <= it; t1 <- max(t[i1])
       i2 <- t > it; t2 <- min(t[i2])
       dt <- t2 - t1
       y1 <- (it-t1)/dt*x[is.element(t,t1),]
       y2 <- (t2 - it)/dt*x[is.element(t,t2),]
       y <- zoo(coredata(y1)+coredata(y2),order.by=it)
     }
  }
  y <- attrcp(x,y)
  if (!is.null(attr(y,'standard.error'))) {
    attr(y,'standard.error') <- matchdate(attr(y,'standard.error'),y)
  }
  nt <- index(y)
  attr(y,'history') <- history.stamp(x)
  #print(index(y)); print(class(index(y)))
  class(y) <- class(x)
  if (inherits(y,'field')) attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],length(nt))
  if (!is-null(attr(y,'count')) attr(y,'count') <- c(attr(y,'count')[1:2],length(nt))
  invisible(y) 
}
