sametimescale <- function(y,X,FUN='mean',verbose=FALSE) {
  ### Function to ensure that station y has the same time scales as X
  
  if (verbose) print('sametimescale')
  tsx <- class(X)[length(class(X))-1]
  tsy <- class(y)[length(class(y))-1]
  if (verbose) print(c(tsx,tsy))
  if (tsx==tsy) return(y)
  
  if (verbose) print('Need to aggregate')
  if (tsx=="day") {
    agrscly <- as.Date(index(y)) 
  } else if (tsx=="month") {
    agrscly <- as.yearmon(index(y)) 
  } else if (tsx=="annual") {
    agrscly <- year(y) 
  } else if (tsx=="year") {
    agrscly <- year(y)
  }
  if (tsx=="season") {
    y <- as.4seasons(y, FUN=FUN, dateindex=TRUE)
    agrscly <- index(y) # KMP 2018-12-26: agrscly needed if verbose = TRUE
    ## y <- aggregate(y, agrscly, match.fun(FUN)) else ## REB, 2018-11-20 - match.fun caused problems
    ## y <- as.4seasons(y, FUN=match.fun(FUN),dateindex=TRUE) ## REB, 2018-11-20 - match.fun caused problems
  } else {
    y <- aggregate(y, agrscly, FUN) 
  }
  
  if (verbose) {
    str(agrscly)
    print(FUN)
  }
  ##
  if(verbose) print(c(class(index(y)),class(index(X))))
  if ( (class(index(y))=='Date') & (class(index(X))=='numeric') & 
       ((tsx=='year') |  (tsx=='annual')) ) {
    index(y) <- year(index(y))
  }
  if ( (class(index(y))=='numeric') & (class(index(X))=='Date') & 
       ((tsx=='year') |  (tsx=='annual')) ) {
    index(y) <- as.Date(paste(index(y),'01-01',sep='-'))
  }
  invisible(y)
}
