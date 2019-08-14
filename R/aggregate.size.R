#' aggregate
#' 
#' The aggregation functions are based on the S3 method for \code{zoo} objects,
#' but takes care of extra house keeping, such as attributes with meta data.
#'
#' \code{aggregate.size} is similar to \code{aggregate.area}, but returns the size statistics (square
#' meters) for individual events (defined as gridboxes touching each other).
#'
#' The function \code{aggregateSize} is exactly the same as \code{aggregate.size}.
#'
#' @aliases aggregateSize aggregateSize.matrix aggregateSize.field
#' @seealso aggregate.area aggregate
#' 
#' @param x A \code{\link{station}} object
#' @param x0 threshold defining an event
#' @param plot a boolean; if TRUE display results as a plot
#' @param a radius of earth (unit: km)
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#'
#' @return The call returns a station object
#'
#' @author R.E. Benestad
#' @keywords utilities
#'
#' @export aggregate.size
aggregate.size <- function(x, ...) {
  aggregateSize(x, ...)
}

#' @export aggregateSize
aggregateSize <- function(x, ...) UseMethod("aggregateSize")

# Aggregate size of events - S3 method for matrix
#' @export aggregateSize.matrix
aggregateSize.matrix <- function(x,x0,plot=FALSE,verbose=FALSE,a=6378,...) {#,a=6.378e06,...) {

    ## Select all grid boxes with values exceeding x0
    if (verbose) print('aggregateSize.matrix')
    ## Copy data that can be masked
    mask <- (x > x0)
    if (sum(mask)==0) return(list(events=x*NA,number=0,statistic=NULL))
    y <- x[mask];  ## Space mask
    z <- y*NA; ns <- sum(mask)
    if (verbose) print(paste('Detected',ns,'grid points from total of',length(x)))
    ## Determine the dimensions of the matrix
    d <- dim(x)
    ## Generate same size matrix with i-indexes
    I <- matrix(rep(1:d[1],d[2]),d[1],d[2])[mask]
    J <- matrix(sort(rep(1:d[2],d[1])),d[1],d[2])[mask]
    if (plot) {
        plot(I,J,pch=19,col='grey70')
        points(I[1],J[1],pch='1')
    }
    ##  Loop: distance (i, j) =1. Store and mask points that are sorted.
    ## Zero, one, two, three or four cases for each point.
    eno <- 1; z[1] <- 1 ## First selected gridbox always belongs to event-number 1
    for (ij in 1:ns) {
        
        ## If the event-number has not been assigned: see if it is adjacent to a gridpoint with an event-number
        if (is.na(z[ij])) {
            ## Distance between the first and all the other grid boxes:
            nbr <- c(0,sqrt( (I[-ij] - I[ij])^2 + (J[-ij] - J[ij])^2 ))
            ## The tested gridpoint has distance to itsel nbr=0
            if (ij > 1) {nbr[1:(ij-1)] <- nbr[2:ij]; nbr[ij] <- 0}
            ## Select the gridboxes attached to the fist grid box: nbr==1
            sel <- (nbr==1)
            if (sum(is.na(sel)) >0) {
	      print("Something is wrong!")
	      browser()
	    }
            if (sum(sel)>0) {
                if (verbose) print(paste('more adjacent points',sum(sel,na.rm=TRUE),eno))
                getnumber <- z[sel]
                assigned <- getnumber[is.finite(getnumber)]
                if (length(assigned)>0) {
                    z[ij] <- assigned[1]; z[sel] <- assigned[1] 
                } else {
                    ## None of the adjacent gridboxes have an event-number, the gridbox gets a new event number
                    eno <- eno + 1
                    z[ij] <- eno
                    z[sel] <- eno
                }
            } else {
                ## If there are no adjacent events, then the gridbox gets a new event number
                eno <- eno + 1
                z[ij] <- eno
            }
            if (plot) points(I[ij],J[ij],pch=as.character(z[ij]))
        }
    }
    ## Synthesise the results in the  
    Z <- x*NA; Z[mask] <- z
    if (plot) image(Z)
    if ( (!is.null(attr(z,'longitude'))) & (!is.null(attr(z,'latitude'))) ) {
      ## Estimate the area of the gridboxes
      W <- Z; for (j in 1:d[2]) W[,j] <- a*cos(pi*lat(z)[j]/180)
      stats <- rep(NA,eno)
      for (i in 1:eno) {stats[i] <- sum(W[Z==eno])}  
    } else stats <- table(z)
    ## The biggest event is always the first
    stats <- sort(stats)
    if (verbose) print(stats)
    return(list(events=Z,number=eno,statistic=stats))
}

# Aggregate size of events - S3 method for field
#' @export aggregateSize.field
aggregateSize.field <- function(x,x0,plot=FALSE,verbose=FALSE,...) {
  if (verbose) print('aggregateSize.field')  
  nt <- length(index(x))
  d <- attr(x,'dimensions')
  sizestats <- list()
  for (it in 1:nt) {
    z <- matrix(x[it,],d[1],d[2])
    attr(z,'longitude') <- lon(x)
    attr(z,'latitude') <- lat(x)
    stat.it <- aggregateSize(z,x0=x0)$statistic
    if(is.null(stat.it)) {
      sizestats[[it]] <- 0
    } else {
      sizestats[[it]] <- as.numeric(stat.it)
    }
  }
  names(sizestats) <- index(x)[1:length(sizestats)]
  ne <- max(unlist(lapply(sizestats,length)))
  size <-  unlist(lapply(sizestats,
       function(x) {
         y <- rep(0,ne)
	 if(length(x)>0) y[1:length(x)] <- x
    	 return(y)
       } ))
  dim(size) <- c(ne,length(sizestats))
  size <- zoo(t(size),order.by=index(x))
  names(size) <- paste('event',1:ne)
  attr(size,'unit') <- 'm^2' 
  if (plot) plot(size)
  return(size)
}

# Test function - do not export
test.aggregate.size <- function(n=62, m=78, verbose=TRUE) {
    x <- matrix(rep(sin(seq(0,4*pi,length=n)),m)*
                sort(rep(cos(seq(0,pi,length=m)),n)),n,m)
    par(mfcol=c(3,1))
    image(x)
    test.results <- aggregateSize(x,x0=0,plot=TRUE,verbose=verbose)
    
}
