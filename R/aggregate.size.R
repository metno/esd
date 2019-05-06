## Function that calculates the size of events
## Matrix lon & lat indices: (i, j).
## @RasmusBenestad, 2018-10-29

aggregate.size <- function(x, ...) UseMethod("aggregate.size")

aggregate.size.matrix <- function(x,x0,plot=FALSE,verbose=FALSE,a=6.378e06,...) {

    ## Select all grid boxes with values exceeding x0
    if (verbose) print('aggregate.size.matrix')
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
            if (sum(is.na(sel)) >0) {print("Something is wrong!"); browser()}
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


aggregate.size.field <- function(x,x0,plot=FALSE,verbose=FALSE,...) {
  if (verbose) print('aggregate.size.field')  
  nt <- length(index(x))
  d <- attr(x,'dimensions')
  sizestats <- list()
  for (it in 1:nt) {
    z <- matrix(x[it,],d[1],d[2])
    attr(z,'longitude') <- lon(x); attr(z,'latitude') <- lat(x)              
    sizestats[[it]] <- aggregate.size(z,x0=x0)$statistic
    if (verbose) print(sizestats[[it]])  
  }
  names(sizestats) <- index(x)
  ne <- max(unlist(lapply(sizestats,length)))
  size <-  unlist(lapply(sizestats,function(x) {y <- rep(0,ne); y[1:length(x)] <- x[]; return(y)}))
  dim(size) <- c(ne,length(sizestats))
  size <- zoo(t(size),order.by=index(x))
  names(size) <- paste('event',1:ne)
  attr(size,'unit') <- 'm^2' 
  if (plot) plot(size)
  return(size)
}


test.events <- function(n=62, m=78, verbose=TRUE) {
    x <- matrix(rep(sin(seq(0,4*pi,length=n)),m)*
                sort(rep(cos(seq(0,pi,length=m)),n)),n,m)
    par(mfcol=c(3,1))
    image(x)
    test.results <- aggregate.size.matrix(x,x0=0,plot=TRUE,verbose=verbose)
    
}
