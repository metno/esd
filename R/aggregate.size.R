## Function that calculates the size of events
## Matrix lon & lat indices: (i, j).
## @RasmusBenestad, 2018-10-29

aggregate.size <- function(x, ...) UseMethod("aggregate.size")

aggregate.size.matrix <- function(x,x0,plot=FALSE,verbose=FALSE,a=6.378e06,...) {
  
  ## Select all grid boxes with values exceeding x0
  if (verbose) print('aggregate.size.matrix')
  ## Copy data that can be masked
  mask <- (x > x0); mask[is.na(mask)] <- FALSE
  if (sum(mask,na.rm=TRUE)==0) return(list(events=x*NA,number=0,statistic=NULL))
  y <- x[mask];  ## Space mask
  z <- y*NA; ns <- sum(mask,na.rm=TRUE)
  if (verbose) print(paste('Detected',ns,'grid points from total of',length(x)))
  ## Determine the dimensions of the matrix
  d <- dim(x)
  ## Generate same size matrix with i-indexes
  I <- matrix(rep(1:d[1],d[2]),d[1],d[2])[mask]
  J <- matrix(sort(rep(1:d[2],d[1])),d[1],d[2])[mask]
  if (verbose) {print(d); print(summary(I)); print(summary(J)); print(table(mask))}
  if (plot) {
    plot(I,J,pch=19,col='grey70',cex=0.5)
    points(I[1],J[1],pch='1')
  }
  ##  Loop: distance (i, j) =1. Store and mask points that are sorted.
  ## Zero, one, two, three or four cases for each point.
  eno <- 1; ij <- 1; z[ij] <- 1 ## First selected gridbox always belongs to event-number 1
  if (verbose) print('Loop to estimate area for each event')
  while (sum(is.na(z))>0) {
    ## If the event-number has not been assigned: 
    ## see if it is adjacent to a gridpoint with an event-number
    
    ## Distance between the first and all the other grid boxes:
    nbr <- sqrt( (I - I[ij])^2 + (J - J[ij])^2 )
    if (verbose) print(min(nbr[-ij]))
    ## Select the gridboxes attached to the fist grid box: nbr==1
    sel <- (nbr==1)
    
    if (is.na(z[ij])) {
      if (sum(is.na(sel)) >0) {print("No adjacent boxes - Something is wrong!"); browser()}
      if (verbose) print(paste('N adjacent points=',sum(sel,na.rm=TRUE),
                               'Cluster no=',eno))
      if (sum(sel)>0) {
        getnumber <- z[sel]
        assigned <- getnumber[is.finite(getnumber)]
        if (verbose) print(assigned)
        if (length(assigned)>0) {
          if (max(assigned,na.rm=TRUE) - min(assigned,na.rm=TRUE)) {
            print("Numbers not the same - something is wrong!"); print(assigned);
          }
          z[ij] <- min(assigned); z[sel] <- min(assigned) 
          eno <- min(assigned)
        } else {
          if (verbose) print(paste('No valid assigned number',sum(sel)))
          ## None of the adjacent gridboxes have an event-number, the gridbox gets a new event number
          eno <- eno + 1
          z[ij] <- eno; z[sel] <- eno
        }
        ij <- min((1:ns)[sel])
      } else {
        if (verbose) print('New cluster') 
        ## If there are no adjacent events, then the gridbox gets a new event number
        ij <- min((1:ns)[is.na(z)])
        eno <- eno + 1
        z[ij] <- eno
      }
      sel <- (z == eno)
      if (plot) points(I[sel],J[sel],pch=as.character(eno),
                       cex=0.5,col=rainbow(21)[(eno-1) %/% 5 + 1])
    } else {
      if (verbose) cat('.')
      if (ij != min((1:ns)[sel])) ij <- min((1:ns)[sel]) else ij <- ij +1
    }
    #rm("nbr"); gc(reset=TRUE)
  }
  ## Synthesise the results in the  
  Z <- x*NA; Z[mask] <- z
  if (plot) image(Z)
  if ( (!is.null(attr(z,'longitude'))) & (!is.null(attr(z,'latitude'))) ) {
    ## Estimate the area of the gridboxes
    W <- Z; for (j in 1:d[2]) W[,j] <- a*cos(pi*lat(z)[j]/180)
    stats <- rep(NA,eno)
    for (i in 1:eno) {stats[i] <- sum(W[Z==eno],na.rm=TRUE)}  
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
