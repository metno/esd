## Function that calculates the size of events
## Matrix lon & lat indices: (i, j).
## @RasmusBenestad, 2018-10-29


aggregate.size.matrix <- function(x,x0,plot=FALSE,verbose=FALSE,...) {

    recursive.search <- function(z,y,I,J,eno) {
      nbr <- (I[-1] - I[1])^2 + (J[-1] - J[1])^2
      sel <- (nbr==1)
      while (sum(sel,na.rm=TRUE)>0) {
          z[1] <- eno; I[1] <- NA; J[1] <- NA
          if (length(sel)> 0) {
              psel <- recursive.search(z[sel],y[sel],I[sel],J[sel],eno)
              
          } else return(NULL)
      }   
      list(z,y,I,J,eno)
    }
    
    if (verbose) print('aggregate.size.matrix')
    ## Copy data that can be masked
    mask <- (x > x0); y <- x[mask];  ## Space mask
    z <- y*0
    if (verbose) print(paste('Detected',sum(mask),'grid points'))
    ## Determine the dimensions of the matrix
    d <- dim(x)
    ## Generate same size matrix with i-indexes
    I <- matrix(rep(1:d[1],d[2]),d[1],d[2])[mask]
    J <- matrix(sort(rep(1:d[2],d[1])),d[1],d[2])[mask]
    if (plot) plot(I,J)
    
    ##  Loop: distance (i, j) =1. Store and mask points that are sorted.
    ## Zero, one, two, three or four cases for each point.
    eno <- 1
    recursive.search(z,y,I,J,eno)
    
}


test.events <- function(n=100, m=80) {
    x <- matrix(rep(sin(seq(0,4*pi,length=n)),m)*
                sort(rep(cos(seq(0,pi,length=m)),n)),n,m)
    aggregate.size.matrix(x,x0=0,plot=TRUE,verbose=TRUE)
}
