## Function that calculates the size of events
## Matrix lon & lat indices: (i, j).
## @RasmusBenestad, 2018-10-29


aggregate.size.matrix <- function(x,x0,plot=FALSE,verbose=FALSE,...) {

    
    recursive.search <- function(z,y,I,J,eno) {
        if (verbose) print(paste('recursive.search:',sum(is.finite(y)),eno))
        if (sum(is.finite(y))==0) browser()
        ## Only keep the points where y has valid data
        valgp <- !is.na(y); y <- y[valgp]; I <- I[valgp]; J <- J[valgp]
        if (sum(valgp,na.rm=TRUE)==0) {
            if (verbose) print('finished')
            return(list(z,y,I,J,eno))
        }

  
        ## Distance between the first and all the other grid boxes:
        nbr <- c(0,sqrt( (I[-1] - I[1])^2 + (J[-1] - J[1])^2 ))
        
        ## Sort the gridpoints according to distance
        srtgp <- order(nbr); nbr <- nbr[srtgp]; y <- y[srtgp]; I <- I[srtgp]; J <- J[srtgp]
        ## Select the girboxes which lie next to the fist grid box:
        sel <- (nbr==1)
        while (sum(sel,na.rm=TRUE)>0) {
            if (verbose) print(paste('more adjacen points',(sum(sel,na.rm=TRUE))))
            
            if (plot) points(I[1],J[1],pch=19,col=rgb(1,0,0,0.2))
            ## If a grid lying next to grid one is identified, then mask grid one, keep a record,
            ## (z[i1] <- eno) and continue search
            z[I[1],J[1]] <- eno; y[1] <- NA; I[1] <- NA; J[1] <- NA
            is <- (1:length(sel))[sel]  
            if ( (length(sel)>=1) & (sum(is.finite(y))==0) ) {
                ## If one adjacent gridbox is detected
                psel <- recursive.search(z,y,I,J,eno)
            }
            ## If more than one neighbouring gridpoints: 2 or more
            if ( (length(sel) >= 2) & (sum(is.finite(y))==0) ) {
                if (verbose) print(paste('>=2',is[1]))
                y[is[1]] <- NA; I[is[1]] <- NA; J[is[1]] <- NA
                psel <- recursive.search(z,y,I,J,eno)
            }
            ## If more than one neighbouring gridpoints: 3 or more
            if ( (length(sel) >= 3)  & (sum(is.finite(y))==0) ) {
                if (verbose) print(paste('>=3',is[1]))
                y[is[2]] <- NA; I[is[2]] <- NA; J[is[2]] <- NA
                psel <- recursive.search(z,y,I,J,eno)
            }
            ## If more than one neighbouring gridpoints: 4
            if ( (length(sel) == 4)  & (sum(is.finite(y))==0) ) {
                if (verbose) print(paste('>=4',is[1]))
                y[is[3]] <- NA; I[is[3]] <- NA; J[is[3]] <- NA
                psel <- recursive.search(z,y,I,J,eno)
            }
            ## Distance between the first and all the other grid boxes:
            nbr <- c(0,sqrt( (I[-1] - I[1])^2 + (J[-1] - J[1])^2 ))
            
            ## Sort the gridpoints according to distance
            srtgp <- order(nbr); nbr <- nbr[srtgp]; y <- y[srtgp]; I <- I[srtgp]; J <- J[srtgp]
            ## Select the girboxes which lie next to the fist grid box:
            sel <- (nbr==1)
            
        }
        ## If there are no more grid boxes next to the ones identified, then return the
        ## data and increase the event number (evo):
        #browser()
        if (verbose) print('...')
        return(list(z,y,I,J,eno))
    }

    if (verbose) print('aggregate.size.matrix')
    ## Copy data that can be masked
    mask <- (x > x0); y <- x[mask];  ## Space mask
    z <- x*NA
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

    ## Loop over unmasked points which have been selected
    while (sum(is.finite(y))>0) {
        Z <- recursive.search(z,y,I,J,eno)
        z <- Z$z; y <- Z$y; I <- Z$I; J <- Z$J; eno <- Z$eno + 1
        browser()
    }
    
    return(Z)
}

test.events <- function(n=40, m=30, verbose=TRUE) {
    x <- matrix(rep(sin(seq(0,4*pi,length=n)),m)*
                sort(rep(cos(seq(0,pi,length=m)),n)),n,m)
    par(mfcol=c(2,1))
    image(x)
    test.results <- aggregate.size.matrix(x,x0=0,plot=TRUE,verbose=verbose)
    
}
