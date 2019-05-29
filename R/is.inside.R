#' A test to see if a point is inside a polygon
#'
#' \code{is.inside} checks whether a point or a set of points is inside a polygon, e.g., borders of a country.
#' 
#' @aliases is.inside test.is.inside
#'
#' @param x an esd-object or a list/data.frame with the elements \code{x$x} and
#' \code{x$y} containing the coordinates.
#' @param y A polygon in the shape of a list/data.frame with the elements
#' \code{x$x} and \code{x$y} containing the coordinates.
#' @param verbose \code{TRUE} prints out diagnostics for the code.
#' @param plot \code{TRUE} provides a graphical disgnostic.
#' @param N Number of tests with random coordinates
#' 
#' @return a boolean; \code{TRUE} if the point(s) \code{x} is/are inside the polygon \code{y}.
#' 
#' @examples
#' 
#' \dontrun{ 
#' library(readINAMdata)
#' data("Moz")
#' 
#' osmoz <- !is.inside(data.frame(x=lons,y=lats),Moz)
#' plot(lons,lats)
#' points(lons[osmoz],lats[osmoz],pch=19,col='red')
#' points(lons[!osmoz],lats[!osmoz],pch=19,col='green')
#' points(lons[is.na(osmoz)],lats[is.na(osmoz)],pch=19,col='black')
#' lines(Moz,type='b')
#' }
#' 
#' @export is.inside
is.inside <- function(x,y,verbose=FALSE,plot=FALSE) {
  if (verbose) print('is.inside')
  if (plot) plot(y,type='l')
  ## x is the location of the esd-object that we want to see is within the polygon y.
  ## If y is a data.frame or list, check if it contains the elements 'x' and 'y'
  ## Assign attributes so we can use lon() and lat()
  if ((is.list(x)) | is.data.frame(x)) {
    els <- names(x); nn <- length(els)
    ixy <- (1:nn)[is.element(substr(tolower(els),1,3),c('x','lon'))][1]
    iyy <- (1:nn)[is.element(substr(tolower(els),1,3),c('y','lat'))][1]
    if (verbose) print(c(els[ixy],els[iyy]))
    X <- data.frame(x=x[[els[ixy]]],y=x[[els[ixy]]])
  } else if ( (!is.null(lon(x))) & (!is.null(lat(x))) ) X <- data.frame(x=lon(x),y=lat(x))
  if ((is.list(y)) | is.data.frame(y)) {
    els <- names(y); nn <- length(els)
    ixy <- (1:nn)[is.element(substr(tolower(els),1,3),c('x','lon'))][1]
    iyy <- (1:nn)[is.element(substr(tolower(els),1,3),c('y','lat'))][1]
    if (verbose) print(c(els[ixy],els[iyy]))
    Y <- data.frame(x=y[[els[ixy]]],y=y[[els[ixy]]])
  } else if ( (!is.null(lon(y))) & (!is.null(lat(y))) ) Y <- data.frame(x=lon(y),y=lat(y))
  
  d <- dim(x)
  if (verbose) print(d)
  
  inside <- rep(FALSE,d[1])
  for (i in 1:d[1]) {
    if (verbose) print(i)
    inside[i] <- is.inside.1(x[i,],y,verbose=verbose,plot=plot)
  }
  if (verbose) print(table(inside))
  if (plot) { 
    col <- rep('lightblue',dim(x)[1])
    col[inside] <- 'darkgreen'
    points(x,pch=19,cex=1.2,col=col)
    #text(x$x,x$y,1:length(x$x),sub=1,cex=0.6,font=2)
  }
  invisible(inside)
}

is.inside.1 <- function(x,y,verbose=FALSE,plot=FALSE) {
  if (verbose) print('is.inside.1')
  ## Check if the point shares the same coordinates as the polygon -if so, move it slightly
  if (sum(x$x == y$x)>0) {
    x$x <- x$x + 0.01 * min(diff(y$x),na.rm=TRUE)
  }
  if (sum(x$y == y$y)>0) {
    x$y <- x$y + 0.01 * min(diff(y$y),na.rm=TRUE)
  }
  
  ## First do the quick simple test: is the coordinate of x within the range of the coordinates of y?
  if ( (x$x < min(y$x,na.rm=TRUE)) |
       (x$x > max(y$x,na.rm=TRUE)) |
       (x$y < min(y$y,na.rm=TRUE)) |
       (x$y > max(y$y,na.rm=TRUE)) ) {
    if (verbose) {print('outside range of the shape'); print(x)}
    return(FALSE)
  }
  
  # if (plot) {
  #   lines(rep(x$x,2),range(y$y),lty=2,col=rgb(1,0.7,0.5))
  #   lines(range(y$x),rep(x$y,2),lty=2,col=rgb(1,0.7,0.5))
  # }
  
  ## Compare the x and y coordinates of x and y
  dx <- (y$x - x$x); nx <- length(dx)
  dy <- (y$y - x$y); ny <- length(dy)
  ## Check if the location has common x/y-points:
  x0 <- (dx == 0); y0 <- (dy == 0)
  ## Check where the coordinates falls between the adjacent points on the border
  x1 <- c(dx[2:nx]*dx[1:(nx-1)],dx[1]*dx[nx])
  y1 <- c(dy[2:ny]*dy[1:(ny-1)],dy[1]*dy[ny]) 
  ## Count the number of intesects between the two axes and the border
  i1 <- (x0 | (x1 < 0)) & (y$y < x$y)
  i2 <- (x0 | (x1 < 0)) & (y$y >= x$y)
  j1 <- (y0 | (y1 < 0)) & (y$x < x$x)
  j2 <- (y0 | (y1 < 0)) & (y$x >= x$x)
  
  ## If all the numbers of intersect are odd numbers, then x is inside y
  #if (plot) points(y[i1 | i2 | j1 | j2,],pch='x',cex=0.75,col='red',lwd=2)
  nx1 <- sum(i1,na.rm=TRUE); nx2 <- sum(i2,na.rm=TRUE)
  ny1 <- sum(j1,na.rm=TRUE); ny2 <- sum(j2,na.rm=TRUE)
  ## Extra check: there should be an even number of line crossings in total for each axis
 
  inside <- (nx1%%2==1) & (nx2%%2==1) & (ny1%%2==1) & (ny2%%2==1)
  
  if (verbose) {
    print(c(nx1,nx2,ny1,ny2))
    print(x)
    print(inside)
  }
  
  if  ( ((nx1 + nx2)%%2 !=0) | ((ny1 + ny2)%%2 !=0) ) {
    if (plot) points(x$x,x$y,col='red')
    print(c(nx1,nx2,ny1,ny2))
    print(x)
    print(inside)
    inside <- NA
  }
  return(inside)
}

## Estimate the distance between the adjacent points in the data frame 
pdist <- function(x,y) {
  d <- (x[,1]-y[,1])^2 + (x[,2]-y[,2])^2 
  return(d)
}

#' @export
test.is.inside <- function(N=5000,verbose=FALSE,plot=TRUE) {
  require(esd)
  data(geoborders)
  ii <- ((geoborders$x > 112) & (geoborders$x < 155) &
           (geoborders$y > -40) & (geoborders$y < -11.5))
  Aus <- geoborders[ii,]
  ok <- is.finite(Aus$x) & is.finite(Aus$y)
  Aus <- Aus[ok,]
  aus <- Aus
  print(Aus[1,]); print(dim(Aus))
  Aus[-1,] <- NA; aus[1,] <- NA
  n <- length(Aus[,1])
  for (i in 2:n) {
    d <- pdist(Aus[i-1,],aus)
    if (sum(is.finite(d))>0) {
      ii <- (1:n)[(d == min(d,na.rm=TRUE))]
      ii <- ii[is.finite(ii)][1]
      if (length(ii)>0) {
        Aus[i,] <- aus[ii,]
        aus[ii,] <- NA
      } else cat(paste(' NULL ii',i))
    } else cat(paste(' zero-d',i))
  }
  D <- c(pdist(Aus[2:n,],Aus[1:(n-1),]),0)
  il <- D < 0.001
  Aus <- Aus[il,]; D <- D[il]
  Aus <- rbind(Aus,Aus[1,])
  n <- length(Aus[,1])
  if (verbose) {print(summary(D)); print(n); print(dim(Aus))}
  if (verbose) print(summary(Aus))
  Aus <- Aus[1:671,]
  lon <- runif(N)*50 + 100
  lat <- runif(N)*40 - 40
  x <- data.frame(x=lon,y=lat)
  inside <- is.inside(x,y=Aus,verbose=verbose,plot=plot)
}
