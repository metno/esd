#' Draw a gridbox with provided colour
#'
#' Draw a gridbox of color \code{col} in a location specified by \code{x}.
#'
#' @param x location of gridbox: c(x0, x1, x2, x3, y1, y2, y3, y4, i)
#'        where x0,x1,x2,x3 and y0,y1,y2,y3 are the four corners of the box on the x- and y-axes,
#'        and i is an index specifying which element of \code{col} to use
#' @param col colors
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @export
gridbox <- function(x,col,verbose=FALSE) {
  if(verbose) print("gridbox")
  i <- round(x[9])
  brightness <- x[10]
  alpha <- x[11]
  
  ## REB 2020-01-26
  if ((brightness < 1) | (alpha < 1)) {
    ## Add shadow effect/
    col <- col2rgb(col)/255 * abs(brightness)
    #print(dim(col))
    col <- rgb(col[1,],col[2,],col[3,],alpha)
  }
  
  if(i==0) {
    border <- rgb(0,0,0,0.1)
  } else {
    border <- col[i]
  }
 
  polygon(c(x[1:4],x[1]),c(x[5:8],x[5]),
          col=col[i],border=border)
}
