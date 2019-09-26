#' Extention of points
#'
#' @param x input vector
#' @param y input vector of same length as x
#' @param col vector of colors
#' @param cex.max maximum size of balls
#' @param n length of color scale
#'
#' @export
balls <- function(x,y=NULL,col=NULL,cex.max=2,n=20) {
  for (i in 1:n) {
    if ((is.null(col)) | length(col)==1) {
      cols <- rgb(i/n,i/n,i/n) 
    } else if (is.vector(col)) {
      cols <- col/i
      cols[cols==0] <- i/n
      cols <- rgb(cols[1],cols[2],cols[3])
    } else 
      if (is.vector(character)) cols <- col[i]    
      
      points(x,y,cex=seq(cex.max,0.1,length=n)[i],
             col=cols)
  }
}
