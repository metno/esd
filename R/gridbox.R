gridbox <- function(x,col,verbose=FALSE) {
  if(verbose) print("gridbox")
  i <- round(x[9])
  if(i==0) {
    border <- rgb(0,0,0,0.1)
  } else {
    border <- col[i]
  }
  polygon(c(x[1:4],x[1]),c(x[5:8],x[5]),
          col=col[i],border=border)
}
