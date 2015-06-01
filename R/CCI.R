# K Parding, 29.05.2015

library(esd)
slp <- slp.ERAINT()

CCI <- function(Z,m=14,nsim=10,it=NULL,is=NULL,
                cyclones=TRUE,accuracy=NULL,verbose=FALSE) {

  stopifnot(inherits(Z,'field'))
  Z <- subset(Z,it=it,is=is)
  resx <- dX(Z,m=m,accuracy=accuracy,verbose=verbose)
  resy <- dY(Z,m=m,accuracy=accuracy,verbose=verbose)

  dslpdx <- resx$dZ
  dslpdx2 <- resx$dZ2
  dslpdy <- resy$dZ
  dslpdy2 <- resy$dZ2
  wind <- dslpdx*0
}
