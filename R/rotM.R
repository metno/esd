#' @export
rotM <- function(x=0,y=0,z=0) {
  X <- -pi*x/180; Y <- -pi*y/180; Z <- -pi*z/180
  cosX <- cos(X); sinX <- sin(X); cosY <- cos(Y)
  sinY <- sin(Y); cosZ <- cos(Z); sinZ <- sin(Z)
  rotM <- rbind( c( cosY*cosZ, cosY*sinZ,-sinY*cosZ ),
                 c(-cosX*sinZ, cosX*cosZ,-sinX*cosZ ),
                 c( sinY*cosX, sinX*cosY, cosX*cosY ) )
  return(rotM)
}
