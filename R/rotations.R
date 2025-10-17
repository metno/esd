# Helper function for X-axis rotation (Pure R_X)
#' @export
rotX <- function(angle_rad) {
  cosA <- cos(angle_rad)
  sinA <- sin(angle_rad)
  rbind(
    c(1, 0, 0),
    c(0, cosA, -sinA),
    c(0, sinA, cosA)
  )
}

# Helper function for Y-axis rotation (Pure R_Y)
#' @export
rotY <- function(angle_rad) {
  cosY <- cos(angle_rad)
  sinY <- sin(angle_rad)
  rbind(
    c(cosY, 0, sinY),
    c(0, 1, 0),
    c(-sinY, 0, cosY)
  )
}

# Helper function for Z-axis rotation (Pure R_Z)
#' @export
rotZ <- function(angle_rad) {
  cosZ <- cos(angle_rad)
  sinZ <- sin(angle_rad)
  rbind(
    c(cosZ, -sinZ, 0),
    c(sinZ, cosZ, 0),
    c(0, 0, 1)
  )
}