validate <- function(x, ...)  UseMethod("validate")

validate.default <- function(x, ...) {
}

validate.eof <- function(x, ...) {
}

validate.pca <- function(x, ...) {
}

validate.eof.field <- function(x, ...) {
}


validate.eof.comb <- function(x, ...) {
  plot(x)
  zz <- attr(x,'appendix.1')
  
  dev.new()
  plot(attr(x,'clim'),type="l",main="Mean values")
  lines(attr(zz,'clim'),col="red")
}

 
validate.ds <- function(x, ...) {
  
}

validate.cca <- function(x, ...) {
}


test.EOF <- function() {
  
}
