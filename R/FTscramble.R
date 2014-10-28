# A phase scrambling of data - to model noise.
# Rasmus Benestad

FTscramble <- function(x) {
  X <- fft(x); n <- length(x)
  Z <- Mod(X)
  ReX <- Re(X)
  ImX <- Im(X) 
  ReY <- runif(n,min=-Z,max=Z)
  ImY <- sqrt(Z^2 - ReY^2)
  ReY[1] <- ReX[1]; ImY[1] <- ImX[1]
#  ReY[n] <- ReX[n]; ImY[n] <- ImX[n]
  Y <- complex(real=ReY, imaginary=ImY)
  y <- fft(Y,inverse=TRUE)
  invisible(Re(y)/n)
}

