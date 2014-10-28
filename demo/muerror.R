# Rasmus Benestad (Athens, 26.09.2014)
# A script meant to test the intrinsic sampling uncertainty associated
# with the mean for a variable that is exponentially distributed through
# a set of Monte-Carlo simulations.
# The implications are that the value for mu estimated from observations
# are associated with an irreducable errorbar that is a function of the
# wet-day frequency.
N <- 10000
x <- 6+3*rnorm(N)
dim(x) <- c(N,1)

emu <- function(mu,n=100) {
  y <- mean(rexp(n,rate=1/mu))
  y
}

y30 <- apply(x,1,emu,n=30)    # 30 rainy days per year
y50 <- apply(x,1,emu,n=50)    # 50 rainy days per year
y100 <- apply(x,1,emu,n=100)  # 100 rainy days per year (f ~ 0.27)

par(bty='n')
plot(x,y30,pch=19,col=rgb(1,0,0,0.2),
     main='Sampling uncertainty for exponential distribution',
     xlab=expression(paste('original ',mu)),
     ylab=expression(paste('sampled ',mu)))
grid()

points(x,y50,pch=19,col=rgb(0,0,1,0.2))
points(x,y100,pch=19,col=rgb(0,1,0,0.2))

lines(range(x),range(x))
# Confidence: make use of the property that var = mean^2
# standard error: +- 2 sigma/sqrt(N-1)
lines(range(x),range(x)*(1+2/sqrt(29)),col='darkred')
lines(range(x),range(x)*(1-2/sqrt(29)),col='darkred')
lines(range(x),range(x)*(1+2/sqrt(49)),col='darkblue')
lines(range(x),range(x)*(1-2/sqrt(49)),col='darkblue')
lines(range(x),range(x)*(1+2/sqrt(99)),col='darkgreen')
lines(range(x),range(x)*(1-2/sqrt(99)),col='darkgreen')


