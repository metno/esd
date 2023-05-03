## This script tests the extent daily and multi-day rainfall amounts follow the exponential distribution
## by plotting qq-plots against fitted distributions. The best normal distribution is also shown for comparison (grey)
## as the sum over larger daily samples is expected to resemble the normal distribution according to the 
## Central Limit Theorem. This test provides a follow-up on DOI: 10.1088/1748-9326/abd4ab
## Rasmus Benestad, 2023-05-03

library(esd)
## Use high-quality daily rain gauge data from Bjørnholt outside Oslo, Norway
data(bjornholt)
print(range(index(bjornholt)))
layout(matrix(1:9,3,3),widths=rep(1,2),heights = rep(1,3))

for (tau in c(1,2,3,5,10,15,30,50,90)) {
  y <- coredata(filter(bjornholt,rep(1,tau)/tau))
  y <- tau*y[seq(1,length(y),by=tau)]
  y <- y[y > 1]
  y <- y[is.finite(y)]
  mu <- mean(y)
  sigma <- sd(y)
  p <- seq(0,1,length=length(y))
  plot(sort(y),qexp(p=p,rate=1/mu),main=paste('Total over',tau,'days'))
  lines(c(0,1000),c(0,1000),lty=2,col='red')
  points(sort(y),qnorm(p=p,mean=mu,sd=sigma),pch=19,col='grey',cex=0.5)
  grid()
}

