## The rain equation. rasmus.benestad@met.no 2017-11-20
## Based on the assumption that the 24-hr accumulated precipitation approximately
## follows an exponential distribution for days with precipitation. The rain equation
## does not provide an accurate description of the extremes in the upper tail, but is
## more analogous to the normal distribution for seasonal temperature.

rainequation <- function(x,x0 = 10,threshold=NULL) {
  fw <- annual(x,FUN='wetfreq',threshold=threshold)
  mu <- annual(x,FUN='wetmean',threshold=threshold)
  pr.gt.x0 <- zoo(coredata(fw)*exp(-x0/coredata(mu)),order.by=index(fw))
  attr(pr.gt.x0,'variable') <- 'Pr(X>x)'
  attr(pr.gt.x0,'unit') <- 'probability'
  return(pr.gt.x0)
}

fract.gt.x <- function(x,threshold) {sum(x > threshold,na.rm=TRUE)/sum(is.finite(x))}

## To test the rain equation
test.rainequation <- function(loc='DE BILT',src='ecad',nmin=150,threshold=20) {
  
  if (is.null(loc)) {
    ss <- select.station(param='precip',nmin=150,src='ecad')
    Y <- station(ss)
    ## Pick a random case
    d <- dim(Y); pick <- order(rnorm(d[2]))[1]
    y <- subset(Y,is=pick)
  } else if (is.character(loc)) y <- station(param='precip',loc=loc,src=src) else
    if (inherits(loc,'station')) y <- loc
  
  d <- dim(y)
  if (!is.null(d)) y <- subset(y,is=1)
  pr <- rainequation(y,threshold=threshold)
  par(bty='n',xpd=TRUE)
  plot(pr,main=paste('The "rain equation" for',loc(y)),lwd=3,
       ylab=paste('fraction of days with more than',threshold,'mm'),xlab='Year')
  obsfrac <- annual(y,FUN='fract.gt.x',threshold=threshold)
  lines(obsfrac,col=rgb(1,0,0,0.7),lwd=2)
  grid()
  legend(year(pr)[1],1.1*max(pr,na.rm=TRUE),
         c(expression(Pr(X>x)==f[w]*e^{-x/mu}),expression(sum(H(X-x))/n)),lty=1,lwd=c(3,2),
       col=c('black','red'),bty='n')
}

## Use a scatter plot to evaluate the rain equation for a selection of rain gauge records.
## Select time series from e.g. ECA&D with a minimum number (e.g. 150) of years with data
scatterplot.rainequation <- function(src='ecad',nmin=150,threshold=c(10,20,30,40)) {
  
  ss <- select.station(param='precip',nmin=150,src='ecad')
  precip <- station(ss)

  d <- dim(precip)
  firstplot <- TRUE
  X <- c(); Y <- X
  for (is in 1:d[2]) {
    y <- subset(precip,is=is)
    print(loc(y))
    for (itr in threshold) {
      pr <- rainequation(y,threshold=itr)
      obsfrac <- annual(y,FUN='fract.gt.x',threshold=itr)
      pr <- matchdate(pr,it=obsfrac); obsfrac <- matchdate(obsfrac,it=pr)
      X <- c(X,coredata(obsfrac)); Y <- c(Y,coredata(pr))
      rng <- range(c(X,Y),na.rm=TRUE)
      par(bty='n',xpd=TRUE)
      plot(X,Y,main='Test the "rain equation"',
           xlim=rng,ylim=rng,pch=19,cex=1.25,col=rgb(0,0,0.7,0.15),
           xlab=expression(sum(H(X-x))/n),ylab=expression(Pr(X>x)==f[w]*e^{-x/mu}))
      grid()
      lines(c(0,0),rep(max(c(X,Y),na.rm=TRUE),2),lty=2,col=rgb(0.5,0.5,0.5,0.3))
    }
  }
  ok <- is.finite(X) & is.finite(Y)
  r <- round(cor(X[ok],Y[ok]),3)
  title(sub=paste('Correlation=',r))
  test.results <- list(x=X,y=Y)
  attr(test.results,'station') <- precip
  attr(test.results,'threshold') <- threshold
  return(test.results)
}
