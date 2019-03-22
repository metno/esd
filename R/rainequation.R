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

fract.gt.x <- function(x,x0) {sum(x > x0,na.rm=TRUE)/sum(is.finite(x))}

rainvar <- function(x,na.rm=FALSE) {
  ## The variance estimated from the integral of the pdf
  sigma2 <- 2*wetfreq(x)*wetmean(x)^3
  return(sigma2)
}

rainvartrend <- function(x,na.rm=TRUE,mean=TRUE,nmin=NULL) {
  ## The rate of change estimated as the first derivative from the analytic expression for sigma^2.
  mu <- annual(x,FUN='wetmean',nmin=nmin)
  fw <- annual(x,FUN='wetfreq',nmin=nmin)
  ds2.dt <- 2*mu^3 * trend.coef(fw) + 6*fw*mu^2*trend.coef(mu)
  if (mean) ds2.dt <- mean(ds2.dt,na.rm=na.rm)
  return(ds2.dt)
}

## To test the rain equation
test.rainequation <- function(loc='DE BILT',src='ecad',nmin=150,x0=20,
                              threshold=1,verbose=FALSE,plot=TRUE) {
  
  if (verbose) {print('test.rainequation'); print(c(x0,threshold))}
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
    if (verbose) print(loc(y))
    pr <- rainequation(y,x0=x0,threshold=threshold)
    obsfrac <- annual(y,FUN='fract.gt.x',x0=x0)
    counts <- annual(y,FUN='count',x0=x0)
    if (plot) {
      par(bty='n',xpd=TRUE)
      plot(pr,main=paste('The "rain equation" for',loc(y)),lwd=3,
           ylab=paste('fraction of days with more than',x0,'mm'),xlab='Year')
      
      lines(obsfrac,col=rgb(1,0,0,0.7),lwd=2)
      grid()
      legend(year(pr)[1],1.1*max(pr,na.rm=TRUE),
             c(expression(Pr(X>x)==f[w]*e^{-x/mu}),expression(sum(H(X-x))/n)),lty=1,lwd=c(3,2),
             col=c('black','red'),bty='n')
    }
    results <- merge(pr,obsfrac,counts)
    results <- attrcp(y,results)
    attr(results,'variable') <- c('probability','frequency','events')
    attr(results,'unit') <- c('fraction','fraction','count')
    return(results)
}

## Use a scatter plot to evaluate the rain equation for a selection of rain gauge records.
## Select time series from e.g. ECA&D with a minimum number (e.g. 150) of years with data
scatterplot.rainequation <- function(src='ecad',nmin=150,x0=c(10,20,30,40),
                                     threshold=1,colour.by='x0',col=NULL,verbose=FALSE) {
  if (verbose) print('scatterplot.rainequation')
  if (is.character(src)) {
    ss <- select.station(param='precip',nmin=150,src='ecad') 
    precip <- station(ss)
  } else if (is.station(src) & is.precip(src)) precip <- src
  
  d <- dim(precip)
  if (is.null(d)) d <- c(length(precip),1)
  if (verbose) print(d)
  
  if (!is.null(colour.by)) {
      if (is.character(('colour.by'))) {
         nc <- switch(tolower(colour.by),'x0'=length(x0),
                                         'stid'=dim(precip)[2],
                                         'lat'=length(lat(precip)),
                                         'lon'=length(lon(precip)),
                                         'alt'=length(alt(precip)))
         if (is.null(col)) cols <- colscal(nc,alpha=0.2)
      }
    }
  if (!is.null(colour.by))
        if (tolower(colour.by)=='x0') col <- cols
  if (is.null(col)) col <- rgb(0,0,0.7,0.15)
  
  firstplot <- TRUE
  X <- c(); Y <- X; COL <- NULL
  
  par(bty='n',xpd=TRUE,mar=par()$mar + c(0,1,0,0))
  
  for (is in 1:d[2]) {
    y <- subset(precip,is=is)
    print(loc(y))
    if (!is.null(colour.by)) {
      if (tolower(colour.by)=='stid') col <- cols[is]
      if (tolower(colour.by)=='lon') col <- cols[(1:d[2])[is.element(order(lon(precip)),lon(y))]]
      if (tolower(colour.by)=='lat') col <- cols[(1:d[2])[is.element(order(lat(precip)),lat(y))]]
      if (tolower(colour.by)=='att') col <- cols[(1:d[2])[is.element(order(alt(precip)),alt(y))]]
    }
    for (itr in x0) {
      if (!is.null(colour.by)) 
        if (tolower(colour.by)=='x0') col <- cols[(1:length(x0))[is.element(x0,itr)]]
      if (is.null(COL)) COL <-col else COL <- c(COL,col)  
      pr <- rainequation(y,x0=itr)
      obsfrac <- annual(y,FUN='fract.gt.x',x0=itr)
      pr <- matchdate(pr,it=obsfrac); obsfrac <- matchdate(obsfrac,it=pr)
      X <- c(X,coredata(obsfrac)); Y <- c(Y,coredata(pr))
      rng <- range(c(X,Y),na.rm=TRUE)
      if (verbose) print(rng)
      if (sum(is.finite(rng))==2) {
        plot(X,Y,main='Test the "rain equation"',
             xlim=rng,ylim=rng,pch=19,cex=1.25,col=COL,
             xlab=expression(sum(H(X-x))/n),ylab=expression(Pr(X>x)==f[w]*e^{-x/mu}))
        grid()
        lines(c(0,0),rep(max(c(X,Y),na.rm=TRUE),2),lty=2,col=rgb(0.5,0.5,0.5,0.3))
      }
    }
  }
  ok <- is.finite(X) & is.finite(Y)
  r <- round(cor(X[ok],Y[ok]),3)
  title(sub=paste('Correlation=',r))
  test.results <- list(x=X,y=Y)
  attr(test.results,'station') <- precip
  attr(test.results,'threshold') <- threshold
  attr(test.results,'col') <- COL
  return(test.results)
}
