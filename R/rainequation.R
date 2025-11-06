#' Rain equation
#' 
#' The rain equation \eqn{Pr(X>x)=fw exp(-x/mu)}{Pr(X>x)==f_w exp(-x/\mu)}
#' estimates the likelihood of 24-hr precipitation exceedint a threshold value
#' x. It is analogous to the normal distribution used to describe the
#' statistical distribution of e.g. daily temperature over a season, but
#' applied to precipitation. It has two parameters, the wet-day frequency fw
#' and the wet-day mean mu. It assumes that the distribution for wet-day
#' precipitation can be approximated with an exponential distribution, and has
#' one tail.
#' 
#' The function \code{rainvar} returns the dail variance of the 24-hr
#' precipitation according to \eqn{sigma^2 = 2 fw * mu^3}{\sigma^2 = 2 f_w
#' \mu^3} and \code{rainvartrend} calculates the first derivative accroding to
#' \eqn{d sigma^2/dt = 2 mu^3 dfw/dt + 6 fw mu^2 dmu/dt}{d \frac{\sigma^2}{dt}=
#' 2 \mu^3 \frac{d f_w}{dt} + 6 f_w * \mu^2 * \frac{d \mu}{dt}}. %% ~~ A
#' concise (1-5 lines) description of what the function does. ~~
#' 
#' 
#' @aliases rainequation fract.gt.x test.rainequation scatterplot.rainequation
#' rainvar rainvartrend
#'
#' @param x A station object - single station
#' @param x0 The threshold value defining an event.
#' @param threshold The threshold defining a 'wet day'.
#' @param src Data source
#' @param nmin Minimum number of years with data
#' @param verbose TRUE for printing out diagnostics
#' @param colour.by Used to plot the data points with differentcolours
#' according to e.g. 'x0', 'stid', 'alt', 'lon', or 'lat'
#' @param col colour palette
#'
#' @return a station object
#'
#' @references Benestad R. and A. Mezghani (2015), On downscaling probabilities
#' for heavy 24-hr precipitation events at seasonal-to-decadal scales, Tellus A
#' 2015, 67, 25954, http://dx.doi.org/10.3402/tellusa.v67.25954
#'
#' @examples
#' data(bjornholt)
#' plot(rainequation(bjornholt))
#' test.rainequation(bjornholt,threshold=30)
#' \dontrun{scatterplot.rainequation()}
#' 
#' @export
rainequation <- function(x,x0 = 10,threshold=NULL) {
  fw <- annual(x,FUN='wetfreq',threshold=threshold)
  mu <- annual(x,FUN='wetmean',threshold=threshold)
  pr.gt.x0 <- zoo(coredata(fw)*exp(-x0/coredata(mu)),order.by=index(fw))
  attr(pr.gt.x0,'variable') <- 'Pr(X>x)'
  attr(pr.gt.x0,'unit') <- 'probability'
  return(pr.gt.x0)
}

#' @export
fract.gt.x <- function(x,x0) {sum(x > x0,na.rm=TRUE)/sum(is.finite(x))}

#' @export
rainvar <- function(x,x0=1,na.rm=FALSE) {
  ## The variance estimated from the integral of the pdf assuming a threshold x0
  ## of 1 mm/day
  mu <- wetmean(x,threshold=x0)
  fw <- wetfreq(x,threshold=x0)
  sigma2 <- fw*(  mu^2 + (x0^2 + 2*x0*mu + 2*mu)*exp(-x0/mu) )
  return(sigma2)
}

#' @export
rainvartrend <- function(x,x0=1,na.rm=TRUE,nmin=NULL,verbose=FALSE) {
  ## The rate of change estimated as the first derivative from the analytic expression for sigma^2.
  if (verbose) {print('rainvartrend'); print(class(x))}
  if (verbose) print('wetmean')
  mu <- annual(x,FUN='wetmean',nmin=nmin,threshold=x0)
  if (verbose) print('wetfreq')
  fw <- annual(x,FUN='wetfreq',nmin=nmin,threshold=x0)
  if (verbose) print('first derivative')
  if (is.null(dim(x))) {
    m <- mean(mu,na.rm=TRUE)
    f <- mean(fw,na.rm=TRUE)
  } else { 
    m <- colMeans(mu,na.rm=TRUE)
    f <- colMeans(fw,na.rm=TRUE)
  }
  if (verbose) {print(summary(m)); print(summary(f))}
  ds2.dt <- (  m^2 + (x0^2 + 2*x0*m + 2*m)*exp(-x0/m) ) * trend.coef(fw) +
            f*( 2*m + (4*x0 + 4*m + x0^3/m^2 + 2*x0/m)*exp(-x0/m) ) * trend.coef(mu)
  if (verbose) {print('Final sigma2-trend estimate'); print(summary(ds2.dt))}
  return(ds2.dt)
}

## To test the rain equation
#' @export
test.rainequation <- function(loc='DE BILT',src='ecad',nmin=150,x0=20,
                              threshold=1,verbose=FALSE,plot=TRUE,new=TRUE) {
  
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
    ok <- is.finite(pr) & is.finite(obsfrac)
    r <- cor(pr[ok],obsfrac[ok])
    if (plot) {
      par(bty='n',xpd=TRUE)
      plot(pr,main=paste('The "rain equation" for',loc(y)),lwd=3,
           ylab=paste('fraction of days with more than',x0,'mm'),xlab='Year',new=new)
      
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
    attr(results,'correlation') <- r
    return(results)
}

## Use a scatter plot to evaluate the rain equation for a selection of rain gauge records.
## Select time series from e.g. ECA&D with a minimum number (e.g. 150) of years with data
#' @export
scatterplot.rainequation <- function(x='ecad',nmin=150,x0=c(10,20,30,40),
                                     type=c('image','contour','points'),plot=TRUE,
                                     threshold=1,colour.by='x0',col=NULL,verbose=FALSE) {
  if (verbose) print('scatterplot.rainequation')
  if (is.character(x)) {
    ss <- select.station(param='precip',nmin=150,x='ecad') 
    precip <- station(ss)
  } else if (is.station(x) & is.precip(x)) precip <- x
  
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
  X <- c(); Y <- X; r <- X; COL <- NULL
  
  par(bty='n',xpd=TRUE,mar=par()$mar + c(0,1,0,0))
  
  for (is in 1:d[2]) {
    y <- subset(precip,is=is)
    cat(loc(y),' ')
    # if (!is.null(colour.by)) {
    #   if (tolower(colour.by)=='stid') col <- cols[is]
    #   if (tolower(colour.by)=='lon') col <- cols[(1:d[2])[is.element(order(lon(precip)),lon(y))]]
    #   if (tolower(colour.by)=='lat') col <- cols[(1:d[2])[is.element(order(lat(precip)),lat(y))]]
    #   if (tolower(colour.by)=='att') col <- cols[(1:d[2])[is.element(order(alt(precip)),alt(y))]]
    # }
    for (itr in x0) {
    #   if (!is.null(colour.by)) 
    #     if (tolower(colour.by)=='x0') col <- cols[(1:length(x0))[is.element(x0,itr)]]
    #   if (is.null(COL)) COL <-col else COL <- c(COL,col)  
      pr <- rainequation(y,x0=itr)
      obsfrac <- annual(y,FUN='fract.gt.x',x0=itr)
      pr <- matchdate(pr,it=obsfrac); obsfrac <- matchdate(obsfrac,it=pr)
      X <- c(X,coredata(obsfrac)); Y <- c(Y,coredata(pr))
      ok <- is.finite(pr) & is.finite(obsfrac)
      r <- c(r,cor(pr[ok],obsfrac[ok]))
      rng <- range(c(X,Y),na.rm=TRUE)
      #if (verbose) cat(rng,' ')
      # if (sum(is.finite(rng))==2) {
      #   plot(X,Y,main='Test the "rain equation"',
      #        xlim=rng,ylim=rng,pch=19,cex=1.25,col=COL,
      #        xlab=expression(sum(H(X-x))/n),ylab=expression(Pr(X>x)==f[w]*e^{-x/mu}))
      #   grid()
      #   lines(c(0,0),rep(max(c(X,Y),na.rm=TRUE),2),lty=2,col=rgb(0.5,0.5,0.5,0.3))
      # }
    }
  }
  
  ## Make 2D bins and count events in each
  if (verbose) cat('2D bins \n')
  x <- X; y <- Y
  nx = 30; ny <- nx
  max.xy <- 1.1*max(c(x,y),na.rm=TRUE)
  seqx <- seq(0,max.xy,length=nx+1)
  seqy <- seq(0,max.xy,length=ny+1)
  seqX <- seq(0,max.xy,length=10*nx)
  seqY <- seq(0,max.xy,length=10*ny)
  if (verbose) { cat(seqx, '\n'); cat(seqy, '\n') }
  Z <- matrix(rep(0,nx*ny),nx,ny)
  for (i in 1:nx)
    for (j in 1:ny) Z[i,j] <- sum( (x >= seqx[i]) & (x < seqx[i+1]) &
                                   (y >= seqy[j]) & (y < seqy[j+1]), na.rm=TRUE )
  ## Make a smoother 2D surface with higher resolution
  ## Start with x dimension
  #image(Z)
  Zx <- matrix(rep(NA,10*nx*ny),length(seqX),ny)
  for (j in 1:ny) Zx[,j] <- approx(1:nx,Z[,j],xout=seq(1,nx,length=10*nx))$y
  Zxy <- matrix(rep(NA,100*nx*ny),length(seqX),length(seqY))
  for (i in seq(1,10*nx,by=1)) Zxy[i,] <- approx(1:ny,Zx[i,],xout=seq(1,ny,length=10*ny))$y
  
  ok <- is.finite(X) & is.finite(Y)
  r.comb <- round(cor(X[ok],Y[ok]),3)
  
  test.results <- list(Zxy=Zxy,x=X,y=Y)
  attr(test.results,'x') <- seqX
  attr(test.results,'y') <- seqY
  attr(test.results,'station') <- precip
  attr(test.results,'threshold') <- threshold
  attr(test.results,'x0') <- x0
  attr(test.results,'max.xy') <- max.xy
  attr(test.results,'col') <- COL
  attr(test.results,'colour.by') <- cols
  attr(test.results,'site_correlation') <- r
  attr(test.results,'combined_correlation') <- r.comb
  class(test.results) <- c('scatterplotrainequation','list')
  
  if (plot) plot(test.results)
  
  invisible(test.results)
}
