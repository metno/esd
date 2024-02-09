## Variability
## Show magnitude of year-to-year differences for each calendar month
require(esd)
require(RColorBrewer)

monthabsdiff <- function(y,mon=1:12,FUN=NULL,verbose=FALSE,col=NULL,sort=TRUE) { 
  
  if (is.null(col)) col <- c(brewer.pal(11,'Spectral'),'grey')
  if (is.null(FUN)) {
    if (is.precip(y)) FUN <- 'sum' else FUN <- 'mean'
  }
  x <- list()
  for (im in mon) {
    ym <- subset(as.monthly(y,FUN=FUN),it=month.abb[im])
    x[[month.abb[im]]] <- zoo(abs(diff(ym)),order.by=year(ym)[-1])
    if (im == mon[1]) x[['annual']] <- x[[month.abb[im]]] else
      x[['annual']] <- x[['annual']] + x[[month.abb[im]]]
  }
  
  if (verbose) str(x)
  
  plot(x[['annual']],type='n',xlab='',ylab=paste('magnitude of interannual monthly',varid(y)),
       main=paste(loc(y),stid(y),varid(y)),ylim=c(0,1.2*max(x[['annual']],na.rm=TRUE)))
  grid()
  t <- year(x[['annual']])
  
  legend(quantile(t,0.01),1.2*max(x[['annual']],na.rm=TRUE),x.intersp = 0.45,seg.len=1.1,
         month.abb[mon],lty=1,lwd=5,col=col[mon],bty='n',cex=0.75,horiz = TRUE,ncol=2)
  for (it in 1:length(t)) {
    y0 <- 0
    ## If sort, show the smallest magnitudes at the bottom and the largest at the top
    if (sort) {
      monyr <- c();
      for (im in mon) monyr <- c(monyr,round(coredata(x[[month.abb[im]]]),2)[it])
      mon <- mon[order(monyr)]
    }
    for (im in mon) {
      xm <- round(coredata(x[[month.abb[im]]]),2)
      if (verbose) print(c(it,im,t[it]-1,y0,t[it],y0+xm[it]))
      rect(t[it]-1,y0,t[it],y0+xm[it],col=col[im],border='black')
      y0 <- y0 + round(xm[it],2)
    }
  }
  #cal <- data.frame(y=coredata(x[['annual']]),x=t)
  #abline(lm(y ~ x,data=cal),lty=2)
}

y <- station(param='t2m',stid=18700,src='metnod.thredds')
monthabsdiff(y)
