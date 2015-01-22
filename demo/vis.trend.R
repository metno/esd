# vis.trend.R
# K. Parding, January 2015
# functions for visual trend analysis

library(esd)
library(zoo)
library(Kendall)
library(ggplot2)
library(latticeExtra)

alltrends <- function(x,minlen=10){
  # Calculate trends of time series x
  # Can this be done faster using apply?
  year <- index(x)
  firstyear <- min(year):(max(year)-minlen+1)
  lastyear <- firstyear+minlen-1
  n <- length(firstyear)
  trends <- matrix(NA,n,n)
  colnames(trends) <- lastyear
  rownames(trends) <- firstyear
  p <- trends
  for (i in firstyear) {
    jvec <- (i+minlen-1):(max(year)+1)
    for (j in jvec) {
      ij <- which(year %in% i:j)
      ij.linfit <- lm(x[ij]~year[ij])
      ij.kendall <- Kendall(x[ij],year[ij])
      iout <- as.numeric(colnames(trends))==j
      jout <- as.numeric(rownames(trends))==i
      trends[jout,iout] <- ij.linfit$coefficients[2]
      p[jout,iout] <- ij.kendall$sl[1]
    }
  }  
  return(list("trends"=trends,"p"=p))
}

plotalltrends <- function(x,unitlabel="unit",varlabel="",
 alpha=0.05,minlen=10,lwd=NA,cticks=NA,new=TRUE) {
  T <- alltrends(x,minlen=minlen)
  trends <- T$trends*10
  p <- T$p
  cols <- colnames(trends)
  rows <- rownames(trends)
  significant <- ifelse(p<alpha,trends,NA)
  
  cmin <- rgb(239,138,98,max=255) # blue
  cmid <- rgb(247,247,247,max=255) # white
  cmax <- rgb(103,169,207,max=255) # red
  rgb.palette <- colorRampPalette(c(cmax,cmid,cmin),space="rgb")
  lattice.options(axis.padding=list(factor=0.5))

  ticks <- seq(1,length(cols),signif(length(cols)/10,1))
  if (is.na(lwd)) lwd <- max(3-0.05*length(cols),0.2)
    
  vmax <- max(abs(trends),na.rm=TRUE)
  vq <- q995(abs(trends))
  dv <- signif(vq/5,1)
  v0 <- 0#signif(dv/2,1)
  vstep <- seq(v0,signif(vq,1),dv)
  if (vmax>max(vstep)) vstep <- c(vstep,vmax)
  vstep <- unique(c(-1*vstep,vstep))
  vstep <- vstep[order(vstep)]

  if (new) dev.new()
  levelplot(t(trends),
    main=paste(c(varlabel," trend (",unitlabel,"/decade)"),collapse=""),
    xlab="end year",ylab="start year",
    scales=list(x=list(rot=45,at=ticks),y=list(at=ticks)),
    col.regions=rgb.palette(120),at=vstep,border.lwd=0) + 
    as.layer(
  levelplot(t(significant),col.regions=rgb.palette(120),at=vstep,
    border.lwd=lwd,border="black"))
 }

# EXAMPLES
# Temperature data from Oslo
data(Oslo)
plotalltrends(annual(Oslo),minlen=20,unitlabel='oC',
              varlabel='T2m Oslo',alpha=0.001)

djf <- annual(subset(Oslo,it='djf'))
mam <- annual(subset(Oslo,it='mam'))
jja <- annual(subset(Oslo,it='jja'))
son <- annual(subset(Oslo,it='son'))
plotalltrends(djf,minlen=20,unitlabel='oC',varlabel='T2m DJF',alpha=0.001)
plotalltrends(mam,minlen=20,unitlabel='oC',varlabel='T2m MAM',alpha=0.001)
plotalltrends(jja,minlen=20,unitlabel='oC',varlabel='T2m JJA',alpha=0.001)
plotalltrends(son,minlen=20,unitlabel='oC',varlabel='T2m SON',alpha=0.001)

# Random time series
N <- 40
y <- runif(N,0,1)
y <- zoo(y,order.by=1800:(1800+N))
plotalltrends(y,alpha=0.05)

