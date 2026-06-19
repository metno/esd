visENSO <- function(x, verbose=TRUE) {
  if (verbose) cat('visENSO ',varid(x),'\n')
  x <- as.monthly(x)
  nino3.4 <- anomaly(NINO3.4()) 
  x <- matchdate(x,nino3.4)
  nino3.4 <- matchdate(nino3.4,x)
  if (verbose) {print(range(index(x))); print(range(index(nino3.4)))}
  yrs <- unique(year(x))
  breaks <- c(-Inf, pretty(coredata(nino3.4), n = 4), Inf)
  if (verbose) print(breaks)
  col <- cut(coredata(nino3.4), breaks = breaks,labels = FALSE)
  if (verbose) print(table(col))
  pal <- colscal(n=length(unique(col)),pal=paste(varid(x),'ipcc',sep='.'))
  if (verbose) print(pal)
  pal[2] <- '#888800'
  main <- switch(varid(x),'t2m'=paste(loc(x),'2m temperature'),
                 'precip'=paste(loc(x),'precipitation'))
  plot(month(x),coredata(x),col=pal[col],pch=19,
       main=main,ylab=esd::unit(x),xlab='',xaxt = "n")
  axis(side = 1, at = 1:12, labels = month.abb)
  legend('topleft',c('El Nino','warmish','Neutral','coolish','La Nina'),
         col=rev(pal),pch=19,bty='n')
  grid()
  for (yr in yrs) {
    z <- subset(x,it=rep(yr,2))
    lines(1:12,coredata(z),lty=3,col=pal[col[year(x)%in%yr]])
  }
}

#x <- as.monthly(station(stid=18700, param="t2m", src='metnod.thredds')) 
visENSO(x)