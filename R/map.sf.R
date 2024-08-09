## Map projections - a wrapper function that uses the sf-package and its projection facilities
## Map projections are: Goode projection, Universal Transverse Mercator projection, Mollweide projection,
## and the Robinson projection

map.sf <- function(x,...,FUN='mean',it=NULL,is=NULL,new=FALSE,
                   projection="+proj=moll",xlim=NULL,ylim=NULL,zlim=NULL,lab='default',
                   colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                                         show=TRUE,type="p",cex=2,h=0.6,v=1),
                   type=c("fill","contour"),gridlines=FALSE,cex=2,
                   lonR=NULL,latR=NULL,axiR=NULL,style='plain',
                   verbose=FALSE,plot=TRUE,add=FALSE) {
  if (verbose) print(paste0('map.sf: ',projection))
  require('sf'); require('oce')
  args <- list(...)
  main <- args$main
  data(coastlineWorld)
  par(mar=c(2.5, 1, 1.5, 1),bty='n')

  lon <- lon(x)
  lat <- lat(x)
  projection <- paste0(projection,' +lat_1=',min(lat),' +lat_2=',max(lat),
                       ' +lon_1=',min(lon),' +lon_2=',max(lon))
  if (verbose) print(projection)
  if (inherits(x,'field')) z <- apply(x,2,FUN=FUN) else z <- x
  if (!is.null(colbar)) {
    if (verbose) print('set colour scheme')
    colbar <- colbar.ini(x,FUN=FUN,colbar=colbar,verbose=FALSE)
  }
  if (verbose) print(c(length(z),length(lon),length(lat)))
  dim(z) <- c(length(lon),length(lat))
  cm <- colormap(z=z,col=colbar$col,breaks=colbar$breaks)
  drawPalette(colormap=cm,cex=0.75,plot=FALSE)
  mapPlot(coastlineWorld, projection=projection, grid=gridlines, col="lightgray",main=main)
  if (length(grep('fill',type))>0) mapImage(lon, lat, z, colormap=cm)
  if (length(grep('contour',type))>0) mapContour(lon, lat, z, col='black')
  mapLines(coastlineWorld, col="lightgray")
  
  if (verbose) print('colour legend')
  par(new=TRUE)
  image.plot(breaks=colbar$breaks, 
             lab.breaks=colbar$breaks,horizontal = TRUE,
             legend.only = TRUE, zlim = range(colbar$breaks),
             col = colbar$col, legend.width = 1,
             axis.args = list(cex.axis = 1,hadj = 0.5,mgp = c(0, 0.5, 0)), 
             border = FALSE, verbose=FALSE)
  title(sub = unit(x), line = 1, cex.sub = 0.7)
  par(new=FALSE)
  ## Return the data in case it should be processed further
  attr(z,'lon') <- lon
  attr(z,'lat') <- lat
  attr(z,'colbar') <- colbar
  invisible(z)
}