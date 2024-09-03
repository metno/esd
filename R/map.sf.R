## Map projections - a wrapper function that uses the sf-package and its projection facilities
## Map projections are: Goode projection, Universal Transverse Mercator projection, Mollweide projection,
## and the Robinson projection

## KMP 2024-08-19 This needs to be exported. Otherwise map.default can't find it. 
#' @export map.sf
map.sf <- function(x,...,FUN='mean',it=NULL,is=NULL,new=FALSE,
                   projection="+proj=moll",xlim=NULL,ylim=NULL,zlim=NULL,lab='default',
                   colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                                         show=TRUE,type="p",cex=2,h=0.6,v=1),
                   type=c("fill","contour"),gridlines=FALSE,cex=2,
                   lonR=NULL,latR=NULL,axiR=NULL,style='plain',
                   verbose=FALSE,plot=TRUE,add=FALSE) {
  if (verbose) print(paste0('map.sf: ',projection))
  ## KMP 2024-08-19: Packages should not be loaded within a function!
  ## It is better to check if the necessary packages are installed and stop with a warning if they are not.
  ## Then you can call the packages explicitly at use. That way it is also easier to see what the external packages are used for.
  ## Is sf actually used here?
  #require('sf'); require('oce')
  if(!requireNamespace("sf",quietly=TRUE)) {
    stop(paste0("Package \"sf\" needed to use the map projection", projection,". Please install it."))
  }
  if(!requireNamespace("oce",quietly=TRUE)) {
    stop(paste0("Package \"oce\" needed to use the map projection", projection,". Please install it."))
  }
  args <- list(...)
  main <- args$main
  data(coastlineWorld, package="oce")
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
  cm <- oce::colormap(z=z,col=colbar$col,breaks=colbar$breaks)
  oce::drawPalette(colormap=cm,cex=0.75,plot=FALSE)
  if(is.null(xlim) | is.null(ylim)) oce::mapPlot(coastlineWorld, projection=projection, 
                                                 grid=gridlines, col="lightgray", main=main) else 
    oce::mapPlot(coastlineWorld, projection=projection, grid=gridlines, col="lightgray", 
                 longitudelim=xlim, latitudelim=ylim, main=main)
  if (length(grep('fill',type))>0) oce::mapImage(lon, lat, z, colormap=cm)
  if (length(grep('contour',type))>0) oce::mapContour(lon, lat, z, col='black')
  oce::mapLines(coastlineWorld, col="lightgray")
  
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