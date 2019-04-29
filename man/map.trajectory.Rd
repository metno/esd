\name{map.trajectory}
\alias{map.trajectory}
\alias{sphere.trajectory}
\alias{lonlat.trajectory}

\title{Plot trajectory maps}
\description{Make different types of trajectory maps. Individual trajectories 
are mapped with map.trajectory. The number density can be visualised
  with map.hexbin.trajectory and map.sunflower.trajectory which are versions of
  \code{\link{scatter.hexbin}} and \code{\link{scatter.sunflower}}
  adapted to show trajectory density.}
\usage{
map.trajectory(x,it=NULL,is=NULL,type="paths",projection="sphere",
                  lonR=10,latR=90,
                  col='red',colmap='rainbow',alpha=0.3,pfit=FALSE,
                  main=NULL,xlim=NULL,ylim=NULL,new=TRUE)
lonlat.trajectory <- function(x,xlim=NULL,ylim=NULL,col='blue',alpha=0.1,
                  lty=1,lwd=1,main=NULL,new=TRUE)
sphere.trajectory <- function(x,xlim=NULL,ylim=NULL,col='blue',alpha=0.1,
                  lty=1,lwd=1,lonR=0,latR=90,main=NULL,new=TRUE)
map.density.trajectory(x,dx=4,dy=2,it=NULL,is=NULL,
      colbar=list(pal='precip',rev=TRUE,breaks=NULL,cex=2,h=0.6,v=1),
      projection='sphere',latR=90,lonR=10,gridlines=FALSE,...)
}
\arguments{
  \item{x}{the trajectory object to be plotted.}
  \item{it}{A list or data.frame providing time index, e.g. month}
  \item{is}{A list or data.frame providing space index, e.g. station record}
  \item{type}{type of map: 'paths' shows trajectories; 'density' shows
    the spatial density of the trajectories)}
  \item{col}{color of trajectories}
  \item{colmap}{Colour scales, either as an output from \code{\link{rbg}}
    or a single character string 'bwr' (blue-white-red) or 'rwb'
    ('red-white-blue')}
  \item{new}{TRUE: create a new graphics device}
  \item{projection}{Projections: c("lonlat","sphere","np","sp") - the
    latter gives stereographic views from the North and south poles.}
  \item{xlim}{see  \code{\link{plot}} - only used for 'lonlat' projection}
  \item{ylim}{see  \code{\link{plot}} - only used for 'lonlat' projection}
  \item{main}{an overall title for the plot}
  \item{lonR}{Only for the spherical projection - see \code{\link{map2sphere}}}
  \item{latR}{Only for the spherical projection - see
    \code{\link{map2sphere}}}
  \item{leg}{logical. If TRUE, legend is shown.}
  \item{alpha}{factor modifying the opacity alpha; typically in [0,1]}
}
\seealso{
  \code{\link{map} \link{map.events}}
}

\author{K. Parding}
\examples{
# plot storm tracks zoomed in on the north Atlantic and northern Europe
data(imilast.M03)
map.trajectory(imilast.M03,col="blue",alpha=0.1,
          projection='latlon',xlim=c(-60,60),ylim=c(30,90),
          new=FALSE)

# spherical projection
map.trajectory(imilast.M03,col="blue",alpha=0.1,projection='sphere',new=FALSE)

# plot number density for grid boxes of width 2 degrees and height 1 degree
map.hexbin.trajectory(imilast.M03,xlim=c(-60,60),ylim=c(30,90),dx=2,dy=1,new=FALSE)
map.sunflower.trajectory(imilast.M03,xlim=c(-60,60),ylim=c(30,90),dx=2,dy=1,new=FALSE)

\dontrun{
# calculate cyclone density, takes a little while
cdens <- as.field(imilast.M03)
map(cdens,new=FALSE)
}
}

\keyword{map}
\keyword{trajectory}


