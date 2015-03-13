\name{map.storm}
\alias{map.storm}
\alias{map.hexbin.storm}
\alias{map.sunflower.storm}
\alias{sphere.storm}
\alias{lonlat.storm}

\title{Plot storm track maps}
\description{Make different types of storm track maps. Individual storm
  tracks are mapped with map.storm. The number density can be visualised
  with map.hexbin.storm and map.sunflower.storm which are versions of
  \code{\link{scatter.hexbin}} and \code{\link{scatter.sunflower}}
  adapted to show storm tracks.}
\usage{
map.storm(x,it=NULL,is=NULL,projection="sphere",lonR=10,latR=90,
                  col='red',colmap='rainbow',alpha=0.3,pfit=FALSE,
                  main=NULL,xlim=NULL,ylim=NULL,new=TRUE)
lonlat.storm <- function(x,xlim=NULL,ylim=NULL,col='blue',alpha=0.1,
                  lty=1,lwd=1,main=NULL,new=TRUE)
sphere.storm <- function(x,xlim=NULL,ylim=NULL,col='blue',alpha=0.1,
                  lty=1,lwd=1,lonR=0,latR=90,main=NULL,new=TRUE)
map.hexbin.storm <- function(x,dx=6,dy=2,Nmax=NULL,
                  xgrid=NULL,ygrid=NULL,add=FALSE,leg=TRUE,
                  xlim=NULL,ylim=NULL,col='red',border='firebrick4',
                  colmap='heat.colors',scale.col=TRUE,scale.size=TRUE,
                  main=NULL,new=TRUE)
map.sunflower.storm <- function(x,dx=6,dy=2,petalsize=7,
                  xgrid=NULL,ygrid=NULL,leg=TRUE,leg.loc=2,
                  xlim=NULL,ylim=NULL,rotate=TRUE,alpha=0.6,
                  main=NULL,new=TRUE)
}
\arguments{
  \item{x}{the storm track object to be plotted.}
  \item{x}{col}
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
  \code{\link{map}},\code{\link{scatter.hexbin}},\code{\link{scatter.sunflower}}
}

\author{K. Parding}
\examples{
# plot storm tracks zoomed in on the north Atlantic and northern Europe
data(imilast.M03)
map.storm(imilast.M03,col="blue",alpha=0.1,
          projection='latlon',xlim=c(-60,60),ylim=c(30,90))

# spherical projection
map.storm(imilast.M03,col="blue",alpha=0.1,projection='sphere')

# plot number density for grid boxes of width 2 degrees and height 1 degree
map.hexbin.storm(imilast.M03,xlim=c(-60,60),ylim=c(30,90),dx=2,dy=1)
map.sunflower.storm(imilast.M03,xlim=c(-60,60),ylim=c(30,90),dx=2,dy=1)
}

\keyword{map}
\keyword{storm}


