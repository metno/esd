\name{map}
\alias{map}
\alias{map.default}
\alias{map.matrix}
\alias{map.station}
\alias{map.comb}
\alias{map.eof}
\alias{map.ds}
\alias{map.field}
\alias{map.corfield}
\alias{map.cca}
\alias{map.mvr}
\alias{map.pca}
\alias{map.trend}
\alias{lonlatprojection}
\alias{map.googleearth}
\alias{rotM}
\alias{gridbox}
\alias{map2sphere}
\alias{vec}

\title{Plot maps for esd objects}
\description{Make map of geophysical data. These plot functions are S3 methods for esd objects.
}
\usage{
map(x,it=NULL,is=NULL,new=FALSE,...)
map.default(x,FUN='mean',it=NULL,is=NULL,new=FALSE,
                        projection="lonlat",
                        xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                        colbar= list(palette=NULL,rev=FALSE,n=10,
                            breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                        type=c("fill","contour"),gridlines=FALSE,
                        lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,...)
map.matrix(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                     xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                     colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                     type=c("fill","contour"),gridlines=FALSE,
                     lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,
                     pattern=1,...)
map.comb(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                     xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                     colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                     type=c("fill","contour"),gridlines=FALSE,
                     lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,
                     pattern=1,...)
map.eof(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                    xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                    colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                    type=c("fill","contour"),gridlines=FALSE,
                    lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,
                    pattern=1,...)
map.ds(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                   colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                   type=c("fill","contour"),gridlines=FALSE,
                   lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,...)
map.field(x,FUN='mean',it=NULL,is=NULL,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                   colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                   type=c("fill","contour"),gridlines=FALSE,
                   lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,na.rm=TRUE,...)
map.corfield(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                   colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                   type=c("fill","contour"),gridlines=FALSE,
                   lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,...)
map.trend(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                   colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                   type=c("fill","contour"),gridlines=FALSE,
                   lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,...)
map.pca(x,it=NULL,is=NULL,pattern=1,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                   colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1),
                   type=c("fill","contour"),gridlines=FALSE,verbose=FALSE,...)
map.mvr(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                   colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                   type=c("fill","contour"),gridlines=FALSE,verbose=FALSE,...)
map.cca(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                   colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1),
                   type=c("fill","contour"),gridlines=FALSE,verbose=FALSE,...)
lonlatprojection(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                   colbar= list(palette=NULL,rev=FALSE,n=10,
                             breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.05),
                   type=c("fill","contour"),gridlines=FALSE,verbose=FALSE,
                   geography=TRUE,fancy=FALSE,...)
map.googleearth(x)
rotM(x=0,y=0,z=0)
gridbox(x,col,density = NULL, angle = 45)
map2sphere(x,it=NULL,is=NULL,lonR=NULL,latR=NULL,axiR=0,new=FALSE,
                       type=c("fill","contour"),colorbar=TRUE,breaks=NULL,
                       gridlines=TRUE,col=NULL,...)
vec(x,y,it=10,a=1,r=1,ix=NULL,iy=NULL,
                projection='lonlat',lonR=NULL,latR=NULL,axiR=0,...)
}
\arguments{
  \item{x}{the object to be plotted; in \code{rotM} x holds a vector
    of x-coordinates.}
  \item{pattern}{Which EOF pattern (mode) to plot}
  \item{col}{Colour scales, either as an output from \code{\link{rbg}}
    or a single character string 'bwr' (blue-white-red) or 'rwb'
    ('red-white-blue')}
  \item{it}{see \code{\link{subset}}}
  \item{is}{see \code{\link{subset}}}
  \item{new}{TRUE: create a new graphics device}
  \item{projection}{Projections: c("lonlat","sphere","np","sp") - the
    latter gives stereographic views from the North and south poles.}
  \item{xlim}{see  \code{\link{plot}} - only used for 'lonlat' projection}
  \item{ylim}{see  \code{\link{plot}} - only used for 'lonlat' projection}
  \item{n}{graphics setting - number of colour breaks}
  \item{breaks}{graphics setting - see \code{\link{image}}}
  \item{type}{graphics setting - colour shading or contour}
  \item{gridlines}{Only for the lon-lat projection}
  \item{lonR}{Only for the spherical projection - see \code{\link{map2sphere}}}
  \item{latR}{Only for the spherical projection - see \code{\link{map2sphere}}}
  \item{axiR}{Only for the spherical projection - see\code{\link{map2sphere}}}
  \item{density}{}
  \item{y}{a vector of y coordinates}
  \item{z}{a vector of z coordinates}
  \item{pattern}{Selects which pattern (see \code{\link{EOF}},  \code{\link{CCA}}) to plot}
  \item{geography}{TRUE: plot geographical features}
  \item{angle}{for hatching}
  \item{a}{used in \code{\link{vec}} to scale the length of the arrows}
  \item{r}{used in \code{\link{vec}} to make a 3D effect of plotting the
    arrows up in the air.}
  \item{ix}{used to subset points for plotting errors}
  \item{iy}{used to subset points for plotting errors}
  \item{colorbar}{TRUE: plot colorbar}
 }
\value{A field object}

\seealso{
  \code{\link{plotstation}}
}

\author{R.E. Benestad}
\examples{
# select stations in ss and map the geographical location of the selected stations with a zoom on Norway.
ss <- select.station(cntr="NORWAY",param="precip",src="GHCND")
map(ss, col="blue",bg="lightblue",xlim = c(-10,30) , ylim = c(50,70))

## Get NACD data and map the mean values
y <- station.nacd()
map(y,FUN='mean',cex=2)
}
\keyword{map}


