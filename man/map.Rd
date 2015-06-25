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
map(x,it=NULL,is=NULL,new=TRUE,...)
map.default(x,it=NULL,is=NULL,new=TRUE,projection="lonlat",
              xlim=NULL,ylim=NULL,zlim=NULL,n=15,
              col=NULL,breaks=NULL,
              type=NULL,gridlines=FALSE,
              lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,...)

map.matrix(x,new=TRUE,projection="lonlat",...)

map.station(x,it=NULL,is=NULL,new=TRUE,projection = "lonlat",
              xlim = NULL, ylim = NULL, zlim = NULL,
              col = "darkred", bg = "orange",
              type = NULL,gridlines= TRUE, 
              lonR = NULL,latR=45,axiR=NULL,verbose=FALSE,
              FUN = NULL, 
              colbar = list(palette="t2m",rev=FALSE,n=10,breaks=NULL,
                            type="p",cex=2,h=0.6,v=1,pos=0.05),
              cex=0.8,zexpr = "alt",
              cex.subset = 1, add.text.subset = FALSE, showall = FALSE, 
              add.text = FALSE, height = NULL, width = NULL, cex.axis = 1, 
              cex.lab = 0.6, pch = 21, from = NULL, to = NULL,
              showaxis = FALSE, border = FALSE,
              full.names = FALSE, full.names.subset = FALSE, 
              text = FALSE, fancy = FALSE, na.rm = TRUE, show.val = FALSE, 
              colorbar = TRUE, legend.shrink = 1, ...)

map.comb(x,it=NULL,is=NULL,new=TRUE,projection = "lonlat",
           xlim=NULL,ylim=NULL,zlim=NULL,n=15,
           colbar=list(palette="t2m",rev=FALSE,n=10,breaks=NULL,type="p",
                       cex=2,h=0.6,v=1,pos=0.05),
           type=c("fill", "contour"),gridlines=FALSE,
           lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,pattern=1,...)

map.eof(x,it=NULL,is=NULL,new=TRUE,projection="lonlat",
           xlim= NULL,ylim=NULL,zlim=NULL,n=15,
           colbar=list(palette="t2m",rev=FALSE,n=10,breaks=NULL,type="p",
                       cex=2,h=0.6,v=1,pos=0.05),
           type=c("fill","contour"),gridlines=FALSE,
           lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,pattern=1,...)

map.ds(x,it=NULL,is=NULL,new=TRUE,projection="lonlat", 
         xlim=NULL,ylim=NULL,zlim=NULL,n=15,
         colbar=list(palette = "t2m",rev=FALSE,n=10,breaks=NULL,
                     type="p",cex = 2,h=0.6,v=1,pos=0.05),
         type= c("fill", "contour"),gridlines=FALSE,
         lonR=NULL,latR=-90,axiR = NULL,verbose = FALSE,...) 

map.field(x,FUN="mean",it=NULL,is=NULL,new=TRUE,projection = "lonlat",
            xlim=NULL,ylim=NULL,zlim=NULL, 
            colbar= list(palette="t2m",rev=FALSE,n=10,breaks = NULL,
                         type = "p", cex = 2, h = 0.6, v = 1, pos = 0.05), 
            type = c("fill", "contour"), gridlines = FALSE,
            lonR = NULL,latR = -90, axiR = NULL, verbose = FALSE,
            na.rm = TRUE, ...)

map.corfield(x,it=NULL,is=NULL,new=TRUE,projection= "lonlat", 
               xlim=NULL,ylim=NULL,zlim=NULL,n = 15,
               colbar=list(palette="t2m",rev=FALSE,n=10,breaks=NULL,
                           type="p",cex=2,h=0.6,v=1,pos=0.05),
               type=c("fill", "contour"),gridlines=FALSE,
               lonR=NULL,latR=-90,axiR = NULL,verbose = FALSE, ...)
 
map.trend(x,it=NULL,is=NULL,new=TRUE,projection="lonlat", 
            xlim=NULL,ylim=NULL,zlim=NULL,n=15,
            colbar=list(palette="t2m",rev=FALSE,n=10,breaks=NULL,type="p",
                        cex=2,h=0.6,v=1,pos=0.05),
            type=c("fill", "contour"),gridlines=FALSE,
            lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,...)

map.pca(x,it=NULL,is=NULL,pattern=1,new=TRUE,projection="lonlat", 
          xlim=NULL,ylim=NULL,zlim=NULL,n=15,
          colbar=list(palette=NULL,rev=FALSE,n=10,breaks=NULL,type="p",
                      cex=2,h=0.6,v=1,pos=0.05),
          type=c("fill","contour"),gridlines=FALSE,verbose = FALSE, ...)


map.mvr(x,it=NULL,is=NULL,new=TRUE,projection="lonlat", 
          xlim=NULL,ylim=NULL,zlim=NULL,n=15,
          colbar=list(palette="t2m",rev=FALSE,n=10,breaks=NULL,type="p",
                      cex=2,h=0.6,v=1,pos=0.05),
          type=c("fill", "contour"),gridlines = FALSE, verbose = FALSE, ...)


map.cca(x,it=NULL,is=NULL,new=TRUE,projection="lonlat",
          xlim=NULL,ylim=NULL,zlim=NULL,n=15,
          colbar=list(palette="t2m",rev=FALSE,n=10,breaks=NULL,type="p",
                      cex = 2,h=0.6,v=1,pos=0.05),
          type = c("fill","contour"),gridlines = FALSE,verbose=FALSE,...)

lonlatprojection(x,it=NULL,is=NULL,new=FALSE,projection="lonlat", 
                   xlim=NULL,ylim=NULL,zlim=NULL,n =15,
                   colbar=list(palette="t2m",rev=FALSE,n=10,breaks=NULL,
                               type="p",cex=2,h=0.6,v=1,pos=0.05),
                   type = c("fill", "contour"),geography=TRUE,
                   gridlines = FALSE, verbose = FALSE,fancy=FALSE,...) 

map.googleearth(x)

rotM(x=0,y=0,z=0)

gridbox(x,col,density = NULL, angle = 45)

map2sphere(x,it=NULL,is=NULL,new = TRUE,
             xlim=NULL,ylim=NULL,zlim = NULL,
             colbar= list(palette="t2m",rev=FALSE,n=10,breaks=NULL,
                          type="p",cex=2,h=0.6,v=1,pos=0.05),
             lonR = NULL, latR = NULL, axiR = 0, 
             type = c("fill", "contour"),
             gridlines = TRUE, fancy = FALSE, verbose = FALSE,...) 

vec(x,y,it=10,a=1,r=1,ix=NULL,iy=NULL,
                projection='lonlat',lonR=NULL,latR=NULL,axiR=0,...)
}
\arguments{
  \item{x}{the object to be plotted; in \code{rotM}, x holds a vector
    of x-coordinates.}
  \item{FUN}{The function to be applied on x before mapping (e.g. \code{mean})}
  \item{pattern}{Which EOF pattern (mode) to plot}
  \item{colbar}{The colour scales defined through \code{colscal}. Users can
                specify the colour palette (`palette'), the number of
                breaks (`n'), values of the `breaks', and the position of
                the color bar from the big plot (`pos'). The `rev' argument,
                will produce a reversed color bar if set to TRUE. The other
                arguments are more specific to \code{col.bar} and are
                used only if argument `fancy' is set to TRUE (not yet finished).}
  \item{it}{see \code{\link{subset}}}
  \item{is}{see \code{\link{subset}}}
  \item{new}{TRUE: create a new graphic device.}
  \item{projection}{Projections: c("lonlat","sphere","np","sp") - the
    latter gives stereographic views from the North and south poles.}
  \item{xlim}{see  \code{\link{plot}} - only used for 'lonlat' projection.}
  \item{ylim}{see  \code{\link{plot}} - only used for 'lonlat' projection.}
  \item{n}{The number of colour breaks in the color bar}
  \item{breaks}{graphics setting - see \code{\link{image}}}
  \item{type}{graphics setting - colour shading or contour}
  \item{gridlines}{Only for the lon-lat projection}
  \item{lonR}{Only for the spherical projection - see \code{\link{map2sphere}}}
  \item{latR}{Only for the spherical projection - see \code{\link{map2sphere}}}
  \item{axiR}{Only for the spherical projection - see \code{\link{map2sphere}}}
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
  \item{colorbar}{Show the color bar in the map (default TRUE). If
    FALSE, the colorbar is not added into the map (ignored).}
  
  \item{cex.subset}{}
  \item{add.text}{Add abbreviated location names.}
  \item{full.names}{Show the full name of the location.}
  \item{showall}{Default is set to FALSE}
  \item{showaxis}{}
  show.val = FALSE, 
  colorbar = TRUE, legend.shrink = 1
  \item{fancy}{If set to true, will use \code{col.bar} instead of
    \code{image} to produce the colour bar} 
  \item{text}{If TRUE, display text info on the map.The default is set
  to FALSE}
  \item{show.val}{Display the values of `x' or `FUN(x)' on top of the
  coloured map.}
\item {...}{further arguments passed to or from other methods.}	 
}
\value{A field object}

\seealso{
  \code{\link{plot.station}}
}

\author{R.E. Benestad}
\examples{
# Select stations in ss and map the geographical location of the selected stations with a zoom on Norway.
ss <- select.station(cntr="NORWAY",param="precip",src="GHCND")
map(ss, col="blue",bg="lightblue",xlim = c(-10,30) , ylim = c(50,70))

## Get NACD data and map the mean values
y <- station.nacd()
map(y,FUN='mean',colbar=list(palette=varid(y),n=10),cex=2)
}
\keyword{map}


