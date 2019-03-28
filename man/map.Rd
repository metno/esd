\name{map}
\alias{map}
\alias{map.default}
\alias{map.matrix}
\alias{map.data.frame}
\alias{map.station}
\alias{map.stationmeta}
\alias{map.stationsummary}
\alias{map.comb}
\alias{map.eof}
\alias{map.ds}
\alias{map.dsensemble}
\alias{map.field}
\alias{map.corfield}
\alias{map.cca}
\alias{map.events}
\alias{map.trajectory}
\alias{map.mvr}
\alias{map.pca}
\alias{map.trend}
\alias{lonlatprojection}
\alias{map.googleearth}
\alias{rotM}
\alias{gridbox}
\alias{map2sphere}
\alias{vec}
\alias{mask}

\title{Plot maps for esd objects}
\description{Make map of geophysical data. These plot functions are S3 methods for esd objects.
}
\usage{
map(x,it=NULL,is=NULL,new=TRUE,...)

map.default(x,FUN="mean",it=NULL,is=NULL,new=FALSE,
            projection="lonlat",xlim=NULL,ylim=NULL,zlim=NULL,n = 15,  
            colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05, 
            type="p",cex=2,show=TRUE,h=0.6,v=1,pos=0.05), 
            type=c("fill","contour"),gridlines=FALSE,lonR=NULL, 
            latR=NULL,axiR=NULL,verbose=FALSE,plot=TRUE, ...)

map.matrix(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
           xlim=NULL,ylim=NULL,zlim=NULL,n=15,
           colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                        pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
           type=c("fill","contour"),gridlines=FALSE,
           lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
           ip=1,plot=TRUE,...)
map.data.frame(x,...)
map.stationmeta(x,...)
map.stationsummary(x,FUN=NULL,cex=1,col='red',pal='t2m',pch=19,nbins=15,
                 new=TRUE,verbose=FALSE,fig=c(0.2,0.25,0.6,0.8),...)
           
map.station(x=NULL,FUN=NULL, it=NULL,is=NULL,new=FALSE,
                         add=FALSE,projection="lonlat",
                         xlim = NULL, ylim = NULL,zlim=NULL,n=15,
                         col='darkred',bg='orange',
                         colbar= list(pal='t2m',col=NULL,rev=FALSE,n=10,
                                      breaks=NULL,type="p",cex=2,h=0.6, v=1,
                                      pos=0.1,show=TRUE),
                         # col=NULL replaced by palette
                         type=NULL,gridlines=TRUE,
                         lonR=NULL,latR=45,axiR=NULL,verbose=FALSE,
                         cex=2,zexpr="alt",cex.subset=1,
                         add.text.subset=FALSE,showall=FALSE,
                         add.text=FALSE,main=NULL,sub=NULL,
                         height=NULL,width=NULL,
                         cex.main=1,cex.sub=0.75,cex.axis=1,cex.lab=0.6,
                         col.main="black",col.sub="grey",
                         font.main=1,font.sub=4,
                         pch=19, from=NULL,to=NULL,showaxis=FALSE,
                         border=FALSE,full.names=FALSE,
                         full.names.subset=FALSE, use.old=FALSE,
                         text=FALSE, fancy=FALSE, 
                         na.rm=TRUE,show.val=FALSE,usegooglemap=FALSE,
                         ##colorbar=TRUE,
                         xlab="lon",ylab="lat",
                         legend.shrink=1,fig=c(0,1,0.05,0.95),
                         mar=rep(2,4),mgp=c(3,1,0),plot=TRUE,...)

map.comb(x,it=NULL,is=NULL,new=TRUE,projection = "lonlat",
         xlim=NULL,ylim=NULL,zlim=NULL,n=15,
         colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                     type="p",cex=2,h=0.6,v=1),
         type=c("fill", "contour"),gridlines=FALSE,
         lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,ip=1,plot=TRUE,...)

map.eof(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",what="eof",
        xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
        colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                    pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
        type=c("fill","contour"),gridlines=FALSE,
        lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
        ip=1,cex=1,plot=TRUE,...)

map.ds(x,it=NULL,is=NULL,new=TRUE,projection="lonlat", 
       xlim=NULL,ylim=NULL,zlim=NULL,n=15,
       colbar=list(pal = NULL,rev=FALSE,n=10,breaks=NULL,
                   type="p",cex = 2,h=0.6,v=1,pos=0.05),
       type= c("fill", "contour"),gridlines=FALSE,
       lonR=NULL,latR=NULL,axiR = NULL,verbose = FALSE,plot=TRUE,...)

map.dsensemble(x,it=c(2000,2099),is=NULL,im=NULL,ip=NULL,
               colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                           show=TRUE,type="p",cex=2,h=0.6,v=1),
               FUN='mean',FUNX='mean',verbose=FALSE,anomaly=FALSE,test=FALSE,
               plot=TRUE,...)

map.field(x,FUN="mean",it=NULL,is=NULL,new=TRUE,projection = "lonlat",
            xlim=NULL,ylim=NULL,zlim=NULL, 
            colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                         type="p",cex=2,h=0.6,v=1), 
            type = c("fill", "contour"), gridlines = FALSE,
            lonR = NULL,latR = NULL, axiR = NULL, verbose = FALSE,
            na.rm = TRUE, plot=TRUE,...)

map.corfield(x,it=NULL,is=NULL,new=TRUE,projection= "lonlat", 
             xlim=NULL,ylim=NULL,zlim=NULL,n = 15,
             colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                         type="p",cex=2,h=0.6,v=1),
             type=c("fill", "contour"),gridlines=FALSE,
             lonR=NULL,latR=NULL,axiR = NULL,verbose = FALSE, ...)
 
map.trend(x,it=NULL,is=NULL,new=TRUE,projection="lonlat", 
          xlim=NULL,ylim=NULL,zlim=NULL,n=15,
          colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                      type="p",cex=2,h=0.6,v=1,),
          type=c("fill", "contour"),gridlines=FALSE,
          lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,plot=TRUE,...)

map.pca(x,it=NULL,is=NULL,ip=1,new=FALSE,projection="lonlat",
        xlim=NULL,ylim=NULL,zlim=NULL,FUN='mean',
        colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                    pos=0.05,show=TRUE,type="p",cex=1,h=0.6,v=1),
        type=c("fill","contour"),gridlines=FALSE,fig=c(0,1,0.05,0.95),
        lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,plot=TRUE,...)


map.mvr(x,it=NULL,is=NULL,new=TRUE,projection="lonlat", 
        xlim=NULL,ylim=NULL,zlim=NULL,n=15,
        colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                    type="p",cex=2,h=0.6,v=1),
        type=c("fill", "contour"),gridlines = FALSE, verbose = FALSE, plot=TRUE,...)


map.cca(x,icca=1,it=NULL,is=NULL,new=FALSE,projection="lonlat",
        xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
        colbar1=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,type="p",
                     cex=2,show=TRUE,h=0.6, v=1,pos=0.05), colbar2= NULL,
        type=c("fill","contour"),gridlines=FALSE,
        lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,cex=2,plot=TRUE,...)

map.trajectory(x,it=NULL,is=NULL,type="paths",
               projection="sphere",verbose=TRUE,plot=TRUE,...)

map.events(x,Y=NULL,it=NULL,is=NULL,xlim=NULL,ylim=NULL,main=NULL,
           param=NA,alpha=0.3,lwd=3,col="black",bg="white",pch=21,cex=1,
           colbar=list(pal="budrd",rev=FALSE,n=10,breaks=NULL,
                       pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
           showaxis=TRUE,fig=c(0,1,0.05,0.95),mgp=c(2,0.5,0),mar=rep(2,4),
           lty=1,type=c("points","trajectory","start","end"),border=FALSE,
           projection="lonlat",latR=NULL,lonR=NULL,new=TRUE,add=FALSE,
           verbose=FALSE,...)

lonlatprojection(x,it=NULL,is=NULL,new=FALSE,projection="lonlat", 
                 xlim=NULL,ylim=NULL,zlim=NULL,n =15,
                 colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                             type="p",cex=2,h=0.6,v=1),
                 type = c("fill", "contour"),geography=TRUE,
                 gridlines = FALSE, verbose = FALSE,fancy=FALSE,plot=TRUE,...) 

map.googleearth(x)

rotM(x=0,y=0,z=0)

gridbox(x,col,density = NULL, angle = 45)

map2sphere(x,it=NULL,is=NULL,new = TRUE,
           xlim=NULL,ylim=NULL,zlim = NULL,
           colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                        type="p",cex=2,h=0.6,v=1,pos=0.05),
           lonR = NULL, latR = NULL, axiR = 0, 
           type = c("fill", "contour"),colorbar=TRUE,
           gridlines = TRUE, fancy = FALSE, verbose = FALSE,...) 

vec(x,y,it=NULL,a=1,r=1,ix=NULL,iy=NULL,
    projection='lonlat',lonR=NULL,latR=NULL,axiR=0,verbose=FALSE,...)
mask(x,land=FALSE)
}
\arguments{
  \item{x}{the object to be plotted; in \code{rotM}, x holds a vector
    of x-coordinates.}
  \item{FUN}{The function to be applied on x before mapping (e.g. \code{mean})}
  \item{colbar}{The colour scales defined through \code{colscal}. Users can
                specify the colour `pal'*ette (`pal'), the number of
                breaks (`n'), values of `breaks', and the position of
                the color bar from the main plot (`pos'). The `rev' argument,
                will produce a reversed color bar if set to TRUE. The other
                arguments (`type',`h' and `v') are more specific to
		\code{col.bar} and are used only if argument `fancy' is
		set to TRUE (not yet finished).}
  \item{it}{see \code{\link{subset}}}
  \item{is}{see \code{\link{subset}}}
  \item{new}{TRUE: create a new graphic device.}
  \item{projection}{Projections: c("lonlat","sphere","np","sp") - the
    latter gives stereographic views from the North and south poles.}
  \item{xlim}{see  \code{\link{plot}} - only used for 'lonlat' and 'sphere' projections.}
  \item{ylim}{see  \code{\link{plot}} - only used for 'lonlat' and 'sphere' projections.}
  \item{n}{The number of colour breaks in the color bar}
  \item{breaks}{graphics setting - see \code{\link{image}}}
  \item{type}{graphics setting - colour shading or contour}
  \item{gridlines}{Only for the lon-lat projection}
  \item{lonR}{Only for the spherical projection used by \code{map2sphere} to change viewing angle}
  \item{latR}{Only for the spherical projection used by \code{map2sphere} to change viewing angle}
  \item{axiR}{Only for the spherical projection used by \code{map2sphere} to change viewing angle}
  \item{density}{}
  \item{y}{a vector of y coordinates}
  \item{z}{a vector of z coordinates}
  \item{ip}{Selects which pattern (see \code{\link{EOF}},  \code{\link{CCA}}) to plot}
  \item{geography}{TRUE: plot geographical features}
  \item{angle}{for hatching}
  \item{a}{used in \code{\link{vec}} to scale the length of the arrows}
  \item{r}{used in \code{\link{vec}} to make a 3D effect of plotting the
    arrows up in the air.}
  \item{ix}{used to subset points for plotting errors}
  \item{iy}{used to subset points for plotting errors}
  \item{colorbar}{Show the color bar in the map (default TRUE). If
    FALSE, the colorbar is not added into the map (ignored).}
  
  \item{cex.subset}{...}
  \item{add.text}{Add abbreviated location names.}
  \item{full.names}{Show the full name of the location.}
  \item{showall}{Default is set to FALSE}
  \item{showaxis}{If set to FALSE, the axis are not displayed in the plot.}
  \item{fancy}{If set to true, will use \code{\link{col.bar}} instead of
    \code{image} to produce the colour bar} 
  \item{text}{If TRUE, display text info on the map.The default is set
    to FALSE}
  \item{show.val}{Display the values of 'x' or 'FUN(x)' on top of the
    coloured map.}
  \item{legend.shrink}{If set, the size of the color bar is shrinked
    (values between 0 and 1)}
  \item{...}{further arguments passed to or from other methods.}
  \item{land}{if TRUE mask land, else mask ocean}
  \item{what}{What to map: ['eof','field] for EOF pattern or the field
  recovered from the EOFs.}
  \item{fig}{see \code{\link{par}}}
  \item{nbins}{number of bins/colour categories}
}
\value{A field object}

\seealso{
  \code{\link{plot.station}}
}

\author{R.E. Benestad}
\examples{
# Select stations in ss and map the geographical location 
# of the selected stations with a zoom on Norway.
ss <- select.station(cntr="NORWAY",param="precip",src="GHCND")
map(ss, col="blue",bg="lightblue",xlim = c(-10,30) , ylim = c(50,70))

## Get NACD data and map the mean values
y <- station.nacd()
map(y,FUN='mean',colbar=list(pal=varid(y),n=10),cex=2)

# Examples of cyclone maps (map.events)
data(storms)
# Subset cyclones from the start of January 2016 lasting a minimum of 10 time steps 
x <- subset(storms,it=c("2016-01-01","2016-01-07"),ic=list(param="trackcount",pmin=10))
# Map with points and lines showing the cyclone centers and trajectories
map(x, type=c("trajectory","points"), col="blue")
## Map with only the trajectory and start and end points
map(x, type=c("trajectory","start","end"), col="red")
## Map showing the cyclone depth (slp at center) as a color scale (rd = red scale)
map(x, param="pcent", type=c('trajectory','start'), 
    colbar=list(pal="rd", rev=TRUE, breaks=seq(980,1000,2)), alpha=0.9)

}
\keyword{map}


