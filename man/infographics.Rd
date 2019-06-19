\name{InfoGraphics}
\alias{wheel}
\alias{wheel.station}
\alias{wheel.spell}
\alias{diagram}
\alias{diagram.dsensemble}
\alias{diagram.station}
\alias{seasevol}
\alias{seasevol.station}
\alias{vis}
\alias{vis.station}
\alias{diagram}
\alias{cumugram}
\alias{climvar}
\title{InfoGraphics.}
\alias{colscal}
\alias{visprob}
\alias{balls}
\alias{graph}
\alias{graph.dsensemble}
\alias{graph.zoo}
\alias{graph.list}
\alias{monthabsdiff}
%\alias{qq.test}
\description{
  Wheel

  Risk

  visprob - boxes with forseen outcomes - area proportional to probability

  conf - confidence intervals and uncertainty - clouds...

  vis

  diagram

  cumugram
  
  graph
  
  \code{monthabsdiff} - presents the magnitude of the interannual variations in monthly time series.
  
  
  \code{balls} draws 3D ball symbols.
}
\usage{
vis(x,...)
vis.station(x,new=FALSE,col=NULL,n=NULL,main=NULL,log.precip=TRUE,plot=TRUE,...)
diagram(x,...)
diagram.dsensemble(x,it=0,...)
diagram.station(x,it=NULL,new=TRUE,plot=TRUE,...)
wheel(x,...)
wheel.station(x,y=NULL,new=TRUE,lwd=2,col=NULL,type=NULL,bg="grey90",verbose=FALSE,...)
wheel.spell(x,y=NULL,new=TRUE,lwd=2,col=NULL,verbose=FALSE,...)
cumugram(x,it==NULL,start='-01-01',prog=FALSE,verbose=FALSE,FUN='mean',...)
climvar(x,FUN='sd',plot=TRUE,...)
colscal(n=30,col="bwr",alpha=NULL,test=FALSE)
seasevol(x,...)
seasevol.station(x,it=NULL,nv=25,...)
visprob(x,...) 
visprob.default(x,...)
visprob.station(x,y=NULL,is=1,...)
visprob.station.precip(x,y=NULL,is=1,threshold=1,...)
graph(x,...)
graph.default(x,img=NULL,pch='fancy',it=NULL,col=rgb(0.5,0.5,0.5,0.5),lwd=5,
                          xlim=NULL,ylim=NULL,new=TRUE,col.obs='black',...)
graph.dsensemble(x,img=NULL,pch='fancy',it=NULL,col=rgb(1,0.7,0.7,0.1),
                 lwd=5,xlim=NULL,ylim=NULL,add=FALSE,new=TRUE,ensmean=FALSE,
                 col.obs='black')
graph.list(x,img=NULL,pch='fancy',it=NULL,
           col=c(rgb(1,1,0.5,0.05),rgb(1,0.5,0.5,0.05),rgb(0.5,1,0.5,0.05)),
           lwd=5,xlim=NULL,ylim=NULL,add=FALSE,new=TRUE,ensmean=FALSE,
           col.obs='black')
graph.zoo(x,img=NULL,pch='fancy',it=NULL,col=rgb(1,0.7,0.7,0.1),
          lwd=5,xlim=NULL,ylim=NULL,xlab='',ylab='',add=FALSE,new=TRUE,ensmean=FALSE,
          col.obs='black')
balls(x,y=NULL,col=NULL,cex.max=2,n=20)
qq.test(x,...)
qp.test.station(x,...)
monthabsdiff(x,mon=1:12,FUN=NULL,verbose=FALSE,col=NULL)
}
\arguments{
 \item{x}{a data object}
 \item{new}{TRUE, invoke \code{dev.new()}}
 \item{col}{Colour palette}
 \item{col.obs}{Colour for the observations}
 \item{n}{number of colour scales}
 \item{main}{main title}
 \item{log.precip}{Use logarithmig scale for precipitation}
 \item{plot}{FALSE, just return the results}
 \item{verbose}{If TRUE, print out diagnostics}
 \item{img}{image from jpeg}
 \item{start}{starting date of the year}
 \item{ensmean}{Ensemble mean}
 \item{pch}{Plot symbol style}
 \item{it}{index time}
 \item{xlim}{see \code{\link{plot}}}
 \item{ylim}{see \code{\link{plot}}}
 \item{xlab}{see \code{\link{plot}}}
 \item{ylab}{see \code{\link{plot}}}
 \item{lwd}{see \code{\link{plot}}}
 \item{add}{TRUE: update existing plot}
 }
\value{A field object}

\seealso{
  \code{\link{map}}, \code{\link{plot.station}}, \code{\link{hist.spell}}  
}

\author{R.E. Benestad and A.  Mezghanil}
\examples{
data(bjornholt)
wheel(bjornholt)

z <- spell(bjornholt,threshold=1)
wheel(z)

}
\keyword{utilities}
