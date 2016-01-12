\name{plot}
\alias{plot.station}
\alias{vis.map}
\alias{plot.eof}
\alias{plot.eof.field}
\alias{plot.eof.var}
\alias{plot.eof.comb}
\alias{plot.field}
\alias{plot.spell}
\alias{plot.cca}
\alias{plot.ds}
\alias{plot.dsx}
\alias{plot.dsensemble}
\alias{plot.diagnose}
\alias{plot.xval}
\alias{plot.diagnose.comb.eof}
\alias{plot.diagnose.matrix}
\alias{plot.diagnose.dsensemble}
\alias{plot.nevents}
\alias{nam2expr}
\title{Plot esd objects}
\description{
  These plot functions are S3 methods for esd objects, based on \code{plot}.
}
\usage{
nam2expr(x)
plot.station(x, plot.type = "single", new = TRUE, lwd = 3, type = "l", 
    pch = 0, main = NULL, col = NULL, xlim = NULL, ylim = NULL, 
    xlab = "", ylab = NULL, errorbar = TRUE, legend.show = FALSE, 
    map.show = TRUE, map.type = "points", mar = c(4.5, 4.5, 0.75, 
        0.5), alpha = 0.3, verbose = FALSE, ...) 
plot.eof(x,new=TRUE,xlim=NULL,ylim=NULL,pattern=1,what=c("pc","eof","var"),...)
plot.eof.field(x,new=TRUE,xlim=NULL,ylim=NULL,pattern=1,what=c("pc","eof","var"),...)
plot.eof.comb(x,new=TRUE,xlim=NULL,ylim=NULL,
                          pattern=1,col=c("red"),what=c("pc","eof","var"),...)
plot.eof.var(x,new=TRUE,xlim=NULL,ylim=NULL,pattern=20,...)
plot.ds(x,plot.type="single",what=c("map","ts"),new=TRUE,
        lwd=3,type='b',pch=0,main=NULL,col=NULL,
        xlim=NULL,ylim=NULL,xlab="",ylab=NULL,...)
plot.cca(x,icca=1,colbar1=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,type="p",cex=2,show=TRUE,
                        h=0.6, v=1,pos=0.05),colbar2=NULL,verbose=FALSE,...) 
plot.field(x,is=NULL,it=NULL,FUN="mean",...)
plot.spell(x)
plot.diagnose(x,...) 
plot.diagnose.comb.eof(x,xlim=NULL,ylim=NULL,verbose=FALSE,add=FALSE,...)
plot.diagnose.matrix(x,xlim=NULL,ylim=NULL,verbose=FALSE,...)
plot.diagnose.dsensemble(x,map=TRUE,verbose=FALSE,new=TRUE,cex=NULL)
plot.xval(x,...)
plot.dsensemble(x,pts=FALSE,showci=TRUE,showtrend=TRUE,legend=TRUE,it=0,
                envcol=c(1,0,0,0.2),...)
plot.nevents(x,verbose=FALSE,...)
vis.map(x,col='red',map.type='points')
}
\arguments{
  \item{x}{the object to be plotted}
  \item{pattern}{Which EOF pattern (mode) to plot}
  \item{col}{Colour}
  \item{icca}{Which CCA pattern to plot}
  \item{is}{For subsetting in space - See \code{link{subset}}, but can
    also be a station value and 
    if provided, the plotting will involve an interploation to the same
    coordinates as defined by \code{is}.}
  \item{it}{For subsetting in time - See \code{link{subset}}.'it=0'
    returns the annual means (mean of DJF + MAM + JJA + SON)}
  \item{FUN}{function}
  \item{map.show}{Show a small map of the location}
  \item{map.type}{}
 }
\value{A field object}

\seealso{\code{\link{plot}}
}

\author{R.E. Benestad}
\examples{
# Example: use aggregate to compute annual mean temperature for Svalbard:
data(Svalbard)
year <- as.numeric( format(index(Svalbard), '\%Y') ) 
y <- aggregate(Svalbard, by=year, FUN='mean', na.rm = FALSE) 
plot(y)


# Example with downscaling:
lon <- c(-12,37)
lat <- c(52,72)
t2m <- t2m.ERA40(lon=lon,lat=lat)
data(Oslo)
ds <- DS(Oslo,t2m)

# Plot the results for January month
# plot(subset(ds,it='Jan'))

# Plot the residuals:
residual <- as.residual(ds)
obs <- as.anomaly(as.calibrationdata(ds))

plot.zoo(obs,lwd=2)
lines(residual,col="red")

print("Global climate model simulation NorESM")
T2m <- t2m.NorESM.M(lon=lon,lat=lat)

# Plot the global mean of the field:
plot(T2m)
# Plot area mean of a sub region
plot(T2m,is=list(lon=c(0,10),lat=c(60,70)))

# Plot interpolated results corresponding to ferder
data(ferder)
plot(T2m,ferder)

# Plot Hovmuller diagram: Not working ...
## plot(T2m,is=list(lon=0)) 

print("Extract a subset - the January month")
x <- subset(t2m,it=1)
X <- subset(T2m,it=1)

print("Combine the fields for computing common EOFs:")
XX <- combine(x,X)

print("Compute common EOFs")
eofxx <- EOF(XX)
plot(eofxx)

print("Downscale the January mean temperature") 
ds.jan <- DS(Oslo,eofxx)
plot(ds.jan)

}
\keyword{hplot}


