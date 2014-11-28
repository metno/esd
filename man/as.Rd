\name{as}
\alias{as.station.data}
\alias{as.station.data.frame}
\alias{as.station.list}
\alias{as.station.ds}
\alias{as.station.pca}
\alias{as.station.field}
\alias{as.station.zoo}
\alias{as.station.spell}
\alias{as.station.eof}
\alias{as.ds}
\alias{as.field}
\alias{as.field.default}
\alias{as.field.zoo}
\alias{as.field.comb}
\alias{as.field.eof}
\alias{as.field.field}
\alias{as.field.station}
\alias{as.anomaly}
\alias{as.anomaly.default}
\alias{as.anomaly.zoo}
\alias{as.anomaly.station}
\alias{as.anomaly.field}
\alias{as.climatology}
\alias{as.eof}
\alias{as.eof.zoo}
\alias{as.eof.eof}
\alias{as.eof.comb}
\alias{as.eof.field}
\alias{as.eof.appendix}
\alias{as.observed}
\alias{as.annual}
\alias{as.annual.default}
\alias{as.annual.numeric}
\alias{as.annual.integer}
\alias{as.annual.yearqtr}
\alias{as.annual.station}
\alias{as.annual.spell}
\alias{as.monthly}
\alias{as.4seasons}
\alias{as.4seasons.default}
\alias{as.4seasons.station}
\alias{as.4seasons.day}
\alias{as.4seasons.field}
\alias{as.4seasons.spell}
\alias{as.seasons}
\alias{as.residual}
\alias{as.residual.ds}
\alias{as.residual.station}
\alias{as.original.data}
\alias{as.original.data.ds}
\alias{as.original.station}
\alias{as.appended}
\alias{as.appended.ds.comb}
\alias{as.appended.eof.comb}
\alias{as.appended.field.comb}
\alias{as.calibrationdata}
\alias{as.calibrationdata.ds}
\alias{as.calibrationdata.station}
\alias{as.fitted.values}
\alias{as.fitted.values.ds}
\alias{as.fitted.values.station}
\alias{as.pattern}
\alias{as.pattern.ds}
\alias{as.pattern.eof}
\alias{as.pattern.mvr}
\alias{as.pattern.cca}
\alias{as.pattern.trend}
\alias{as.pattern.field}
\alias{as.pattern.corfield}
\alias{as.pca}
\alias{as.pca.ds}
\alias{as.pca.station}
\alias{as.comb}
\alias{as.comb.eof}
\title{Conversion to esd objects.}
\description{
  Various methods for converting objects from one shape to
  another. These methods do the house keeping, keeping track of
  attributes and metadata.

  \code{as.field.station} uses \code{regrid} to generate a field based
  on bi-linear interpolation of station values and their coordinates. Unfinished...
}
\usage{
as.4seasons(x,FUN='mean',...)
as.4seasons.default(x,FUN='mean',...)
as.4seasons.station(x,FUN='mean',...)
as.4seasons.day(x,FUN='mean',na.rm=TRUE,dateindex=TRUE,...)
as.4seasons.field(x,FUN='mean',...)
as.4seasons.spell(x,FUN="mean",...)
as.seasons(x,start='01-01',end='12-31',FUN='mean', ...)
as.annual(x,...)
as.annual.default(x, ...)
as.annual.numeric(x, ...)
as.annual.integer(x, ...)
as.annual.yearqtr(x, frac = 0, ...)
as.annual.spell(x, ...)
as.annual.station(x, ...)
as.anomaly(x,...)
as.anomaly.default(x,ref=NULL,monthly=NULL,na.rm=TRUE)
as.anomaly.station(x,ref=NULL,monthly=NULL,na.rm=TRUE)
as.anomaly.field(x,ref=NULL,monthly=NULL,na.rm=TRUE)
as.anomaly.zoo(x,ref=NULL,monthly=NULL,na.rm=TRUE)
as.appended(x,...)
as.appended.ds.comb(x,iapp=1)
as.appended.eof.comb(x,iapp=1)
as.appended.field.comb(x,iapp=1)
as.eof(x,...)
as.eof.zoo(x,...)
as.eof.eof(x,iapp=NULL)
as.eof.comb(x,iapp=NULL)
as.eof.field(x,iapp=NULL,...)
as.eof.appendix(x,iapp=1)
as.calibrationdata(x)
as.calibrationdata.ds(x)
as.calibrationdata.station(x)
as.climatology(x,...)
as.comb(x,...)
as.comb.eof(x,...)
as.ds(x)
as.field(x,...)
as.field.zoo(x,lon,lat,param,unit,
             longname=NA,quality=NA,src=NA,url=NA,
             reference=NA,info=NA,calendar='gregorian',
             greenwich=TRUE, method= NA,type=NA,aspect=NA)
as.field.default(x,index,lon,lat,param,unit,
                 longname=NA,quality=NA,src=NA,url=NA,
                 reference=NA,info=NA,calendar='gregorian',
                 greenwich=TRUE, method= NA,type=NA,aspect=NA)
as.field.eof(x,...)
as.field.comb(x,iapp=NULL,...)
as.field.station(x,...)
as.fitted.values(x)
as.fitted.values.ds(x)
as.fitted.values.station(x)
as.monthly(x,fun='mean')
as.observed.station(x)
as.original.data(x)
as.original.data.ds(x)
as.original.data.station(x)
as.pattern(x)
as.pattern.ds(x)
as.pattern.eof(x)
as.pattern.cca(x)
as.pattern.mvr(x)
as.pattern.field(x)
as.pattern.trend(x)
as.pattern.corfield(x)
as.pca(x)
as.pca.ds(x)
as.pca.station(x)
as.residual(x)
as.residual.ds(x)
as.residual.station(x)
as.station(x,...)
as.station.zoo(x,loc=NA,param=NA,unit=NA,lon=NA,lat=NA,alt=NA,
                      cntr=NA,longname=NA,stid=NA,quality=NA,src=NA,url=NA,
                      reference=NA,info=NA, method= NA)
as.station.data.frame(x,location=NA,param=NA,unit=NA,lon=NA,lat=NA,alt=NA,
                      cntr=NA,longname=NA,stid=NA,quality=NA,src=NA,url=NA,
                      reference=NA,info=NA, method= NA)
as.station.zoo(x,location=NA,param=NA,unit=NA,lon=NA,lat=NA,alt=NA,
                      cntr=NA,longname=NA,stid=NA,quality=NA,src=NA,url=NA,
                      reference=NA,info=NA, method= NA,type=NA,aspect=NA)
as.station.list(x)
as.station.ds(x)
as.station.pca(x)
as.station.field(x)
as.station.spell(x)
as.station.eof(x,pattern=1:10)
as.pca.ds(x)
as.pca.station(x)
}
\arguments{
  \item{x}{Data object}
  \item{location}{define location attribute \code{attr(x,'location')}}
  \item{param}{define variable attribute \code{attr(x,'variable')}}
  \item{unit}{define unit attribute \code{attr(x,'unit')}}
  \item{lon}{define longitude attribute \code{attr(x,'longitude')}}
  \item{lat}{define latitude attribute \code{attr(x,'latitude')}}
  \item{alt}{define altitude attribute \code{attr(x,'altitude')}}
  \item{cntr}{define country attribute \code{attr(x,'country')}}
  \item{longname}{define long-name attribute \code{attr(x,'loongname')}}
  \item{stid}{define station ID attribute \code{attr(x,'station_id')}}
  \item{quality}{define quality attribute \code{attr(x,'quality')}}
  \item{src}{define source attribute \code{attr(x,'source')}}
  \item{url}{define URL attribute \code{attr(x,'URL')}}
  \item{reference}{define reference attribute \code{attr(x,'reference')}}
  \item{info}{define info attribute \code{attr(x,'info')}}
  \item{method}{define method attribute \code{attr(x,'method')}}
  \item{FUN}{function}
  \item{na.rm}{TRUE: ignore NA values}
  \item{dateindex}{}
  \item{monthly}{}
  \item{aspect}{}
  \item{iapp}{For values greater than 1, select the corresponding
    appended field in 'comb' objects (e.g. 1 gives
    \code{attr(x,'appendix.1')})}
  \item{pattern}{Which EOF pattern (mode) to extract as a series for PC}
 }
\value{A field object}

\seealso{
}

\author{R.E. Benestad and A.  Mezghanil}
\examples{
# Example: how to generate a new station object.
data <- round(matrix(rnorm(20*12),20,12),2); colnames(data) <- month.abb
x <- data.frame(year=1981:2000,data)
X <- as.station.data.frame(x,loc="",param="noise",unit="none")

# Example: how to generate a new field object.
year <- sort(rep(1991:2000,12))
month <- rep(1:12,length(1991:2000))
n <-length(year)
lon <- seq(-30,40,by=5); nx <- length(lon)
lat <- seq(40,70,by=5); ny <- length(lat)
# Time dimension should come first, space second.
y <- matrix(rnorm(nx*ny*n),n,nx*ny)
index <- as.Date(paste(year,month,1,sep="-"))
Y <- as.field(y,index,lon,lat,param="noise",unit="none")
map(Y)
plot(EOF(Y))

data(Oslo)
plot(as.anomaly(Oslo))

data(ferder)
plot(annual(ferder,FUN="min"))
plot(annual(ferder,FUN="IQR",na.rm=TRUE))
plot(as.4seasons(ferder))

data(bjornholt)
plot(annual(bjornholt,FUN="exceedance",fun="counts"))
plot(annual(bjornholt,FUN="exceedance",fun="freq"))
plot(annual(bjornholt,FUN="exceedance"))

# Test the as.4seasons function:
data(ferder)
#Daily data:
yd <- ferder
# Monthly data:
ym <- aggregate(ferder,as.yearmon)
ym <- zoo(coredata(ym),as.Date(index(ym)))
ym <- attrcp(ferder,ym)
plot(ym)

#Monthly reanalyses:
t2m <- t2m.ERAINT(lon=c(-30,40),lat=c(50,70))
T2m <- as.4seasons(t2m)
#Extract the grid point with location corresponding to that of the station:
x <- regrid(t2m,is=ferder)
x4s <- as.4seasons(x)
X4s <- regrid(T2m,is=ferder)
y4s1 <- as.4seasons(yd)
y4s2 <- as.4seasons(ym)
plot.zoo(y4s1,lwd=2,xlim=as.Date(c("1980-01-01","2000-01-01")),ylim=c(-10,20))
lines(y4s2,col="red",lty=2)
lines(x4s,col="darkblue",lwd=2)
lines(X4s,col="lightblue",lty=2)

# Select a random season
data(bjornholt)
data(ferder)
plot(as.seasons(ferder,FUN='CDD'))
plot(as.seasons(ferder,start='05-17',end='11-11',FUN='HDD'))
}
\keyword{utilities}
