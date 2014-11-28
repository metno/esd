\name{annual}
\alias{annual}
\alias{annual.zoo}
\alias{annual.default}
\alias{annual.dsensemble}
\alias{annual.station}
\alias{annual.spell}
\alias{annual.field}
\alias{year}
\alias{month}
\alias{day}
\alias{season}
\alias{seasonal.name}
\alias{pentad}
\title{Conversion to esd objects.}
\description{
  \code{annual} aggregates time series into annual values
  (e.g. means). \code{year},  \code{month}, \code{season} return the
  years, months, and days associated with the data.
}
\usage{
annual(x,...)
annual.zoo(x,FUN="mean",na.rm=TRUE,nmin=NULL, ...)
annual.default(x,FUN="mean",na.rm=TRUE,nmin=NULL, ...)
annual.dsensemble(x,FUN="mean")
annual.station(x,FUN="mean",nmin=NULL, ...)
annual.spell(x,FUN="mean",nmin=NULL,...)
annual.field(x,FUN="mean",na.rm=TRUE, nmin=NULL, ...)
year(x, ...)
month(x, ...)
day(x, ...)
season(x,format="numeric", ...)
season.name()
pentad(x,l=5,...)
}
\arguments{
  \item{x}{a station, field object, or a date}
  \item{FUN}{see \code{\link{aggregate.zoo}}}
  \item{nmin}{Minimum number of data points (e.g. days or months) with
    valid data accepted for annual estimate. NULL demands complete years.}
  \item{format}{'numeric' or 'character'}
  \item{na.rm}{TRUE: ignore NA - see see \code{\link{mean}}}
  \item{l}{length of window}
 }
\value{Same as x, or a numeric for \code{year}, \code{month},
  \code{day}, or \code{pentad}.
}

\seealso{
  \code{\link{as.annual}},\code{\link{aggregate.station}}
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
}
\keyword{utilities}
