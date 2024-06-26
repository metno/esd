% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset.R
\name{subset.Default}
\alias{subset.Default}
\alias{subset}
\alias{subset.station}
\alias{subset.eof}
\alias{subset.pca}
\alias{subset.cca}
\alias{subset.events}
\alias{subset.trajectory}
\alias{subset.trend}
\alias{subset.corfield}
\alias{subset.ds}
\alias{subset.dsensemble}
\alias{subset.comb}
\alias{subset.field}
\alias{subset.spell}
\alias{subset.zoo}
\alias{station.subset}
\alias{subset.stationmeta}
\alias{subset.matrix}
\alias{subset.pattern}
\alias{subset.mvcomb}
\title{Subsetting esd objects}
\usage{
\method{subset}{Default}(x, it = NULL, is = NULL, ip = NULL, verbose = FALSE)
}
\arguments{
\item{x}{Data object from which the subset is taken}

\item{it}{A list or data.frame providing time index, e.g. a range of years like c(1979,2010), a season ('djf'), or a month ('dec' or 'december').}

\item{is}{A list or data.frame providing space index, e.g. a list of longitude and latitude range like list(lon=c(0,60), lat=c(35,60)).}

\item{ip}{selection of patterns in PCA or EOF (used for e.g. filtering the data)}

\item{verbose}{If TRUE, print out diagnosics}

\item{ensemble.aggregate}{If TRUE, call \code{subset.dsensemble.multi} if
appropriate.}

\item{ic}{Argument of \code{subset.events}: A list providing criteria for selection of cyclones, 
ic = list(param, pmax, pmin, FUN), where param is a parameter or element type,
pmax and pmin are the upper and lower limit of the parameter.  If FUN is
"any" (default setting), subset selects cyclones or trajectories that are
within the chosen range at any point during their lifetime.  If FUN is "all"
and x is an 'events' object, subset selects all individual cyclones within
the range (pmin, pmax).  If FUN is "all" and x is a 'trajectory' object,
subset selects cyclone trajectories that are within the chosen range at all
points during their lifetime.)}

\item{\dots}{additional arguments}
}
\value{
An object of the same class as the input object
}
\description{
The subset method tries to be 'intelligent', and if the list has no names,
then the list contains two vectors of length 2, then this is interpreted as
a region, e.g. argument \code{is = list(
c(lon.min,lon.max),c(lat.min,lat.max) )}. If, on the other hand, \code{is =
list( lon=1:50,lat=55:65 )}, then the function picks the longitudes and
latitudes which match these. This makes it flexible so that one can pick any
irregular sequence.
}
\examples{

data(Oslo)
# January months:
jan <- subset(Oslo,it="jan")
# The last 10 years:
recent <- subset(Oslo,it=c(2003,2012))
# JJA season
jja <- subset(Oslo,it="jja")
# Seasonl values for MAM
mam <- subset(as.4seasons(Oslo),it="mam")

data(ferder)
# Aggregated values for May
may <- subset(as.monthly(Oslo),it="may")
# The last 10 aggregated annual values 
recent.ann <- subset(as.annual(Oslo),it=2004:2013)

gcm <- t2m.NorESM.M()
# Extract July months from a field:
gcm.jul <- subset(gcm,it="jul")

# Extract a period from a field:
gcm.short <- subset(gcm.jul,it=c(1950,2030))

# Extract data for the region 0-50E/55-65N
X <-  subset(gcm,is=list(c(0,50),c(55,65)))

# Extract data for a specific set of longitudes and latitudes
Z <-  subset(gcm,is=list(lon=c(1,30),lat=c(58,63)))

t2m <- t2m.NCEP(lon=c(-10,30),lat=c(50,70))
cal <- subset(t2m,it=c("1948-01-01","1980-12-31"))

# Example on how to split the data into two parts for
# split-sample test...

T2M <- as.annual(t2m.NCEP(lon=c(-10,30),lat=c(50,70)))
cal <- subset(T2M,it=c(1948,1980))
pre <- subset(T2M,it=c(1981,2012))
comb <- combine(cal,pre) 
X <- EOF(comb)
plot(X)

data(ferder)
y <- as.annual(ferder)
z <- DS(y,X)
plot(z, new=FALSE)

# Test of subset the commutative property of subset and combine: 
T2M <- as.4seasons(t2m.NCEP(lon=c(-10,30),lat=c(50,70)))
GCM <- as.4seasons(t2m.NorESM.M(lon = range(lon(T2M))+c(-2,2), lat = range(lat(T2M))+c(-2,2)))
XY <- combine(T2M,GCM)
X1 <- subset(XY,it="mam")
X2 <- combine(subset(T2M,it="mam"),subset(GCM,it="mam"))
eof1 <- EOF(X1)
eof2 <- EOF(X2)
plot(merge(eof1[,1],eof2[,1]),plot.type='single',
     col=c('red','blue'),lty=c(1,2),lwd=c(4,2), new=FALSE)
# OK - identical results

# Extract storm tracks for specific periods, regions and characteristics
# from the sample 'events' object \code{storms} (North Atlantic storms identified from ERA5 data)
data(storms)

# Subset deep cyclones...
x <- subset(storms, ic=list(param="pcent", pmax=970, FUN="any"))
# ... and trajectories with a lifetime of at least 12 time steps (72 hours)
x <- subset(x, ic=list(param="trackcount", pmin=12))

# Subset cyclones in the region 10W-10E/55-65N
x.is <- subset(x,is=list(lat=c(55,65),lon=c(-10,10)))
# ...and all cyclones passing going through the region
x.is2 <- subset(x,it=which(x$trajectory \%in\% x.is$trajectory))

# Subset cyclones in the spring season (march, april, may)
x.mam <- subset(x, it="mam")
# Subset cyclones in december 2016
x.201612 <- subset(x,it=c("2016-12-01","2016-12-31")) 
map(x.201612, new=FALSE)
 
}
\seealso{
matchdate sort.station
}
\author{
R.E. Benestad and A.  Mezghani
}
\keyword{utilities}
