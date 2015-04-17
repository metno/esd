\name{spell}
\alias{spell}
\alias{spell.default}
\alias{spell.station}
\alias{hist.spell}
\alias{count}
\alias{wetfreq}
\alias{wetmean}
\alias{nevents}
\alias{exceedance}
\alias{exceedance.default}
\alias{exceedance.station}
\alias{exceedance.field}
\alias{HDD}
\alias{CDD}
\alias{GDD}
\alias{coldwinterdays}
\title{Spell statistics}
\description{Statistics of spell durations (consecutive wet and dry
  days), e.g. dry and wet periods or duration of extremes.

\code{exceedance} estimates statistics for peak-over-treshold, and
\code{nevents} returns the number of events with exceeding values
(e.g. the number of rainy days X > 1 mm/day). \code{wetfreq} resturs n
the wet-day frequency (a fraction) and \code{wetmean} wet-day mean.}
\usage{
spell(x,threshold,...)
spell.default(x,threshold,upper=150,...)
spell.station(x,threshold,upper=150,...)
hist.spell(x, ...)
count(x,threshold=1,fraction=FALSE) 
exceedance(x,threshold=1,fun='mean',...)
exceedance.default(x,threshold=1,fun='mean',...)
exceedance.station(x,threshold=1,fun='mean',...)
exceedance.field(x,threshold=1,fun='mean',...)
nevents(x,threshold=1)
wetfreq(x,threshold=1)
wetmean(x,threshold=1)
HDD(x,x0=18,na.rm=TRUE)
CDD(x,x0=22,na.rm=TRUE)
GDD(x,x0=10,na.rm=TRUE)
coldwinterdays(x,dse=NULL,threshold=0,verbose=FALSE,plot=TRUE)
}

\arguments{
  \item{x}{station or field object}
  \item{threshold}{threshold value}
  \item{upper}{upper limit for maximum length - ignore any above this
  because they are likely errornous}
  \item{fraction}{TRUE: divide the number of counts by number of samples}
  \item{fun}{function}
}
\value{Station or field objects}

\seealso{plot}

\author{R.E. Benestad and A.  Mezghanil}
\examples{
# Example 1 : 
precip <- station.metnod(stid="18700",param="precip")
x <- spell(precip,threshold=.1)
x.ann <- annual(x,FUN="max")
plot(x.ann,plot.type="multiple")
# Example 2 :
x11() ; plot(x)

# Growing degree days:
data(ferder)
plot(as.seasons(ferder,FUN='GDD'))

# Mild winter days - number of days in the winter season with
# above freezing temperatures
data(ferder)
try(coldwinterdays(ferder))
}
\keyword{utilities}
