\name{corfield}
\alias{corfield}
\alias{corfield.default}
\alias{corfield.zoo}
\alias{corfield.field}
\alias{corfield.field.station}
\alias{corfield.station}
\alias{corfield.eof}
\alias{corfield.trajectory}
\title{Correlation}
\description{
Compute the correlation between field objects and station/field.
}
\usage{
corfield(x,y,...)
corfield.zoo(x,y,plot=TRUE,use='pairwise.complete.obs',verbose=FALSE,new=TRUE,
             colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE))
corfield.field(x,y,plot=TRUE,use='pairwise.complete.obs',verbose=FALSE,new=TRUE,
               colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE),...)
corfield.field.station(x,y,plot=TRUE,verbose=FALSE,new=TRUE,
                       colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE),
                       use='pairwise.complete.obs',...)
corfield.station(x,y,plot=TRUE,verbose=FALSE,
                 use='pairwise.complete.obs',na.action='new=TRUE,na.omit',
                 colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE),...)
corfield.eof(x,y,ip=1,plot=TRUE,new=TRUE,
             colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE),
             use='pairwise.complete.obs',na.action='na.omit',...)
corfield.trajectory(x,y,it=NULL,is=NULL,param=NULL,FUN="count",
                    unit=NULL,longname=NULL,loc=NULL,
                    use="pairwise.complete.obs",
                    colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE),
                    method="pearson")
}
\arguments{
    \item{x}{data object}
    \item{y}{data object}
    \item{plot}{TRUE: plot the results}
    \item{use}{see \code{\link{cor}}}.
    \item{new}{see \code{link{map}}}.
    \item{ip}{index EOF pattern}
 }
\value{Map of correlation}
\author{R.E. Benestad and A. Mezghani}
\examples{
x <- t2m.ERAINT(lon=c(-40,30),lat=c(0,50))
y <- t2m.NCEP(lon=c(-40,30),lat=c(0,50))
r <- corfield(annual(x),annual(y))

data(Oslo)
t2m <- t2m.ERAINT()
x <- subset(Oslo,it='january')
y <- subset(t2m,it='january')
r <- corfield(x,y)
}
\keyword{manip}
