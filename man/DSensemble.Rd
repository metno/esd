\name{DSensemble}
\alias{DSensemble}
\alias{DSensemble.default}
\alias{DSensemble.t2m}
\alias{DSensemble.precip}
\alias{DSensemble.annual}
\alias{DSensemble.mu}
\alias{DSensemble.mu.worstcase}
\title{Downscale ensemble runs}
\description{
Downscales an ensemble of climate model runs, e.g. CMIP5, taking the
results to be seasonal climate statistics. For temperature, the result
hold the seasonal mean and standard deviation, whereas for
precipitation, the results hold the wet-day mean, the wet-day frequency,
and the wet/dry-spell statistics. The call assumes that netCDF files
containing the climate model ensemble runs are stores in a file
structure, linked to the path argument and the rcp argument.

These methods are based on \code{\link{DS}}, and \code{DSensemble} is
designed to make a number of checks and evaluations in addition to
performing the DS on an ensemble of models. It is based on a similar
philosophy as the old R-package '\code{clim.pact}', but there is a new
default way of handling the predictors. In order to attempt to ensure a
degree of consistency between the downscaled results and those of the
GCMs, a fist covariate is introduced before the principal components
(PCs) describing the \code{\link{EOF}s}. The argument
\code{area.mean.expl=TRUE} will take the time series describing area
mean value for the selected predictor domain as the first covariate,
followed by the PCs. These are then used in the regression analysis. 

The argument \code{non.stationarity.check} is used to conduct an
additional test, taking the GCM results as 'pseudo-reality' where the
predictand is replaced by GCM results interpolated to the same location
as the provided predictand. The time series with interpolated values are
then used as predictor in calibrating the model, and used to predict
future values. This set of prediction is then compared with the
interpolated value itself to see if the dependency between the large and
small scales in the model world is non-stationary.

Other chekch include cross-validation (\code{\link{crossval}}) and
diagnostics comparing the sample of ensemble results with the
observations: number of observations outside the predicted 90-percent
conf. int and comparing trends for the past. 
}
\usage{
DSensemble(y,...)
DSensemble.default(y,path='CMIP5.monthly/',rcp='rcp45',...) 
DSensemble.t2m(y,plot=TRUE,path="CMIP5.monthly/",
                           predictor="~/ERA40_t2m_mon.nc",
                           rcp="rcp45",biascorrect=FALSE,
                           area.mean.expl=FALSE,
                           non.stationarity.check=FALSE,
                           eofs=1:6,lon=c(-15,15),lat=c(-10,10),
                           select=NULL,FUN="mean",FUNX="mean",
                           pattern="tas_Amon_ens_",verbose=FALSE)
DSensemble.precip(y,plot=TRUE,path="CMIP5.monthly/",
                  rcp="rcp45",biascorrect=FALSE,
                  predictor="~/ERA40_pr_mon.nc",
                  area.mean.expl=FALSE,
                  non.stationarity.check=FALSE,
                  eofs=1:6,lon=c(-15,15),lat=c(-10,10),
                  select=NULL,FUN="exceedance",
                  FUNX="sum",threshold=1,
                  pattern="pr_Amon_ens_",verbose=FALSE)
DSensemble.annual(y,plot=TRUE,path="CMIP5.monthly/",
                  rcp="rcp45",biascorrect=FALSE,
                  predictor="ERA40_t2m_mon.nc",
                  non.stationarity.check=FALSE,
                  area.mean.expl=FALSE,type='ncdf',
                  eofs=1:6,lon=c(-10,10),lat=c(-10,10),
                  abscoords=FALSE,select=NULL,FUN=NULL,
                  FUNX="mean",threshold=1,
                  pattern="pr_Amon_ens_",verbose=FALSE,nmin=NULL)
DSensemble.mu(y,plot=TRUE,path="CMIP5.monthly/",
              rcp="rcp45",biascorrect=FALSE,
              predictor="~/ERA40_t2m_mon.nc",
              non.stationarity.check=FALSE,
              eofs=1:16,lon=c(-30,20),lat=c(-20,10),
              select=NULL,FUN="wetmean",
              FUNX="C.C.eq",threshold=1,
              pattern="tas_Amon_ens_",verbose=FALSE)
DSensemble.mu.worstcase(y,plot=TRUE,path="~/CMIP5.monthly/",
                        predictor="~/ERA40_t2m_mon.nc",
                        rcp="rcp45",biascorrect=FALSE,
                        eofs=1:6,lon=c(-20,20),lat=c(-10,10),
                        select=NULL,FUN="wetmean",
                        pattern="tas_Amon_ens_",verbose=FALSE)
}
\arguments{
    \item{y}{A station object.}
    \item{plot}{Plot intermediate results if TRUE.}
    \item{path}{The path where the GCM results are stored.}
    \item{rcp}{Which (RCP) scenario}
    \item{area.mean.expl}{When TRUE, subtract the area mean for the domain
  and use as a the first co-variate before the PCs from the EOF analysis.}
 }
\value{A 'dsensembele' object - a list object holding DS-results.}
\author{R.E. Benestad and A. Mezghani}
\examples{
data(Oslo)
## Download NorESM1-M from 'climexp.knmi.nl' in default directory
## (home directory for linux/mac users)
url <-"http://climexp.knmi.nl/CMIP5/monthly/tas/"
## Download NorESM1-ME from the rcp45
noresm <- "tas_Amon_NorESM1-M_rcp45_000.nc"
if (!file.exists(noresm))
    download.file(url=file.path(url,noresm), destfile=noresm,
                 method = "auto",quiet=FALSE,mode = "w",cacheOK = TRUE)
fioesm <- "tas_Amon_FIO-ESM_rcp45_000.nc"
if (!file.exists(fioesm))
    download.file(url=file.path(url,fioesm), destfile=fioesm,
                  method = "auto", quiet = FALSE, mode = "w",cacheOK = TRUE)
## Downscale 2m temperature             
rcp4.5 <- DSensemble.t2m(Oslo,path='~',rcp='',pattern = "tas_Amon_",biascorrect=TRUE, area.mean.expl=TRUE, predictor = "~/air.2m.mon.mean.nc",plot=TRUE)

# Evaluation: (1) combare the past trend with downscaled trends for same
# interval by ranking and by fitting a Gaussian to the model ensemble;
# (2) estimate the probabilty for the counts outside the 90
# percent confidence interval according to a binomial distribution.

diagnose(rcp4.5, plot = TRUE, type = "target")
}
\keyword{manip}
