\name{DSensemble}
\alias{DSensemble}
\alias{DSensemble.default}
\alias{DSensemble.t2m}
\alias{DSensemble.precip}
\alias{DSensemble.annual}
\alias{DSensemble.season}
\alias{DSensemble.mu}
\alias{DSensemble.mu.worstcase}
\alias{DSensemble.pca}
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
(PCs) describing the \code{\link{EOF}s}. 

\code{DSensemble.pca} is used to downscale a predictor represented in
terms of PCA, and can reduce the computation time significantly. See
Benestad et al. (2015) \url{http://dx.doi.org/10.3402/tellusa.v67.28326}. 

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
                           predictor="ERA40_t2m_mon.nc",
                           rcp="rcp45",biascorrect=FALSE,
                           non.stationarity.check=FALSE,type='ncdf',
                           eofs=1:6,lon=c(-20,20),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                           select=NULL,FUN="mean",FUNX="mean",
                           pattern="tas_Amon_ens_",
                           path.ds=NULL,file.ds="DSensemble.rda",
                           nmin=NULL,verbose=FALSE)
DSensemble.precip((y,plot=TRUE,path="CMIP5.monthly/",
                              rcp="rcp45",biascorrect=FALSE,
                              predictor="ERA40_pr_mon.nc",
                              non.stationarity.check=FALSE,
                              type='ncdf',
                              eofs=1:6,lon=c(-10,10),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                              select=NULL,FUN="wetmean",
                              FUNX="sum",threshold=1,
                              pattern="pr_Amon_ens_",verbose=FALSE,nmin=NULL)
DSensemble.annual(y,plot=TRUE,path="CMIP5.monthly/",
                              rcp="rcp45",biascorrect=FALSE,
                              predictor="ERA40_t2m_mon.nc",
                              non.stationarity.check=FALSE,type='ncdf',
                              eofs=1:6,lon=c(-10,10),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                              abscoords=FALSE,select=NULL,FUN=NULL,
                              FUNX="mean",threshold=1,
                              pattern="tas_Amon_ens_",verbose=FALSE,nmin=NULL)
DSensemble.season(y,season="djf",plot=TRUE,path="CMIP5.monthly/",
                           predictor="slp.mon.mean.nc",
                           rcp="rcp45",biascorrect=FALSE,
                           non.stationarity.check=FALSE,type='ncdf',
                           eofs=1:6,lon=c(-20,20),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                           select=NULL,FUN="mean",FUNX="mean",
                           pattern="psl_Amon_ens_",
                           path.ds=NULL,file.ds=NULL,
                           nmin=NULL,verbose=FALSE)
DSensemble.mu(y,plot=TRUE,path="CMIP5.monthly/",
                          rcp="rcp45",biascorrect=FALSE,
                          predictor=list(t2m="data/ncep/air.mon.mean.nc",
                                         olr="data/ncep/OLR.mon.mean.nc",
                                         slp="data/ncep/slp.mon.mean.nc"),
                          non.stationarity.check=FALSE,type='ncdf',
                          eofs=1:16,lon=c(-30,20),lat=c(-20,10),it=NULL,rel.cord=TRUE,
                          select=NULL,FUN="wetmean",threshold=1,
                          pattern=c("tas_Amon_ens_","slp_Amon_ens_"),verbose=FALSE,nmin=365)
DSensemble.mu.worstcase(y,plot=TRUE,path="CMIP5.monthly/",
                                    predictor="ERA40_t2m_mon.nc",
                                    rcp="rcp45",biascorrect=FALSE,
                                    lon=c(-20,20),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                                    select=NULL,FUN="wetmean",type='ncdf',
                                    pattern="tas_Amon_ens_",verbose=FALSE)
DSensemble.pca(y,plot=TRUE,path="CMIP5.monthly/",
                          rcp="rcp45",biascorrect=FALSE,
                          predictor="ERA40_t2m_mon.nc",
                          non.stationarity.check=FALSE,
                          eofs=1:16,lon=c(-30,20),lat=c(-20,10), it=NULL,rel.cord=TRUE,
                          select=NULL,FUN="mean",rmtrend=TRUE,
                          FUNX="mean",threshold=1,type='ncdf',
                          pattern="tas_Amon_ens_",verbose=FALSE,
                          file.ds="DSensemble.rda",path.ds=NULL,nmin=NULL)
}
\arguments{
    \item{y}{A station object.}
    \item{plot}{Plot intermediate results if TRUE.}
    \item{path}{The path where the GCM results are stored.}
    \item{rcp}{Which (RCP) scenario}
    \item{biascorrect}{TRUE, apply a bias adjustment using \code{\link{biasfix}}}
    \item{predictor}{}
    \item{non.stationarity.check}{}
    \item{eofs}{Which EOFs to include in the step-wise multiple regresion}
    \item{rmtrend}{TRUE: detrend before calibrating the regression model.}
    \item{lon}{Longitude range for predictor}
    \item{lat}{Latitude range for predictor}
    \item{rel.cord}{TRUE: use the range relative to predictand; FALSE
use absolute range}
    \item{it}{Used to extract months or a time period. See \code{\link{subset}}.}
    \item{select}{GCMs to select, e.g .subsample the ensemble (1:3
selects the three first GCMs)}
    \item{FUN}{Function for aggregating the predictand (daily),
e.g. 'mean', 'wetmean'}
    \item{treshold}{Used together with FUN for some functions ('wetmean').}
    \item{nmin}{Minimum number of day used in \code{\link{annual}} used
for aggregating the predictand/predictor}
    \item{FUNX}{Function for transforming the predictor, e.g. '\code{\link{C.C.eq}}'
to estimate the saturation water vapout}
    \item{type}{Type of netCDF used in \code{\link{retrieve}} for
reading GCM data.}
    \item{pattern}{File name pattern for GCM data.}
    \item{verbose}{TRUE for checking and debugging the functions.}
    \item{file.ds}{Name of file savming the results.}
    \item{path.ds}{Path of file savming the results.}
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
rcp4.5 <- DSensemble.t2m(Oslo,path='~',rcp='',pattern = "tas_Amon_",
biascorrect=TRUE, predictor = "~/air.2m.mon.mean.nc",plot=TRUE)

# Evaluation: (1) combare the past trend with downscaled trends for same
# interval by ranking and by fitting a Gaussian to the model ensemble;
# (2) estimate the probabilty for the counts outside the 90
# percent confidence interval according to a binomial distribution.

diagnose(rcp4.5, plot = TRUE, type = "target")
}
\keyword{manip}
