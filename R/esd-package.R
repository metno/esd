# Conversion to esd objects.
# 
# Various methods for converting objects from one shape to another. These
# methods do the house keeping, keeping track of attributes and metadata.
# 
# \code{as.field.station} uses \code{regrid} to generate a field based on
# bi-linear interpolation of station values and their coordinates.
# Unfinished...
# 
# `field' object must be a two dimensional object (time,nlon*nlat). Time index
# must be the first dimension. The second dimension is a product of the number
# of longitudes by the number of latitudes. For a different format, users can
# use `aperm' function (Transpose an array by permuting its dimensions) to
# rearrange the object to meet the required format and `dim' function can be
# used to transform a 3D (d1,d2,d3) object into 2D (d1,d2*d3) object as
# dim(3D) <- c(d1,d2*d3)
# 
# @aliases as.station.data as.station.data.frame as.station.list as.station.ds
# as.station.pca as.station.field as.station.zoo as.station.spell
# as.station.eof as.station.dsensemble as.station.dsensemble.pca as.ds
# as.field as.field.default as.field.zoo as.field.comb as.field.eof
# as.field.ds as.field.field as.field.station as.field.dsensemble.eof
# as.anomaly as.anomaly.default as.anomaly.zoo as.anomaly.station
# as.anomaly.field as.climatology as.eof as.eof.zoo as.eof.eof as.eof.ds
# as.eof.comb as.eof.field as.eof.appendix as.annual as.annual.default
# as.annual.numeric as.annual.integer as.annual.yearqtr as.annual.station
# as.annual.spell as.monthly as.4seasons as.4seasons.default
# as.4seasons.station as.4seasons.day as.4seasons.field as.4seasons.spell
# as.seasons as.residual as.residual.ds as.residual.station as.original.data
# as.original.data.ds as.original.station as.appended as.appended.ds.comb
# as.appended.eof.comb as.appended.field.comb as.calibrationdata
# as.calibrationdata.ds as.calibrationdata.station as.fitted.values
# as.fitted.values.ds as.fitted.values.station as.pattern as.pattern.ds
# as.pattern.eof as.pattern.mvr as.pattern.cca as.pattern.trend
# as.pattern.field as.pattern.corfield as.pca as.pca.ds as.pca.station as.comb
# as.comb.eof as.trajectory
# @param x Data object
# @param location define location attribute \code{attr(x,'location')}
# @param param define variable attribute \code{attr(x,'variable')}
# @param unit define unit attribute \code{attr(x,'unit')}
# @param lon define longitude attribute \code{attr(x,'longitude')}
# @param lat define latitude attribute \code{attr(x,'latitude')}
# @param alt define altitude attribute \code{attr(x,'altitude')}
# @param cntr define country attribute \code{attr(x,'country')}
# @param longname define long-name attribute \code{attr(x,'loongname')}
# @param stid define station ID attribute \code{attr(x,'station_id')}
# @param quality define quality attribute \code{attr(x,'quality')}
# @param src define source attribute \code{attr(x,'source')}
# @param url define URL attribute \code{attr(x,'URL')}
# @param reference define reference attribute \code{attr(x,'reference')}
# @param info define info attribute \code{attr(x,'info')}
# @param method define method attribute \code{attr(x,'method')}
# @param FUN function
# @param na.rm TRUE: ignore NA values
# @param nmin Minimum number of valid data points. Default value is set to
# nmin=85 for daily values, nmin=3 for monthly values, and nmin=300 for annual
# values.
# @param dateindex if TRUE, convert format of index from year quarter to date
# @param aspect decription of data type, e.g., original, anomaly, fitted, or
# predicted
# @param iapp For values greater than 1, select the corresponding appended
# field in 'comb' objects (e.g. 1 gives \code{attr(x,'appendix.1')})
# @param ip Which EOF pattern (mode) to extract as a series for PC
# @return A field object
# @author R.E. Benestad and A.  Mezghani
# @keywords utilities
# @examples
# 
# # Example: how to generate a new station object.
# data <- round(matrix(rnorm(20*12),20,12),2)
# colnames(data) <- month.abb
# x <- data.frame(year=1981:2000,data)
# X <- as.station.data.frame(x,loc="",param="noise",unit="none")
# 
# # Example: how to generate a new field object.
# year <- sort(rep(1991:2000,12))
# month <- rep(1:12,length(1991:2000))
# n <-length(year)
# lon <- seq(-30,40,by=5)
# nx <- length(lon)
# lat <- seq(40,70,by=5)
# ny <- length(lat)
# # Time dimension should come first, space second.
# y <- matrix(rnorm(nx*ny*n),n,nx*ny)
# index <- as.Date(paste(year,month,1,sep="-"))
# Y <- as.field(y,index=index,lon=lon,lat=lat,param="noise",unit="none")
# map(Y)
# plot(EOF(Y))
# 
# data("Oslo")
# plot(as.anomaly(Oslo))
# 
# data(ferder)
# plot(annual(ferder,FUN="min"))
# plot(annual(ferder,FUN="IQR",na.rm=TRUE))
# plot(as.4seasons(ferder))
# 
# data(bjornholt)
# plot(annual(bjornholt,FUN="count",threshold=1))
# plot(annual(bjornholt,FUN="wetfreq",threshold=1))
# plot(annual(bjornholt,FUN="wetmean",threshold=1))
# 
# ## Test the as.4seasons function:
# data("ferder")
# ## Daily data:
# yd <- ferder
# ##  Monthly data:
# ym <- as.monthly(ferder, FUN="mean")
# plot(ym)
# 
# ## Monthly reanalyses:
# t2m <- t2m.NCEP(lon=c(-30,40),lat=c(50,70))
# T2m <- as.4seasons(t2m)
# # Extract the grid point with location corresponding to that of the station:
# x <- regrid(t2m,is=ferder)
# x4s <- as.4seasons(x)
# X4s <- regrid(T2m,is=ferder)
# y4s1 <- as.4seasons(yd)
# y4s2 <- as.4seasons(ym)
# plot(y4s1,lwd=2,xlim=as.Date(c("1980-01-01","2000-01-01")),ylim=c(-10,20))
# lines(y4s2,col="red",lty=2)
# lines(x4s,col="darkblue",lwd=2)
# lines(X4s,col="lightblue",lty=2)
# 
# ## Select a random season
# data("bjornholt")
# data("ferder")
# plot(as.seasons(ferder,FUN='mean'))
# plot(as.seasons(ferder,start='05-17',end='11-11',FUN='mean'))
# 
# \dontrun{
# ## PCA-example - switch back and forth between station and PCAs
# t2m <- station.narp()
# pca.t2m <- as.pca(annual(t2m))
# y <- as.station(pca.t2m)
# 
# ## Quick example of downscaling PCA and then extract the single station series
# predictor <- annual(t2m.NCEP(lon=c(-20,30),lat=c(50,80)))
# pca.t2m <- subset(pca.t2m,ip = 1:3)
# z.pca <- DSensemble.pca(pca.t2m,predictor=predictor,select=1:4)
# z <- as.station.dsensemble.pca(z.pca)
# plot(z[[1]])
# }
# 
# ## Example of transforming cyclone data from an 'events' object (individual events)
# ## into a 'trajectory' object (interpolated trajectories) which takes less space
# data(storms)
# x <- as.trajectory(storms)
# map(x, alpha=0.1)





# Sample data.
# 
# Different data sets: station data from northern Europe (NACD, NARP) and
# historic reconstructions (Oslo, Svalbard) from Dr. Nordli, Met Norway.
# 
# The object \code{station.meta} contains station information, used in the
# methods \code{\link{station}}.
# 
# Also reduced representation of re-analyses, where the data have been sampled
# by skipping grid points to reduce the spatial dimensions and stored as 20
# EOFS (30 for precipitation). The data compression facilitated by the EOFs
# can provide 80-90\% of the variance in the data. ESD uses the large-scale
# features from these reanalyses, and hence this information loss may be
# acceptable for downscaling work.
# 
# A reduced copy of the NorESM (M RCP 4.5) is also provided for the examples
# and demonstrations on how the downscaling can be implemented. Note:
# downscaling for end-users should never be based on one GCM simulation alone.
# 
# The object \code{geoborders} contains data on coastlines and borders, used
# in the methods \code{\link{map}}.
# 
# \code{glossstations} contains META-data for GLOSS stations taken from the
# table in
# \url{http://www.gloss-sealevel.org/station_handbook/stations/#.Vqtc6kL4phg}
# 
# Some data sets, such as NINO3.4 and NAOI come with a 'frozen' version in the
# package, but there are also functions that read the most recent version of
# these indeces from the Internet. GSL reads the global mean sea level.
# 
# 
# @aliases station.meta NACD NARP Oslo Svalbard t2m.NORDKLIM precip.NORDKLIM
# geoborders t2m.NCEP sst.NCEP precip.ERAINT slp.NCEP slp.atlantic.ERA5
# t2m.NorESM.M t2m.DNMI sst.DNMI slp.DNMI eof.t2m.NCEP eof.sst.NCEP
# eof.slp.NCEP eof.precip.ERAINT eof.t2m.NorESM.M eof.t2m.DNMI eof.sst.DNMI
# eof.slp.DNMI NAOI sunspots NINO3.4 SOI GSL GSL.nasa AMO bjornholt vardo
# ferder dse.ferder HadCRUT4 NASAgiss dse.Oslo glossstations
# IPCC.AR5.Table.9.A.1 cmipgcmresolution global.t2m.gcm QBO CET CO2
# read.hurdat2 read.best.track storms
# @param lon longitude range c(lin.min,lon.max)
# @param lat latitude range
# @param anomaly TRUE: return anomaly
# @param url source of data
# @param plot TRUE:plot
# @return Numeric vectors/matrices with a set of attributes describing the
# data.
# @author R.E. Benestad
# @seealso \code{\link{aggregate.area}} \code{\link{as.4seasons}},
# \code{\link{annual}}
# @keywords datasets
# @examples
# 
# data(Oslo)
# year <- as.numeric( format(index(Oslo), '%Y') ) 
# plot(aggregate(Oslo, by=year,FUN='mean', na.rm = FALSE), new=FALSE)
# 
# data(etopo5)
# z <- subset(etopo5,is=list(lon=c(-10,30),lat=c(40,60)))
# map(z, new=FALSE)
# 
# 





# InfoGraphics.
# 
# Wheel
# 
# Risk
# 
# visprob - boxes with forseen outcomes - area proportional to probability
# 
# conf - confidence intervals and uncertainty - clouds...
# 
# vis
# 
# diagram
# 
# cumugram
# 
# graph
# 
# \code{balls} draws 3D ball symbols.
# 
# 
# @aliases wheel wheel.station wheel.spell diagram diagram.dsensemble
# diagram.station seasevol seasevol.station vis vis.station diagram cumugram
# climvar colscal visprob balls graph graph.dsensemble graph.zoo graph.list
# @param x a data object
# @param new TRUE, invoke \code{dev.new()}
# @param col Colour palette
# @param col.obs Colour for the observations
# @param n number of colour scales
# @param main main title
# @param log.precip Use logarithmig scale for precipitation
# @param plot FALSE, just return the results
# @param verbose If TRUE, print out diagnostics
# @param img image from jpeg
# @param start starting date of the year
# @param ensmean Ensemble mean
# @param pch Plot symbol style
# @param it index time
# @param xlim see \code{\link{plot}}
# @param ylim see \code{\link{plot}}
# @param xlab see \code{\link{plot}}
# @param ylab see \code{\link{plot}}
# @param lwd see \code{\link{plot}}
# @param add TRUE: update existing plot
# @return A field object
# @author R.E. Benestad and A.  Mezghanil
# @seealso \code{\link{map}}, \code{\link{plot.station}},
# \code{\link{hist.spell}}
# @keywords utilities
# @examples
# 
# data(bjornholt)
# wheel(bjornholt)
# 
# z <- spell(bjornholt,threshold=1)
# wheel(z)
# 
# 






#' Manual for esd
#' 
#' Help, assistance and manuals for \code{esd} and empirical-statistical
#' downsclaing (ESD).
#' 
#' 
#' @name manual
#' @aliases manual ABC4ESD esd.tips esd.issues downscaling.about
#' @docType data
#' @keywords datasets
#' @examples
#' 
#' manual()
#' 
NULL





#' Calibrated models.
#' 
#' \code{mu.eq.f.tx} contains a regression model \eqn{mu=f(tx)}{\mu==f(t_x)},
#' relating the wet-day mean (mu) to the mean daily maximum temperature -
#' however, the input used in this model is the saturation vapour pressure
#' according to the Clausius-Clapeyron equation, and 'f()' here is 'beta *
#' C.C.eq()'.
#' 
#' 
#' @aliases mu.eq.f.tx
#' @return lm() object.
#' @author R.E. Benestad
#' @keywords datasets
#' @examples
#' 
#' ## Retrieve the model
#' data(mu.eq.f.tx)
#' ## Sample data - temperature from Ferder Lighthouse
#' data(ferder)
#' pre <- data.frame(x=mean(C.C.eq(ferder),na.rm=TRUE))
#' ## Predict the wet-day mean based on mean temperatures
#' predict(mu.eq.f.tx,newdata=pre)
#' 
NULL





#' %% ~~ data name/kind ... ~~ Oslo monthly mean temperature time series
#' 
#' %% ~~ A concise (1-5 lines) description of the dataset. ~~ Oslo temperature
#' monthy record from 1837 up to now.
#' 
#' %% ~~ If necessary, more details than the __description__ above ~~ Oslo
#' surface temperature recorded on a monthly basis from 1837 up to 2012. It
#' corresponds to a blended recnostruction (1837-1936) and instrumental data
#' (1937-2012). An homogenisation procedure has been carried out on the data by
#' Øyvind Nordli at MET Norway.
#' 
#' @name Oslo
#' @docType data
#' @format The format is: 'zoo' series from 1837-01-01 to 2014-02-01 \cr Data:
#' atomic [1:2126] NaN NaN NaN 3.1 8.9 13.5 16 15.7 11.2 7.6 ... \cr - attr(*,
#' "location")= chr "Oslo" \cr - attr(*, "variable")= chr "T2m" \cr - attr(*,
#' "unit")= chr "deg C" \cr - attr(*, "longitude")= num 10.7 \cr - attr(*,
#' "latitude")= num 59.9 \cr - attr(*, "altitude")= num 94 \cr - attr(*,
#' "country")= chr "Norway" \cr - attr(*, "longname")= chr "temperature at 2m"
#' \cr - attr(*, "station_id")= num 18700 \cr - attr(*, "quality")= chr
#' "homogenised" \cr - attr(*, "calendar")= chr "gregorian" \cr - attr(*,
#' "source")= chr "Dr. Nordli, 2013, met.no" \cr - attr(*, "URL")= logi NA \cr
#' - attr(*, "type")= chr "observation"\cr - attr(*, "aspect")= chr "original"
#' \cr - attr(*, "reference")= chr "Nordli et al. (in progress). 'The Oslo
#' Temperature series 1837-2012: Homogeneity testing and Climate Analysis'" \cr
#' - attr(*, "info")= logi NA \cr - attr(*, "method")= chr "Blended
#' recnostruction (1877-1936) and instrumental data (1937-)" \cr - attr(*,
#' "history")=List of 3 \cr ..$ call :length 19 as.station.data.frame(oslo, loc
#' = "Oslo", param = "T2m", unit = "deg C", lon = 10.7, lat = 59.9, alt = 94,
#' cntr = "Norway", longname = "temperature at 2m", ... \cr .. ..- attr(*,
#' "srcref")=Class 'srcref' atomic [1:8] 85 1 90 104 1 104 85 90 \cr .. .. ..
#' ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment:
#' 0x23b1980> \cr ..$ timestamp : chr "Tue Dec 10 14:53:51 2013" \cr ..$
#' sessioninfo:List of 3 \cr .. ..$ R.version : chr "R version 3.0.2
#' (2013-09-25)" \cr .. ..$ esd.version: chr "esd_0.2-1" \cr .. ..$ platform :
#' chr "x86_64-pc-linux-gnu (64-bit)" \cr Index: Date[1:2126], format:
#' "1837-01-01" "1837-02-01" "1837-03-01" "1837-04-01" ... \cr
#' @references %% ~~ possibly secondary sources and usages ~~ Nordli et al. (in
#' progress). 'The Oslo Temperature series 1837-2012: Homogeneity testing and
#' Climate Analysis'
#' @source %% ~~ reference to a publication or URL from which the data were
#' obtained ~~ MET Norway
#' @keywords datasets
#' @examples
#' 
#' data(Oslo)
#' ## maybe str(Oslo) ; plot(Oslo) ...
#' 
NULL





#' %% ~~function to do ... ~~ Simple and handy functions.
#' 
#' \code{lag.station} and \code{lag.station} are wrap-around functions for
#' \code{\link[zoo]{lag.zoo}} that maintains all the attributes.
#' \code{attrcp(x,y)} passes on attributes from \code{x} to \code{y} and
#' returns the \code{y} with new attributes. \code{ensemblemean} returns the
#' ensemble mean for dsensemble objects. The argument \code{FUN} can also be
#' used to estiamte other statistics, e.g \code{FUN='q95'} where
#' \code{q95=function(x) apply(x,1,quantile,probs=0.95)} %% ~~ A concise (1-5
#' lines) description of what the function does. ~~ \code{TGW} uses
#' triangulation of pressure measurements to estimate geostrophic wind based on
#' Alexandersson et al. (1998), Glob. Atm. and Oce. Sys. \code{stand} gives a
#' standardised time series.
#' 
#' \code{arec} compares the number of record-breaking events to the number of
#' events for a stime series of iid data (\code{sum(1/1:n)}) \code{strstrip}
#' strips off spaces before and after text in a string.
#' 
#' 'as.decimal' converts between degree-minute-second into decimal value.
#' 
#' 'cv' computes the coefficient of variation.
#' 
#' 'nv' count the number of valid data points.
#' 
#' 'q5','q95' and 'q995' are shortcuts to the 5\%, 95\%, and 99.5\%
#' percentiles.
#' 
#' 'trend.coef' and 'trend.pval' return the coefficient and the p-value of the
#' linear trend.
#' 
#' 'exit' is a handy function for exiting the R session without saving.
#' 
#' 'figlab' is a handy function for labelling figures (e.g. 'Figure 1')
#' 
#' 'ndig' estimates the number of digits for round(x,ndig), e.g. in scales for
#' plotting.
#' 
#' @aliases as.decimal nv cv q5 q95 q995 lag.station lag.field filt
#' filt.default exit figlab ndig attrcp ensemblemean propchange stand rmse RMSE
#' firstyear lastyear eofvar test.ds.field test.num.predictors arec
#' arec.default arec.station lastrains lastelementrecord strstrip bin
#' @param x A data.frame or a coredata zoo object.
#' @param na.rm If TRUE, remove NA's from data
#' @param type 'ma' for moving average (box-car), 'gauss' for Gaussian, 'binom'
#' for binomial filter,' parzen' for Parzen filter, 'hanning' for Hanning
#' filter, or 'welch' for Welch filter.
#' @param lowpass True for smoothing, otherwise the highpass results is
#' returned
#' @param triangle a group of three stations with sea-level pressure', e.g.
#' from ECA\&D.
#' @param nbins number of bins/categories
#' @return \item{as.decimal }{Decimal value} \item{trend.coef }{Linear trend
#' per decade}
#' @author A. Mezghani
#' @keywords rtools ~kwd2
#' @examples
#' 
#' ## Monthly mean temperature at 'Oslo - Blindern' station 
#' data(Oslo)
#' ## Compute the linear trend and the p-value on annual aggregated values 
#' tr <- trend.coef(coredata(annual(Oslo)))
#' pval <- trend.pval(coredata(annual(Oslo)))
#' \dontrun{
#' pp <- station(param='slp',cntr='Denmark',src='ecad')
#' wind <- TGW(subset(pp,is=c(1,3,10))
#' plot(wind)
#' ws <- sqrt(wind[,1]^2 + wind[,2]^2)
#' plot(ws)
#' hist(ws)
#' 
#' ## Estimate wind for a larger group of stations
#' wind <- geostrophicwind(pp,nmax=10)
#' u <- subset(wind,param='u')
#' v <- subset(wind,param='u')
#' ws <- sqrt(u^2+v^2)
#' ws <- attrcp(v,ws)
#' class(ws) <- class(v)
#' attr(ws,'variable')='windspeed'
#' attr(ws,'longname')='geostrophic wind speed'
#' map(ws,FUN='quantile',probs=0.98)
#' 
#' ## Test firstyears on HadCRUT4
#' if (!file.exists('~/Downloads/HadCRUT.4.6.0.0.median.nc')) {
#'   print('Download HadCRUT4')
#'   download.file('https://crudata.uea.ac.uk/cru/data/temperature/HadCRUT.4.6.0.0.median.nc',
#'                 dest='~/Downloads/HadCRUT.4.6.0.0.median.nc') 
#' }
#' 
#' Obs <- annual(retrieve('~/Downloads/HadCRUT.4.6.0.0.median.nc',param='temperature_anomaly'))
#' lons <- rep(lon(Obs),length(lat(Obs)))
#' lats <- sort(rep(lat(Obs),length(lon(Obs))))
#' fy <- firstyear(Obs)
#' map(subset(Obs,it=1))
#' points(lons[fy==1850],lats[fy==1850])
#' map(Obs,FUN='firstyear')
#' }
#' 
NULL





#' Weather station metadata
#' 
#' Meta datasets based on different data sources or datasets included in esd
#' package. Mainly weather stations' coordinates and other meta data.
#' 
#' 
#' @name station.meta
#' @docType data
#' @format The format is: \cr List of 12 \cr $ station_id: chr [1:284239]
#' "6447" "6193" "21100" "25140" ... \cr $ location : chr [1:284239] "UCCLE"
#' "HAMMERODDE_FYR" "VESTERVIG" "NORDBY" ... \cr $ country : chr [1:284239]
#' "BELGIA" "DENMARK" "DENMARK" "DENMARK" ... \cr $ longitude : num [1:284239]
#' 4.35 14.78 8.32 8.4 10.6 ... \cr $ latitude : num [1:284239] 50.8 55.3 56.8
#' 55.4 55.9 ... \cr $ altitude : num [1:284239] 100 11 18 5 11 9 4 51 85 105
#' ... \cr $ element : chr [1:284239] "101" "101" "101" "101" ... \cr $ start :
#' chr [1:284239] "1833" "1853" "1873" "1872" ... \cr $ end : chr [1:284239] NA
#' NA NA NA ... \cr $ source : chr [1:284239] "NACD" "NACD" "NACD" "NACD" ...
#' \cr $ wmo : num [1:284239] 6447 6193 -999 -999 -999 ... \cr $ quality : int
#' [1:284239] 2 2 2 2 2 1 1 2 5 5 ... \cr - attr(*, "history")= chr [1:7] \cr
#' "meta2esd.R - data taken from the clim.pact package and consolidated for
#' NACD and NARP" "nordklim.meta.rda" "ecad.meta.rda" "ghcnd.meta.rda" ... \cr
#' - attr(*, "date")=function () \cr - attr(*, "call")= language meta2esd(save
#' = TRUE) \cr - attr(*, "author")= chr "R.E. Benestad & A. Mezghani" \cr -
#' attr(*, "URLs")= chr [1:6] "www.dmi.dk/dmi/sr96-1.pdf" \cr
#' "http://www.norden.org/en/publications/publikationer/2005-450" \cr
#' "http://www.smhi.se/hfa_coord/nordklim/" "http://eca.knmi.nl/" ... \cr
#' @seealso \code{\link{select.station}},\code{\link{station}}
#' @keywords datasets
#' @examples
#' 
#' data(station.meta)
#' str(station.meta)
#' map(station.meta)
#' 
NULL





#' Various formulas, equations and transforms.
#' 
#' Computes different formulas.
#' 
#' \code{C.C.eq}: Clapeyron-Clausius equation (saturation evaporation pressure)
#' where \code{x} is a data object holding the temperature.
#' 
#' \code{precip.vul}: and index for the vulerability to precipitation defined
#' as wetmean(x)/wetfreq(x). High when the mean intensity is high and/or the
#' frequency is low (it rains seldom, but when it rains, it really pours down).
#' 
#' \code{t2m.vul}: and index for the vulerability to temperature defined as the
#' mean spell length for heat waves with temperatures exceeding 30C (default).
#' 
#' \code{precip.rv}: a rough estimate of the return value for precipitation
#' under the assumption that it is exponentially distributed. Gives apprximate
#' answers for low return levels (less than 20 years). Advantage, can be
#' predicted given wet-day mean and frequency.
#' 
#' \code{nv}: number of valid data points.
#' 
#' \code{precip.Pr}: rough estimate of the probability of more than x0 of rain
#' based on an exponential distribution.
#' 
#' \code{t2m.Pr}: rough estimate of the probability of more than x0 of rain
#' based on a normal distribution.
#' 
#' \code{NE}: predicts the number of events given the probability Pr.
#' 
#' 
#' @aliases C.C.eq precip.vul t2m.vul precip.rv nv precip.Pr t2m.Pr NE
#' @param x a data object
#' @param p a probability
#' @param x0 a threshold value
#' @param tau time scale (years)
#' @param is which of the spell results [1,2]
#' @param na.rm See \code{\link{mean}}.
#' @return The right hand side of the equation
#' @author R. Benestad, MET Norway
#' @keywords parameter,element
#' @examples
#' 
#' t2m <- t2m.DNMI(lon=c(-70,-10),lat=c(20,60))
#' es <- C.C.eq(t2m)
#' map(es)
#' 
NULL



