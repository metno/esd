# not exported
ar1 <- function(x,...) acf(x,plot=FALSE,na.action = na.pass)$acf[2]

# not exported
ltp <- function(x,type='exponential',...) {
  # Rybski et al. (2006), i:10.1029/2005GL025591
  ar <- acf(x,plot=FALSE)$acf
  n <- min(11,(1:length(ar))[ar <= 0])-1
  gamma <- log(ar[1:n])
  l <- 1:length(gamma)
#  qfit <- lm(gamma ~ l + I(l^2))
  qfit <- lm(gamma ~ l)
  gamma <- qfit$coefficients[2]
  if (type=='exponential') y <- gamma else
  if (type=='hurst') y <- 1 - gamma/2
  return(y)
}




#' Downscale ensemble runs
#' 
#' Downscales an ensemble of climate model runs, e.g. CMIP5, taking the results
#' to be seasonal climate statistics. For temperature, the result hold the
#' seasonal mean and standard deviation, whereas for precipitation, the results
#' hold the wet-day mean, the wet-day frequency, and the wet/dry-spell
#' statistics. The call assumes that netCDF files containing the climate model
#' ensemble runs are stores in a file structure, linked to the path argument
#' and the rcp argument.
#' 
#' These methods are based on \code{\link{DS}}, and \code{DSensemble} is
#' designed to make a number of checks and evaluations in addition to
#' performing the DS on an ensemble of models. It is based on a similar
#' philosophy as the old R-package '\code{clim.pact}', but there is a new
#' default way of handling the predictors. In order to attempt to ensure a
#' degree of consistency between the downscaled results and those of the GCMs,
#' a fist covariate is introduced before the principal components (PCs)
#' describing the \code{\link{EOF}s}.
#' 
#' \code{DSensemble.pca} is used to downscale a predictor represented in terms
#' of PCA, and can reduce the computation time significantly. See Benestad et
#' al. (2015) \url{http://dx.doi.org/10.3402/tellusa.v67.28326}.
#' 
#' The argument \code{non.stationarity.check} is used to conduct an additional
#' test, taking the GCM results as 'pseudo-reality' where the predictand is
#' replaced by GCM results interpolated to the same location as the provided
#' predictand. The time series with interpolated values are then used as
#' predictor in calibrating the model, and used to predict future values. This
#' set of prediction is then compared with the interpolated value itself to see
#' if the dependency between the large and small scales in the model world is
#' non-stationary.
#' 
#' Other chekch include cross-validation (\code{\link{crossval}}) and
#' diagnostics comparing the sample of ensemble results with the observations:
#' number of observations outside the predicted 90-percent conf. int and
#' comparing trends for the past.
#' 
#' The 'bias correction' is described in Imbert and Benestad (2005),
#' \emph{Theor. Appl. Clim.} \url{http://dx.doi.org/10.1007/s00704-005-0133-4}.
#' 
#' 
#' @aliases DSensemble DSensemble.default DSensemble.station DSensemble.t2m
#' DSensemble.precip DSensemble.annual DSensemble.season DSensemble.field
#' DSensemble.mu.worstcase DSensemble.pca DSensemble.eof
#'
#' @importFrom graphics par grid segments text axis abline
#' @importFrom stats ks.test pnorm acf sd na.pass lm rnorm window start end
#' cor.test
#'
#' @param y A station object.
#' @param plot Plot intermediate results if TRUE.
#' @param path The path where the GCM results are stored.
#' @param rcp Which (RCP) scenario
#' @param biascorrect TRUE, apply a bias adjustment using \code{\link{biasfix}}
#' @param predictor The predictor, a field or EOF object
#' @param non.stationarity.check If TRUE perform stationarity test - work in
#' progress
#' @param ip Which EOFs to include in the step-wise multiple regression.
#' @param rmtrend TRUE: detrend before calibrating the regression model.
#' @param lon Longitude range for predictor
#' @param lat Latitude range for predictor
#' @param rel.cord TRUE: use the range relative to predictand; FALSE use
#' absolute range
#' @param it Used to extract months or a time period. See \code{\link{subset}}.
#' @param select GCMs to select, e.g .subsample the ensemble (1:3 selects the
#' three first GCMs)
#' @param FUN Function for aggregating the predictand (daily), e.g. 'mean',
#' 'wetmean'
#' @param threshold Used together with FUN for some functions ('wetmean').
#' @param nmin Minimum number of day used in \code{\link{annual}} used for
#' aggregating the predictand/predictor
#' @param FUNX Function for transforming the predictor, e.g.
#' '\code{\link{C.C.eq}}' to estimate the saturation water vapout
#' @param type Type of netCDF used in \code{\link{retrieve}} for reading GCM
#' data.
#' @param pattern File name pattern for GCM data.
#' @param verbose TRUE for checking and debugging the functions.
#' @param file.ds Name of file saving the results.
#' @param path.ds Path of file saving the results.
#' @param xfuns Names of functions which do not work in \code{annual(x,FUN=f)}.
#' These functions are used using the following code
#' \code{annual(f(x),FUN="mean")}
#' @param mask TRUE mask out land
#' @param ds.1900.2099 Default, only downscale for the period 1900-2099
#' @param \dots additional arguments
#'
#' @return A 'dsensembele' object - a list object holding DS-results.
#' 
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' # Import historical temperature data from Oslo
#' data(Oslo)
#' 
#' ## Download NorESM1-M from 'climexp.knmi.nl' in default directory
#' ## (home directory for linux/mac users)
#' url <-"http://climexp.knmi.nl/CMIP5/monthly/tas"
#' ## Download NorESM1-ME for the emission scenario RCP4.5
#' noresm <- "tas_Amon_NorESM1-M_rcp45_000.nc"
#' if (!file.exists(noresm)) {
#'   download.file(url=file.path(url,noresm), destfile=noresm,
#'                 method="auto", quiet=FALSE, mode="w", cacheOK=TRUE)
#' }
#' 
#' ## Download FIO-ESM for the emission scenario RCP4.5
#' fioesm <- "tas_Amon_FIO-ESM_rcp45_000.nc"
#' if (!file.exists(fioesm)) {
#'   download.file(url=file.path(url,fioesm), destfile=fioesm,
#'                 method="auto", quiet=FALSE, mode="w", cacheOK=TRUE)
#' }
#' 
#' ## Downscale the predictor (ERA-interim reanalysis 2m temperature)
#' predictor <- "air.2m.mon.mean.nc"
#' if (!file.exists(predictor)) {
#'   url <-"http://climexp.knmi.nl/ERA-interim/erai_t2m.nc"
#'   download.file(url=url, destfile=predictor,
#'                 method="auto", quiet=FALSE, mode="w", cacheOK=TRUE)
#' }
#' 
#' # Downscale the temperature in Oslo
#' rcp4.5 <- DSensemble.t2m(Oslo, path='~', rcp='', pattern="tas_Amon_",
#'                          biascorrect=TRUE, predictor = predictor,
#'                          plot=TRUE, verbose=TRUE)
#' 
#' ## Evaluation: 
#' ## (1) combare the past trend with downscaled trends for same
#' ## interval by ranking and by fitting a Gaussian to the model ensemble;
#' ## (2) estimate the probabilty for the counts outside the 90
#' ## percent confidence interval according to a binomial distribution.
#' diagnose(rcp4.5, plot = TRUE, type = "target")
#' }
#' 
#' @export DSensemble
DSensemble<-function(y,...) UseMethod("DSensemble")

#' @export DSensemble.default
DSensemble.default <- function(y,...,path='CMIP5.monthly/',rcp='rcp45') {
   ## 
  stopifnot(!missing(y),inherits(y,"station"),
            file.exists(paste(file.path(path,rcp,fsep = .Platform$file.sep))))
  
  if (is.null(attr(y,'aspect'))) attr(y,'aspect') <- "original"
  
  if (inherits(y,'pca')) {
    z <- DSensemble.pca(y,path=path,rcp=rcp,...) 
  } else if (inherits(y,c('eof','field'))) {
    z <- DSensemble.eof(y,path=path,rcp=rcp,...)
  } else if (inherits(y,'annual')) {
    z <- DSensemble.annual(y,path=path,rcp=rcp,threshold=1,...) 
  } else if (inherits(y,'season')) {
    z <- DSensemble.season(y,path=path,rcp=rcp,...) 
  } else if (is.T(y)) {
    z <- DSensemble.t2m(y,path=path,rcp=rcp,...) 
  } else if (is.precip(y)) {
    z <- DSensemble.precip(y,path=path,rcp=rcp,threshold=1,...) 
  } else {
    z <- NULL 
  }
  return(z)
}

#' @export DSensemble.t2m
DSensemble.t2m <- function(y,...,plot=TRUE,path="CMIP5.monthly/",predictor="ERA40_t2m_mon.nc",
                           rcp="rcp45",biascorrect=FALSE,non.stationarity.check=FALSE,
                           ip=1:6,lon=c(-20,20),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                           select=NULL,FUN="mean",FUNX="mean",xfuns='C.C.eq',
                           pattern="tas_Amon_ens_",path.ds=NULL,file.ds="DSensemble.rda",
                           nmin=NULL,verbose=FALSE,ds.1900.2099=TRUE) {

  if (!inherits(y,'day')) warning('station is not daily data')
  if (verbose) print("predictand")
  #if ( (deparse(substitute(FUN))=='sd') | (deparse(substitute(FUN))=='ar1') )
  if(verbose) print("DSensemble.t2m")
  if ((FUN=='sd') | (FUN =='ar1')) {
    y <- anomaly(y)
    attr(y,'aspect') <- 'original'
  }

  if(verbose) print("aggregate time series")
  if (is.null(nmin)) nmin4 <- nmin else nmin4 <- nmin*4
  ya <- annual(y,FUN=FUN,nmin=nmin4)
  y <- as.4seasons(y,FUN=FUN,nmin=nmin)

  if(verbose) print("Set lon/lat predictor range")
  if ( !is.na(attr(y,'longitude'))[1] & (any(lon>0) & any(lon<0)) & rel.cord)
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if ( !is.na(attr(y,'latitude'))[1] & (any(lat>0) & any(lat<0)) & rel.cord)
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )
  
  if (sum(!is.finite(lon))>0) 
    warning(paste('Bad longitude range provided: ',paste(lon,collapse='-')))
  if (sum(!is.finite(lat))>0) 
    warning(paste('Bad latitude range provided: ',paste(lat,collapse='-')))
  
  # The units
  if(verbose) print("Arrange units")
  if ( (unit(y)[1] == "deg C") | (unit(y)[1] == "degree Celsius") |
       (unit(y)[1] == "degC") | (unit(y)[1] == "Celsius") )
        unit <- expression(degree*C) else
        unit <- attr(y,'unit')
  if (verbose) print(paste('Units:',unit))

  if (plot) {
    if(verbose) print("Plot station data (predictand)")
    ylim <- c(0,0)
    ylim <- switch(FUN,'mean'=c(-2,8),'sd'=c(-0.5,1),'ar1'=c(-0.5,0.7)) # assuming y is the temperature?
    if (verbose) print(paste('set ylim based on "',FUN,'" -> c(',ylim[1],', ',ylim[2],')',sep=''))
    par(bty="n")
    plot.zoo(ya,type="b",pch=19,main=attr(y,'location'),
             xlab="year",ylab=unit,
             sub=paste('Station: ',attr(y,'station_id'),'; coordinates: ',
             round(attr(y,'longitude'),4),'E/',
             round(attr(y,'latitude'),4),'N; ',
             attr(y,'altitude'),'m.a.s.l',sep=''),
             ylim=ylim + range(coredata(ya),na.rm=TRUE),xlim=c(1900,2100))
    grid()
  }
  if(verbose) print("Retrieve predictor data")
  if (is.character(predictor))
    t2m <- retrieve(file=predictor,lon=lon,lat=lat,
                    verbose=verbose) else
  if (inherits(predictor,'field'))
    t2m <- subset(predictor,is=list(lon=lon,lat=lat))

  if(verbose) print("Aggregate seasonal values")
  T2M <- as.4seasons(t2m,FUN=FUNX,nmin=nmin)
  DJF <- subset(as.4seasons(t2m,FUN=FUNX,nmin=nmin),it='djf')
  MAM <- subset(as.4seasons(t2m,FUN=FUNX,nmin=nmin),it='mam')
  JJA <- subset(as.4seasons(t2m,FUN=FUNX,nmin=nmin),it='jja')
  SON <- subset(as.4seasons(t2m,FUN=FUNX,nmin=nmin),it='son')

  ok1 <- is.finite(rowSums(DJF))
  DJF <- subset(DJF,it=range(year(DJF)[ok1]))
  ok2 <- is.finite(rowSums(MAM))
  MAM <- subset(MAM,it=range(year(MAM)[ok2]))
  ok3 <- is.finite(rowSums(JJA))
  JJA <- subset(JJA,it=range(year(JJA)[ok3]))
  ok4 <- is.finite(rowSums(SON))
  SON <- subset(SON,it=range(year(SON)[ok4]))

  rm("t2m"); gc(reset=TRUE)
  
  # Ensemble GCMs
  if(verbose) print("Retrieve & arrange GCMs")
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.names=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) select <- 1:N else
      select <- select[select<=N]; N <- length(select)
  if (verbose) {print('GCMs:'); print(path); print(ncfiles[select])}
  
  if(verbose) print("Set up results matrix & table of diagnostics")
  years <- sort(rep(1900:2100,4))
  months <- rep(c(1,4,7,10),length(1900:2100))
  m <- length(years)
  X <- matrix(rep(NA,N*m),N,m)
  gcmnm <- rep("",N)
  scorestats <- matrix(rep(NA,N*9),N,9)

  colnames(scorestats) <- c("1-r.xval","mean.diff","1-sd.ratio",
                            "1-autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration',
                            '1-R2')

  t <- as.Date(paste(years,months,'01',sep='-'))

  if(verbose) print("Quick test")  
  #load('control.ds.1.rda'); T2M.ctl -> T2M
  flog <- file("DSensemble.t2m-log.txt","at")
  
  if (verbose) print("loop...")
  for (i in 1:N) {
    if (verbose) print(ncfiles[select[i]])
    gcm <- retrieve(file = ncfiles[select[i]],
                          lon=range(lon(T2M))+c(-2,2),
                          lat=range(lat(T2M))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099)) else
                      gcm <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
    if (length(index(gcm))<=1) print(paste('Problem selecting GCM results in period',
                                           min(year(y),na.rm=TRUE),'2099'))
    #gcmnm[i] <- attr(gcm,'model_id'))
    #gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
    # KMP: 10.03.2017 - pass on additional information about GCM runs (gcm + rip - realization, initialization, physics version)
    ## REB: 2020-12-23: These lines cause some problems sometimes:
    ## Error in names(attributes(gcm))[grep("realization", names(attributes(gcm)))][[1]] : 
    ##  subscript out of bounds
    ## Included checks to avoid unnecessary stops 
    nmattsgcm <- names(attributes(gcm))
    if (length(grep("realization",nmattsgcm)) > 0) nm.r <- attributes(gcm)[grep("realization",nmattsgcm)][[1]] else nm.r <- '_'
    if (length(grep("initialization",nmattsgcm)) > 0) nm.i <- attributes(gcm)[grep("initialization",nmattsgcm)][[1]] else nm.i <- '_'
    if (length(grep("physics",nmattsgcm)) > 0) nm.p <- attributes(gcm)[grep("physics",nmattsgcm)][[1]] else nm.p <- '_'
    rip <- paste0("r",attr(gcm,nm.r),"i",attr(gcm,nm.i),"p",attr(gcm,nm.p))
    gcmnm.i <- paste0(attr(gcm,'model_id'),".",rip)
    gcmnm[i] <- gcmnm.i
    #gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'parent_experiment_rip'),sep="-")
    # REB: 30.04.2014 - new lines...
    DJFGCM <- subset(as.4seasons(gcm,FUN=FUNX,nmin=nmin),it='djf')
    MAMGCM <- subset(as.4seasons(gcm,FUN=FUNX,nmin=nmin),it='mam')
    JJAGCM <- subset(as.4seasons(gcm,FUN=FUNX,nmin=nmin),it='jja')
    SONGCM <- subset(as.4seasons(gcm,FUN=FUNX,nmin=nmin),it='son')
    rm("gcm"); gc(reset=TRUE)
    X.DJF <- combine(DJF,DJFGCM)
    X.MAM <- combine(MAM,MAMGCM)
    X.JJA <- combine(JJA,JJAGCM)
    X.SON <- combine(SON,SONGCM)
    if (verbose) print("- - - > EOFs")
    Z1 <- try(EOF(X.DJF))
    if (inherits(Z1,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(Z1[[1]],con=flog)
    }
    Z2 <- try(EOF(X.MAM))
    if (inherits(Z2,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(Z2[[1]],con=flog)
    }
    Z3 <- try(EOF(X.JJA))
    if (inherits(Z3,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(Z3[[1]],con=flog)
    }
    Z4 <- try(EOF(X.SON))
    if (inherits(Z4,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(Z4[[1]],con=flog)
    }
    # The test lines are included to assess for non-stationarity
    if (non.stationarity.check) {
      ## KMP 2018-11-12: Not sure if this test works.
      ## It's only done for DJF here but could be expanded to other seasons
      ## if it is a useful diagnostic.
      testGCM <- subset(DJFGCM,it=range(year(DJF))) # REB 29.04.2014
      testy <- as.station(regrid(testGCM,is=DJF))  # REB 29.04.2014
      attr(testGCM,'source') <- 'testGCM'        # REB 29.04.2014
      testZ <- combine(testGCM,DJFGCM)              # REB 29.04.2014
      rm("testGCM"); gc(reset=TRUE)
    }
    ##
    # REB: 30.04.2014 - new lines...
    if (verbose) print("- - - > DS (seasonal)")
    if (verbose) print(class(attr(Z1,'appendix.1')))
    if (biascorrect) try(Z1 <- biasfix(Z1))
    ds1 <- try(DS(subset(y,it='djf'),Z1,ip=ip))
    if (inherits(ds1,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds1[[1]],con=flog)
    }
    if (biascorrect) try(Z2 <- biasfix(Z2))
    ds2 <- try(DS(subset(y,it='mam'),Z2,ip=ip))
    if (inherits(ds2,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds2[[1]],con=flog)
    }
    if (biascorrect) try(Z3 <- biasfix(Z3))
    ds3 <- try(DS(subset(y,it='jja'),Z3,ip=ip))
    if (inherits(ds3,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds3[[1]],con=flog)
    }
    if (biascorrect) try(Z4 <- biasfix(Z4))
    ds4 <- try(DS(subset(y,it='son'),Z4,ip=ip))
    if (inherits(ds4,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds4[[1]],con=flog)
    }
    if (verbose) print("Combine the 4 seasons")
    ds <- try(combine(list(ds1,ds2,ds3,ds4)))
    rm("Z1","Z2","Z3","Z4")
    gc(reset=TRUE) ##rm("ds1","ds2","ds3","ds4")
    if (inherits(ds,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds[[1]],con=flog)
    } else {
      x.val <- c(crossval(ds1),crossval(ds2),crossval(ds3),crossval(ds4))
      attr(ds,'evaluation') <- x.val

      if (verbose) print("post-processing")
      z <- attr(ds,'appendix.1')
    #save(file='inside.dsens.1.rda',ds,y,Z1)

      if (non.stationarity.check) {
        testds <- DS(testy,testZ,biascorrect=biascorrect,ip=ip)   # REB 29.04.2014
        testz <- attr(testds,'appendix.1')                      # REB 29.04.2014
        difference.z <- testy - testz                           # REB 29.04.2014
      }
      i1 <- is.element(paste(years,months,sep='-'),
                       paste(year(z),month(z),sep='-'))
      i2 <- is.element(paste(year(z),month(z),sep='-'),
                       paste(years,months,sep='-'))
    #
      X[i,i1] <- z[i2]

    # Diagnose the residual: ACF, pdf, trend. These will together with the
    # cross-validation and the common EOF diagnostics provide a set of
    # quality indicators.
      cal <- coredata(attr(ds,"original_data"))
      fit <- coredata(attr(ds,"fitted_values"))
      res <- as.residual(ds)
      res.trend <- 10*diff(range(trend(res)))/diff(range(year(res)))
      ks <- round(ks.test(coredata(res),pnorm)$p.value,4)
      ar <- as.numeric(acf(trend(cal-fit,result="residual"),plot=FALSE)[[1]][2]) ##ar <- ar1(coredata(res))
      if (verbose) print(paste("Examine residuals: trend=",
                               round(res.trend,3),'D/decade; K.S. p-val',
                               round(ks,2),'; AR(1)=',round(ar,2)))

    # Evaluation: here are lots of different aspects...
    # Get the diagnostics: this is based on the analysis of common EOFs...

      xval <- attr(ds,'evaluation')
      r.xval <- cor(xval[,1],xval[,2])
      if (verbose) print(paste("x-validation r=",r.xval))
    
      dsa <- annual(ds)                     # annual mean value

      xy <- merge.zoo(annual(z),ya)
      ds.ratio <- sd(xy[,1],na.rm=TRUE)/sd(xy[,2],na.rm=TRUE)
      if (verbose) print(paste("sd ratio=",ds.ratio))

    #print(names(attributes(ds)))
      if (biascorrect) {
        diag <- attr(ds,'diagnose')
        if ( (verbose) & !is.null(diag)) str(diag)
      } else diag <- NULL
    
    # diagnose for ds-objects
      ##
      if (verbose) print('...')
      if (is.null(diag)) {
        ##diag <- diagnose(z,plot=FALSE)
        scorestats[i,] <- c(1-r.xval,NA,NA,NA,res.trend,ks,1-ar,1-ds.ratio,
                            1-round(var(xval[,2])/var(xval[,1]),2))
        mdiff <- (mean(subset(ya,it=range(year(dsa))),na.rm=TRUE)-
                  mean(subset(dsa,it=range(year(ya))),na.rm=TRUE))/sd(ya,na.rm=TRUE)
        srati <- sd(subset(dsa,it=range(year(ya))),na.rm=TRUE)/
                 sd(subset(ya,it=range(year(dsa))),na.rm=TRUE)
        arati <- ar1(dsa)/ar1(ya)
      } else {

    # Extract the mean score for leading EOF from the 4 seasons:
        mdiff <- mean(c(diag$s.1$mean.diff[1]/diag$s.1$sd0[1],
                        diag$s.2$mean.diff[1]/diag$s.2$sd0[1],
                        diag$s.3$mean.diff[1]/diag$s.3$sd0[1],
                        diag$s.4$mean.diff[1]/diag$s.4$sd0[1]))
        srati <- mean(c(diag$s.1$sd.ratio[1],diag$s.2$sd.ratio[1],
                            diag$s.3$sd.ratio[1],diag$s.4$sd.ratio[1]))
        arati <- mean(c(diag$s.1$autocorr.ratio[1],diag$s.2$autocorr.ratio[1],
                            diag$s.3$autocorr.ratio[1],diag$s.4$autocorr.ratio[1]))
      }

      scorestats[i,] <- c(1-r.xval,mdiff,1-srati,1-arati,res.trend,ks,ar,
                          1-ds.ratio,
                          1- var(xval[,2])/var(xval[,1]))
      if (verbose) print('scorestats')
      if (verbose) print(scorestats[i,])

      quality <- 100*(1-mean(abs(scorestats[i,]),na.rm=TRUE))
     
      if (plot) {
        qcol <- quality
        qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
        cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)
        lines(annual(z),lwd=2,col=cols[qcol])
        lines(ya,type="b",pch=19)
        lines(dsa,lwd=2,col="grey")
      }
      R2 <- round(100*sd(xval[,2])/sd(xval[,1]),2)
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(100*r.xval,2),
                  "R2=",R2,'% ','Common EOF: bias=',round(mdiff,2),
                  ' sd1/sd2=',round(srati,3),
                  "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',round(quality)))
    }
  }
  if(verbose) print("Done with downscaling!")
  rm("DJFGCM")

  X <- zoo(t(X),order.by=t)
  colnames(X) <- gcmnm
  attr(X,"model_id") <- gcmnm
  #X <- attrcp(y,X)
  attr(X,'station') <- y
  attr(X,'predictor') <- attr(T2M,'source')
  attr(X,'domain') <- list(lon=lon,lat=lat)
  attr(X,'scorestats') <- scorestats
  attr(X,'path') <- path
  attr(X,'scenario') <- rcp
  attr(X,'history') <- history.stamp(y)
  if (non.stationarity.check)
    attr(X,'on.stationarity.check') <- difference.z else
    attr(X,'on.stationarity.check') <- NULL
  class(X) <- c("dsensemble","zoo")
  if (!is.null(path.ds)) file.ds <- file.path(path.ds,file.ds)
  save(file=file.ds,X)#"DSensemble.rda",X)
  print("---")
  invisible(X)
} 
#save(file=paste("dscmip5_",attr(y,'location'),"_",N,"_rcp4.5.rda",sep=""),rcp4.5)

#' @export DSensemble.precip
DSensemble.precip <- function(y,...,plot=TRUE,path="CMIP5.monthly/",rcp="rcp45",biascorrect=FALSE,
                              predictor="ERA40_pr_mon.nc",non.stationarity.check=FALSE,
                              ip=1:6,lon=c(-10,10),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                              select=NULL,FUN="wetmean",FUNX="sum",xfuns='C.C.eq',threshold=1,
                              pattern="pr_Amon_ens_",verbose=FALSE,nmin=NULL,ds.1900.2099=TRUE) {
  # FUN: exceedance, wetfreq, wet, dry

  if (verbose) print('DSensemble.precip')
#  if (deparse(substitute(FUN))=='spell') {
  if ( (FUN=='wet') | (FUN=='dry')) {
    y <- spell(y,threshold=threshold)
    y <- annual(y,nmin=nmin)
    #plot(y); 
    y <- switch(FUN,'wet'=subset(y,is=1),'dry'=subset(y,is=2))
  } else
#    y <- annual(y,FUN=FUN,threshold=threshold)
  # REB: the two next lines give errors... :(
#    if (sum(is.element(names(formals(FUN)),'threshold')==1))
#        y <- annual(y,FUN=FUN,threshold=threshold,nmin=nmin) else
        y <- annual(y,FUN=FUN,nmin=nmin)
  #
  index(y) <- year(y)
  
  if (!is.na(attr(y,'longitude')) & rel.cord)
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude')) & rel.cord)
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )

  if (sum(!is.finite(lon))>0) 
    warning(paste('Bad longitude range provided: ',paste(lon,collapse='-')))
  if (sum(!is.finite(lat))>0) 
    warning(paste('Bad latitude range provided: ',paste(lat,collapse='-')))
  
  # Get the predictor: ERA40
  if (verbose) print("predictor")
  if (is.character(predictor))
    pre <- retrieve(file=predictor,lon=lon,lat=lat,
                    verbose=verbose) else
  if (inherits(predictor,'field')) pre <- predictor
  rm("predictor"); gc(reset=TRUE)
  attr(pre,"source") <- "ERA40"

  # Use proportional variaions
  if (verbose) print("Annual mean")
  if (!is.annual(pre)) {
    if (sum(is.element(FUNX,xfuns))==0)  
      PREX <- annual(pre,FUN=FUNX,nmin=nmin) else
      eval(parse(text=paste('PREX <- annual(',FUNX,'(pre),FUN="mean",nmin=nmin)',sep="")))
  }
  #print("estimate %-changes")
  PRE.ref <- subset(PREX,it=1961:1990)
#  PRE <- zoo(100*coredata(PREX)/colMeans(coredata(PRE.ref)),order.by=year(PREX))
  if (is.precip(PREX)) {
    PRE <- zoo(100*coredata(PREX)/mean(c(coredata(subset(PREX,it=1961:1990)))),
             order.by=year(PREX))
    PRE <- attrcp(PREX,PRE)
    attr(PRE, "unit" ) <- "%"
    attr(PRE, "dimensions" ) <- attr(PREX, "dimensions" )
    class(PRE) <- class(PREX)
  } else  PRE <- PREX
  rm("PREX"); gc(reset=TRUE)
  
  if (verbose) print("graphics")
  unit <- attr(y,'unit')

  if (plot) {
    ylim <- switch(FUN,
                   'exceedance'=c(0,10),'wetmean'=c(0,10),
                   'wetfreq'=c(0,0),'spell'=c(0,0),
                   'mean'=c(-10,50),'sd'=c(-5,10),'ar1'=c(-0.5,0.7),
                   'HDD'=c(0,5000),'CDD'=c(0,500),'GDD'=c(0,2000))
    if (is.null(ylim)) ylim <- c(0,0)
    par(bty="n")
    plot.zoo(y,type="b",pch=19,main=attr(y,'location'),
             xlab="year",ylab=unit,
             sub=paste('Station: ',attr(y,'station_id'),'; coordinates: ',
             round(attr(y,'longitude'),4),'E/',
             round(attr(y,'latitude'),4),'N; ',
             attr(y,'altitude'),'m.a.s.l',sep=''),
             ylim=ylim + range(coredata(y),na.rm=TRUE),xlim=c(1900,2100))
    grid()
  }

  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.names=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) select <- 1:N else
                       N <- length(select)
  if (verbose) print(ncfiles[select])

  # set up results matrix and tables of diagnostics:
  years <- sort(1900:2100)
  m <- length(years)
  X <- matrix(rep(NA,N*m),N,m)
  gcmnm <- rep("",N)
  scorestats <- matrix(rep(NA,N*9),N,9)
  colnames(scorestats) <- c("1-r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration','1-R2')

  flog <- file("DSensemble.precip-log.txt","at")
  for (i in 1:N) {
    #
    gcm <- retrieve(file = ncfiles[select[i]],
                    lon=range(lon(PRE))+c(-2,2),lat=range(lat(PRE))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099)) else
                      gcm <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
    if (length(index(gcm))<=1) print(paste('Problem selecting GCM results in period',
                                           min(year(y),na.rm=TRUE),'2099'))
    # KMP: 10.03.2017 - pass on additional information about GCM runs (gcm + rip - realization, initialization, physics version)
    rip <- NULL
    if(any(grepl("rip", names(attributes(gcm))))) {
      nm.rip <- names(attributes(gcm))[grepl("rip",names(attributes(gcm)))][[1]]
      if(!is.null(attr(gcm, nm.rip))) {
        rip <- attr(gcm, nm.rip)
        if(!grepl("r[0-9]{1,2}i[0-9]{1,2}p[0-9]{1,2}", rip)) rip <- NULL
      }
    }
    if(is.null(rip)) {
      nm.r <- names(attributes(gcm))[grep("realization",names(attributes(gcm)))][[1]]
      nm.i <- names(attributes(gcm))[grep("initialization",names(attributes(gcm)))][[1]]
      nm.p <- names(attributes(gcm))[grep("physics",names(attributes(gcm)))][[1]]
      rip <- paste0("r",attr(gcm,nm.r),"i",attr(gcm,nm.i),"p",attr(gcm,nm.p))
    }
    gcmnm.i <- paste0(attr(gcm,'model_id'),".",rip)
    #gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
    #gcmnm[i] <- attr(gcm,'model_id')
    if (verbose) print(varid(gcm))
    
    if (sum(is.element(FUNX,xfuns))==0)
      GCMX <- annual(gcm,FUN=FUNX,nmin=nmin) else
      eval(parse(text=paste('GCMX <- annual(',FUNX,'(gcm),FUN="mean",nmin=nmin)',sep="")))
    if (is.precip(GCMX)) {
          GCM <- zoo(100*coredata(GCMX)/mean(c(coredata(subset(GCMX,it=1961:1990)))),
                     order.by=year(GCMX))
      GCM <- attrcp(GCMX,GCM)
      attr(GCM, "unit" ) <- "%"
      attr(GCM, "dimensions" ) <- attr(GCMX, "dimensions" )
      class(GCM) <- class(GCMX)
    } else
          GCM <- GCMX
    #
    #str(GCM)
    model.id <- attr(gcm,'model_id')
    rm("gcm","GCMX"); gc(reset=TRUE)
    if (verbose) print("combine")
    #
    PREGCM <- combine(PRE,GCM)
    if (verbose) print("EOF")
    Z <- EOF(PREGCM)
    
    # The test lines are included to assess for non-stationarity
    if (non.stationarity.check) {
      testGCM <- subset(GCM,it=range(year(PRE))) # REB 29.04.2014
      testy <- as.station(regrid(testGCM,is=y))  # REB 29.04.2014
      attr(testGCM,'source') <- 'testGCM'        # REB 29.04.2014
      testZ <- EOF(combine(testGCM,GCM))         # REB 29.04.2014
      rm("testGCM"); gc(reset=TRUE)
    }
    rm("GCM"); gc(reset=TRUE)
     
    # The test lines are included to assess for non-stationarity
    if (non.stationarity.check) {
      testds <- DS(testy,testZ,biascorrect=biascorrect,ip=ip)  # REB 29.04.2014
      testz <- attr(testds,'appendix.1')                     # REB 29.04.2014
      difference.z <- testy - testz                          # REB 29.04.2014
    }
    
    if (verbose) print("diagnose")
    diag <- diagnose(Z)
    if (biascorrect) Z <- biasfix(Z)
    if (verbose) print("- - - > DS (precip)")
    if (verbose) print(class(attr(Z,'appendix.1')))
    ds <- try(DS(y,Z,ip=ip,verbose=verbose))
    if (inherits(ds,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds[[1]],con=flog)
    } else {
      if (verbose) print("post-processing")
      z <- attr(ds,'appendix.1')
      i1 <- is.element(years,year(z))
      i2 <- is.element(year(z),years)
    #
      X[i,i1] <- z[i2]

    # Diagnose the residual: ACF, pdf, trend. These will together with the
    # cross-validation and the common EOF diagnostics provide a set of
    # quality indicators.
      cal <- coredata(attr(ds,"original_data"))
      fit <- coredata(attr(ds,"fitted_values"))
      res <- as.residual(ds)
      res.trend <- 10*diff(range(trend(res)))/diff(range(year(res)))
      ks <- ks.test(coredata(res),pnorm)$p.value
      ar <- as.numeric(acf(trend(cal-fit,result="residual"),plot=FALSE)[[1]][2])
    ## ar <- ar1(coredata(res))

    # Evaluation: here are lots of different aspects...
    # Get the diagnostics: this is based on the analysis of common EOFs...

      xval <- attr(ds,'evaluation')
      r.xval <- cor(xval[,1],xval[,2])

    #
      xy <- merge.zoo(z,y)
      ds.ratio <- sd(xy[,1],na.rm=TRUE)/sd(xy[,2],na.rm=TRUE)
    
    # Extract the mean score for leading EOF from the 4 seasons:
      mdiff <- diag$mean.diff[1]/diag$sd0[1]
      srati <- 1 - diag$sd.ratio[1]
      arati <- 1 - diag$autocorr.ratio[1]
      scorestats[i,] <- c(1-r.xval,mdiff,srati,arati,res.trend,ks,ar,1-ds.ratio,
      1-var(xval[,2])/var(xval[,1]))
      if (verbose) print(scorestats[i,])
      quality <- 100*(1-mean(abs(scorestats[i,]),na.rm=TRUE))

      index(z) <- year(z); index(ds) <- year(ds)
      if (plot) {
        qcol <- quality
        qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
        cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)
        lines(z,lwd=2,col=cols[qcol])
        lines(y,type="b",pch=19)
        lines(ds,lwd=2,col="grey")
     }
      #
      R2 <- round(100*sd(xval[,2])/sd(xval[,1]),2)
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(100*r.xval,2),
                  "R2=",R2,'% ', 'Common EOF: bias=',round(mdiff,2),
                  ' sd1/sd2=',round(srati,3),
                  "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',round(quality)))
    }
  }

  #
  X <- zoo(t(X),order.by=years)
  colnames(X) <- gcmnm
  attr(X,"model_id") <- gcmnm
  #X <- attrcp(y,X)
  attr(X,'station') <- y
  attr(X,'predictor') <- attr(PRE,'source')
  attr(X,'domain') <- list(lon=lon,lat=lat)
  attr(X,'scorestats') <- scorestats
  attr(X,'path') <- path
  attr(X,'scenario') <- rcp
  if (non.stationarity.check)
    attr(X,'on.stationarity.check') <- difference.z else
    attr(X,'on.stationarity.check') <- NULL
  attr(X,'history') <- history.stamp(y)
  class(X) <- c("dsensemble","zoo")
  save(file="DSensemble.rda",X)
  print("---")
  invisible(X)
}

#' @export DSensemble.annual
DSensemble.annual <- function(y,...,plot=TRUE,path="CMIP5.monthly/",rcp="rcp45",biascorrect=FALSE,
                              predictor="ERA40_t2m_mon.nc",non.stationarity.check=FALSE,
                              ip=1:6,lon=c(-10,10),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                              abscoords=FALSE,select=NULL,FUN=NULL,FUNX="mean",xfuns='C.C.eq',threshold=1,
                              pattern="tas_Amon_ens_",verbose=FALSE,nmin=NULL,ds.1900.2099=TRUE) {
  # FUN: exceedance, wetfreq, wet, dry
  
  if (verbose) print('DSensemble.annual')
#  if (deparse(substitute(FUN))=='spell') {
  index(y) <- year(y)
  
  if (!abscoords) {
    ## If relative coordinates:
    if (!is.na(attr(y,'longitude')) & rel.cord)
      lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
    if (!is.na(attr(y,'latitude')) & rel.cord)
      lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )
  }

  if (sum(!is.finite(lon))>0) 
    warning(paste('Bad longitude range provided: ',paste(lon,collapse='-')))
  if (sum(!is.finite(lat))>0) 
    warning(paste('Bad latitude range provided: ',paste(lat,collapse='-')))
  
  # Get the predictor: ERA40
  if (verbose) print("predictor")
  if (is.character(predictor))
    pre <- retrieve(file=predictor,lon=lon,lat=lat,
                    verbose=verbose) else
  if (inherits(predictor,'field')) pre <- predictor
  rm("predictor"); gc(reset=TRUE)
  attr(pre,"source") <- "ERA40"
  
  ## KMP 2017-09-12: don't use subset if pre is annual data and it is months or season!
  ## if it is character, then then extraction of months reduces number of
  ## months per year.
  if ((is.null(nmin)) & (is.character(it))) nmin <- length(it)
  if (!is.null(it) & !(is.character(it) & inherits(pre,"annual"))) {
    if (verbose) print('Extract some months or a time period')
    if (verbose) print(it)
    pre <- subset(pre,it=it)
  }

  # Use proportional variations
  ## KMP 2017-09-12: don't calculate the annual mean if pre is already annual
  if (verbose) print("Annual mean")
  if(inherits(pre,"annual")) {
    PRE <- pre
  } else {
    if (sum(is.element(FUNX,xfuns))==0) {
      PRE <- annual(pre,FUN=FUNX,nmin=nmin) 
    } else {
      eval(parse(text=paste('PRE <- annual(',FUNX,'(pre),FUN="mean",nmin=nmin)',sep="")))
    }
  }
  
  if (verbose) print("graphics")
  unit <- attr(y,'unit')

  if (plot) {
    ylim <- c(0,10)
    par(bty="n")
    plot.zoo(y,type="b",pch=19,main=attr(y,'location'),
             xlab="year",ylab=unit,
             sub=paste('Station: ',attr(y,'station_id'),'; coordinates: ',
             round(attr(y,'longitude'),4),'E/',
             round(attr(y,'latitude'),4),'N; ',
             attr(y,'altitude'),'m.a.s.l',sep=''),
             ylim=ylim + range(coredata(y),na.rm=TRUE),xlim=c(1900,2100))
    grid()
  }

  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.names=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) select <- 1:N else
                       N <- length(select)
  if (verbose) print(ncfiles[select])

  # set up results matrix and tables of diagnostics:
  years <- sort(1900:2100)
  m <- length(years)
  X <- matrix(rep(NA,N*m),N,m)
  gcmnm <- rep("",N)
  scorestats <- matrix(rep(NA,N*9),N,9)
  colnames(scorestats) <- c("1-r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration','1-R2')

  flog <- file("DSensemble.precip-log.txt","at")
  for (i in 1:N) {
    #
    gcm <- try(retrieve(file = ncfiles[select[i]],
      lon=range(lon(PRE))+c(-2,2),lat=range(lat(PRE))+c(-2,2),verbose=verbose))
    if(inherits(gcm,"try-error")) {
      writeLines(ncfiles[select[i]],con=flog)
    } else {
      if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099)) else
                        gcm <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
      if (length(index(gcm))<=1) print(paste('Problem selecting GCM results in period',
                                             min(year(y),na.rm=TRUE),'2099'))
      # KMP: 10.03.2017 - pass on additional information about GCM runs (gcm + rip - realization, initialization, physics version)
      rip <- NULL
      if(any(grepl("rip", names(attributes(gcm))))) {
        nm.rip <- names(attributes(gcm))[grepl("rip",names(attributes(gcm)))][[1]]
        if(!is.null(attr(gcm, nm.rip))) {
          rip <- attr(gcm, nm.rip)
          if(!grepl("r[0-9]{1,2}i[0-9]{1,2}p[0-9]{1,2}", rip)) rip <- NULL
        }
      }
      if(is.null(rip)) {
        nm.r <- names(attributes(gcm))[grep("realization",names(attributes(gcm)))][[1]]
        nm.i <- names(attributes(gcm))[grep("initialization",names(attributes(gcm)))][[1]]
        nm.p <- names(attributes(gcm))[grep("physics",names(attributes(gcm)))][[1]]
        rip <- paste0("r",attr(gcm,nm.r),"i",attr(gcm,nm.i),"p",attr(gcm,nm.p))
      }
      gcmnm[i] <- paste0(attr(gcm,'model_id'),".",rip)
      #gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
      #gcmnm[i] <- attr(gcm,'model_id')
      if (verbose) print(varid(gcm))
      
      if (!is.null(it)) {
        if (verbose) print('Extract some months or a time period')
        if (verbose) print(it)
        gcm <- subset(gcm,it=it)
      }
      
      if (sum(is.element(FUNX,xfuns))==0) {
        GCM <- annual(gcm,FUN=FUNX,nmin=nmin)
      } else {
        eval(parse(text=paste('GCM <- annual(',FUNX,'(gcm),FUN="mean",nmin=nmin)',sep="")))
      }
  
      model.id <- attr(gcm,'model_id')
  
      if (verbose) print("combine")
  
      PREGCM <- combine(PRE,GCM)
      if (verbose) print("EOF")
      Z <- EOF(PREGCM)
      
      # The test lines are included to assess for non-stationarity
      if (non.stationarity.check) {
        testGCM <- subset(GCM,it=range(year(PRE))) # REB 29.04.2014
        testy <- as.station(regrid(testGCM,is=y))  # REB 29.04.2014
        attr(testGCM,'source') <- 'testGCM'        # REB 29.04.2014
        testZ <- EOF(combine(testGCM,GCM))         # REB 29.04.2014
        rm("testGCM"); gc(reset=TRUE)
      }
      rm("GCM"); gc(reset=TRUE)
       
      # The test lines are included to assess for non-stationarity
      if (non.stationarity.check) {
        testds <- DS(testy,testZ,biascorrect=biascorrect,ip=ip)  # REB 29.04.2014
        testz <- attr(testds,'appendix.1')                     # REB 29.04.2014
        difference.z <- testy - testz                          # REB 29.04.2014
      }
      
      if (verbose) print("diagnose")
      diag <- diagnose(Z)
      if (biascorrect) Z <- biasfix(Z)
      if (verbose) print("- - - > DS (annual)")
      if (verbose) print(class(attr(Z,'appendix.1')))
      ds <- try(DS(y,Z,ip=ip,verbose=verbose))
      
      if (inherits(ds,"try-error")) {    
        writeLines(gcmnm[i],con=flog)
        writeLines(ds[[1]],con=flog)
      } else {
        if (verbose) print("post-processing")
        z <- attr(ds,'appendix.1')
        i1 <- is.element(years,year(z))
        i2 <- is.element(year(z),years)
        #
        X[i,i1] <- z[i2]
  
        # Diagnose the residual: ACF, pdf, trend. These will together with the
        # cross-validation and the common EOF diagnostics provide a set of
        # quality indicators.
        cal <- coredata(attr(ds,"original_data"))
        fit <- coredata(attr(ds,"fitted_values"))
        res <- as.residual(ds)
        res.trend <- 10*diff(range(trend(res)))/diff(range(year(res)))
        ks <- ks.test(coredata(res),pnorm)$p.value
        ar <- as.numeric(acf(trend(cal-fit,result="residual"),plot=FALSE)[[1]][2])
        ## ar <- ar1(coredata(res))

        # Evaluation: here are lots of different aspects...
        # Get the diagnostics: this is based on the analysis of common EOFs...
 
        xval <- attr(ds,'evaluation')
        r.xval <- cor(xval[,1],xval[,2])

        #
        xy <- merge.zoo(z,y)
        ds.ratio <- sd(xy[,1],na.rm=TRUE)/sd(xy[,2],na.rm=TRUE)
    
        # Extract the mean score for leading EOF from the 4 seasons:
        mdiff <- diag$mean.diff[1]/diag$sd0[1]
        srati <- 1 - diag$sd.ratio[1]
        arati <- 1 - diag$autocorr.ratio[1]
        scorestats[i,] <- c(1-r.xval,mdiff,srati,arati,res.trend,ks,ar,1-ds.ratio,
        1-var(xval[,2])/var(xval[,1]))
        if (verbose) print(scorestats[i,])
        quality <- 100*(1-mean(abs(scorestats[i,]),na.rm=TRUE))
  
        index(z) <- year(z); index(ds) <- year(ds)
        if (plot) {
          qcol <- quality
          qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
          cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)
          lines(z,lwd=2,col=cols[qcol])
          lines(y,type="b",pch=19)
          lines(ds,lwd=2,col="grey")
        }
        #
        R2 <- round(100*sd(xval[,2])/sd(xval[,1]),2)
        print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(100*r.xval,2),
                  "R2=",R2,'% ','Common EOF: bias=',
                  round(mdiff,2),' sd1/sd2=',round(srati,3),
                  "mean=",round(mean(coredata(y),na.rm=TRUE),2),
                  'quality=',round(quality)))
      }
    }
  }
  #
  X <- zoo(t(X),order.by=years)
  colnames(X) <- gcmnm
  attr(X,"model_id") <- gcmnm
  #X <- attrcp(y,X)
  attr(X,'station') <- y
  attr(X,"lon") <- attr(y,"lon")
  attr(X,"lat") <- attr(y,"lat")
  attr(X,"longname") <- attr(y,"longname")
  attr(X,'predictor') <- attr(PRE,'source')
  attr(X,'domain') <- list(lon=lon,lat=lat)
  attr(X,'scorestats') <- scorestats
  attr(X,'path') <- path
  attr(X,'scenario') <- rcp
  if (non.stationarity.check) {
    attr(X,'on.stationarity.check') <- difference.z 
  } else {
    attr(X,'on.stationarity.check') <- NULL
  }
  attr(X,'history') <- history.stamp(y)
  class(X) <- c("dsensemble","zoo")
  
  
  save(file="DSensemble.rda",X)
  print("---")
  invisible(X)
}

#' @export DSensemble.season
DSensemble.season <- function(y,...,season=NULL,plot=TRUE,path="CMIP5.monthly/",predictor="slp.mon.mean.nc",
                           rcp="rcp45",biascorrect=FALSE,non.stationarity.check=FALSE,
                           ip=1:6,lon=c(-20,20),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                           select=NULL,FUN="mean",FUNX="mean",xfuns='C.C.eq',
                           pattern="psl_Amon_ens_",lev=NULL,levgcm=NULL,path.ds=NULL,file.ds=NULL,
                           nmin=NULL,verbose=FALSE,ds.1900.2099=TRUE) {

  if(verbose) print("DSensemble.season")

  if(is.null(season)) {
    if(inherits(y,"season")) {
      season <- unique(season(y))
    } else {
      season <- "djf"
    }
  }
  
  if ((FUN=='sd') | (FUN =='ar1')) {
    y <- anomaly(y)
    attr(y,'aspect') <- 'original'
  }
       
  if(verbose) print("Set lon/lat predictor range")
  if ( !is.na(attr(y,'longitude'))[1] & (any(lon>0) & any(lon<0)) & rel.cord)
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if ( !is.na(attr(y,'latitude'))[1] & (any(lat>0) & any(lat<0)) & rel.cord)
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )

  if (sum(!is.finite(lon))>0) 
    warning(paste('Bad longitude range provided: ',paste(lon,collapse='-')))
  if (sum(!is.finite(lat))>0) 
    warning(paste('Bad latitude range provided: ',paste(lat,collapse='-')))
  
  if(verbose) print("Arrange units")
  if ( (unit(y)[1] == "deg C") | (unit(y)[1] == "degree Celsius") |
       (unit(y)[1] == "degC") | (unit(y)[1] == "Celsius") )
        unit <- expression(degree*C) else
        unit <- attr(y,'unit')
  if (verbose) print(paste('Units:',unit))
  
  if(verbose) print("aggregate time series")
  sm <- eval(parse(text=paste("season.abb()$",season,sep="")))
  s1 <- as.character(sm[1])
  s2 <- as.character(sm[length(sm)])
  if (nchar(s1)==1) s1 <- paste("0",s1,sep="")
  if (nchar(s2)==1) s2 <- paste("0",s2,sep="")
  if (is.null(nmin)) nmin <- length(sm)
  ys <- as.4seasons(y,start=paste(s1,"01",sep="-"),
                            end=paste(s2,"28",sep="-"),FUN=FUN)
  if(!is.null(attr(ys,"n.valid"))) ys <- ys[attr(ys,"n.valid")>=nmin]
  if (FUN=="sum" & grepl("month",attr(ys,"unit"))) {
    attr(ys,"unit") <- gsub("month","season",attr(ys,"unit"))
  }

  if (plot) {
    if(verbose) print("Plot station data (predictand)")
    ylim <- c(0,0)
    ylim <- switch(FUN,'mean'=c(-2,8),'sd'=c(-0.5,1),'ar1'=c(-0.5,0.7),
                   'sum'=c(-6,12))
    if (verbose) print(paste('set ylim based on "',FUN,
                             '" -> c(',ylim[1],', ',ylim[2],')',sep=''))
    par(bty="n")
    plot.zoo(ys,type="b",pch=19,main=attr(y,'location'),
             xlab="year",ylab=unit,
             sub=paste('Station: ',attr(y,'station_id'),'; coordinates: ',
             round(attr(ys,'longitude'),4),'E/',
             round(attr(ys,'latitude'),4),'N; ',
             attr(ys,'altitude'),'m.a.s.l',sep=''),
             ylim=ylim + range(coredata(ys),na.rm=TRUE),
             xlim=as.Date(c("1900-01-01","2100-12-31")))
    grid()
  }

  if(verbose) print("Retrieve predictor data")
  if (is.character(predictor))
    slp <- retrieve(file=predictor,lon=lon,lat=lat,lev=lev,
                    verbose=verbose) else
  if (inherits(predictor,'field'))
    slp <- subset(predictor,is=list(lon=lon,lat=lat))

  if(verbose) print("Aggregate seasonal values")
  SLP <- subset(as.4seasons(slp,FUN=FUNX,nmin=nmin),it=season)
  ok <- is.finite(rowSums(SLP))
  SLP <- subset(SLP,it=range(year(SLP)[ok]))
  rm("slp"); gc(reset=TRUE)
  
  # Ensemble GCMs
  if(verbose) print("Retrieve & arrange GCMs")
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.names=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) {
    select <- 1:N 
  } else {
    select <- select[select<=N]
    N <- length(select)
  }
  if (verbose) {print('GCMs:'); print(path); print(ncfiles[select])}

  if(verbose) print("Set up results matrix & table of diagnostics")
  ## KMP 2017-10-19: Don't need to keep all seasons in X
  months <- switch(season, "djf"=1, "mam"=4, "jja"=7, "son"=10)
  years <- sort(rep(1900:2100,length(months)))
  months <- rep(months,length(1900:2100))
  #years <- sort(rep(1900:2100,4))
  #months <- rep(c(1,4,7,10),length(1900:2100))
  m <- length(years)
  X <- matrix(rep(NA,N*m),N,m)
  gcmnm <- rep("",N)
  scorestats <- matrix(rep(NA,N*9),N,9)
  colnames(scorestats) <- c("1-r.xval","mean.diff","1-sd.ratio",
                            "1-autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration',
                            '1-R2')
  t <- as.Date(paste(years,months,'01',sep='-'))

  if(verbose) print("Quick test")  
  flog <- file("DSensemble.season-log.txt","at")

  if (verbose) print("loop...") 
  for (i in 1:N) {
    if (verbose) print(ncfiles[select[i]])
    gcm <- retrieve(file = ncfiles[select[i]],
                          lon=range(lon(SLP))+c(-2,2),lev=levgcm,
                          lat=range(lat(SLP))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099)) else
                      gcm <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
    if (length(index(gcm))<=1) print(paste('Problem selecting GCM results in period',
                                           min(year(y),na.rm=TRUE),'2099'))
    # KMP: 10.03.2017 - pass on additional information about GCM runs (gcm + rip - realization, initialization, physics version)
    rip <- NULL
    if(any(grepl("rip", names(attributes(gcm))))) {
      nm.rip <- names(attributes(gcm))[grepl("rip",names(attributes(gcm)))][[1]]
      if(!is.null(attr(gcm, nm.rip))) {
        rip <- attr(gcm, nm.rip)
        if(!grepl("r[0-9]{1,2}i[0-9]{1,2}p[0-9]{1,2}", rip)) rip <- NULL
      }
    }
    if(is.null(rip)) {
      nm.r <- names(attributes(gcm))[grep("realization",names(attributes(gcm)))][[1]]
      nm.i <- names(attributes(gcm))[grep("initialization",names(attributes(gcm)))][[1]]
      nm.p <- names(attributes(gcm))[grep("physics",names(attributes(gcm)))][[1]]
      rip <- paste0("r",attr(gcm,nm.r),"i",attr(gcm,nm.i),"p",attr(gcm,nm.p))
    }
    gcmnm[i] <- paste0(attr(gcm,'model_id'),".",rip)
    #gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
    GCM <- subset(as.4seasons(gcm,FUN=FUNX,nmin=nmin),it='djf')
    rm("gcm"); gc(reset=TRUE)
    SLPGCM <- combine(SLP,GCM)
    if (verbose) print("- - - > EOFs")
    Z <- EOF(SLPGCM)

    # The test lines are included to assess for non-stationarity
    if (non.stationarity.check) {
      testGCM <- subset(GCM,it=range(year(SLP))) # REB 29.04.2014
      testy <- as.station(regrid(testGCM,is=ys))  # REB 29.04.2014
      attr(testGCM,'source') <- 'testGCM'        # REB 29.04.2014
      testZ <- combine(testGCM,GCM)              # REB 29.04.2014
      rm("testGCM"); gc(reset=TRUE)
    }

    if (verbose) print("- - - > DS (seasonal)")
    if (verbose) print(class(attr(Z,'appendix.1')))
    if (biascorrect) try(Z <- biasfix(Z))
    ds <- try(DS(ys,Z,ip=ip))
    if (inherits(ds,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds[[1]],con=flog)
    } else {
      attr(ds,'evaluation') <- crossval(ds)
      if (verbose) print("post-processing")
      z <- attr(ds,'appendix.1')
      if (non.stationarity.check) {
        testds <- DS(testy,testZ,biascorrect=biascorrect,ip=ip)   # REB 29.04.2014
        testz <- attr(testds,'appendix.1')                      # REB 29.04.2014
        difference.z <- testy - testz                           # REB 29.04.2014
      }
      i1 <- is.element(paste(years,months,sep='-'),
                       paste(year(z),month(z),sep='-'))
      i2 <- is.element(paste(year(z),month(z),sep='-'),
                       paste(years,months,sep='-'))
      X[i,i1] <- z[i2]

    # Diagnose the residual: ACF, pdf, trend. These will together with the
    # cross-validation and the common EOF diagnostics provide a set of
    # quality indicators.
      cal <- coredata(attr(ds,"original_data"))
      fit <- coredata(attr(ds,"fitted_values"))
      res <- as.residual(ds)
      res.trend <- 10*diff(range(trend(res)))/diff(range(year(res)))
      ks <- round(ks.test(coredata(res),pnorm)$p.value,4)
      ar <- as.numeric(acf(trend(cal-fit,result="residual"),plot=FALSE)[[1]][2])

      if (verbose) print(paste("Examine residuals: trend=",
                               round(res.trend,3),'D/decade; K.S. p-val',
                               round(ks,2),'; AR(1)=',round(ar,2)))
    # Evaluation: here are lots of different aspects...
    # Get the diagnostics: this is based on the analysis of common EOFs...
      xval <- attr(ds,'evaluation')
      r.xval <- cor(xval[,1],xval[,2])
      if (verbose) print(paste("x-validation r=",r.xval))
    
      xy <- merge.zoo(z,ys)
      ds.ratio <- sd(xy[,1],na.rm=TRUE)/sd(xy[,2],na.rm=TRUE)
      if (verbose) print(paste("sd ratio=",ds.ratio))
  
    #print(names(attributes(ds)))
      if (biascorrect) {
        diag <- attr(ds,'diagnose')
        if ( (verbose) & !is.null(diag)) str(diag)
      } else diag <- NULL
    
    # diagnose for ds-objects
      ##
      if (verbose) print('...')
       if (is.null(diag)) {
        ##diag <- diagnose(z,plot=FALSE)
        scorestats[i,] <- c(1-r.xval,NA,NA,NA,res.trend,ks,1-ar,1-ds.ratio,
                            1-round(var(xval[,2])/var(xval[,1]),2))
        mdiff <- (mean(subset(ys,it=range(year(ds))),na.rm=TRUE)-
                  mean(subset(ds,it=range(year(ys))),na.rm=TRUE))/
                    sd(ys,na.rm=TRUE)
        srati <- sd(subset(ds,it=range(year(ys))),na.rm=TRUE)/
                 sd(subset(ys,it=range(year(ds))),na.rm=TRUE)
        arati <- ar1(zoo(ds,order.by=year(ds)))/ar1(zoo(ys,order.by=year(ys)))
      } else {

    # Extract the mean score for leading EOF from the 4 seasons:
        mdiff <- mean(c(diag$s.1$mean.diff[1]/diag$s.1$sd0[1],
                        diag$s.2$mean.diff[1]/diag$s.2$sd0[1],
                        diag$s.3$mean.diff[1]/diag$s.3$sd0[1],
                        diag$s.4$mean.diff[1]/diag$s.4$sd0[1]))
        srati <- mean(c(diag$s.1$sd.ratio[1],diag$s.2$sd.ratio[1],
                            diag$s.3$sd.ratio[1],diag$s.4$sd.ratio[1]))
        arati <- mean(c(diag$s.1$autocorr.ratio[1],diag$s.2$autocorr.ratio[1],
                        diag$s.3$autocorr.ratio[1],
                        diag$s.4$autocorr.ratio[1]))
      }

      scorestats[i,] <- c(1-r.xval,mdiff,1-srati,1-arati,res.trend,ks,ar,
                          1-ds.ratio,
                          1- var(xval[,2])/var(xval[,1]))
      if (verbose) print('scorestats')
      if (verbose) print(scorestats[i,])
      quality <- 100*(1-mean(abs(scorestats[i,]),na.rm=TRUE))
     
      if (plot) {
        qcol <- quality
        qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
        cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)
        lines(z,lwd=2,col=cols[qcol])
        lines(ys,type="b",pch=19)
        lines(ds,lwd=2,col="grey")
      }
      R2 <- round(100*sd(xval[,2])/sd(xval[,1]),2)
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(100*r.xval,2),
                  "R2=",R2,'% ','Common EOF: bias=',round(mdiff,2),
                  ' sd1/sd2=',round(srati,3),
                  "mean=",round(mean(coredata(ys),na.rm=TRUE),2),
                  'quality=',round(quality)))
    }
  }
  if(verbose) print("Done with downscaling!")
  rm("GCM")
  
  X <- zoo(t(X),order.by=t)
  colnames(X) <- gcmnm
  attr(X,"model_id") <- gcmnm
  #X <- attrcp(ys,X)
  attr(X,"season") <- season
  attr(X,"longname") <- attr(y,"longname")
  attr(X,'station') <- ys
  attr(X,'predictor') <- attr(SLP,'source')
  attr(X,'domain') <- list(lon=lon,lat=lat)
  attr(X,'scorestats') <- scorestats
  attr(X,'path') <- path
  attr(X,'scenario') <- rcp
  attr(X,'history') <- history.stamp(y)
  if (non.stationarity.check) {
    attr(X,'non.stationarity.check') <- difference.z
  } else {
    attr(X,'non.stationarity.check') <- NULL
  }
  class(X) <- c("dsensemble","season","zoo")
  if (is.null(file.ds)) {
    file.ds <- paste("DSensemble",rcp,N,attr(y,"variable"),
                     season,"rda",sep="")
  }
  if (!is.null(path.ds)) file.ds <- file.path(path.ds,file.ds)
  save(file=file.ds,X)
  print("---")
  invisible(X)
}


# DSensemble.mu <- function(y,plot=TRUE,path="CMIP5.monthly/",
#                           rcp="rcp45",biascorrect=FALSE,
#                           predictor=list(t2m="data/ncep/air.mon.mean.nc",
#                                          olr="data/ncep/OLR.mon.mean.nc",
#                                          slp="data/ncep/slp.mon.mean.nc"),
#                           non.stationarity.check=FALSE,
#                           ip=1:16,lon=c(-30,20),lat=c(-20,10),it=NULL,rel.cord=TRUE,
#                           select=NULL,FUN="wetmean",threshold=1,
#                           pattern=c("tas_Amon_ens_","olr_Amon_ens_","slp_Amon_ens_"),
#                           verbose=FALSE,nmin=365,ds.1900.2099=TRUE) {
# 
# # This function is for downscaling wet-day mean using a combination of predictors
# 
#   # Get the global mean temeprature: pentads
# 
#   # Or a combination of OLR + C.C.eq(t2m) + SLP
# 
# 
#     # FUN: exceedance, wetfreq, wet, dry
# 
#   if (verbose) print('DSensemble.mu')
#   print("DSensemble.mu is not finished. Try something else.")
#   # if (verbose) print(paste('The predictor: annual',FUN))
#   # if (!inherits(y,'annual')) y <- annual(y,FUN=FUN,threshold=threshold,nmin=nmin)
#   # index(y) <- year(y)
#   # 
#   # if (!is.na(attr(y,'longitude')) & rel.cord)
#   #   lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
#   # if (!is.na(attr(y,'latitude')) & rel.cord)
#   #   lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )
#   # 
#   # if (sum(!is.finite(lon))>0) 
#   #   warning(paste('Bad longitude range provided: ',paste(lon,collapse='-')))
#   # if (sum(!is.finite(lat))>0) 
#   #   warning(paste('Bad latitude range provided: ',paste(lat,collapse='-')))
#   # 
#   # # Get the predictor: NCEP/NCAR
#   # if (verbose) print(paste("Get the set of predictors:",names(predictor),collapse=' '))
#   # if (is.character(predictor[[1]])) pre1 <- retrieve(file=predictor[[1]],
#   #                                                    lon=lon,lat=lat,
#   #                                                    verbose=verbose)
#   # if (is.character(predictor[[2]])) pre2 <- retrieve(file=predictor[[2]],
#   #                                                    lon=lon,lat=lat,
#   #                                                    verbose=verbose)
#   # if (is.character(predictor[[3]])) pre3 <- retrieve(file=predictor[[3]],
#   #                                                    lon=lon,lat=lat,
#   #                                                    verbose=verbose)
#   # 
#   # # Combine the predictors
#   # if (verbose) print("Annual mean - predictors")
#   # PREX1 <- annual(C.C.eq(pre1),FUN='mean') # Clausius-Claperyron eq. -> sat. vapour pressure
#   # PREX2 <- annual(pre2,FUN='mean')
#   # PREX3 <- annual(pre3,FUN='mean')
#   # 
#   # if (verbose) print("graphics")
#   # unit <- attr(y,'unit')
#   # 
#   # if (plot) {
#   #   #ylim <- switch(deparse(substitute(FUN)),
#   #   ylim <- switch(FUN,
#   #                  'exceedance'=c(0,10),'wetmean'=c(0,10),
#   #                  'wetfreq'=c(0,0),'spell'=c(0,0),
#   #                  'mean'=c(-10,50),'sd'=c(-5,10),'ar1'=c(-0.5,0.7),
#   #                  'HDD'=c(0,5000),'CDD'=c(0,500),'GDD'=c(0,2000))
#   #   if (is.null(ylim)) ylim <- c(0,0)
#   #   par(bty="n")
#   #   plot.zoo(y,type="b",pch=19,main=attr(y,'location'),
#   #            xlab="year",ylab=unit,
#   #            sub=paste('Station: ',attr(y,'station_id'),'; coordinates: ',
#   #            round(attr(y,'longitude'),4),'E/',
#   #            round(attr(y,'latitude'),4),'N; ',
#   #            attr(y,'altitude'),'m.a.s.l',sep=''),
#   #            ylim=ylim + range(coredata(y),na.rm=TRUE),xlim=c(1900,2100))
#   #   grid()
#   # }
#   # 
#   # # Ensemble GCMs
#   # path <- file.path(path,rcp,fsep = .Platform$file.sep)
#   # ## KMP 2018-11-02: pattern1, pattern2, pattern3 have not been defined
#   # ## but ncfiles1, ncfiles2, ncfiles3 are used later
#   # #ncfiles <- list.files(path=path,pattern=pattern,full.names=TRUE)
#   # pattern1 <- pattern[1]
#   # pattern2 <- pattern[2]
#   # pattern3 <- pattern[3]
#   # ncfiles1 <- list.files(path=path,pattern=pattern1,full.names=TRUE)
#   # ncfiles2 <- list.files(path=path,pattern=pattern2,full.names=TRUE)
#   # ncfiles3 <- list.files(path=path,pattern=pattern3,full.names=TRUE)
#   # ncfiles <- ncfiles[1]
#   # 
#   # N <- length(ncfiles)
#   # if (is.null(select)) select <- 1:N else
#   #                      N <- length(select)
#   # if (verbose) print(ncfiles[select])
#   # 
#   # # set up results matrix and tables of diagnostics:
#   # years <- sort(1900:2100)
#   # m <- length(years)
#   # X <- matrix(rep(NA,N*m),N,m)
#   # gcmnm <- rep("",N)
#   # scorestats <- matrix(rep(NA,N*9),N,9)
#   # colnames(scorestats) <- c("r.xval","mean.diff","sd.ratio","autocorr.ratio",
#   #                           "res.trend","res.K-S","res.ar1",
#   #                           'amplitude.ration','1-R2')
#   # 
#   # flog <- file("DSensemble.precip-log.txt","at")
#   # dse <- list(description='DSensemble.mu')
#   # 
#   # for (i in 1:N) {
#   #   ## Need to ensure that the different predictor files match...
#   #   print(paste(i,N,ncfiles1[select[i]],ncfiles2[select[i]],ncfiles3[select[i]]))
#   #   gcm1 <- retrieve(file = ncfiles1[select[i]],
#   #                   lon=range(lon(PRE1))+c(-2,2),lat=range(lat(PRE1))+c(-2,2),verbose=verbose)
#   #   if (ds.1900.2099) gcm1 <- subset(gcm1,it=c(1900,2099)) else
#   #                     gcm1 <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
#   #   if (length(index(gcm1))<=1) print(paste('Problem selecting GCM results in period',
#   #                                          min(year(y),na.rm=TRUE),'2099'))
#   #   gcm2 <- retrieve(file = ncfiles2[select[i]],
#   #                   lon=range(lon(PRE2))+c(-2,2),lat=range(lat(PRE2))+c(-2,2),verbose=verbose)
#   #   if (ds.1900.2099) gcm2 <- subset(gcm2,it=c(1900,2099)) else
#   #                     gcm2 <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
#   #   gcm3 <- retrieve(file = ncfiles3[select[i]],
#   #                   lon=range(lon(PRE3))+c(-2,2),lat=range(lat(PRE3))+c(-2,2),verbose=verbose)
#   #   if (ds.1900.2099) gcm3 <- subset(gcm3,it=c(1900,2099)) else
#   #                     gcm3 <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
#   #   # KMP: 10.03.2017 - pass on additional information about GCM runs (gcm + rip - realization, initialization, physics version)
#   #   gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'parent_experiment_rip'),sep="-")
#   #   #gcmnm[i] <- paste(attr(gcm1,'model_id'),attr(gcm,'realization'),sep="-")
#   #   #gcmnm[i] <- attr(gcm,'model_id')
#   #   if (verbose) print(varid(gcm1))
#   #   
#   #   GCM1 <- annual(C.C.eq(gcm1),FUN='mean')
#   #   GCM2 <- annual(gcm2,FUN='mean')
#   #   GCM3 <- annual(gcm3,FUN='mean')
#   # 
#   #   model.id <- attr(gcm1,'model_id')
#   #   rm("gcm","GCMX"); gc(reset=TRUE)
#   #   if (verbose) print("combine the three predictors")
#   #   #
#   #   PREGCM1 <- combine(PRE1,GCM1)
#   #   PREGCM2 <- combine(PRE2,GCM2)
#   #   PREGCM3 <- combine(PRE3,GCM3)
#   #   if (verbose) print("EOF")
#   #   Z1 <- EOF(PREGCM1)
#   #   Z2 <- EOF(PREGCM2)
#   #   Z3 <- EOF(PREGCM3)
#   #   
#   #   if (verbose) print("diagnose")
#   #   diag1 <- diagnose(Z1)
#   #   diag2 <- diagnose(Z2)
#   #   diag3 <- diagnose(Z3)
#   # 
#   #   if (biascorrect) 
#   #     X <- list(Z1=biasfix(Z1),
#   #               Z2=biasfix(Z2),
#   #               Z3=biasfix(Z3)) else
#   #     X < list(Z1=Z1,
#   #              Z2=Z2,
#   #              Z3=Z3)
#   #   x <- zoo(X[[1]])
#   #   np <- length(names(x))
#   #   for (i in 2:np) {
#   #       x <- merge(x,zoo(X[[i]]),all=TRUE)
#   #       w <- c(w,attr(X[[i]],'eigenvalues')/sum(attr(X[[i]],'eigenvalues')))
#   #       id <- c(id,rep(i,length(attr(X[[i]],'eigenvalues'))))
#   #     }
#   #   if (verbose) print(c(dim(x),length(w)))
#   #   t <- index(x)
#   #   ## apply the weights
#   #   x <- x %*% diag(w)
#   #   xm <- rowMeans(x)
#   #   x <- x[is.finite(xm),]; t <- t[is.finite(xm)]
#   # 
#   #   ## Apply an SVD to the combined PCs to extract the common signal in the
#   #   ## different predictors - these are more likely to contain real physics
#   #   ## and be related to the predictand.
#   #   if (verbose) print('svd')
#   #   udv <- svd(coredata(x))
#   #   if (verbose) print(summary(udv))
#   # 
#   #   ## If the predictor is a common EOF, then also combine the appended fields
#   #   ## the same way as the original flield.
#   #   if (inherits(X[[1]],'comb')) {
#   #       z <- z %*% diag(w)
#   #       udvz <- svd(coredata(z))
#   #   }
#   # 
#   #   eof <- zoo(udv$u[,1:20],order.by=t)
#   #   ## Let the pattern contain the weights for the EOFs in the combined
#   #   ## PC matrix, rather than spatial patterns. The spatial patterns are
#   #   ## then reconstructed from these.
#   #   pattern <- matrix(rep(1,length(udv$v[,1:20])),dim(udv$v[,1:20]))
#   # 
#   #   ## Do a little 'fake': here the pattern is not a geographical map but weight
#   #   ## for the EOFs.
#   #   dim(pattern) <- c(1,dim(pattern))
#   #   if (verbose) str(pattern)
#   #   attr(eof,'eigenvalues') <- udv$d[1:20]
#   #   attr(eof,'pattern') <- rep(1,20)
#   #   names(eof) <- paste("X.",1:20,sep="")
#   #   
#   #   class(eof) <- class(X[[1]])
#   # 
#   #   ## Downscale the results:
#   #   if (verbose) print("- - - > DS")
#   #   ds <- try(DS(y,eof,ip=ip,verbose=verbose))
#   #   if (inherits(ds,"try-error")) {    
#   #     writeLines(gcmnm[i],con=flog)
#   #     writeLines(ds[[1]],con=flog)
#   #   } else {
#   #     if (verbose) print("post-processing")
#   #     z <- attr(ds,'appendix.1')
#   #     i1 <- is.element(years,year(z))
#   #     i2 <- is.element(year(z),years)
#   #   #
#   # 
#   #   # Diagnose the residual: ACF, pdf, trend. These will together with the
#   #   # cross-validation and the common EOF diagnostics provide a set of
#   #   # quality indicators.
#   #     cal <- coredata(attr(ds,"original_data"))
#   #     fit <- coredata(attr(ds,"fitted_values"))
#   #     res <- as.residual(ds)
#   #     res.trend <- 10*diff(range(trend(res)))/diff(range(year(res)))
#   #     ks <- ks.test(coredata(res),pnorm)$p.value
#   #     ar <- as.numeric(acf(trend(cal-fit,result="residual"),plot=FALSE)[[1]][2])
#   #   ## ar <- ar1(coredata(res))
#   # 
#   #   # Evaluation: here are lots of different aspects...
#   #   # Get the diagnostics: this is based on the analysis of common EOFs...
#   # 
#   #     xval <- attr(ds,'evaluation')
#   #     r.xval <- cor(xval[,1],xval[,2])
#   # 
#   #   #
#   #     xy <- merge.zoo(z,y)
#   #     ds.ratio <- sd(xy[,1],na.rm=TRUE)/sd(xy[,2],na.rm=TRUE)
#   #   
#   #   # Extract the mean score for leading EOF from the 4 seasons:
#   #     mdiff <- diag$mean.diff[1]/diag$sd0[1]
#   #     srati <- 1 - diag$sd.ratio[1]
#   #     arati <- 1 - diag$autocorr.ratio[1]
#   #     
#   #     attr(z,'scorestats') <- c(1-r.xval,mdiff,srati,arati,res.trend,ks,ar,
#   #                               1-ds.ratio,1-var(xval[,2])/var(xval[,1]))
#   #     if (verbose) print('scorestats')
#   #     if (verbose) print(attr(z,'scorestats'))
#   #     dse[[i]] <- z
#   #     
#   #     quality <- 100*(1-mean(abs(scorestats[i,]),na.rm=TRUE))
#   #     index(z) <- year(z); index(ds) <- year(ds)
#   # 
#   #     if (plot) {
#   #       qcol <- quality
#   #       qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
#   #       cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)
#   #       lines(z,lwd=2,col=cols[qcol])
#   #       lines(y,type="b",pch=19)
#   #       lines(ds,lwd=2,col="grey")
#   #    }
#   #     #
#   #     R2 <- round(100*sd(xval[,2])/sd(xval[,1]),2)
#   #     print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(100*r.xval,2),
#   #                 "R2=",R2,'% ','Common EOF: bias=',round(mdiff,2),
#   #                 ' sd1/sd2=',round(srati,3),
#   #                 "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',round(quality)))
#   #   }
#   # }
#   # 
#   # #
#   # names(dse) <- gcmnm
#   # #X <- attrcp(y,X)
#   # attr(dse,'station') <- y
#   # attr(dse,'predictor') <- c(attr(PRE1,'source'),attr(PRE2,'source'),attr(PRE3,'source'))
#   # attr(dse,'domain') <- list(lon=lon,lat=lat)
#   # attr(dse,'scenario') <- rcp
#   # attr(dse,'history') <- history.stamp(y)
#   # class(X) <- c("dsensemble","zoo")
#   # save(file="DSensemble.rda",X)
#   # print("---")
#   # invisible(dse)
# }

#' @export DSensemble.mu.worstcase
DSensemble.mu.worstcase <- function(y,...,plot=TRUE,path="CMIP5.monthly/",predictor="ERA40_t2m_mon.nc",
                                    rcp="rcp45",biascorrect=FALSE,n=6,lon=c(-20,20),lat=c(-10,10),
                                    it=NULL,rel.cord=TRUE,select=NULL,FUN="wetmean",
                                    pattern="tas_Amon_ens_",mask=FALSE,verbose=FALSE,ds.1900.2099=TRUE) {
  if (verbose) print('DSensemble.mu.worstcase')

  ## The predictor is based on the seasonal variations and assumes that the seasnoal cycle in the
  ## wet-day mean mu is follows a systematic dependency to the seasonal variations in the temperature
  ## - the calibration uses the Clausius Clapeiron equation to estimate the saturation water vapour
  ## rather than using the temeprature directly.
  if (verbose) print(paste('The predictand: seasonal',FUN))

  ys <- aggregate(y,by=month,FUN=FUN)
  ya <- aggregate(y,by=year,FUN=FUN)
  ns <- 1
  
  ## A group of stations - use PCA
  if (!is.null(dim(ys))) {
    if (verbose) print('Use PCA for multiple stations')
    pca.ys <- PCA(ys,n=n)
    pca.ya <- PCA(ya)
    ns <- n
  }
  
  if (!is.na(attr(y,'longitude')) & rel.cord)
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude')) & rel.cord)
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )

  if (sum(!is.finite(lon))>0) 
    warning(paste('Bad longitude range provided: ',paste(lon,collapse='-')))
  if (sum(!is.finite(lat))>0) 
    warning(paste('Bad latitude range provided: ',paste(lat,collapse='-')))
  
  if (verbose) print("predictor")
  if (is.character(predictor)) {
    pre <- retrieve(file=predictor,lon=lon,lat=lat,verbose=verbose)
    if (mask) pre=mask(pre,land=TRUE)
    # KMP 2019-05-28: replaced spatial.avg.field with aggregate.area
    #pre <- spatial.avg.field(C.C.eq(pre)) 
    pre <- aggregate.area(C.C.eq(pre), FUN="mean")
  } else if (inherits(predictor,'field')) {
    # KMP 2019-05-28: replaced spatial.avg.field with aggregate.area
    #pre <- spatial.avg.field(predictor)
    pre <- aggregate.area(predictor, FUN="mean")
  }
  rm("predictor"); gc(reset=TRUE)
  ## Estimate the reference level
  normal61.90 <- mean(coredata(subset(pre,it=c(1961,1990))))
  x <- aggregate(pre,by=month,FUN="mean")

  if (verbose) print('Prepare the calibration data')
  ## Loop over PCAs if multipple stations

  if (n>1) results <- list()

  ## Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.names=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) select <- 1:N else
                                 N <- length(select)
  if (verbose) print(ncfiles[select])

  ## set up results matrix and tables of diagnostics:
  years <- sort(1900:2100)
  m <- length(years)
  gcmnm <- rep("",N)
  
  for (is in 1:n) {
    X <- matrix(rep(NA,N*m),N,m)
    if (verbose) print(paste('is=',is,'n=',n))
    if (n==1) cal <- data.frame(y=coredata(ys),x=coredata(x)) else
              cal <- data.frame(y=coredata(ys)[,is],x=coredata(x))
    attributes(cal$y) <- NULL
                        
    if (is.null(dim(ys))) stats <- cor.test(cal$y,cal$x) else
    stats <- cor.test(as.matrix(cal)[,1],cal$x)
    wc.model <- lm(y ~ x, data=cal)
    if (plot) {
      par(bty='n',cex.sub=0.7,col.sub='grey40')
      ylim <- range(cal$y,na.rm=TRUE); xlim=range(cal$x,na.rm=TRUE); dy <- diff(ylim)/25
      plot(cal$x,cal$y,pch=19,cex=1.5,col='grey',
           ylab=expression(paste(mu,' (mm/day)')),
           xlab=expression(paste(e[s],' (Pa)')),
           ylim=ylim,xlim=xlim,
           main='Worst-case based on seasonal variations',
           sub=paste(loc(y),' (',round(lon(y),2),'E/',round(lat(y),2),'N; ',alt(y),'m.a.s.l.)',sep=''))
      segments(x0=cal$x,y0=cal$y,x1=cal$x,y1=cal$y+2*attr(ys,'standard.error'),col='grey')
      segments(x0=cal$x,y0=cal$y,x1=cal$x,y1=cal$y-2*attr(ys,'standard.error'),col='grey')
      segments(x0=cal$x,y0=cal$y,x1=cal$x+2*attr(x,'standard.error'),y1=cal$y,col='grey')
      segments(x0=cal$x,y0=cal$y,x1=cal$x-2*attr(x,'standard.error'),y1=cal$y,col='grey')
      points(cal$x,cal$y,pch=19,cex=1.5,col='grey')
      grid()
      abline(wc.model)
      text(xlim[1],ylim[2],paste('Correlation=',round(stats$estimate,2),
                                 '(','p-value=',100*round(stats$p.value,4),'%)'),
           pos=4,cex=0.7,col='grey')
      text(xlim[1],ylim[2]-dy,paste('Regression: y=',round(wc.model$coeff[1],4), '+',
                                    round(wc.model$coeff[2],4), 'x (R2=',
                                    round(summary(wc.model)$r.squared,2),')'),
           pos=4,cex=0.7,col='grey')
      par(new=TRUE,fig=c(0.5,0.97,0.1,0.5),yaxt='n',xpd=TRUE,cex.axis=0.7,col.axis='grey')
      plot((cal$x - mean(cal$x))/sd(cal$x),type='l',lwd=2,ylab='',xlab='',col=rgb(0.6,0.3,0))
      axis(1,col='grey')
      lines((cal$y - mean(cal$y))/sd(cal$y),type='l',lwd=2,col=rgb(0,0.3,0.6))
      #dev.copy2eps(file='DSensemble.mu.worstcase.cal.eps')
    }

    if (plot) {
      dev.new()
      if (n>1) ya <- pca.ya[,is]
      plot(ya,xlim=c(1900,2100),
           ylim=range(aggregate(y,by=year,FUN='wetmean'),na.rm=TRUE)*c(0.75,1.5))
      grid()
    }
  
    sd.noise <- max(attr(ys,'standard.error'),sd(wc.model$residuals))
    for (i in 1:N) {
      if (verbose) print(ncfiles[select[i]])
      gcm <- retrieve(file = ncfiles[select[i]],lon=lon,lat=lat,
                      verbose=FALSE)
      if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099)) else
                        gcm <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
      if (length(index(gcm))<=1) print(paste('Problem selecting GCM results in period',
                                             min(year(y),na.rm=TRUE),'2099'))
      if (verbose) print(paste('mask=',mask))
      if (mask) gcm <- mask(gcm,land=TRUE)
      # KMP: 10.03.2017 - pass on additional information about GCM runs (gcm + rip - realization, initialization, physics version)
      rip <- NULL
      if(any(grepl("rip", names(attributes(gcm))))) {
        nm.rip <- names(attributes(gcm))[grepl("rip",names(attributes(gcm)))][[1]]
        if(!is.null(attr(gcm, nm.rip))) {
          rip <- attr(gcm, nm.rip)
          if(!grepl("r[0-9]{1,2}i[0-9]{1,2}p[0-9]{1,2}", rip)) rip <- NULL
        }
      }
      if(is.null(rip)) {
        nm.r <- names(attributes(gcm))[grep("realization",names(attributes(gcm)))][[1]]
        nm.i <- names(attributes(gcm))[grep("initialization",names(attributes(gcm)))][[1]]
        nm.p <- names(attributes(gcm))[grep("physics",names(attributes(gcm)))][[1]]
        rip <- paste0("r",attr(gcm,nm.r),"i",attr(gcm,nm.i),"p",attr(gcm,nm.p))
      }
      gcmnm[i] <- paste0(attr(gcm,'model_id'),".",rip)
      #gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
      if (verbose) print('spatial average')
      # KMP 2019-05-28: replaced spatial.avg.field with aggregate.area
      #GCM <- spatial.avg.field(C.C.eq(gcm))
      GCM <- aggregate.area(C.C.eq(gcm), FUN='mean')
      z <- annual(GCM,FUN="max")
      z <- z - mean(coredata(subset(z,it=c(1961,1990)))) + normal61.90
      i1 <- is.element(year(z),years)
      i2 <- is.element(years,year(z))
      if (verbose) print(c(i,sum(i1),sum(i2)))
      prex <- data.frame(x=coredata(z[i1]))
      if (verbose) print(summary(prex))
      if (verbose) print('prediction')
      z.predict <- predict(wc.model, newdata=prex) +
                         rnorm(n=sum(i1),sd=sd.noise)
      if (length(z.predict) != sum(i2)) {
        print('problem discovered')
        browser()
      }
      X[i,i2] <- z.predict
      if (verbose) print(paste("i=",i,"GCM=",gcmnm[i],sum(i2)))
      #if (sum(i2) != length(years)) 
      if (plot) lines(years[i2],X[i,i2],col=rgb(0,0.3,0.6,0.2))
    }
    if (plot) {
      if (n==1) lines(aggregate(y,by=year,FUN='wetmean'),col='red',lwd=3) else
                lines(PCA(aggregate(y,by=year,FUN='wetmean'))[,is],col='red',lwd=3)
    }
  
    X <- zoo(t(X),order.by=years)
    colnames(X) <- gcmnm
    attr(X,"model_id") <- gcmnm
    attr(X,'station') <- aggregate(y,by=year,FUN='wetmean')
    attr(X,'predictor') <- attr(pre,'source')
    attr(X,'domain') <- list(lon=lon,lat=lat)
    attr(X,'path') <- path
    attr(X,'scenario') <- rcp
    attr(X,'history') <- history.stamp(y)
    class(X) <- c("dsensemble","zoo")
    save(file="DSensemble.rda",X)
    if (verbose) print("--- end of iteration")
    if (n>1) eval(parse(text=paste('results$pca.',is,' <- X',sep='')))
  }
  if (n>1) X <- results
  if (verbose) print(names(X))
  if (verbose) print("--- Exit DSensemble.mu.worstcase")
  invisible(X)   
}

#' @export DSensemble.pca
DSensemble.pca <- function(y,...,plot=TRUE,path="CMIP5.monthly/",rcp="rcp45",biascorrect=FALSE,
                           predictor="ERA40_t2m_mon.nc",non.stationarity.check=FALSE,
                           ip=1:16,lon=c(-30,20),lat=c(-20,10), it=NULL,rel.cord=TRUE,
                           select=NULL,FUN="mean",rmtrend=TRUE,FUNX="mean",xfuns='C.C.eq',
                           threshold=1,pattern="tas_Amon_ens_",verbose=FALSE,
                           file.ds="DSensemble.rda",path.ds=NULL,nmin=NULL,ds.1900.2099=TRUE,test=FALSE) {

  if (verbose) print('DSensemble.pca')
  cls <- class(y)

  # This function is for downscaling PCA to represent a group of stations
  if (!is.na(attr(y,'longitude'))[1] & (any(lon>0) & any(lon<0)) & rel.cord)
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude'))[1] & (any(lat>0) & any(lat<0)) & rel.cord)
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )

  if (sum(!is.finite(lon))>0) 
    warning(paste('Bad longitude range provided: ',paste(lon,collapse='-')))
  if (sum(!is.finite(lat))>0) 
    warning(paste('Bad latitude range provided: ',paste(lat,collapse='-')))
  
  if (is.character(predictor)) {
    if (verbose) print('retrieve the predictor from netCDF file')
    t2m <- retrieve(file=predictor,lon=lon,lat=lat,
                    verbose=verbose)
      if (!is.null(it)) {
        if (verbose) print('Extract some months or a time period')
        if (verbose) print(it)
        t2m <- subset(t2m,it=it,verbose=verbose)
        ## if it is character, then then extraction of months reduces number of
        ## days per year.
      }
  } else if (inherits(predictor,'field')) {
    if (verbose) print('use the predictor provided as an argument')
    t2m <- predictor
    lon <- range(lon(t2m))
    lat <- range(lat(t2m))
    if (!is.annual(t2m) & !is.null(it)) {
      if (verbose) print('Extract months/time period:')
      if (verbose) print(it)
      t2m <- subset(t2m,it=it,verbose=verbose)
    }
  }
  ## If some months are selected, make sure that the minimum number of months
  ## required in the annual aggregation is updated
  if ((is.null(nmin)) & (is.character(it))) nmin <- length(it)
  
  if (inherits(y,'season')) {
    if (verbose) print('seasonal data found in the predictand')
    if (FUNX !='C.C.eq') {
      if (verbose) print(paste('apply',FUNX,'to the predictor'))
      T2M <- as.4seasons(t2m,FUN=FUNX,nmin=nmin) 
    } else {
      if (verbose) print('apply C.C.eq to the predictor:')
      eval(parse(text=paste('T2M <- as.4seasons(',FUNX,'(t2m),FUN="mean",nmin=nmin)',sep="")))
    }
    T2M <- matchdate(T2M,y,verbose=verbose)

    # Recursive: do each season seperately if there are more than one season
    if (length(table(season(y)))>1) {
      if (verbose) print('--- Apply DS to seasons seperately ---')
      Z <- list(info=paste('DSensemble.pca for different seasons: ',
                           paste(lon,collapse='-'),'E/',paste(lat,collapse='-'),'N',sep=''))
      for (season in names(table(season(y)))) {
        if (verbose) print(paste('Select',season))
        z <- DSensemble.pca(subset(y,it=season),plot=plot,path=path,
                            rcp=rcp,biascorrect=biascorrect,predictor=T2M,
                            non.stationarity.check=non.stationarity.check,
                            ip=ip,lon=lon,lat=lat,rel.cord=FALSE,
                            select=select,FUN=FUN,rmtrend=rmtrend,
                            FUNX=FUNX,xfuns=xfuns,threshold=threshold,
                            pattern=pattern,verbose=verbose,nmin=nmin)
        eval(parse(text=paste('Z$',season,' <- z',sep='')))
      }
      if (verbose) print('--- Results returned as a list ---')
      return(Z)
    }

  } else if (inherits(y,'annual')) {
    if (verbose) print('annual data')
    if (!inherits(predictor,'field')) {
      T2M <- t2m 
    } else if (!is.annual(t2m)) {
      if (verbose) print(paste('Aggregate annually',FUNX,'for calibration'))
      if (sum(is.element(FUNX,xfuns))==0) {
        T2M <- annual(t2m,FUN=FUNX,nmin=nmin) 
      } else {
        eval(parse(text=paste('T2M <- annual(',FUNX,'(t2m),FUN="mean",nmin=nmin)',sep="")))
      }
    } else {
      T2M <- t2m
    }
    ## Match the date
    T2M <- matchdate(T2M,y)
  } else if (inherits(y,'month')) {
    if (verbose) print('monthly data')
    T2M <- matchdate(t2m,y)
    ##if (FUNX=='C.C.eq') 
    ##  T2M <- mask(T2M,land=TRUE)
  }
  if (inherits(T2M,"eof")) T2M <- as.field(T2M)
  rm("predictor","t2m"); gc(reset=TRUE)
  if (verbose) {print('Check T2M:'); print(class(T2M)); print(index(T2M))}

  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.names=TRUE)
  N <- length(ncfiles)

  if (!is.null(select)) {
    select <- select[select<=N]
    N <- length(select)
  } else {
    select <- 1:N
  }
  if (verbose) {print('GCMs:'); print(path); print(ncfiles[select])}

  d.y <- dim(y)
  years <- 1900:2100
  m <- length(years)
  months <- rep(month(y)[1],m)
  X <- matrix(rep(NA,N*m*d.y[2]),N,m*d.y[2])
  dim(X) <- c(d.y[2],N,m)
  gcmnm <- rep("",N)
  scorestats <- matrix(rep(NA,N*9),N,9)
  colnames(scorestats) <- c("1-r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration',
                            "1-R2")

  t <- as.Date(paste(years,months,'01',sep='-'))

  if (plot) {
    par(bty='n')
    index(y) <- year(y)
    plot.zoo(y[,1],lwd=3,main='PC1',ylab='',xlab='',xlim=range(years))
  }
  
  flog <- file("DSensemble.pca-log.txt","at")
  
  ## Set up a list environment to keep all the results
  dse.pca <- list(info=paste('DSensemble.pca for different seasons: ',
                             paste(lon,collapse='-'),'E/',paste(lat,collapse='-'),'N',sep=''),
                  pca=y) ## KMP 06-08-2015
  if (verbose) print("loop...") 

  for (i in 1:N) {
    if (verbose) print(ncfiles[select[i]])
    gcm <- retrieve(file = ncfiles[select[i]],
                          lon=range(lon(T2M))+c(-2,2),
                          lat=range(lat(T2M))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099)) else
                      gcm <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
    if (length(index(gcm))<=1) print(paste('Problem selecting GCM results in period',
                                           min(year(y),na.rm=TRUE),'2099'))
    
    if (!is.null(it)) {
      if (verbose) print('Extract some months or a time period')
      if ( is.null(nmin) & is.character(it[1]) ) warning(paste("The argument 'it' is set but not 'nmin'; it=",
                                       paste(it,collapse="-")))
      if (verbose) print(it)
      gcm <- subset(gcm,it=it,verbose=verbose)
    }
    rip <- NULL
    if(any(grepl("rip", names(attributes(gcm))))) {
      nm.rip <- names(attributes(gcm))[grepl("rip",names(attributes(gcm)))][[1]]
      if(!is.null(attr(gcm, nm.rip))) {
        rip <- attr(gcm, nm.rip)
        if(!grepl("r[0-9]{1,2}i[0-9]{1,2}p[0-9]{1,2}", rip)) rip <- NULL
      }
    }
    if(is.null(rip)) {
      if(any(grepl("realization", names(attributes(gcm)))))
        nm.r <- names(attributes(gcm))[grep("realization",names(attributes(gcm)))][[1]] else
        nm.r <- 'NA'
      if(any(grepl("initialization", names(attributes(gcm)))))
        nm.i <- names(attributes(gcm))[grep("initialization",names(attributes(gcm)))][[1]] else
        nm.i <- 'NA'
      if(any(grepl("physics", names(attributes(gcm)))))
        nm.p <- names(attributes(gcm))[grep("physics",names(attributes(gcm)))][[1]] else
        nm.p <- 'NA'
      rip <- paste0("r",attr(gcm,nm.r),"i",attr(gcm,nm.i),"p",attr(gcm,nm.p))
    }
    gcmnm.i <- paste0(attr(gcm,'model_id'),".",rip)
    if (verbose) {
        print(paste('Extract month/season/annual data nmin=',nmin))
        print(class(y))
        print(FUNX)
    }
    if (inherits(y,'season')) {
      if (sum(is.element(FUNX,xfuns))==0) {
          if (verbose) print(paste('No special transformation (PCA)',FUNX,nmin)) 
          GCM <- as.4seasons(gcm,FUN=FUNX,nmin=nmin)
       } else {
           if (verbose) print('Need to aggregate FUNX(gcm)')
           eval(parse(text=paste('GCM <- as.4seasons(',FUNX,'(gcm),FUN="mean",nmin=nmin)',sep="")))
       }
      if (verbose) {print('Check: index(T2M)'); print(index(T2M))}
      GCM <- subset(GCM,it=season(T2M)[1])
    } else if (inherits(y,'annual')) {
      if (verbose) print(paste('Annualy aggregated',FUNX,'for GCM'))
      if (sum(is.element(FUNX,xfuns))==0)
          GCM <- annual(gcm,FUN=FUNX,nmin=nmin) else
          eval(parse(text=paste('GCM <- annual(',FUNX,'(gcm),FUN="mean",nmin=nmin)',sep="")))
    } else if (inherits(y,'month')) {
      if (length(table(month(y)))==1) 
        GCM <- subset(gcm,it=month.abb[month(y)[1]]) 
      else
        GCM <- gcm
      if (!is.null(FUNX)) {
        GCM <- do.call(FUNX,list(GCM))
      }
    }
    
    if (verbose) print('Estimate commne EOFs - combine fields')	
    if (is.null(src(T2M))) attr(T2M,'source') <- 'reanalysis'
    T2MGCM <- combine(T2M,GCM)
    if (verbose) print("- - - > EOFs")
    Z <- try(EOF(T2MGCM,verbose=verbose))
    
    ## The test lines are included to assess for non-stationarity
    if (non.stationarity.check) {
      testGCM <- subset(GCM,it=range(year(T2M))) # REB 29.04.2014
      testy <- as.station(regrid(testGCM,is=y))  # REB 29.04.2014
      attr(testGCM,'source') <- 'testGCM'        # REB 29.04.2014
      testZ <- combine(testGCM,GCM)              # REB 29.04.2014
      rm("testGCM"); gc(reset=TRUE)
    }
    rm("gcm","GCM"); gc(reset=TRUE)

    if (verbose) print("- - - > DS (pca)")
    Z0 <- Z
    if (verbose) print(class(attr(Z,'appendix.1')))
    if (biascorrect) Z <- biasfix(Z)
    ds <- try(DS(y,Z,ip=ip,rmtrend=rmtrend,verbose=verbose))
    if(inherits(ds,"try-error")) {
      print(paste("esd failed for",gcmnm.i))
    } else {
      if (verbose) print("post-processing")
      gcmnm[i] <- gcmnm.i
   
      ## Keep the results for the projections:
      if (verbose) print('Extract the downscaled projection')
      z <- attr(ds,'appendix.1') ## KMP 09.08.2015
      
      ## REB: 2016-11-29
      if (test) {
        ## model takes up too much space! can it be stored more efficiently?
        ## REB 2016-11-29: remove most of the contents and keep only a small part
      if (verbose) print('Add reduced model information')
        for (iii in 1:dim(ds)[2]) {
          print(names(attr(ds,'model')[[iii]]))
          attr(ds,'model')[[iii]]$residuals <- NULL
          attr(ds,'model')[[iii]]$effects <- NULL
          attr(ds,'model')[[iii]]$rank <- NULL
          attr(ds,'model')[[iii]]$fitted.values <- NULL
          attr(ds,'model')[[iii]]$assign <- NULL
          attr(ds,'model')[[iii]]$qr <- NULL
          attr(ds,'model')[[iii]]$df.residual <- NULL
          attr(ds,'model')[[iii]]$xlevels <- NULL
          attr(ds,'model')[[iii]]$model <- NULL
          attr(ds,'model')[[iii]]$terms <- NULL
          print(names(attr(ds,'model')[[iii]]))
        }
      }
      
      attr(z,'predictor.pattern') <- attr(ds,'predictor.pattern')
      attr(z,'evaluation') <- attr(ds,'evaluation')

      # REB 2016-11-28: adjust results to have same mean as observations in overlapping period:
      if (verbose) print('adjust offset of predicted PCs for overlapping period')
      index(z) <- year(z)
      zolp <- window(zoo(z),start=start(y),end=end(y))
      coredata(z) <- t(t(coredata(z)) - mean(coredata(zolp)) + colMeans(coredata(y)))
      ## y is a pca with no missing values; z has no NAs.
      if (verbose) print(round(colMeans(y),2))             
      
      cl <- paste('dse.pca$i',i,'_',gsub('-','.',gcmnm[i]),' <- z',sep='')
      eval(parse(text=cl))
      if (verbose) {
        print('Test to see if as.station has all information needed')
        test.stations.ds <- as.station(ds)
        a <- attrcp(y,z);  class(a) <- c("ds",class(y))
        test.stations.z <- as.station(a)
      }
       
      # Diagnose the residual: ACF, pdf, trend. These will together with the
      # cross-validation and the common EOF diagnostics provide a set of
      # quality indicators.

      ## REB 2016-10-20 revised code
      if (verbose) print('examine residuals...')
      cal <- attr(ds,"original_data")
      fit <- attr(ds,"fitted_values")
      res <- cal - fit 
      #REBres <- as.residual(ds)
      res.trend <- 10*diff(range(trend(res)))/diff(range(year(res)))
      ks <- round(ks.test(coredata(res),pnorm)$p.value,4)
#      ar <- as.numeric(acf(trend(cal-fit,result="residual"),
#                           plot=FALSE)[[1]][2]) 
      ar <- as.numeric(acf(coredata(trend(res,result="residual")[,1]),
                         plot=FALSE)[[1]][2])
      if (verbose) print(paste("Residual trend=",
                             round(res.trend,3),'D/decade; K.S. p-val',
                             round(ks,2),'; AR(1)=',round(ar,2)))

      # Evaluation: here are lots of different aspects...
      # Get the diagnostics: this is based on the analysis of common EOFs...

      xval <- attr(ds,'evaluation')
      r.xval <- round(cor(xval[,1],xval[,2]),3)
      if (verbose) print(paste("x-validation r=",r.xval))
      ds.ratio <- round(sd(ds[,1],na.rm=TRUE)/sd(y[,1],na.rm=TRUE),4)
      
      if (verbose) print(paste("sd ratio=",ds.ratio))
      if (verbose) print(names(attributes(ds)))
      if (biascorrect) {
        if (verbose) print('biascorrect')
        diag <- attr(ds,'diagnose')
        if ( (verbose) & !is.null(diag)) str(diag)
      } else diag <- NULL
      
      # diagnose for ds-objects
      if (verbose) print('...')
      #
      if (is.null(diag)) {
        if (verbose) print('no diag')
        ##diag <- diagnose(ds,plot=FALSE)
        scorestats[i,] <- c(1-r.xval,NA,NA,NA,res.trend,ks,ar,1-ds.ratio,
        1-round(var(xval[,2])/var(xval[,1]),2))
        mdiff <- (mean(subset(y,it=range(year(ds))),na.rm=TRUE)-
                  mean(subset(ds,it=range(year(y))),na.rm=TRUE))/
                    sd(y,na.rm=TRUE)
        srati <- sd(subset(ds,it=range(year(y))),na.rm=TRUE)/
                 sd(subset(y,it=range(year(ds))),na.rm=TRUE)
        arati <- ar1(ds)/ar1(y)
      } else {
        if (verbose) print('diag ok')
     # Extract the mean score for leading EOF from the 4 seasons:
        mdiff <- mean(c(diag$s.1$mean.diff[1]/diag$s.1$sd0[1],
                        diag$s.2$mean.diff[1]/diag$s.2$sd0[1],
                        diag$s.3$mean.diff[1]/diag$s.3$sd0[1],
                        diag$s.4$mean.diff[1]/diag$s.4$sd0[1]))
        srati <- mean(1 - c(diag$s.1$sd.ratio[1],diag$s.2$sd.ratio[1],
                            diag$s.3$sd.ratio[1],diag$s.4$sd.ratio[1]))
        arati <- mean(1 - c(diag$s.1$autocorr.ratio[1],
                            diag$s.2$autocorr.ratio[1],
                            diag$s.3$autocorr.ratio[1],
                            diag$s.4$autocorr.ratio[1]))
        
      }
      scorestats[i,] <- c(1-r.xval,mdiff,srati,arati,res.trend,ks,ar,1-ds.ratio,
                          1-round(var(xval[,2])/var(xval[,1]),2))
      if (verbose) print('scorestats')
      if (verbose) print(scorestats[i,])
      quality <- 100*(1-mean(abs(scorestats[i,]),na.rm=TRUE))
      R2 <- round(100*sd(xval[,2])/sd(xval[,1]),2)
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(100*r.xval,2),
                  "R2=",R2,'% ','Common EOF: bias=',round(mdiff,2),
                  ' sd1/sd2=',round(srati,3),
                  "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',
                  round(quality)))

      index(y) <- year(y); index(z) <- year(z)
      if (plot) {
        qcol <- quality
        qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
        cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)
        lines(z[,1],lwd=2,col=cols[qcol])
        lines(y[,1],lwd=3,main='PC1')
      }
    }
   
    if (verbose) print('Downscaling finished')
  }

  ## Unpacking the information tangled up in GCMs, PCs and stations:
  ## Save GCM-wise in the form of PCAs
  gcmnm <- gsub('-','.',gcmnm)

  #Z <- attrcp(y,Z)
  if (verbose) print('Set attributes')
  if (test) {
    attr(dse.pca,'model') <- attr(ds,'model') ## KMP 09-08-2015
    attr(dse.pca,'ceof0') <- Z0
    attr(dse.pca,'ceof') <- Z
  }
  attr(dse.pca,'predictor') <- attr(T2M,'source')
  attr(dse.pca,"longname") <- attr(y,"longname")
  attr(dse.pca,'domain') <- list(lon=lon,lat=lat)
  attr(dse.pca,'scorestats') <- scorestats
  attr(dse.pca,'path') <- path
  attr(dse.pca,'scenario') <- rcp
  attr(dse.pca,'model_id') <- gcmnm
  attr(dse.pca,'variable') <- attr(y,"variable")[1]
  attr(dse.pca,'unit') <- attr(y,"unit")[1]
  attr(dse.pca,'history') <- history.stamp(y)
  # KMP 2018-11-12: difference.z is not defined in this function
  #if (non.stationarity.check) {
  #  attr(dse.pca,'on.stationarity.check') <- difference.z
  #} else {
  #  attr(dse.pca,'on.stationarity.check') <- NULL
  #}
  class(dse.pca) <- c("dsensemble","pca","list")

  if(!is.null(path.ds)) file.ds <- file.path(path.ds,file.ds)
  if (verbose) print(file.ds)
  save(file=file.ds,dse.pca)
  if (verbose) print("---")
  invisible(dse.pca)
}

#' @export DSensemble.eof
DSensemble.eof <- function(y,...,plot=TRUE,path="CMIP5.monthly",rcp="rcp45",biascorrect=FALSE,
                           predictor="ERA40_slp_mon.nc",non.stationarity.check=FALSE,
                           ip=1:5,lon=c(-30,20),lat=c(-20,10),it=NULL,rel.cord=TRUE,nmin=NULL,
                           lev=NULL,levgcm=NULL,select=NULL,FUN="mean",rmtrend=TRUE,FUNX="mean",
                           xfuns='C.C.eq',threshold=1,pattern="psl_Amon_ens_",verbose=FALSE,
                           file.ds="DSensemble.eof.rda",path.ds=NULL,ds.1900.2099=TRUE,test=FALSE) {

  if(verbose) print("DSensemble.eof")
  stopifnot(inherits(y,c("EOF","field")))
  cls <- class(y)

  if (!is.na(attr(y,'longitude'))[1] & (any(lon>0) & any(lon<0)) & rel.cord)
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude'))[1] & (any(lat>0) & any(lat<0)) & rel.cord)
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )

  if (sum(!is.finite(lon))>0) 
    warning(paste('Bad longitude range provided: ',paste(lon,collapse='-')))
  if (sum(!is.finite(lat))>0) 
    warning(paste('Bad latitude range provided: ',paste(lat,collapse='-')))

  if (is.character(predictor))
    slp <- retrieve(file=predictor,lon=lon,lat=lat,lev=lev,
                    verbose=verbose) else
  if (inherits(predictor,'field'))
    slp <- subset(predictor,is=list(lon=lon,lat=lat))

  if (inherits(y,'season')) {
    if (verbose) print('seasonal data')
    if (FUNX!='C.C.eq') SLP <- as.4seasons(slp,FUN=FUNX) else
                        eval(parse(text=paste('SLP <- as.4seasons(',FUNX,'(slp),FUN="mean",nmin=nmin)',sep="")))
    SLP <- matchdate(SLP,y)

    if (length(table(season(y)))>1) {
      if (verbose) print('--- Apply DS to seasons seperately ---')
      Z <- list(info=paste('DSensemble.pca for different seasons: ',
                                paste(lon,collapse='-'),'E/',paste(lat,collapse='-'),'N',sep=''))
      ## KMP 2016-10-25: Looping over seasons will not work if y is an eof object.
      ##   I added a temporary fix, turning the multi-season eof object into a field
      ##   before selecting a season. This could result in a serious loss of information.
      if(inherits(y,"eof")) y <- as.field(y)
      for (season in names(table(season(y)))) {
        if (verbose) print(paste('Select',season))
        z <- DSensemble.eof(subset(y,it=season),plot=plot,path=path,
                            rcp=rcp,biascorrect=biascorrect,predictor=SLP,
                            non.stationarity.check=non.stationarity.check,
                            ip=ip,lon=lon,lat=lat,rel.cord=FALSE,
                            select=select,FUN=FUN,rmtrend=rmtrend,
                            FUNX=FUNX,threshold=threshold,
                            pattern=pattern,verbose=verbose,nmin=nmin)
        eval(parse(text=paste('Z$',season,' <- z',sep='')))
      }
      if (verbose) print('--- Results returned as a list ---')
      return(Z)
    }

  } else if (inherits(y,'annual')) {
    if (verbose) print('annual data')
    if (!is.null(it)) {
      if (verbose) print('Extract some months or a time period')
      if (verbose) print(it)
      slp <- subset(slp,it=it)
      if (is.null(nmin) & is.character(it)) nmin <- length(it)
      if (is.null(nmin)) nmin <- 1
    }
    if (!is.annual(slp)) {
      if (FUNX!='C.C.eq')
        SLP <- annual(slp,FUN=FUNX,nmin=nmin) else
        SLP <- annual(C.C.eq(slp),FUN='mean',nmin=nmin)
  } else SLP <- slp
    ## Synchronise
    if (verbose) {print(index(SLP)); print(index(y))}
    SLP <- matchdate(SLP,y)
  } else if (inherits(y,'month')) {
    SLP <- matchdate(slp,y)
  }
  if (inherits(SLP,"eof")) SLP <- as.field(SLP)
  rm("predictor","slp"); gc(reset=TRUE)

  if(!inherits(y,"eof")) y <- EOF(y)
      
  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.names=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) select <- 1:N else
                       N <- length(select)
  if (verbose) {print('GCMs:'); print(path); print(ncfiles[select])}

  d.y <- dim(y)
  years <- 1900:2100
  m <- length(years)
  months <- rep(month(y)[1],m)
  X <- matrix(rep(NA,N*m*d.y[2]),N,m*d.y[2])
  dim(X) <- c(d.y[2],N,m)
  gcmnm <- rep("",N)
  scorestats <- matrix(rep(NA,N*9),N,9)
  colnames(scorestats) <- c("1-r.xval","mean.diff","sd.ratio","autocorr.diff",
                            "res.trend","res.K-S","res.ar1",'amplitude.ratio',
                            "1-R2")

  t <- as.Date(paste(years,months,'01',sep='-'))

  flog <- file("DSensemble.eof-log.txt","at")

  ## Set up a list environment to keep all the results
  dse.eof <- list(info=paste('DSensemble.pca for different seasons: ',
                             paste(lon,collapse='-'),'E/',paste(lat,collapse='-'),'N',sep=''),eof=y) 
  if (verbose) print("loop...") 
  for (i in 1:N) {
    if (verbose) print(ncfiles[select[i]])
    gcm <- try(retrieve(file = ncfiles[select[i]],
                        lon=range(lon(SLP))+c(-2,2),
                        lat=range(lat(SLP))+c(-2,2),
                        lev=levgcm,verbose=verbose))
    if(inherits(gcm,"try-error")) {
      print(paste("retrieve failed for",ncfiles[select[i]]))
    } else {
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099)) else
                      gcm <- subset(gcm,it=c(min(year(y),na.rm=TRUE),2099))
    if (length(index(gcm))<=1) print(paste('Problem selecting GCM results in period',
                                           min(year(y),na.rm=TRUE),'2099'))
    ## KMP 2016-08-09 added separate level input for slp and gcm
    ##                because they can have levels of different units
    if(is.null(levgcm) & !is.null(attr(gcm,"level")))
      levgcm <- attr(gcm,"level")
    if (!is.null(it)) {
      if ((nmin==0) & is.character(it)) nmin <- length(it) 
      if (verbose) print('Extract some months or a time period')
      if (verbose) {print(it); print(nmin)}
      gcm <- subset(gcm,it=it)
    }
    #gcmnm[i] <- attr(gcm,'model_id')
    #gcmnm.i <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-r")
    # KMP: 10.03.2017 - pass on additional information about GCM runs (gcm + rip - realization, initialization, physics version)
    #gcmnm.i <- paste(attr(gcm,'model_id'),attr(gcm,'parent_experiment_rip'),sep="-")   
    rip <- NULL
    if(any(grepl("rip", names(attributes(gcm))))) {
      nm.rip <- names(attributes(gcm))[grepl("rip",names(attributes(gcm)))][[1]]
      if(!is.null(attr(gcm, nm.rip))) {
        rip <- attr(gcm, nm.rip)
        if(!grepl("r[0-9]{1,2}i[0-9]{1,2}p[0-9]{1,2}", rip)) rip <- NULL
      }
    }
    if(is.null(rip)) {
      nm.r <- names(attributes(gcm))[grep("realization",names(attributes(gcm)))][[1]]
      nm.i <- names(attributes(gcm))[grep("initialization",names(attributes(gcm)))][[1]]
      nm.p <- names(attributes(gcm))[grep("physics",names(attributes(gcm)))][[1]]
      rip <- paste0("r",attr(gcm,nm.r),"i",attr(gcm,nm.i),"p",attr(gcm,nm.p))
    }
    gcmnm.i <- paste0(attr(gcm,'model_id'),".",rip)
    if (verbose) print(class(y))

    ## REB 2016-10-25
    if (inherits(y,'season')) {
      if (verbose) print(paste('Seasonally aggregated',FUNX,'for GCM'))
      if (sum(is.element(FUNX,xfuns))==0) {
          if (verbose) print('No special transformation (EOF)') 
          GCM <- as.4seasons(gcm,FUN=FUNX,nmin=nmin)
       } else {
           if (verbose) print('Need to aggregate FUNX(gcm)')
           eval(parse(text=paste('GCM <- as.4seasons(',FUNX,'(gcm),FUN="mean",nmin=nmin)',sep="")))
       }
      if (verbose) {print(season(SLP)[1]); print(dim(GCM))}
      GCM <- subset(GCM,it=season(SLP)[1])
    } else if (inherits(y,'annual')) {
      if (verbose) print(paste('Annualy aggregated',FUNX,'for GCM'))
      if (sum(is.element(FUNX,xfuns))==0)
          GCM <- annual(gcm,FUN=FUNX,nmin=nmin) else
          eval(parse(text=paste('GCM <- annual(',FUNX,'(gcm),FUN="mean",nmin=nmin)',sep="")))
    } else if (inherits(y,'month')) {
      if (length(table(month(y)))==1) 
        GCM <- subset(gcm,it=month.abb[month(y)[1]]) 
      else
        GCM <- gcm
      if (!is.null(FUNX)) {
        GCM <- do.call(FUNX,list(GCM))
      }
    }
    
    #ds.parts <- TRUE
    #if(ds.parts) {
    #  GCM <- subset(GCM,it=c(seq(1900,1920),year(SLP),seq(1981,2100)))
    #}
      
    ## REB 2016-10-25
    #if (inherits(y,'season')) {
    #  if (FUNX!='C.C.eq') GCM <- as.4seasons(gcm,FUN=FUNX,nmin=nmin) else
    #                      GCM <- as.4seasons(C.C.eq(gcm),FUN='mean',nmin=nmin)
    #  GCM <- subset(GCM,it=season(SLP)[1])
    #} else if (inherits(y,'annual')) {
    #  if (FUNX!='C.C.eq') GCM <- annual(gcm,FUN=FUNX,nmin=nmin) else
    #                      GCM <- annual(C.C.eq(gcm),FUN='mean',nmin=nmin)
    #} else if (inherits(y,'month')) {
    #  if (length(table(month(y)))==1)
    #    GCM <- subset(gcm,it=month.abb[month(y)[1]]) else
    #    GCM <- gcm
    #}
    
    if (verbose) {str(SLP); str(GCM)}
    SLPGCM <- combine(SLP,GCM)
    if (verbose) print("- - - > EOFs")
    
    Z <- try(EOF(SLPGCM))
    
    ## The test lines are included to assess for non-stationarity
    ## KMP 2018-11-02: Does the nonstationarity test work for DSensemble.eof?
    if (non.stationarity.check) {
      testGCM <- subset(GCM,it=range(year(SLP)))      
      testy <- regrid(testGCM,is=as.field(y)) 
      attr(testGCM,'source') <- 'testGCM'
      testZ <- combine(testGCM,GCM)
      difference.z <- testy - testZ
      rm("testGCM"); gc(reset=TRUE)
    }
    rm("gcm","GCM"); gc(reset=TRUE)

    if (verbose) print("- - - > DS (eof)")
    Z0 <- Z
    if (verbose) print(class(attr(Z,'appendix.1')))
    if (biascorrect) Z <- biasfix(Z)
    
    diag <- diagnose(Z)
    
    ds <- try(DS(y,Z,ip=ip,verbose=verbose))
    if(inherits(ds,"try-error")) {
      print(paste("esd failed for",gcmnm.i))
    } else {
      if (verbose) print("post-processing")
      gcmnm[i] <- gcmnm.i
      ## Unpacking the information tangled up in GCMs, PCs and stations:
      ## Save GCM-wise in the form of PCAs
      gcmnm[i] <- gsub('-','.',gcmnm[i])
      gcmnm[i] <- gsub('/','.',gcmnm[i])
      
      ## Keep the results for the projections:
      if (verbose) print('Extract the downscaled projection')
      z <- attr(ds,'appendix.1') ## KMP 09.08.2015
      ##attr(z,'model') <- attr(ds,'model') ## KMP 09-08-2015
      ## model takes up too much space! can it be stored more efficiently?
      
      ## REB: 2016-11-29
      if (test) {
        ## model takes up too much space! can it be stored more efficiently?
        ## REB 2016-11-29: remove most of the contents and keep only a small part
        if (verbose) print('Reduced model information')
        for (iii in 1:dim(ds)[2]) {
          print(names(attr(ds,'model')[[iii]]))
          attr(ds,'model')[[iii]]$residuals <- NULL
          attr(ds,'model')[[iii]]$effects <- NULL
          attr(ds,'model')[[iii]]$rank <- NULL
          attr(ds,'model')[[iii]]$fitted.values <- NULL
          attr(ds,'model')[[iii]]$assign <- NULL
          attr(ds,'model')[[iii]]$qr <- NULL
          attr(ds,'model')[[iii]]$df.residual <- NULL
          attr(ds,'model')[[iii]]$xlevels <- NULL
          attr(ds,'model')[[iii]]$model <- NULL
          attr(ds,'model')[[iii]]$terms <- NULL
          print(names(attr(ds,'model')[[iii]]))
        }
      }
      
      attr(z,'predictor.pattern') <- attr(ds,'predictor.pattern')
      attr(z,'evaluation') <- attr(ds,'evaluation')
    
      ## Store the results in a list element
      cl <- paste('dse.eof$i',i,'_',gcmnm[i],' <- z',sep='')
      eval(parse(text=cl))
      if (verbose) {
        print('Test to see if as.field has all information needed')
        test.field.ds <- as.field(ds)
        a <- attrcp(y,z);  class(a) <- class(y)
        test.field.z <- as.field(a)
      }
       
      # Diagnose the residual: ACF, pdf, trend. These will together with the
      # cross-validation and the common EOF diagnostics provide a set of
      # quality indicators.
      ## REB 2016-10-20
      cal <- attr(ds,"original_data")
      fit <- attr(ds,"fitted_values")
      if (verbose) print('examine residuals...')
      res <- cal - fit ## REB 2016-10-20 test PCs
      res.trend <- 10*diff(range(trend(res)))/diff(range(year(res)))
      ks <- round(ks.test(coredata(res),pnorm)$p.value,4)
      ar <- as.numeric(acf(coredata(trend(res,result="residual")[,1],
                         plot=FALSE))[[1]][2])
      ## REB 2016-10-20 - end of revised code

      if (verbose) print(paste("Residual trend=",
                             round(res.trend,3),'D/decade; K.S. p-val',
                             round(ks,2),'; AR(1)=',round(ar,2)))

      # Evaluation: here are lots of different aspects...
      # Get the diagnostics: this is based on the analysis of common EOFs...

      xval <- attr(ds,'evaluation')
      r.xval <- round(sapply(seq(1,2*ncol(ds),2),function(i) cor(xval[,i],xval[,i+1])),3)
                                        #round(cor(xval[,1],xval[,2]),3)
      if (verbose) print(paste("x-validation r=",paste(r.xval,collapse=", ")))
      ds.ratio <- round(sapply(seq(1,ncol(ds)),
                               function(i) sd(ds[,i],na.rm=TRUE)/sd(y[,i],na.rm=TRUE)),3)
                                        #round(sd(ds[,1],na.rm=TRUE)/sd(y[,1],na.rm=TRUE),4)
      
      if (verbose) print(paste("sd ratio=",paste(ds.ratio,collapse=", ")))
      if (verbose) print(names(attributes(ds)))
      if (biascorrect) {
        if (verbose) print('biascorrect')
        diag <- attr(ds,'diagnose')
        if ( (verbose) & !is.null(diag)) str(diag)
      } else diag <- NULL
      
      # diagnose for ds-objects
      if (verbose) print('...')
      #
      if (is.null(diag)) {
        if (verbose) print('no diag')
        mdiff <- (mean(subset(y,it=range(year(y))),na.rm=TRUE)-
                  mean(subset(ds,it=range(year(ds))),na.rm=TRUE))/
                    sd(y,na.rm=TRUE)
        srati <- 1 - sd(subset(ds,it=range(year(ds))),na.rm=TRUE)/
                    sd(subset(y,it=range(year(y))),na.rm=TRUE)
        adiff <- as.numeric(acf(coredata(y)[,1],plot=FALSE,na.action=na.pass)[[1]][2])-
                 as.numeric(acf(coredata(ds)[,1],plot=FALSE,na.action=na.pass)[[1]][2])
       } else {
        if (verbose) print('diag ok')
     # Extract the mean score for leading EOF from the 4 seasons:
        mdiff <- mean(c(diag$s.1$mean.diff[1]/diag$s.1$sd0[1],
                        diag$s.2$mean.diff[1]/diag$s.2$sd0[1],
                        diag$s.3$mean.diff[1]/diag$s.3$sd0[1],
                        diag$s.4$mean.diff[1]/diag$s.4$sd0[1]))
        srati <- mean(1 - c(diag$s.1$sd.ratio[1],diag$s.2$sd.ratio[1],
                            diag$s.3$sd.ratio[1],diag$s.4$sd.ratio[1]))
        adiff <- mean(1 - c(diag$s.1$autocorr.ratio[1],
                            diag$s.2$autocorr.ratio[1],
                            diag$s.3$autocorr.ratio[1],
                            diag$s.4$autocorr.ratio[1]))
        
      }
      scorestats[i,] <- c(1-r.xval[1],mdiff,srati,adiff,res.trend,ks,ar,1-ds.ratio[1],
                          1-round(var(xval[,2])/var(xval[,1]),2))
      if (verbose) print('scorestats')
      if (verbose) print(scorestats[i,])
      quality <- 100*(1-mean(abs(scorestats[i,]),na.rm=TRUE))
      R2 <- round(100*sd(xval[,2])/sd(xval[,1]),2)
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(100*r.xval[1],2),
                  "R2=",R2,'% ','Common EOF: bias=',round(mdiff,2),
                  ' sd1/sd2=',round(srati,3),
                  "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',
                  round(quality)))
    }
    }
    
    if (verbose) print('Downscaling finished')

  }

  #Z <- attrcp(y,Z)
  if (verbose) print('Set attributes')
  if (test) {
    attr(dse.eof,'model') <- attr(ds,'model') ## KMP 09-08-2015
    attr(dse.eof,'ceof0') <- Z0
    attr(dse.eof,'ceof') <- Z
  }
  attr(dse.eof,'predictor') <- attr(SLP,'source')
  attr(dse.eof,'domain') <- list(lon=lon,lat=lat)
  attr(dse.eof,'scorestats') <- scorestats
  attr(dse.eof,'path') <- path
  attr(dse.eof,'scenario') <- rcp
  attr(dse.eof,'variable') <- attr(y,"variable")[1]
  attr(dse.eof,'unit') <- attr(y,"unit")[1]
  attr(dse.eof,'unitarea') <- attr(y,"unitarea")
  attr(dse.eof,'history') <- history.stamp(y)
  if (non.stationarity.check)
    attr(dse.eof,'on.stationarity.check') <- difference.z else
    attr(dse.eof,'on.stationarity.check') <- NULL
  class(dse.eof) <- c("dsensemble","eof","list")

  if(!is.null(path.ds)) file.ds <- file.path(path.ds,file.ds)
  save(file=file.ds,dse.eof)
  if (verbose) print("---")
  invisible(dse.eof)
}

#' @export DSensemble.field
DSensemble.field <- function(y,...,plot=TRUE,path="CMIP5.monthly/",rcp="rcp45",biascorrect=FALSE,
                           predictor="ERA40_t2m_mon.nc",non.stationarity.check=FALSE,
                           ip=1:16,lon=c(-30,20),lat=c(-20,10),it=c('djf','mam','jja','son'),
                           rel.cord=TRUE,select=NULL,FUN="mean",rmtrend=TRUE,FUNX="mean",
                           xfuns='C.C.eq',threshold=1,pattern="tas_Amon_ens_",verbose=FALSE,
                           file.ds="DSensemble.rda",path.ds=NULL,nmin=NULL,ds.1900.2099=TRUE) {
  ## For downscaling gridded predictand. This is a wrap-around which extracts the season or aggregates
  ## to annual values and then calls the other types for the downscaling.
  ## KMP 2016-10-25: Redirect to DSensemble.eof
  if(verbose) print("DSensemble.field")
  dse.eof <- DSensemble.eof(y,plot=plot,path=path,rcp=rcp,biascorrect=biascorrect,
                           predictor=predictor,non.stationarity.check=non.stationarity.check,
                           ip=ip,lon=lon,lat=lat,it=it,rel.cord=rel.cord,select=select,
                           FUN=FUN,rmtrend=rmtrend,FUNX=FUNX,xfuns=xfuns,threshold=threshold,
                           pattern=pattern,verbose=verbose,
                           file.ds=file.ds,path.ds=path.ds,nmin=nmin)
  invisible(dse.eof)
}

#' @export DSensemble.station
DSensemble.station <- function(y,...,verbose=FALSE) {
  if(verbose) print("DSensemble.station")
  dse <- DSensemble.default(y=y,...,verbose=verbose)
  return(dse)
}
