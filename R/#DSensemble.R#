# R.E. Benestad - replacement for old ds.one to downscale the CMIP 5 ensemble
# seasonal mean and standard deviation fortemperature
#

ar1 <- function(x,...) acf(x,plot=FALSE,na.action = na.pass)$acf[2]

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


DSensemble<-function(y,...) UseMethod("DSensemble")

DSensemble.default <- function(y,path='CMIP5.monthly/',rcp='rcp45',...) {
   ## 
  stopifnot(!missing(y),inherits(y,"station"),
            file.exists(paste(file.path(path,rcp,fsep = .Platform$file.sep))))
  
  if (is.null(attr(y,'aspect'))) attr(y,'aspect') <- "original"
  
  if (inherits(y,'annual'))
    z <- DSensemble.annual(y,path=path,rcp=rcp,threshold=1,...) else
  if (is.T(y))
    z <- DSensemble.t2m(y,path=path,rcp=rcp,...) else
  if (is.precip(y))
    z <- DSensemble.precip(y,path=path,rcp=rcp,threshold=1,...) else
    z <- NULL 
    
  return(z)
}

DSensemble.t2m <- function(y,plot=TRUE,path="CMIP5.monthly/",
                           predictor="ERA40_t2m_mon.nc",
                           rcp="rcp45",biascorrect=FALSE,
                           non.stationarity.check=FALSE,type='ncdf4',
                           ip=1:6,lon=c(-20,20),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                           select=NULL,FUN="mean",FUNX="mean",xfuns='C.C.eq',
                           pattern="tas_Amon_ens_",
                           path.ds=NULL,file.ds="DSensemble.rda",
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

  ylim <- c(0,0)
  ylim <- switch(FUN,'mean'=c(-2,8),'sd'=c(-0.5,1),'ar1'=c(-0.5,0.7))
  if (verbose) print(paste('set ylim based on "',FUN,'" -> c(',ylim[1],', ',ylim[2],')',sep=''))
  
  if (plot) {
    if(verbose) print("Plot station data (predictand)")
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
    t2m <- retrieve(ncfile=predictor,lon=lon,lat=lat,
                    type=type,verbose=verbose) else
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
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
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

  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)

  if(verbose) print("Quick test")  
  #load('control.ds.1.rda'); T2M.ctl -> T2M
  flog <- file("DSensemble.t2m-log.txt","at")
  
  if (verbose) print("loop...")
  for (i in 1:N) {
    if (verbose) print(ncfiles[select[i]])
    gcm <- retrieve(ncfile = ncfiles[select[i]],type=type,
                          lon=range(lon(T2M))+c(-2,2),
                          lat=range(lat(T2M))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099))
    #gcmnm[i] <- attr(gcm,'model_id'))
    gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
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
    #save(file='inside.dsens.rda',T2M,GCM)
    #rm("GCM"); gc(reset=TRUE)
    # The test lines are included to assess for non-stationarity
    if (non.stationarity.check) {
      testGCM <- subset(GCM,it=range(year(T2M))) # REB 29.04.2014
      testy <- as.station(regrid(testGCM,is=y))  # REB 29.04.2014
      attr(testGCM,'source') <- 'testGCM'        # REB 29.04.2014
      testZ <- combine(testGCM,GCM)              # REB 29.04.2014
      rm("testGCM"); gc(reset=TRUE)
    }
    ##
    # REB: 30.04.2014 - new lines...
    if (verbose) print("- - - > DS")
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
    #save(file='inside.dsens.1.rda',ds,y,Z)

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
      qcol <- quality
      qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
     
      if (plot) {
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
  rm("GCM")

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

DSensemble.precip <- function(y,plot=TRUE,path="CMIP5.monthly/",
                              rcp="rcp45",biascorrect=FALSE,
                              predictor="ERA40_pr_mon.nc",
                              non.stationarity.check=FALSE,
                              type='ncdf4',
                              ip=1:6,lon=c(-10,10),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                              select=NULL,FUN="wetmean",
                              FUNX="sum",xfuns='C.C.eq',threshold=1,
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
    pre <- retrieve(ncfile=predictor,lon=lon,lat=lat,
                    type=type,verbose=verbose) else
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
  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)
  unit <- attr(y,'unit')

  ylim <- switch(FUN,
                 'exceedance'=c(0,10),'wetmean'=c(0,10),
                 'wetfreq'=c(0,0),'spell'=c(0,0),
                 'mean'=c(-10,50),'sd'=c(-5,10),'ar1'=c(-0.5,0.7),
                 'HDD'=c(0,5000),'CDD'=c(0,500),'GDD'=c(0,2000))
  if (is.null(ylim)) ylim <- c(0,0)
  
  if (plot) {
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
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
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
    gcm <- retrieve(ncfile = ncfiles[select[i]],type=type,
                    lon=range(lon(PRE))+c(-2,2),lat=range(lat(PRE))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099))
    gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
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
    if (verbose) print("- - - > DS")
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
      qcol <- quality
      qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100

      index(z) <- year(z); index(ds) <- year(ds)
      if (plot) {
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

DSensemble.annual <- function(y,plot=TRUE,path="CMIP5.monthly/",
                              rcp="rcp45",biascorrect=FALSE,
                              predictor="ERA40_t2m_mon.nc",
                              non.stationarity.check=FALSE,type='ncdf4',
                              ip=1:6,lon=c(-10,10),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                              abscoords=FALSE,select=NULL,FUN=NULL,
                              FUNX="mean",xfuns='C.C.eq',threshold=1,
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
    pre <- retrieve(ncfile=predictor,lon=lon,lat=lat,
                    type=type,verbose=verbose) else
  if (inherits(predictor,'field')) pre <- predictor
  rm("predictor"); gc(reset=TRUE)
  attr(pre,"source") <- "ERA40"

  if (!is.null(it)) {
    if (verbose) print('Extract some months or a time period')
    if (verbose) print(it)
    pre <- subset(pre,it=it)
      ## if it is character, then then extraction of months reduces number of
      ## months per year.
    if ((is.null(nmin)) & (is.character(it))) nmin <- length(it)
  }

  # Use proportional variations
  if (verbose) print("Annual mean")
  if (sum(is.element(FUNX,xfuns))==0) PRE <- annual(pre,FUN=FUNX,nmin=nmin) else
  eval(parse(text=paste('PRE <- annual(',FUNX,'(pre),FUN="mean",nmin=nmin)',sep="")))
  if (verbose) print("graphics")
  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)
  unit <- attr(y,'unit')
  ylim <- c(0,10)

  if (plot) {
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
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
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
    gcm <- retrieve(ncfile = ncfiles[select[i]],type=type,
                    lon=range(lon(PRE))+c(-2,2),lat=range(lat(PRE))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099))
    gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realisation'),sep="-")
    #gcmnm[i] <- attr(gcm,'model_id')
    if (verbose) print(varid(gcm))

    if (!is.null(it)) {
      if (verbose) print('Extract some months or a time period')
      if (verbose) print(it)
      gcm <- subset(gcm,it=it)
    }
    
    if (sum(is.element(FUNX,xfuns))==0) GCM <- annual(gcm,FUN=FUNX,nmin=nmin) else
    eval(parse(text=paste('GCM <- annual(',FUNX,'(gcm),FUN="mean",nmin=nmin)',sep="")))

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
    if (verbose) print("- - - > DS")
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
      qcol <- quality
      qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100

      index(z) <- year(z); index(ds) <- year(ds)
      if (plot) {
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
  if (non.stationarity.check)
    attr(X,'on.stationarity.check') <- difference.z else
    attr(X,'on.stationarity.check') <- NULL
  attr(X,'history') <- history.stamp(y)
  class(X) <- c("dsensemble","zoo")
  save(file="DSensemble.rda",X)
  print("---")
  invisible(X)
}

DSensemble.season <- function(y,season="djf",plot=TRUE,path="CMIP5.monthly/",
                           predictor="slp.mon.mean.nc",
                           rcp="rcp45",biascorrect=FALSE,
                           non.stationarity.check=FALSE,type='ncdf4',
                           ip=1:6,lon=c(-20,20),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                           select=NULL,FUN="mean",FUNX="mean",xfuns='C.C.eq',
                           pattern="psl_Amon_ens_",
                           path.ds=NULL,file.ds=NULL,
                           nmin=NULL,verbose=FALSE,ds.1900.2099=TRUE) {

  if(verbose) print("DSensemble.season")
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
  ys <- ys[attr(ys,"n.valid")>=nmin]
  if (FUN=="sum" & grepl("month",attr(ys,"unit"))) {
    attr(ys,"unit") <- gsub("month","season",attr(ys,"unit"))
  }

  ylim <- c(0,0)
  ylim <- switch(FUN,'mean'=c(-2,8),'sd'=c(-0.5,1),'ar1'=c(-0.5,0.7),
                 'sum'=c(-6,12))
  if (verbose) print(paste('set ylim based on "',FUN,
                           '" -> c(',ylim[1],', ',ylim[2],')',sep=''))
  
  if (plot) {
    if(verbose) print("Plot station data (predictand)")
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
    slp <- retrieve(ncfile=predictor,lon=lon,lat=lat,
                    type=type,verbose=verbose) else
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
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
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
  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)

  if(verbose) print("Quick test")  
  flog <- file("DSensemble.season-log.txt","at")

  if (verbose) print("loop...") 
  for (i in 1:N) {
    if (verbose) print(ncfiles[select[i]])
    gcm <- retrieve(ncfile = ncfiles[select[i]],type=type,
                          lon=range(lon(SLP))+c(-2,2),
                          lat=range(lat(SLP))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099))
    gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realisation'),sep="-")
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

    if (verbose) print("- - - > DS")
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
      qcol <- quality
      qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
     
      if (plot) {
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
  attr(X,'predictor') <- attr(T2M,'source')
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


DSensemble.mu <- function(y,plot=TRUE,path="CMIP5.monthly/",
                          rcp="rcp45",biascorrect=FALSE,
                          predictor=list(t2m="data/ncep/air.mon.mean.nc",
                                         olr="data/ncep/OLR.mon.mean.nc",
                                         slp="data/ncep/slp.mon.mean.nc"),
                          non.stationarity.check=FALSE,type='ncdf4',
                          ip=1:16,lon=c(-30,20),lat=c(-20,10),it=NULL,rel.cord=TRUE,
                          select=NULL,FUN="wetmean",threshold=1,
                          pattern=c("tas_Amon_ens_","slp_Amon_ens_"),
                          verbose=FALSE,nmin=365,ds.1900.2099=TRUE) {

# This function is for downscaling wet-day mean using a combination of predictors

  
  
  # Get the global mean temeprature: pentads

  # Or a combination of OLR + C.C.eq(t2m) + SLP


    # FUN: exceedance, wetfreq, wet, dry

  if (verbose) print('DSensemble.mu')

  if (verbose) print(paste('The predictor: annual',FUN))
  if (!inherits(y,'annual')) y <- annual(y,FUN=FUN,threshold=threshold,nmin=nmin)
  index(y) <- year(y)
  
  if (!is.na(attr(y,'longitude')) & rel.cord)
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude')) & rel.cord)
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )

  if (sum(!is.finite(lon))>0) 
    warning(paste('Bad longitude range provided: ',paste(lon,collapse='-')))
  if (sum(!is.finite(lat))>0) 
    warning(paste('Bad latitude range provided: ',paste(lat,collapse='-')))
  
  # Get the predictor: NCEP/NCAR
  if (verbose) print(paste("Get the set of predictors:",names(predictor),collapse=' '))
  if (is.character(predictor[[1]])) pre1 <- retrieve(ncfile=predictor[[1]],
                                                     type=type,
                                                     lon=lon,lat=lat,
                                                     verbose=verbose)
  if (is.character(predictor[[2]])) pre2 <- retrieve(ncfile=predictor[[2]],
                                                     type=type,
                                                     lon=lon,lat=lat,
                                                     verbose=verbose)
  if (is.character(predictor[[3]])) pre3 <- retrieve(ncfile=predictor[[3]],
                                                     type=type,
                                                     lon=lon,lat=lat,
                                                     verbose=verbose)
  
  # Combine the predictors
  if (verbose) print("Annual mean - predictors")
  PREX1 <- annual(C.C.eq(pre1),FUN='mean') # Clausius-Claperyron eq. -> sat. vapour pressure
  PREX2 <- annual(pre2,FUN='mean')
  PREX3 <- annual(pre3,FUN='mean')
  
  if (verbose) print("graphics")
  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)
  unit <- attr(y,'unit')
  #ylim <- switch(deparse(substitute(FUN)),

  ylim <- switch(FUN,
                 'exceedance'=c(0,10),'wetmean'=c(0,10),
                 'wetfreq'=c(0,0),'spell'=c(0,0),
                 'mean'=c(-10,50),'sd'=c(-5,10),'ar1'=c(-0.5,0.7),
                 'HDD'=c(0,5000),'CDD'=c(0,500),'GDD'=c(0,2000))
  if (is.null(ylim)) ylim <- c(0,0)
  
  if (plot) {
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
  ncfiles1 <- list.files(path=path,pattern=pattern1,full.name=TRUE)
  ncfiles2 <- list.files(path=path,pattern=pattern2,full.name=TRUE)
  ncfiles3 <- list.files(path=path,pattern=pattern3,full.name=TRUE)
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
  colnames(scorestats) <- c("r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",
                            'amplitude.ration','1-R2')

  flog <- file("DSensemble.precip-log.txt","at")
  dse <- list(description='DSensemble.mu')
  
  for (i in 1:N) {
    ## Need to ensure that the different predictor files match...
    print(paste(i,N,ncfiles1[select[i]],ncfiles2[select[i]],ncfiles3[select[i]]))
    gcm1 <- retrieve(ncfile = ncfiles1[select[i]],type=type,
                    lon=range(lon(PRE1))+c(-2,2),lat=range(lat(PRE1))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm1 <- subset(gcm1,it=c(1900,2099))
    gcm2 <- retrieve(ncfile = ncfiles2[select[i]],type=type,
                    lon=range(lon(PRE2))+c(-2,2),lat=range(lat(PRE2))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm2 <- subset(gcm2,it=c(1900,2099))
    gcm3 <- retrieve(ncfile = ncfiles3[select[i]],type=type,
                    lon=range(lon(PRE3))+c(-2,2),lat=range(lat(PRE3))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm3 <- subset(gcm3,it=c(1900,2099))
    gcmnm[i] <- paste(attr(gcm1,'model_id'),attr(gcm,'realization'),sep="-")
    #gcmnm[i] <- attr(gcm,'model_id')
    if (verbose) print(varid(gcm1))
    
    GCM1 <- annual(C.C.eq(gcm1),FUN='mean')
    GCM2 <- annual(gcm2,FUN='mean')
    GCM3 <- annual(gcm3,FUN='mean')

    model.id <- attr(gcm1,'model_id')
    rm("gcm","GCMX"); gc(reset=TRUE)
    if (verbose) print("combine the three predictors")
    #
    PREGCM1 <- combine(PRE1,GCM1)
    PREGCM2 <- combine(PRE2,GCM2)
    PREGCM3 <- combine(PRE3,GCM3)
    if (verbose) print("EOF")
    Z1 <- EOF(PREGCM1)
    Z2 <- EOF(PREGCM2)
    Z3 <- EOF(PREGCM3)
    
    if (verbose) print("diagnose")
    diag1 <- diagnose(Z1)
    diag2 <- diagnose(Z2)
    diag3 <- diagnose(Z3)

    if (biascorrect) 
      X <- list(Z1=biasfix(Z1),
                Z2=biasfix(Z2),
                Z3=biasfix(Z3)) else
      X < list(Z1=Z1,
               Z2=Z2,
               Z3=Z3)
    x <- zoo(X[[1]])
    np <- length(names(x))
    for (i in 2:np) {
        x <- merge(x,zoo(X[[i]]),all=TRUE)
        w <- c(w,attr(X[[i]],'eigenvalues')/sum(attr(X[[i]],'eigenvalues')))
        id <- c(id,rep(i,length(attr(X[[i]],'eigenvalues'))))
      }
    if (verbose) print(c(dim(x),length(w)))
    t <- index(x)
    ## apply the weights
    x <- x %*% diag(w)
    xm <- rowMeans(x)
    x <- x[is.finite(xm),]; t <- t[is.finite(xm)]

    ## Apply an SVD to the combined PCs to extract the common signal in the
    ## different predictors - these are more likely to contain real physics
    ## and be related to the predictand.
    if (verbose) print('svd')
    udv <- svd(coredata(x))
    if (verbose) print(summary(udv))

    ## If the predictor is a common EOF, then also combine the appended fields
    ## the same way as the original flield.
    if (inherits(X[[1]],'comb')) {
        z <- z %*% diag(w)
        udvz <- svd(coredata(z))
    }

    eof <- zoo(udv$u[,1:20],order.by=t)
    ## Let the pattern contain the weights for the EOFs in the combined
    ## PC matrix, rather than spatial patterns. The spatial patterns are
    ## then reconstructed from these.
    pattern <- matrix(rep(1,length(udv$v[,1:20])),dim(udv$v[,1:20]))

    ## Do a little 'fake': here the pattern is not a geographical map but weight
    ## for the EOFs.
    dim(pattern) <- c(1,dim(pattern))
    if (verbose) str(pattern)
    attr(eof,'eigenvalues') <- udv$d[1:20]
    attr(eof,'pattern') <- rep(1,20)
    names(eof) <- paste("X.",1:20,sep="")
    
    class(eof) <- class(X[[1]])

    ## Downscale the results:
    if (verbose) print("- - - > DS")
    ds <- try(DS(y,eof,ip=ip,verbose=verbose))
    if (inherits(ds,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds[[1]],con=flog)
    } else {
      if (verbose) print("post-processing")
      z <- attr(ds,'appendix.1')
      i1 <- is.element(years,year(z))
      i2 <- is.element(year(z),years)
    #

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
      
      attr(z,'scorestats') <- c(1-r.xval,mdiff,srati,arati,res.trend,ks,ar,
                                1-ds.ratio,1-var(xval[,2])/var(xval[,1]))
      if (verbose) print('scorestats')
      if (verbose) print(attr(z,'scorestats'))
      dse[[i]] <- z
      
      quality <- 100*(1-mean(abs(scorestats[i,]),na.rm=TRUE))
      qcol <- quality
      qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100

      index(z) <- year(z); index(ds) <- year(ds)
      if (plot) {
        lines(z,lwd=2,col=cols[qcol])
        lines(y,type="b",pch=19)
        lines(ds,lwd=2,col="grey")
     }
      #
      R2 <- round(100*sd(xval[,2])/sd(xval[,1]),2)
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(100*r.xval,2),
                  "R2=",R2,'% ','Common EOF: bias=',round(mdiff,2),
                  ' sd1/sd2=',round(srati,3),
                  "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',round(quality)))
    }
  }

  #
  names(dse) <- gcmnm
  #X <- attrcp(y,X)
  attr(dse,'station') <- y
  attr(dse,'predictor') <- c(attr(PRE1,'source'),attr(PRE2,'source'),attr(PRE3,'source'))
  attr(dse,'domain') <- list(lon=lon,lat=lat)
  attr(dse,'scenario') <- rcp
  attr(dse,'history') <- history.stamp(y)
  class(X) <- c("dsensemble","zoo")
  save(file="DSensemble.rda",X)
  print("---")
  invisible(dse)
}


DSensemble.mu.worstcase <- function(y,plot=TRUE,path="CMIP5.monthly/",
                                    predictor="ERA40_t2m_mon.nc",
                                    rcp="rcp45",biascorrect=FALSE,n=6,
                                    lon=c(-20,20),lat=c(-10,10),it=NULL,rel.cord=TRUE,
                                    select=NULL,FUN="wetmean",type='ncdf4',
                                    pattern="tas_Amon_ens_",mask=FALSE,
                                    verbose=FALSE,ds.1900.2099=TRUE) {
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
    pre <- retrieve(ncfile=predictor,lon=lon,lat=lat,type=type,verbose=verbose)
    if (mask) pre=mask(pre,land=TRUE)
    pre <- spatial.avg.field(C.C.eq(pre))
  } else if (inherits(predictor,'field')) pre <- spatial.avg.field(predictor)
  rm("predictor"); gc(reset=TRUE)
  ## Estimate the reference level
  normal61.90 <- mean(coredata(subset(pre,it=c(1961,1990))))
  x <- aggregate(pre,by=month,FUN="mean")

  if (verbose) print('Prepare the calibration data')
  ## Loop over PCAs if multipple stations

  if (n>1) results <- list()

  ## Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
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
      gcm <- retrieve(ncfile = ncfiles[select[i]],lon=lon,lat=lat,
                      type=type,verbose=FALSE)
      if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099))
      if (verbose) print(paste('mask=',mask))
      if (mask) gcm <- mask(gcm,land=TRUE)
      gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
      if (verbose) print('spatial average')
      GCM <- spatial.avg.field(C.C.eq(gcm))
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


DSensemble.pca <- function(y,plot=TRUE,path="CMIP5.monthly/",
                           rcp="rcp45",biascorrect=FALSE,
                           predictor="ERA40_t2m_mon.nc",
                           non.stationarity.check=FALSE,
                           ip=1:16,lon=c(-30,20),lat=c(-20,10), it=NULL,
                           rel.cord=TRUE,
                           select=NULL,FUN="mean",rmtrend=TRUE,
                           FUNX="mean",xfuns='C.C.eq',threshold=1,type='ncdf4',
                           pattern="tas_Amon_ens_",verbose=FALSE,
                           file.ds="DSensemble.rda",path.ds=NULL,nmin=NULL,ds.1900.2099=TRUE) {

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
    t2m <- retrieve(ncfile=predictor,lon=lon,lat=lat,
                    type=type,verbose=verbose)
      if (!is.null(it)) {
        if (verbose) print('Extract some months or a time period')
        if (verbose) print(it)
        t2m <- subset(t2m,it=it,verbose=verbose)
        ## if it is character, then then extraction of months reduces number of
        ## days per year.
      }
  } else if (inherits(predictor,'field')) {
    t2m <- predictor
    lon <- range(lon(t2m))
    lat <- range(lat(t2m))
  }
  ## If some months are selected, make sure that the minimum number of months
  ## requiired in the annual aggregation is updated
  if ((is.null(nmin)) & (is.character(it))) nmin <- length(it)

  if (inherits(y,'season')) {
    if (verbose) print('seasonal data')
    T2M <- as.4seasons(t2m,FUN=FUNX,nmin=nmin)
    T2M <- matchdate(T2M,y)

    # Recursive: do each season seperately if there are more than one season
    if (length(table(season(y)))>1) {
      if (verbose) print('--- Apply DS to seasons seperately ---')
      Z <- list(info='DSensemble.pca for different seasons')
      for (season in names(table(season(y)))) {
        if (verbose) print(paste('Select',season))
        z <- DSensemble.pca(subset(y,it=season),plot=plot,path=path,
                            rcp=rcp,biascorrect=biascorrect,predictor=T2M,
                            non.stationarity.check=non.stationarity.check,
                            ip=ip,lon=lon,lat=lat,rel.cord=FALSE,
                            select=select,FUN=FUN,rmtrend=rmtrend,
                            FUNX=FUNX,xfuns=xfuns,threshold=threshold,type=type,
                            pattern=pattern,verbose=verbose,nmin=nmin)
        eval(parse(text=paste('Z$',season,' <- z',sep='')))
      }
      if (verbose) print('--- Results returned as a list ---')
      return(Z)
    }

  } else if (inherits(y,'annual')) {
    if (verbose) print('annual data')
    if  (!inherits(predictor,'field')) T2M <- t2m else if (!is.annual(t2m)) {
      if (verbose) print(paste('Aggregate annually',FUNX,'for calibration'))
        if (sum(is.element(FUNX,xfuns))==0)
          T2M <- annual(t2m,FUN=FUNX,nmin=nmin) else
          eval(parse(text=paste('T2M <- annual(',FUNX,'(t2m),FUN="mean",nmin=nmin)',sep="")))
      } else T2M <- t2m
    ## Match the date
    T2M <- matchdate(T2M,y)
  } else if (inherits(y,'month')) {
    if (verbose) print('monthly data')
    T2M <- matchdate(t2m,y)
    ## browser()
    ##if (FUNX=='C.C.eq') 
    ##  T2M <- mask(T2M,land=TRUE)
  }
  if (inherits(T2M,"eof")) T2M <- as.field(T2M)
  rm("predictor","t2m"); gc(reset=TRUE)
  
  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
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
  colnames(scorestats) <- c("1-r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration',
                            "1-R2")

  t <- as.Date(paste(years,months,'01',sep='-'))

  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)

  if (plot) {
    par(bty='n')
    index(y) <- year(y)
    plot.zoo(y[,1],lwd=3,main='PC1',ylab='',xlab='',xlim=range(years))
  }
  
  flog <- file("DSensemble.pca-log.txt","at")
  
  ## Set up a list environment to keep all the results
  dse.pca <- list(info='DSensemble.pca',pca=y) ## KMP 06-08-2015
  if (verbose) print("loop...") 
  for (i in 1:N) {
    if (verbose) print(ncfiles[select[i]])
    gcm <- retrieve(ncfile = ncfiles[select[i]],type=type,
                          lon=range(lon(T2M))+c(-2,2),
                          lat=range(lat(T2M))+c(-2,2),verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099))
    if (!is.null(it)) {
      if (verbose) print('Extract some months or a time period')
      if (is.null(nmin)) warning(paste("The argument 'it' is set but not 'nmin'; it=",
                                       paste(it,collapse="-")))
      if (verbose) print(it)
      gcm <- subset(gcm,it=it,verbose=verbose)
    }
    gcmnm.i <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-r")

    if (verbose) {
        print(paste('Extract month/season/annual data nmin=',nmin))
        print(class(y))
        print(FUNX)
    }
    if (inherits(y,'season')) {
      if (sum(is.element(FUNX,xfuns))==0) {
          if (verbose) print('No special transformation') 
          GCM <- as.4seasons(gcm,FUN=FUNX,nmin=nmin)
       } else {
           if (verbose) print('Need to aggregate FUNX(gcm)')
           eval(parse(text=paste('GCM <- as.4seasons(',FUNX,'(gcm),FUN="mean",nmin=nmin)',sep="")))
       }
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

    if (verbose) print("- - - > DS")
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
      ##attr(z,'model') <- attr(ds,'model') ## KMP 09-08-2015
      ## model takes up too much space! can it be stored more efficiently?
      attr(z,'predictor.pattern') <- attr(ds,'predictor.pattern')
      attr(z,'evaluation') <- attr(ds,'evaluation')
    
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
      qcol <- quality
      qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
      R2 <- round(100*sd(xval[,2])/sd(xval[,1]),2)
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(100*r.xval,2),
                  "R2=",R2,'% ','Common EOF: bias=',round(mdiff,2),
                  ' sd1/sd2=',round(srati,3),
                  "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',
                  round(quality)))

      if (plot) {
        index(y) <- year(y); index(z) <- year(z)
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
  attr(dse.pca,'predictor') <- attr(T2M,'source')
  attr(dse.pca,"longname") <- attr(y,"longname")
  attr(dse.pca,'domain') <- list(lon=lon,lat=lat)
  attr(dse.pca,'scorestats') <- scorestats
  attr(dse.pca,'path') <- path
  attr(dse.pca,'scenario') <- rcp
  attr(dse.pca,'variable') <- attr(y,"variable")[1]
   attr(dse.pca,'unit') <- attr(y,"unit")[1]
  attr(dse.pca,'history') <- history.stamp(y)
  if (non.stationarity.check)
    attr(dse.pca,'on.stationarity.check') <- difference.z else
    attr(dse.pca,'on.stationarity.check') <- NULL
  class(dse.pca) <- c("dsensemble","pca","list")

  if(!is.null(path.ds)) file.ds <- file.path(path.ds,file.ds)
  save(file=file.ds,dse.pca)
  if (verbose) print("---")
  invisible(dse.pca)
}


DSensemble.eof <- function(y,lplot=TRUE,path="CMIP5.monthly",
                           rcp="rcp45",biascorrect=FALSE,
                           predictor="ERA40_slp_mon.nc",
                           non.stationarity.check=FALSE,
                           ip=1:5,lon=c(-30,20),lat=c(-20,10),it=NULL,
                           rel.cord=TRUE,nmin=NULL,lev=NULL,levgcm=NULL,
                           select=NULL,FUN="mean",rmtrend=TRUE,
                           FUNX="mean",xfuns='C.C.eq',threshold=1,type='ncdf4',
                           pattern="psl_Amon_ens_",verbose=FALSE,
                           file.ds="DSensemble.eof.rda",path.ds=NULL,ds.1900.2099=TRUE) {

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
    slp <- retrieve(ncfile=predictor,lon=lon,lat=lat,lev=lev,
                    type=type,verbose=verbose) else
  if (inherits(predictor,'field'))
    slp <- subset(predictor,is=list(lon=lon,lat=lat))

  if (inherits(y,'season')) {
    if (verbose) print('seasonal data')
    SLP <- as.4seasons(slp,FUN=FUNX)
    SLP <- matchdate(SLP,y)

    if (length(table(season(y)))>1) {
      if (verbose) print('--- Apply DS to seasons seperately ---')
      Z <- list(info='DSensemble.eof for different seasons')
      for (season in names(table(season(y)))) {
        if (verbose) print(paste('Select',season))
        z <- DSensemble.eof(subset(y,it=season),lplot=lplot,path=path,
                            rcp=rcp,biascorrect=biascorrect,predictor=SLP,
                            non.stationarity.check=non.stationarity.check,
                            ip=ip,lon=lon,lat=lat,rel.cord=FALSE,
                            select=select,FUN=FUN,rmtrend=rmtrend,
                            FUNX=FUNX,threshold=threshold,type=type,
                            pattern=pattern,verbose=verbose,nmin=nmin)
        eval(parse(text=paste('Z$',season,' <- z',sep='')))
      }
      if (verbose) print('--- Results returned as a list ---')
      return(Z)
    }

  } else if (inherits(y,'annual')) {
    if (verbose) print('annual data')
    if (!is.null(it)) {
      if (verbose) print('Extract some months of a time period')
      if (verbose) print(it)
      slp <- subset(slp,it=it)
      if ((is.null(nmin)) & (is.character(it))) nmin <- length(it)
    }
    if (!is.annual(slp)) {
      if (FUNX!='C.C.eq')
        SLP <- annual(slp,FUN=FUNX,nmin=nmin) else
        SLP <- annual(C.C.eq(slp),FUN='mean',nmin=nmin)
  } else SLP <- slp
    SLP <- matchdate(SLP,y)
  } else if (inherits(y,'month')) {
    SLP <- matchdate(slp,y)
  }
  if (inherits(SLP,"eof")) SLP <- as.field(SLP)
  rm("predictor","slp"); gc(reset=TRUE)

  if(!inherits(y,"eof")) y <- EOF(y)
      
  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
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

  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)

  flog <- file("DSensemble.eof-log.txt","at")

  ## Set up a list environment to keep all the results
  dse.eof <- list(info='DSensemble.eof',eof=y) 
  if (verbose) print("loop...") 
  for (i in 1:N) {
    if (verbose) print(ncfiles[select[i]])
    gcm <- retrieve(ncfile = ncfiles[select[i]],type=type,
                          lon=range(lon(SLP))+c(-2,2),
                          lat=range(lat(SLP))+c(-2,2),
                          lev=levgcm,verbose=verbose)
    if (ds.1900.2099) gcm <- subset(gcm,it=c(1900,2099))
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
    gcmnm.i <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-r")
    if (verbose) print(class(y))

    ## REB 2016-10-25
    if (inherits(y,'season')) {
      if (verbose) print(paste('Seasonally aggregated',FUNX,'for GCM'))
      if (sum(is.element(FUNX,xfuns))==0) {
          if (verbose) print('No special transformation') 
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
    if (non.stationarity.check) {
      testGCM <- subset(GCM,it=range(year(SLP)))      
      testy <- regrid(testGCM,is=as.field(y)) 
      attr(testGCM,'source') <- 'testGCM'
      testZ <- combine(testGCM,GCM)            
      rm("testGCM"); gc(reset=TRUE)
    }
    rm("gcm","GCM"); gc(reset=TRUE)

    if (verbose) print("- - - > DS")
    if (biascorrect) Z <- biasfix(Z)
    
    diag <- diagnose(Z)
    
    ds <- try(DS(y,Z,ip=ip,verbose=verbose))
    if(inherits(ds,"try-error")) {
      print(paste("esd failed for",gcmnm.i))
    } else {
      if (verbose) print("post-processing")
      gcmnm[i] <- gcmnm.i
   
      ## Keep the results for the projections:
      if (verbose) print('Extract the downscaled projection')
      z <- attr(ds,'appendix.1') ## KMP 09.08.2015
      ##attr(z,'model') <- attr(ds,'model') ## KMP 09-08-2015
      ## model takes up too much space! can it be stored more efficiently?
      attr(z,'predictor.pattern') <- attr(ds,'predictor.pattern')
      attr(z,'evaluation') <- attr(ds,'evaluation')
    
      cl <- paste('dse.eof$i',i,'_',gsub('-','.',gcmnm[i]),' <- z',sep='')
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
   
    if (verbose) print('Downscaling finished')
  
  }
  
  ## Unpacking the information tangled up in GCMs, PCs and stations:
  ## Save GCM-wise in the form of PCAs
  gcmnm <- gsub('-','.',gcmnm)

  #Z <- attrcp(y,Z)
  if (verbose) print('Set attributes')
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


DSensemble.field <- function(y,plot=TRUE,path="CMIP5.monthly/",
                           rcp="rcp45",biascorrect=FALSE,
                           predictor="ERA40_t2m_mon.nc",
                           non.stationarity.check=FALSE,
                           ip=1:16,lon=c(-30,20),lat=c(-20,10),
                           it=c('djf','mam','jja','son'),
                           rel.cord=TRUE,
                           select=NULL,FUN="mean",rmtrend=TRUE,
                           FUNX="mean",xfuns='C.C.eq',threshold=1,type='ncdf4',
                           pattern="tas_Amon_ens_",verbose=FALSE,
                           file.ds="DSensemble.rda",path.ds=NULL,nmin=NULL,ds.1900.2099=TRUE) {
  ## For downscaling gridded predictand. This is a wrap-around which extracts the season or aggregates
  ## to annual values and then calls the other types for the downscaling.
  
}

DSensemble.station <- function(y,plot=TRUE,path="CMIP5.monthly/",
                           rcp="rcp45",biascorrect=FALSE,
                           predictor="ERA40_t2m_mon.nc",
                           non.stationarity.check=FALSE,
                           ip=1:16,lon=c(-30,20),lat=c(-20,10),
                           it=c('djf','mam','jja','son'),
                           rel.cord=TRUE,
                           select=NULL,FUN="mean",rmtrend=TRUE,
                           FUNX="mean",xfuns='C.C.eq',threshold=1,type='ncdf4',
                           pattern="tas_Amon_ens_",verbose=FALSE,
                           file.ds="DSensemble.rda",path.ds=NULL,nmin=NULL,ds.1900.2099=TRUE) {
}
