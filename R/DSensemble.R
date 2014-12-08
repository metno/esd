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

  stopifnot(!missing(y),inherits(y,"station"),file.exists(paste(path,rcp,sep="")))
  
  if (is.null(attr(y,'aspect'))) attr(y,'aspect') <- "original"
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
                           non.stationarity.check=FALSE,
                           area.mean.expl=FALSE,
                           eofs=1:6,lon=c(-20,20),lat=c(-10,10),
                           select=NULL,FUN="mean",FUNX="mean",
                           pattern="tas_Amon_ens_",verbose=FALSE) {
  
  #print("predictand")
  #if ( (deparse(substitute(FUN))=='sd') | (deparse(substitute(FUN))=='ar1') )
  if ((FUN=='sd') | (FUN =='ar1')) {
    y <- anomaly(y)
    attr(y,'aspect') <- 'original'
  }
  
  ya <- annual(y,FUN=FUN)
  y <- as.4seasons(y,FUN=FUN)
  

#  yr <- as.4seasons(ya,FUN=ltp)
  if (!is.na(attr(y,'longitude')))
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude')))
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )
  #lon <- round( attr(y,'longitude') + lon )
  #lat <- round( attr(y,'latitude') + lat )

  # The units
  if ( (attr(y,'unit') == "deg C") | (attr(y,'unit') == "degree Celsius") )
        unit <- expression(degree*C) else
        unit <- attr(y,'unit')

  #ylim <- switch(deparse(substitute(FUN)),
  ylim <- switch(FUN,
                 'mean'=c(-2,8),'sd'=c(-0.5,1),'ar1'=c(-0.5,0.7))
  if (plot) {
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

  # Get the predictor: ERA40
  #print("predictor")
  #data(ferder)
  #lon <- lon(ferder) + c(-20,20)
  #lat <- lat(ferder) + c(-10,10)
  
  if (is.character(predictor))
    t2m <- retrieve.ncdf4(ncfile=predictor,lon=lon,lat=lat) else
  if (inherits(predictor,'field'))
    t2m <- subset(predictor,is=list(lon=lon,lat=lat))

  #rm("predictor"); gc(reset=TRUE)
  #t2m <- t2m.ERA40(lon=lon,lat=lat)
  T2M <- as.4seasons(t2m,FUN=FUNX)
  # Fix - there is a bug with 'T2M <- as.4seasons(t2m,FUN=FUNX)' - the date is not correct
  DJF <- subset(as.4seasons(t2m,FUN=FUNX),it='djf')
  MAM <- subset(as.4seasons(t2m,FUN=FUNX),it='mam')
  JJA <- subset(as.4seasons(t2m,FUN=FUNX),it='jja')
  SON <- subset(as.4seasons(t2m,FUN=FUNX),it='son')

  ok1 <- is.finite(rowSums(DJF))
  DJF <- subset(DJF,it=range(year(DJF)[ok1]))
  ok2 <- is.finite(rowSums(MAM))
  MAM <- subset(MAM,it=range(year(MAM)[ok2]))
  ok3 <- is.finite(rowSums(JJA))
  JJA <- subset(JJA,it=range(year(JJA)[ok3]))
  ok4 <- is.finite(rowSums(SON))
  SON <- subset(SON,it=range(year(SON)[ok4]))

  #ok <- is.finite(rowSums(T2M))
  #T2M <- subset(T2M,it=range(year(T2M)[ok]))
  #eofjja <- EOF(subset(T2M,it=3)); save(file='eof.rda',eofjja)
  #load('T2M.rda')
  rm("t2m"); gc(reset=TRUE)
  
  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) select <- 1:N else
                       N <- length(select)
  if (verbose) print(ncfiles[select])

  # set up results matrix and tables of diagnostics:
  years <- sort(rep(1900:2100,4))
  months <- rep(c(1,4,7,10),length(1900:2100))
  m <- length(years)
  X <- matrix(rep(NA,N*m),N,m)
  gcmnm <- rep("",N)
  scorestats <- matrix(rep(NA,N*8),N,8)
  colnames(scorestats) <- c("r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration')

  t <- as.Date(paste(years,months,'01',sep='-'))

  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)

  # Quick test:
  #load('control.ds.1.rda'); T2M.ctl -> T2M
  flog <- file("DSensemble.t2m-log.txt","at")
  
  if (verbose) print("loop...") 
  for (i in 1:N) {
    gcm <- retrieve.ncdf4(ncfile = ncfiles[select[i]],
                          lon=range(lon(T2M)),lat=range(lat(T2M)))
    #gcmnm[i] <- attr(gcm,'model_id')
    gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
    # REB: 30.04.2014 - new lines...
    DJFGCM <- subset(as.4seasons(gcm,FUN=FUNX),it='djf')
    MAMGCM <- subset(as.4seasons(gcm,FUN=FUNX),it='mam')
    JJAGCM <- subset(as.4seasons(gcm,FUN=FUNX),it='jja')
    SONGCM <- subset(as.4seasons(gcm,FUN=FUNX),it='son')
    rm("gcm"); gc(reset=TRUE)
    X.DJF <- combine(DJF,DJFGCM)
    X.MAM <- combine(MAM,MAMGCM)
    X.JJA <- combine(JJA,JJAGCM)
    X.SON <- combine(SON,SONGCM)

    #load('control.ds.1.rda')
    #X.JJA <- PREGCM
    # REB: 30.04.2014 - new lines...
    if (verbose) print("- - - > EOFs")
    
    Z1 <- try(EOF(X.DJF,area.mean.expl=area.mean.expl))
    if (inherits(Z1,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(Z1[[1]],con=flog)
    }
    Z2 <- try(EOF(X.MAM,area.mean.expl=area.mean.expl))
    if (inherits(Z2,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(Z2[[1]],con=flog)
    }
    Z3 <- try(EOF(X.JJA,area.mean.expl=area.mean.expl))
    if (inherits(Z3,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(Z3[[1]],con=flog)
    }
    Z4 <- try(EOF(X.SON,area.mean.expl=area.mean.expl))
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

    # REB: 30.04.2014 - new lines...
    if (verbose) print("- - - > DS")
    if (biascorrect) Z1 <- biasfix(Z1)
    ds1 <- try(DS(subset(y,it='djf'),Z1,eofs=eofs))
    if (inherits(ds1,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds1[[1]],con=flog)
    }
    if (biascorrect) Z2 <- biasfix(Z2)
    ds2 <- try(DS(subset(y,it='mam'),Z2,eofs=eofs))
    if (inherits(ds2,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds2[[1]],con=flog)
    }
    if (biascorrect) Z3 <- biasfix(Z3)
    ds3 <- try(DS(subset(y,it='jja'),Z3,eofs=eofs))
    if (inherits(ds3,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds3[[1]],con=flog)
    }
    if (biascorrect) Z4 <- biasfix(Z4)
    ds4 <- try(DS(subset(y,it='son'),Z4,eofs=eofs))
    if (inherits(ds4,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds4[[1]],con=flog)
    }
    if (verbose) print("Combine the 4 seasons")
    ## browser()
    ds <- try(combine(list(ds1,ds2,ds3,ds4)))
    ##ds <- c(zoo(ds1),zoo(ds2),zoo(ds3),zoo(ds4))
    ##ds <- attrcp(y,ds)
    ##attr(ds,'appendix.1') <- c(attr(ds1,'appendix.1'),
    ##                           attr(ds2,'appendix.1'),
    ##                           attr(ds3,'appendix.1'),
    ##                           attr(ds4,'appendix.1'))
    ##class(ds) <- class(y)
    
    if (inherits(ds,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds[[1]],con=flog)
    } else {
      x.val <- c(crossval(ds1),crossval(ds2),crossval(ds3),crossval(ds4))
      attr(ds,'evaluation') <- x.val
#    browser()
#    ds <- DS(y,Z,biascorrect=biascorrect,eofs=eofs,
#             area.mean.expl=area.mean.expl,verbose=verbose)
#    ds <- DS.t2m.season.field(y,Z,biascorrect=biascorrect,eofs=eofs,
#                              area.mean.expl=area.mean.expl,verbose=verbose)
      if (verbose) print("post-processing")
      z <- attr(ds,'appendix.1')
    #save(file='inside.dsens.1.rda',ds,y,Z)
    ##browser()
    # The test lines are included to assess for non-stationarity
      if (non.stationarity.check) {
        testds <- DS(testy,testZ,biascorrect=biascorrect,
                     area.mean.expl=area.mean.expl,eofs=eofs)   # REB 29.04.2014
        testz <- attr(testds,'appendix.1')                      # REB 29.04.2014
        difference.z <- testy - testz                           # REB 29.04.2014
      }
      i1 <- is.element(paste(years,months,sep='-'),
                       paste(year(z),month(z),sep='-'))
      i2 <- is.element(paste(year(z),month(z),sep='-'),
                       paste(years,months,sep='-'))
    #browser()
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
                               res.trend,'D/decade; K.S. p-val',
                               ks,'; AR(1)=',ar))

    # Evaluation: here are lots of different aspects...
    # Get the diagnostics: this is based on the analysis of common EOFs...

      xval <- attr(ds,'evaluation')
      r.xval <- round(cor(xval[,1],xval[,2]),3)
      if (verbose) print(paste("x-validation r=",r.xval))
    
      dsa <- annual(ds)                     # annual mean value

      xy <- merge.zoo(annual(z),ya)
      ds.ratio <- round(sd(xy[,1],na.rm=TRUE)/sd(xy[,2],na.rm=TRUE),4)
      if (verbose) print(paste("sd ratio=",ds.ratio))

    #print(names(attributes(ds)))
      if (biascorrect) {
        diag <- attr(ds,'diagnose')
        if ( (verbose) & !is.null(diag)) str(diag)
      } else diag <- NULL
    
    # diagnose for ds-objects
      
      if (verbose) print('...')
      if (is.null(diag)) {
        diag <- diagnose(z,plot=FALSE)
        scorestats[i,] <- c(1-r.xval,NA,NA,NA,res.trend,ks,ar,ds.ratio)
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
        srati <- mean(1 - c(diag$s.1$sd.ratio[1],diag$s.2$sd.ratio[1],
                            diag$s.3$sd.ratio[1],diag$s.4$sd.ratio[1]))
        arati <- mean(1 - c(diag$s.1$autocorr.ratio[1],diag$s.2$autocorr.ratio[1],
                            diag$s.3$autocorr.ratio[1],diag$s.4$autocorr.ratio[1]))
      }
      scorestats[i,] <- c(1-r.xval,mdiff,srati,arati,res.trend,ks,ar,ds.ratio)
      if (verbose) print(scorestats[i,])

      quality <- 100*(1-mean(scorestats[i,]))
      qcol <- quality
      qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100
     
      if (plot) {
        lines(annual(z),lwd=2,col=cols[qcol])
        lines(ya,type="b",pch=19)
        lines(dsa,lwd=2,col="grey")
      #zoo(subset(y,it=3),order.by=year(subset(y,it=3)))
#      jja <- zoo(z[is.element(month(z),7)],order.by=year(z[is.element(month(z),7)]))
#      lines(jja,lwd=2,col=cols[quality])
#      lines(zoo(subset(y,it=3),order.by=year(subset(y,it=3))),type="b",pch=19)
#      lines(zoo(subset(ds,it=3),order.by=year(subset(ds,it=3))),lwd=2,col="grey")
      #browser()
      }
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(r.xval,2),
                  "R2=",round(100*sd(xval[,2])/sd(xval[,1]),2),'% ',
                'Common EOF: bias=',round(mdiff,2),' 1- sd1/sd2=',round(srati,3),
                "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',round(quality)))
    }
  }

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
  if (area.mean.expl) attr(X,'area.mean.expl') <- TRUE else
                      attr(X,'area.mean.expl') <- FALSE
  class(X) <- c("dsensemble","zoo")
  save(file="DSensemble.rda",X)
  print("---")
  invisible(X)
} 
#save(file=paste("dscmip5_",attr(y,'location'),"_",N,"_rcp4.5.rda",sep=""),rcp4.5)

DSensemble.precip <- function(y,plot=TRUE,path="CMIP5.monthly/",
                              rcp="rcp45",biascorrect=FALSE,
                              predictor="ERA40_pr_mon.nc",
                              non.stationarity.check=FALSE,
                              area.mean.expl=FALSE,
                              eofs=1:6,lon=c(-10,10),lat=c(-10,10),
                              select=NULL,FUN="exceedance",
                              FUNX="sum",threshold=1,
                              pattern="pr_Amon_ens_",verbose=FALSE) {
  # FUN: exceedance, wetfreq, wet, dry

  if (verbose) print('DSensemble.precip')
#  if (deparse(substitute(FUN))=='spell') {
  if ( (FUN=='wet') | (FUN=='dry')) {
    y <- spell(y,threshold=threshold)
    y <- annual(y)
    #plot(y); browser()
    y <- switch(FUN,'wet'=subset(y,is=1),'dry'=subset(y,is=2))
  } else
#    y <- annual(y,FUN=FUN,threshold=threshold)
    if (sum(is.element(names(formals(FUN)),'threshold')==1))
        y <- annual(y,FUN=FUN,threshold=threshold) else
        y <- annual(y,FUN=FUN)
  #browser()
  index(y) <- year(y)
  
  if (!is.na(attr(y,'longitude')))
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude')))
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )
  
  # Get the predictor: ERA40
  if (verbose) print("predictor")
  if (is.character(predictor))
    pre <- retrieve(ncfile=predictor,lon=lon,lat=lat) else
  if (inherits(predictor,'field')) pre <- predictor
  rm("predictor"); gc(reset=TRUE)
  attr(pre,"source") <- "ERA40"

  # Use proportional variaions
  if (verbose) print("Annual mean")
  if (FUNX!='C.C.eq')
    PREX <- annual(pre,FUN=FUNX) else
    PREX <- annual(C.C.eq(pre),FUN='mean')
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
  scorestats <- matrix(rep(NA,N*8),N,8)
  colnames(scorestats) <- c("r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration')

  flog <- file("DSensemble.precip-log.txt","at")
  for (i in 1:N) {
    #browser()
    gcm <- retrieve(ncfile = ncfiles[select[i]],
                    lon=range(lon(PRE)),lat=range(lat(PRE)))
    gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
    #gcmnm[i] <- attr(gcm,'model_id')
    if (verbose) print(varid(gcm))
    
    if (FUNX!='C.C.eq')
      GCMX <- annual(gcm,FUN=FUNX) else
      GCMX <- annual(C.C.eq(gcm),FUN='mean')
#    GCM <- zoo(100*coredata(GCMX)/colMeans(coredata(subset(GCMX,it=1961:1990))),order.by=year(GCMX))
    if (is.precip(GCMX)) {
          GCM <- zoo(100*coredata(GCMX)/mean(c(coredata(subset(GCMX,it=1961:1990)))),
                     order.by=year(GCMX))
      GCM <- attrcp(GCMX,GCM)
      attr(GCM, "unit" ) <- "%"
      attr(GCM, "dimensions" ) <- attr(GCMX, "dimensions" )
      class(GCM) <- class(GCMX)
    } else
          GCM <- GCMX
    #browser()
    #str(GCM)
    model.id <- attr(gcm,'model_id')
    rm("gcm","GCMX"); gc(reset=TRUE)
    if (verbose) print("combine")
    #browser()
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
      testds <- DS(testy,testZ,biascorrect=biascorrect,
                   area.mean.expl=area.mean.expl,eofs=eofs)  # REB 29.04.2014
      testz <- attr(testds,'appendix.1')                     # REB 29.04.2014
      difference.z <- testy - testz                          # REB 29.04.2014
    }
    
    if (verbose) print("diagnose")
    diag <- diagnose(Z)
    if (biascorrect) Z <- biasfix(Z)
    if (verbose) print("- - - > DS")
    ds <- try(DS(y,Z,eofs=eofs,area.mean.expl=area.mean.expl,verbose=verbose))
    if (inherits(ds,"try-error")) {    
      writeLines(gcmnm[i],con=flog)
      writeLines(ds[[1]],con=flog)
    } else {
      if (verbose) print("post-processing")
      z <- attr(ds,'appendix.1')
      i1 <- is.element(years,year(z))
      i2 <- is.element(year(z),years)
    #browser()
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

    #browser()
      xy <- merge.zoo(z,y)
      ds.ratio <- sd(xy[,1],na.rm=TRUE)/sd(xy[,2],na.rm=TRUE)
    
    # Extract the mean score for leading EOF from the 4 seasons:
      mdiff <- diag$mean.diff[1]/diag$sd0[1]
      srati <- 1 - diag$sd.ratio[1]
      arati <- 1 - diag$autocorr.ratio[1]
      scorestats[i,] <- c(1-r.xval,mdiff,srati,arati,res.trend,ks,ar,ds.ratio)
      
      quality <- 100*(1-mean(scorestats[i,]))
      qcol <- quality
      qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100

      index(z) <- year(z); index(ds) <- year(ds)
      if (plot) {
        lines(z,lwd=2,col=cols[qcol])
        lines(y,type="b",pch=19)
        lines(ds,lwd=2,col="grey")
     }
      #browser()
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(r.xval,2),
                  "R2=",round(100*sd(xval[,2])/sd(xval[,1]),2),'% ',
                  'Common EOF: bias=',round(mdiff,2),' 1- sd1/sd2=',round(srati,3),
                  "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',round(quality)))
    }
  }

  #browser()
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
  if (area.mean.expl) attr(X,'area.mean.expl') <- TRUE else
                      attr(X,'area.mean.expl') <- FALSE
  attr(X,'history') <- history.stamp(y)
  class(X) <- c("dsensemble","zoo")
  save(file="DSensemble.rda",X)
  print("---")
  invisible(X)
}



DSensemble.mu <- function(y,plot=TRUE,path="CMIP5.monthly/",
                          rcp="rcp45",biascorrect=FALSE,
                          predictor="ERA40_t2m_mon.nc",
                          non.stationarity.check=FALSE,
                          eofs=1:16,lon=c(-30,20),lat=c(-20,10),
                          select=NULL,FUN="wetmean",
                          FUNX="C.C.eq",threshold=1,
                          pattern="tas_Amon_ens_",verbose=FALSE) {

# This function is for downscaling wet-day mean using SLP to
# remove the internal variations and then the global mean temperature
# to downscale the slow changes due to a global warming.

  # Get the global mean temeprature: pentads

  # Or a combination of OLR + C.C.eq(t2m) + SLP

  Z <- DSensemble.precip(y=y,plot=plot,path=path,rcp=rcp,
                         biascorrect=biascorrect,predictor=predictor,
                         non.stationarity.check=non.stationarity.check,
                         eofs=eofs,lon=lon,lat=lat,select=select,
                         FUN=FUN,FUNX=FUNX,threshold=threshold,
                         pattern=pattern,verbose=verbose)
  return(Z)
}


DSensemble.pca <- function(y,plot=TRUE,path="CMIP5.monthly/",
                          rcp="rcp45",biascorrect=FALSE,
                          predictor="ERA40_t2m_mon.nc",
                          non.stationarity.check=FALSE,
                          eofs=1:16,lon=c(-30,20),lat=c(-20,10),
                          select=NULL,FUN="wetmean",
                          FUNX="C.C.eq",threshold=1,
                          pattern="tas_Amon_ens_",verbose=FALSE) {
  
  # This function is for downscaling PCA to represent a group of stations
  

}
