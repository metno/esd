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

DSensemble.t2m <- function(y,plot=TRUE,path="~/CMIP5.monthly/",
                           predictor="~/ERA40_t2m_mon.nc",
                           rcp="rcp45",biascorrect=FALSE,
                           non.stationarity.check=FALSE,
                           area.mean.expl=FALSE,
                           eofs=1:6,lon=c(-20,20),lat=c(-10,10),
                           select=NULL,FUN="mean",FUNX="mean",
                           pattern="tas_Amon_ens_",verbose=FALSE,nmin=NULL) {
  
  #print("predictand")
  #if ( (deparse(substitute(FUN))=='sd') | (deparse(substitute(FUN))=='ar1') )
  if ((FUN=='sd') | (FUN =='ar1')) {
    y <- anomaly(y)
    attr(y,'aspect') <- 'original'
  }
  
  ya <- annual(y,FUN=FUN,nmin=nmin*4)
  y <- as.4seasons(y,FUN=FUN,nmin=nmin)
  

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
    t2m <- retrieve(ncfile=predictor,lon=lon,lat=lat) else
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
    gcm <- retrieve(ncfile = ncfiles[select[i]],
                          lon=range(lon(T2M))+c(-2,2),
                          lat=range(lat(T2M))+c(-2,2))
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
    ## 
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
#    
#    ds <- DS(y,Z,biascorrect=biascorrect,eofs=eofs,
#             area.mean.expl=area.mean.expl,verbose=verbose)
#    ds <- DS.t2m.season.field(y,Z,biascorrect=biascorrect,eofs=eofs,
#                              area.mean.expl=area.mean.expl,verbose=verbose)
      if (verbose) print("post-processing")
      z <- attr(ds,'appendix.1')
    #save(file='inside.dsens.1.rda',ds,y,Z)
    ##
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
      #
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
                              pattern="pr_Amon_ens_",verbose=FALSE,nmin=NULL) {
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
    #
    gcm <- retrieve(ncfile = ncfiles[select[i]],
                    lon=range(lon(PRE))+c(-2,2),lat=range(lat(PRE))+c(-2,2))
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
      #
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(r.xval,2),
                  "R2=",round(100*sd(xval[,2])/sd(xval[,1]),2),'% ',
                  'Common EOF: bias=',round(mdiff,2),' 1- sd1/sd2=',round(srati,3),
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
                          predictor=list(t2m="data/ncep/air.mon.mean.nc",
                                         olr="data/ncep/OLR.mon.mean.nc",
                                         slp="data/ncep/slp.mon.mean.nc"),
                          non.stationarity.check=FALSE,
                          eofs=1:16,lon=c(-30,20),lat=c(-20,10),
                          select=NULL,FUN="wetmean",
                          threshold=1,
                          pattern=c("tas_Amon_ens_","slp_Amon_ens_"),verbose=FALSE,nmin=nmin) {

# This function is for downscaling wet-day mean using a combination of predictors

  
  
  # Get the global mean temeprature: pentads

  # Or a combination of OLR + C.C.eq(t2m) + SLP


    # FUN: exceedance, wetfreq, wet, dry

  if (verbose) print('DSensemble.mu')

  if (verbose) print(paste('The predictor: annual',FUN))
  if (!inherits(y,'annual')) y <- annual(y,FUN=FUN,threshold=threshold,nmin=nmin)
  index(y) <- year(y)
  
  if (!is.na(attr(y,'longitude')))
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude')))
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )
  
  # Get the predictor: NCEP/NCAR
  if (verbose) print(paste("Get the set of predictors:",names(predictor),collapse=' '))
  if (is.character(predictor[[1]])) pre1 <- retrieve(ncfile=predictor[[1]],lon=lon,lat=lat)
  if (is.character(predictor[[2]])) pre2 <- retrieve(ncfile=predictor[[2]],lon=lon,lat=lat)
  if (is.character(predictor[[3]])) pre3 <- retrieve(ncfile=predictor[[3]],lon=lon,lat=lat)

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
  scorestats <- matrix(rep(NA,N*8),N,8)
  colnames(scorestats) <- c("r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration')

  flog <- file("DSensemble.precip-log.txt","at")
  dse <- list(description='DSensemble.mu')
  
  for (i in 1:N) {
    ## Need to ensure that the different predictor files match...
    print(paste(i,N,ncfiles1[select[i]],ncfiles2[select[i]],ncfiles3[select[i]]))
    gcm1 <- retrieve(ncfile = ncfiles1[select[i]],
                    lon=range(lon(PRE1))+c(-2,2),lat=range(lat(PRE1))+c(-2,2))
    gcm2 <- retrieve(ncfile = ncfiles2[select[i]],
                    lon=range(lon(PRE2))+c(-2,2),lat=range(lat(PRE2))+c(-2,2))
    gcm3 <- retrieve(ncfile = ncfiles3[select[i]],
                    lon=range(lon(PRE3))+c(-2,2),lat=range(lat(PRE3))+c(-2,2))
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
    ds <- try(DS(y,eof,eofs=eofs,verbose=verbose))
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
      attr(z,'scorestats') <- c(1-r.xval,mdiff,srati,arati,res.trend,ks,ar,ds.ratio)
      dse[[i]] <- z
      
      quality <- 100*(1-mean(scorestats[i,]))
      qcol <- quality
      qcol[qcol < 1] <- 1;qcol[qcol > 100] <- 100

      index(z) <- year(z); index(ds) <- year(ds)
      if (plot) {
        lines(z,lwd=2,col=cols[qcol])
        lines(y,type="b",pch=19)
        lines(ds,lwd=2,col="grey")
     }
      #
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(r.xval,2),
                  "R2=",round(100*sd(xval[,2])/sd(xval[,1]),2),'% ',
                  'Common EOF: bias=',round(mdiff,2),' 1- sd1/sd2=',round(srati,3),
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
                                    rcp="rcp45",biascorrect=FALSE,
                                    lon=c(-20,20),lat=c(-10,10),
                                    select=NULL,FUN="wetmean",
                                    pattern="tas_Amon_ens_",verbose=FALSE) {
  if (verbose) print('DSensemble.mu.worstcase')

  ## The predictor is based on the seasonal variations and assumes that the seasnoal cycle in the
  ## wet-day mean mu is follows a systematic dependency to the seasonal variations in the temperature
  ## - the calibration uses the Clausius Clapeiron equation to estimate the saturation water vapour
  ## rather than using the temeprature directly.
  if (verbose) print(paste('The predictor: seasonal',FUN))
  ys <- aggregate(y,by=month,FUN=FUN)

  if (!is.na(attr(y,'longitude')))
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude')))
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )

  if (verbose) print("predictor")
  if (is.character(predictor))
    pre <- spatial.avg.field(C.C.eq(retrieve(ncfile=predictor,lon=lon,lat=lat))) else
  if (inherits(predictor,'field')) pre <- predictor
  rm("predictor"); gc(reset=TRUE)
  normal61.90 <- mean(coredata(subset(pre,it=c(1961,1990))))
  x <- aggregate(pre,by=month,FUN="mean")
  cal <- data.frame(y=coredata(ys),x=coredata(x))
  stats <- cor.test(cal$y,cal$x)
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
    dev.copy2eps(file='DSensemble.mu.worstcase.cal.eps')
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

  if (plot) {
    dev.new()
    plot(aggregate(y,by=year,FUN='wetmean'),xlim=c(1900,2100),
         ylim=range(aggregate(y,by=year,FUN='wetmean'),na.rm=TRUE)*c(0.75,1.5))
    grid()
    
  }
  
  for (i in 1:N) {
    #
      gcm <- retrieve(ncfile = ncfiles[select[i]],lon=lon,lat=lat)
      gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-")
      GCM <- spatial.avg.field(C.C.eq(gcm))
      z <- annual(GCM,FUN="max")
      z <- z - mean(coredata(subset(z,it=c(1961,1990)))) + normal61.90
      i1 <- is.element(year(z),years)
      i2 <- is.element(years,year(z))
      prex <- data.frame(x=coredata(z[i1]))
      X[i,i2] <- predict(wc.model, newdata=prex) +
        rnorm(n=sum(i1),sd=max(attr(ys,'standard.error'))) 
      if (plot) lines(years,X[i,i2],col=rgb(0,0.3,0.6,0.2))
      print(paste("i=",i,"GCM=",gcmnm[i]))
    }
  if (plot) lines(aggregate(y,by=year,FUN='wetmean'),col='red',lwd=3)
  
  X <- zoo(t(X),order.by=years)
  colnames(X) <- gcmnm
  attr(X,"model_id") <- gcmnm
  #X <- attrcp(y,X)
  attr(X,'station') <- aggregate(y,by=year,FUN='wetmean')
  attr(X,'predictor') <- attr(pre,'source')
  attr(X,'domain') <- list(lon=lon,lat=lat)
  attr(X,'path') <- path
  attr(X,'scenario') <- rcp
  attr(X,'history') <- history.stamp(y)
  class(X) <- c("dsensemble","zoo")
  save(file="DSensemble.rda",X)
  print("---")
  invisible(X)   
}


DSensemble.pca <- function(y,plot=TRUE,path="CMIP5.monthly/",
                          rcp="rcp45",biascorrect=FALSE,
                          predictor="ERA40_t2m_mon.nc",
                          non.stationarity.check=FALSE,
                          eofs=1:16,lon=c(-30,20),lat=c(-20,10),
                          select=NULL,FUN="mean",rmtrend=TRUE,
                          FUNX="mean",threshold=1,
                          pattern="tas_Amon_ens_",verbose=FALSE,nmin=nmin) {

  if (verbose) print('DSensemble.pca')
  # This function is for downscaling PCA to represent a group of stations
  if (!is.na(attr(y,'longitude'))[1])
    lon <- round( range(attr(y,'longitude'),na.rm=TRUE) + lon )
  if (!is.na(attr(y,'latitude'))[1])
    lat <- round( range(attr(y,'latitude'),na.rm=TRUE) + lat )

  if (is.character(predictor))
    t2m <- retrieve(ncfile=predictor,lon=lon,lat=lat) else
  if (inherits(predictor,'field'))
    t2m <- subset(predictor,is=list(lon=lon,lat=lat))

  if (inherits(y,'season')) {
    T2M <- as.4seasons(t2m,FUN=FUNX,nmin=nmin)
    T2M <- matchdate(T2M,y)

    # Recursive: do each season seperately if there are more than one season
    if (length(table(season(y)))>1) {
      if (verbose) print('--- Apply DS to seasons seperately ---')
      Z <- list(info='DSensembe.pca for different seasons')
      for (season in names(table(season(y)))) {
        if (verbose) print(paste('Select',season))
        z <- DSensemble.pca(subset(y,it=season),plot=plot,path=path,
                            rcp=rcp,biascorrect=biascorrect,predictor=T2M,
                            non.stationarity.check=non.stationarity.check,
                            eofs=eofs,lon=lon,lat=lat,
                            select=select,FUN=FUN,rmtrend=rmtrend,
                            FUNX=FUNX,threshold=threshold,
                            pattern=pattern,verbose=verbose,nmin=nmin)
        eval(parse(text=paste('Z$',season,' <- z',sep='')))
      }
      if (verbose) print('--- Results returned as a list ---')
      return(Z)
    }

  } else if (inherits(y,'annual')) {
    T2M <- annual(t2m,FUN=FUNX,nmin=nmin)
    T2M <- matchdate(T2M,y)
  }
  
  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) select <- 1:N else
                       N <- length(select)
  if (verbose) print(ncfiles[select])

  d.y <- dim(y)
  years <- 1900:2100
  m <- length(years)
  months <- rep(month(y)[1],m)
  X <- matrix(rep(NA,N*m*d.y[2]),N,m*d.y[2])
  dim(X) <- c(d.y[2],N,m)
  gcmnm <- rep("",N)
  scorestats <- matrix(rep(NA,N*8),N,8)
  colnames(scorestats) <- c("r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration')

  #
  t <- as.Date(paste(years,months,'01',sep='-'))

  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)

  flog <- file("DSensemble.pca-log.txt","at")

  if (verbose) print("loop...") 
  for (i in 1:N) {
    gcm <- retrieve(ncfile = ncfiles[select[i]],
                          lon=range(lon(T2M))+c(-2,2),
                          lat=range(lat(T2M))+c(-2,2))
    #gcmnm[i] <- attr(gcm,'model_id')
    gcmnm[i] <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-r")
    if (inherits(y,'season')) {
      GCM <- as.4seasons(gcm,FUN=FUNX)
      GCM <- subset(gcm,it=season(T2M)[1])
    } else if (inherits(y,'annual')) {
      GCM <- annual(gcm,FUN=FUNX)

    }
    rm("gcm"); gc(reset=TRUE)
    T2MGCM <- combine(T2M,GCM)
    if (verbose) print("- - - > EOFs")
    Z <- try(EOF(T2MGCM))

        # The test lines are included to assess for non-stationarity
    if (non.stationarity.check) {
      testGCM <- subset(GCM,it=range(year(T2M))) # REB 29.04.2014
      testy <- as.station(regrid(testGCM,is=y))  # REB 29.04.2014
      attr(testGCM,'source') <- 'testGCM'        # REB 29.04.2014
      testZ <- combine(testGCM,GCM)              # REB 29.04.2014
      rm("testGCM"); gc(reset=TRUE)
    }

    if (verbose) print("- - - > DS")
    if (biascorrect) Z <- biasfix(Z)
    ds <- try(DS(y,Z,eofs=eofs,verbose=verbose))
    if (verbose) print("post-processing")
    z <- attr(ds,'appendix.1')

    # The test lines are included to assess for non-stationarity
    if (non.stationarity.check) {
        if (verbose) print('non.stationarity.check')
        testds <- DS(testy,testZ,biascorrect=biascorrect,eofs=eofs)
                                                                # REB 29.04.2014
        testz <- attr(testds,'appendix.1')                      # REB 29.04.2014
        difference.z <- testy - testz                           # REB 29.04.2014
    }
        
    i1 <- is.element(paste(years,months,sep='-'),
                     paste(year(z),month(z),sep='-'))
    i2 <- is.element(paste(year(z),month(z),sep='-'),
                     paste(years,months,sep='-'))
    
    if (verbose) print(paste('i=',i,gcmnm[i],'data points',
                             sum(i1),'=',sum(i2)))
    X[,i,i1] <- z[i2,]

    # Diagnose the residual: ACF, pdf, trend. These will together with the
    # cross-validation and the common EOF diagnostics provide a set of
    # quality indicators.
    cal <- coredata(attr(ds,"original_data"))
    fit <- coredata(attr(ds,"fitted_values"))
    if (verbose) print('examine residuals...')
    res <- as.residual(ds)
    res.trend <- 10*diff(range(trend(res)))/diff(range(year(res)))
    ks <- round(ks.test(coredata(res),pnorm)$p.value,4)
#      ar <- as.numeric(acf(trend(cal-fit,result="residual"),
#                           plot=FALSE)[[1]][2]) 
    ar <- as.numeric(acf(coredata(trend(res,result="residual")[,1]),
                         plot=FALSE)[[1]][2])
    if (verbose) print(paste("Residual trend=",
                             res.trend,'D/decade; K.S. p-val',
                             ks,'; AR(1)=',ar))

    # Evaluation: here are lots of different aspects...
    # Get the diagnostics: this is based on the analysis of common EOFs...

    xval <- attr(ds,'evaluation')
    r.xval <- round(cor(xval[,1],xval[,2]),3)
    if (verbose) print(paste("x-validation r=",r.xval))
    ds.ratio <- round(sd(ds[,1],na.rm=TRUE)/sd(y[,1],na.rm=TRUE),4)
      
    if (verbose) print(paste("sd ratio=",ds.ratio))

    #print(names(attributes(ds)))
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
        diag <- diagnose(ds,plot=FALSE)
        scorestats[i,] <- c(1-r.xval,NA,NA,NA,res.trend,ks,ar,ds.ratio)
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
      scorestats[i,] <- c(1-r.xval,mdiff,srati,arati,res.trend,ks,ar,ds.ratio)
      if (verbose) print(scorestats[i,])
      quality <- 100*(1-mean(scorestats[i,],na.rm=TRUE))
      print(paste("i=",i,"GCM=",gcmnm[i],' x-valid cor=',round(r.xval,2),
            "R2=",round(100*sd(xval[,2])/sd(xval[,1]),2),'% ',
            'Common EOF: bias=',round(mdiff,2),' 1- sd1/sd2=',round(srati,3),
            "mean=",round(mean(coredata(y),na.rm=TRUE),2),'quality=',
                    round(quality)))
  }

  if (verbose) print('Downscaling finished - need to convert PCAs to stations')

  # Unpacking the information tangled up in GCMs, PCs and stations:
  # Save GCM-wise
  gcmnm <- gsub('-','.',gcmnm)
  Y <- as.station(y)
  Z <- list(info='DSensemble.pca')
  if (verbose) print('Unscramble the results')
  for (i in 1:N) {
    pca1 <- zoo(t(X[,i,]),order.by=t)
    pca1 <- attrcp(y,pca1)
    class(pca1) <- class(y) 
    z.pca <- as.station(pca1)
    class(z.pca) <- c("dsensemble","zoo")
    attr(z.pca,'GCM') <- gcmnm[i]
    cl <- paste('Z$i',i,'_',gcmnm[i],' <- z.pca',sep='')
    eval(parse(text=cl))
    if (verbose) print(paste(i,cl))
  }

  # Save the information station-wise
  if (verbose) print('Save the results station-wise')
  x <- matrix(rep(NA,dim(X)[3]*dim(X)[2]),dim(X)[3],dim(X)[2])
  for (i in 1:(length(names(Z[[2]])))) {
    x[,] <- NA
    for (j in 1:N){
      x[,j] <- Z[[j+1]][,i]
    }
    yloc <- loc(Z[[j+1]])[i]
    yloc <- gsub('-','.',yloc)
    yloc <- gsub(' ','.',yloc)
    yloc <- gsub('/','.',yloc)
    attr(z.pca,'station') <- subset(Y,is=i)
    eval(parse(text=paste('Z$',yloc,' <- z.pca',sep='')))
  }

  #Z <- attrcp(y,Z)
  attr(Z,'predictor') <- attr(T2M,'source')
  attr(Z,'domain') <- list(lon=lon,lat=lat)
  attr(Z,'scorestats') <- scorestats
  attr(Z,'path') <- path
  attr(Z,'scenario') <- rcp
  attr(Z,'history') <- history.stamp(y)
  if (non.stationarity.check)
    attr(Z,'on.stationarity.check') <- difference.z else
    attr(Z,'on.stationarity.check') <- NULL
  attr(Z,'area.mean.expl') <- FALSE

  save(file="DSensemble.rda",Z)
  print("---")
  invisible(Z)
}
