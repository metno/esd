# R.E. Benestad - replacement for old ds.one to downscale the CMIP 5 ensemble
# monthly mean temperature
# 30.01.2013
# ReinClim

source("esd/R/DS.R",local=environment())
path="CMIP5.monthly/"
rcp="rcp45"
biascorrect=FALSE
plot=TRUE

  si <- list(R.version=sessionInfo()$R.version$version.string,
             esd.version=paste(sessionInfo()$otherPkgs$esd$Package,
                               sessionInfo()$otherPkgs$esd$Version,sep="_"),
             platform=sessionInfo()$platform)

  path <- paste(path,rcp,"/",sep="")
  ncfiles <- list.files(path=path,pattern="tas_Amon_ens_rcp45",full.name=TRUE)
  print(ncfiles)
  N <- length(ncfiles)
  y <- station.metno("Oslo - Blinder")
  attr(y,'aspect') <- "original"
# Need to fix the variable name:
  #attr(y,'variable') <- 't2m'
  ym <- as.4seasons(y,FUN="mean")
  ys <- as.4seasons(y,FUN="sd")
  lon <- round( attr(y,'longitude') + c(-10,30) )
  lat <- round( attr(y,'latitude') + c(-15,10) )
  t2m <- t2m.ERA40(lon=lon,lat=lat)
  T2M <- as.4seasons(t2m)

  if ( (attr(y,'unit') == "deg C") | (attr(y,'unit') == "degree Celsius") )
        unit <- expression(degree*C) else
        unit <- attr(y,'unit')

 if (plot) {
    par(bty="n")
    plot.zoo(annual(ym),type="b",pch=19,main=attr(y,'location'),
             xlab="year",ylab=unit,
             sub=paste('Station: ',attr(y,'station_id'),'; coordinates: ',
             attr(y,'longitude'),'E/',attr(y,'latitude'),'N; ',
               attr(y,'altitude'),'m.a.s.l',sep=''),
             ylim=c(-2,8)+ range(coredata(annual(ym)),na.rm=TRUE),xlim=c(1900,2100))
  }
  X <- list(station=y,domain=list(lon=lon,lat=lat),
            scenario=rcp,path=path)

  for (i in 1:N) {
    gcm <- retrieve.ncdf4(ncfile = ncfiles[i],lon=lon,lat=lat)
  # Need to fix the source attribute:
    attr(gcm,'source') <-  attr(gcm,'model_id')
    print(paste("i=",i,"GCM=",attr(gcm,'model_id')))

  # Need to fix the history attribute
    attr(gcm,'history') <- list(call="retrieve.ncdf4(ncfile = ncfiles[i],lon=lon,lat=lat)",
                                timestamp=date(),session=si)
  
    GCM <- as.4seasons(gcm)
    Z <- combine(T2M,GCM)

    print("- - - > DS")
    dsm <- DS(ym,Z,biascorrect=biascorrect)
    print("HERE - done the DS-call...")
    dss <- DS(ym,Z,biascorrect=biascorrect)
    print("HERE - done the DS-call...")
    
    eval(parse(text=paste("X$dsm.gcm.",i," <- dsm",sep="")))
    eval(parse(text=paste("X$dss.gcm.",i," <- dss",sep="")))
    if (plot) {
      lines(annual(dsm),lwd=1,col="grey")
      lines(annual(attr(dsm,'appendix.1')),lwd=2,col="steelblue")
      lines(annual(ym),type="b",pch=19)
    }
  }
  X <- attrcp(y,X)
  attr(diag,'history') <- history.stamp(y)
  class(X) <- dsensemble
  save(file="DSensemble.rda",X)
  invisible(X)

#save(file=paste("dscmip5_",attr(y,'location'),"_",N,"_rcp4.5.rda",sep=""),rcp4.5)

  
           


