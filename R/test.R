## Name		: test.retrieve.ncdf4
## Description	: This routine contains a series of test functions which compute the global mean 2m-temperature anomalies and other predifined regions based on both CMIP3 and CMIP5 experiments.
## Author 	: Abdelkader Mezghani, METNO
## contact 	: abdelkaderm@met.no
## Created	: 21-03-2013
## Last Updates	: 25-04-2013 ; 22-10-2013 ; 14-05-2014

## require	: ncdf4,zoo,retrieve.ncdf4,spatial.ave.field


## test <- function(path,...) UseMethod("test")

## Comparing cmip3 to cmip5 over the Arctic
plot.cmip35.global <- function(saveplot=TRUE) { 
  data(global.t2m.cmip5)
  data(global.t2m.cmip3)
  dev.new()
  fig1.zoo(z=list(z1=arctic.t2m.cmip3,z2=arctic.t2m.cmip5),select="noresm")
  if (saveplot) dev.copy2pdf(file="Global_mean_t2m_anomaly_CMIP3-5_1986-2005.pdf")
}

## Comparing cmip3 to cmip5 over Scandinavia
plot.cmip35.global <- function(saveplot=TRUE) { 
  data(scandinavia.t2m.cmip5)
  data(scandinavia.t2m.cmip3)
  dev.new()
  fig1.zoo(z=list(z1=arctic.t2m.cmip3,z2=arctic.t2m.cmip5),select="noresm")
  if (saveplot) dev.copy2pdf(file="Scandinavia_mean_t2m_anomaly_CMIP3-5_1986-2005.pdf")
}
## Comparing cmip3 to cmip5 over the Arctic
plot.cmip35.arctic <- function(saveplot=TRUE) { 
  data(arctic.t2m.cmip5)
  data(arctic.t2m.cmip3)
  dev.new()
  fig1.zoo(z=list(z1=arctic.t2m.cmip3,z2=arctic.t2m.cmip5),select="noresm")
  if (saveplot) dev.copy2pdf(file="Arctic_mean_t2m_anomaly_CMIP3-5_1986-2005.pdf")
}

## test for cmip3 over Scandinavia
test.cmip3.scandinavia <- function(...) {
  scandinavia.t2m.cmip3 <- test.retrieve.ncdf4(path="CMIP3.monthly/20c3m-sresa1b",param="tas",cntr = "Scandinavia", lon=c(-10,30),lat=c(50,70),experiment="CMIP3",climatology=c(1986,2005),fig=TRUE,saveinfile =NULL,verbose=FALSE) 
  save(scandinavia.t2m.cmip3,file="scandinavia.t2m.cmip3.rda")
}

## test for cmip5 over Scandinavia
test.cmip5.scandinavia <- function(...) {
  test.retrieve.ncdf4(path="CMIP5.monthly/rcp45/",param="tas",cntr = "Scandinavia", lon=c(-10,30),lat=c(50,70),experiment="CMIP5",climatology=c(1986,2005),fig=TRUE,saveinfile =NULL,verbose=FALSE) 
  scandinavia.t2m.cmip5 <- z
  save(scandinavia.t2m.cmip5,file="scandinavia.t2m.cmip5.rda")
}

## test cmip3 Global
test.cmip3.global <- function(...) {
  global.t2m.cmip3 <- test.retrieve.ncdf4(path="CMIP3.monthly/20c3m-sresa1b",param="tas",cntr = "Global", lon=NULL,lat=NULL,experiment="CMIP3",climatology=c(1986,2005),fig=TRUE,saveinfile =NULL,verbose=FALSE) 
  save(global.t2m.cmip3,file="global.t2m.cmip3.rda")
}

## test cmip5 Global
test.cmip5.global <- function(...) {
   global.t2m.cmip5 <- test.retrieve.ncdf4(path="CMIP5.monthly/rcp45/",param="tas",cntr = "Global", lon=NULL,lat=NULL,experiment="CMIP5",climatology=c(1986,2005),fig=TRUE,saveinfile =NULL,verbose=FALSE) 
   save(global.t2m.cmip5,file="global.t2m.cmip5.rda")
 }

## test for cmip5 over Arctic
test.cmip5.arctic <- function(...) {
  arctic.t2m.cmip5 <- test.retrieve.ncdf4(path="~/CMIP5.monthly/rcp45/",param="tas",cntr = "Arctic", lon=NULL,lat=c(60,90),experiment="CMIP5",climatology=c(1986,2005),fig=TRUE,saveinfile =TRUE,verbose=FALSE) 
  save(arctic.t2m.cmip5,file="arctic.t2m.cmip5.rda")
}

## test for cmip3 over Arctic
test.cmip3.arctic <- function(...) {
  arctic.t2m.cmip3 <- test.retrieve.ncdf4(path="CMIP3.monthly/20c3m-sresa1b",param="auto",cntr = "Arctic", lon=NULL,lat=c(60,90),experiment="CMIP3",climatology=c(1986,2005),fig=TRUE,saveinfile =TRUE,outfile=NULL) 
  save(arctic.t2m.cmip3,file="arctic.t2m.cmip3.rda")
}

## Function test # 
test.retrieve.ncdf4 <- function(path="CMIP3.monthly/20c3m-sresa1b",param="auto",cntr = "Global", lon=NULL,lat=NULL,experiment="CMIP3",climatology=c(1986,2005),fig=TRUE,saveinfile=TRUE,outfile=NULL,verbose=FALSE) {
  ncfiles <- list.files(path=path,pattern=paste("tas","_",sep=""),full.names=TRUE) # the dash is added to avoid selecting tasmin and tasmax !
  if (length(ncfiles)<1) ncfiles <- path
  
  ## Initialize
  vmodel <- rep("-",length(ncfiles))
  for (k in 1:length(ncfiles)) {
    if (substr(ncfiles[k],nchar(ncfiles[k])-1,nchar(ncfiles[k])) != "nc")
      stop("Netcdf file has not been found or the path has not been set correctly !")
    ncid <- nc_open(ncfiles[k])
    print("-------------------------------------------")
    print(ncid$filename)
    gcm <- retrieve.ncdf4(ncfile=ncid,param="auto", lon=lon,lat=lat, lev = NULL,time = NULL,ncdf.check=TRUE,miss2na = TRUE,verbose = verbose)
    project_id <- attr(gcm,"project_id")
    model_id <- attr(gcm,"model_id")
    realization <- attr(gcm,"realization")
    experiment_id <- attr(gcm,"experiment_id")
    title <- attr(gcm,"title")
    ## Compute the spatial average along lon and lat
    gcm2 <- spatial.avg.field(gcm)
    ## Compute the annual mean 
    gcm.am <- aggregate(x=gcm2,by=format(time(gcm2),"%Y"),FUN=mean)
    year <- time(gcm.am)   
    ## Compute the anomalies relative to the climatology
    agcm.ave <- gcm.am-mean(gcm.am[is.element(as.numeric(year),c(climatology[1]:climatology[2]))])
    if (k == 1) {
      X <- matrix(NA,length(year),length(ncfiles)) 
      if (fig)
        plot(agcm.ave,col="steelblue" ,xlim=c(1900,2100),ylim=c(-5,5)) 
    } else if (fig) lines(agcm.ave,col="steelblue")
    if (!is.null(model_id)) { 
      vmodel[k] <- paste(model_id,realization,sep="-")
    }
    X[,k] <- agcm.ave 
  }
  ## output format
  z <- zoo(x=X,order.by=year)
  attr(z,"country") <- cntr
  attr(z,"variable") <- "t2m"
  attr(z,"longname") <- "2m-temperature" 
  attr(z,"timeunit") <- "year"
  attr(z,"climatology") <- paste(climatology,collapse="-")
  attr(z,"model_id") <- vmodel
  attr(z,"project_id") <- project_id 
  attr(z,"experiment") <- experiment
  attr(z,"title") <- title
  attr(z,"author") <- "A. Mezghani, MET Norway , Email : abdelkaderm@met.no"
  attr(z,"Call") <- match.call()
  attr(z,"date") <- date()
  attr(z, "anomaly")        <- TRUE
  attr(z, "time:method")    <- "annual:mean"
  attr(z, "spatial:method") <- "cell:mean"
  attr(z, "title")	    <- "Spatial average of the annual 2m-temperature anomalies"
  attr(z, "URL") 	    <- "http://climexp.knmi.nl/"
  if (!is.null(saveinfile) | (is.logical(saveinfile))) {
    if (is.logical(saveinfile)) {
      outvar <- tolower(paste(cntr,param,project_id,sep="."))
      outfile <-paste(outvar,"rda",sep=".")
    }
    else if (is.character(saveinfile)) outfile <- saveinfile
    attr(z,"file") <- tolower(outfile)
    eval(parse(text=paste(outvar," <- z",sep="")))
    eval(parse(text=paste("save(",outvar,",file='",outfile,"')",sep="")))
  }
  invisible(z)
}
## e.g.
## test.retrieve.ncdf4 <- function(path="CMIP5.monthly/rcp45/",file="Norway_mean_t2m_anomaly_cmip3_1986-2005.rda")

## compute trends in Tmin
trend <- function(x) {
    tr <- rep(NA,dim(tmin)[2])
    for (i in 1:dim(tmin)[2]) {
        print(i)
        tr[i] <- as.numeric(lm(coredata(subset(tmin,is=i)) ~ index(tmin))$coefficients[[2]])
    }
    return(tr)
}

sumisnotna <- function(x) sum(!is.na(x))/length(x)

## compute scaling factor based on availability
## scale <- apply(tmin,2,sumisnotna)
## add scale in attributes
## attr(tmin,"na") <- scale
## map(tmin,FUN="trend",cex.bar=2,showall=TRUE,add.text=TRUE,showaxis=TRUE)
## dev.copy2eps(file="SHARED/data/nordic/fig6_tmin_map_nmin50.eps")
