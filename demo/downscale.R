print('>--- R-job for downscaling Nordic temperatures for KSS ---<')



## Function definitions:
downscale <- function(Y,predictor,it='djf',param='t2m',FUN='mean',FUNX='mean',
                      period=c(1950,2020),plot=plot,rcp='rcp45',verbose=FALSE,
                      lon=c(-20,40),lat=c(50,80),ip=1:7,n=5,
                      pattern='tas_Amon_',ds.interval=c(1950,2100),
                      rel.cord=FALSE,select=NULL,path=NULL,label=NULL) {
  
  if (verbose) print(paste('downscale',param,it,FUN,'     (',FUNX,')'))
  
  ## Use a time and space window:  
  if (verbose) print(paste('subset',paste(period,collapse='-')))
  Y <- subset(Y,it=period)
  nmin <- 360
  nmin <- switch(class(Y)[2],'month'=12,'annual'=1,'season'=3)
  
  ## Estimate seasonal means & weed out stations with little data
  if (verbose) print('season/annual')
  if(is.list(it)) {
    Y4 <- as.seasons(Y,FUN=FUN,start=it$start,end=it$end)
  } else if(is.numeric(it)) {
    Y4 <- as.seasons(Y,FUN=FUN,start=it[1],end=it[length(it)])    
  } else if ( (it!='annual') & (length(it)==1) ) {
    print(paste("Season:",it))
    Y4 <- subset(as.4seasons(Y,FUN=FUN),it=it,nmin=nmin)
  } else if(it=="annual") {
    Y4 <- as.annual(Y,FUN=FUN,nmin=nmin)
  } else Y4 <- Y
  print(table(month(Y4)))
  
  if (verbose) print('weed out missing data')
  ok <- apply(coredata(Y4),1,nv)
  Y4 <- subset(Y4,it=ok>0)
  nok <- apply(coredata(Y4),2,nv)
  Y4 <- subset(Y4,is=nok>15)
  
  print(paste(round(100*sum(!is.finite(Y4))/length(!is.finite(Y4))),'% missing ',
              ' in',paste(range(index(Y4)),collapse=' - ') ,sep=''))
  if (plot) map(Y,FUN=FUN,cex=-2)
  
  nmiss <- round(100*sum(!is.finite(Y4))/length(!is.finite(Y4)))
  
  ## Fill missing data using PCA-based regression
  if (verbose) print('pcafill')
  Z <- pcafill(Y4)
  ## Negative precipitation is impossible - clip to zero
  if (is.precip(Z)) coredata(Z)[Z<0]<- 0
  if(is.precip(Z) & FUN=='wetfreq') {
    attr(Z, "variable") <- "wetfreq"
    attr(Z, "unit") <- "fraction of days"
    attr(Z, "is.precip(Z) & FUN==longname") <- "wet-day frequency"
  } else if(is.precip(Z) & FUN=='wetfreq') {
    attr(Z, "variable") <- "wetmean"
    attr(Z, "is.precip(Z) & FUN==longname") <- "wet-day mean precipitation"
  }   
  
  if (verbose) print('pca')
  pca <- PCA(Z,n=n)
  if (plot) plot(pca)
  
  ## Downscale results
  print(paste('DSensemble',varid(predictor),it))
  
  ## Select files with data starting in the 19th or 20th century
  ## and ends at the end of the 21st century or later
  pattern <- paste0(pattern,
                    c(paste0(".*._1[8,9].*.-2", 
                             c("100","099","200","300","400","500")),
                      ".*._rcp[0-9]{2}_[0-9]{1,3}"),
                    ".nc", collapse="|")
  ## Downscale ensemble
  dse.pca <- DSensemble(pca,predictor=predictor,FUNX=FUNX,verbose=verbose,
                        biascorrect=TRUE,rcp=rcp,ip=ip,select=select,path=path,
                        ds.interval=ds.interval,nmin=1,pattern=pattern,nmin.fit=20,
                        lon=lon,lat=lat,rel.cord=rel.cord,it=it,plot=plot) #else
  attr(dse.pca,'N.missing') <- nmiss
  gcmnm <- sapply(attr(dse.pca, "model_id"), function(x) {
    gcm <- gsub(".r[0-9]{1,3}i[0-9]{1,3}p[0-9]{1,3}.*","",x)
    gsub(gcm, gsub("[.]", "_", gcm), x) } )
  dsenm <- names(dse.pca)
  i.gcm <- grep("i[0-9]{1,3}", dsenm)
  dsenm[i.gcm] <- paste0("i", seq(length(i.gcm)), "_", gcmnm)
  names(dse.pca) <- dsenm
  invisible(dse.pca)
}
