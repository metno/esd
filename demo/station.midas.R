## EU-Circle.
## This script reads the rain gauge data from the MIDAS data base (UK Met Office)
## downloaded from CEDA.
## Metadata in http://browse.ceda.ac.uk/browse/badc/ukmo-midas/metadata/CPAS/CPAS.DATA
## Data in: http://browse.ceda.ac.uk/browse/badc/ukmo-midas/data/RD/yearly_files
## All stations are stored in single ASCII files for each year (not an efficient way).
## The data is not open, but an application must be sent to CEDA. This script has been
## made to read downloaded ASCII files with data.
## Rasmus.Benestad@met.no, 2016-01-29, Meteorologisk institutt, Oslo, Norway

require(esd)

station.midas <- function(stid=NULL,loc=NULL,lon=c(-10,4),lat=c(50,60),alt=NULL,county=NULL,cntr=NULL,it=NULL,nmin=30,
                          path='data/midas/',pattern='midas_raindrnl',metaid='CPAS.DATA',verbose=TRUE,plot=TRUE) {
  metacolnames <- c('SRC_ID','SRC_NAME','ID_TYPE','ID','MET_DOMAIN_NAME','SRC_CAP_BGN_DATE','SRC_CAP_END_DATE',
                    'PRIME_CAPABILITY_FLAG','RCPT_METHOD_NAME','DB_SEGMENT_NAME','DATA_RETENTION_PERIOD',
                    'LOC_GEOG_AREA_ID','POST_CODE','WMO_REGION_CODE','HIGH_PRCN_LAT','HIGH_PRCN_LON',
                    'GRID_REF_TYPE','EAST_GRID_REF','NORTH_GRID_REF','ELEVATION','HYDR_AREA_ID',
                    'DRAINAGE_STREAM_ID','ZONE_TIME','SRC_BGN_DATE','SRC_END_DATE')

  ## Reading the meta-data is tricky because of inconsistent use of commas.
  ## Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  : 
  ## line 79261 did not have 25 elements

  if (verbose) print('station.midas')
  ref <- paste('Met Office (2012): Met Office Integrated Data Archive System (MIDAS) Land and Marine Surface Stations Data (1853-current). ',
               'NCAS British Atmospheric Data Centre, date of citation. http://catalogue.ceda.ac.uk/uuid/220a65615218d5c9cc9e4785a3234bd0')
  meta.file <- list.files(path=path,pattern=metaid,full.names=TRUE)
  metatest <- readLines(meta.file,n=1)
  metawidths=gregexpr(',',metatest)[[1]]
  metawidths <- c(metawidths[1],diff(metawidths))
  metawidths <- c(metawidths,nchar(metatest) - sum(metawidths))
  metadata <- read.fwf(meta.file,widths=metawidths,col.names=metacolnames,nrows=79000)
  stids <- as.character(metadata$SRC_ID); stids <- as.integer(substr(stids,1,nchar(stids)-1))
  locs <- as.character(metadata$SRC_NAME); locs <- substr(locs,1,nchar(locs)-1); locs <- gsub("  ","",locs)
  alts <- as.character(metadata$ELEVATION); alts <- as.numeric(substr(alts,1,nchar(alts)-1))
  lats <- as.character(metadata$HIGH_PRCN_LAT); lats <- as.numeric(substr(lats,1,nchar(lats)-1))
  lons <- as.character(metadata$HIGH_PRCN_LON); lons <- as.numeric(substr(lons,1,nchar(lons)-1))
  cntrs <- as.character(metadata$LOC_GEOG_AREA_ID); cntrs <- substr(cntr,1,nchar(cntrs)-1); cntrs <- gsub("  ","",cntrs)
  t1 <- as.character(metadata$SRC_CAP_BGN_DATE); t1 <- substr(t1,1,nchar(t1)-1); t1 <- as.Date(t1)
  t2 <- as.character(metadata$SRC_CAP_END_DATE); t2 <- substr(t2,1,nchar(t2)-1); t2 <- as.Date(t2)
  nyrs <- year(t2) - year(t1) + 1 
  ii <- is.finite(stids)
  if (!is.null(stid)) i1 <- is.element(stids,stid) else i1 <- ii
  if (!is.null(loc)) i2 <- is.element(substr(toupper(locs),1,nchar(loc)),toupper(loc)) else i2 <- ii
  if (!is.null(alt)) {
    if (alt >= 0) i3 <- (alts >= alt) else i3 <- (alts < alt)
  } else i3 <- ii
  if (!is.null(lat)) i4 <- (lats <= max(lat)) & (lats >= min(lat))
  if (!is.null(lon)) i5 <- (lons <= max(lon)) & (lons >= min(lon))
  if (!is.null(cntr)) i6 <- is.element(substr(toupper(cntrs),1,nchar(cntr)),toupper(cntr)) else i6 <- ii
  if (!is.null(nmin)) i7 <-(nmin >= nyrs) else i7 <- ii
  if (!is.null(it)) i8 <- ( (year(t1) <= min(it)) & (year(t2) >= max(it)) ) else i8 <- ii

  if (verbose) print(paste('#maching ID=',sum(i1),' #maching loc=',sum(i2),' #matching alt=',sum(i3),' #machning lat=',sum(i4),
                           ' #matchin lon=',sum(i5),' #matching coyntry=',sum(i6),' #matching length=',sum(i7),' #matching interval=',sum(i8)))
  ii <- i1 & i2 & i3 & i4 & i5 & i6 & i7 & i8
  stids <- stids[ii]
  locs <- locs[ii]
  lats <- lats[ii]
  lons <- lons[ii]
  alts <- alts[ii]
  cntrs <- cntrs[ii]
  t1 <- t1[ii]
  t2 <- t2[ii]
  nyrs <- nyrs[ii]
  if ( (plot) & (sum(ii)>0) ) {
    nt <- nyrs; nt[nt > 50] <- 50
    z <- abs(alts/max(alts,na.rm=TRUE))
    z[!is.finite(z)] <- 0
    colt <- rgb(1-z,1,1-z,nt/50)
    plot(lons,lats,cex=0.5,col='grey')
    points(lons,lats,col=colt,pch=19,cex=0.5)
    data(geoborders)
    lines(geoborders)
  }
  
  if (sum(ii)==0) { print('station.midas: could not find any station'); return(NULL) }
  ##

  if (verbose) print(paste('read the data:',sum(ii),'stations'))
  data.files <- list.files(path=path,pattern=pattern,full.names=TRUE)
  Y <- list()
  for (i in 1:length(data.files)) {
    if (verbose) print(paste(i,length(data.files),data.files[i]))
    x <- try(read.table(data.files[i],sep=','))
    if (inherits(x,"try-error")) print('failed') else {
      im <- is.element(x$V1,stids)
      i0 <- is.element(stids,x$V1)
      s <- as.integer(rownames(table(x$V1[im]))); ns <- length(s)
      t <- as.Date(x$V3[im])
      X <- x$V10[im]
      time <- t[!duplicated(t)]; nt <- length(time)
      if (verbose) {print(range(time)); print(length(time))}
      y <- matrix(rep(NA,ns*nt),ns,nt)
      for (i in 1:ns) {
        if (verbose) print(paste('station ID',s[i],' location=',locs[i0][i]))
        if (plot) points(lons[i0][i],lats[i0][i])
        iii <- is.element(x$V1[im],s[i])
        precip <- X[iii]
        i1 <- is.element(time,t[iii])
        if (verbose) print(c(sum(i1),sum(iii)))
      ## Quality check!
        ok <- TRUE
        if ( (sum(duplicated(t[iii]))>0) | (sum(iii) != sum(i1)) ) {
          print('------- Some errors were detected ---------')
          if (sum(duplicated(t[iii]))>0) {
            idup <- is.element(t[iii],t[iii][duplicated(t[iii])])
            print('duplicated dates'); print(t[iii][idup])
            print('duplicated precipitation'); print(precip[idup])
            precip <- precip[!duplicated(t[iii])]
            tnd <- t[iii][!duplicated(t[iii])]
            i1 <- is.element(time,tnd)
          }
          if (sum(iii) != sum(i1)) {
            print(paste('Still, not matching number of elements!',sum(iii),sum(i1)))
          ok <- FALSE
          }
        }
        if (ok) y[i,i1] <- precip
      }
    }

    ## Create a time seies object
    y <- zoo(t(y),order.by=time)
    y <- as.station(y,param='precip',unit='mm/day',
                  loc=locs[i0],lon=lons[i0],lat=lats[i0],alt=alts[i0],stid=stids[i0],
                  src='MIDAS (UK MetOffice)',ref=ref,
                  url='http://browse.ceda.ac.uk/browse/badc/ukmo-midas/data/RD/yearly_files')
    eval(parse(text=paste('Y$x.',year(time)[1],' <- y',sep='')))
  }

}

y <- station.midas(lon=c(-6,-2),lat=c(50,51.5))
