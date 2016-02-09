## EU-Circle.
## This script reads the rain gauge data from the MIDAS data base (UK Met Office)
## downloaded from CEDA.
## Metadata in http://browse.ceda.ac.uk/browse/badc/ukmo-midas/metadata/CPAS/CPAS.DATA
## Data in: http://browse.ceda.ac.uk/browse/badc/ukmo-midas/data/RD/yearly_files
## File description in: http://badc.nerc.ac.uk/data/ukmo-midas/RD_Table.html
## All stations are stored in single ASCII files for each year (not an efficient way).
## The data is not open, but an application must be sent to CEDA. This script has been
## made to read downloaded ASCII files with data.
## Rasmus.Benestad@met.no, 2016-01-29, Meteorologisk institutt, Oslo, Norway

require(esd)

dealwithduplicates <- function(precip,t1,action="remove",verbose=FALSE) {
  print('dealwithduplicates')
  idup <- is.element(t1,t1[duplicated(t1)])
  if (verbose) {print('duplicated dates'); print(t1[idup])}
  if (verbose) {print('duplicated precipitation'); print(precip[idup])}
  if (action=="remove") x <- precip[!duplicated(t1)]
  if (action=="NA") {
    precip[idup] <- NA
    x <- precip[!duplicated(t1)]
  }
  if (action=="sum") {
    for (it in t1[duplicated(t1)]) {
      iii <- is.element(t1,it)
      precip[iii] <- sum(precip[iii])
    }
    x <- precip[!duplicated(t1)]
  }  
  return(x)
}

station.midas <- function(stid=NULL,loc=NULL,lon=c(-10,4),lat=c(50,60),alt=NULL,county=NULL,
                          cntr=NULL,it=NULL,nmin=30,
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
  stids <- as.character(metadata$ID); stids <- as.integer(substr(stids,1,nchar(stids)-1))
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
  ## Get the station values

  colnames <- c('ID','ID_TYPE','OB_DATE','VERSION_NUM','MET_DOMAIN_NAME','OB_END_CTIME','OB_DAY_CNT',
                 'SRC_ID','REC_ST_IND','PRCP_AMT','OB_DAY_CNT_Q','PRCP_AMT_Q','METO_STMP_TIME',
                 'MIDAS_STMP_ETIME','PRCP_AMT_J')
  
  if (verbose) print(paste('read the data:',sum(ii),'stations'))
  data.files <- list.files(path=path,pattern=pattern,full.names=TRUE)
  data.files <- data.files[1:(length(data.files)-2)] # fudge
  Y <- list() ## Set up a list object for storing the results temporarily
  for (i in 1:length(data.files)) {
    if (verbose) print(paste(i,length(data.files),data.files[i]))
    x <- try(read.table(data.files[i],sep=',',
             col.names=colnames)) # Use try - some files contain inconsistencies leading to errors
    if (inherits(x,"try-error")) print('failed') else {
      im <- is.element(x$ID,stids)              # select the stations according to given criterea
                                                # im gives both stations and different times
      i0 <- is.element(stids,x$ID)              # double check to pick only metadata for the selected stations
      s <- as.integer(rownames(table(x$ID[im]))); ns <- length(s) ## list of stations selected
      t <- as.POSIXlt(as.character(x$OB_DATE[im]))                    # time index for the selected stations
      X <- x$PRCP_AMT[im]                       # X contains all the selected stations for all available times
      time <- t[!duplicated(t)]; nt <- length(time)  ## common time index for all selected stations in this file
      if (verbose) {
        print(paste('Number of stations=',sum(im),sum(i0),'total=',length(rownames(table(x$ID)))))
        print(range(time)); print(length(time))
      }
      y <- matrix(rep(NA,ns*nt),ns,nt)          # Set up a matrix for a zoo object for the stations data
      for (i in 1:ns) {                         # Loop over all the single stations
        if (verbose) print(paste('station ID',s[i],' location=',locs[i0][i]))
        if (plot) points(lons[i0][i],lats[i0][i])
        iii <- is.element(x$ID[im],s[i]) &      # Point to all the times for one respective station
               (x$OB_DAY_CNT==1)
        precip <- X[iii]                        # Extract precipitation for single station
        info <- paste(x$ID_TYPE[iii][1],x$VERSION_NUM[iii][1],x$MET_DOMAIN_NAME[iii][1],
                      x$SRC_ID[iii][1],x$REC_ST_IND[iii][1])
        t1 <- t[iii]                            # Time index for single station
        i1 <- is.element(time,t1)               # Synchonise the times for the single station with that of the station group
        i2 <- is.element(t1,time)  
        if (verbose) print(c(sum(i1),sum(iii)))
      ## Quality check!
        ok <- TRUE
        if ( (sum(duplicated(t1))>0) | (sum(iii) != sum(i1)) ) { # Check for duplicated times for single station
          print('------- Some errors were detected ---------')
          if (sum(duplicated(t1))>0) {                           # See if it looks OK the duplicated data are removed
            print('Duplicated times!')
            #browser()
            precip <- dealwithduplicates(precip,t1,action="remove")
            tnd <- t1[!duplicated(t1)]
            i1 <- is.element(time,tnd); iii <- is.finite(precip)
            i2 <- is.element(tnd,time)
          }
          if (sum(iii) != sum(i1)) {
            print(paste('Still, not matching number of elements!',sum(iii),sum(i1)))
            print(table(x$OB_DAY_CNT[iii]))
            ok <- FALSE
            #browser()
          }
        }
        if (ok & (sum(i1)==sum(i2))) y[i,i1] <- precip[i2] else {
          print(paste('Something went wrong! sum(i1)=',sum(i1),'sum(i2)=',sum(i2)))
          #browser()
        }
      }
    }

    ## Create a time seies object
    y <- zoo(t(y),order.by=as.Date(time))
    y <- as.station(y,param='precip',unit='mm/day',
                  loc=locs[i0],lon=lons[i0],lat=lats[i0],alt=alts[i0],stid=stids[i0],
                  src='MIDAS (UK MetOffice)',ref=ref,info=info,
                  url='http://browse.ceda.ac.uk/browse/badc/ukmo-midas/data/RD/yearly_files')
    nok <- apply(coredata(y),2,'nv')
    y <- subset(y,is=(nok>=300))
    eval(parse(text=paste('Y$x.',year(time)[1],' <- y',sep='')))
  }
  invisible(Y)
}

#y <- station.midas(lon=c(-6,-2),lat=c(50,51.5))
y <- station.midas(lon=c(-5,-3),lat=c(50,51))
