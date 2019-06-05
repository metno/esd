## Read the sea levels for tidal gauges (EU-Circle & eSACP)
## Rasmus.Benestad@met.no, 2016-01-29, Meteorologisk institutt, Oslo, Norway

# French tidal stations: http://www.sonel.org/-Tide-gauges,29-.html?lang=en
# Daily means
#' @export
station.sonel <- function(...,urls=c('http://www.sonel.org/msl/Demerliac/VALIDATED/dCHERB.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dRSCOF.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dLCONQ.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dBREST.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dHAVRE.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dDIEPP.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dBOULO.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dCALAI.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dDUNKE.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dCONCA.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dPTUDY.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dSNAZA.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dBOUCA.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dSJLUZ.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dPBLOC.slv',
                             'http://www.sonel.org/msl/Demerliac/VALIDATED/dLROCH.slv'),
                      verbose=FALSE) {

  param='sea-level'; unit='mm'; src='sonel'; cntr='France'

  loc <- c('CHERBOURG','ROSCOFF','LE_CONQUET','BREST','LE_HAVRE','DIEPPE',
           'BOULOGNE-SUR-MER','CALAIS','DUNKERQUE','CONCARNEAU','PORT TUDY',
           'SAINT-NAZAIRE','BOUCAU-BAYONNE','SAINT JEAN-DE-LUZ','PORT_BLOC',
           'LA ROCHELLE')
  lon <- c(-1.63563001,-3.96585989,-4.78082991,-4.49483800,0.10599000,1.08444000,1.57746005,
           1.86772001,2.36664009,-3.90738010,-3.44585200,-2.20155000,-1.51482999,-1.68162300,
           -1.06157005,-1.22073600)
  lat <- c(49.65129852,48.71839905,48.35940170,48.38285000,49.48189926,49.92918000,50.72750092,
           50.96939850,51.04809952,47.87369919,47.64427400,47.26686200,43.52730179,43.39523900,
           45.56850052,46.15847800)

  for (i in 1:length(urls)) {
    if (verbose) print(paste(i,loc[i],urls[i]))
    sl <- read.table(urls[i])
    x <- zoo(sl[[2]],order.by=as.Date(sl[[1]]))
    y <- as.station(x,param=param,unit=unit,loc=loc[i],lon=lon[i],lat=lat[i],
                    src=src,cntr=cntr,url=urls[i])
    if (i==1) Y <- y else Y <- combine(Y,y)
  }
  return(Y)
}

#' @export
station.gloss <- function(...,url='https://www.psmsl.org/data/obtaining/rlr.monthly.data/rlr_monthly.zip',is=NULL,verbose=TRUE) {
  if (!file.exists('rlr_monthly.zip')) download.file(url,'rlr_monthly.zip')
  con1 <- unzip('rlr_monthly.zip', files="rlr_monthly/filelist.txt")
  meta <- read.table('rlr_monthly/filelist.txt',sep=';')
  if (verbose) print(paste('station.gloss:',dim(meta)[1],'stations'))
  if (is.null(is)) is <- 1:length(meta$V1) else {
    if (is.character(is)) {
      for (i in 1:length(is)) is[i] <- grep(toupper(is[i]),toupper(as.character(meta$V4)))
      is <- as.numeric(is)
    } else if (is.list(is)) {
      il <- is; is <- 1:length(meta$V1)
      if (!is.null(il$lon)) i1 <- meta$V3>= min(il$lon) & meta$V3<= max(il$lon) else i1 <- rep(TRUE,length(is))
      if (!is.null(il$lat)) i2 <- meta$V2>= min(il$lat) & meta$V2<= max(il$lat) else i2 <- rep(TRUE,length(is))
      #print(c(sum(i1),sum(i2)))
      is <- is[i1 & i2]
      #print(is)
      if (verbose) {print(range(meta$V3[is])); print(range(meta$V2[is]))}
    } else if (is.numeric(is))
      is <- (1:length(meta$V1))[is.element(meta$V1,is)]
  }
  if (verbose) print(paste('Reading',length(is),'stations'))
  iv <- 1
  for (i in meta$V1[is]) {
    filename <- paste0("rlr_monthly/data/",i,".rlrdata")
    con1 <- unzip('rlr_monthly.zip', files=filename)
    x <- read.table(filename,sep=';')
    x$V2[x$V2 < -999] <- NA; yr <- trunc(x$V1)
    y <- zoo(x$V2,order.by=as.Date(paste(yr,round(12*(x$V1 - yr)+0.5),'01',sep='-')))
    ii <- is.element(meta$V1,i)
    loc <- strstrip(as.character(meta[ii,4]))
    if (verbose) print(paste(iv,i,loc)); iv <- iv + 1
    y <- as.station(y,loc=loc,lon=meta[ii,3],lat=meta[ii,2],alt=0,src='GLOSS',url=url,stid=meta[ii,6],
                    param = 'sea-level',unit='mm')
    if (i==meta$V1[is][1]) Y <- y else Y <- combine.stations(Y,y)
  }
  return(Y)
}

#' @export
station.newlyn <- function(...,path='data/gloss-241_Newlyn',verbose=TRUE) {
  if (!file.exists(path)) {
    download.file('http://www.gloss-sealevel.org/extlink/https%3A//www.bodc.ac.uk/data/online_delivery/international_sea_level/gloss/ascii/g241.zip',destfile='newlyn.zip')
    dir.create(path)
    system(paste('unzip newlyn.zip d',path))
  }
  metadata <- read.table(paste(path,'ig241.txt',sep='/'),skip=5,nrows=1)
  files <- list.files(path=path,pattern='.lst',full.names=TRUE)
  for (i in 1:length(files)) {
    if (verbose) print(files[i])
    testline <- readLines(files[i],n=1)
    if (substr(testline,1,4)=='BODC') xin <- try(read.table(files[i],skip=13)) else {
      xxin <- readLines(files[i])
      keeplines <- grep('/',xxin)
      xxin <- xxin[keeplines]
      keeplines <- grep(')',xxin)
      writeLines(con='newlyn.txt',xxin[keeplines])
      xin <- try(read.table('newlyn.txt'))
    }
    if (inherits(xin,"try-error")) print('failed') else {
      yymmdd <- gsub('/','-',as.character(xin$V2))
      if (verbose) print(c(yymmdd[1],yymmdd[length(yymmdd)]))
      hr <- gsub('.',':',as.character(xin$V3),fixed=TRUE)
      if (i==1) z <- zoo(xin$V4*1000,order.by=as.POSIXlt(paste(yymmdd,hr))) else {
        zz <- zoo(xin$V4*1000,order.by=as.POSIXlt(paste(yymmdd,hr)))
        it1 <- !is.element(index(z),index(zz))
        t <- c(index(z)[it1],index(zz))
        z <- zoo(c(coredata(z),coredata(zz)),order.by=t)
      }
    }
  }
  
  z <- as.station(z,loc=metadata$V1,lon=metadata$V4,lat=metadata$V3,alt=0,src='GLOSS',
                  cntr='UK',param='sea-level',unit='mm',
                  url='http://www.gloss-sealevel.org/station_handbook/stations/241/#.VqnZSkL4phh')

  class(z) <- c('station','hour','zoo')
  return(z)

 
}
