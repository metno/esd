#' 
#' A function that lists files on thredds catalogue.
#' @seealso retrieve
#'  
#' @param caturl url of the thredds ctalogue containing the files
#' @param extension the extension ot the netCDF files (default is ".nc")
#' @param pattern text pattern to filter file names
#' @param verbose if true promt test lines
#'
#' @examples
#' 
#' list_thredds()
#' list_thredds(pattern=c("rcp45","r1i1p1","RACMO"))
#' list_thredds(caturl="http://esgf3.dkrz.de/thredds/catalog/esgcet/12/CMIP6.ScenarioMIP.DKRZ.MPI-ESM1-2-HR.ssp585.r1i1p1f1.day.pr.gn.v20190710.html")
#' caturl <- "https://thredds.met.no/thredds/catalog/KSS/Klima_i_Norge_2100/bias_corrected/3DBC/cross-validation/noresm-r1i1p1-remo/tasmin/catalog.html"
#' extension <- '.nc4'
#' list_thredds(caturl,extension)
#' @export

list_thredds <- function(caturl="https://thredds.met.no/thredds/catalog/KSS/Klima_i_Norge_2100/seasonal_RCM/catalog.html", 
                         extension=".nc", pattern = NULL, verbose=FALSE) {
  
  if (verbose) print(paste('list_thredds',caturl))
  
  URL <- url(caturl)
  readLines(URL) -> xmlcode
  if (verbose) print(xmlcode)
  close(URL)
  files<- xmlcode[grep(extension,xmlcode,fixed=TRUE)]
  files <- gsub('.*<tt>','',files)
  files <- gsub('.*<code>','',files)
  files <- gsub('</tt></a></td>','',files)
  files <- gsub('</code></a>','',files)
  ## Remove rubbish
  stupid <- c("<!DOCTYPE html PUBLIC '-//W3C//DTD HTML 4.01 Transitional//EN'","Surface Observations Archive",
              "<tr>","&nbsp;<tt><td>","--<tt><td>","<table>","<body>","<html>","<head>",
              "</tr>","&nbsp;</tt></td>","--</tt></td>","</table>","</body>","</html>","</head>",
              "<th align","http://www.w3.org","<meta http-equiv","Catalog https://thredds.met.no",
              "Get messages","If you have questions or problems","<HR size","<a href=","<link rel")
  for (stu in stupid) if (length(grep(stu,files)>0)) files <- files[-grep(stu,files)]
  if (!is.null(pattern)) for (i in 1:length(pattern)) files <- files[grep(pattern[i],files,ignore.case = TRUE)]
  attr(files,'caturl') <- caturl
  attr(files,'pattern') <- pattern
  attr(files,'history') <- history.stamp()
  return(files)
}


#' A function that provides metadata for S-ENDA data on thredds.
#' @seealso retrieve list_thredds station.s_enda
#'  
#' @param caturl url of the thredds ctalogue containing the files
#' @param verbose if true promt test lines
#'
#' @examples
#' 
#' meta <- meta.s_enda()
#' @export
meta.s_enda <- function(caturl="https://thredds.met.no/thredds/catalog/met.no/observations/surface/catalog.html",verbose=FALSE) {
  if (verbose) print('meta.s-enda')
  stids <- list_thredds(caturl=caturl,extension = '/')
  if (verbose) print(stids)
  ns <- length(stids); lons <- rep(NA,ns); lats <-lons; 
  locs <- rep("?",ns); longnames <- locs; srcs <- locs; params <- locs; urls <- locs; 
  
  for (is in 1:ns) {
    url <- sub("surface/catalog.html",paste0("surface/",stids[is],"catalog.html"),caturl)
    if (verbose) print(url)
    vars <- list_thredds(caturl=url)
    if (verbose) print(vars)
    if (length(vars)>0) { 
      urlf <- sub("catalog.html","",paste0(url,vars[1]))
      ## https://thredds.met.no/thredds/catalog/met.no/observations/surface/99950/wind_speed_10m_st_99950.nc
      ## https://thredds.met.no/thredds/dodsC/met.no/observations/surface/99950/wind_speed_10m_st_99950.nc.html
      urlf <- sub('catalog','dodsC',urlf)
      if (verbose) print(urlf)
      ncid <- nc_open(urlf)
      lon <- ncvar_get(ncid,'longitude')
      lat <- ncvar_get(ncid,'latitude')
      #if (verbose) print(ncid$)
      loc <- ncatt_get(ncid,0,'platform')$value
      lonm <- ncatt_get(ncid,0,'title')$value
      src <- paste(ncatt_get(ncid,0,'source')$value,ncatt_get(ncid,0,'institution')$value)
      nc_close(ncid)
      lons[is] <- lon; lats[is] <- lat; urls[is] <- urlf; locs[is] <- loc; longnames[is] <- lonm
      srcs[is] <- src
      l2 <- regexpr('_st_',vars)-1; l1 <- rep(1,length(vars))
      if (length(vars)>1) params[is] <- paste(substr(vars,l1,l2),collapse=', ') else params[is] <- substr(vars,l1,l2)
    }
  }
  
  ## att <- c("station_id","location","country","longitude","latitude","altitude","element","start","end",
  ## "source","wmo","quality")
  alts <- rep(NA,ns); start <- alts; end <- alts; wmo <- alts
  ok <- is.finite(lons) & is.finite(lats)
  if (verbose) print(sum(ok))
  meta <- data.frame(station_id = sub('/','',stids)[ok],location=locs[ok],longitude=lons[ok],latitude=lats[ok],
                     altitude=alts[ok],element=params[ok],start=start[ok],end=end[ok], source=srcs[ok], 
                     wmo=wmo[ok],url=urls[ok])
  class(meta) <- c('stationmeta','S-ENDA','data.frame')
  return(meta)
}

#' A function that provides reads S-ENDA data from thredds.
#' @seealso retrieve list_thredds meta.s_enda
#'  
#' @param stid station ID or the S-ENDA metadata object
#' @param param parameter ('precip'|'t2m')
#' @param freq frequency of the data 
#' @param FUN function for aggregating data (e.g. sub-daily for generating daily data)
#' @param start.precip Hour of daily start of the accumulated rainfall
#' @param caturl url of the thredds ctalogue containing the files
#' @param verbose if true promt test lines
#'
#' @examples
#' 
#' x <- station.s_enda()
#' x <- station.s_enda(param='t2m')
#' ## Using the S-ENDA metadata (it takes some time to fetch themetadata...)
#' meta <- meta.s_enda()
#' x <- station.s_enda(meta[1:2,])
#' @export
station.s_enda <- function(stid=18700,param='precip',freq='day',FUN='default',start.precip=7,...,
                           caturl="https://thredds.met.no/thredds/dodsC/met.no/observations/surface/",verbose=FALSE) {
  if (verbose) print('station.s_enda')
  if (inherits(stid,'S-ENDA')) { 
    if (verbose) print('Using S-ENDA metadata')
    #param <- stid$element
    stid <- stid$station_id
    if (verbose) print(table(param))
  }
  ## Get thevariable names
  ns <- length(param)
  vars <- rep('?',ns)
  for (is in 1:ns) { 
    vars[is] <- switch(param[is],'precip'='precipitation_amount','t2m'='air_temperature_2m','fx'='wind_speed_10m')
    if (is.na(vars[is])) vars[is] <- param
  }
  
  ## get station data based on station_id
  urls <- paste0(caturl,stid,'/',vars,'_st_',stid,'.nc')
  if (verbose) print(urls)
  ns <- length(urls)
  lons <- rep(NA,ns); lats <- lons; alts <- lons; 
  locs <- rep('?',ns); units <- locs; params <- locs; lonms <- locs; srcs <- locs
  
  
  is <- 0
  for (url in urls) {
    ncid <- try(nc_open(url))
    if (!inherits(ncid,'try-error')) { 
      is <- is + 1
      x <- ncvar_get(ncid,'obs')
      units[is] <- ncatt_get(ncid,'obs','units')$value
      t <- ncvar_get(ncid,'time')
      torg <- ncatt_get(ncid,'time','units')$value
      torg <- sub('seconds','days',torg)
      tim <- as.POSIXlt(t,origin=sub('days since','',torg))
      lons[is] <- ncvar_get(ncid,'longitude')
      lats[is] <- ncvar_get(ncid,'latitude')
      #if (verbose) print(ncid$)
      locs[is] <- ncatt_get(ncid,0,'platform')$value
      if (verbose) print(locs[is])
      lonms[is] <- ncatt_get(ncid,0,'title')$value
      srcs[is] <- paste(ncatt_get(ncid,0,'source')$value,ncatt_get(ncid,0,'institution')$value)
      nc_close(ncid)
      X <- zoo(x,order.by=tim)
      if (freq=='day') time <- as.Date(tim) else 
        if (freq=='hour') time <- round(tim/3600)*3600 else
          if (freq=='month') time <- as.yearmon(tim)
      if (verbose) print(time[1:50])
      if (FUN=='default') {
        if (param=='precip') {
          ## The standard observations start 07:00 in the morning
          FUN <- 'sum' 
          index(X) <- index(X) - start.precip*3600
          if (verbose) print(c(length(X),length(tim),length(t),length(x),length(time)))
        } else FUN <- 'mean'
      }
      if (is==1) Y <- aggregate(X,by=time,FUN=FUN) else
        Y <- merge(Y,aggregate(X,by=time,FUN=FUN),...)
    }
  }
  Y <- as.station(Y,param=param,unit=units,lon=lons,lat=lats,alt=alts,cntr=rep('Norway',ns),
                  longname=lonms,stid=stid,url=urls)
  return(Y)
}



