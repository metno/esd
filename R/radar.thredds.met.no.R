#' Function to read radar data from thredds.met.no
#' @aliases radar
#'
#' @param url URL for the data on thredds.met.no
#' @param lons Longitude selection - if NULL read all
#' @param lats Latitude selection - if NULL read all
#' @param param Variable name
#' @param FUN Function for daily aggregation. =NULL gives raw data
#' @param it Intex time - the years to select
#' @param verbose write out diagnostics
#' @param plot plot the results while reading. 
#' 
#' @examples 
#' \dontrun{
#' Z <- radar(lons = c(10,12), lats = c(59,61), it=2015)
#' map(Z)
#' 
#' y <- station.thredds(stid=18700,param='precip')
#' x <- regrid(Z,is=y)
#' xy <- combine.stations(subset(y,it=x),x)
#' plot(xy,new=FALSE)
#' plot(as.monthly(xy,na.rm=TRUE,FUN='mean'),new=FALSE)
#' }
#' @seealso station.thredds, meta.thredds
#' 
#' @export
radar <- function(url='https://thredds.met.no/thredds/catalog/remotesensingradaraccr/',
                  lons = c(9.5,11.5), lats = c(59,61),
                  param='lwe_precipitation_rate',FUN='sum',it=2010:2019,
                  verbose=FALSE,plot=FALSE) {
  if (verbose) print('esd::radar')
  results <- list()
  if(is.dates(it)) {
    dates.ym <- seq(as.Date(strftime(min(it),"%Y-%m-01")), as.Date(max(it)), by="month") 
    dates <- seq(as.Date(min(it)), as.Date(max(it)), by="day")
  } else {
    dates.ym <- seq(as.Date(paste0(min(it),"-01-01")), as.Date(paste0(max(it),"-12-31")), by="month")
    dates <- seq(as.Date(paste0(min(it),"-01-01")), as.Date(paste0(max(it),"-12-31")), by="day")
  }
  #for (yr in it) {
  #  for (mo in as.character(1:12)) {
  for(ym in dates.ym) {
    yr <- strftime(as.Date(ym), "%Y")
    mo <- strftime(as.Date(ym), "%m")
      #if (nchar(mo)==1) mo <- paste0('0',mo)
      contents <- readLines(paste0(url,'/',yr,'/',mo,'/catalog.html'))
      contents <- contents[grep('dataset',contents)]
      contents <- gsub("<a href='catalog.html?dataset=remotesensingradaraccr","",contents,fixed=TRUE)
      contents <- gsub("</tt></a></td>","",contents,fixed=TRUE)
      contents <- gsub("'","#",contents,fixed=TRUE)
      contents <- contents[grep('.nc',contents,fixed=TRUE)]
      #for (i in 1:length(contents)) {
      ivec <- grep(paste(strftime(dates,"%Y%m%d"),collapse="|"), contents)
      for (i in ivec) {
        eol <- regexpr("#",contents[i])[1]-1
        filename <- gsub('catalog/','dodsC/',paste0(url,substr(contents[i],1,eol)))
        if (verbose) print(filename)
        ncid <- try(nc_open(filename=filename))
        if(inherits(ncid,"try-error")) browser()
	lon <- ncvar_get(ncid,'lon')
        lat <- ncvar_get(ncid,'lat')
        dim0 <- dim(lon)
        if (verbose) print(dim0)
        if (!is.null(lons)) {
          if (verbose) print(lons)
          lon[lon > max(lons)] <- NA; lon[lon < min(lons)] <- NA
        }
        if (!is.null(lons)) {
          if (verbose) print(lats)
          lat[lat > max(lats)] <- NA; lat[lat < min(lats)] <- NA
        }
        keep <- !is.na(lon) & !is.na(lat)
        dim(keep) <- dim0
        cxl <- cumsum(rowSums(keep,na.rm=TRUE))
        cxr <- cumsum(rev(rowSums(keep,na.rm=TRUE)))
        cyl <- cumsum(colSums(keep,na.rm=TRUE))
        cyr <- cumsum(rev(colSums(keep,na.rm=TRUE)))
        x0 <- max(c(1,sum(cxl <= 0) - 1))
        y0 <- max(c(1,sum(cyl <= 0) - 1))
        xn <- max(c(1,dim0[1]-x0-sum(cxr <= 0) - 1))
        yn <- max(c(1,dim0[2]-y0-sum(cyr <= 0) - 1))
        start <- c(x0,y0,1); count<- c(xn,yn,-1)
        
        ## Test
        # image(keep)
        # contour(lon,add=TRUE)
        # contour(lat,add=TRUE,lty=2)
        # lines(cyl/max(cyl),seq(0,1,length=length(cyl)),col='blue',lwd=2)
        # lines(cyr/max(cyr),seq(1,0,length=length(cyr)),col='blue',lwd=2,lty=2)
        # lines(seq(0,1,length=length(cxl)),cxl/max(cxl),col='red',lwd=2)
        # lines(seq(1,0,length=length(cxr)),cxr/max(cxr),col='red',lwd=2,lty=2)
        # rect(x0/dim0[1],y0/dim0[2],(x0+xn)/dim0[1],(y0+yn)/dim0[2],
        #      col=rgb(0,0,0,0.1),border=rgb(0,0,0,0.2),lty=2)
        # 
        
        if (verbose) print(c(start,NA,count))
        ## Sanity check!
        if (start[1]+count[1] > dim0[1] | start[2]+count[2]>dim0[2]) {
          print('Problem detection!')
          print(start + count); print(dim0)
          browser()
        } 
          
        lon <- ncvar_get(ncid,'lon',start = start[1:2], count=count[1:2])
        lat <- ncvar_get(ncid,'lat',start = start[1:2], count=count[1:2])
        if (is.null(results[['lon']])) { 
          results[['lon']] <- lon
          if (verbose) print(summary(c(lon)))
        }
        if (is.null(results[['lat']])) { 
          results[['lat']] <- lat
          if (verbose) print(summary(c(lat)))
        }
        
        z <- ncvar_get(ncid,param,start=start,count=count)
        t <- ncvar_get(ncid,'time') # units: seconds since 1970-01-01 00:00:00 +00:00
        t <- as.POSIXlt(t,origin='1970-01-01 00:00:00')
        d <- dim(z)
        #if(is.null(d)) browser()
        #if(length(d)<3) browser()
        ## If FUN is provided, do the daily aggregation
        if (!is.null(FUN)) {
          t <- as.Date(t[1])
          dim(z) <- c(d[1]*d[2],d[3])
          z <- apply(z,1,FUN=FUN)
          dim(z) <- c(d[1],d[2])
        }
        if (verbose) print(t)
        if (is.null(results[['lat']])) 
          results[['time']] <- t else results[['time']] <- c(results[['time']],t)
        results[[paste(param,t[1],sep='.')]] <- z
        if (verbose) print(dim(results[[paste(param,t[1],sep='.')]]))
        if (plot) {
          if (!is.null(FUN)) image(z) else image(z[,,1])
          contour(lon,add=TRUE)
          contour(lat,add=TRUE,lty=2)
        }
        nc_close(ncid)
        if (!verbose) cat('.')
      #}
    }
  }
  if (verbose) print('Read all the data')
  ## Extract and remove the coordinates from the list
  lons <- c(results[['lon']]); results[['lon']] <- NULL
  lats <- c(results[['lat']]); results[['lat']] <- NULL
  time <- results[['time']]; results[['time']] <- NULL
  if (!is.null(FUN)) time <- as.Date(time) else time <- as.POSIXlt(time)
  if (verbose) str(results)
  
  ## Ensure that the times match the elements
  iit <- names(results)
  nct <- nchar(iit)
  iit <- as.Date(substr(iit,nct-9,nct))
  time <- time[is.element(time,iit)]
  time <- time[!duplicated(time)]
  ## convert the elements in list to matrix
  nt <- length(time); nr <- length(results[[1]])
  #Z <- do.call(rbind,lapply(results,matrix,ncol=length(lon),byrow=TRUE))
  radarZ <- matrix(rep(NA,nt*nr),nt,nr)
  if (nt != length(results)) {
    print('Potential problem found - number of elements != length of time')
    str(results)
    print(nt)
    browser()
  }
  for (i in 1:nt) radarZ[i,] <- c(results[[i]])
  if (verbose) str(results)
  radarZ <- zoo(radarZ,order.by=time)
  attr(radarZ,'longitude') <- t(lons)
  attr(radarZ,'latitude') <- t(lats)
  attr(radarZ,'url') <- url
  attr(radarZ,'variable') <- param
  attr(radarZ,'unit') <- 'mm'
  attr(radarZ,'aspect') <- 'radar reflectivity'
  attr(radarZ,'greenwich') <- TRUE
  attr(radarZ,'source') <- 'The Norwegian Meteorological Institute'
  class(radarZ) <- c('station','day','zoo')
  radarZ <- as.field(radarZ)
  invisible(radarZ)
}