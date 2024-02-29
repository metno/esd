#' InfoGraphics
#'
#' Various functions for visual display of data and statistics
#'
#' \code{cumugram} shows the running cumulative mean or sum of a time series
#'
#' @seealso wheel graph visprob conf vis diagram scatter plot map
#' 
#' @param x an input object of class 'station'
#' @param it A list or data.frame providing time index, e.g. month
#' @param start year and month, e.g., '-01-01' to start in january
#' @param prog a boolean; if TRUE show prognosis for end of year in cumugram
#' @param plot a boolean; if TRUE show the plot
#' @param FUN a function
#' @param verbose a boolean; if TRUE print information about progress
#' @param main main title
#' @param \dots additional arguments
#'
#' @examples
#' data(bjornholt)
#' cumugram(bjornholt)
#'
#' @export
cumugram <- function(x,it=NULL,start='-01-01',prog=FALSE,plot=TRUE,verbose=FALSE,FUN='mean',main=NULL,...) {
  stopifnot(!missing(x),inherits(x,"station"))
  
  if(verbose) print("cumugram")
  yrs <- as.numeric(rownames(table(year(x))))
  today <- Sys.Date(); yesterday <- seq(today, length.out=2, by=-1)[2]
  
  #print(yrs)
  ny <- length(yrs)
  j <- 1:ny
  col <- rgb(j/ny,abs(sin(pi*j/ny)),(1-j/ny),0.3)
  class(x) <- "zoo"
  
  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") )
      unit <- expression(degree*C) else
      unit <- attr(x,'unit')
  titletext <- paste('Running cumulative',FUN,'of')
  if (is.null(main)) 
    eval(parse(text=paste("main <- paste('",titletext,"',
                          tolower(attr(x,'longname')),sep=' ')")))
  
   if (plot) {dev.new(); par(bty="n")}
  
  z <- coredata(x)
  ylim <- c(NA,NA)
  
  #print('Find the y-range')
  y.rest <- rep(NA,ny); y2n <- y.rest
  ylim <- max(coredata(x),na.rm=TRUE) # to avoid getting warnings with empty vectors.
  md.today <- format(Sys.Date(), '-%m-%d')
  md.yesterday <- format(Sys.Date()-1, '-%m-%d')
  for (i in 1:ny) {
    #if(verbose) print(paste("Year", i, "of", ny))
    y <- window(x,start=as.Date(paste(yrs[i],start,sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    if("-02-29" %in% c(md.yesterday, md.today) & !leapyear(yrs[i])) {
      if(md.today=="-02-29") md.today <- "-03-01"
      if(md.yesterday=="02-29") {
        md.yesterday <- "-03-01"
	md.today <- "-03-02"
      }
    }
    y.rest[i] <- mean(coredata(window(x,start=as.Date(paste(yrs[i], md.today, sep='')),
                                      end=as.Date(paste(yrs[i], '-12-31', sep='')))))
    y2n[i] <- mean(coredata(window(x,end=as.Date(paste(yrs[i], md.yesterday, sep='')),
                                   start=as.Date(paste(yrs[i], '-01-01', sep='')))))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    if (FUN=='mean') z <- cumsum(coredata(y))/1:length(y) else
    if (FUN=='sum') z <- cumsum(coredata(y))
    ok <- is.finite(z)
    #rint(c(i,yrs[i],range(z[ok],na.rm=TRUE),ylim))
    ylim[!is.finite(ylim)] <- NA
    ylim[1] <- min(c(ylim,y[ok]),na.rm=TRUE)
    ylim[2] <- max(c(ylim,y[ok]),na.rm=TRUE)
  }
  #print(ylim)
  names(y2n) <- yrs
  y2n <- round(sort(y2n,decreasing=TRUE),2)
  
  if (plot) {
    par0 <- par()
    plot(c(0,length(y)),ylim,
         type="n",xlab="",
         main=main,sub=attr(x,'location'),ylab=ylab(x),...)
    grid()
  }
  cm <- rep(NA,ny)
  
  #browser()

  mm <- format(yesterday, "%m")
  dd <- format(yesterday, "%d")
  period <- paste('YYYY',start,' to YYYY-',paste(mm,dd,sep='-'),sep='')
  if (verbose) {print(yesterday); print(mm); print(dd); print(period)}
  
  if (verbose) print('No. year min max ylim[1] ylim[2]')
  Y <- matrix(rep(NA,ny*366),ny,366)
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],start,sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],start,sep='')))
    if (FUN=='mean') z <- cumsum(coredata(y))/1:length(y) else
    if (FUN=='sum') z <- cumsum(coredata(y))
    Y[i,1:length(z)] <- z
    
    if (FUN=='mean') cm[i] <- mean(coredata(window(x,
                                   start=as.Date(paste(yrs[i],start,sep='')),
                                   end=as.Date(paste(yrs[i],mm,dd,sep='-'))))) else 
                     cm[i] <- sum(coredata(window(x,
                                    start=as.Date(paste(yrs[i],start,sep='')),
                                    end=as.Date(paste(yrs[i],mm,dd,sep='-')))))
    if (plot) lines(t,z,lwd=2,col=col[i])
    if (verbose) print(c(i,yrs[i],cm[i],range(z[ok],na.rm=TRUE),ylim))
  }
  Y <- t(Y); colnames(Y) <- yrs
  if (is.null(it)) {
    if (plot) {
      lines(t,z,lwd=5,col="black")
      lines(t,z,lwd=2,col=col[i])
    }
  } else {
    y <- window(x,start=as.Date(paste(it,start,sep='')),
                    end=as.Date(paste(it,'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(it,start,sep='')))
    if (FUN=='mean') z <- cumsum(coredata(y))/1:length(y)  else
    if (FUN=='sum') z <- cumsum(coredata(y))  
    if (plot) { 
      lines(t,z,lwd=5,col="black")
      lines(t,z,lwd=2,col=col[i])
    }
  }
  tn <- t[length(t)]; 

  ## Here is some difference between monthly and daily data:
  if (!is.na(coredata(z[length(z)]))) zn <- coredata(z[length(z)]) else
                                      zn <- coredata(z[length(z)-1])
  n <- max(table(year(x)))
  if (n>=365) n <- as.numeric(diff(as.Date(c(paste(yrs[ny],start,sep=''),
                                             paste(yrs[ny],'-12-31',sep=''))))+1)
  if (n>=365) tm <- julian(as.Date('1900-12-31')) - julian(as.Date('1900-01-01')) else
              tm <- julian(as.Date('1900-12-01')) - julian(as.Date('1900-01-01'))
  #browser()
  zp <- length(z)/n * zn + (n-length(z))/n * quantile(y.rest,0.95,na.rm=TRUE)
  zm <- length(z)/n * zn + (n-length(z))/n * quantile(y.rest,0.05,na.rm=TRUE)
  zz <- length(z)/n * zn + (n-length(z))/n * mean(y.rest,na.rm=TRUE)
  if (prog & plot) {
    polygon(c(tn,rep(tm,2),tn),c(zn,zp,zm,zn),
            col=rgb(0.5,0.5,0.5,0.1),border=rgb(0.5,0.5,0.5,0.2),lwd=2)
    lines(c(tn,tm),c(zn,zz),col=rgb(0.3,0.3,0.3,0.1),lwd=3)
    text(tm,zp,round(zp,1),pos=3,cex=0.5,col='grey40')
    text(tm,zm,round(zm,1),pos=1,cex=0.5,col='grey40')
    text(tm,zz,round(zz,1),pos=4,cex=0.75)
    print(paste('Prognosis for end-of-year: ',round(zz,1),' (',round(zm,1),',',round(zp,1),')',sep=''))
  }

  if (plot) { 
    if (!is.precip(x))
      par(new=TRUE,fig=c(0.70,0.85,0.20,0.35),mar=c(0,3,0,0),
          cex.axis=0.7,yaxt="s",xaxt="n",las=1)
    else
      par(new=TRUE,fig=c(0.70,0.85,0.70,0.85),mar=c(0,3,0,0),
          cex.axis=0.7,yaxt="s",xaxt="n",las=1)
    colbar <- rbind(1:ny,1:ny)
    image(1:2,yrs,colbar,col=col)
    par(fig=par0$fig,mar=par0$mar,cex.axis=par0$cex.axis,
        yaxt=par0$yaxt,xaxt=par0$xaxt,las=par0$las)
  }
  srt <- order(cm,decreasing=TRUE)
  if (verbose) print(y2n)
  result <- cbind(yrs[srt],cm[srt])
  if (verbose) print(round(t(result)))
  colnames(result) <- c('year','cumulated')
  attr(result,'period')  <- period
  attr(result,'Y') <- Y
  invisible(result)
  
}
