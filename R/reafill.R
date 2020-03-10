#' Fill in gaps of missing data.
#' \code{reafill} is an alternative to \code{pcafill}, using ERA5 to fill in missing data and to evaluate
#' the station data based on an ordinary linear regression between data from the
#' reanalysis interpolated to the coordinates of the station data. 
#' 
#' The reanalysis is usually more complete than the station data and the final result 
#' keeps the original data wherever valid and fills in the gaps and extends the coverage
#' with informatino from the reanalysis wherever appropriate. This function is an alternative
#' to \code{pcafill} and using DS to downscale local data. Whereas \code{pcafill} is more suited 
#' for aggregated (monthly/seasonal/annual) data for a group of stations within a region with common
#' variabiliet, the \code{reafill} function is more geared to daily data. \code{pcafill} does not make
#' use of additional infromation other than assuming a stable spatio-temporal covariance structure whereas
#' \code{reafill} makes use of addtional information from reanalyses.  
#' 
#' \code{test.reafill} provides a testing routine for \code{reafill} on sample data 
#' (Oslo monthly temperature) where gaps of missing data have been introduced. The test 
#' consists of comparing with the data that has been removed before applying \code{reafill}.
#' 
#' @param x the station data with gaps that need interpolation
#' @param file Name of the reanalysis data file (netCDF). NB use daily data if x contains daily data. 
#' @param anomaly (Not yet working) subtract the mean annual cycle before interpolation and then add it back for recovering original form.
#' @param verbose Print out checks for diagnosing
#' @param plot Graphical diagnostics
#'  
#' @aliases pcafill reafill test.reafill 
#' @author R.E. Benestad
#' @export reafill
reafill <- function(x,file,anomaly=TRUE,verbose=FALSE,plot=FALSE,delta=0.3) {
  ## Read only the reanalysis data for the local region surrounding the station data 
  lon <- round(range(lon(x)) + delta*c(-1,1),3)
  lat <- round(range(lat(x)) + delta*c(-1,1),3)
  if (verbose) print(paste('reanfill',min(lon),max(lon),min(lat),max(lat)))
  Y <- retrieve(file,lon=lon,lat=lat)
  if (verbose) print(class(x))
  if (inherits(x,'month')) Y <- as.monthly(Y) else 
    if (inherits(x,'seasonal')) Y <- as.4seasons(Y) else
      if (inherits(x,'annual')) Y <- as.annual(Y)
  if (verbose) print(class(Y))
  if (verbose) print(loc(x))
  ## Extract bilienarly interpolated series from reanalysis corresponding to the station
  y <- regrid(Y,is=x)
  ## If anomaly==TRUE, then subtract the mean annual cycle before fitting. Use the climatology to 
  ## recover the original data.
  if (anomaly) {
    x0 <- x; y0 <- y
    if (inherits(y,'day')) index(y) <- as.Date(index(y))
    x <- anomaly(x); y <- anomaly(y)
  }
  if (inherits(x,'day')) index(y) <- as.Date(index(y))
  z <- y
  ## Take into acount the possibility of multiple (n) stations
  n <- dim(x)[2]; if (is.null(n)) n <- 1
  if (verbose) print(paste(n,'station(s)'))
  evaluation <- list()
  for (is in 1:n) {
    x1 <- subset(x,is=is); z1 <- subset(z, is=is)
    if (verbose) print(loc(x1))
    xz <- merge(zoo(x1),zoo(z1))
    cal <- data.frame(x=coredata(xz[,1],z=coredata(xz[,2])))
    fill <- data.frame(z=coredata(z1))
    ## fit is the interpolated data based on ordinary linear regression
    fit <- zoo(predict(lm(x ~ z, data=cal),newdata=fill),order.by=index(z1))
    evaluation[[is]] <- summary(lm(x ~ z, data=cal))
    ## combine the fitted and original data
    trange <- range(c(index(x1),index(z1)))
    if (verbose) print(trange)
    ## Common times
    tc <- c(index(x1),index(z1))
    tc[duplicated(tc)] <- NA; tc <- sort(tc[is.finite(tc)])
    ## Set up an empty time seris
    s1 <- rep(NA,length(tc))
    ## Add the interpolated data 
    if (verbose) print('Interpolated data')
    s1[is.element(tc,index(z1))] <- coredata(fit)
    ## Add the original data (overwrite interpolated data wher they overlap)
    if (verbose) print('original valid data')
    good <- is.finite(coredata(x1))
    s1[is.element(tc,index(x1))][good] <- coredata(x1)[good]
    s1 <- zoo(s1,order.by=tc)
    if (plot) {
      plot(s1,lwd=2,type='l',col='grey')
      lines(zoo(x1),col='black')
      lines(fit,col='red',lty=2)
    }
    ## The final 
    if (is > 1) z <- combine.stations(s1,fit) else z <- s1
  }
  if (anomaly) { 
    ## The climatology needs to be added as a repeated seqment
    clim <- as.climatology(x0)
    tyrs <- table(year(z))
    nyrs <- length(rownames(tyrs))
    if ((max(tyrs)==12) & (length(clim)==12)) {
      ## Monthly data
      for (im in 1:12) {
        it <- is.element(month(z),im)
        coredata(z)[it] <- coredata(z)[it] + rep(clim[im],sum(it))
      } 
      #plot(zoo(x0),lwd=3,col='grey'); lines(z,col='red')
    } else {
      ## Daily data
      wt <- 2*pi/365.25 * c(as.numeric(index(clim)))
      cal <- data.frame(y=coredata(clim),s1=sin(wt),c1=cos(wt),
                        s2=sin(2*wt),c2=cos(2*wt),s3=sin(3*wt),c3=cos(3*wt))
      wt <- 2*pi/365.25 * as.numeric(index(z))
      pre <- data.frame(s1=sin(wt),c1=cos(wt),
                        s2=sin(2*wt),c2=cos(2*wt),s3=sin(3*wt),c3=cos(3*wt))
      clim <- predict(lm(y ~c1 + s1 + c2 + s2 + c3 + s3, data=cal),newdata=pre)
      z <- z + clim
    }
    
  }
  
  attr(z,'evaluation') <- evaluation
  invisible(z) 
}

## Test the fill-in function reafill by inserting NAs into series
#' @export test.reafill
test.reafill <- function(x=NULL,file='~/data/ERA5/ERA5_t2m_mon.nc',n=3,nmiss=100,anomaly=TRUE) {
  print('test.reafill')
  if (is.null(x)) data("Oslo") else Oslo <- x
  #data("ferder"); Oslo <- ferder
  z0 <- matrix(rep(NA,n*nmiss),n,nmiss); z <- z0
  for (i in 1:n) {
    ## copy time series
    y <- subset(Oslo,it=c(1979,2018))
    ## Insert random NAs into the series
    seed <- order(rnorm(length(y)))[1:nmiss]
    cy <- coredata(y); z0[i,] <- cy[seed]; cy[seed] <- NA; c(cy) -> coredata(y)
    x <- reafill(y,file=file,verbose=TRUE,plot=FALSE,delta=0.3,anomaly=anomaly)
    z[i,] <- coredata(x)[seed]
  }
  plot(c(z0),c(z),main='Test reafill',xlab='original data',ylab='interpolated data',
       sub=paste('Correlation=',round(cor(c(z0),c(z)),3)),pch=19,col=rgb(0,0,0,0.2))
  lines(c(-100,100),c(-100,100),lty=2,col='red')
  grid()
}