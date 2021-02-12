#' Simple and handy functions
#'
#' \code{attrcp(x,y)} passes on attributes from \code{x} to \code{y} and
#' returns the \code{y} with new attributes.
#'
#' \code{ensemblemean} returns the ensemble mean for dsensemble objects. The argument \code{FUN} can also be
#' used to estimate other statistics, e.g \code{FUN='q95'} where
#' \code{q95=function(x) apply(x,1,quantile,probs=0.95)}
#'
#' \code{TGW} uses triangulation of pressure measurements to estimate geostrophic wind based on
#' Alexandersson et al. (1998), Glob. Atm. and Oce. Sys.
#'
#' \code{stand} gives a standardised time series.
#'
#' \code{zeros} counts the occurrence of zero values in a vector.
#'
#' \code{nv} returns the number of valid points.
#'
#' \code{missval} computes the percentage of missing data.
#'
#' \code{arec} compares the number of record-breaking events to the number of
#' events for a time series of iid data (\code{sum(1/1:n)})
#' 
#' \code{strstrip} strips off spaces before and after text in a string.
#' 
#' \code{as.decimal} converts between degree-minute-second into decimal value.
#' 
#' \code{cv} computes the coefficient of variation.
#' 
#' \code{nv} count the number of valid data points.
#' 
#' \code{q5}, \code{q95}, \code{q975} and \code{q995} are shortcuts to the 5\%, 95\%, 97.5\% and 99.5\% percentiles.
#' 
#' \code{trend.coef} and \code{trend.pval} return the coefficient and the p-value of the linear trend.
#' 
#' \code{exit} is a handy function for exiting the R session without saving.
#' 
#' \code{figlab} is a handy function for labelling figures (e.g. 'Figure 1')
#' 
#' \code{ndig} estimates the number of digits for round(x,ndig), e.g. in scales for plotting.
#' 
#' \code{factor2numeric} transforms a factor to a numeric object
#' 
#' \code{rmse} and \code{RMSE} calculate the root-mean-square error
#' 
#' @aliases as.decimal nv cv q5 q95 q975 q995 filt
#' filt.default exit figlab ndig ensemblemean propchange stand rmse RMSE
#' firstyear lastyear eofvar test.num.predictors arec
#' arec.default arec.station lastrains lastdry lastelementrecord strstrip bin
#' factor2numeric zeros missval
#' @seealso attrcp
#'
#' @importFrom stats quantile qgeom qpois dnorm filter 
#'
#' @param x A data.frame or a coredata zoo object.
#' @param na.rm If TRUE, remove NA's from data
#' @param type 'ma' for moving average (box-car), 'gauss' for Gaussian, 'binom'
#' for binomial filter,' parzen' for Parzen filter, 'hanning' for Hanning
#' filter, or 'welch' for Welch filter.
#' @param lowpass True for smoothing, otherwise the highpass results is
#' returned
#' @param triangle a group of three stations with sea-level pressure', e.g.
#' from ECA\&D.
#' @param nbins number of bins/categories
#' @return \item{as.decimal }{Decimal value} \item{trend.coef }{Linear trend
#' per decade}
#' @author A. Mezghani
#' @keywords rtools ~kwd2
#' @examples
#' 
#' ## Monthly mean temperature at 'Oslo - Blindern' station 
#' data(Oslo)
#' ## Compute the linear trend and the p-value on annual aggregated values 
#' tr <- trend.coef(coredata(annual(Oslo)))
#' pval <- trend.pval(coredata(annual(Oslo)))
#' \dontrun{
#' pp <- station(param='slp',cntr='Denmark',src='ecad')
#' wind <- TGW(subset(pp,is=c(1,3,10))
#' plot(wind)
#' ws <- sqrt(wind[,1]^2 + wind[,2]^2)
#' plot(ws)
#' hist(ws)
#' 
#' ## Estimate wind for a larger group of stations
#' wind <- geostrophicwind(pp,nmax=10)
#' u <- subset(wind,param='u')
#' v <- subset(wind,param='u')
#' ws <- sqrt(u^2+v^2)
#' ws <- attrcp(v,ws)
#' class(ws) <- class(v)
#' attr(ws,'variable')='windspeed'
#' attr(ws,'longname')='geostrophic wind speed'
#' map(ws,FUN='quantile',probs=0.98)
#' 
#' ## Test firstyears on HadCRUT4
#' if (!file.exists('~/Downloads/HadCRUT.4.6.0.0.median.nc')) {
#'   print('Download HadCRUT4')
#'   download.file('https://crudata.uea.ac.uk/cru/data/temperature/HadCRUT.4.6.0.0.median.nc',
#'                 dest='~/Downloads/HadCRUT.4.6.0.0.median.nc') 
#' }
#' 
#' Obs <- annual(retrieve('~/Downloads/HadCRUT.4.6.0.0.median.nc',param='temperature_anomaly'))
#' lons <- rep(lon(Obs),length(lat(Obs)))
#' lats <- sort(rep(lat(Obs),length(lon(Obs))))
#' fy <- firstyear(Obs)
#' map(subset(Obs,it=1))
#' points(lons[fy==1850],lats[fy==1850])
#' map(Obs,FUN='firstyear')
#' }
#'
#' @export as.decimal
as.decimal <- function(x=NULL) {
    ## converts from degree min sec format to degrees ...
    ##x is in the form "49 deg 17' 38''"
  if (!is.null(x)) {
        deg <- as.numeric(substr(x,1,2)) 
        min <- as.numeric(substr(x,4,5))
        sec <- as.numeric(substr(x,7,8))     
        x <- deg + min/60 + sec/3600
    }
    return(x)
}

#' @export
nv <- function(x) sum(is.finite(x))

#' @export
ndig <- function(x) {
  i0 <- (x==0) & !is.finite(x)
  if (sum(i0)>0) x[i0] <- 1
  y <- -trunc(log(abs(x))/log(10))
  if (sum(i0)>0) y[i0] <- 0
  return(y)
}

#' @export
strstrip <- function(x) {
  if (is.na(x)) return(NA)
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(NA)
  while (substr(x,1,1)==' ') x <- substr(x,2,nchar(x))
  while (substr(x,nchar(x),nchar(x))==' ') x <- substr(x,1,nchar(x)-1)
  return(x)
}

#' @export
eofvar <- function(x) {
  if (inherits(x,c('eof','pca'))) {
    attr(x,'eigenvalues')^2/attr(x,'tot.var')*100
  } else {
    NULL
  }
}

#' @export
firstyear <- function(x,na.rm=FALSE,verbose=FALSE) {
  if (verbose) print('firstyear')
  yrs <- year(x)
  if (verbose) print(range(as.numeric(yrs)))
  if (is.null(dim(x))) y <- min(yrs[is.finite(x)]) else { 
    nv <- apply(x,2,'nv')
    y <- rep(NA,length(nv))
    ok <- (1:length(nv))[nv > 0]
    y[ok] <- apply(x[,ok],2,function(x,yrs=yrs) min(yrs[is.finite(x)]),yrs)
    #for (i in ok) y[i] <- min(yrs[is.finite(x[,i])])
  }
  y[!is.finite(y)] <- NA
  if (verbose) print(table(as.numeric(y)))
  return(y)
}

#' @export
lastyear <- function(x,na.rm=FALSE,verbose=FALSE) {
  if (verbose) print('lastyear')
  yrs <- year(x)
  if (verbose) print(range(as.numeric(yrs)))
  if (is.null(dim(x))) y <- max(yrs[is.finite(x)]) else { 
    nv <- apply(x,2,'nv')
    y <- rep(NA,length(nv))
    ok <- (1:length(nv))[nv > 0]
    y[ok] <- apply(x[,ok],2,function(x,yrs=yrs) max(yrs[is.finite(x)]),yrs)
  }     
  y[!is.finite(y)] <- NA
  if (verbose) print(table(as.numeric(y)))
  return(y)
}

#' @export difftime.month
difftime.month <- function(enddate, startdate, verbose=FALSE) {
  if(verbose) print("difftime.month")
  yend <- as.integer(strftime(enddate,"%Y"))
  ystart <- as.integer(strftime(startdate,"%Y"))
  mend <- as.integer(strftime(enddate,"%m"))
  mstart <- as.integer(strftime(startdate,"%m"))
  return((yend-ystart)*12 + mend-mstart)
}

## Iterate using n number of predictands in the downscaling and retrive the cross-val given the number of predictands
#' @export test.num.predictors
test.num.predictors <- function(x=NA,y=NA,nmax.x=6,nmin.x=3,nmax.y=4,nam.x='NA', nam.y.res='NA', nam.y='NA', nam.x.dom='NA',nam.t='NA',verbose=FALSE) {
  predictor_field <- x
  predictand_field <- y
  max_EOFs_predictor <- nmax.x
  min_EOFs_predictor <- nmin.x
  max_EOFs_predictand <- nmax.y
  min_EOFs_predictand <- 1
  if (verbose) cat('The input data should be a field. \n The field is iteratively downscaling with a different number of predictor EOFs. \n A dataframe with a summary of cross-validation correlation coefficients (Q2) and, if provided, meta-data is returned.\n')
  stopifnot (inherits(predictor_field,'field') & inherits(predictor_field,'zoo'))
  stopifnot (inherits(predictand_field,'field') & inherits(predictand_field,'zoo')) #add not eof not ds
  
  training<-setNames(data.frame(matrix(ncol =max_EOFs_predictand-min_EOFs_predictand+2, nrow = max_EOFs_predictor-min_EOFs_predictor+2), stringsAsFactors = FALSE), c(paste('y.EOF',seq(min_EOFs_predictand,max_EOFs_predictand)),"wQ2_lim0.2"))
  rownames(training)<- c('y R2:', paste(seq(min_EOFs_predictor,max_EOFs_predictor),'x.EOFs Q2:'))

  n_predictand_EOFs <- max_EOFs_predictand 
  if (verbose) cat('\n Downscaling',n_predictand_EOFs,'EOFs from the predictand field. \n',fill = TRUE)
  predictand <- EOF(predictand_field,n=max_EOFs_predictand)
  training[1,1:(ncol(training)-1)]<-round(eofvar(predictand),3)
  
  for (n_predictor_EOFs in rev(seq(min_EOFs_predictor,max_EOFs_predictor))) {
    if (verbose) cat('\n Using',n_predictor_EOFs,'EOFs from the coarse input as predictors.\n',fill = TRUE)
    predictor <- EOF(predictor_field,n=n_predictor_EOFs)
    ds_obj <- DS(predictand,predictor)
  
    #Since ortho one can reduce the DS model later
    xc=cor(attr(ds_obj,'evaluation'),attr(ds_obj,'evaluation'))
    crossvals<-diag(xc[seq(2,nrow(xc),2),seq(1,ncol(xc),2)])
    training[n_predictor_EOFs-1,1:(ncol(training)-1)]<-round(crossvals,3)
    training[n_predictor_EOFs-1,ncol(training)]<-round(sum(eofvar(predictand)[crossvals>0.2]/sum(eofvar(predictand))*crossvals[crossvals>0.2]),3)
  }
  training$y.name=nam.y
  training$x.name=nam.x
  training$y.res=nam.y.res
  training$t.name=nam.t
  training$dom.name=nam.x.dom
  tm<-cbind(setNames(data.frame(rownames(training)),c('names')), data.frame(training, row.names=NULL))
  tmn<-tm[,c(ncol(tm)-3, ncol(tm)-2, ncol(tm)-1,  ncol(tm),1:(ncol(tm)-4))]
  if (verbose) {cat('\n .')  
    print(tmn) }
  invisible(tmn)
}

## compute the percentage of missing data in x
#' @export
missval <- function(x) sum(is.na(coredata(x)))/length(coredata(x))

## compute the quantile 95% of x
#' @export
q95 <- function(x,na.rm=TRUE) quantile(x,probs=.95,na.rm=na.rm)

## compute the quantile 5% of x
#' @export
q5 <- function(x,na.rm=TRUE) quantile(x,probs=.05,na.rm=na.rm)

## compute the quantile 5% of x
#' @export
q995 <- function(x,na.rm=TRUE) quantile(x,probs=.995,na.rm=na.rm)

## compute the quantile 5% of x
#' @export
q975 <- function(x,na.rm=TRUE) quantile(x,probs=.975,na.rm=na.rm)

## count the number of valid data points
#' @export
nv <- function(x,...) sum(is.finite(x))

## Compute the coefficient of variation of x
#' @export
cv <- function(x,na.rm=TRUE) {sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm)}

#' @export
stand <- function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

#' @export
zeros <- function(x) (sum(is.infinite(1/x)) > 0)

## Estimate the root-mean-squared-error
#' @export
rmse <- function(x,y,na.rm=TRUE) {
  z <- sqrt( (x - y)^2 )
  z <- sum(z,na.rm=na.rm)/sum(is.finite(z),na.rm=na.rm)
  return(z)             
}

#' @export
RMSE <- function(x,y,...) return(rmse(x,y,...))

#' Wrap-around for lag.zoo to work on station and field objects
#'
#' @aliases lag.field
#'
#' @importFrom stats lag
#'
#' @export lag.station
lag.station <- function(x,...) {
  y <- lag(zoo(x),...)
  y <- attrcp(x,y)
  class(y) <- class(x)
  invisible(y)
}

#' @export lag.field
lag.field <- function(x,...) lag.station(x,...)

#' @export
exit <- function() q(save="no")

#' @export
filt <- function(x,n,type='ma',lowpass=TRUE) UseMethod("filt")

#' @exportS3Method
#' @export filt.default
filt.default <- function(x,n,type='ma',lowpass=TRUE) {
  
  # A number of different filters using different window
  # shapes.
  #
  # R.E. Benestad, July, 2002, met.no.
  #
  # ref: Press et al. (1989), Numerical Recipes in Pascal, pp. 466
  #library(ts)

  # Moving-average (box-car) filter
  ma.filt <- function(x,n) {
    y <- filter(x,rep(1,n)/n)
    y
  }

  # Gaussian filter with cut-off at 0.025 and 0.975 of the area.
  gauss.filt <- function(x,n) {
    i <- seq(0,qnorm(0.975),length=n/2)
    win <- dnorm(c(sort(-i),i))
    win <- win/sum(win)
    y <- filter(x,win)
    y
  }

  # Binomial filter
  binom.filt <- function(x,n) {
    win <- choose(n-1,0:(n-1))
    win <- win/max(win,na.rm=T)
    win[is.na(win)] <- 1
    win <- win/sum(win,na.rm=T)
    y <- filter(x,win)
    y
  }

  # Parzen filter (Press,et al. (1989))
  parzen.filt  <-  function(x,n) {
    j <- 0:(n-1)
    win <- 1 - abs((j - 0.5*(n-1))/(0.5*(n+1)))
    win <- win/sum(win)
    y <- filter(x,win)
    y
  }

  # Hanning filter (Press,et al. (1989))
  hanning.filt  <-  function(x,n) {
    j <- 0:(n-1)
    win <- 0.5*(1-cos(2*pi*j/(n-1)))
    win <- win/sum(win)
    y <- filter(x,win)
    y
  }

  # Welch filter (Press,et al. (1989))
  welch.filt  <-  function(x,n) {
    j <- 0:(n-1)
    win <- 1 - ((j - 0.5*(n-1))/(0.5*(n+1)))^2
    win <- win/sum(win)
    y <- filter(x,win)
    y
  }

  y <- coredata(x)
  z <- eval(parse(text=paste(type,'filt(y,n)',sep='.')))
  hp <- as.numeric(y - coredata(z))
  if (!is.null(dim(x))) dim(hp) <- dim(x)
  if (lowpass) coredata(x) <- coredata(z) else
               coredata(x) <- hp
  attr(x,'history') <- history.stamp(x)
  return(x)
}

#' @export
figlab <- function(x,xpos=0.001,ypos=0.001) {
  par(new=TRUE,fig=c(0,1,0,1),xaxt='n',yaxt='n',bty='n',mar=rep(0,4))
  plot(c(0,1),c(0,1),type='n')
  text(xpos,ypos,x,cex=1.2,pos=4,col='grey30')
}

#' @export
ensemblemean <- function(x,FUN='rowMeans') {
  if (inherits(x,'pca')) z <- as.station(x) else z <- x
  ## Estimate the ensemble mean
  zz <- unlist(lapply(coredata(z),FUN=FUN))
  zm <- matrix(zz,length(index(z[[1]])),length(z))
  zm <- zoo(zm,order.by=index(z[[1]]))
  zm <- as.station(zm,param=varid(z),unit=unit(z),
                   loc=unlist(lapply(z,loc)),lon=unlist(lapply(z,lon)),
                   lat=unlist(lapply(z,lat)),alt=unlist(lapply(z,alt)),
                   longname=attr(x,'longname'),aspect=attr(x,'aspect'),
                   info='Ensemble mean ESD')
  invisible(zm)
}

#' @export
propchange <- function(x,it0=c(1979,2013)) {
  z <- coredata(x)
  if (is.null(dim(z)))
      z <- 100*(z/mean(coredata(subset(x,it=it0)),na.rm=TRUE)) else
      z <- 100*t(t(z)/apply(coredata(subset(x,it=it0)),2,'mean',na.rm=TRUE))
  attributes(z) <- NULL
  z -> coredata(x)  
  x
}

#' @export
arec <- function(x,...) UseMethod("arec")

#' @exportS3Method
#' @export arec.default
arec.default <- function(x,...) {
  y <- length(records(x))/sum(1/(1:nv(x)))
  return(y)
}

#' @exportS3Method
arec.station <- function(x,...) {
  y <- unlist(lapply(records(x),length))/apply(x,2,function(x) sum(1/(1:nv(x))))
  return(y)
}

#' @export
lastrains <- function(x,x0=1,uptodate=TRUE,verbose=FALSE) {
  if (verbose) print('lastrains')
  ## Clean up missing values
  x <- x[is.finite(x)]
  y <- cumsum(rev(coredata(x)))
  z <- sum(y < x0,na.rm=TRUE)
  if (uptodate) if (Sys.Date() - end(x) > 1) z <- NA 
  return(z)
}

#' @export
lastdry <- function(x,x0=1,uptodate=TRUE,verbose=FALSE) {
  if (verbose) print('lastdry')
  ## Clean up missing values
  x <- x[is.finite(x)]
  y <- rev(coredata(x))
  z <- (1:length(y))[y < x0][1] - 1
  if (uptodate) if (Sys.Date() - end(x) > 1) z <- NA 
  return(z)
}

#' @export
lastelementrecord <- function(x,verbose=FALSE) {
  ## Checks last element of the record to see if they are the highest - a record
  if (verbose) print('lastelementrecord')
  ## If minimum, then multiply x with -1
  y <- coredata(x)
  nt <- length(index(x))
  if (length(dim(y)) == 2) {
    z <- rep(0,dim(y)[2])
    validlast <- is.finite(y[nt,])
    if (verbose) print(sum(validlast))
    if (sum(validlast)>1)
      z[validlast] <- apply(y[,validlast],2,function(x) if (x[length(x)] == max(x,na.rm=TRUE)) 1 else 0) else
    if (sum(validlast)>1) if (y[nt,validlast]==max(y[,validlast],na.rm=TRUE)) z <- 1    
  } else {
    z <- 0
    if(is.finite(y[nt])) 
      if (y[nt]==max(y,na.rm=TRUE)) z <- 1
  }
  return(z)
}

#' @export
bin <- function(x,nbins=5,labels=NULL,na.rm=TRUE) {
  if (na.rm) good <- is.finite(x) else good <- rep(TRUE,length(x))
  rank <- order(x[good])
  y <- (rank -1) %/% trunc(length(rank)/nbins) + 1
  z <- rep(NA,length(x))
  z[good] <- y
  if (is.null(labels)) names(z) <- quantile(x,probs = (1:nbins -0.5)/nbins)
  return(z)
}

#' @export
factor2numeric <- function(f) {
  if(!is.null(levels(f))) {return(as.numeric(levels(f))[f])
  } else return(as.numeric(f))
}

