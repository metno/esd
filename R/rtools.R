## Author A. Mezghani
## Description Contains some rtools ...
## Created 14.11.2014

as.decimal <- function(x=NULL) {
    ## converts from degree min sec format to degrees ...
    ##x is in the form "49 deg 17' 38''"
  if (!is.null(x)) {
        deg <-as.numeric(substr(x,1,2)) 
        min <- as.numeric(substr(x,4,5))
        sec <- as.numeric(substr(x,7,8))     
        x <- deg + min/60 + sec/3600
    }
    return(x)
}



## Function to compare downscaled field with original field: subtracts the original data
## from the downscaled data (REB + HBE)
test.ds.field <- function(x,what='xval',verbose=FALSE) {
  if (verbose) print('fieldtest')
  stopifnot (inherits(x,'eof') & inherits(x,'ds'))
  if (verbose) print(colnames(attr(x,'evaluation')))
  isel <- is.element(colnames(attr(x,'evaluation')),paste('X.PCA',1:dim(x)[2],sep='.'))
  eof1 <- attr(x,'evaluation')[,isel]
  isel <- is.element(colnames(attr(x,'evaluation')),paste('Z.PCA',1:dim(x)[2],sep='.'))
  eof2 <- attr(x,'evaluation')[,isel]
  if (verbose) {str(eof1); print(names(attributes(x)))}
  attr(eof1,'variable') <- varid(x)
  attr(eof1,'unit') <- unit(x)
  attr(eof1,'longname') <- attr(x,'longname')
  attr(eof1,'pattern') <- attr(x,'pattern')
  attr(eof1,'eigenvalues') <- attr(x,'eigenvalues')
  attr(eof1,'dimensions') <- c(dim(attr(x,'evaluate'))[1],dim(attr(x,'pattern'))[3])
  attr(eof1,'mean') <- attr(x,'mean')
  attr(eof1,'longitude') <- lon(x)
  attr(eof1,'latitude') <- lat(x)
  attr(eof1,'max.autocor') <- attr(x,'max.autocor')
  attr(eof1,'eigenvalues') <- attr(x,'eigenvalues')
  attr(eof1,'sum.eigenv') <- attr(x,'sum.eigenv')
  attr(eof1,'tot.var') <- attr(x,'tot.var')
  attr(eof1,'aspect') <- 'anomaly'
  attr(eof1,'dimnames') <- NULL   # REB 2016-03-04
  class(eof1) <- c("eof",class(x))[3:6]
  eof2 <- attrcp(eof1,eof2)
  class(eof2) <- c("eof",class(x))[3:6]
  if (verbose) print(c(dim(eof1),dim(eof2)))
  if (verbose) print('estimated EOFs - now get the fields..')
  x1 <- as.field(eof1)
  x2 <- as.field(eof2)
  if (verbose) print(c(dim(x1),dim(x2)))
  coredata(x1) <- coredata(x1 - x2)
  attr(x,'history') <- history.stamp()
  attr(x,'info') <- 'orginale - downscaled'
  invisible(x1)
}

## Simplify life to extract the variance of the EOFs
eofvar <- function(x) if (inherits(x,c('eof','pca'))) 
                          attr(x,'eigenvalues')^2/attr(x,'tot.var')*100 else NULL

## Iterate using n number of predictands in the downscaling and retrive the cross-val given the number of predictands   
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
missval <- function(x) sum(is.na(coredata(x)))/length(coredata(x))

## compute the quantile 95% of x
q95 <- function(x,na.rm=TRUE) quantile(x,probs=.95,na.rm=na.rm)

## compute the quantile 5% of x
q5 <- function(x,na.rm=TRUE) quantile(x,probs=.05,na.rm=na.rm)

## compute the quantile 5% of x
q995 <- function(x,na.rm=TRUE) quantile(x,probs=.995,na.rm=na.rm)

## compute the quantile 5% of x
q975 <- function(x,na.rm=TRUE) quantile(x,probs=.975,na.rm=na.rm)

## count the number of valid data points
nv <- function(x,...) sum(is.finite(x))

## Compute the coefficient of variation of x
cv <- function(x,na.rm=TRUE) {sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm)}

stand <- function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

## Estimate the root-mean-squared-error
rmse <- function(x,y,na.rm=TRUE) {
  z <- sqrt( (x - y)^2 )
  z <- sum(z,na.rm=na.rm)/sum(is.finite(z),na.rm=na.rm)
  return(z)             
}
RMSE <- function(x,y,...) return(rmse(x,y,...))

# Wrap-around for lag.zoo to work on station and field objects:
lag.station <- function(x,...) {
  y <- lag(zoo(x),...)
  y <- attrcp(x,y)
  class(y) <- class(x)
  invisible(y)
}

lag.field <- function(x,...) lag.station(x,...)
  
exit <- function() q(save="no")

filt <- function(x,...) UseMethod("filt")

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
  
figlab <- function(x,xpos=0.001,ypos=0.001) {
  par(new=TRUE,fig=c(0,1,0,1),xaxt='n',yaxt='n',bty='n',mar=rep(0,4))
  plot(c(0,1),c(0,1),type='n')
  text(xpos,ypos,x,cex=1.2,pos=4,col='grey30')
}

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


propchange <- function(x,it0=c(1979,2013)) {
  z <- coredata(x)
  if (is.null(dim(z)))
      z <- 100*(z/mean(coredata(subset(x,it=it0)),na.rm=TRUE)) else
      z <- 100*t(t(z)/apply(coredata(subset(x,it=it0)),2,'mean',na.rm=TRUE))
  attributes(z) <- NULL
  z -> coredata(x)  
  x
}

arec <- function(x,...) UseMethod("arec")

arec.default <- function(x,...) {
  y <- length(records(x))/sum(1/(1:nv(x)))
  return(y)
}

arec.station <- function(x,...) {
  y <- unlist(lapply(records(x),length))/apply(x,2,function(x) sum(1/(1:nv(x))))
  return(y)
}

lastrains <- function(x,x0=1,uptodate=TRUE,verbose=FALSE) {
  if (verbose) print('lastrains')
  ## Clean up missing values
  x <- x[is.finite(x)]
  y <- cumsum(rev(coredata(x)))
  z <- sum(y < x0,na.rm=TRUE)
  if (uptodate) if (Sys.Date() - end(x) > 1) z <- NA 
  return(z)
}

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

bin <- function(x,nbins=5,labels=NULL,na.rm=TRUE) {
  if (na.rm) good <- is.finite(x) else good <- rep(TRUE,length(x))
  rank <- order(x[good])
  y <- (rank -1) %/% trunc(length(rank)/nbins) + 1
  z <- rep(NA,length(x))
  z[good] <- y
  if (is.null(labels)) names(z) <- quantile(x,probs = (1:nbins -0.5)/nbins)
  return(z)
}
