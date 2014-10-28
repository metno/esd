# R.E. Benestad, met.no, Oslo, Norway 12.04.2013
# rasmus.benestad@met.no
#------------------------------------------------------------------------


trend<-function(x,result="trend",model="y ~ t",...) UseMethod("trend")

trend.default <- function(x,result="trend",model="y ~ t",...) {
  trendx <- data.frame(t=1:length(x),y=x)
  eval(parse(text=paste("xt <- lm(",model,",data=trendx)")))
  y <- switch(result,"trend"=zoo(predict(xt,newdata=trendx),order.by=index(x)),
                     "residual"=zoo(xt$residuals,order.by=index(x)))
  attr(y,'history') <- history.stamp(x)
  return(y)
}

trend.one.station <- function(x,result="trend",model="y ~ t",...) {
  #print("trend.station")
  #print(class(index(x)))
  if (class(index(x))=="Date") {
    #print("HERE")
    year <- as.numeric( format(index(x), '%Y') ) 
    month <- as.numeric( format(index(x), '%m') )
    day <- as.numeric( format(index(x), '%d') )
    trendx <- data.frame(t=year + (month - 0.5)/12 + (day - 1)/365.25,
                         y=coredata(x))
  } else if ( (class(index(x))=="numeric") |
              (class(index(x))=="integer") )
                trendx <- data.frame(t=index(x),y=coredata(x))
  #print(summary(trendx))
  eval(parse(text=paste("xt <- lm(",model,",data=trendx)")))
  #print(summary(summary(xt)))
  #print(result)
  y <- switch(result,"trend"=zoo(predict(xt,newdata=trendx),order.by=index(x)),
                     "residual"=zoo(xt$residuals,order.by=index(x)))
  #str(y)
  #mostattributes(y) <- attributes(x)
#  nattr <- softattr(x)
#  print(nattr)
#  for (i in 1:length(nattr))
#    attr(y,nattr[i]) <- attr(x,nattr[i])
  attr(y,'coefficients') <- xt$coefficients
  attr(y,'original data') <-  x
  attr(y,'aspect') <- result
  #attr(y,'call') <- match.call()
  y <- attrcp(x,y,ignore='aspect')
  attr(y,'history') <- history.stamp(x)
  #print(".")
  return(y)
}

trend.station <- function(x,result="trend",model="y ~ t",...) {
  # Allow for a set of stations.
  d <- dim(x)
  if (is.null(d)) y <- trend.one.station(x,result=result,model=model) else {
      Y <- x*NA
      for (i in 1:d[2]) {
        y <- trend.one.station(x[,i],result=result,model=model)
        Y[,i] <- y
      }
      y <- Y
      #nattr <- softattr(x)
      #print("+")
      #for (i in 1:length(nattr))
      #  attr(y,nattr[i]) <- attr(x,nattr[i])
      y <- attrcp(x,y)
    }
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  return(y)
}

trend.eof <- function(x,result="trend",model="y ~ t",...) {
  class(x) -> cls
  if (class(index(x))=="Date") {
    year <- as.numeric( format(index(x), '%Y') ) 
    month <- as.numeric( format(index(x), '%m') )
    day <- as.numeric( format(index(x), '%d') )
    t <- year + (month - 0.5)/12 + (day - 1)/365.25
  } else if ( (class(index(x))=="numeric") |
              (class(index(x))=="integer") )
    t=index(x)
  #print("detrend.eof")
  d <- dim(x)
  Y <- x
  #print(dim(y))
  nc <- sum(is.element(strsplit(model,"")[[1]],"t")) + 1
  coefficients <- matrix(rep(NA,nc*d[2]),nc,d[2])
  for (i in 1:d[2]) {
    trendx <- data.frame(t=t,y=coredata(x[,i]))
    eval(parse(text=paste("xt <- lm(",model,",data=trendx)")))
  #print(summary(xt))
    Y[,i] <- switch(result,"trend"=zoo(predict(xt,newdata=trendx),order.by=index(x)),
                           "residual"=zoo(xt$residuals,order.by=index(x)))
    coefficients[,i] <- xt$coefficients
    #xt <- lm(y ~ t,data=trendx)
    #z <- predict(xt); print(length(z))
    #Y[,i] <- xt$residual
    #print(summary(Y[,i]))
  }
  #X <- zoo(Y,order.by=index(x))
  #nattr <- softattr(x)
  #for (i in 1:length(nattr))
  #  attr(Y,nattr[i]) <- attr(x,nattr[i])
  #mostattributes(Y) <- attributes(x)
  attr(Y,'coefficients') <- coefficients
  attr(Y,'original data') <-  x
  #attr(Y,'call') <- match.call()
  Y <- attrcp(x,Y)
  attr(Y,'history') <- history.stamp(x)
  
  #plot(x[,1]); lines(X[,1],col="red",lty=2)
  return(Y)
}

trend.field <- function(x,result="trend",model="y ~ t",...) {

  gettrend <- function(x,model="y ~ t") {
    #browser()
    #print(match.call()); str(x); print(model)
    #print(table(year))
    trendx <- data.frame(y=x,t=1:length(x))
    eval(parse(text=paste("trendfit <- lm(",model,",data=trendx)")))
    trend <- predict(trendfit,newdata=trendx)
    return(trend)
  }
  
  class(x) -> cls
  
  print("detrend.field")
  d <- dim(x)
  Y <- x
  #print(dim(Y))
#  for (i in 1:d[2]) {
#    #browser()
#    
#    
#    #print(summary(xt))
#    #str(xt)

  t <- index(x)
  datetype <- class(t)

  if ( (datetype=="numeric") |
                ( (datetype=="character") & (nchar(t[1]==4)) ) ) {
    index(x) <- as.Date(paste(t,'-01-01',sep=''))
  }
  year <- as.numeric( format(index(x), '%Y') ) 
  month <- as.numeric( format(index(x), '%m') )
  day <- as.numeric( format(index(x), '%d') )
  t <- year + (month - 0.5)/12 + (day - 1)/365.25
  nt <- length(t)
  
  #print("test"); ttt <- gettrend(x[,1]); plot(ttt,main="test: trend")
  xt <- apply(x,2,gettrend,model)
  trend <- zoo(xt,order.by=index(x))
  #print("residual"); print(dim(x)); print(dim(trend))
  Y <- x - trend
  #print("attributes")
  Y <- attrcp(x,Y,ignore="aspect")
  #print("pattern"); print(nt)
  pattern <- 10*(coredata(trend)[nt,] - coredata(trend)[1,])/(t[nt]-t[1])
  #print("dimensions")
  dim(pattern) <- attr(x,'dimension')[1:2]
  attr(Y,'dimension') <- attr(x,'dimension')
  attr(Y,'pattern') <- pattern
  attr(Y,'aspect') <- result
  attr(Y,'history') <- history.stamp(x)
  if (result=="trend") {
    class(Y) <- c("trend",class(x))
    
  } else {
    class(Y) <- class(x)
  }
  return(Y)
}

trend.zoo <- function(x,result="trend",model="y ~ t",...) {
  trendx <- data.frame(t=year(x),y=coredata(x))
  eval(parse(text=paste("xt <- lm(",model,",data=trendx)")))
  y <- switch(result,"trend"=zoo(predict(xt,newdata=trendx),order.by=index(x)),
                     "residual"=zoo(xt$residuals,order.by=index(x)))
  attr(y,'history') <- history.stamp(x)
  if (!is.null(attr(x,'unit'))) unit <- attr(x,'unit') else
                                unit <- 'x'
  attr(y,'unit') <- paste(unit,'/ year')
  return(y)
  invisible(y)
}
