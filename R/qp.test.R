#' quantile-quantile plot
#' 
#' @param x input object
#' @param p percentiles to show in plot, a numeric vector
#' @param threshold threshold used only if input is precipitation data - only values above threshold are included in analysis 
#' 
#' @export qp.test
qp.test <- function(x,p=c(seq(0.1,0.95,0.05),0.97,0.98,0.99),threshold=1,...) UseMethod("qp.test")

#' @export qp.test.station
qp.test.station <- function(x,p=c(seq(0.1,0.95,0.05),0.97,0.98,0.99),threshold=1,...) {
  if (is.precip(x)) qp.test.precip(x,...) else
  if (is.T(x)) qp.test.t2m(x,...)
}

qp.test.precip <- function(x,p=c(seq(0.1,0.95,0.05),0.97,0.98,0.99),threshold=1,...) {
  # From qqplotter:
  x[x < threshold] <- NA
  if (is.null(dim(x))) {
    qp <- quantile(x,prob=p,na.rm=TRUE)
    qmu <- -log(1-p)*mean(x,na.rm=TRUE)
  } else {
    qp <- apply(coredata(x),2,quantile,prob=p,na.rm=TRUE)
    qmu <- -log(1-p)%o%apply(coredata(x),2,wetmean,na.rm=TRUE)
    #fw <- round(100*apply(coredata(x),2,wetfreq))
  }
  plot(qp,qmu,pch=19,col=rgb(0,0,1,0.2),
       xlab=expression(q[p]),ylab=expression(-log(1-p)*mu))
  lines(range(qp,qmu),range(qp,qmu))
  grid()
}

qp.test.t2m <- function(x,p=c(0.01,0.02,0.03,0.04,seq(0.1,0.95,0.05),
                                0.97,0.98,0.99),...) {
  d <- dim(x); if (is.null(d)) d <- c(length(x),1)
  if (is.null(dim(x))) {
    qp <- quantile(x,prob=p,na.rm=TRUE)
    qmu <-  qnorm(p=p,mean=mean(coredata(x),na.rm=TRUE),sd=sd(coredata(x),na.rm=TRUE))
  } else {
    qp <- apply(coredata(x),2,quantile,prob=p,na.rm=TRUE)
    qmu <- qp + NA
    for (i in 1:length(p))
      qmu[i,] <- qnorm(p=p[i],mean=apply(coredata(x),2,mean,na.rm=TRUE),
                       sd=apply(coredata(x),2,sd,na.rm=TRUE))
  }
  plot(qp,qmu,pch=19,col=rgb(1,0,0,0.2),
       xlab=expression(q[p]),ylab='qnorm(p,mean(x),sd(x))')
  lines(range(qp,qmu),range(qp,qmu))
  grid() 
}
