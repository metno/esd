plot.ds.pca <- function(x,ip=1,
                        colbar1=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                     type="p",cex=1,show=TRUE,h=0.6, v=1,pos=0.05),
                        colbar2=NULL,mar=c(3,2,2,0.5),mgp=c(1,0.5,0),
                        verbose=FALSE,...) {
  y <- x # quick fix
  
  if (verbose) print('plot.ds.pca')
  if (is.null(colbar2)) colbar2 <- colbar1
  attr(y,'longname') <- attr(y,'longname')[1]
  #par(fig=c(0,0.45,0.5,0.975),new=TRUE)
  par(fig=c(0,0.5,0.5,0.975), mar=mar, mgp=mgp) #par(fig=c(0,0.45,0.5,0.975))
  
  if (verbose) print('PCA ip')
  map.pca(y,ip=ip,new=FALSE,colbar=colbar1,
          fig=c(0,0.5,0.5,0.975),
          main=paste("PCA Pattern #",ip," (unitless)",sep=""),
          verbose=verbose,...)
  
  #title(paste("PCA Pattern # ",ip,sep=""))
  par(fig=c(0.55,0.975,0.5,0.975),new=TRUE)
  if (verbose) print('Predictor pattern')
  pp <- attr(y,'predictor.pattern')
  attr(pp,"variable") <- paste("EOF.Pattern.",ip,sep="")
  attr(pp,"unit") <- "unitless"
  map(pp,ip=ip,new=FALSE,
      colbar=colbar2,verbose=verbose,
      main="")#paste("EOF Pattern # ",ip,sep=""))
  title(paste("EOF Pattern #",ip,sep=""))
  if (!is.null(attr(y,'evaluation'))) {
    if (verbose) print('Evaluation results')
    par(fig=c(0.05,0.45,0.05,0.475),new=TRUE)
    ## Get the right pattern
    xvp <- (ip-1)*2 +1
    xok <- is.finite(attr(y,'evaluation')[,xvp]) & is.finite(attr(y,'evaluation')[,xvp+1])
    xcor <- cor(attr(y,'evaluation')[xok,xvp],attr(y,'evaluation')[xok,xvp+1])
    
    par(mgp=c(2,0.5,0), mar=c(3,3.5,1.5,1.5))
    plot(attr(y,'evaluation')[,xvp],attr(y,'evaluation')[,xvp+1],
         main=paste('Cross-validation: r=',round(xcor,2)),
         xlab='original data',ylab='prediction',pch=19,col="grey")
    lines(range(c(attr(y,'evaluation')),na.rm=TRUE),
          range(c(attr(y,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(attr(y,'evaluation')[,1]),
                      x=coredata(attr(y,'evaluation')[,2]))
    xvalfit <- lm(y ~ x, data = cal)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    #browser()
    #legend("bottomleft", )
    par(fig=c(0.55,0.975,0.05,0.475),new=TRUE)
    xlim <- range(index(attr(y,'original_data')),index(y))
    ylim <- range(attr(y,'original_data')[,ip],y[,ip],na.rm=TRUE) + 
      diff(range(attr(y,'original_data')[,ip],y[,ip],na.rm=TRUE))*c(-0.1,0.2)
    y0 <- attr(y,'original_data')
    plot(y0[,ip], lwd=2, type='b', pch=19, xlim=xlim, ylim=ylim, 
         xlab="Date",ylab=paste("PC #",ip,sep=""))
    if ( (class(index(y))=='Date') & (class(index(y0))=='numeric') & inherits(y,'annual') ) 
      index(y) <- year(index(y))
    if ( (class(index(y))=='numeric') & (class(index(y0))=='Date') & inherits(y,'annual') ) 
      index(y) <- as.Date(paste(index(y),'01-01',sep='-'))
    
    lines(zoo(y[,ip]),lwd=2,col='red',type='b')
    legend(x=index(attr(y,'original_data')[,ip])[1],
           y=max(attr(y,'original_data')[,ip],na.rm=TRUE)+
             0.2*diff(range(attr(y,'original_data')[,ip])),
           legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
           lwd=c(2,2),pch=c(21,19),bty="n")
  } else {
    par(fig=c(0.05,0.975,0.05,0.475),new=TRUE)
    plot(attr(y,'original_data')[,ip],lwd=2,type='b',pch=19)
    lines(zoo(y[,ip]),lwd=2,col='red',type='b')
    xvalfit <- NULL
  }
}
