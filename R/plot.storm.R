## Author 	 Kajsa Parding
## Last update   08.04.2015

plot.storm <- function(x,it=NULL,is=NULL,
      main=NULL,xlim=NULL,ylim=NULL,
      col=NULL,pch=0,type='l',lwd=3,
      xlab="",ylab=NULL,new=TRUE,verbose=FALSE) {
  y <- subset.storm(x,it=it,is=is)
  n <- count.storm(y,by='year')
  plot.station(n,main=main,new=new,col=col,
               xlim=xlim,ylim=ylim,type=type,
               lwd=lwd,pch=pch,xlab=xlab,ylab=ylab,
               legend.show=FALSE)
  if (is.null(col)) col <- rainbow(length(y[1,]))
  if (sum(is.na(attr(y,'longitude')))==0 |
      sum(is.na(attr(y,'latitude')))==0) {
    if(verbose) print("done plotting time series. let's try the legend")
    par(xpd=TRUE)
    leg <- ''
    if (!any(is.na(attr(y,'longitude')))) leg <- paste(leg,
              round(attr(y,'longitude'),2)[1],'–',
              round(attr(y,'longitude'),2)[2],"E/",sep="")
    if (!any(is.na(attr(y,'latitude')))) leg <- paste(leg,
           round(attr(y,'latitude'),2)[1],'–',
           round(attr(y,'latitude'),2)[2],"N",sep="")
    if (is.null(col)) col <- rainbow(1)
    if(verbose) print(paste('legend',leg,'color',col[1]))
    legend("bottomleft",inset=c(0,-0.25),legend=leg,bty="n",cex=0.6,ncol=3,
         text.col="grey40",lty=1,col=col)
  }
}
