## Author 	 Kajsa Parding
## Last update   24.04.2015

plot.trajectory <- function(x,it=NULL,is=NULL,
      main=NULL,xlim=NULL,ylim=NULL,
      col=NULL,pch=0,type='l',lwd=3,
      xlab="",ylab=NULL,new=TRUE,verbose=FALSE) {
  y <- subset(x,it=it,is=is)
  n <- count.trajectory(y,by='year')
  if(new) dev.new()
  plot.station(n,main=main,new=new,col=col,
               xlim=xlim,ylim=ylim,type=type,
               lwd=lwd,pch=pch,xlab=xlab,ylab=ylab,
               legend.show=FALSE)
  if (is.null(col)) col <- rainbow(length(y[1,]))
  if (sum(is.na(lon(y)))==0 |
      sum(is.na(lat(y)))==0) {
    if(verbose) print("adding legend")
    par(xpd=TRUE)
    leg <- ''
    if (!any(is.na(lon(y))) & !is.null(lon(y))) 
      leg <- paste(leg,round(lon(y),2)[1],'–',round(lon(y),2)[2],"E/",sep="")
    if (!any(is.na(lat(y))) & !is.null(lat(y))) 
      leg <- paste(leg,round(lat(y),2)[1],'–',round(lat(y),2)[2],"N",sep="")
    if (is.null(col)) col <- rainbow(1)
    if(verbose) print(paste('legend:',leg,', color:',col[1]))
    legend("bottomleft",inset=c(0,-0.25),legend=leg,bty="n",cex=0.6,ncol=3,
         text.col="grey40",lty=1,col=col)
  }
  invisible(n)
}
