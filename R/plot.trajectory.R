## Author 	 Kajsa Parding
## Rasmus: small fix
## Last update   24.04.2015, 05.04.2017

plot.trajectory <- function(x,it=NULL,is=NULL,
      main=NULL,xlim=NULL,ylim=NULL,
      col=NULL,pch=0,type='l',lwd=3,
      xlab="",ylab=NULL,new=TRUE,verbose=FALSE,...) {
  y <- subset(x,it=it,is=is)
  n <- count.trajectory(y,by='year')
  if(new) dev.new()
  plot.station(n,main=main,new=new,col=col,
               xlim=xlim,ylim=ylim,type=type,
               lwd=lwd,pch=pch,xlab=xlab,ylab=ylab,
               legend.show=FALSE,...)
  if (is.null(col)) col <- rainbow(length(y[1,]))
  if (!is.null(attr(y,"longitude")) |
      !is.null(attr(y,"latitude"))) {
    if(verbose) print("adding legend")
    par(xpd=TRUE)
    leg <- ""
    if (!any(is.na(attr(y,"longitude"))) & !is.null(attr(y,"longitude"))) {
      leg <- paste(leg,paste(round(range(attr(y,"longitude")),2),sep="-"),"E/",sep="")
    }
    if (!any(is.na(attr(y,"latitude"))) & !is.null(attr(y,"latitude"))) {
      leg <- paste(leg,paste(round(range(attr(y,"latitude")),2),sep="-"),"N",sep="")
    }
    if (is.null(col)) col <- rainbow(1)
    if(verbose) print(paste('legend:',leg,', color:',col[1]))
    legend("bottomleft",inset=c(0,-0.25),legend=leg,bty="n",cex=0.6,ncol=3,
         text.col="grey40",lty=1,col=col)
  }
  invisible(n)
}
