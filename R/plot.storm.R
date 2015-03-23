## Author 	 Kajsa Parding
## Last update   23.03.2015

plot.storm <- function(x,it=NULL,is=NULL,
      col='black',type='b',lwd=1,lty=1,pch=19,
      main=NULL,xlim=NULL,ylim=NULL,new=TRUE) {

  y <- subset.storm(x,it=it,is=is)
  n <- count.storm(y,by='month')
  i <- c(min(which(min(month(n))==month(n))),
         max(which(max(month(n))==month(n))))
  ny <- aggregate(n[i[1]:i[2]],by=year,FUN=sum)
    
  if (new | dev.cur()==1) {
    dev.new()
    plot(ny,type='n',main=main,xlim=xlim,ylim=ylim,
         xlab='Time',ylab='Storms count  (events/year)')
  }
  lines(ny,type=type,col=col,lwd=lwd,lty=lty,pch=pch)
  
}


