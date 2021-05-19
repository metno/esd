#' Diagnosis plots
#'
#' Produce plots to diagnose and examine different data types and checks for consistency. The method
#' plot.diagnose.matrix has an additional agrument 'calibrate' that can make the target plot presentation
#' stricter or leaner (default is 100).  
#'
#' @aliases plot.diagnose.comb.eof plot.diagnose.dsensemble plot.diagnose.matrix
#' @seealso diagnose
#'
#' @exportS3Method
#' @export plot.diagnose
plot.diagnose <- function(x,...) {
  if ( (inherits(x,"eof")) & (inherits(x,"comb")) ) plot.diagnose.comb.eof(x,...) else
    if (inherits(x,"dsensembles")) plot.diagnose.dsensemble(x,...)
}

#' @export plot.diagnose.comb.eof
plot.diagnose.comb.eof <- function(x,xlim=NULL,ylim=NULL,add=FALSE,new=TRUE,
                                   alpha=0.5,lwd=2,verbose=FALSE,...) {
  if (verbose) print('plot.diagnose.comb.eof')
  stopifnot(!missing(x), inherits(x,"diagnose"),
            inherits(x,"eof"),inherits(x,"comb"))
  
  n <- length(x$mean.diff)
  j <- 1:n
  col <- rgb(j/n,abs(sin(pi*j/n)),(1-j/n),alpha)
  
  ## KMP 2016-11-02 xlim changed because of changes in diagnose.comb.eof
  #if (is.null(xlim)) xlim <- range(abs(c(0,1,x$mean.diff)),na.rm=TRUE)
  if (is.null(xlim)) xlim <- range(c(-1,1,x$mean.diff),na.rm=TRUE)
  if (is.null(ylim)) ylim <- range(c(-1,1,x$sd.ratio),na.rm=TRUE)
  
  if (!add) {
    if (new) dev.new()
    par(bty="n")
    par0 <- par()
    wt <- 0:360
    plot(cos(pi*wt/180),sin(pi*wt/180),type="l",
         ## KMP 2016-11-02 changed definition of the sd ratio and diff (see diagnose.comb.eof)
         #xlab="mean difference",
         #ylab=expression(paste("|",sigma[pre] - sigma[ref],"|/",sigma[pre])),
         #xlab=expression(paste("(mean difference)/",sigma[ref])),
         xlab=expression(paste("(",mean[ref] - mean[pre],")/",sigma[ref])),
         ylab=expression(paste("(",sigma[ref] - sigma[pre],")/",sigma[ref])),
         main=paste("Diagnostics: common EOFs",attr(x,'variable')),
         xlim=xlim,ylim=ylim,col="grey",
         sub=paste(x$calibrationdata," - ",rownames(x$mean.diff),collapse = "/"))
    lines(c(-10,10),rep(0,2))
    lines(rep(0,2),c(-10,10))
    grid()
    xpos <- xlim[2] - 0.2*diff(xlim)
    ypos <- ylim[1] + 0.1*diff(ylim)
    #    legend(xlim[1],ylim[2],c("same sign","different sign"),
    #           pch=c(19,21),bty="n",col="grey")
    legend(xpos,ypos,c("same sign","different sign"),
           pch=c(19,21),bty="n",col="grey")
    par(xpd=TRUE)
    #    text(xlim[1],ylim[2],'AR(1) - symbol size',col='grey40',pos=3)
    text(xlim[2]-0.05*diff(xlim),ypos+0.025*diff(ylim),
         'AR(1) - symbol size',col='grey40',pos=3,cex=0.8)
    text(xlim[2],ylim[2],'EOF #',col='grey40',cex=0.8,pos=3)
    
    par(new=TRUE,fig=c(0.85,0.95,0.70,0.85),mar=c(0,3,0,0),
        cex.axis=0.7,yaxt="s",xaxt="n",las=1)
    colbar <- rbind(1:n,1:n)
    image(1:2,1:n,colbar,col=col)
    par(fig=par0$fig,mar=par0$mar,cex.axis=par0$cex.axis,
        yaxt=par0$yaxt,xaxt=par0$xaxt,las=par0$las,new=TRUE)
    plot(cos(pi*wt/180),sin(pi*wt/180),type="n",xlim=xlim,ylim=ylim,
         xlab='',ylab='',main='',sub='')
    par(par0$new)
  }
  cex <- x$autocorr.ratio*1.5;
  pch <- rep(19,n); pch[cex < 0] <- 21
  cex <- abs(cex); cex[cex > 2.5] <- 2.5
  if (verbose) {
    print('Mean difference:');print(x$mean.diff)
    print('Ration of standard deviation');print(x$sd.ratio)
    print('Size');print(cex)
    print('col');print(col)
  }
  points(x$mean.diff,x$sd.ratio,pch=pch,col=col,lwd=lwd,cex=cex)
}

#' @export plot.diagnose.matrix
plot.diagnose.matrix <- function(x,xlim=NULL,ylim=NULL,verbose=FALSE,new=TRUE,...) {
  if (verbose) print('plot.diagnose.matrix')
  x <- as.data.frame(x)
  par(bty="n")
  if (is.null(xlim)) xlim <- range(abs(c(0,1,x$mean.diff)),na.rm=TRUE)
  if (is.null(ylim)) ylim <- range(c(-1,1,1-x$sd.ratio),na.rm=TRUE)
  wt <- 0:360
  if (new) dev.new()
  plot(cos(pi*wt/180),sin(pi*wt/180),type="l",
       xlab="mean difference",ylab=expression(1- sigma[p*r*e]/sigma[r*e*f]),
       main=paste("Diagnostics: common EOFs",attr(x,'variable')),
       xlim=xlim,ylim=ylim,col="grey",
       sub=paste(x$calibrationdata," - ",rownames(x$mean.diff),collapse = "/"))
  lines(c(0,10),rep(0,2))
  lines(rep(0,2),c(0,10))
  n <- length(x$mean.diff)
  j <- 1:n
  col <- rgb(j/n,abs(sin(pi*j/n)),(1-j/n))
  cex <- x$autocorr.ratio;
  cex[!is.finite(cex)] <- 1
  pch <- rep(19,n); pch[cex < 0] <- 21
  cex <- abs(cex); cex[cex > 2] <- 2
  if (verbose) {
    print('Mean difference:');print(x$mean.diff)
    print('Ration of standard deviation');print(x$sd.ratio)
    print('Size');print(cex)
    print('col');print(col)
    points(x$mean.diff,1-x$sd.ratio,pch=pch,col='grey75',cex=1)
  }
  
  points(abs(x$mean.diff),1-x$sd.ratio,pch=pch,col=col,cex=cex)
  legend(xlim[1],ylim[2],c("same sign","different sign"),
         pch=c(19,21),bty="n",col="grey")
  par(xpd=TRUE)
  text(xlim[1],ylim[2],'AR(1) - symbol size',col='grey40',pos=3)
  
  text(xlim[2],ylim[2],'EOF #',col='grey40',cex=0.8,pos=3)
  par(new=TRUE,fig=c(0.85,0.95,0.70,0.85),mar=c(0,3,0,0),
      cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:n,1:n)
  image(1:2,1:n,colbar,col=col)  
}

#' @exportS3Method
#' @export plot.diagnose.dsensemble
plot.diagnose.dsensemble <- function(x,new=TRUE,mgp=c(2,1,0),cex=NULL,map.show=TRUE,
                                     map.type=NULL,verbose=FALSE,main=NULL,...) {
  if (verbose) {
    print('plot.diagnose.dsensemble')
    for (i in 1:length(x$outside)) { 
      print(c(x$outside[i],qbinom(c(0.05,0.95),size=x$N,prob=0.1)))
      print(x$deltaobs[i])
      print(summary(x$deltagcm[i,]))
    }
  }
  ## REB 2021-05-10. The …operator, technically known as the ellipsis, allows a function to take arguments that are not predefined in 
  ## its definition. The ellipsis is useful when we don’t know how many arguments a function may take. 
  ## Use this to calibrate the target plot. The argument 'calibrate' 
  l <- list(...)
  ic <- match('calibrate',names(l))
  if (!is.na(ic)) {
    if (verbose) print(paste('calibrate=',l$calibrate))
    r0 <- l$calibrate
  } else r0 <- 200
  if (verbose) print(paste('r0=',r0))
  
  if (is.null(map.type)) {
    if( inherits(x,"field") | length(lon(x))!=length(lat(x)) |
        (length(lon(x))==2 & length(lat(x))==2) ) {
      map.type <- "rectangle"
    } else {
      map.type <- "points"
    }
  }
  
  Y <- -c(r0*(0.5-pbinom(x$outside,size=x$N,prob=0.1)))
  if (verbose) print(Y)
  ## KMP 2020-08-04: Changed X (probability of observed trends) because the 
  ## calculations did not seem right. Trend distributions for different stations 
  ## are mixed up when applying mean(deltagcm) and sd(deltagcm). 
  #X <- -round(200*(0.5-pnorm(x$deltaobs,mean=mean(x$deltagcm),
  #                           sd=sd(x$deltagcm))),2)
  ## REB 2020-11-12: the changes caused a problem - probably used in different ways..
  if (!is.null(dim(x$deltagcm))) 
    X <- -r0*sapply(seq_along(x$deltaobs), function(i) 0.5 - pnorm(x$deltaobs[i], 
                                                                      mean=mean(x$deltagcm[i,]), sd=sd(x$deltagcm[i,]))) else
    X <- -r0*(0.5-pnorm(x$deltaobs,mean=mean(x$deltagcm),  sd=sd(x$deltagcm)))                               
  if (verbose) print(X)            
  if (new) {
    dev.new()
    par(bty="n",fig=c(0.05,0.95,0,0.95),mgp=mgp)
  }
  if (!is.null(cex)) par(cex=cex)
  plot(c(-100,100),c(-100,100),type="n",axes=FALSE,mgp=mgp,
       ylab="",xlab="",main=main)
  par(las=0)
  mtext("trend",side=1,line=1.5,cex=par("cex"))
  mtext("standard deviation",side=2,line=2,srt=180,cex=par("cex"))
  u <- par("usr")
  dx <- (u[2]-u[1])/20
  dy <- (u[4]-u[3])/20
  arrows(u[1]+dx,u[3],u[2]-dx,u[3],
         lwd=0.75,length=0.1,angle=20,code=2,xpd=NA)
  arrows(u[2]-dx,u[3],u[1]+dx,u[3],
         lwd=0.75,length=0.1,angle=20,code=2,xpd=NA)
  arrows(u[1],u[4]-dy,u[1],u[3]+dy,
         lwd=0.75,length=0.1,angle=20,code=2,xpd=NA)
  arrows(u[1],u[3]+dy,u[1],u[4]-dy,
         lwd=0.75,length=0.1,angle=20,code=2,xpd=NA)
  mtext("ensemble > obs",side=1,line=0,adj=0,cex=par("cex")*0.65)
  mtext("ensemble < obs",side=1,line=0,adj=1,cex=par("cex")*0.65)
  mtext("ensemble > obs",side=2,line=0.5,adj=0,cex=par("cex")*0.65)
  mtext("ensemble < obs",side=2,line=0.5,adj=1,cex=par("cex")*0.65)  
  bcol=c("grey95","grey40")
  for (i in 1:10) {
    r <- (11-i)*10
    polygon(r*cos(pi*seq(0,2,length=360)),
            r*sin(pi*seq(0,2,length=360)),
            col=bcol[i %% 2 + 1],border="grey15")
  }
  for (i in seq(0,90,by=1))
    points(X,Y,pch=19,cex=2 - i/50,col=rgb(i/90,0,0))
  if (map.show) {
    if(verbose) print("add map") 
    if(is.null(x$xrange) & !is.null(attr(x$z,"lon"))) {
      x$xrange <- c(attr(x$z,"lon")-15,attr(x$z,"lon")+15)
    }
    if(is.null(x$yrange) & !is.null(attr(x$z,"lat"))) {
      x$yrange <- c(attr(x$z,"lat")-10,attr(x$z,"lat")+10)
    }
    if (!is.null(x$xrange) & !is.null(x$xrange)) {
      data("geoborders", envir = environment())
      lon <- geoborders$x
      lat <- geoborders$y
      ok <- lon>(min(x$xrange)-1) & lon<(max(x$xrange)+1) &
        lat>(min(x$yrange)-1) & lat<(max(x$yrange)+1)
      lon2 <- attr(geoborders,"borders")$x
      lat2 <- attr(geoborders,"borders")$y
      ok2 <- lon2>(min(x$xrange)-1) & lon2<(max(x$xrange)+1) &
        lat2>(min(x$yrange)-1) & lat2<(max(x$yrange)+1)
      par(fig=c(0.78,0.97,0.78,0.97),new=TRUE, mar=c(0,0,0,0),
          cex.main=0.75,xpd=NA,col.main="grey",bty="n")
      plot(lon[ok],lat[ok],lwd=1,col="black",type='l',xlab=NA,ylab=NA,
           axes=FALSE)
      axis(1,mgp=c(3,.5,0))
      axis(2,mgp=c(2,.5,0))
      lines(lon2[ok2],lat2[ok2],col = "pink",lwd=1)
      if("points" %in% map.type) {
        points(attr(x$z,"lon"),attr(x$z,"lat"),pch=21,cex=1,col='black',bg='red',lwd=1)
      }
      if("rectangle" %in% map.type) {
        rect(min(attr(x$z,"lon")),min(attr(x$z,"lat")),max(attr(x$z,"lon")),
             max(attr(x$z,"lat")),lwd=1,col=NA,border='red',lty=2)
      }
    } 
  }  
}