plot.station <- function(x,plot.type="single",new=TRUE,
                         lwd=3,type='l',pch=0,main=NULL,col=NULL,
                         xlim=NULL,ylim=NULL,xlab="",ylab=NULL,
                         errorbar=TRUE,legend.show=FALSE,
                         map.show=TRUE,map.type=NULL,map.insert=TRUE,
                         usegooglemap=TRUE,zoom=NULL,
                         cex.axis=1.2,cex.lab=1.2,cex.main=1.2,
                         mar=c(4.5,4.5,0.75,0.5),fig=NULL,
                         alpha=0.5,alpha.map=0.7,
                         verbose=FALSE,...) {
  
  if (verbose) print('plot.station')
  par(las=1)
  if (!is.numeric(lon(x)) | !is.numeric(lat(x))) {
    map.show <- FALSE
  }
  if(map.show) {
    if (verbose) print('show map')
    if (is.null(map.type)) {
      if( inherits(x,"field") | length(lon(x))!=length(lat(x)) |
          (length(lon(x))==2 & length(lat(x))==2) ) {
        map.type <- "rectangle"
      } else {
        map.type <- "points"
      }
    }
    if (verbose) print(map.type)
  }
  
  if(is.null(fig) & new) {
    fig <- c(0,1,0,0.95)
    if (map.show & map.insert) fig[4] <- 0.8
  }
  if (legend.show) fig[3] <- 0.05  
  ## if (is.null(ylim))
  ##     if (is.null(dim(x)))
  ##         ylim <- pretty(x)
  ##     else              
  ##         ylim <- apply(x,2,pretty,n=5)
  if (is.null(xlim))
    xlim <- range(index(x))
  #if (is.null(ylim))
  #  ylim <- pretty(as.numeric(x))
  if (verbose) {print(xlim); print(ylim)}
  
  if (plot.type=="single") {
    if (is.null(ylab)) {
      ylab <- esd::ylab(x) # ggplot2 ylab can interfere with esd
    }
    if (inherits(ylab,"try-error")) ylab <- attr(x,'unit')
  } else {
    if (is.null(ylab)) { 
      if ((length(levels(factor(stid(x))))>1) & (length(levels(factor(varid(x))))<=1)) {
        ylab <- stid(x)
      } else 
        ylab <- varid(x)
    } else {
      if (is.null(main)) {
        if ((length(levels(factor(stid(x))))>1) & (length(levels(factor(varid(x))))<=1)) {
          main <- levels(factor((attr(x,'longname'))))[1]
        } else {
          main <- levels(factor(loc(x)))[1]
        }
      }
    }  
  }
  #if (is.null(main)) main <- attr(x,'longname')[1] 
  if (is.null(col)) {
    if (is.null(dim(x))) {
      col <- "blue"
    } else if (!is.null(lon(x)) & !is.null(lat(x)) &
               length(lon(x))==dim(x)[2] &
               length(lat(x))==dim(x)[2]) {
      nx <- (lon(x)-min(lon(x)))/diff(range(lon(x)))
      ny <- (lat(x)-min(lat(x)))/diff(range(lat(x)))
      if ( all(is.finite(nx) & is.finite(ny)) ) {
        col <- rgb(1-ny,nx,ny,1)
      } else {
        col <- rainbow(dim(x)[2])
      }
    } else {
      col <- rainbow(length(x[1,]))  
    }
  }
  if(is.null(alpha.map)) alpha.map <- alpha
  col.map <- adjustcolor(col,alpha.f=alpha.map)
  col <- adjustcolor(col,alpha.f=alpha)
  
  ns <- length(stid(x))
  #  if ( (ns > 1) & (plot.type=="multiple") ) {
  #    for (i in 1:ns) {
  #        z <- try(eval(parse(text=paste("ylab[",i,"] <- expression(",ylab[i],
  #                      "*phantom(0)*(",unit[i],"))"))),silent=TRUE)
  #        if (inherits(z,"try-error")) ylab[i] <- unit[i]
  #      }
  #  }
  
  errorbar <- errorbar & !is.null(err(x))
  
  if(map.show & !map.insert) {
    vis.map(x,col.map,map.type,add.text=FALSE,map.insert=map.insert,
            cex.axis=cex.axis,cex=1.8,usegooglemap=usegooglemap,
            zoom=zoom,verbose=verbose)
    new <- TRUE
  }
  
  #print(ylab)
  cls <- class(x)
  if("seasonalcycle" %in% cls) xaxt <- "n" else  xaxt <- NULL
  class(x) <- "zoo"
  if(new) dev.new()
  if(!is.null(fig)) par(cex.axis=1,fig=fig,mar=mar)
  par(bty="n",xaxt="s",yaxt="s",xpd=FALSE)
  plot.zoo(x,plot.type=plot.type,xlab=xlab,ylab=ylab,
           col=col,xlim=xlim,ylim=ylim,lwd=lwd,type=type,pch=pch,
           cex.axis=cex.axis,cex.lab=cex.lab,cex.main=cex.main,
           xaxt=xaxt,main=main,...)
  #mtext(main,side=3,line=1,adj=0,cex=cex.main)
  if("seasonalcycle" %in% cls) {
    axis(1,at=seq(1,12),labels=month.abb,cex.axis=cex.axis,las=2)
  }
  par0 <- par()
  if (plot.type=="single") {
    if (errorbar) {
      # REB 2014-10-03: add an errorbar to the plots.
      segments(index(x),x-err(x),index(x),x+err(x),
               lwd=3,col=rgb(0.5,0.5,0.5,0.25))
      #      d.err <- dim(err(x))
      #      dt <- 0.3*diff(index(x))[1]
      #      if (is.null(d.err)) d.err <- c(length(err(x),1)
      #      for (i in 1:d.err[2]) {
      #        for (j in 1:d.err[splot.dse1])
      #          lines(rep(index(x)[j],2),rep(x[j],2) + err(x)[j]*c(-1,1),
      #                lwd=3,col=rgb(1,0.5,0.5,0.25))
      #          lines(rep(index(x)[j],2) + dt*c(-1,1),rep(x[j],2) + err(x)[j],
      #                lwd=1,col=rgb(1,0.5,0.5,0.25))
      #          lines(rep(index(x)[j],2) + dt*c(-1,1),rep(x[j],2) - err(x)[j],
      #                lwd=1,col=rgb(1,0.5,0.5,0.25))
      #       }
    }
    
    par(fig=c(0,1,0,0.1),new=TRUE, mar=c(0,0,0,0),xaxt="s",yaxt="s",bty="n")
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    
    #if(legend.show) legend(0.01,0.75,loc(x),bty='n',ncol=4,
    #                       text.col=col,cex=0.75)
    #title(main=loc(x),cex=1)
    
    if(legend.show) {
      legend(0.01,0.95,paste(attr(x,'location'),": ",
                             #attr(x,'aspect'),
                             #attr(x,'longname')," - ",
                             round(attr(x,'longitude'),2),"E/",
                             round(attr(x,'latitude'),2),"N (",
                             attr(x,'altitude')," masl)",sep=""),
             bty="n",cex=0.6,ncol=3,text.col="grey40",lty=1,col=col)
    }
    if (map.show & map.insert) vis.map(x,col.map,map.type=map.type,cex=1,
                                       cex.axis=0.65,add.text=FALSE,
                                       map.insert=map.insert,usegooglemap=usegooglemap,
                                       zoom=zoom,verbose=verbose)
    par(fig=par0$fig,mar=par0$mar,new=TRUE)
    plot.zoo(x,plot.type=plot.type,type="n",xlab="",ylab="",
             xaxt="n",yaxt="n",xlim=xlim,ylim=ylim,new=FALSE)
    
  }
}
