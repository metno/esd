## Plots one time series
plot.dsensemble.one <-  function(x,pts=FALSE,it=0,
                                 envcol=rgb(1,0,0,0.2),legend.show=TRUE,ylab=NULL,
                                 obs.show=TRUE,target.show=TRUE,map.show=TRUE,map.type=NULL,map.insert=TRUE,
                                 usegooglemap=TRUE,new=FALSE,xrange=NULL,yrange=NULL,
                                 alpha=0.5,alpha.map=0.7,mar=c(5.1,4.5,4.1,2.1),
                                 cex.axis=1, cex.lab=1.2, cex.main=1.2, 
                                 verbose=FALSE,...) {
  if(verbose) print("plot.dsensemble.one")
  stopifnot(inherits(x,'dsensemble'))
  
  if (is.null(map.type)) {
    if (verbose) print(class(x))
    if( inherits(x,"field") | length(lon(x))!=length(lat(x)) |
        (length(lon(x))==2 & length(lat(x))==2) ) {
      map.type <- "rectangle"
    } else {
      map.type <- "points"
    }
  }
  
  if (verbose) {print(map.type); print(attr(x,'station'))}
  if (!is.null(attr(x,'station')) & !inherits(attr(x,'station'),c('annual','season'))) {
    z <- subset(x,it=it,verbose=verbose) 
  } else {
    z <- x
  }
  
  if (verbose) {
    print("diagnose")
    class(z)
  }
  diag <- diagnose(z,plot=FALSE,verbose=verbose)
  
  y <- attr(z,'station')
  attr(y,'standard.error') <- NULL
  if (verbose) print(paste('lon=',lon(y),'lat=',lat(y))) 
  
  d <- dim(z)
  index(y) <- year(y)
  
  if(map.show | target.show) {
    pscl <- c(0.9,1.3)
  } else {
    pscl <- c(0.9,1.1)
  }
  
  if (max(coredata(z),na.rm=TRUE) < 0) pscl <- rev(pscl)
  args <- list(...)
  if (verbose) print(names(args))
  ixl <- grep('xlim',names(args))
  if (length(ixl)==0) xlim <- range(year(z)) else
    xlim <- args[[ixl]]
  iyl <- grep('ylim',names(args))
  if (length(iyl)==0) ylim <- pscl*range(coredata(z),na.rm=TRUE) else
    ylim <- args[[iyl]]  
  #print("...")
  if(new) dev.new()
  index(y) <- year(y)
  if(!is.null(mar)) par(mar=mar)
  par0 <- par()
  if (obs.show) obscol <- 'black' else obscol='white'
  plot(y,type="b",pch=19,xlim=xlim,ylim=ylim,col=obscol,main='',
       cex.axis=cex.axis,cex.lab=cex.lab,
       ylab=ylab,map.show=FALSE,new=new, verbose=verbose)
  grid()
  usr <- par()$usr; mar <- par()$mar; fig <- par()$fig
  t <- index(z)
  
  if (pts) for (i in 1:d[2]) {
    points(year(t),coredata(z[,i]),pch=19,col="red",cex=0.3)
  }
  
  # Produce a transparent envelope
  nt <- length(index(z))
  t2 <- c(year(t),rev(year(t)))
  
  col <- rgb(rep(1,49),seq(0.95,0.1,length=49),seq(0.95,0.1,length=49),0.1)
  ## REB 2016-11-25
  if(is.null(alpha.map)) alpha.map <- alpha
  col.map <- adjustcolor(col,alpha.f=alpha.map)
  col <- adjustcolor(col,alpha.f=alpha)
  
  mu <- apply(coredata(z),1,mean,na.rm=TRUE)
  si <- apply(coredata(z),1,sd,na.rm=TRUE)
  for (ii in 1:49) {
    qp1 <- qnorm(1-ii/50,mean=coredata(mu),sd=coredata(si))
    qp2 <- qnorm(ii/50,mean=coredata(mu),sd=coredata(si))
    ci <- c(qp1,rev(qp2))
    polygon(t2[!is.na(ci)],ci[!is.na(ci)], col= envcol, border=NA)
  }
  q05 <- qnorm(0.05,mean=mu,sd=si)
  q95 <- qnorm(0.95,mean=mu,sd=si)
  
  lcol <- adjustcolor(envcol,offset=c(0.5,0.5,0.5,0.2))
  lines(zoo(mu,order.by=year(z)),lwd=3,col=lcol)
  lines(zoo(q05,order.by=year(z)),lty=2,col=lcol)  
  lines(zoo(q95,order.by=year(z)),lty=2,col=lcol)
  if (obs.show) lines(y,type="b",pch=19)
  
  if (!is.null(diag)) {
    index(diag$y) <- year(diag$y)
    outside <- diag$above | diag$below
    points(zoo(coredata(diag$y)[which(outside)],
               order.by=year(diag$y)[which(outside)]),col="grey")
  }
  
  title(main=toupper(loc(x)),cex.main=cex.main)
  if ((target.show) & (!is.null(diag))) {
    if (verbose) print('add target diagnostic')
    par(fig=c(0.23,0.45,0.75,0.95),new=TRUE, mar=c(0,0,0,0),xaxt="s",yaxt="n",bty="n",
        cex.main=0.75,xpd=NA,col.main="grey30")
    plot(diag,map.show=FALSE,new=FALSE,cex=0.75)
  } 
  
  if(map.show & !map.insert) {
    vis.map(x,"red",map.type,add.text=FALSE,map.insert=map.insert,
            cex.axis=cex.axis,cex=1.5,usegooglemap=usegooglemap,
            xrange=xrange,yrange=yrange,
            verbose=verbose,...)
    new <- TRUE
  }
  # REB 2016-11-25
  #if (map.show) {
  #  if(verbose) print("add map")
  #  if(is.null(xrange) & !is.null(lon(y))) {
  #    xrange <- range(lon(y)) + c(-15,15)
  #  }
  #  if(is.null(yrange) & !is.null(lat(y))) {
  #    yrange <- range(lat(y)) + c(-10,10)
  #  }
  #  if (!is.null(xrange) & !is.null(xrange)) {
  #    data("geoborders", envir = environment())
  #    lon <- geoborders$x
  #    lat <- geoborders$y
  #    lon2 <- attr(geoborders,"borders")$x
  #    lat2 <- attr(geoborders,"borders")$y
  #    par(fig=c(0.7,0.95,0.78,0.98),new=TRUE, mar=c(0,0,0,0),
  #      cex.main=0.75,xpd=FALSE,col.main="grey",bty="n")
  #    plot(lon,lat,lwd=1,col="black",type='l',xlab=NA,ylab=NA,
  #       axes=FALSE,xlim=xrange,ylim=yrange)
  #    axis(1,mgp=c(3,.5,0))
  #    axis(2,mgp=c(2,.5,0))
  #    lines(lon2,lat2,col = "pink",lwd=1)
  #    if("points" %in% map.type) {
  #      points(lon(y),lat(y),pch=21,cex=1,col='black',bg='red',lwd=1)
  #  }
  #    if("rectangle" %in% map.type) {
  #      rect(min(lon(y)),min(lat(y)),max(lon(y)),max(lat(y)),lwd=1,col=NA,border='red',lty=2)
  #    }
  #  } else if (verbose) print(paste('lon=',lon(y),'lat=',lat(y)))
  #}  
  # finished plotting
  
  if (legend.show) {
    par(fig=c(0.1,0.5,0.2,0.25),new=TRUE,mar=c(0,0,0,0),xaxt="n",yaxt="n",bty="n",xpd=NA)
    #par(fig=c(0.1,0.5,0.65,0.70),new=TRUE, mar=c(0,0,0,0),xaxt="n",yaxt="n",bty="n")
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    ## KMP 2017-09-13: move this text!!! vvv
    if(!is.null(diag)) {
      legend(0.05,0.90,c(paste("Past trend:",round(diag$deltaobs,2)),
                         paste(diag$robs,'%'),
                         paste(diag$outside,"observations"),
                         "p-value: "),
             bty="n",cex=0.7,text.col="grey40")
      #legend(0.5,0.90,c(expression(paste(levels(factor(attr(x,'unit')))[1]/d*e*c*a*d*e)),
      #                  "ensemble trends > obs.",
      #                  "outside ensemble 90% conf.int.",
      #                  paste(round(100*pbinom(diag$outside,size=diag$N,prob=0.1)),"%")),
      #        bty="n",cex=0.7,text.col="grey40")
      legend(0.5,0.90,c(paste(levels(factor(attr(y,'unit')))[1],"/decade",sep=""),
                        "ensemble trends > obs.",
                        "outside ensemble 90% conf.int.",
                        paste(round(100*pbinom(diag$outside,size=diag$N,prob=0.1)),"%")),
             bty="n",cex=0.7,text.col="grey40")
    }
  }
  if (map.show & map.insert) vis.map(x,"red",map.type=map.type,cex=1.5,
                                     cex.axis=cex.axis*0.65,add.text=FALSE,
                                     map.insert=map.insert,usegooglemap=usegooglemap,
                                     xrange=xrange,yrange=yrange,
                                     verbose=verbose,...)
  par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
      fig=c(0,1,0.1,1),new=TRUE)
  par(fig=fig,new=TRUE, mar=mar)
  plot(usr[1:2],usr[3:4],
       type="n",ylab="",xlab="",xlim=usr[1:2],ylim=usr[3:4])
  
  # target: perfect score is bull's eye
  # from diagnostics?
  par(fig=par0$fig,new=TRUE, mar=par0$mar)
  if(verbose) print("exit plot.dsensemble.one")
  invisible(z)
}
