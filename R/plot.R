plot <- function(x,y, ...)  UseMethod("plot")

plot.list <- function(x,is=NULL,
                      col=c(rgb(1,1,0.5,0.05),rgb(1,0.5,0.5,0.05),rgb(0.5,1,0.5,0.05)),
                      lwd=3,xlim=NULL,ylim=NULL,...) {
  if (!is.null(is)) y <- subset(x,it=is) else y <- x[[1]]
  plot(y,img=img,col=col[1],lwd=lwd,xlim=xlim,ylim=ylim)
  for (j in c(2:length(x),1)) {
    if (!is.null(it)) y <- subset(x[[j]],it=it) else y <- x[[j]]
    for (i in 1:dim(y)[2]) lines(y[,i],lwd=7,col=col[j])
    lines(attr(y,'station'),lwd=3,col=rgb(0.5,0.5,0.5,0.25))
  }
}

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
    if (is.null(ylab))
      ylab <- ylab(x)
    if (inherits(ylab,"try-error")) ylab <- unit(x)
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
  ##browser()
  plot.zoo(x,plot.type=plot.type,xlab=xlab,ylab=ylab,
           col=col,xlim=xlim,ylim=ylim,lwd=lwd,type=type,pch=pch,
           cex.axis=cex.axis,cex.lab=cex.lab,xaxt=xaxt,main=main,...)
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


vis.map <- function(x,col='red',map.type=NULL,
                    xrange=NULL,yrange=NULL,cex=1,
                    add.text=FALSE,cex.axis=NULL,
                    map.insert=TRUE,verbose=FALSE,
                    usegooglemap=TRUE,zoom=NULL,...) {
  if(verbose) {print('vis.map'); print(lon(x)); print(lat(x)); print(zoom)}
  ## KMP 2017-06-07 Weird problem: cex.axis is not found even though it is an argument to the function.
  ## It looks like cex.axis exists but when applying 'print' the following error message shows up: 
  ## 'Warning: restarting interrupted promise evaluation. Error in print(cex.axis) : object 'cex.axis' not found'
  cex.axis <- 0.7  # Temporary fix
  
  ## REB 2016-11-25: for dsensemble object
  if (is.null(lon(x))) attr(x,'longitude') <- lon(attr(x,'station'))
  if (is.null(lat(x))) attr(x,'latitude') <- lat(attr(x,'station'))
  if(is.null(xrange)) xrange <- range(lon(x)) + c(-5,5)
  if(is.null(yrange)) yrange <- range(lat(x)) + c(-2,2)
  if(!map.insert) new <- TRUE else new <- FALSE
  
  if (is.null(map.type)) {
    if( inherits(x,"field") | length(lon(x))!=length(lat(x)) |
        (length(lon(x))==2 & length(lat(x))==2) ) {
      map.type <- "rectangle"
    } else {
      map.type <- "points"
    }
  }
  
  ## REB: 2016-10-12 - add the possibility to use google maps
  ## KMP 2018-10-31: Don't use require inside the esd package. 
  ## Instead check if it the external package is installed and then 
  ## call it explicitly, e.g., RgoogleMaps::GetMap().
  ## Also add the package under 'Suggested' in the DESCRIPTION file.
  if (!requireNamespace("RgoogleMaps", quietly = TRUE)) {
    usegooglemap <- FALSE
  }
  
  if(usegooglemap) {
    if (is.null(zoom)) {
      if (verbose) print('zoom not defined')
      if (length(lon(x))==1) {
        zoom <- 5 
      } else {
        ## zoom = 12 is very local, zoom = 1 is the world
        mxdst <- max(diff(range(lat(x))),diff(range(lon(x))))
        zoom <- 1 - floor(0.75*log(mxdst/360))
      }
    }
    if (!is.finite(zoom)) zoom <- 5
    if (verbose) print(paste('zoom=',zoom))
    bgmap <- try(RgoogleMaps::GetMap(center=c(lat=mean(lat(x)),lon=mean(lon(x))),
                  destfile = "map.station.esd.png",
                  maptype = "mobile", zoom=zoom))
    if(inherits(bgmap,"try-error")) {
      usegooglemap <- FALSE
    } else {
      if(map.insert) {
        par(fig=c(0.75,0.95,0.75,0.95),new=TRUE,
             mar=c(0,0,0,0),xpd=NA,col.main="grey",bty="n")
      }
      if(map.type=="rectangle") {
        xx <- c(rep(max(lat(x)),2), rep(min(lat(x)),2), max(lat(x)))
        yy <- c(range(lon(x)), rev(range(lon(x))), min(lon(x)))
        RgoogleMaps::plotmap(xx, yy, bgmap, pch=19, col=col, cex=0.25)
        RgoogleMaps::PlotOnStaticMap(bgmap, lat=xx, lon=yy, lwd=1, col=col, FUN=lines, add=TRUE)
      } else {
        RgoogleMaps::plotmap(lat(x), lon(x), bgmap, pch=19, col=col, cex=2)
      }
    }
  }
  
  if(!usegooglemap) {
    if (verbose) {
      print('basic map')
      print(cex.axis)
    }
    data("geoborders", envir = environment())
    lon <- geoborders$x
    lat <- geoborders$y
    ok <- lon>(min(xrange)-1) & lon<(max(xrange)+1) &
          lat>(min(yrange)-1) & lat<(max(yrange)+1) &
          is.finite(lon) & is.finite(lat)
    lon2 <- attr(geoborders,"borders")$x
    lat2 <- attr(geoborders,"borders")$y
    ok2 <- lon2>(min(xrange)-1) & lon2<(max(xrange)+1) &
           lat2>(min(yrange)-1) & lat2<(max(yrange)+1) &
          is.finite(lon2) & is.finite(lat2)
    if (verbose) {print(sum(ok)); print(range(lon[ok])); print(range(lat[ok]))}
    if(map.insert) {
      par(fig=c(0.76,0.97,0.76,0.97),new=TRUE,
          mar=c(0,0,0,0),xpd=NA,col.main="grey",bty="n")
    } else {
      dev.new()
    }
    plot(lon[ok],lat[ok],lwd=1,col="black",type="p",pch='.',cex=2,
        #type='l', KMP 2016-03-16 problem with lines in map
        xlab=NA,ylab=NA,axes=FALSE,new=new,
        xlim=xrange,ylim=yrange)
      #xlim=range(c(lon[ok],lon2[ok2]),na.rm=TRUE),
      #ylim=range(c(lat[ok],lat2[ok2]),na.rm=TRUE))
    par(xpd=FALSE)
    lines(lon,lat) ## REB: 2016-11-25 need more solid lines.
    axis(1,mgp=c(3,0.5,0.3),cex.axis=cex.axis)
    axis(2,mgp=c(2,0.5,0.3),cex.axis=cex.axis)
    lines(lon2,lat2,col = "pink",lwd=1)
    #lines(lon2[ok2],lat2[ok2],col = "pink",lwd=1)
    if (verbose) print(map.type)
    if (map.type=="points") {
      if (verbose) {print(c(lon(x),lat(x),cex)); print(col)}
      points(lon(x),lat(x),pch=21,cex=cex,col=col,bg=col,lwd=1)
      if (add.text) text(lon(x),lat(x),labels=loc(x),col=col) 
    } else if (map.type=="rectangle") {
      rect(min(lon(x)),min(lat(x)),max(lon(x)),max(lat(x)),
            border="black",lwd=1,lty=2)
    }
  }
  if(verbose) print("exit vis.map")
}

plot.eof <- function(x,new=FALSE,xlim=NULL,ylim=NULL,
                     ip=1,what=c("pc","eof","var"),
                     colbar=list(pal=NULL,rev=FALSE,n=10,alpha=0.8,
                         breaks=NULL,type="p",cex=2,show=TRUE,
                         h=0.6,v=1,pos=0.05),
                     verbose=FALSE,is=NULL,it=NULL,...) {
  if (verbose) print(paste('plot.eof',paste(what,collapse=',')))
  if (inherits(x,"comb"))
    plot.eof.comb(x,new=new,xlim=xlim,ylim=ylim,
                  ip=ip,what=what,colbar=colbar,verbose=verbose,...) else
  if (inherits(x,"field"))
    plot.eof.field(x,new=new,xlim=xlim,ylim=ylim,
                   ip=ip,what=what,colbar=colbar,verbose=verbose,
                   it=it,is=is,...) else
    print("x does not have 'comb' or 'field' aspects...")
}




plot.eof.field <- function(x,new=FALSE,xlim=NULL,ylim=NULL,ip=1,
                           what=c("pc","eof","var"),## colbar=NULL,
                           cex.axis=0.9,cex.main=0.9,cex.lab=0.9,
                           verbose=FALSE,it=NULL,is=NULL,cex=1,...) {
  if (verbose) print(paste('plot.eof.field',paste(what,collapse=',')))
  n <- ip
  what <- tolower(what)
  if ('field' %in% what) {
      ## Expand EOF to original field before plotting
      if (verbose) print('Transform eof to field before plot')
      x <- subset(x,it=it,is=is)
      y <- as.field(x)
      z <- plot(y,xlim=xlim,ylim=ylim,new=new,...)
      invisible(z)
  }
  #str(ip); stop("HERE")
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100* D^2/tot.var
  if (length(what)==3) {
      mfrow <- c(2,2)
      
  } else
  if (length(what)==2) mfrow <- c(2,1) else
  if (length(what)==1) mfrow <- c(1,1)
  if (new) dev.new()
  ## par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)
  par(mfrow=mfrow)##,mar=c(1,1,1,2)) ##,bty="n",xaxt="n",yaxt="n")
  if (length(grep('eof',what))>0) {
      if (verbose) {print('Show map'); print(class(x))}
      if (inherits(x,'eof')) {## inherits(x,'pca') |
          par(fig=c(0,0.5,0.5,1))
          ## par(fig=c(0.025,0.5,0.5,0.975)) ## c(0,0.45,0.5,0.975) c(0.05,0.5,0.55,0.95)
          map(x,ip=ip,verbose=verbose,
              cex.main=cex.main,cex.axis=cex.axis,
              cex.lab=cex.lab,cex=cex,...) ## AM formely new=FALSE colbar=colbar,
      } else if (inherits(x,'pca')) {
          fig <- c(0,0.5,0.5,1)
          par(fig=fig)
          main1 <- paste('Leading EOF#',ip, ' (',
                         round(var.eof[ip],digits=2),"%)",sep='')
          map(x,ip=ip,verbose=verbose,
              cex.main=cex.main,cex.axis=cex.axis,
              cex.lab=cex.lab,cex=cex,fig=fig,...) ## colbar=colbar,
          title(main=src(x)[1],cex.main=cex.main*0.8,
                col.main="grey40",adj=0,line=0)
          title(main=main1,cex.main=cex.main)
      }
  }
  ##  if (length(grep('pc',what))>0) result <- as.station(x) else
#  if (length(grep('var',what))>0) result <- attr(x,'tot.var')
    
  ylab <- paste("PC",n)
  main <- paste('First',n,"leading EOFs: ", ## attr(x,'longname')
                 round(sum(var.eof[1:n]),1),"% of variance")
  
  if (length(grep('var',what))>0) {
    par(new=TRUE,fig=c(0.5,1,0.5,1))##,xaxt="s",yaxt="s")fig=c(0.5,0.95,0.5,0.975) 
    plot.eof.var(x,ip=ip,new=FALSE,cex.main=cex.main,
                 cex.axis=cex.axis,bty="n",cex=cex)
  }
  
  #print(main)
  if (length(grep('pc',what))>0) {
    ##par(bty="n", ##,xaxt="s",yaxt="s",xpd=FALSE,
      par(fig=c(0.05,1,0.025,0.475),new=TRUE) ##,cex.axis=0.9,cex.lab=1) ##(0.05,0.95,0.02,0.45)
      main <- paste('Leading PC#',ip,' of ',attr(x,'longname'),
                 " - Explained variance = ",round(var.eof[ip],digits=2),
                    "%",sep='')

      if(inherits(x,"seasonalcycle")) xaxt <- "n" else  xaxt <- NULL
      plot.zoo(x[,n],lwd=2,ylab=ylab,main=main,xlim=xlim,ylim=ylim,
               cex.main=cex.main,bty="n",cex.axis=cex.axis,
               cex.lab=cex.lab,xaxt=xaxt)
      if(inherits(x,"seasonalcycle")) axis(1,at=seq(1,12),labels=month.abb,
                                           cex.axis=cex.axis,las=2)
      #axis(1,at=pretty(index(x),n=10),labels=,cex.axis=0.9)
      grid()
  }
 
  par(fig=c(0,1,0,0.55),new=TRUE, mar=c(1,1,1,1),xaxt="n",yaxt="n",bty="n")
  plot(c(0,1),c(0,1),type="n",xlab="",ylab="")

  varnm <- varid(x)[1]
  legend(0,0.83,varnm,bty="n",cex=0.8,ncol=2,text.col="grey40")
  
  par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
      fig=c(0,1,0.1,1),new=FALSE)
  #par(fig=c(0,1,0,0.1),new=NEW, mar=c(0,0,0,0))  
}


plot.eof.comb <- function(x,new=FALSE,xlim=NULL,ylim=NULL,
                          ip=1,col=c("red"),alpha=1,
                          what=c("pc","eof","var"),colbar=NULL,verbose=FALSE,...) {
  if (verbose) print("plot.eof.comb")
  n <- ip
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100* D^2/tot.var

  if (length(what)==3) mfrow <- c(2,2) else
  if (length(what)==2) mfrow <- c(2,1) else
                       mfrow <- NULL
  
  if (new) dev.new()
  #par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)
  if (!is.null(mfrow)) par(mfrow=mfrow)

  if (length(grep('eof',what))>0) {
    if (!is.null(mfrow)) par(fig=c(0,0.5,0.5,1))
    map(x,ip=ip,verbose=verbose,colbar=colbar,...)
  }

  n.app <- attr(x,'n.apps')
  col <- rep(col,n.app)
  src <- rep("",n.app+1)
  src[1] <- attr(x,'source')
  ylab <- paste("PC",n)
  main <- paste("EOF: ",n,"accounts for",
                round(var.eof[n],1),"% of variance")
  
  if (length(grep('var',what))>0)  {
#    par(xaxt="s",yaxt="s")
#    plot.eof.var(x,new=FALSE,cex.main=0.7)
    if (!is.null(mfrow)) par(new=TRUE,fig=c(0.5,1,0.5,1))##,xaxt="s",yaxt="s")fig=c(0.5,0.95,0.5,0.975) 
    plot.eof.var(x,ip=ip,new=FALSE,cex.main=0.8,cex.axis=0.9,bty="n")
  }

  if (is.null(ylim)) {
    ylim <- range(coredata(x[,n]))
    for (i in 1:n.app) {
      z <- attr(x,paste('appendix.',i,sep=""))
      ylim <- range(c(ylim,coredata(z[,n])),na.rm=TRUE)
    }  
  }

  if (is.null(xlim)) {
    xlim <- range(index(x))
    for (i in 1:n.app) {
      z <- attr(x,paste('appendix.',i,sep=""))
      xlim <- range(xlim,index(z))
    }
  }
  if(is.character(xlim)) xlim <- as.Date(xlim)

  if (length(grep('pc',what))>0) {
    if (verbose) {print('time series');print(index(x)); print(index(attr(x,'appendix.1')))}
    if ( (sd(coredata(x)[,n])/sd(coredata(attr(x,'appendix.1'))[,n]) > 100) |
         (sd(coredata(attr(x,'appendix.1'))[,n])/sd(coredata(x)[,n]) > 100) )
       warning('plot.comb.eof: PCs have very different scales')
#    par(bty="n",xaxt="s",yaxt="s",xpd=FALSE,
#      fig=c(0.1,0.9,0.1,0.5),new=TRUE,cex.axis=0.6,cex.lab=0.6)
#    plot.zoo(x[,n],lwd=2,ylab=ylab,main=main,sub=attr(x,'longname'),
#                                          xlim=xlim,ylim=ylim)
    if (!is.null(mfrow)) par(fig=c(0.025,1,0.025,0.475),new=TRUE) ##,cex.axis=0.9,cex.lab=1) ##(0.05,0.95,0.02,0.45)
      main <- paste('Leading PC#',ip,'of ',attr(x,'longname'),
                 " - Explained variance = ",round(var.eof[ip],digits=2),
                    "%",sep='')
      
      plot.zoo(x[,n],lwd=2,ylab=ylab,main=main,xlim=xlim,ylim=ylim,
               cex.main=0.8,bty="n",cex.axis=0.9,cex.lab=1,xaxt="n")
      taxis <- pretty(index(x[,n]),n=10)              # REB 2016-03-03
      if (min(diff(taxis))> 360) taxisl <- year(taxis)  else
                                 taxisl <- taxis      # REB 2016-03-03
      if (verbose) print(taxisl)
      axis(1,at=taxis,labels=taxisl,cex.axis=0.9)      # REB 2016-03-03
      grid()

      ## Plot the common PCs
      for (i in 1:n.app) {
        z <- attr(x,paste('appendix.',i,sep=""))
        lines(z[,n],col=adjustcolor(col[i],alpha.f=alpha),lwd=2)
        if (verbose) print(attr(z,'source'))
        if (!is.null(attr(z,'source'))) src[i+1] <- attr(z,'source') else
                                        src[i+1] <- paste('x',i,sep='.')
      }

      lines(x[,n],lwd=2,col="black")
    }
#    par(xaxt="n",yaxt="n",bty="n",fig=c(0,1,0,0.1),
#        mar=rep(0,4),new=TRUE)
#    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
#    legend(0,1,src,col=c("black",col),lwd=2,ncol=4,bty="n",cex=0.7)
#    par(xaxt="n",yaxt="n",bty="n",fig=par0$fig,mar=par0$mar,new=TRUE)
#    plot.zoo(x[,n],type="n",xlab="",ylab="")
  
  par(fig=c(0,1,0,0.55),new=TRUE, mar=c(0,0,0,0),xaxt="n",yaxt="n",bty="n")
  plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
  varnm <- varid(x)
  legend(0,0.83,varnm,bty="n",cex=0.8,ncol=2,text.col="grey40")
  
  par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
      fig=c(0,1,0.1,1),new=TRUE)
  par(fig=c(0,1,0,0.1),new=TRUE, mar=c(0,0,0,0))  
}


plot.ds <- function(x,plot.type="multiple",what=c("map","ts",'xval'),new=TRUE,
                    lwd=1,type='l',pch=0,main=NULL,col=NULL,
                    colbar=list(pal=NULL,rev=FALSE,n=10,
                        breaks=NULL,type="p",cex=2,show=TRUE,
                        h=0.6, v=1,pos=0.05),
                    xlim=NULL,ylim=NULL,xlab="",ylab=NULL,verbose=FALSE,...) {
  if (verbose) print(paste('plot.ds',paste(what,collapse=',')))
  if (inherits(x,'pca')) {
    plot.ds.pca(x,verbose=verbose,...)
    return()
  }
  if (inherits(x,'eof')) {
    plot.ds.eof(x,verbose=verbose,...)
    return()
  }
  
  unit <- attr(x,'unit')
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
  for (i in 1:length(unit)) {
    if ((unit[i]=='degree Celsius') | (unit[i]=='deg C') | (unit[i]=='degC'))
         unit[i] <- 'degree*C'
  }
  
  if (is.null(ylab))
   ylab <- ylab(x)
  
  if (verbose)  print(ylab)
  if (is.null(main)) main <- attr(x,'longname')[1]               
  if (is.null(col)) col <- rainbow(length(x[1,]))  
  
  cols <- rep("blue",100)
  model <- attr(x,'model')
  
  if (length(what)==2) mfrow <- c(2,1) else
  if (length(what)==1) mfrow <- c(1,1)
  
  if (new) dev.new()
  if (plot.type=="single") new <- TRUE
  par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)

  if (sum(is.element(what,'map'))>0) {
    if (verbose) print('Show map...')
    par(fig=c(0,0.5,0.5,1))
    map(x,new=FALSE,colbar=list(show=FALSE),verbose=verbose,...)
    points(lon(x),lat(x),lwd=3,cex=1.5)
  }

  if ( (sum(is.element(what,'xval'))>0)  & (!is.null(attr(x,'evaluation'))) ){
    par(new=TRUE,fig=c(0.5,1,0.5,1)) 
     
    plot(attr(x,'evaluation')[,1],attr(x,'evaluation')[,2],
         main='Cross-validation',xlab='original data',
         ylab='prediction',pch=19,col="grey")
    lines(range(c(attr(x,'evaluation')),na.rm=TRUE),
          range(c(attr(x,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(attr(x,'evaluation')[,1]),
                      x=coredata(attr(x,'evaluation')[,2]))
    xvalfit <- lm(y ~ x, data = cal)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    par(bty="n",fig=c(0.6,0.95,0.48,0.52),mar=c(0,0,0,0),new=TRUE,
        xaxt='n',yaxt='n',cex.sub=0.7)
    plot(c(0,1),c(0,1),type='n',xlab='',ylab='')
    ok <- is.finite(attr(x,'evaluation')[,1]) &
          is.finite(attr(x,'evaluation')[,2])
    text(par()$xaxp[1],mean(par()$yaxp[1:2]),
          paste('correlation=',
          round(cor(attr(x,'evaluation')[ok,1],attr(x,'evaluation')[ok,2]),2)),
         pos=4,cex=0.8,col='grey')
  }  else {
    xvalfit <- NULL
  }

  Y0 <- as.original.data(x)
  #print(index(Y0)); print(index(x))
  yX <- merge.zoo(Y0,x,all=FALSE)
  #print(summary(yX))
  y0 <- yX$Y0

  if (!is.null(attr(x,'n.apps'))) ns <- attr(x,'n.apps') else
                                  ns <- 0
  y.rng <- NA; x.rng <- NA
  if (ns > 0) {
    #print("Add other DS results")
    for (i in 1:ns)
      eval(parse(text=paste("y <- attr(x,'appendix.",i,"')",sep="")))
      #print(summary(y))
      y.rng <- range(y.rng,y,na.rm=TRUE)
      x.rng <- range(x.rng,index(y),na.rm=TRUE)
  }
  
  if (is.null(ylim)) {
    ylim <- range(coredata(x),coredata(y0),y.rng,na.rm=TRUE)
  }
  if (is.null(xlim)) {
    xlim <- range(index(x),index(y0),x.rng,na.rm=TRUE)
  }

  par(fig=c(0.025,1,0.025,0.475),new=TRUE)
  par(bty="n",fig=c(0,1,0.1,0.5),mar=c(1,4.5,1,1),new=TRUE, xaxt='s',yaxt='s')
  ds <- list(obs=y0)
  plot.zoo(y0,plot.type=plot.type,ylab=ylab,xlab=xlab,
           main=main,xlim=xlim,ylim=ylim,lwd=1,type='b',pch=19)
  par0 <- par()
  grid()
  if (verbose) print(c(class(index(x)),class(index(y0))))
  if ( (class(index(x))=='Date') & (class(index(y0))=='numeric') & inherits(x,'annual') ) 
    index(x) <- year(index(x))
  if ( (class(index(x))=='numeric') & (class(index(y0))=='Date') & inherits(x,'annual') ) 
    index(x) <- as.Date(paste(index(x),'01-01',sep='-'))
  lines(x,col="red",type="l",lwd=lwd)
  
  cal0 <- data.frame(y=coredata(y0),t=year(y0))
  cal1 <- data.frame(y=coredata(x),t=year(x))
 
  trend0 <- lm(y ~ t, data=cal0)
  trend1 <- lm(y ~ t, data=cal1)
  lines(zoo(predict(trend0),order.by=index(y0)),lty=2)
  lines(zoo(predict(trend1),order.by=index(x)),lty=2,col='red')

  st0 <- summary(trend0); st1 <- summary(trend1)
  obstrend <- paste('obs. trend: ', round(st0$coefficients[2],2),' (',
                    round(st0$coefficients[2]-2*st0$coefficients[4],2),',',
                    round(st0$coefficients[2]+2*st0$coefficients[4],2),')',
                    unit(x),'/decade',sep='')
  dstrend <- paste('obs. trend: ', round(st1$coefficients[2],2),' (',
                    round(st1$coefficients[2]-2*st1$coefficients[4],2),',',
                    round(st1$coefficients[2]+2*st1$coefficients[4],2),')',
                    unit(x),'/decade',sep='')
  
  if (is.null(attr(x,'source'))) attr(x,'source') <- 'ESD'
  if (is.na(attr(x,'source'))) attr(x,'source') <- 'ESD'
  legtext <- c("Observations",attr(x,'source')) 
  legcol <- c("black","red")

  if (ns > 0) {
    #print("Add other DS results")
    for (i in 1:ns) {
      eval(parse(text=paste("y <- attr(x,'appendix.",i,"')",sep="")))
      lines(zoo(coredata(y),order.by=index(y)),col=cols[i],lwd=lwd)
      legcol <- c(legcol,cols[i])
      legtext <- c(legtext,attr(y,'source'))
      eval(parse(text=paste("ds$result.",i,
                   " <- attr(x,'appendix.",i,"')",sep="")))
    }
    
  } else  legcol <- c("black","red")

  ## Replot observations and prediction for calibration period

  lines(y0,lwd=1,type='b',pch=19)
  lines(x,col="red",type="l",lwd=lwd)
  #print(legcol)
  if (!is.null(attr(x,'appendix.1'))) legend <- c("Obs.","Cal.","Proj") else
                                      legend <- c("Obs.","Cal.")
  legend(x="topleft",legend=legend,bty="n",horiz=TRUE,
         col=c("black","red","blue"),lwd=c(1,1,1),pch=c(19,1,1))
  
  if (plot.type=="single") {
    par(fig=c(0,1,0,0.1),new=TRUE, mar=c(0,0,0,0),xaxt="s",yaxt="s",bty="n")
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    legend(0.01,0.95,c(paste(attr(x,'location'),": ",
                           #attr(x,'aspect'),
                           #attr(x,'longname')," - ",
                           round(attr(x,'longitude'),2),"E/",
                           round(attr(x,'latitude'),2),"N (",
                           attr(x,'altitude')," masl)",sep=""),
                       obstrend,dstrend),ncol=3,
           bty="n",cex=0.6,text.col="grey40",
           lwd=c(1,rep(2,ns),1,1),lty=c(1,rep(1,ns),2,2),
           col=c(col,'black',col[1]))
    
    par(fig=c(0,1,0.05,0.95),new=TRUE,mar=par0$mar,xaxt="n",yaxt="n",bty="n")
    plot.zoo(x,plot.type=plot.type,type="n",ylab="",xlab="",xlim=xlim,ylim=ylim)
  }
  invisible(list(trend0=trend0,trend1=trend1,xvalfit=xvalfit))
}



plot.eof.var <- function(x,ip=1,new=TRUE,xlim=NULL,ylim=NULL,n=20,...) {
  n <- min(c(n,length(attr(x,'eigenvalues'))))
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100* D^2/tot.var
  nt <- length(index(x)) 
  n.eff <- round(nt * (1.0-attr(x,'max.autocor'))/
                      (1.0+attr(x,'max.autocor')))  
  dD <- D*sqrt(2.0/n.eff)
  main <- paste(n,"leading EOFs: ", ## attr(x,'longname')
                 round(sum(var.eof[1:n]),digits=1),"% of variance")
  ##main <- paste(attr(x,'longname'),n,"leading EOFs: ",
  ##               round(sum(var.eof[1:n]),1),"% of variance")
  if (is.null(xlim)) xlim <- c(0,n) ##c(0.7,n+0.3)
  if (is.null(ylim)) ylim <- c(0,100)
  if (new) dev.new()
  ##par(bty="n") ##,xaxt="n")
  plot(var.eof,type="n",
       main=main,ylab="Variance (%)",xlab="EOF order",
       xlim=xlim,ylim=ylim,...)
  lines(cumsum(var.eof),lwd=3,col=rgb(1,0.8,0.8))
  grid() ## nx=21,ny=12)
  lines(var.eof,type="b",pch=19)
  for (i in 1:length(var.eof)) {
    lines(rep(i,2),100*c((D[i]+dD[i])^2/tot.var,(D[i]-dD[i])^2/tot.var),
           lty=1,col="darkgrey")
    lines(c(i-0.25,i+0.25),100*rep((D[i]+dD[i])^2/tot.var,2),
          lwd=1,col="darkgrey")
    lines(c(i-0.25,i+0.25),100*rep((D[i]-dD[i])^2/tot.var,2),
          lwd=1,col="darkgrey")
  }
  points(var.eof,cex=1.5)
  points(var.eof,pch=20,cex=1.2,col="darkgrey")
  points(ip,var.eof[ip],pch=20,col="red",cex=1.2)
  attr(var.eof,'errorbar') <- cbind(100*(D-dD)^2/tot.var,100*(D+dD)^2/tot.var)
  invisible(var.eof)
}



plot.field <- function(x,is=NULL,it=NULL,FUN="mean",map.type='rectangle',verbose=FALSE,...) {
  if (verbose) print("plot.field")
  stopifnot(!missing(x),inherits(x,'field'))

  d <- dim(x)
  if (d[2]==1) {
    if (verbose) print('one grid point')
    plot.station(x,verbose=verbose,...)
    return()
  }

  # To distinguish between tpyes of plots - movmuller or are mean
  twocoord <- (length(is)==2)
  if (twocoord) {
    l1 <- length(is[[1]]); l2 <- length(is[[2]])
  } else {
    l1 <- 0; l2 <- 0
  }
  #print(c(l1,l2))
  if (inherits(is,'station')) {
    # If a station object is provided - extract a time series for a
    # the corresponding location
      lon <- attr(is,'longitude')
      lat <- attr(is,'latitude')
  } else if ((inherits(x,'list')) & twocoord & (l1==1) & (l2==1) ) {
    # Is spatial coordinates  
      lon <- is[[1]]
      lat <- is[[2]]
      z <- NULL
  } else if (!is.null(is)) {
      nms <- names(is)
      lon <- attr(x,'longitude')
      lat <- attr(x,'latitude')
                                        #print(nms)
      if (length(nms)==2) {
          lon <- is[[1]]
          lat <- is[[2]]
          y <- subset(x,is=is,it=it)
          z <- aggregate.area(y,FUN=FUN)
      } else if ( (length(nms)==1) & (tolower(nms)=="lon") ) {
      # Hovmuller diagram along latitude
      #print(is[[1]]); print(lon); print(d)
          if (length(is[[1]])== 1) {
              picklon <- lon[max( (1:length(lon))[lon <= is[[1]]] )]
        #print(picklon)
              xy <- rep(lon,length(lat))
              yx <- sort(rep(lat,length(lon)))
              ix <- is.element(xy,picklon)
              z <- x[,ix]
              z <- attrcp(x,z)
        #image(z)
        attr(z,'longitude') <- picklon
        attr(z,'latitude') <- attr(x,'latitude')
        #print(attr(x,'dimensions'))
        attr(z,'dimensions') <- c(1,attr(x,'dimensions')[2:3])
        class(z) <- c('xsection',class(x)[-1])
      } else if (length(is[[1]])== 2) {
        y <- subset(x,is=list(lon=is[[1]],lat=range(lat)))
        xy <- rep(attr(x,'longitude'),d[2])
        yx <- sort(rep(attr(x,'latitude'),d[1]))
        X <- as.data.frame(coredata(x)); colnames(yx)
        Z <- aggregate(x,by=yx,FUN=FUN)
        z <- zoo(Z,order.by=index(x))
        z <- attrcp(x,z)
        class(z) <- class(x)
      }
    } else if ( (length(nms)==1) & (tolower(nms)=="lat") ) {
      # Hovmuller diagram along longitude
      if (length(is[[1]])== 1) {
        picklat <- lat[max( (1:length(lat))[lat <= is[[1]]] )]
        xy <- rep(attr(x,'longitude'),length(lat))
        yx <- sort(rep(attr(x,'latitude'),length(lon)))
        iy <- is.element(yx,picklat)
        z <- x[,iy]
        z <- attrcp(x,z)
        attr(z,'latitude') <- picklat
        attr(z,'longitude') <- attr(x,'longitude')
        attr(z,'dimensions') <- c(attr(x,'dimensions')[1],1,attr(x,'dimensions')[3])
        class(z) <- c('xsection',class(x)[-1])
      } else if (length(is[[1]])== 2) {
          y <- subset(x,is=list(lon=range(lon),lat=is[[1]]))
        xy <- rep(attr(x,'longitude'),d[2])
        yx <- sort(rep(attr(x,'latitude'),d[1]))
        X <- as.data.frame(coredata(x)); colnames(xy)
        Z <- aggregate(x,by=xy,FUN=FUN)
        z <- zoo(Z,order.by=index(x))
        class(z) <- class(x)
      }
    }
    
  } else {lon <- NULL; lat <- NULL}

  if ( is.null(lon) & is.null(lat) ) {
    #print("aggregate")
    z <- aggregate.area(x,is=is,FUN=FUN)
    class(z) <- c('station',class(x)[-1])
  } else if ( (length(lon)==1) & (length(lat)==1) ) {
    #print("regrid")
    z <- regrid(x,list(x=lon,y=lat))
    class(z) <- c('station',class(x)[-1])
  }
  #print("plot")
  plot(z,map.type=map.type,...)
  z <- attrcp(x,z,ignore=c("longitude","latitude"))
  attr(z,'history') <- history.stamp(x)
  if (inherits(x,'station')) lines(y,col="red",lwd=2)
  invisible(z)
}

plot.pca <- function(y,cex=1,verbose=FALSE,new=TRUE,...) {
  if (verbose) print('plot.pca')
  if(inherits(y,"trajectory")) {
    plot.pca.trajectory(y,cex=cex,new=new,verbose=verbose,...)
  } else {
    attr(y,'longname') <- attr(y,'longname')[1]
    plot.eof.field(y,verbose=verbose,new=new,cex=cex,...)
  }
}


plot.ds.pca <- function(x,ip=1,verbose=FALSE,
                        colbar1=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                            type="p",cex=1,show=TRUE,
                            h=0.6, v=1,pos=0.05),colbar2=NULL,...) {
  y <- x # quick fix

  if (verbose) print('plot.ds.pca')
  if (is.null(colbar2)) colbar2 <- colbar1
  attr(y,'longname') <- attr(y,'longname')[1]
  #par(fig=c(0,0.45,0.5,0.975),new=TRUE)
  par(fig=c(0,0.5,0.5,0.975)) #par(fig=c(0,0.45,0.5,0.975))

  if (verbose) print('PCA ip')
  map.pca(y,ip=ip,verbose=verbose,new=FALSE,colbar=colbar1,fig=c(0,0.5,0.5,0.975),...)

  title(paste("PCA Pattern # ",ip,sep=""))
  par(fig=c(0.55,0.975,0.5,0.975),new=TRUE)
  if (verbose) print('Predictor pattern')
  map(attr(y,'predictor.pattern'),ip=ip,new=FALSE,
      colbar=colbar2,verbose=verbose,
      main=paste("EOF Pattern # ",ip,sep=""))
  #title(paste("EOF Ip # ",ip,sep=""))
  if (!is.null(attr(y,'evaluation'))) {
    if (verbose) print('Evaluation results')
    par(fig=c(0.05,0.45,0.05,0.475),new=TRUE)
    ## Get the right pattern
    xvp <- (ip-1)*2 +1
    xok <- is.finite(attr(y,'evaluation')[,xvp]) & is.finite(attr(y,'evaluation')[,xvp+1])
    xcor <- cor(attr(y,'evaluation')[xok,xvp],attr(y,'evaluation')[xok,xvp+1])
    plot(attr(y,'evaluation')[,xvp],attr(y,'evaluation')[,xvp+1],
         main=paste('Cross-validation: r=',round(xcor,2)),
         xlab='original data',ylab='prediction',pch=19,col="grey")
    lines(range(c(attr(y,'evaluation')),na.rm=TRUE),
          range(c(attr(y,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(attr(y,'evaluation')[,1]),
                      x=coredata(attr(y,'evaluation')[,2]))
    xvalfit <- lm(y ~ x, data = cal)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    par(fig=c(0.55,0.975,0.05,0.475),new=TRUE)
    xlim <- range(index(attr(y,'original_data')),index(y))
    ylim <- range(attr(y,'original_data')[,ip],y[,ip],na.rm=TRUE)
    y0 <- attr(y,'original_data')
    plot(y0[,ip],lwd=2,type='b',pch=19,xlim=xlim,ylim=ylim)
    if ( (class(index(y))=='Date') & (class(index(y0))=='numeric') & inherits(y,'annual') ) 
      index(y) <- year(index(y))
    if ( (class(index(y))=='numeric') & (class(index(y0))=='Date') & inherits(y,'annual') ) 
      index(y) <- as.Date(paste(index(y),'01-01',sep='-'))
    
    lines(zoo(y[,ip]),lwd=2,col='red',type='b')
    legend(x=index(attr(y,'original_data')[,ip])[1],
           y=max(attr(y,'original_data')[,ip],na.rm=TRUE)+
               diff(range(attr(y,'original_data')[,ip]))/10,
           legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
           lwd=c(2,2),pch=c(21,19),bty="n")
  } else {
    par(fig=c(0.05,0.975,0.05,0.475),new=TRUE)
    plot(attr(y,'original_data')[,ip],lwd=2,type='b',pch=19)
    lines(zoo(y[,ip]),lwd=2,col='red',type='b')
    xvalfit <- NULL
  }
}

plot.ds.eof <- function(x,ip=1,
                        colbar1=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,type="p",cex=2,show=TRUE,
                        h=0.6, v=1,pos=0.05),colbar2=NULL,verbose=FALSE,...) {
  y <- x # quick fix
  if (verbose) print('plot.ds.eof')
  if (is.null(colbar2)) colbar2 <- colbar1 
  attr(y,'longname') <- attr(y,'longname')[1]
  par(fig=c(0,0.5,0.5,1),mar=c(3,5,4.2,1),mgp=c(3,0.5,0.5))
  map.eof(y,ip=ip,verbose=verbose,new=FALSE,colbar=colbar1,
          main=paste("Predictand EOF pattern # ",ip,sep=""),...)
  par(fig=c(0.5,1,0.5,1),mar=c(3,4,4.2,1),new=TRUE)
  map(attr(y,'predictor.pattern'),ip=ip,new=FALSE,
      colbar=colbar2,verbose=verbose,
      main=paste("Predictor EOF pattern # ",ip,sep=""))
  #title(paste("EOF Pattern # ",ip,sep=""))
  if (!is.null(attr(y,'evaluation'))) {
    par(fig=c(0,0.5,0,0.48),mar=c(3,4.5,3,1),new=TRUE)
    pc.obs <- attr(y,'evaluation')[,1+2*(ip-1)]
    pc.ds <- attr(y,'evaluation')[,2+2*(ip-1)]
    plot(pc.obs,pc.ds,main='Cross-validation',xlab='original data',
         ylab='prediction',pch=19,col="grey",
         xlim=range(pc.obs,pc.ds),ylim=range(pc.obs,pc.ds))
    lines(range(c(attr(y,'evaluation')),na.rm=TRUE),
          range(c(attr(y,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(pc.obs),x=coredata(pc.ds))
    xvalfit <- lm(y ~ x, data = cal)
    r.xval <- round(cor(pc.obs,pc.ds),2)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    text(min(pc.obs)+diff(range(pc.obs))/12,max(pc.ds),
         paste("r =",r.xval),#round(xvalfit$coefficients[2],digits=2)),
         pos=4,cex=0.9)
    par(fig=c(0.5,1,0,0.48),mar=c(3,4.5,3,1),new=TRUE)
    xlim <- range(index(attr(y,'original_data')),index(y))
    ylim <- range(attr(y,'original_data')[,ip],y[,ip],na.rm=TRUE)
    y0 <- attr(y,'original_data')
    plot(y0[,ip],
         main=paste("PC",ip),ylab="",
         ylim=ylim*c(1.2,1.2),xlim=xlim,
         lwd=2,type='b',pch=19)
    if ( (class(index(y))=='Date') & (class(index(y0))=='numeric') & inherits(y,'annual') ) 
      index(y) <- year(index(y))
    if ( (class(index(y))=='numeric') & (class(index(y0))=='Date') & inherits(y,'annual') ) 
      index(y) <- as.Date(paste(index(y),'01-01',sep='-'))
    lines(zoo(y[,ip]),lwd=2,col='red',type='b')
    legend(x=index(attr(y,'original_data')[,ip])[1],
           y=max(attr(y,'original_data')[,ip],na.rm=TRUE)+
               diff(range(attr(y,'original_data')[,ip]))/3,
           legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
           lwd=c(2,2),pch=c(21,19),bty="n")
  } else {
    par(fig=c(0,1,0,0.48),mar=c(3,4.5,3,1),new=TRUE)
    plot(attr(y,'original_data')[,ip],
         main="PC1",ylab="",
         ylim=range(attr(y,'original_data')[,ip])*c(1.6,1.6),
         lwd=2,type='b',pch=19)
    lines(zoo(y[,ip]),lwd=2,col='red',type='b')
    legend(x=index(attr(y,'original_data')[,ip])[1],
           y=max(attr(y,'original_data')[,ip],na.rm=TRUE)+
               diff(range(attr(y,'original_data')[,ip]))/3,
           legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
           lwd=c(2,2),pch=c(21,19),bty="n")
    xvalfit <- NULL
  }  
}

vis.pca <- function(x,cex=1.5,new=TRUE) {

  y <- x # quick fix
  col <- colscal(col=varid(y)); nc <- length(col)
  #if (is.precip(y)) col <- rev(col)
  lon <- attr(y,'longitude') 
  lat <- attr(y,'latitude') 
  N <- length(lon)
  R2 <- round(100*attr(y,'eigenvalues')^2/attr(y,'tot.var'),2)

  #print(N); print(length(attr(y,'mean')))
  m <- min(3,dim(attr(y,'pattern'))[2])
  
  # Set scale for colour scheme
  #str(y)
  a.T <- matrix(rep(NA,4*N),4,N)
  ax <- quantile(abs(attr(y,'mean')),0.99,na.rm=TRUE)
  if (min(attr(y,'mean'))<0) scale0 <- seq(-ax,ax,length=nc) else
                             scale0 <- seq(0,ax,length=nc)
  ax <- quantile(abs(attr(y,'pattern')),0.99,na.rm=TRUE)
  scale <- seq(-ax,ax,length=nc)

  #print("here")
  for (i in 1:N) {
    a.T[1,i] <-  sum(attr(y,'mean')[i] > scale0)
    for (j in 1:m) 
      a.T[j+1,i] <-  sum(attr(y,'pattern')[i,j] > scale)
  }
  a.T[a.T < 1] <- 1; a.T[a.T > 100] <- 100

  if (new) dev.new(width=5,height=7)
  par(mfrow=c(3,2),mar=c(3.5,3,3.5,3),bty="n",xaxt="n",yaxt="n")
  
  plot(lon,lat,
       main="Climatology",
       col=col[a.T[1,]],pch=19,xlab="",ylab="",cex=cex)
  points(lon,lat,cex=cex)
  data("geoborders",envir=environment())
  lines(geoborders,col='grey40')
  lines(geoborders$x - 360,geoborders$y,col='grey40')
  points(lon,lat,cex=cex,col=col[a.T[1,]],pch=19)

  plot(lon,lat,
       main=paste("EOF #1:",R2[1],"% of variance"),
       col=col[a.T[2,]],pch=19,xlab="",ylab="",cex=cex)
  points(lon,lat,cex=cex)
  lines(geoborders)
  lines(geoborders$x - 360,geoborders$y)
  points(lon,lat,cex=cex,col=col[a.T[2,]],pch=19)

  plot(lon,lat,
       main=paste("EOF #2:",R2[2],"% of variance"),
       col=col[a.T[3,]],pch=19,xlab="",ylab="",cex=cex)
  points(lon,lat,cex=cex)
  lines(geoborders,col='grey40')
  lines(geoborders$x - 360,geoborders$y,col='grey40')
  points(lon,lat,cex=cex,col=col[a.T[3,]],pch=19)

  plot(lon,lat,
       main=paste("EOF #3:",R2[3],"% of variance"),
       col=col[a.T[4,]],pch=19,xlab="",ylab="",cex=cex)
  points(lon,lat,cex=cex)
  lines(geoborders,col='grey40')
  lines(geoborders$x - 360,geoborders$y,col='grey40')
  points(lon,lat,cex=cex,col=col[a.T[4,]],pch=19)

  par(mar=c(1,0,0,0),fig=c(0.1,0.3,0.665,0.695),new=TRUE,cex.axis=0.6)
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale0)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),labels=nl)

  par(mar=c(1,0,0,0),fig=c(0.1,0.3,0.32,0.35),new=TRUE,cex.axis=0.6,xaxt="n")
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),labels=nl)

  par(mar=c(1,0,0,0),fig=c(0.6,0.8,0.665,0.695),new=TRUE,cex.axis=0.6,xaxt="n")
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),labels=nl)

  par(mar=c(1,0,0,0),fig=c(0.6,0.8,0.32,0.35),new=TRUE,cex.axis=0.6,xaxt="n")
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),labels=nl)
  
  par(mfcol=c(1,1),fig=c(0,1,0,0.33),new=TRUE,xaxt="s",yaxt="n",bty="n",
      mar=c(2,2,1,1))
  ylim <- 2*range(coredata(y[,1:m]),na.rm=TRUE)
  plot(y[,1]+0.5*ylim[2],lwd=2,ylim=ylim)
  grid()
  col <- c("red","blue")
  for (j in 1:m) lines(y[,j+1]+(1-j)*0.5*ylim[2],lwd=2,col=col[j])
  legend(index(y)[1],ylim[1],c('PC 1','PC 2','PC 3'),
         col=c('black','red','blue'),bty='n',lwd=2)
  invisible(a.T)
}

plot.mvr <- function(x) {
  plot(x$fitted.values)
}


plot.cca <- function(x,icca=1,
                     colbar1=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,type="p",cex=2,show=TRUE,
                        h=0.6, v=1,pos=0.05),colbar2=NULL,verbose=FALSE,new=TRUE,...) {
  if (verbose) print("plot.cca")
  if (new) dev.new()
  par(mfrow=c(2,2),bty="n",xaxt="n",yaxt="n")
  if (is.null(colbar2)) colbar2 <- colbar1
  map.cca(x,icca=icca,colbar1=colbar1,colbar2=colbar2,verbose=verbose,...)

  w.m <- zoo((x$w.m[,icca]-mean(x$w.m[,icca],na.rm=TRUE))/
             sd(x$w.m[,icca],na.rm=TRUE),order.by=x$index)
  v.m <- zoo((x$v.m[,icca]-mean(x$v.m[,icca],na.rm=TRUE))/
             sd(x$v.m[,icca],na.rm=TRUE),order.by=x$index)
  r <- cor(x$w.m[,icca],x$v.m[,icca])
  par(bty="n",xaxt="s",yaxt="s",xpd=FALSE,mar=c(2,1.5,1.5,0.5),
      fig=c(0.02,1,0.1,0.45),new=TRUE,cex.axis=0.8,cex.lab=0.8)
  plot(w.m,col="blue",lwd=2,
       main=paste("CCA pattern ",icca," for ",varid(x),
         "; r= ",round(r,2),sep=""),xlab="",ylab="")
  lines(v.m,col="red",lwd=2)

  par(fig=c(0,1,0,0.1),new=TRUE, xaxt="n",yaxt="n",bty="n",
      mar=c(0,0,0,0))
  plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
  legend(0.01,0.90,c(paste(attr(x$X,'source')[1],attr(x$X,'variable')[1]),
                     paste(attr(x$Y,'source')[1],attr(x$Y,'variable')[1])),
         col=c("red","blue"),lwd=2,lty=1,
         bty="n",cex=0.5,ncol=2,text.col="grey40")
  
  ##par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
  ##    fig=c(0,1,0.1,1),new=TRUE)
  ## par(fig=c(0,1,0,0.1),new=TRUE, mar=c(0,0,0,0))
}


plot.list <- function(x,...) {
  plot(combine.ds(x),...)
}


plot.diagnose <- function(x,...) {
  if ( (inherits(x,"eof")) & (inherits(x,"comb")) ) plot.diagnose.comb.eof(x,...) else
  if (inherits(x,"dsensembles")) plot.diagnose.dsensemble(x,...)
}

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


plot.diagnose.dsensemble <- function(x,new=TRUE,mgp=c(2,1,0),cex=NULL,map.show=TRUE,
                                     map.type=NULL,verbose=FALSE,main=NULL,...) {
  if (verbose) print('plot.diagnose.dsensemble')

  if (is.null(map.type)) {
    if( inherits(x,"field") | length(lon(x))!=length(lat(x)) |
        (length(lon(x))==2 & length(lat(x))==2) ) {
      map.type <- "rectangle"
    } else {
      map.type <- "points"
    }
  }
  
  Y <- -round(200*(0.5-pbinom(x$outside,size=x$N,prob=0.1)),2)
  X <- -round(200*(0.5-pnorm(x$deltaobs,mean=mean(x$deltagcm),
                             sd=sd(x$deltagcm))),2)
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

nam2expr <- function(x) {
  y <- x
  for (i in 1:length(y)) {
    z <- switch(tolower(x[i]),
                  't2m'=expression(T[2 * m]),
                  'tmax'=expression(T[x]),
                  'tmin'=expression(T[n]),
                  'tas'=expression(T[2 * m]))
    if (is.null(z)) y[i] <- x[i] else y[i] <- z
  }
  return(y)
}

 
plot.xval <- function(x,new=TRUE,...) {
  if (new) dev.new()
  par(bty="n")
  unit <- attr(x,'unit')
  cols <- rgb(seq(0,1,length=20),rep(0,20),rep(0,20))
  cexs <- seq(1.5,0.5,length=20)^2
  eofindex <- as.integer(gsub('X.','',
                         names(attr(x,'original_model')$coefficients[-1])))
  print(eofindex)
  if (unit=='deg C') unit <- expression(degree*C)
  
  class(x) <- "zoo"
  plot(x[,1],type="b",pch=19,lwd=2,
       main=paste("Cross-validation:",attr(x,'location'),attr(x,'variable')),
       xlab="",ylab= unit,
       sub=paste("r=",attr(x,'correlation'),"rmse=", attr(x,'rmse'),
         unit))
  lines(x[,2],lwd=2,col="red")
  lines(attr(x,'fitted_values_all'),col="red",lty=3,pch="x")
  
  ## par(new=TRUE,fig=c(0.1,1,0.1,0.5),bty="n")
  ## plot(c(0,1),c(0,1),col="white",xlab="",ylab="",axes=F)
  legend(rep(range(index(x))[1],2),rep(range(x)[2],2),c("obs","x-valid","fit to all"), col=c("black","red","red"),lwd=c(2,2,1),lty=c(1,1,3),bty="n")
  
  dev.new()
  par(bty="n")
  boxplot(t(attr(x,'beta')[-1,]),col="grey90",border="grey40",
          xlim=range(eofindex),main="Model coefficients",
          ylab=expression(beta),xlab="EOF number")
  for (i in 1:20) 
    points(eofindex,attr(x,'original_model')$coefficients[-1],
           pch=19,col=cols[i],cex=cexs[i])
  grid()
}

plot.dsensemble.pca <- function(x,pts=FALSE,target.show=TRUE,map.show=TRUE,it=0,ip=1,
                               envcol=rgb(1,0,0,0.2),legend.show=TRUE,verbose=FALSE,
                               ...) {
  if (verbose) print("plot.dsensemble.pca")
  stopifnot(inherits(x,'dsensemble') & inherits(x,'pca'))
  d <- index(x[[3]])
  pc <- x[3:length(x)]
  pc <- array(unlist(pc), dim = c(dim(pc[[1]]), length(pc)))
  pc <- lapply(seq(dim(pc)[2]), function(x) pc[ , x, ])
  fn <- function(x) {
    x <- zoo(x,order.by=d)
    class(x) <- c("dsensemble","station","zoo")
    invisible(x)
  }
  pc <- lapply(pc,fn)
  for (i in 1:length(pc)) {
    attr(pc[[i]],"station") <- as.station(x[[2]][,i],param=attr(x,"variable"),
                              longname=attr(x,"longname"),unit=attr(x,"unit"))
  }
  plot(pc[[ip]],ylab=paste("PC",ip,sep=""))
}

plot.dsensemble <- function(x,verbose=FALSE,plot = TRUE, ...) {
  if(verbose) print("plot.dsensemble")
  if (inherits(x,c('pca','eof'))) {
    y <- plot.dsensemble.multi(x,verbose=verbose, plot = plot, ...) 
  } else if (inherits(x,'zoo')) {
    y <- plot.dsensemble.one(x,verbose=verbose,...) 
  } else if (inherits(x,'station')) {
    x <- as.station(x,verbose=verbose)
    y <- plot(x,verbose=verbose,...)
  } else {
    print(paste('Unknown class - do not know how to plot',
                paste(class(x),collapse=", ")))
    y <- x
  }
  if(verbose) print("exit plot.dsensemble")
  invisible(y)
}

## Plot multiple stations/spatially aggregated field.
plot.dsensemble.multi <- function(x,it=c(2000,2099),FUNX='mean',verbose=FALSE,
                                  anomaly=FALSE,test=FALSE, plot = TRUE, ...) {
  if (verbose) print('plot.dsensemble.multi')
  
  if (inherits(x,c('pca','eof'))) {
    Y <- expandpca(x,it=it,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
    if (plot) plot(Y,verbose=verbose,...)
    #invisible(Y)
  } else {
    #return(NULL)
    Y <- NULL
  }
  if (verbose) print('exit plot.dsensemble.multi')
  invisible(Y)
}

## Plots one time series
plot.dsensemble.one <-  function(x,pts=FALSE,it=0,
                             envcol=rgb(1,0,0,0.2),legend.show=TRUE,ylab=NULL,
                             obs.show=TRUE,target.show=TRUE,map.show=TRUE,map.type=NULL,map.insert=TRUE,
                             usegooglemap=TRUE,new=FALSE,xrange=NULL,yrange=NULL,
                             alpha=0.5,alpha.map=0.7,mar=c(5.1,4.5,4.1,2.1),
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
  
  if (verbose) {print("diagnose"); class(z)}
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
       ylab=ylab,map.show=FALSE,new=new)
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
  
  lines(zoo(mu,order.by=year(z)),col=rgb(1,0.7,0.7),lwd=3)
  lines(zoo(q05,order.by=year(z)),col=rgb(0.5,0.5,0.5),lty=2)  
  lines(zoo(q95,order.by=year(z)),col=rgb(0.5,0.5,0.5),lty=2)
  if (obs.show) lines(y,type="b",pch=19)

  if (!is.null(diag)) {
    index(diag$y) <- year(diag$y)
    outside <- diag$above | diag$below
    points(zoo(coredata(diag$y)[which(outside)],
             order.by=year(diag$y)[which(outside)]),col="grey")
  }

  title(main=toupper(loc(x)),cex.main=1)
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
    legend(0.05,0.90,c(paste("Past trend:",round(diag$deltaobs,2)),
                      paste(diag$robs,'%'),
                      paste(diag$outside,"observations"),
                      "p-value: "),
            bty="n",cex=0.7,text.col="grey40")
    #legend(0.5,0.90,c(expression(paste(levels(factor(unit(x)))[1]/d*e*c*a*d*e)),
    #                  "ensemble trends > obs.",
    #                  "outside ensemble 90% conf.int.",
    #                  paste(round(100*pbinom(diag$outside,size=diag$N,prob=0.1)),"%")),
    #        bty="n",cex=0.7,text.col="grey40")
    legend(0.5,0.90,c(paste(levels(factor(unit(y)))[1],"/decade",sep=""),
                      "ensemble trends > obs.",
                      "outside ensemble 90% conf.int.",
                      paste(round(100*pbinom(diag$outside,size=diag$N,prob=0.1)),"%")),
            bty="n",cex=0.7,text.col="grey40")    
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

plot.xsection <- function(x,...) {
                                        #print("plot.xsection")
    d <- attr(x,'dimensions')
                                        #print(d)
    X <- coredata(x)
    
    if (d[1]==1) {
        attr(X,'longitude') <- index(x)
        attr(X,'latitude') <- attr(x,'latitude')
        attr(X,'dimensions') <- attr(x,'dimensions')[c(3,2)]
        
    } else {
        attr(X,'longitude') <- attr(x,'longitude')
        attr(X,'latitude') <- index(x)
        attr(X,'dimensions') <- attr(x,'dimensions')[c(1,3)]
        X <- t(X)
    }
    attr(X,'variable') <- attr(x,'variable')
    
    attr(X,'unit') <- attr(x,'unit')
    attr(X,'source') <- attr(x,'source')

                                        # print(dim(X)); print(c(length(lon(X)),length(lat(X))))
    lonlatprojection(x=X,what="fill",geography=FALSE,...)
}


# plot wet/dry or cold/warm spells.
plot.spell <- function(x,xlim=NULL,ylim=NULL) {
  bar <- function(x,col) {
    rect(x[1],x[2],x[3],x[4],col=col,border=col)
  }
  t <- index(x)
  h <- coredata(x[,1]); l <- coredata(x[,2])
  ih <- is.finite(h); il <- is.finite(l)
  h <- h[ih]; th1 <- t[ih]
  l <- l[il]; tl1 <- t[il]
  th2 <- th1 + h; tl2 <- tl1 + l
  par(bty="n")

  col <- c('red','blue')
  runs <- c('hot','cold')
  spelltype <- 'hot and cold'
  if (sum(is.element(attr(x,'variable'),c('wet','dry')))>0) {
    col <- c('darkblue','brown')
    runs <- c('wet','dry')
    spelltype <- 'wet and dry' 
  }

  tunit <- attr(x,'threshold.unit')[1]
  for (i in 1:length(tunit)) {
    if ( (is.na(tunit[i]) | is.null(tunit[i])) ) tunit[i] <- " "
    if ((tunit[i]=='degree Celsius') | (tunit[i]=='deg C') | (tunit[i]=='degC'))
         tunit[i] <- 'degree*C'
  }


  plot(range(t),c(-1,1)*max(c(h,l),na.rm=TRUE),type="n",
       xlab="",ylab="Spell length",xlim=xlim,ylim=ylim,
       main=paste(attr(x,'location')[1],": ",spelltype[1],sep=""))
  leg <- eval(parse(text=paste("expression(paste(X > ",
                      attr(x,'threshold'),"*",tunit,"))")))
  text(t[1],0.75*max(c(h,l),na.rm=TRUE),leg,srt=90,cex=0.7,col="grey")
  leg <- eval(parse(text=paste("expression(paste(X <= ",
                      attr(x,'threshold'),"*",tunit,"))")))
  text(t[1],-0.75*max(c(h,l),na.rm=TRUE),leg,srt=90,cex=0.7,col="grey")
  lines(range(t),rep(0,2))
  apply(cbind(th1,rep(0,length(h)),th2,h),1,bar,col[1])
  apply(cbind(tl1,rep(0,length(l)),tl2,-l),1,bar,col[2])
  
}

plot.ssa <- function(ssa,main="SSA analysis",sub="")  {
    if ( (class(ssa)[1]!="SSA") ) stop("Need an 'SSA' object")
    nt <- ssa$nt
    dev.new()
    plot(ssa$d,main=main,sub=sub,ylab="Singular value",pch=20,col="grey50")
    points(ssa$d)
    grid()

    dev.new()
    par(mfcol=c(3,1))
    plot(ssa$v[,1],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA vector: mode 1",lwd=3,col="grey70")
    grid()
    plot(ssa$v[,2],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA vector: mode 2",lwd=3,col="grey70")
    grid()
    plot(ssa$v[,3],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA vector: mode 3",lwd=3,col="grey70")
    grid()


    dev.new()
    par(mfcol=c(3,1))
    if (class(ssa)[3] == "monthly.station.record") {
      yy <- sort(rep(ssa$x$yy,12)); yy <- yy[1:ssa$Nm]
      mm <- rep(1:12,nt); mm <- mm[1:ssa$Nm]
      #print(rbind(yy,mm))
      #print(dim(ssa$v)); print(dim(ssa$u)); print(length(yy));
      #print(length(mm));  print(ssa$Nm); print(nt); print(length(ssa$x$yy))
      plot(yy + (mm-0.5)/12, ssa$u[,1],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(yy + (mm-0.5)/12, ssa$u[,2],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(yy + (mm-0.5)/12, ssa$u[,3],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
    } else if (class(ssa)[3] == "daily.station.record") {
      plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
           ssa$u[,1],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
           ssa$u[,2],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
           ssa$u[,3],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
    } else {
      plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
           ssa$u[,1],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
           ssa$u[,2],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
      plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
           ssa$u[,3],type="l",main=main,sub=sub,
           xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
      grid()
    }
  }

plot.nevents <- function(x,verbose=FALSE,main=NULL,xlab=NULL,ylab=NULL,col=NULL,...) {
  # Plot the results from 
  if (verbose) print('plot.nevents')
  par(bty='n')
  if (is.null(main)) main <- loc(x)
  if (is.null(xlab)) xlab <- ""
  if (is.null(ylab)) ylab <- attr(x,'info')
  if (is.null(col)) {
    if (is.T(attr(x,'observation')))
      col <- c(rgb(0.5,0.5,0.7,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.8,0.5,0.8,0.5),
               rgb(0.3,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.6,0.3,0.6,0.5)) else
      col <- c(rgb(0.3,0.3,0.6,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.5,0.5,0.7,0.5),
               rgb(0.6,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.8,0.5,0.8,0.5))
  }
  plot.zoo(x,plot.type='single',lwd=5,main=main,
       xlab=xlab,ylab=ylab,col=col,...)
  grid()
  points(attr(x,'observation'),pch=19)
  lines(attr(x,'nwd.pre'),col=rgb(0.5,0.5,0.5,0.5))
}

barplot.station <- function(x,threshold=0,...) {
    stopifnot(inherits(x,'station'))
    x.above <- x.below <- x
    x.above[x < threshold] <- NA
    x.below[x > threshold] <- NA
    ylim <- range(pretty(coredata(x)),na.rm=TRUE)
    barplot(as.numeric(x),col='white',ylim=ylim,border=NA)
    barplot(as.numeric(x.above),col='red',names.arg=year(x),
            ylab=paste(varid(x),'[',unit(x),']'),axes=FALSE,
            border=NA,add=TRUE)
    barplot(as.numeric(x.below),col='blue',axes=FALSE,border=NA,add=TRUE)
    title(toupper(loc(x)))
}
