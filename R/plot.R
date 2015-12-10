##plot <- function(x,y, ...)  UseMethod("plot")
     
plot.station <- function(x,plot.type="single",new=TRUE,
                         lwd=3,type='l',pch=0,main=NULL,col=NULL,
                         xlim=NULL,ylim=NULL,xlab="",ylab=NULL,
                         errorbar=TRUE,legend.show=FALSE,
                         map.show=TRUE,map.type="points",
                         mar=c(4.5,4.5,0.75,0.5),
                         alpha=0.3,verbose=FALSE,...) {

  if (verbose) print('plot.station')

  if (!is.numeric(lon(x)) | !is.numeric(lat(x))) {
    map.show <- FALSE
  } else if (length(lon(x))!=length(lat(x))) {
    map.type <- "rectangle"
  }
  
  fig <- c(0,1,0,0.95)
  if (map.show) fig[4] <- 0.8
  if (legend.show) fig[3] <- 0.05  
  par(bty="n",xaxt="s",yaxt="s",xpd=FALSE,cex.axis=1,fig=fig,mar=mar)
  ## browser()
  ## if (is.null(ylim))
  ##     if (is.null(dim(x)))
  ##         ylim <- pretty(x)
  ##     else              
  ##         ylim <- apply(x,2,pretty,n=5)
  if (is.null(xlim))
    xlim <- range(index(x))
  if (is.null(ylim))
    ylim <- pretty(as.numeric(x))
  
  unit <- attr(x,'unit')[1]
  for (i in 1:length(unit)) {
    if ( (is.na(unit[i]) | is.null(unit[i])) ) unit[i] <- " "
    if ((unit[i]=='degree Celsius') | (unit[i]=='deg C') | (unit[i]=='degC'))
         unit[i] <- 'degree*C'
  }
  
  if (plot.type=="single") {
      if (is.null(ylab))
          ylab <- try(eval(parse(text=paste("ylab <- expression(",varid(x),
                                     "*phantom(0)*(",unit,"))"))),silent=TRUE)
      if (inherits(ylab,"try-error")) ylab <- unit(x)
  }
  else if (is.null(ylab) & (length(levels(factor(stid(x))))>1))
      ylab <- stid(x)
  else if (is.null(ylab))
      ylab <- apply(x,1,varid)
  
  if (is.null(main)) main <- attr(x,'longname')[1] 
  if (is.null(col)) {
    if (is.null(dim(x))) {
      col <- "blue"
    } else if (!is.null(lon(x)) & !is.null(lat(x)) &
               length(lon(x))==dim(x)[2] &
               length(lat(x))==dim(x)[2]) {
      nx <- (lon(x)-min(lon(x)))/diff(range(lon(x)))
      ny <- (lat(x)-min(lat(x)))/diff(range(lat(x)))
      col <- rgb(1-ny,nx,ny,alpha)
    } else {
      col <- adjustcolor(rainbow(length(x[1,])),alpha=alpha)
    }
  } else {
    col <- adjustcolor(col,alpha.f=alpha)
  }
  
  ns <- length(stid(x))
#  if ( (ns > 1) & (plot.type=="multiple") ) {
#    for (i in 1:ns) {
#        z <- try(eval(parse(text=paste("ylab[",i,"] <- expression(",ylab[i],
#                      "*phantom(0)*(",unit[i],"))"))),silent=TRUE)
#        if (inherits(z,"try-error")) ylab[i] <- unit[i]
#      }
#  }

  errorbar <- errorbar & !is.null(err(x))
  
  #print(ylab)
  class(x) <- "zoo"
  plot.zoo(x,plot.type=plot.type,xlab=xlab,ylab=ylab,
           col=col,xlim=xlim,ylim=ylim,lwd=lwd,type=type,pch=pch,...)
  mtext(main,side=3,line=1,adj=0,cex=1.1)
  par0 <- par()
  
  if (plot.type=="single") {
    if (errorbar) {
      # REB 2014-10-03: add an errorbar to the plots.
      segments(index(x),x-err(x),index(x),x+err(x),
               lwd=3,col=rgb(0.5,0.5,0.5,0.25))
#      d.err <- dim(std.err)
#      dt <- 0.3*diff(index(x))[1]
#      if (is.null(d.err)) d.err <- c(length(std.err),1)
#      for (i in 1:d.err[2]) {
#        for (j in 1:d.err[splot.dse1])
#          lines(rep(index(x)[j],2),rep(x[j],2) + std.err[j]*c(-1,1),
#                lwd=3,col=rgb(1,0.5,0.5,0.25))
#          lines(rep(index(x)[j],2) + dt*c(-1,1),rep(x[j],2) + std.err[j],
#                lwd=1,col=rgb(1,0.5,0.5,0.25))
#          lines(rep(index(x)[j],2) + dt*c(-1,1),rep(x[j],2) - std.err[j],
#                lwd=1,col=rgb(1,0.5,0.5,0.25))
#       }
    }
    
    par(fig=c(0,1,0,0.1),new=TRUE, mar=c(0,0,0,0),xaxt="s",yaxt="s",bty="n")
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    
    if(legend.show) {
      legend(0.01,0.95,paste(attr(x,'location'),": ",
                           #attr(x,'aspect'),
                           #attr(x,'longname')," - ",
                           round(attr(x,'longitude'),2),"E/",
                           round(attr(x,'latitude'),2),"N (",
                           attr(x,'altitude')," masl)",sep=""),
           bty="n",cex=0.6,ncol=3,text.col="grey40",lty=1,col=col)
    }
   if(map.show) vis.map(x,col,map.type)
   par(fig=par0$fig,mar=par0$mar,bty="n",xaxt="n",yaxt="n",xpd=FALSE,new=TRUE)
   plot.zoo(x,plot.type=plot.type,type="n",xlab="",ylab="",
            xlim=xlim,ylim=ylim,new=FALSE)
  }
}

vis.map <- function(x,col='red',map.type='points') {
#  print('vis.map')
       xrange <- range(lon(x)) + c(-10,10)
       yrange <- range(lat(x)) + c(-5,5)
       data(geoborders)
       lon <- geoborders$x
       lat <- geoborders$y
       ok <- lon>(min(xrange)-1) & lon<(max(xrange)+1) &
             lat>(min(yrange)-1) & lat<(max(yrange)+1)
       lon2 <- attr(geoborders,"borders")$x
       lat2 <- attr(geoborders,"borders")$y
       ok2 <- lon2>(min(xrange)-1) & lon2<(max(xrange)+1) &
              lat2>(min(yrange)-1) & lat2<(max(yrange)+1)
       par(fig=c(0.76,0.97,0.76,0.97),new=TRUE, mar=c(0,0,0,0),
           xpd=NA,col.main="grey",bty="n")
       plot(lon[ok],lat[ok],lwd=1,col="black",type='l',
            xlab=NA,ylab=NA,axes=FALSE)
       axis(1,mgp=c(3,.5,0),cex.axis=0.75)
       axis(2,mgp=c(2,.5,0),cex.axis=0.75)
       lines(lon2[ok2],lat2[ok2],col = "pink",lwd=1)
       if (map.type=="points") {
         points(lon(x),lat(x),pch=21,cex=1,col=col,bg=col,lwd=1)
       } else if (map.type=="rectangle") {
         rect(min(lon(x)),min(lat(x)),max(lon(x)),max(lat(x)),
              border="black",lwd=1,lty=2)
       }
}

plot.eof <- function(x,new=FALSE,xlim=NULL,ylim=NULL,
                     pattern=1,what=c("pc","eof","var"),
                     colbar=list(pal=NULL,rev=FALSE,n=10,
                         breaks=NULL,type="p",cex=2,show=TRUE,
                         h=0.6,v=1,pos=0.05),
                     verbose=FALSE,...) {
  if (verbose) print(paste('plot.eof',paste(what,collapse=',')))
  if (inherits(x,"comb"))
    plot.eof.comb(x,new=new,xlim=xlim,ylim=ylim,
                  pattern=pattern,what=what,colbar=colbar,verbose=verbose,...) else
  if (inherits(x,"field"))
    plot.eof.field(x,new=new,xlim=xlim,ylim=ylim,
                   pattern=pattern,what=what,colbar=colbar,verbose=verbose,...) else
    print("x does not have 'comb' or 'field' aspects...")
}




plot.eof.field <- function(x,new=FALSE,xlim=NULL,ylim=NULL,pattern=1,
                           what=c("pc","eof","var"),## colbar=NULL,
                           verbose=FALSE,...) {
  if (verbose) print(paste('plot.eof.field',paste(what,collapse=',')))
  n <- pattern
  what <- tolower(what)
  #str(pattern); stop("HERE")
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
          map(x,pattern=pattern,verbose=verbose,...) ## AM formely new=FALSE colbar=colbar,
      } else if (inherits(x,'pca')) {
          par(fig=c(0,0.5,0.5,1))
          main1 <- paste('Leading EOF#',pattern, ' (',
                         round(var.eof[pattern],digits=2),"%)",sep='')
          map(x,pattern=pattern,verbose=verbose,...) ## colbar=colbar,
          title(main=src(x)[1],cex.main=0.6,col.main="grey40",adj=0,line=0)
          title(main=main1,cex.main=0.8)
      }
  }
  ##  if (length(grep('pc',what))>0) result <- as.station(x) else
#  if (length(grep('var',what))>0) result <- attr(x,'tot.var')
    
  ylab <- paste("PC",n)
  main <- paste('First',n,"leading EOFs: ", ## attr(x,'longname')
                 round(sum(var.eof[1:n]),1),"% of variance")
  ## browser()
  if (length(grep('var',what))>0) {
    par(new=TRUE,fig=c(0.5,1,0.5,1))##,xaxt="s",yaxt="s")fig=c(0.5,0.95,0.5,0.975) 
    plot.eof.var(x,pattern=pattern,new=FALSE,cex.main=0.8,cex.axis=0.9,bty="n")
  }
  
  #print(main)
  #browser()
  if (length(grep('pc',what))>0) {
    ##par(bty="n", ##,xaxt="s",yaxt="s",xpd=FALSE,
      par(fig=c(0.05,1,0.025,0.475),new=TRUE) ##,cex.axis=0.9,cex.lab=1) ##(0.05,0.95,0.02,0.45)
      main <- paste('Leading PC#',pattern,' of ',attr(x,'longname'),
                 " - Explained variance = ",round(var.eof[pattern],digits=2),
                    "%",sep='')
      
      plot.zoo(x[,n],lwd=2,ylab=ylab,main=main,xlim=xlim,ylim=ylim,
               cex.main=0.8,bty="n",cex.axis=0.9,cex.lab=1)
      #axis(1,at=pretty(index(x),n=10),labels=,cex.axis=0.9)
      grid()
  }
 
  par(fig=c(0,1,0,0.55),new=TRUE, mar=c(1,1,1,1),xaxt="n",yaxt="n",bty="n")
  plot(c(0,1),c(0,1),type="n",xlab="",ylab="")

  varnm <- varid(x)
  legend(0,0.83,varnm,bty="n",cex=0.8,ncol=2,text.col="grey40")
  
  par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
      fig=c(0,1,0.1,1),new=TRUE)
  par(fig=c(0,1,0,0.1),new=TRUE, mar=c(0,0,0,0))  
}


plot.eof.comb <- function(x,new=FALSE,xlim=NULL,ylim=NULL,
                          pattern=1,col=c("red"),
                          what=c("pc","eof","var"),colbar=NULL,verbose=FALSE,...) {
  if (verbose) print("plot.eof.comb")
  n <- pattern
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100* D^2/tot.var

  if (length(what)==3) mfrow <- c(2,2) else
  if (length(what)==2) mfrow <- c(2,1)
  
  if (new) dev.new()
  #par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)
  par(mfrow=mfrow)

  if (length(grep('eof',what))>0) {
    par(fig=c(0,0.5,0.5,1))
    map(x,pattern=pattern,verbose=verbose,colbar=colbar,...)
  }

  n.app <- attr(x,'n.apps')
  col <- rep(col,n.app)
  src <- rep("",n.app+1)
  src[1] <- attr(x,'source')
  ylab <- paste("PC",1:n)
  main <- paste("EOF: ",n,"accounts for",
                round(var.eof[n],1),"% of variance")

  if (length(grep('var',what))>0)  {
#    par(xaxt="s",yaxt="s")
#    plot.eof.var(x,new=FALSE,cex.main=0.7)
    par(new=TRUE,fig=c(0.5,1,0.5,1))##,xaxt="s",yaxt="s")fig=c(0.5,0.95,0.5,0.975) 
    plot.eof.var(x,new=FALSE,cex.main=0.8,cex.axis=0.9,bty="n")
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

  if (length(grep('pc',what))>0) {
#    par(bty="n",xaxt="s",yaxt="s",xpd=FALSE,
#      fig=c(0.1,0.9,0.1,0.5),new=TRUE,cex.axis=0.6,cex.lab=0.6)
#    plot.zoo(x[,n],lwd=2,ylab=ylab,main=main,sub=attr(x,'longname'),
#                                          xlim=xlim,ylim=ylim)
      par(fig=c(0.025,1,0.025,0.475),new=TRUE) ##,cex.axis=0.9,cex.lab=1) ##(0.05,0.95,0.02,0.45)
      main <- paste('Leading PC#',pattern,'of ',attr(x,'longname'),
                 " - Explained variance = ",round(var.eof[pattern],digits=2),
                    "%",sep='')
      
      plot.zoo(x[,n],lwd=2,ylab=ylab,main=main,xlim=xlim,ylim=ylim,
               cex.main=0.8,bty="n",cex.axis=0.9,cex.lab=1,xaxt="n")
      axis(1,at=pretty(index(x[,n]),n=10),cex.axis=0.9)    
      grid()
#    par0 <- par()

      ## Plot the common PCs
      for (i in 1:n.app) {
        z <- attr(x,paste('appendix.',i,sep=""))
        lines(z[,n],col=col[i],lwd=2)
        if (verbose) print(attr(z,'source'))
        src[i+1] <- attr(z,'source')
      }
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
    plot.ds.pca(x,verbose=verbose)
    return()
  }
  if (inherits(x,'eof')) {
    plot.ds.eof(x,verbose=verbose)
    return()
  }
  
  unit <- attr(x,'unit')
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
  for (i in 1:length(unit)) {
    if ((unit[i]=='degree Celsius') | (unit[i]=='deg C') | (unit[i]=='degC'))
         unit[i] <- 'degree*C'
  }
  
  if (is.null(ylab))
   ylab <- try(eval(parse(text=paste("ylab <- expression(",varid(x),
                               "*phantom(0)*(",unit,"))"))),silent=TRUE)
  if (inherits(ylab,"try-error")) ylab <- unit(x)
  
  if (verbose)  print(ylab)
  if (is.null(main)) main <- attr(x,'longname')[1]               
  if (is.null(col)) col <- rainbow(length(x[1,]))  
  
  cols <- rep("blue",100)
  model <- attr(x,'model')
  ## browser()
  if (length(what)==2) mfrow <- c(2,1) else
  if (length(what)==1) mfrow <- c(1,1)
  
  if (new) dev.new()
  if (plot.type=="single") new <- TRUE
  par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)

  if (sum(is.element(what,'map'))>0) {
    if (verbose) print('Show map...')
    par(fig=c(0,0.5,0.5,1)) ## par(bty="n",fig=c(0,0.5,0.5,1),mar=c(1,1,1,1))
    map(x,new=FALSE,colbar=list(show=FALSE),verbose=verbose,...)
    points(lon(x),lat(x),lwd=3,cex=1.5)
  }

  if ( (sum(is.element(what,'xval'))>0)  & (!is.null(attr(x,'evaluation'))) ){
    #if (is.null(attr(x,'evaluation'))) attr(x,'evaluation') <- crossval(x)
     par(new=TRUE,fig=c(0.5,1,0.5,1)) ##par(bty="n",fig=c(0.55,0.95,0.55,0.95),mar=c(4,3,1,1),new=TRUE, xaxt='s',yaxt='s',cex.sub=0.7)
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
    text(0,0.5,paste('correlation=',
           round(cor(attr(x,'evaluation')[ok,1],attr(x,'evaluation')[ok,2]),2)),
         pos=4,cex=0.8,col='grey')
  }  else {
    xvalfit <- NULL
  }
  
 
#  print("predict")
#  y <- predict(x)
#  print(dim(y))
  
  #print("HERE")

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
  
  if (is.null(ylim))
    ylim <- range(coredata(x),coredata(y0),y.rng,na.rm=TRUE)
  if (is.null(xlim))
    xlim <- range(index(x),index(y0),x.rng,na.rm=TRUE)


  par(fig=c(0.025,1,0.025,0.475),new=TRUE)
  par(bty="n",fig=c(0,1,0.1,0.5),mar=c(1,4.5,1,1),new=TRUE, xaxt='s',yaxt='s')
  ds <- list(obs=y0)
  plot.zoo(y0,plot.type=plot.type,ylab=ylab,xlab=xlab,
           main=main,xlim=xlim,ylim=ylim,lwd=1,type='b',pch=19)
  par0 <- par()
  grid()
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

  legend(x="topleft",legend=c("Obs.","Cal.","Proj"),bty="n",horiz=TRUE,
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



plot.eof.var <- function(x,pattern=1,new=TRUE,xlim=NULL,ylim=NULL,n=20,...) {
  n <- min(c(n,length(attr(x,'eigenvalues'))))
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100* D^2/tot.var
  nt <- length(index(x)) 
  n.eff <- round(nt * (1.0-attr(x,'max.autocor'))/
                      (1.0+attr(x,'max.autocor')))  
  dD <- D*sqrt(2.0/n.eff)
  main <- paste('First',n,"leading EOFs: ", ## attr(x,'longname')
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
  points(pattern,var.eof[pattern],pch=20,col="red",cex=1.2)
  attr(var.eof,'errorbar') <- cbind(100*(D-dD)^2/tot.var,100*(D+dD)^2/tot.var)
  invisible(var.eof)
}



plot.field <- function(x,is=NULL,it=NULL,FUN="mean",...) {
  #print("plot.field")
  stopifnot(!missing(x),inherits(x,'field'))

  d <- dim(x)
  if (d[2]==1) {
    plot.station(x,...)
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
  plot(z,...)
  z <- attrcp(x,z,ignore=c("longitude","latitude"))
  attr(z,'history') <- history.stamp(x)
  if (inherits(x,'station')) lines(y,col="red",lwd=2)
  invisible(z)
}

plot.pca <- function(y,verbose=FALSE,...) {
  if (verbose) print('plot.pca')
  attr(y,'longname') <- attr(y,'longname')[1]
  plot.eof.field(y,verbose=verbose,new=TRUE,...)
}

plot.ds.pca <- function(y,pattern=1,verbose=FALSE,colbar=NULL,...) {
  if (verbose) print('plot.ds.pca')
  attr(y,'longname') <- attr(y,'longname')[1]
  #par(fig=c(0,0.45,0.5,0.975),new=TRUE)
  par(fig=c(0,0.5,0.5,0.975)) #par(fig=c(0,0.45,0.5,0.975))
  map.pca(y,pattern=pattern,verbose=verbose,new=FALSE,colbar=FALSE,...)
  title(paste("PCA Pattern # ",pattern,sep=""))
  par(fig=c(0.55,0.975,0.5,0.975),new=TRUE)
  map(attr(y,'predictor.pattern'),it=pattern,new=FALSE,
      colbar=colbar,verbose=verbose,
      main=paste("EOF Pattern # ",pattern,sep=""))
  #title(paste("EOF Pattern # ",pattern,sep=""))
  if (!is.null(attr(y,'evaluation'))) {
    par(fig=c(0.05,0.45,0.05,0.475),new=TRUE)
    plot(attr(y,'evaluation')[,1],attr(y,'evaluation')[,2],
         main='Cross-validation',xlab='original data',
         ylab='prediction',pch=19,col="grey")
    lines(range(c(attr(y,'evaluation')),na.rm=TRUE),
          range(c(attr(y,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(attr(y,'evaluation')[,1]),
                      x=coredata(attr(y,'evaluation')[,2]))
    xvalfit <- lm(y ~ x, data = cal)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    par(fig=c(0.55,0.975,0.05,0.475),new=TRUE)
    plot(attr(y,'original_data')[,pattern],lwd=2,type='b',pch=19)
    lines(zoo(y[,pattern]),lwd=2,col='red',type='b')
    legend(x=index(attr(y,'original_data')[,pattern])[1],
           y=max(attr(y,'original_data')[,pattern],na.rm=TRUE)+
               diff(range(attr(y,'original_data')[,pattern]))/10,
           legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
           lwd=c(2,2),pch=c(21,19),bty="n")
  } else {
    par(fig=c(0.05,0.975,0.05,0.475),new=TRUE)
    plot(attr(y,'original_data')[,pattern],lwd=2,type='b',pch=19)
    lines(zoo(y[,pattern]),lwd=2,col='red',type='b')
    xvalfit <- NULL
  }
}

plot.ds.eof <- function(y,pattern=1,verbose=FALSE,colbar=list(show=FALSE),...) {
  if (verbose) print('plot.ds.eof')
  attr(y,'longname') <- attr(y,'longname')[1]
  par(fig=c(0,0.5,0.5,1),mar=c(3,5,4.2,1),mgp=c(3,0.5,0.5))
  map.eof(y,pattern=pattern,verbose=verbose,new=FALSE,colbar=colbar,
          main=paste("Predictand EOF pattern # ",pattern,sep=""),...)
  par(fig=c(0.5,1,0.5,1),mar=c(3,4,4.2,1),new=TRUE)
  map(attr(y,'predictor.pattern'),it=pattern,new=FALSE,
      colbar=colbar,verbose=verbose,
      main=paste("Predictor EOF pattern # ",pattern,sep=""))
  #title(paste("EOF Pattern # ",pattern,sep=""))
  if (!is.null(attr(y,'evaluation'))) {
    par(fig=c(0,0.5,0,0.48),mar=c(3,4.5,3,1),new=TRUE)
    pc.obs <- attr(y,'evaluation')[,1*pattern]
    pc.ds <- attr(y,'evaluation')[,1*pattern+1]
    plot(pc.obs,pc.ds,main='Cross-validation',xlab='original data',
         ylab='prediction',pch=19,col="grey",
         xlim=range(pc.obs,pc.ds),ylim=range(pc.obs,pc.ds))
    lines(range(c(attr(y,'evaluation')),na.rm=TRUE),
          range(c(attr(y,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(pc.obs),x=coredata(pc.ds))
    xvalfit <- lm(y ~ x, data = cal)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    text(min(pc.obs)+diff(range(pc.obs))/12,max(pc.ds),
         paste("r =",round(xvalfit$coefficients[2],digits=2)),
         pos=4,cex=0.9)
    par(fig=c(0.5,1,0,0.48),mar=c(3,4.5,3,1),new=TRUE)
    plot(attr(y,'original_data')[,pattern],
         main="PC1",ylab="",
         ylim=range(attr(y,'original_data')[,pattern])*c(1.6,1.6),
         lwd=2,type='b',pch=19)
    lines(zoo(y[,pattern]),lwd=2,col='red',type='b')
    legend(x=index(attr(y,'original_data')[,pattern])[1],
           y=max(attr(y,'original_data')[,pattern],na.rm=TRUE)+
               diff(range(attr(y,'original_data')[,pattern]))/3,
           legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
           lwd=c(2,2),pch=c(21,19),bty="n")
  } else {
    par(fig=c(0,1,0,0.48),mar=c(3,4.5,3,1),new=TRUE)
    plot(attr(y,'original_data')[,pattern],
         main="PC1",ylab="",
         ylim=range(attr(y,'original_data')[,pattern])*c(1.6,1.6),
         lwd=2,type='b',pch=19)
    lines(zoo(y[,pattern]),lwd=2,col='red',type='b')
    legend(x=index(attr(y,'original_data')[,pattern])[1],
           y=max(attr(y,'original_data')[,pattern],na.rm=TRUE)+
               diff(range(attr(y,'original_data')[,pattern]))/3,
           legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
           lwd=c(2,2),pch=c(21,19),bty="n")
    xvalfit <- NULL
  }  
}

vis.pca <- function(y,cex=1.5,new=TRUE) {

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
  data(geoborders,envir=environment())
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
  axis(1,at=seq(0,1,length=length(nl)),label=nl)

  par(mar=c(1,0,0,0),fig=c(0.1,0.3,0.32,0.35),new=TRUE,cex.axis=0.6,xaxt="n")
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),label=nl)

  par(mar=c(1,0,0,0),fig=c(0.6,0.8,0.665,0.695),new=TRUE,cex.axis=0.6,xaxt="n")
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),label=nl)

  par(mar=c(1,0,0,0),fig=c(0.6,0.8,0.32,0.35),new=TRUE,cex.axis=0.6,xaxt="n")
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),label=nl)
  
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


plot.cca <- function(x,icca=1,colbar=list(pal=NULL,rev=FALSE,n=10,
                        breaks=NULL,type="p",cex=2,show=TRUE,
                        h=0.6, v=1,pos=0.05),verbose=FALSE,...) {
  if (verbose) print("plot.cca")
  dev.new()
  par(mfrow=c(2,2),bty="n",xaxt="n",yaxt="n")
  map.cca(x,icca=icca,colbar=colbar,verbose=verbose,...)

  w.m <- zoo((x$w.m[,icca]-mean(x$w.m[,icca],na.rm=TRUE))/
             sd(x$w.m[,icca],na.rm=TRUE),order.by=x$index)
  v.m <- zoo((x$v.m[,icca]-mean(x$v.m[,icca],na.rm=TRUE))/
             sd(x$v.m[,icca],na.rm=TRUE),order.by=x$index)
  r <- cor(x$w.m[,icca],x$v.m[,icca])
  par(bty="n",xaxt="s",yaxt="s",xpd=FALSE,mar=c(2,1.5,1.5,0.5),
      fig=c(0.02,1,0.1,0.45),new=TRUE,cex.axis=0.8,cex.lab=0.8)
  plot(w.m,col="blue",lwd=2,
       main=paste("CCA pattern ",icca," for ",varid(x),
         "; r= ",round(r,2),sep=""),
       xlab="",ylab="")
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

plot.diagnose.comb.eof <- function(x,xlim=NULL,ylim=NULL,verbose=FALSE,...) {
  stopifnot(!missing(x), inherits(x,"diagnose"),
            inherits(x,"eof"),inherits(x,"comb"))
  dev.new()
  par(bty="n")
  if (is.null(xlim)) xlim <- range(abs(c(0,1,x$mean.diff)),na.rm=TRUE)
  if (is.null(ylim)) ylim <- range(c(-1,1,1-x$sd.ratio),na.rm=TRUE)
  wt <- 0:360
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

plot.diagnose.matrix <- function(x,xlim=NULL,ylim=NULL,verbose=FALSE,...) {
  if (verbose) print('plot.diagnose.matrix')
  x <- as.data.frame(x)
  par(bty="n")
  if (is.null(xlim)) xlim <- range(abs(c(0,1,x$mean.diff)),na.rm=TRUE)
  if (is.null(ylim)) ylim <- range(c(-1,1,1-x$sd.ratio),na.rm=TRUE)
  wt <- 0:360
  dev.new()
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
                                     map.type="points",verbose=FALSE,...) {
  if (verbose) print('plot.diagnose.dsensemble')
  Y <- -round(200*(0.5-pbinom(x$outside,size=x$N,prob=0.1)),2)
  X <- -round(200*(0.5-pnorm(x$deltaobs,mean=mean(x$deltagcm),
                             sd=sd(x$deltagcm))),2)
  if (new) {
    dev.new()
    par(bty="n",fig=c(0.05,0.95,0,0.95),mgp=mgp)
  }
  if (!is.null(cex)) par(cex=cex)
  plot(c(-100,100),c(-100,100),type="n",axes=FALSE,mgp=mgp,ylab="",xlab="")
  mtext("trend",side=1,line=1.5,cex=par("cex"))
  mtext("magnitude",side=2,line=2,cex=par("cex"))
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
  mtext("ensemble > obs",side=1,line=0,adj=0,cex=par("cex")*0.75)
  mtext("ensemble < obs",side=1,line=0,adj=1,cex=par("cex")*0.75)
  mtext("ensemble > obs",side=2,line=0.5,adj=0,cex=par("cex")*0.75)
  mtext("ensemble < obs",side=2,line=0.5,adj=1,cex=par("cex")*0.75)  
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
      data(geoborders)
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

 
plot.xval <- function(x,...) {
  dev.new()
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

plot.dsensemble.pca <- function(x,pts=FALSE,showci=TRUE,showtrend=TRUE,it=0,
                               envcol=rgb(1,0,0,0.2),legend=TRUE,verbose=FALSE,...) {
  if (verbose) print("plot.dsensemble.pca")
  stopifnot(inherits(x,'dsensemble') & inherits(x,'pca'))
  d <- index(x[[3]])
  pc <- x[3:length(x)]
  pc <- array(unlist(pc), dim = c(dim(pc[[1]]), length(pc)))
  pc <- lapply(split(pc, arrayInd(seq_along(pc),dim(pc))[,2]),array,dim=dim(pc)[-2])
  fn <- function(x) {
    x <- zoo(x,order.by=d)
    class(x) <- c("dsensemble","station","zoo")
    invisible(x)
  }
  pc <- lapply(pc,fn)
  for (i in 1:length(pc)) {
    attr(pc[[i]],"station") <- x[[2]][,i]
  }
  plot(pc[[1]])
}

plot.dsensemble <-  function(x,pts=FALSE,showci=TRUE,showtrend=TRUE,it=0,
                             envcol=rgb(1,0,0,0.2),legend=TRUE,
                             map.show=FALSE,map.type="points",new=TRUE,
                             xrange=NULL,yrange=NULL,verbose=FALSE,...) {
  if(verbose) print("plot.dsensemble")
  stopifnot(inherits(x,'dsensemble'))
  #print("subset") 
  ## browser()
  #if (inherits(x,'pca')) {
  #  plot.dsensemble.pca(x,)
  #} else {
    
  if (!inherits(attr(x,'station'),'annual')) z <- subset(x,it=it) else
    z <- x
  #print("diagnose")
  ##browser()
  diag <- diagnose(z,plot=FALSE)
  
  y <- attr(z,'station')
  attr(y,'standard.error') <- NULL
  
  d <- dim(z)
  #str(z); str(y); print(class(z)); print(year(z))
  #browser()
  index(y) <- year(y)
  #print("---")

  pscl <- c(0.9,1.1)
  if (max(coredata(z),na.rm=TRUE) < 0) pscl <- rev(pscl)
  args <- list(...)
  #print(names(args))
  ixl <- grep('xlim',names(args))
  if (length(ixl)==0) xlim <- range(year(z)) else
                      xlim <- args[[ixl]]
  iyl <- grep('ylim',names(args))
  if (length(iyl)==0) ylim <- pscl*range(coredata(z),na.rm=TRUE) else
                      ylim <- args[[iyl]]  
  #print("...")
  if(new) dev.new()
  plot(y,type="b",pch=19,xlim=xlim,ylim=ylim,col="black",main='')
  grid()
  usr <- par()$usr; mar <- par()$mar; fig <- par()$fig
  t <- index(z)
  
  if (pts) for (i in 1:d[2])
    points(year(t),coredata(z[,i]),pch=19,col="red",cex=0.3)

  # Produce a transparent envelope
  nt <- length(index(z))
  t2 <- c(year(t),rev(year(t)))
  col <- rgb(rep(1,49),seq(0.95,0.1,length=49),seq(0.95,0.1,length=49))
  
  for (ii in 1:49) {
    qp1 <- qnorm(1-ii/50,mean=coredata(diag$mu),sd=coredata(diag$si))
    qp2 <- qnorm(ii/50,mean=coredata(diag$mu),sd=coredata(diag$si))
    ci <- c(qp1,rev(qp2))
    #print(c(length(ci),length(t2)))
    polygon(t2,ci, col= envcol, ,border=NA)
                                        # transparency not good for hard copies
    #polygon(t2,ci,col=col[ii],border=NA)
  }
  #polygon(t2,c(q95,rev(q05)),col=rgb(0,1,0,0.2),border=NA)
  lines(zoo(diag$mu,order.by=year(diag$mu)),col=rgb(1,0.7,0.7),lwd=3)
  lines(zoo(diag$q05,order.by=year(diag$q05)),col=rgb(0.5,0.5,0.5),lty=2)  
  lines(zoo(diag$q95,order.by=year(diag$q95)),col=rgb(0.5,0.5,0.5),lty=2)
  #browser()
  lines(y,type="b",pch=19)

  # statistics: past trends 
#  i1 <- is.element(year(y),year(z))
#  i2 <- is.element(year(z),year(y))
#  obs <- data.frame(y=y[i1],t=year(y)[i1])
#  #print(summary(gcm))
#  deltaobs <- lm(y ~ t,data=obs)$coefficients[2]*10  # deg C/decade
#  deltagcm <- rep(NA,d[2])
#  for (j in 1:d[2]) {
#    gcm <- data.frame(y=z[i2,j],t=year(z)[i2])
#    deltagcm[j] <- lm(y ~ t,data=gcm)$coefficients[2]*10  # deg C/decade
#  }
#  robs <- round(100*sum(deltaobs < deltagcm)/d[2])
  #print(deltaobs); print(deltagcm); print(order(c(deltaobs,deltagcm))[1])
#
#  # number of points outside conf. int. (binom)
#  above <- y[i1] > q95[i2]
#  below <- y[i1] < q05[i2]
#  outside <- sum(above) + sum(below)
#  #print(outside); print(pbinom(outside,size=sum(i1),prob=0.1))
  
  #str(diag$y); browser()
  index(diag$y) <- year(diag$y)
  #points(diag$y,cex=0.5,col="green")
  #points(diag$y[diag$above | diag$below],col="grey")
  outside <- diag$above | diag$below
  points(zoo(coredata(diag$y)[which(outside)],
             order.by=year(diag$y)[which(outside)]),col="grey")
  #browser()
  
#  if (FALSE) {
#  if (showci) {
#    par(fig=c(0.15,0.4,0.7,0.8),new=TRUE, mar=c(0,0,0,0),xaxt="n",yaxt="n",bty="n",
#        cex.main=0.75,xpd=NA,col.main="grey")
#    par(fig=c(0.05,0.25,0.85,0.95),new=TRUE, mar=c(0,0,0,0),
#        xaxt="n",yaxt="n",bty="n",
#        cex.main=0.75,xpd=NA,col.main="grey")
#    plot(0:(diag$N/2),dbinom(x=0:(diag$N/2), size=diag$N, prob=0.1),
#         type="l",col="grey20",
#         main="Obs. outside simulated 90%",xlab="",ylab="",lwd=6,
#         ylim=c(0,1.25*max(dbinom(x=0:(diag$N/2), size=diag$N, prob=0.1))))
#    lines(0:(diag$N/2),dbinom(x=0:(diag$N/2), size=diag$N, prob=0.1),
#          col="grey70",lwd=5)
#    points(diag$outside,dbinom(x=diag$outside, size=diag$N, prob=0.1),pch=19)
#  }

  if (showtrend) {
#    par(fig=c(0.6,0.9,0.25,0.4),new=TRUE, mar=c(0,0,0,0),xaxt="s",yaxt="n",bty="n",
#        cex.main=0.75,xpd=NA,col.main="grey30")
    par(fig=c(0.23,0.45,0.78,0.98),new=TRUE, mar=c(0,0,0,0),xaxt="s",yaxt="n",bty="n",
        cex.main=0.75,xpd=NA,col.main="grey30")
#    h <- hist(diag$deltagcm,plot=FALSE)
#    hist(diag$deltagcm,freq=FALSE,col="grey80",lwd=2,border="grey",
#         ,main="Obs. and simulated trends",
#        xlab="",ylab="",ylim=c(0,1.2*max(h$density)))
#    lines(seq(min(diag$deltagcm),max(diag$deltagcm),length=300),
#         dnorm(seq(min(diag$deltagcm),max(diag$deltagcm),length=300),
#               mean=mean(diag$deltagcm),sd=sd(diag$deltagcm)),
#         col="grey20",lwd=7)
#    lines(seq(min(diag$deltagcm),max(diag$deltagcm),length=300),
#         dnorm(seq(min(diag$deltagcm),max(diag$deltagcm),length=300),
#               mean=mean(diag$deltagcm),sd=sd(diag$deltagcm)),
#         col="pink",lwd=5)
#  #lines(h$mids,h$density)
#    points(diag$deltaobs,dnorm(diag$deltaobs,mean=mean(diag$deltagcm),
#    sd=sd(diag$deltagcm)),pch=19)
    plot(diag,map.show=FALSE,new=FALSE,cex=0.75)
  } 
  
  if (map.show) {
    if(verbose) print("add map") 
    if(is.null(xrange) & !is.null(attr(z,"lon"))) {
      xrange <- range(attr(z,"lon")) + c(-15,15)
    }
    if(is.null(yrange) & !is.null(attr(z,"lat"))) {
      yrange <- range(attr(z,"lat")) + c(-10,10)
    }
    if (!is.null(xrange) & !is.null(xrange)) {
      data(geoborders)
      lon <- geoborders$x
      lat <- geoborders$y
      ok <- lon>(min(xrange)-1) & lon<(max(xrange)+1) &
      lat>(min(yrange)-1) & lat<(max(yrange)+1)
      lon2 <- attr(geoborders,"borders")$x
      lat2 <- attr(geoborders,"borders")$y
      ok2 <- lon2>(min(xrange)-1) & lon2<(max(xrange)+1) &
      lat2>(min(yrange)-1) & lat2<(max(yrange)+1)
      par(fig=c(0.7,0.95,0.78,0.98),new=TRUE, mar=c(0,0,0,0),
        cex.main=0.75,xpd=NA,col.main="grey",bty="n")
      plot(lon[ok],lat[ok],lwd=1,col="black",type='l',xlab=NA,ylab=NA,
         axes=FALSE)
      axis(1,mgp=c(3,.5,0))
      axis(2,mgp=c(2,.5,0))
      lines(lon2[ok2],lat2[ok2],col = "pink",lwd=1)
      if("points" %in% map.type) {
        points(attr(z,"lon"),attr(z,"lat"),pch=21,cex=1,col='black',bg='red',lwd=1)
      }
      if("rectangle" %in% map.type) {
        rect(min(attr(z,"lon")),min(attr(z,"lat")),max(attr(z,"lon")),
             max(attr(z,"lat")),lwd=1,col=NA,border='red',lty=2)
      }
    } 
  }  
  # finished plotting

  if (legend) {
    par(fig=c(0.1,0.5,0.2,0.25),new=TRUE,mar=c(0,0,0,0),xaxt="n",yaxt="n",bty="n")
    #par(fig=c(0.1,0.5,0.65,0.70),new=TRUE, mar=c(0,0,0,0),xaxt="n",yaxt="n",bty="n")
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
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
  par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
      fig=c(0,1,0.1,1),new=TRUE)
  par(fig=fig,new=TRUE, mar=mar)
  plot(usr[1:2],usr[3:4],
       type="n",ylab="",xlab="",xlim=usr[1:2],ylim=usr[3:4])

  # target: perfect score is bull's eye
  # from diagnostics?
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
    newFig()
    plot(ssa$d,main=main,sub=sub,ylab="Singular value",pch=20,col="grey50")
    points(ssa$d)
    grid()

    newFig()
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


    newFig()
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

plot.nevents <- function(x,verbose=FALSE,...) {
  # Plot the results from 
  if (verbose) print('plot.nevents')
  par(bty='n')
  if (is.T(attr(x,'observation')))
    col <- c(rgb(0.5,0.5,0.7,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.8,0.5,0.8,0.5),
         rgb(0.3,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.6,0.3,0.6,0.5)) else
    col <- c(rgb(0.3,0.3,0.6,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.5,0.5,0.7,0.5),
         rgb(0.6,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.8,0.5,0.8,0.5))
  plot.zoo(x,plot.type='single',lwd=5,main=loc(x),
       xlab="",ylab=attr(x,'info'),col=col,...)
  grid()
  points(attr(x,'observation'),pch=19)
  lines(attr(x,'nwd.pre'),col=rgb(0.5,0.5,0.5,0.5))
}

