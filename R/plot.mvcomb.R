#' @export
plot.mvcomb <- function(x, ..., verbose=FALSE) {
  if (verbose) print("plot.mvcomb")
  stopifnot(inherits(x, "mvcomb"))
  if(inherits(x, "ds")) {
    plot.mvcomb.ds(x, ..., verbose=verbose)
  } else if(inherits(x, "eof")) {
    plot.mvcomb.eof(x, ..., verbose=verbose)
  } else {
    plot.mvcomb.field(x, ..., verbose=verbose)
  }
}

plot.mvcomb.field <- function(x,...,is=NULL,it=NULL,FUN="mean",type="l",
                              col=NULL,lty=NULL,pch=NULL,ylab=NULL,lwd=2,
                              ylim=c(-1,1),fig=NULL,mar=c(4.5,4.5,0.75,0.5),
                              map.show=FALSE,new=TRUE,add=FALSE,verbose=FALSE) {
  if(verbose) print("plot.mvcomb.field")
  ## What should this function plot?
  x <- subset(x, it=it, is=is)
  id <- attr(x, "id")
  if(is.null(col)) col <- colscal(pal="cat", n=length(unique(id)))
  if(is.null(pch)) pch <- seq(1, length(unique(id)))
  if(is.null(lty)) lty <- seq(1, length(unique(id)))
  if(is.null(lwd)) lwd <- seq(2, length(unique(id)))
  if(length(col)==1) col <- rep(col, length(unique(id)))
  if(length(pch)==1) pch <- rep(pch, length(unique(id)))
  if(length(lty)==1) lty <- rep(lty, length(unique(id)))
  if(length(lty)==1) lwd <- rep(lwd, length(unique(id)))
  if(is.null(ylab)) ylab <- paste("Normalized",
              paste(unlist(attr(x, "variable")), collapse=" & "))
  if(is.null(fig) & new) {
    fig <- c(0,1,0,0.95)
    if (map.show) fig[4] <- 0.8
  }
  if(new) dev.new()
  if(!is.null(fig)) par(cex.axis=1,fig=fig,mar=mar)
  par(bty="n",xaxt="s",yaxt="s",xpd=FALSE,new=add)
  y <- list(); j <- 1
  for(i in unique(id)) {
    xi <- apply(coredata(x)[,id==i], 1, FUN=FUN)
    if(j==1) {
      plot(xi, col=col[j], pch=pch[j], lty=lty[j], lwd=lwd[j],
           type=type, ylab=ylab, ylim=ylim, ...)
    } else {
      lines(xi, col=col[j], pch=pch[j], lty=lty[j], lwd=lwd[j], 
            type=type, ...)
    }
    y[[paste0("var",i)]] <- xi
    j <- j+1
  }
  if(map.show) {
    for(j in seq_along(unique(id))) {
      attr(xi,"longitude") <- attr(x,"longitude")[[j]]
      attr(xi,"latitude") <- attr(x,"latitude")[[j]]
      if(j==1) {
        vis.map(xi,col=col[[j]],map.type=NULL,
                xrange=NULL,yrange=NULL,cex=1,
                add.text=FALSE,cex.axis=NULL,
                map.insert=TRUE,verbose=verbose)
      } else {
        rect(min(lon(xi)),min(lat(xi)),max(lon(xi)),max(lat(xi)),
             border=col[[j]],lwd=1,lty=2)
      }
    }
  }
  invisible(y)
}

plot.mvcomb.eof <- function(x, ..., new=FALSE, xlim=NULL, ylim=NULL, 
                            ip=1, col=c("red"), alpha=1, cex.axis=0.9, cex.main=0.9, cex.lab=0.9,
                            what=c("pc","eof","var"), colbar=NULL, verbose=FALSE) {
  if (verbose) print("plot.mvcomb.eof")
  
  if(!inherits(x, "eof")) x <- EOF(x)
  ## Save the original graphics settings
  par0 <- par()
  n <- ip
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100*D^2/tot.var
  nv <- length(attr(x, "pattern")) 
  
  #if (length(what)==3) mfrow <- c(2,2) else
  #  if (length(what)==2) mfrow <- c(2,1) else
  #    mfrow <- NULL
  mfrow <- NULL
  add <- FALSE
  if (new) dev.new()
  par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)
  if (!is.null(mfrow)) par(mfrow=mfrow)
  
  if (length(grep('var',what))>0)  {
    if (is.null(mfrow)) par(new=add, fig=c(0,0.4,0.5,1))
    plot.eof.var(x,ip=ip,new=FALSE,cex.main=0.8,cex.axis=0.9,bty="n")
    add <- TRUE
  }
  
  if (length(grep('pc',what))>0) {
    if(is.null(mfrow)) par(fig=c(0.4,1,0.5,1),mar=c(3,3,2,2),new=add)
    ylab <- paste("PC",n)
    main <- paste0('PC#',ip,' of ',paste(attr(x,'longname'), collapse=" & "),
                   "\nExplained variance = ",round(var.eof[ip],digits=2),"%")
    xn <- coredata(x)[,n]
    if (is.null(ylim)) ylim <- range(coredata(xn),na.rm=TRUE) + 
      0.1*c(-1,1)*diff(range(coredata(xn),na.rm=TRUE))
    if(inherits(x,"seasonalcycle")) xaxt <- "n" else  xaxt <- NULL
    if(inherits(index(xn),"PCICt")) {
      caldays <- as.numeric(substr(attr(x,"calendar"),1,3))
      index(xn) <- as.numeric(format(index(x),"%Y")) + 
        (as.numeric(format(index(x),"%j"))+as.numeric(format(index(x),"%H"))/24)/caldays
    }
    plot.zoo(xn, lwd=2, ylab=ylab, main=main, xlim=xlim, ylim=ylim,
             cex.main=cex.main, bty="n", cex.axis=cex.axis,
             cex.lab=cex.lab, xaxt=xaxt)
    if(inherits(x,"seasonalcycle")) axis(1,at=seq(1,12),labels=month.abb,
                                         cex.axis=cex.axis,las=2)
    grid()
    add <- TRUE
  }
  
  if(length(grep('eof',what))>0) {
    for(i in seq(nv)) {
      if(is.null(mfrow)) par(fig=c((i-1)/nv, i/nv, 0, 0.5), new=add)
      map.mvcomb(x, ip=ip, iv=i, new=FALSE, lab=paste("normalized", attr(x, "variable")[i]))
      add <- TRUE
    }
  }
  
  ## Reset the graphics settings that have been changed to original
  par(fig=par0$fig, mar=par0$mar, mgp=par0$mgp, xaxt=par0$xaxt , yaxt=par0$yaxt)
  
}

plot.mvcomb.ds <- function(x,...,plot.type="multiple",what=c("map","ts",'xval'),new=TRUE,
                           ip=1,lwd=1,type='l',pch=0,main=NULL,col=NULL,
                           colbar=list(pal=NULL,rev=FALSE,n=10,
                                       breaks=NULL,type="p",cex=2,show=TRUE,
                                       h=0.6, v=1,pos=0.05),
                           xlim=NULL,ylim=NULL,xlab="",ylab=NULL,verbose=FALSE) {
  if (verbose) print(paste('plot.mvcomb.ds',paste(what,collapse=',')))
  
  if (inherits(x,'pca')) {
    if(!is.null(dim(x))) y <- zoo(coredata(x)[,ip],order.by=index(x)) else y <- x
    y <- attrcp(x,y)
    attr(y, "model") <- attr(x, "model")[[ip]]
    attr(y, "pattern") <- attr(x, "pattern")[[ip]]
    attr(y, "evaluation") <- attr(x, "evaluation")[,seq(2*ip-1,2*ip)]
    attr(y, "variable") <- paste0("PC",ip,".",attr(x, "variable"))
    attr(y, "unit") <- "unitless"
    org <- zoo(attr(x, "original_data")[,ip], 
               order.by=index(attr(x, "original_data")))
    attr(org, "variable") <- paste0("PC",ip,".",attr(x, "variable"))
    attr(org, "unit") <- "unitless"
    attr(y, "original_data") <- org
    attr(y, "fitted_values") <- zoo(attr(x, "fitted_values")[,ip], 
                                    order.by=index(attr(x, "fitted_values")))
    class(y) <- class(x)[!class(x) %in% "pca"]
    z <- plot.mvcomb.ds(y, verbose=verbose, plot.type=plot.type,
                        what=what, new=new, lwd=lwd, type=type, pch=pch, main=main,
                        col=col, colbar=colbar, xlim=xlim, ylim=ylim, 
                        xlab=xlab, ylab=ylab, ...)
    invisible(z)
  } else {
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
      ncomb <- length(attr(x,"pattern"))
      #par(fig=c(0,0.5,0.5,1))
      #map(x,new=FALSE,colbar=list(show=FALSE),verbose=verbose,...)
      #points(lon(x),lat(x),lwd=3,cex=1.5)
      for(i in seq(1,ncomb)) {
        par(fig=c((i-1)/ncomb,i/ncomb,0.5,1), new=i>1)
        map.mvcomb(x,iv=i,ip=ip,new=FALSE,colbar=list(show=FALSE),verbose=verbose,...)
        points(lon(x),lat(x),lwd=3,cex=1.5)
      }
    }
    
    if ( (sum(is.element(what,'xval'))>0)  & (!is.null(attr(x,'evaluation'))) ) {
      #par(new=TRUE,fig=c(0.5,1,0.5,1))
      par(new=TRUE,fig=c(0,0.5,0,0.5))
      plot(attr(x,'evaluation')[,1],attr(x,'evaluation')[,2],
           main='Cross-validation',xlab='original data',
           ylab='prediction',pch=19,col="grey")
      lines(range(c(attr(x,'evaluation')),na.rm=TRUE),
            range(c(attr(x,'evaluation')),na.rm=TRUE),lty=2)
      cal <- data.frame(y=coredata(attr(x,'evaluation')[,1]),
                        x=coredata(attr(x,'evaluation')[,2]))
      xvalfit <- lm(y ~ x, data = cal)
      abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
      #par(bty="n",fig=c(0.6,0.95,0.48,0.52),mar=c(0,0,0,0),new=TRUE,
      #    xaxt='n',yaxt='n',cex.sub=0.7)
      #plot(c(0,1),c(0,1),type='n',xlab='',ylab='')
      ok <- is.finite(attr(x,'evaluation')[,1]) &
        is.finite(attr(x,'evaluation')[,2])
      text(par()$usr[1] + diff(range(par()$usr[1:2]))/24,
           par()$usr[4] - diff(range(par()$usr[3:4]))/12,
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
    ## KMP 2021-02-26: Plot attribute 'fitted value' instead of coredata
    
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
      #ylim <- range(coredata(x),coredata(y0),y.rng,na.rm=TRUE)
      ylim <- range(attr(x,"fitted_values"),coredata(y0),y.rng,na.rm=TRUE)
    }
    if (is.null(xlim)) {
      xlim <- range(index(x),index(y0),x.rng,na.rm=TRUE)
    }
    
    #par(fig=c(0.025,1,0.025,0.475),new=TRUE)
    #par(bty="n",fig=c(0,1,0.1,0.5),mar=c(1,4.5,1,1),new=TRUE, xaxt='s',yaxt='s')
    par(bty="n",fig=c(0.5,1,0.1,0.5),mar=c(1,4.5,1,1),new=TRUE, xaxt='s',yaxt='s')
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
    #lines(x,col="red",type="l",lwd=lwd)
    lines(attr(x,"fitted_values"),col="red",type="l",lwd=lwd)
    
    cal0 <- data.frame(y=coredata(y0),t=year(y0))
    #cal1 <- data.frame(y=coredata(x),t=year(x))
    cal1 <- data.frame(y=attr(x,"fitted_values"),t=year(x))
    
    trend0 <- lm(y ~ t, data=cal0)
    trend1 <- lm(y ~ t, data=cal1)
    lines(zoo(predict(trend0),order.by=index(y0)),lty=2)
    lines(zoo(predict(trend1),order.by=index(x)),lty=2,col='red')
    
    st0 <- summary(trend0); st1 <- summary(trend1)
    obstrend <- paste('obs. trend: ', round(st0$coefficients[2],2),' (',
                      round(st0$coefficients[2]-2*st0$coefficients[4],2),',',
                      round(st0$coefficients[2]+2*st0$coefficients[4],2),')',
                      attr(x,'unit'),'/decade',sep='')
    dstrend <- paste('obs. trend: ', round(st1$coefficients[2],2),' (',
                     round(st1$coefficients[2]-2*st1$coefficients[4],2),',',
                     round(st1$coefficients[2]+2*st1$coefficients[4],2),')',
                     attr(x,'unit'),'/decade',sep='')
    
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
    #lines(x,col="red",type="l",lwd=lwd)
    lines(attr(x,"fitted_values"),col="red",type="l",lwd=lwd)
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
      #plot.zoo(x,plot.type=plot.type,type="n",ylab="",xlab="",xlim=xlim,ylim=ylim)
      plot.zoo(attr(x,"fitted_values"),plot.type=plot.type,type="n",
               ylab="",xlab="",xlim=xlim,ylim=ylim)
    }
    invisible(list(trend0=trend0,trend1=trend1,xvalfit=xvalfit))
  }
}
