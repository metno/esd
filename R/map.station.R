## Author 	 Rasmus E. Bnestad
## Updated 	 by Abdelkader Mezghani and Kajsa Parding
## Rasmus. E. Benestad - attempt to simplify by splitting up
## Last update   27.07.2017
## Includes	 map.station() ; test.map.station()
## Require 	 geoborders.rda

genfun <- function(x,FUN,verbose=FALSE) {
  if (sum(is.element(names(attributes(x)),FUN))>0){
    ## REB 2015-12-17: Use FUN to colour the symbols according to some attribute:
    FUN <- eval(parse(text=paste("attr(x,'",FUN,"')")))
  } else if (sum(is.element(names(x),FUN))>0){
    ## REB 2015-12-17: Use FUN to colour the symbols according to some list element (stationmeta-objects):
    if (verbose) print('FUN refers to a list element')
    FUN <- eval(parse(text=paste("function(x,...) x$",FUN,sep='')))
    return(FUN)
  }
}

## Simplified function for mapping station objects.
#' @exportS3Method
#' @export map.station
map.station <- function(x=NULL,FUN=NULL, it=NULL,is=NULL,new=FALSE,
                        add=FALSE,projection="lonlat",
                        xlim = NULL, ylim = NULL,zlim=NULL,n=15,
                        col='darkred',bg='orange',
                        colbar= list(pal='t2m',col=NULL,rev=FALSE,n=10,
                                     breaks=NULL,type="p",cex=2,h=0.6, v=1,
                                     pos=0.1,show=FALSE),
                        # col=NULL replaced by palette
                        type=NULL,gridlines=TRUE,
                        lonR=NULL,latR=45,axiR=NULL,verbose=FALSE,
                        cex=2,zexpr="alt",cex.subset=1,
                        add.text.subset=FALSE,showall=FALSE,
                        add.text=FALSE,main=NULL,sub=NULL,
                        height=NULL,width=NULL,
                        cex.main=1,cex.sub=0.75,cex.axis=1,cex.lab=0.9,
                        col.main="black",col.sub="grey",col.border="grey",
                        font.main=1,font.sub=4,
                        pch=19, from=NULL,to=NULL,showaxis=FALSE,
                        border=FALSE,full.names=FALSE,
                        full.names.subset=FALSE, use.old=FALSE,
                        text=FALSE, fancy=TRUE, 
                        na.rm=TRUE,show.val=FALSE,#usegooglemap=FALSE,
                        ##colorbar=TRUE,
                        add.significance=FALSE, pval=0.01, col.pval='black', lwd.pval=2,
                        xlab="lon",ylab="lat",
                        legend.shrink=1,fig=c(0,1,0.05,0.95),
                        mar=rep(2,4),mgp=c(3,1,0),plot=TRUE,...) { 
  if ( (inherits(x,"stationmeta")) | (projection != 'lonlat') | use.old) {#| usegooglemap) {
    map.station.old(x=x,FUN=FUN,it=it,is=is,new=new,projection=projection,
                    xlim=xlim,ylim=ylim,zlim=zlim,n=n,col=col,bg=bg,
                    colbar=colbar,xlab=xlab,ylab=ylab,type=type,gridlines=gridlines,
                    lonR=lonR,latR=latR,axiR=axiR,verbose=verbose,
                    cex=cex,zexpr=zexpr,cex.subset=cex.subset,
                    add.text.subset=add.text.subset,showall=showall,add.text=add.text,
                    height=height,width=width,cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,
                    pch=pch,from=from,to=to,showaxis=showaxis,border=border,
                    full.names=full.names,full.names.subset=full.names.subset,
                    text=text,fancy=fancy,na.rm=na.rm,show.val=show.val,
                    legend.shrink=legend.shrink,...)#,usegooglemap=usegooglemap)
  } else {
    if (verbose) print('map.station - new version')
    if (is.null(add)) new <- FALSE
    if (new) dev.new()
    if ( (!is.null(it)) | (!is.null(is)) ) x <- subset(x,it=it,is=is)
    if (!is.null(FUN)) if (FUN=='trend') {
      FUN <- 'trend.coef'
      if(is.null(colbar$pal)) colbar$pal <- 't2m'
      if(is.null(colbar$rev)) if (is.precip(x)) colbar$rev=TRUE
    } else {
      ## KMP 2018-05-31: Allow user to set different pal than default
      if(is.null(colbar$pal)) colbar$pal <- varid(x)[1]
    }
    if (!is.null(FUN)) {
      if (!(FUN %in% names(attributes(x)))) {
        if(is.null(dim(x))) dim(x) <- c(length(x),1)
        ## The function first tries a function that allows the argument 'na.rm'
        y <- try(apply(coredata(x),2,FUN,na.rm=na.rm))
        ## If not, do it with a function that doesn't have the argument 'na.rm'
        if (inherits(y,"try-error")) y <- apply(coredata(x),2,FUN)
      } else {
        y <- attr(x,FUN); FUN <- NULL
      }    
      if (verbose) {print('Contents of y:'); print(summary(y))}
      colbar <- colbar.ini(y,colbar=colbar)
      if (verbose) print('Set colour scheme')
      wr <- round(strtoi(paste('0x',substr(colbar$col,2,3),sep=''))/255,2)
      wg <- round(strtoi(paste('0x',substr(colbar$col,4,5),sep=''))/255,2)
      wb <- round(strtoi(paste('0x',substr(colbar$col,6,7),sep=''))/255,2)
      col <- rep(colbar$col[1],length(y))
      for (i in 1:length(y)) {
        ii <- round(approx(0.5*(colbar$breaks[-1]+colbar$breaks[-length(colbar$breaks)]),1:length(colbar$col),
                           xout=y[i],rule=2)$y)
        if (is.finite(ii)) {
          if (ii < 1) ii <- 1
          if (ii > length(colbar$col)) ii <- length(colbar$col)
          col[i] <- rgb(wr[ii],wg[ii],wb[ii],0.7)
        } else col[i] <- rgb(0.5,0.5,0.5,0.2)
      }
      if(is.null(colbar$show)) {
        show.colbar <- TRUE
      } else {
        show.colbar <- colbar$show
      }
    } else {
      y <- rep(1,length(lon(x)))
      show.colbar <- FALSE
    }
    if (verbose) print(paste('show.colbar=',show.colbar))
    
    par(mar=mar,mgp=mgp,bty='n',xaxt='n',yaxt='n',cex.axis=0.7,
        col.axis='grey30',col.lab='grey30',las=1)
    ## KMP 2017-07-28: fig creates problems when you want to add map.station as a subplot.
    ## With this solution you have to use add=TRUE and set fig to your subplot or to NULL.
    ## KMP 2023-01-31: Only use fig when it is defined. 
    #if(is.null(fig)) fig <- par()$fig
    if (plot) {
      if(!is.null(fig)) {
        if (is.null(add)) add <- FALSE 
        par(fig=fig,new=add)
      }	
      
      ## Avoid errors when plotting the colorbar with small figure windows
      ## fin collects information about the figure size (in inches)
      #fin <- par()$fin
      
      ## For checking & debugging
      if (verbose) {
        #print(paste('window size=',fin,collapse=' '))
        print(summary(lon(x)))
        print(summary(lat(x)))
        str(col)
      }
      
      if (is.null(xlim)) {
        if (length(is$lon)>1) {
          xlim <- range(is$lon,na.rm=TRUE) + c(-1,1)
        } else {
          xlim <- range(lon(x),na.rm=TRUE) + c(-4,4)
        }
      }
      
      if (is.null(ylim)) {
        if (length(is$lat)>1) {
          ylim <- range(is$lat,na.rm=TRUE) + c(-1.5,1)
        } else {
          ylim <- range(lat(x),na.rm=TRUE) + c(-2.5,2)
        }
      }
      
      plot(lon(x),lat(x),xlim=xlim,ylim=ylim,col=col,pch=pch,cex=cex,new=FALSE,
           cex.lab=cex.lab,xlab=xlab,ylab=ylab)
      if (add.text) text(lon(x),lat(x),substr(loc(x),1,6),cex=cex.lab,col='grey',pos=1)
      
      if(showaxis | gridlines) {
        par(xaxt="s",yaxt="s",las=1,col.axis='grey',col.lab='grey',
            cex.lab=cex.lab,cex.axis=cex.axis)
        if(diff(range(xlim))>10) {
          dlon <- round(round(diff(range(xlim))/4)/5)*5
        } else if(diff(range(xlim))<2) {
          dlon <- round(diff(range(xlim))/3, 2)
        } else {
          dlon <- round(diff(range(xlim))/3)
        }
        if(diff(range(ylim))>10) {
          dlat <- round(round(diff(range(ylim))/4)/5)*5
        } else if(diff(range(ylim))<2) {
          dlat <- round(diff(range(ylim))/3, 2)
        } else {
          dlat <- round(diff(range(ylim))/3)
        }
        axis(3,seq(floor(par("xaxp")[1]/dlon)*dlon,par("xaxp")[2],by=dlon),col='grey')
        axis(4,seq(floor(par("yaxp")[1]/dlat)*dlat,par("yaxp")[2],by=dlat),col='grey')
        if (gridlines) grid()
      }
      
      data("geoborders", envir = environment())
      lines(geoborders$x,geoborders$y,col=col.border, lwd=1.5)
      if (border) lines(attr(geoborders,'border')$x,attr(geoborders,'border')$y,
                        col=adjustcolor(col.border, alpha.f=0.7), lwd=0.75)#'grey')
      
      if (show.colbar) {
        par(xpd=TRUE)
        if(fancy) {
          if (verbose) print('show.colorbar')
          #image.plot(breaks=colbar$breaks,
          #           lab.breaks=colbar$breaks,horizontal = TRUE,
          #           legend.only = T, zlim = range(colbar$breaks),
          #           col = colbar$col, legend.width = 1,
          #           axis.args = list(cex.axis = cex.axis, hadj = 0.5,mgp = c(0, 0.5, 0)), 
          #           border = FALSE)
          #image(colbar$breaks,1:2,cbind(colbar$breaks,colbar$breaks),
          #      col=colbar$col,axes=FALSE)
          #par(mar=c(2,1,2,1),mgp=c(2,0.4,0),cex.axis=cex.axis,col.axis='grey')
          #axis(1,colbar$breaks)
          ## KMP 2023-02-16: testing alternative colorbar
          dy <- diff(ylim)*0.1
          below <- c(min(xlim), min(ylim)-dy/2, max(xlim), min(ylim)+dy/2)
          rect(below[1], below[2], below[3], below[4], 
               col = "white", border = "white")
          col.bar(below[1],below[2],below[3],below[4],
                  colbar$breaks,horiz=TRUE,pch=15,v=1,h=1,
                  col=colbar$col,cex=2,cex.lab=colbar$cex.lab,
                  type=colbar$type,verbose=FALSE,vl=1,border=FALSE)
        } else {
          ## REB 2023-01-31
          ## Make a simple legend for colour scales based on points along the lower part of the figure
          if (verbose) {print('Alternative colour bar'); print(colbar$breaks)}
          nb <- length(colbar$breaks) - 1
          xs <- seq(xlim[1],xlim[2],length=nb)
          par(xpd = TRUE)
          rect(xlim[1]-1,ylim[1]-1.5,xlim[2]+1,ylim[1]+0.5,col=rgb(1,1,1,0.7),border=rgb(0.5,0.5,0.5,0.7),lwd=2)
          points(xs,rep(ylim[1],nb)-0.5,col=colbar$col,pch=19,cex=3)
          if (nb > 15) ib <- c(1,(1:nb)[(1:nb)%%3 == 0],nb) else ib <- 1:nb
          levels <- round(0.5*(colbar$breaks[-1] + colbar$breaks[-(nb+1)])[ib],2)
          text(xs[ib],rep(ylim[1],sum(ib))-0.5,levels,cex=0.7,col='grey40')
        }
      }
      
      if(!is.null(FUN)) if(add.significance & FUN %in% c("trend","trend.coef")) {
        pval.x <- apply(x, 2, trend.pval) # pval.x = the p-values of the trends at each station in the object x
        if(is.null(col.pval)) col.pval <- colscal("gray.colors", n=length(pval)) # col.pval = colors for the thresholds pval
        if(is.null(lwd.pval)) lwd.pval <- pretty(c(2,1), n=length(pval)) # lwd.pval = line width for the thresholds pval
        if(length(col.pval)<length(pmax)) col.pval <- rep(col.pval, length(pval)) # if only one color has been defined but several thresholds
        if(length(lwd.pval)<length(pmax)) lwd.pval <- rep(lwd.pval, length(pval)) # if only line width has been defined but several thresholds
        pmin <- 0
        i <- 1
        for(p in sort(pval)) { # loop through all the thresholds and add points
          points(lon(x)[pval.x>=pmin & pval.x<p], lat(x)[pval.x>=pmin & pval.x<p], 
                 col=col.pval[i], lwd=lwd.pval[i], pch=21, cex=cex)
          pmin <- p
          i <- i+1
        }
      }
      ## Add a title
      if(is.null(main)) {
        if(!is.null(FUN)) {
          main <- paste(esd::varid(x)[1]," (",attr(x,'unit')[1],")",sep="")
          if(is.function(FUN)) {
            main <- paste(as.character(quote(FUN)),"of",main)
          } else {
            main <- paste(as.character(FUN),"of",main)
          }
        } else main <- esd::varid(x)[1]
      }
      title(main=main,sub=sub,line=-1,adj=0,cex.main=cex.main,cex.sub=cex.sub,
            col.main=col.main,col.sub=col.sub,font.main=font.main,font.sub=font.sub)
    }
    if (verbose) print('Organise output')
    if (inherits(x,'station')) {
      dim(y) <- c(1,length(y))
      y <- zoo(y,order.by=1)
      if (verbose) print(dim(y))
      class(y) <- class(x)
      y <- attrcp(x,y)
      attr(y,'period') <- paste(range(index(x)))
    }
    attr(y,'history') <- history.stamp(x)
    invisible(y)
  }
}

###
# Internal function - no need to export
map.station.old <- function (x=NULL,FUN=NULL, it=NULL,is=NULL,new=FALSE,
                             projection="lonlat",
                             xlim = NULL, ylim = NULL,zlim=NULL,n=15,
                             col='darkred',bg='orange',
                             colbar= list(pal='t2m',col=NULL,rev=FALSE,n=10,
                                          breaks=NULL,type="p",cex=2,h=0.6, v=1,
                                          pos=0.1,show=TRUE),
                             # col=NULL replaced by palette
                             xlab=NULL,ylab=NULL,
                             type=NULL,gridlines=TRUE,
                             lonR=NULL,latR=45,axiR=NULL,verbose=FALSE,
                             cex=2,zexpr="alt",cex.subset=1,
                             add.text.subset=FALSE,showall=FALSE,
                             add.text=FALSE,
                             height=NULL,width=NULL,
                             cex.main=1,cex.axis=1,cex.lab=0.6,
                             pch=21, from=NULL,to=NULL,showaxis=FALSE,
                             border=FALSE,full.names=FALSE,
                             full.names.subset=FALSE, 
                             text=FALSE, fancy=TRUE, 
                             na.rm=TRUE,show.val=FALSE,
                             col.subset="red", bg.subset="red",
                             #usegooglemap=FALSE,
                             ##colorbar=TRUE,
                             legend.shrink=1,...) { 
  ##
  if (verbose) {
    print(paste('map.station',FUN))
    print(class(x))
  }
  arg <- list(...)
  if (inherits(x,'station')) {
    attr(x,'unit') <- as.character(attr(x,'unit'))
    attr(x,'variable') <- as.character(varid(x))
    
    ## REB 2016-11-28: some objects contain the attribute 'mean' which gets in the way.
    if (!is.null(FUN)) {
      if ((FUN=='mean') & (!is.null(attr(x,'mean')))) attr(x,'mean') <- NULL
    }
    if (inherits(x,"stationmeta")) {
      x$years <- as.numeric(x$end) - as.numeric(x$start) + 1
      if (!is.null(FUN)) if (FUN=='alt') FUN <- 'altitude'
      if (verbose) print(names(x))
    }
    
    if (!is.null(FUN)) {
      if (is.character(FUN)) {
        if (FUN=="NULL") FUN <- NULL else FUN <- genfun(x,FUN)
      } else if (!is.function(FUN)) {
        x <- FUN
        FUN <- NULL
      }
    }
    if (verbose) print(FUN)
  }
  fig0 <- c(0,1,0,1); mar0 <- rep(2,4)
  par0 <- par(fig=fig0,mar=mar0)
  if ((par()$mfcol[1]> 1) | (par()$mfcol[2]> 1)) new <- FALSE
  
  if (verbose) print(paste("List of arguments in the three-dots listed below ",arg,sep=""))
  
  if (sum(is.element(type,c('fill','contour')))) {
    x0 <- x
    x <- as.field(x)
  }
  
  if ((!is.null(FUN)) & is.character(FUN)) if (FUN=='trend') FUN <- 'trend.coef'
  
  if (verbose) print(paste(projection,'projection'))
  if (projection=="sphere") {
    sphere(x,lonR=lonR,latR=latR,axiR=axiR,
           gridlines=gridlines,xlim=xlim,ylim=ylim,
           col=colbar$col,new=new,FUN=FUN,cex=cex,
           cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,
           verbose=verbose,
           xlab=xlab,ylab=ylab,...)
  } else if (projection=="np") {
    sphere(x,lonR=lonR,latR=90,axiR=axiR,
           gridlines=gridlines,xlim=xlim,ylim=ylim,
           col=colbar$col,new=new,FUN=FUN,
           cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,...) 
  } else if (projection=="sp") {
    sphere(x,lonR=lonR,latR=-90,axiR=axiR,
           gridlines=gridlines,xlim=xlim,ylim=ylim,
           col=colbar$col,new=new,FUN=FUN,
           cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,
           verbose=verbose,...)
    ## else if (projection=="lonlat")
    ##    lonlatprojection(x=X,xlim=xlim,ylim=ylim, n=colbar$n,col=colbar$col,breaks=colbar$breaks,new=new,
    ##                     type=type,gridlines=gridlines,...)
  } else if (projection=="lonlat") {
    data("geoborders", envir = environment())
    if (zexpr == "alt") zexpr <- "sqrt( station.meta$alt/max(station.meta$alt,na.rm=TRUE) )"
    if (verbose & !is.null(x)) print(class(x)) else if (verbose) print('x is null')
    if (!is.null(x)) { 
      if (inherits(x,"stationmeta")) {
        ss <- x
        ss$variable <- apply(as.matrix(ss$element),1,esd2ele)
      } else if (inherits(x,"station")) {
        if (is.null(dim(x))) dim(x) <- c(length(x),1)
        ss <- list(station_id=attr(x,"station_id"),
                   location=attr(x,'location'),
                   country=attr(x,'country'),
                   longitude=attr(x,"longitude"),
                   latitude=attr(x,'latitude'),
                   altitude=attr(x,'altitude'),
                   variable=attr(x,"variable"),
                   longname=attr(x,"longname"),
                   start=rep(start(x),dim(x)[2]),
                   end=rep(end(x),dim(x)[2]),
                   source=attr(x,"source"))
      }
    } else {
      ss <- select.station()
    }
    if (verbose) {
      print('The station metadata:')
      str(ss)
    }  
    
    ## REB 2021-07-021: bugfix for new GHCND metadata
    if (verbose) print('Set element if missing...')
    if ( (is.null(ss$variable)) & !is.null(x$variable) ) ss$variable <- x$variable
    if (is.null(attr(ss,"element"))) ss$element <-apply(as.matrix(ss$variable),1,esd2ele)   
    
    if (verbose) str(ss)
    
    if (!is.null(unlist(is))) { ## highlight a subset of station
      if (verbose) print('HERE')
      if (is.null(is$x)) {
        highlight <- select.station(x=is$x,loc=is$loc,stid = is$stid,
                                    param = is$param, lon = is$lon,
                                    lat = is$lat, alt = is$alt,
                                    cntr = is$cntr, src = is$src,
                                    nmin = is$nmin , it = is$it)
        highlight$variable <- apply(as.matrix(highlight$element),1,
                                    esd2ele)
      } else if (inherits(is$x,"station")) {
        if (is.null(dim(is$x))) dim(is$x) <- c(length(is$x),1)
        highlight <- list(station_id=attr(is$x,"station_id"),
                          location=attr(is$x,'location'),
                          country=attr(is$x,'country'),
                          longitude=attr(is$x,"longitude"),
                          latitude=attr(is$x,'latitude'),
                          altitude=attr(is$x,'altitude'),
                          variable=attr(is$x,"variable"),
                          longname=attr(is$x,"longname"),
                          start=rep(start(is$x),dim(is$x)[2]),
                          end=rep(end(is$x),dim(is$x)[2]),
                          source=attr(is$x,"source"))
      }
      highlight$element <- apply(as.matrix(highlight$variable),1,esd2ele)
    }  else highlight <-  NULL
    
    ## Set negative altitude to NA
    ss$altitude[ss$altitude < 0] <- NA
    
    tte <- "rwb"
    
    ## An attempt to set the size of the symbols automatically,
    ## but this fails if ss is a list.
    if (inherits(x,'station')) {
      nok <- apply(coredata(x),2,FUN='nv')
      if (verbose) {
        print('nok:')
        print(nok)
      }
      if ((is.null(cex)) & !is.null(dim(ss))) {
        cex <- 5/log(dim(ss)[1]) 
      } else {
        if (is.null(cex)) {
          cex <- 5/log(length(ss[[1]])) 
        } else if (cex==0) {
          cex <- 1.25*nok/max(nok,na.rm=TRUE)
          if (cex<0) cex <- abs(cex)*nok/max(nok,na.rm=TRUE)
        }
      }
    }
    
    ## Select a subdomain in the x-axis
    if (verbose) {
      print('cex:')
      print(cex)
    }
    if (is.null(xlim)) {
      if ((is.null(highlight) | showall)) {
        if(length(is$lon) > 1) {
          xlim <- range(ss$longitude, na.rm = TRUE) + c(-2,2)
        } else {
          xlim <- range(ss$longitude, na.rm = TRUE) + c(-5,5)
        }
      } else {
        if (length(is$lon) > 1) {
          xlim <- range(highlight$longitude, na.rm = TRUE) + c(-2,2)
        } else {
          xlim <- range(highlight$longitude, na.rm = TRUE) + c(-5,5)
        }
      }
    }
    ## Select a subdomain in the y-axis
    if (is.null(ylim)) {
      if ((is.null(highlight) | showall)) {
        if (length(is$lat) > 1) {
          ylim <- range(ss$latitude, na.rm = TRUE) + c(-1,1)
        } else {
          ylim <- range(ss$latitude, na.rm = TRUE) + c(-2,2)
        } 
      } else {
        if (length(is$lat) > 1) {
          ylim <- range(highlight$latitude, na.rm = TRUE) + c(-1,1)
        } else {
          ylim <- range(highlight$latitude, na.rm = TRUE) + c(-2,2)
        }
      }
    }
    
    ## scaling factor to apply on cex ...
    if (verbose) print('scale:')
    if (!inherits(x,"stationmeta") & !is.null(attr(x,'na'))) {
      scale <- attr(x,'na')
    } else {
      scale <- 1
    }
    ##print(par()$fig)
    ## Transform x using FUN and insert color bar
    ##
    
    if (verbose) {
      print('FUN:')
      print(FUN)
    }
    if (!is.null(FUN)) {
      if (is.function(FUN)) {
        if (verbose) print('function')
        y <- as.numeric(FUN(x))
        colbar <- colbar.ini(y,colbar=colbar,verbose=verbose)
      } else if (is.character(FUN)) {
        if (verbose) print('string')
        if (is.element(FUN,c('lon','lat','alt'))) {
          if (verbose) print(paste('FUN=',FUN,'(lon/lat/alt)'))
          if (is.character(FUN)) eval(parse(text=paste('y <-',FUN,'(x)',sep="")))               
        } else {
          if (verbose) print(FUN)
          if (is.element("na.rm",names(formals(FUN))) |
              is.element("...",names(formals(FUN))) |
              (is.element(FUN,c("max","min","sum","mean","sd")))) {
            y <- apply(coredata(x),2,FUN=FUN,na.rm=TRUE)
          } else if (FUN=="trend") {
            if (verbose) print("trend")
            y <- apply(x,2,FUN=FUN,na.omit=FALSE)
          } else {
            if (verbose) print('other')
            y <- apply(coredata(x),2,FUN=FUN) ## ,na.rm=TRUE)
          }
        }
        y <- attrcp(x,y)
        class(y) <- class(x)
        ## AM 30-06-2015 ...
        if (is.logical(colbar)) {
          ## If colbar set to FALSE, treat it as set to NULL
          if (!colbar) colbar <- NULL else
            colbar <- list(pal='t2m',rev=FALSE,n=10,
                           breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.1,show=TRUE)
        }
        
        ##                       
        ##if (!is.null(colbar)) {
        colbar <- colbar.ini(y,FUN=FUN,colbar=colbar,verbose=verbose)
      }
      if (verbose) print("length(col) =",length(colbar$col))
      
      ## Range of scale
      y.rng <- range(y,na.rm=TRUE)
      if (verbose) print(paste("range of mapped values",paste(y.rng,
                                                              collapse="/")))          
      ## find color index in colbar
      icol <- apply(as.matrix(y),2,findInterval,colbar$breaks)
      
      ## if (is.null(col)) col <- colbar$col[icol]
      ## if (is.null(bg)) bg <- colbar$col[icol]           
      if (verbose) print(range(y,na.rm=TRUE))
    }
    
    if (!is.null(FUN)) col <- NULL
    if (is.null(FUN)) bg.all <- bg
    
    ##scale <- apply(y,2,function(x) sum(!is.na(x))/length(x))
    if (!is.null(attr(x,'na')) & (!inherits(x,"stationmeta"))) {
      scale <- attr(x,'na')
    } else {
      scale <- 1
    }
    
    ##points(ss$longitude, ss$latitude, pch = pch, bg=bg , col=col,
    ##       cex = cex*scale, xlab = "", ylab = "", xlim = xlim, ylim = ylim,...)       
    
    if(is.null(colbar$pos)) pos <- 0.05
    
    ##fig0 <- par0$fig
    if (is.null(colbar$show)) colbar$show <- TRUE ## quick fix REB
    if (!is.null(FUN) & (!is.null(colbar)) & colbar$show) {
      if (showaxis) {
        fig0[3] <- par0$fig[3] + colbar$pos
        ## (par0$fig[4]-par0$fig[3])/150 ##0.075
      } else {
        fig0[3] <- par0$fig[3] + colbar$pos
        ## (par0$fig[4]-par0$fig[3])/120 ##0.05
      }
    } else {
      fig0 <- par0$fig
    }
    #par(fig=fig0)
    
    if(new) dev.new()
    par(fig=fig0,mar=mar0)
    
    if (!is.null(highlight)) {
      plot(highlight$longitude, highlight$latitude, pch = pch, col = col,
           bg = bg.all, cex = cex*scale, xlab = "", ylab = "",
           xlim = xlim, ylim = ylim , axes =FALSE , frame.plot = FALSE,
           cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab, new=FALSE)
    } else if (!is.null(ss) & !is.null(FUN)) {
      plot(ss$longitude, ss$latitude, pch = pch, col = "white",
           bg = "white", cex = cex * scale, xlab = "", ylab = "",
           xlim = xlim, ylim = ylim , axes = FALSE ,
           frame.plot = FALSE,
           cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab, new=FALSE)
    } else {
      plot(ss$longitude, ss$latitude, pch = pch, col = col, bg = bg,
           cex = cex*scale, xlab = "", ylab = "", xlim = xlim,
           ylim = ylim , axes = FALSE , frame.plot = FALSE,
           cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab, new=FALSE)
    }
    par(new=FALSE)
    
    ## Add geoborders
    data("geoborders", envir = environment())
    lines(geoborders$x, geoborders$y, col = "black")
    lines(attr(geoborders, "borders")$x, attr(geoborders, "borders")$y,
          col = "pink")##"grey90"
    
    if (showall) {
      ss.all <- select.station(param=is$param)
      points(ss.all$longitude,ss.all$latitude,pch=".",col="grey50",
             bg="grey",cex=cex/2)
    }
    
    if (!is.null(highlight)) {
      points(highlight$longitude, highlight$latitude, pch = pch, col = col,
             bg = bg.all, cex = cex*scale, xlab = "", ylab = "",
             xlim = xlim, ylim = ylim, cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab)
    } else if (!is.null(ss) & !is.null(FUN)) {
      points(ss$longitude, ss$latitude, pch = pch, col = "white",
             bg = "white", cex = cex*scale, xlab = "", ylab = "",
             xlim = xlim, ylim = ylim, cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab)
    } else {
      points(ss$longitude, ss$latitude, pch = pch, col = col, bg = bg,
             cex = cex*scale, xlab = "", ylab = "", xlim = xlim,
             ylim = ylim , cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab)
    }
    ## par(fig=par0$fig)
    ## print(par()$fig)
    ## add search info to plot
    if (text) {
      if (!is.null(highlight)) {
        title(main=paste("SOURCE(S): ",
                         paste(levels(factor(highlight$source)),
                               collapse="/" )), line=3,cex.main=cex.main)
        title(main=paste(length(levels(factor(highlight$location)))),
              line=2,cex.main=cex.main)
        title(main=paste(min(highlight$start,na.rm=TRUE),"/",
                         max(highlight$end,na.rm=TRUE)),
              line=2,cex.main=cex.main,adj=1)
        if (!is.null(FUN)) {
          title(main=paste(paste(toupper(apply(as.matrix(levels(factor(highlight$variable))),
                                               1,esd2ele)),collapse="/"),toupper(FUN),sep="/"),
                line=2,cex.main=cex.main , adj = 0)
        } else {
          title(main=paste(toupper(apply(as.matrix(levels(factor(highlight$variable))),
                                         1,esd2ele)),collapse="/"),line=2,cex.main=cex.main , adj = 0)
        }
      } else {
        title(main=paste("SOURCE(S) : ",
                         paste(levels(factor(ss$source)),collapse="/" )),
              line=3,cex.main=cex.main)
        title(main=paste(length(ss$location)),line=2,cex.main=.8)
        title(main=paste(max(ss$start,na.rm=TRUE),"/",
                         min(ss$end,na.rm=TRUE)),line=2,cex.main=.8,adj=1)
        title(main=paste(paste(toupper(levels(factor(ss$variable))),
                               collapse="/"),toupper(FUN),sep="/"),
              line=2,cex.main=cex.main , adj = 0)
      }
    }
    ## title(main=attr(z,"title"),line=2.2,cex.main=0.7)
    ## add margin text
    if (text) mtext(paste(("ESD package - map.station() - MET Norway 2014"),
                          "(www.met.no)",sep=" "),side=1,line=4,cex=0.6)
    
    par1 <- par()    
    if (!is.null(FUN)) {
      ##if (is.null(col)) colbar$col <- rep(col,length(colbar$col[icol]))
      ## 
      if (!is.null(col)) col <- col else col <- colbar$col[icol]
      points(ss$longitude, ss$latitude, pch = pch,
             bg=colbar$col[icol], col=col, ##col=colbar$col[icol]
             cex = cex*scale, xlab = "", ylab = "",
             xlim = xlim, ylim = ylim,...)
      
      par(fig=fig0,new=TRUE)
      
      ## print(par()$fig)
      
      if (!is.null(highlight)) {
        points(highlight$longitude, highlight$latitude, pch = 21 ,
               col = col.subset,
               bg=bg.subset, cex = cex.subset,...)
        
      }
      
      ## Add geoborders
      data("geoborders", envir = environment())
      lines(geoborders$x, geoborders$y, col = "grey50")
      lines(attr(geoborders, "borders")$x,
            attr(geoborders, "borders")$y, col = "pink") ##"grey90"
      
      if (show.val) text(ss$longitude,ss$latitude,round(y,digits=2),
                         cex=cex/4)
      
      par(fig=fig0,new=TRUE)
      ## print(par()$fig)
      ## add color bar
      if (colbar$show) {
        if (fancy & !is.null(colbar))
          col.bar(colbar$breaks,horiz=TRUE,pch=21,v=1,h=1,
                  col=colbar$col,cex=2,cex.lab=colbar$cex.lab,
                  type=colbar$type,verbose=FALSE,vl=1,border=FALSE)
      } else if (!is.null(colbar)) {
        ##fig1 <- par0$fig
        par(fig=par0$fig,new=TRUE)
        image.plot(lab.breaks=colbar$breaks,horizontal = TRUE,
                   legend.only = T, zlim = range(colbar$breaks),
                   col = colbar$col, legend.width = 1,
                   axis.args = list(cex.axis = cex.axis,
                                    xaxp=c(range(colbar$breaks),n=colbar$n)),
                   border = FALSE,
                   legend.shrink=legend.shrink)
        ##bigplot=c(0,0,1,1),
        ##smallplot=c(0,0,1,1))
      }
    }    
    
    ##par(fig=fig0,new=TRUE) ## AM 18-06-2015 comment
    if (verbose) {
      print(paste('colbar$breaks are ',
                  paste(colbar$breaks,collapse="/"),
                  'of length',length(colbar$breaks)))
      print(paste('colbar$col are ',
                  paste(colbar$col,collapse="/"),'of length',
                  length(colbar$col)))      
      print(paste('colbar$n ', paste(colbar$n,collapse="/")))
    }
    
    ## add text if TRUE
    if (!is.null(unlist(is))) {
      if (add.text.subset & full.names.subset) {
        text(highlight$longitude, highlight$latitude, highlight$location,
             pos=3,cex=cex.subset/2)
      }
    } else if(!is.null(highlight)) {
      text(highlight$longitude, highlight$latitude,
           substr(toupper(highlight$location),1,3),pos=3,cex=cex.subset/2)
    }
    if (add.text) {
      if (full.names) {
        text(ss$longitude, ss$latitude,ss$location,pos=3,cex=cex/2)
      } else {
        text(ss$longitude, ss$latitude,substr(toupper(ss$location),1,3),pos=3,cex=cex/2)
      }
    } 
    
    if (showaxis) title(xlab="Longitude", line=2.2, cex.lab=cex.lab, col="grey30") 
    if (showaxis) title(ylab="Latitude", line=2.2, cex.lab=cex.lab, col="grey30") 
    
    ## format axes
    if (showaxis) axis(1,pretty(seq(xlim[1],xlim[2],by=5),n=5),
                       cex.axis=cex.axis,col="grey30",col.ticks="grey30") # 0.7
    if (showaxis) axis(2,pretty(seq(ylim[1],ylim[2],by=5),n=5),
                       cex.axis=cex.axis,col="grey30",col.ticks="grey30")
    if (showaxis) axis(4,pretty(seq(ylim[1],ylim[2],by=5),n=5),
                       cex.axis=cex.axis,col="grey30",col.ticks="grey30")
    if (showaxis) axis(3,pretty(seq(xlim[1],xlim[2],by=5),n=5),
                       cex.axis=cex.axis,col="grey30",col.ticks="grey30")
    ## add grid
    par(fig=par0$fig,new=TRUE)
    if (gridlines) grid()
    
    ## lines(geoborders$x, geoborders$y, col = "black")
    ## lines(attr(geoborders, "borders")$x, attr(geoborders, "borders")$y, col = "grey90")
  }
  ## return(par1)
}



sphere <- function(x,n=30,FUN="mean",lonR=10,latR=45,axiR=0,xlim=NULL,ylim=NULL,
                   gridlines=TRUE,col="green",bg="darkgreen",cex=0.2,
                   cex.axis=1,cex.lab=1,cex.main=1.5,pch=".",new=TRUE,verbose=FALSE,...) {
  if(verbose) print("sphere")
  x0 <- x
  
  ## KMP 2016-12-21: To handle xlim in greenwich format, e.g., 180-360
  if(!is.null(xlim)) {
    greenwich <- (min(xlim)>0 & max(xlim)>180)
    g2dl(x,greenwich=greenwich,verbose=verbose)
  } else greenwich <- NULL
  
  ## Data to be plotted:
  if (inherits(x,"stationmeta")) {
    lon <- x$longitude
    lat <- x$latitude
    param <- param2ele(x$ele)
    unit <- " "  
  } else if (inherits(x,"station")) {
    lon <- attr(x,'longitude')
    lat <- attr(x,'latitude')
    param <- as.character(levels(factor(attr(x,'parameter'))))
  }
  ## KMP 2016-12-21: To handle xlim in greenwich format, e.g., 180-360
  if(is.null(greenwich)) greenwich <- (min(lon)>0 & max(lon)>180)
  
  ## To deal with grid-conventions going from north-to-south or east-to-west:
  ##srtx <- order(attr(x,'longitude')); lon <- lon[srtx]
  ##srty <- order(attr(x,'latitude')); lat <- lat[srty]
  
  if (!is.null(FUN)) {
    map <- apply(as.matrix(x),2,FUN,na.rm=TRUE) ##map <- x[srtx,srty]
  } else {
    map <- x
  }
  
  # Rotatio:
  # longitudinal rotation
  if(!is.null(xlim)) {
    lonR <- mean(xlim, na.rm=TRUE)
  } else if (is.null(lonR)) {
    lonR <- mean(lon, na.rm=TRUE)  
  }
  # latitudinal rotation
  if(!is.null(ylim)) {
    latR <- mean(ylim, na.rm=TRUE) 
  } else if (is.null(latR)) {
    latR <- mean(lat, na.rm=TRUE)
  }
  # axiR: rotation of Earth's axis
  
  # coastline data:
  data("geoborders",envir=environment())
  ## KMP 10-11-2015: apply xlim and ylim
  gx <- geoborders$x
  gy <- geoborders$y
  ok <- is.finite(gx) & is.finite(gy)
  if(greenwich) gx[gx<0 & ok] <- gx[gx<0 & ok] + 360
  ## KMP 2017-09-19: New xlim/ylim code line 829
  #if (!is.null(xlim)) ok <- ok & gx>=min(xlim) & gx<=max(xlim)
  #if (!is.null(ylim)) ok <- ok & gy>=min(ylim) & gy<=max(ylim)
  theta <- pi*gx[ok]/180
  phi <- pi*gy[ok]/180
  #ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  #theta <- pi*geoborders$x[ok]/180; phi <- pi*geoborders$y[ok]/180
  x <- sin(theta)*cos(phi)
  y <- cos(theta)*cos(phi)
  z <- sin(phi)
  
  # Calculate contour lines if requested...  
  #contourLines
  ##lonxy <- rep(lon,length(lat))
  ##latxy <- sort(rep(lat,length(lon)))
  ##map<- c(map)
  
  # Remove grid boxes with missign data:
  ##ok <- is.finite(map)
  #print(paste(sum(ok)," valid grid point"))
  ##lonxy <- lonxy[ok]; latxy <- latxy[ok]; map <- map[ok]
  
  # Define the grid box boundaries:
  ##dlon <- min(abs(diff(lon))); dlat <- min(abs(diff(lat)))
  ##Lon <- rbind(lonxy - 0.5*dlon,lonxy + 0.5*dlon,
  ##             lonxy + 0.5*dlon,lonxy - 0.5*dlon)
  ##Lat <- rbind(latxy - 0.5*dlat,latxy - 0.5*dlat,
  ##             latxy + 0.5*dlat,latxy + 0.5*dlat)
  Theta <- pi*lon/180; Phi <- pi*lat/180
  
  # Transform -> (X,Y,Z):
  X <- sin(Theta)*cos(Phi)
  Y <- cos(Theta)*cos(Phi)
  Z <- sin(Phi)
  #print(c( min(x),max(x)))
  
  ## Define colour pal:
  if (!is.null(FUN)) {
    if (is.null(col)) col <- colscal(n=n,pal=varid(x)) else
      if (length(col)==1) {
        col <- colscal(pal=col,n=n)
      }
    nc <- length(col)
    index <- round( nc*( map - min(map) )/
                      ( max(map) - min(map) ) )
  } 
  # Rotate coastlines:
  a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
  a <- rotM(x=latR,y=0,z=0) %*% a
  x <- a[1,]; y <- a[2,]; z <- a[3,]
  
  # Grid coordinates:
  #d <- dim(X)
  #print(d)
  
  # Rotate data grid:  
  A <- rotM(x=0,y=0,z=lonR) %*% rbind(c(X),c(Y),c(Z))
  A <- rotM(x=latR,y=0,z=0) %*% A
  X <- A[1,]; Y <- A[2,]; Z <- A[3,]
  #dim(X) <- d; dim(Y) <- d; dim(Z) <- d
  #print(dim(rbind(X,Z)))
  
  # Rotate xlim and ylim
  if(!is.null(xlim) & !is.null(ylim)) {
    thetalim <- pi*xlim/180
    philim <- pi*ylim/180
    Xlim <- sin(thetalim)*cos(philim)
    Ylim <- cos(thetalim)*cos(philim)
    Zlim <- sin(philim)
    Alim <- rotM(x=0,y=0,z=lonR) %*% rbind(c(Xlim),c(Ylim),c(Zlim))
    Alim <- rotM(x=latR,y=0,z=0) %*% Alim
    Xlim <- Alim[1,]; Ylim <- Alim[2,]; Zlim <- Alim[3,]
  } else {
    Xlim <- range(x, na.rm=TRUE)
    Zlim <- range(z, na.rm=TRUE)
  }
  
  # Plot the results:
  if(new) dev.new()
  par(bty="n",xaxt="n",yaxt="n",new=TRUE)
  plot(Xlim,Zlim,pch=".",col="white",xlab="",ylab="",
       cex.axis=cex.axis,cex.lab=cex.lab)
  #plot(x,z,pch=".",col="white",xlab="",ylab="",
  #     cex.axis=cex.axis,cex.lab=cex.lab)
  #par0 <- par()
  
  # plot the grid boxes, but only the gridboxes facing the view point:
  ##Visible <- Y > 0 ##colMeans(Y) > 0
  ##X <- X[,Visible]; Y <- Y[,Visible]; Z <- Z[,Visible]
  ##index <- index[Visible]
  ##apply(rbind(X,Z,index),2,gridbox,cols)
  # c(W,E,S,N, colour)
  # xleft, ybottom, xright, ytop
  
  ##if (!is.null(FUN)) {
  ##    colb <- colscal(n=length(breaks)) 
  ##    col <- colb[findInterval(map,breaks)]
  ##    bg <- col
  ##    nc <- length(colb)
  ##}
  
  ## Initialise colbar
  colbar <- colbar.ini(map,FUN=FUN)
  breaks <- colbar$breaks
  colb <- colbar$col 
  col <- colb[findInterval(map,breaks)]
  bg <- col
  nc <- length(colb)
  visible <- Y > 0
  points(X[visible],Z[visible],cex=cex,pch=pch,col=col,bg=bg)
  
  ## Add contour lines?
  ## Plot the coast lines  
  visible <- y > 0
  points(x[visible],z[visible],pch=".")
  #plot(x[visible],y[visible],type="l",xlab="",ylab="")
  lines(cos(pi/180*1:360),sin(pi/180*1:360),col="black")
  
  ## Add grid ?
  if (!is.null(FUN)) {    
    ## Colourbar:  
    par(fig = c(0.3, 0.7, 0.05, 0.10),mar=rep(0,4),cex=0.8,
        new = TRUE, mar=c(1,0,0,0), xaxt = "s",yaxt = "n",bty = "n")
    #print("colourbar")
    ##breaks <- round( nc*(seq(min(map),max(map),length=nc)- min(map) )/                 ( max(map) - min(map) ) )
    bar <- cbind(breaks,breaks)
    image(seq(breaks[1],breaks[length(breaks)],length=nc),
          c(1,2),bar,col=col,cex.axis=cex.axis)
    
    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
        fig=c(0,1,0,1),new=TRUE)
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    text(0.1,0.95,param,cex=cex.main,pos=4)
    text(0.72,0.002,unit,pos=4)  
  }
  ##result <- data.frame(x=colMeans(Y),y=colMeans(Z),z=c(map))
  if (inherits(x0,"stationmeta")) {
    result <- data.frame(x=Y,y=Z)
  } else if (inherits(x0,"station")) {
    result <- data.frame(x=Y,y=Z,z=map)
  }
  invisible(result)
}

## Perform a series of tests that produce and save maps of the different data sources into the local directory.
# do not export
test.map.station <- function(save=FALSE) {
  
  map.station(src="NACD",col="darkgreen",bg="green")
  if (save) dev.copy2pdf(file="NACD-Network.pdf")
  if (save) dev.off()
  map.station(src="ECAD",col="darkgreen",bg="green")
  if (save) dev.copy2pdf(file="ECAD-Network.pdf")
  if (save) dev.off()
  map.station(src="GHCNM",col="darkgreen",bg="green")
  if (save) dev.copy2pdf(file="GHCNM-Network.pdf")
  if (save) dev.off()
  map.station(src="GHCND",col="darkgreen",bg="green")
  if (save) dev.copy2pdf(file="GHCND-Network.pdf")
  if (save) dev.off()
  map.station(src="NORDKLIM",col="darkgreen",bg="green")
  if (save) dev.copy2pdf(file="NORDKLIM-Network.pdf")
  if (save) dev.off()
  map.station(src="NARP",col="darkgreen",bg="green")
  if (save) dev.copy2pdf(file="NARP-Network.pdf")
  if (save) dev.off()
  
  locs <- map.station(ele=101)
  locs <- map.station(cntr="NORWAY",ele=101,src="ECAD")
}


## The main function to produce map of subset of stations
#' @export map.stationmeta
map.stationmeta <- function(x,...)
  map.stationsummary(x,...)

#' @export map.data.frame
map.data.frame <- function(x,...) {
  
  att <- c("station_id","location","country","longitude","latitude","altitude","element","start","end","source","wmo","quality")
  
  if (sum(is.element(names(x),att))==12) {   
    class(x) <- c("stationmeta","data.frame")
    map.station(x,...)
  } else if (inherits(x,c('stationsummary','stationmeta'))) {
    map.stationsummary(x,...)
  } else print("x is not a stationmeta object")
}

map.stationsummary <- function(x,FUN=NULL,cex=1,cex0=1,col='red',pal='t2m',pch=19,nbins=15,rev=FALSE,
                               new=TRUE,verbose=FALSE,fig=c(0.2,0.25,0.6,0.8),
                               hist=TRUE,lon=NULL,lat=NULL,...) {
  if (verbose) print(match.call())
  if (new) dev.new()
  ok <- rep(TRUE,length(x$longitude)); ok2 <- ok
  if (!is.null(lon)) {
    keep <- (x$longitude >= min(lon)) & (x$longitude <= max(lon))
    x <- x[keep,]
  }
  if (!is.null(lat)) {
    keep <- (x$latitude >= min(lat)) & (x$latitude <= max(lat))
    x <- x[keep,]
  }
  if (!is.null(FUN)) {
    if (verbose) {print(paste('FUN=',FUN)); print(names(x))}
    ## If FUN specified, change the colours
    if (length(grep(FUN,names(x)[is.element(nchar(names(x)),nchar(FUN))]))==1) {
      z <- x[[FUN]] 
      ok <- is.finite(z); ok2 <- ok
      if (verbose) print(summary(z))
      colbar <- colscal(n=nbins,pal=pal,rev=rev)
      breaks <- pretty(z,nbins)
      #ic <- trunc(nbins*(z - min(z,na.rm=TRUE))/(max(z,na.rm=TRUE) - min(z,na.rm=TRUE))) + 1
      #breaks <- round(seq(min(z,na.rm=TRUE),max(z,na.rm=TRUE),length=nbins),2)
      #ic[ic > nbins] <- nbins
      ic <- rep(NA,length(z))
      if (verbose) print(summary(z))
      for (i in 1:length(z)) ic[i] <- sum(breaks <= z[i],na.rm=TRUE) 
      ok <- ok & (ic > 0) & (ic <= length(breaks))
      if (verbose) {print(c(NA,breaks)); print(table(ic))}
      col <- colbar[ic[ok]]
    } else if (verbose) print('No match')
  }
  if (is.character(cex)) {
    if (verbose) print(paste('cex=',cex))
    ## If FUN specified, change the colours
    if (length(grep(cex,names(x)))==1) {
      z2 <- x[[cex]]  
      ok2 <- is.finite(z2)
      if (is.numeric(z2)) {
        cex <- cex0*sqrt( (z2 - min(z2,na.rm=TRUE))/(max(z2,na.rm=TRUE) - min(z2,na.rm=TRUE)) )
      }
      cex <- cex[ok2]
      if (verbose) print(summary(cex))
    }
  }
  
  par(bty='n',mar=c(3,1,3,2),xaxt='n',yaxt='n')
  if (!is.null(FUN)) {  
    nf <- layout(matrix(c(rep(1,56),0,0,rep(2,4),0,0), 8, 8, byrow = TRUE), respect = TRUE)
    nd <- 2
    if (max(abs(z),na.rm=TRUE) > 100) nd <- 0 else if (max(abs(z),na.rm=TRUE) >= 10) nd <- 1
    main=paste0(FUN,' (mean= ',round(mean(z,na.rm=TRUE),nd),', sd=',round(sd(z,na.rm=TRUE),nd),
                ' [',round(min(z,na.rm=TRUE),nd),', ',round(max(z,na.rm=TRUE),nd),'])')
  } else main <- ''
  plot(x$longitude[ok&ok2],x$latitude[ok&ok2],col=col,cex=cex,pch=pch,xlab='',ylab='',
       main=main,...)
  data("geoborders",envir = environment())
  grid()
  par(yaxt='s',xaxt='s')
  axis(1,at=10*round(x$longitude/10),col='grey',col.axis='grey',col.lab='grey',col.ticks='grey')
  axis(2,at=10*round(x$latitude/10),col='grey',col.axis='grey',col.lab='grey',col.ticks='grey')
  lines(geoborders,col='grey')
  lines(attr(geoborders,'borders'),col='lightgreen')
  points(x$longitude[ok&ok2],x$latitude[ok&ok2],col=col,cex=cex,pch=pch)
  if (!is.null(FUN)) {
    par0 <- par()
    par(mar=c(3,0,1,0),cex.lab=0.5,yaxt='n',xaxt='s')
    image(breaks, 1:2, cbind(breaks, breaks), col = colbar, cex.axis = 0.75)
    par(par0)
    
    if (hist) {
      par(new=TRUE,fig=c(0.05,0.30,0.03,0.32))
      hist(z,col='grey',main='')
      par(new=FALSE,fig=c(0,1,0,1))
    }
  }
}
