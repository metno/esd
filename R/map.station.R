## Author 	 Rasmus E. Bnestad
## Updated 	 by Abdelkader Mezghani and Kajsa Parding
## Rasmus. E. Benestad - attempt to simpify by splitting up
## Last update   27.07.2017
## Includes	 map.station() ; test.map.station()
## Require 	 geoborders.rda

genfun <- function(x,FUN) {
  if (sum(is.element(names(attributes(x)),FUN))>0){
    ## REB 2015-12-17: Use FUN to colour the symbols according to some attribute:
    FUN <- eval(parse(text=paste("attr(x,'",FUN,"')")))
  } else if (sum(is.element(names(x),FUN))>0){
    ## REB 2015-12-17: Use FUN to colour the symbols according to some list element (stationmeta-objects):
    if (verbose) print('FUN refers to a list element')
    FUN <- eval(parse(text=paste("function(x,...) x$",FUN,sep='')))
    return(FUN)
  }}

## Simplified function for mapping station objects.
map.station <- function (x=NULL,FUN=NULL, it=NULL,is=NULL,new=FALSE,
                         add=FALSE,projection="lonlat",
                         xlim = NULL, ylim = NULL,zlim=NULL,n=15,
                         col='darkred',bg='orange',
                         colbar= list(pal='t2m',col=NULL,rev=FALSE,n=10,
                                      breaks=NULL,type="p",cex=2,h=0.6, v=1,
                                      pos=0.1,show=TRUE),
                         # col=NULL replaced by palette
                         type=NULL,gridlines=TRUE,
                         lonR=NULL,latR=45,axiR=NULL,verbose=FALSE,
                         cex=2,zexpr="alt",cex.subset=1,
                         add.text.subset=FALSE,showall=FALSE,
                         add.text=FALSE,main=NULL,sub=NULL,
                         height=NULL,width=NULL,
                         cex.main=1,cex.sub=0.75,cex.axis=1,cex.lab=0.6,
                         col.main="black",col.sub="grey",
                         font.main=1,font.sub=4,
                         pch=19, from=NULL,to=NULL,showaxis=FALSE,
                         border=FALSE,full.names=FALSE,
                         full.names.subset=FALSE, use.old=FALSE,
                         text=FALSE, fancy=FALSE, 
                         na.rm=TRUE,show.val=FALSE,usegooglemap=FALSE,
                         ##colorbar=TRUE,
                         xlab="lon",ylab="lat",
                         legend.shrink=1,fig=c(0,1,0.05,0.95),
                         mar=rep(2,4),mgp=c(3,1,0),...) { 
  if ( (inherits(x,"stationmeta")) | (projection != 'lonlat') | usegooglemap | use.old) {
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
                    usegooglemap=usegooglemap,legend.shrink=legend.shrink,...)
  } else {
    if (verbose) print('map.station - new version')
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
        y <- apply(coredata(x),2,FUN,na.rm=na.rm)
      } else {
        y <- attr(x,FUN); FUN <- NULL
      }    
      if (verbose) print(summary(y))
      colbar <- colbar.ini(y,colbar=colbar)
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
      show.colbar <- TRUE
    } else {
      y <- rep(1,length(lon(x)))
      show.colbar <- FALSE
    }
    
    ## KMP 2017-07-28: fig creates problems when you want to add map.station as a subplot.
    ## With this solution you have to use add=TRUE and set fig to your subplot or to NULL.
    if(is.null(fig)) fig <- par(fig)
    if(add) {
      par(fig=fig,mar=mar,mgp=mgp,new=TRUE,bty='n',xaxt='n',yaxt='n',cex.axis=0.7,
          col.axis='grey30',col.lab='grey30',las=1)
    } else {
      par(fig=fig,mar=mar,mgp=mgp,new=FALSE,bty='n',xaxt='n',yaxt='n',cex.axis=0.7,
          col.axis='grey30',col.lab='grey30',las=1)
    }
    
    ## Avoid errors when plotting the colorbar with small figure windows
    ## fin collects information about the figure size (in inches)
    fin <- par()$fin
    
    ## For checking & debugging
    if (verbose) {
      print(paste('window size=',fin,collapse=' '))
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
        ylim <- range(is$lat,na.rm=TRUE) + c(-1,1)
      } else {
        ylim <- range(lat(x),na.rm=TRUE) + c(-2,2)
      }
    }
    
    plot(lon(x),lat(x),xlim=xlim,ylim=ylim,col=col,pch=pch,cex=2,new=FALSE,
         xlab=xlab,ylab=ylab)
    if (add.text) text(lon(x),lat(x),substr(loc(x),1,6),cex=0.6,col='grey',pos=1)
    
    if(showaxis | gridlines) {
      par(xaxt="s",yaxt="s",las=1,col.axis='grey',col.lab='grey',
          cex.lab=0.9,cex.axis=0.9)
      axis(3,seq(floor(par("xaxp")[1]/5)*5,par("xaxp")[2],by=5),col='grey')
      axis(4,seq(floor(par("yaxp")[1]/5)*5,par("yaxp")[2],by=5),col='grey')
      if (gridlines) grid()
    }

    data("geoborders")
    lines(geoborders$x,geoborders$y)
    if (border) lines(attr(geoborders,'border')$x,attr(geoborders,'border')$y,col='grey')
    
    if (show.colbar) {
      ## KMP 2017-07-28: If fig is something other than the default
      ## the colbar may be misplaced relative to the plot.
      if (fin[2] >= 8) {
        cf <- c(0.2,0.8,0,0.1)*fig
        par(new=TRUE,fig=cf,mar=rep(1.5,4),yaxt='n')
        #par(new=TRUE,fig=c(0.2,0.8,0,0.1),mar=rep(1.5,4),yaxt='n') 
      } else {
        cf <- c(0.2,0.8,0,0.15)*fig
        par(new=TRUE,fig=cf,mar=rep(2,4),yaxt='n')
        #par(new=TRUE,fig=c(0.2,0.8,0,0.15),mar=rep(2,4),yaxt='n')
      }
      image(colbar$breaks,1:2,cbind(colbar$breaks,colbar$breaks),col=colbar$col,axes=FALSE)
      par(mar=c(2,1,2,1),mgp=c(2,0.4,0),cex.axis=0.7,col.axis='grey')
      axis(1,colbar$breaks)
      
    } 
    par(new=TRUE,fig=fig,mar=mar,yaxt='n',xaxt='n')
    plot(lon(x),lat(x),type='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
    
    ## Add a title
    if(is.null(main)) {
      main <- paste(attr(x,"param")[1],"  (",attr(x,"unit")[1],")",sep="")
      if(!is.null(FUN)) {
        if(is.function(FUN)) {
          main <- paste(as.character(quote(FUN)),"of",main)
        } else {
          main <- paste(as.character(FUN),"of",main)
        }
      }
    }
    title(main=main,sub=sub,line=-2,adj=0,cex.main=cex.main,cex.sub=cex.sub,
          col.main=col.main,col.sub=col.sub,font.main=font.main,font.sub=font.sub)
  }
  if (verbose) print('Organise output')
  dim(y) <- c(1,length(y))
  y <- zoo(y,order.by=1)
  if (verbose) print(dim(y))
  class(y) <- class(x)
  y <- attrcp(x,y)
  attr(y,'period') <- paste(range(index(x)))
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

###

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
                             text=FALSE, fancy=FALSE, 
                             na.rm=TRUE,show.val=FALSE,usegooglemap=FALSE,
                             ##colorbar=TRUE,
                             legend.shrink=1,...) { 
  ##
  if (verbose) {
    print(paste('map.station',FUN))
    print(class(x))
  }
  arg <- list(...)
  attr(x,'unit') <- as.character(unit(x))
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
           xlab=xlab,ylab=ylab...)
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
      print('The station metadata')
      str(ss)
    }  
    
    if (is.null(attr(ss,"element"))) ss$element <-apply(as.matrix(ss$variable),1,esd2ele)   
    
    if (verbose) str(ss)
    
    if (!is.null(unlist(is))) { ## highlight a subset of station
      
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
          xlim <- range(ss$longitude, na.rm = TRUE) + c(-1,1)
        } else {
          xlim <- range(ss$longitude, na.rm = TRUE) + c(-4,4)
        }
      } else {
        if (length(is$lon) > 1) {
          xlim <- range(highlight$longitude, na.rm = TRUE) + c(-1,1)
        } else {
          xlim <- range(highlight$longitude, na.rm = TRUE) + c(-4,4)
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
    par(fig=par0$fig,mar=mar0)
    
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
    par(fig=fig0)
    
    ## REB: 2016-10-12 - add the possibility to use google maps
    if (("RgoogleMaps" %in% rownames(installed.packages()) == TRUE) &
         (projection=="lonlat") & usegooglemap) {
      require(RgoogleMaps)
      mxdst <- max(diff(range(ss$latitude)),diff(range(ss$longitude)))
      if (!is.finite(mxdst) | mxdst==0) {
        zoom <- 3 
      } else {
        zoom <- 7 - round(log(mxdst))
      }
      bgmap <- GetMap(center=c(lat=mean(ss$latitude),lon=mean(ss$longitude)),
                      destfile = "map.station.esd.png",
                      maptype = "mobile", zoom=zoom)
      plotmap(ss$latitude, ss$longitude, bgmap)
      
      print('Unfinished')
      return()
    }

    if (!is.null(highlight)) {
      plot(highlight$longitude, highlight$latitude, pch = pch, col = col,
           bg = bg.all, cex = cex*scale, xlab = "", ylab = "",
           xlim = xlim, ylim = ylim , axes =FALSE , frame.plot = FALSE,
           cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab)
    } else if (!is.null(ss) & !is.null(FUN)) {
      plot(ss$longitude, ss$latitude, pch = pch, col = "white",
           bg = "white", cex = cex * scale, xlab = "", ylab = "",
           xlim = xlim, ylim = ylim , axes = FALSE ,
           frame.plot = FALSE,
           cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab)
    } else {
      plot(ss$longitude, ss$latitude, pch = pch, col = col, bg = bg,
           cex = cex*scale, xlab = "", ylab = "", xlim = xlim,
           ylim = ylim , axes = FALSE , frame.plot = FALSE,
           cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab)
    }
    #if ( ("RgoogleMaps" %in% rownames(installed.packages()) == TRUE) )
    #     par(new=FALSE) else ## REB: 2016-10-12 - add the possibility to use google maps
    #     par(new=TRUE)
    par(new=FALSE)
    
    ## Add geoborders
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
          title(main=paste(paste(toupper(apply(as.matrix(levels(,
                                                                factor(highlight$variable))),
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
    if (is.null(col)) col <- colscal(n=n,col=varid(x)) else
      if (length(col)==1) {
        pal <- col
        col <- colscal(pal=pal,n=n)
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


## The main function to produce map of subseted stations
map.stationmeta <- function(...)
  map.station(...)

map.data.frame <- function(x,...) {
  
  att <- c("station_id","location","country","longitude","latitude","altitude","element","start","end","source","wmo","quality")
  
  if (sum(is.element(names(x),att))==12) {   
    class(x) <- c("stationmeta","data.frame")
    map.station(x,...)
  }
  else print("x is not a stationmeta object")
}
