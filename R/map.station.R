## Author 	 Rasmus E. Bnestad
## Updated 	 by Abdelkader Mezghani and Kajsa Parding
## Rasmus. E. Benestad - attempt to simplify by splitting up
## Last update   08.10.2025
## Includes	 map.station() ; test.map.station()
## Require 	 geoborders.rda

#genfun <- function(x,FUN,verbose=FALSE) {
#  if (sum(is.element(names(attributes(x)),FUN))>0){
#    ## REB 2015-12-17: Use FUN to colour the symbols according to some attribute:
#    FUN <- eval(parse(text=paste("attr(x,'",FUN,"')")))
#  } else if (sum(is.element(names(x),FUN))>0){
#    ## REB 2015-12-17: Use FUN to colour the symbols according to some list element (stationmeta-objects):
#    if (verbose) print('FUN refers to a list element')
#    FUN <- eval(parse(text=paste("function(x,...) x$",FUN,sep='')))
#    return(FUN)
#  }
#}
## KMP 2025-10-08: Updating genfun so that it works with character input such as
## "mean" or "sd", transforming them to functions using get(FUN)
genfun <- function(x,FUN,verbose=FALSE) {
  if(!is.null(FUN)) {
    if (sum(is.element(names(attributes(x)),FUN))>0){
      ## REB 2015-12-17: Use FUN to colour the symbols according to some attribute:
      FUN <- eval(parse(text=paste("attr(x,'",FUN,"')")))
    } else if (sum(is.element(names(x),FUN))>0){
      ## REB 2015-12-17: Use FUN to colour the symbols according to some list element (stationmeta-objects):
      if (verbose) print('FUN refers to a list element')
      FUN <- eval(parse(text=paste("function(x,...) x$",FUN,sep='')))
    } else if(exists(FUN) & mode(get(FUN))=="function") {
      FUN <- get(FUN)
    }
    return(FUN)
  }
}


## Simplified function for mapping station objects.
#' @exportS3Method
#' @export map.station
map.station <- function(x=NULL,FUN=NULL, it=NULL,is=NULL,new=FALSE,
                        add=FALSE,projection="lonlat",
                        xlim = NULL, ylim = NULL, zlim=NULL, n=15,
                        col='darkred',bg='orange',
                        colbar= list(pal='t2m',col=NULL,rev=FALSE,n=6,
                                     breaks=NULL,type="p",cex=2,h=0.6, v=1,
                                     pos=0.1,show=FALSE),
                        # col=NULL replaced by palette
                        type=NULL,gridlines=TRUE,
                        lonR=NULL,latR=45,axiR=NULL,verbose=FALSE,
                        cex=2,zexpr="alt",cex.subset=1,
                        add.text.subset=FALSE,showall=FALSE,
                        add.text=FALSE,n.text=5,main=NULL,sub=NULL,
                        height=NULL,width=NULL,
                        cex.main=1,cex.sub=0.75,cex.axis=1,cex.lab=0.9,
                        col.main="black",col.sub="grey",col.border="grey",
                        col.text='grey',
                        font.main=1,font.sub=4,
                        pch=19, from=NULL,to=NULL,showaxis=FALSE,
                        border=FALSE,full.names=FALSE,
                        full.names.subset=FALSE, use.old=FALSE,
                        text=FALSE, fancy=TRUE, 
                        na.rm=TRUE,show.val=FALSE,#usegooglemap=FALSE,
                        ##colorbar=TRUE,
                        add.significance=FALSE, pval=0.01, col.pval='black', lwd.pval=2,
                        pch.significance=21, pch.positive=24, pch.negative=25,
                        xlab="lon", ylab="lat",
                        legend.shrink=1,fig=NULL,#fig=c(0,1,0.05,0.95),
                        mar=rep(2,4), mgp=c(3,1,0), 
                        plot=TRUE,...) { 
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
      ## REB 2024-04-30
      ## KMP 2024-09-05: rev is already applied in colbar.ini. No need to do it again. You are just reversing it back. 
      #if (colbar$rev) {wr <- rev(wr); wg <- rev(wg); wb <- rev(wb); colbar$col <- rev(colbar$col)}
      col <- rep(colbar$col[1],length(y))
      for (i in 1:length(y)) {
        ii <- round(approx(0.5*(colbar$breaks[-1]+colbar$breaks[-length(colbar$breaks)]),
                           1:length(colbar$col),xout=y[i],rule=2)$y)
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
    
    if(is.character(FUN)) if(grepl("trend", FUN) & !is.null(pch.positive) & !is.null(pch.negative)) {
      pch <- rep(pch[1], length(y))
      pch[y>0] <- pch.positive
      pch[y<0] <- pch.negative
      pch.significance <- pch
    }
    
    par0 <- par()
    par(mar=mar,mgp=mgp,bty='n',xaxt='n',yaxt='n',cex.axis=0.7,
        xpd=FALSE,col.axis='grey30',col.lab='grey30',las=1)
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
      
      plot(lon(x),lat(x),xlim=xlim,ylim=ylim,col=col,bg=col,pch=pch,cex=cex,new=FALSE,
           cex.lab=cex.lab,xlab=xlab,ylab=ylab)
      if (add.text) text(lon(x),lat(x),substr(loc(x),1,n.text),cex=cex.lab,col=col.text,pos=1)
      
      if(showaxis | gridlines) {
        par(xaxt="s",yaxt="s",las=1,col.axis='grey',col.lab='grey',
            cex.lab=cex.lab,cex.axis=cex.axis)
        if(diff(range(xlim))>10) {
          dlon <- round(round(diff(range(xlim))/4)/5)*5
        } else if(diff(range(xlim))<2) {
          dlon <- round(diff(range(xlim))/3, 2)
        } else {
          dlon <- round(diff(range(xlim))/4)
        } 
        if(diff(range(ylim))>10) {
          dlat <- round(round(diff(range(ylim))/4)/5)*5
        } else if(diff(range(ylim))<2) {
          dlat <- round(diff(range(ylim))/3, 2)
        } else {
          dlat <- round(diff(range(ylim))/4)
        }
        if (dlon==0) dlon <- 1 ## fudge
        if (verbose) {print(xlim); print(ylim); print(paste('dlon=',dlon,'dlat=',dlat))}
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
          rect(below[1], below[2]-2*dy, below[3], below[4], 
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
                 col=col.pval[i], bg=NA, lwd=lwd.pval[i], cex=cex, 
                 pch=pch.significance[pval.x>=pmin & pval.x<p])
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
    ## KMP 2024-10-17: Do NOT return mar to its original value par0$mar if it has been
    ## changed within the function. For some reason this messes with the map so that
    ## additional elements (points, rectangles etc) cannot be added at the edges of the map!
    par(mgp=par0$mgp,bty=par0$bty,xaxt=par0$xaxt, #mar=par0$mar,
        yaxt=par0$yaxt,cex.axis=par0$cex.axis,xpd=par0$xpd,
        col.axis=par0$col.axis,col.lab=par0$col.lab,
        las=par0$las,xpd=par0$xpd)
    
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
    # rectangle used as a test that new elements can be added to the maps
    #rect(-10, 60, 30, 90, col=adjustcolor("red", alpha.f=0.1))
    #par2 <- par()
    #for(nm in names(par1)) if(!identical(par1[[nm]], par2[[nm]])) {
    #    if(verbose) print(paste0("par$", nm, " is different"))
    #    if(verbose) print(paste(paste(par1[[nm]], collapse=" "), 
    #                            paste(par2[[nm]], collapse=" ")))
    #}
    invisible(y)
  }
}

###
# Internal function - no need to export
map.station.old <- function(x=NULL,FUN=NULL, it=NULL,is=NULL,new=FALSE,
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
    print(paste('map.station.old',FUN))
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
        if (FUN=="NULL") FUN <- NULL else {
          if (FUN=='trend') FUN <- 'trend.coef'
          FUN <- genfun(x,FUN)
        }
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
  
  ## KMP 2025-10-08: nothing from arg is used. Why print it?
  #if (verbose) print(paste("List of arguments in the three-dots listed below ",arg,sep=""))
  
  if(projection=="np" ) {
    projection <- "sphere"
    latR <- 90 
  } else if (projection=="sp") {
    projection <- "sphere"
    latR <- -90
  }
  
  if(is.null(type)) type <- "p"
  
  if (verbose) {
    print(paste(projection,'projection'))
    if(projection=="sphere") print(paste('lonR =', lonR, 'latR = ', latR)) 
    if(!is.null(xlim)) print(paste('xlim =', paste(xlim, collapse="-")))
    if(!is.null(ylim)) print(paste('ylim =', paste(ylim, collapse="-")))
  }
  
  if (sum(is.element(type,c('fill','contour')))) {
    ## KMP 2025-10-08: If type is fill or contour, there should be a redirect to map2sphere, 
    ##  right, going through map which will redirect. We could also go on to also 
    ##  add point values for stations, but this has not been implemented yet.
    x0 <- x
    x <- as.field(x)
    z <- map(x,lonR=lonR,latR=latR,projection=projection,xlim=xlim,ylim=ylim,
             col_contour=col_contour, breaks_contour=breaks_contour,
             lab=lab,type=type,gridlines=gridlines,colbar=colbar,new=new,...)
    invisible(z)
  }
  
  if (projection=="sphere") {
    sphere(x,lonR=lonR,latR=latR,axiR=axiR,
           gridlines=gridlines,xlim=xlim,ylim=ylim,
           col=colbar$col,new=new,FUN=FUN,cex=cex,
           cex.main=cex.main,cex.axis=cex.axis,cex.lab=cex.lab,
           verbose=verbose,colbar=colbar,
           xlab=xlab,ylab=ylab,...)
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
                   legend.only = TRUE, zlim = range(colbar$breaks),
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
                   cex.axis=1,cex.lab=1,cex.main=1.5,pch=21,
                   colbar= list(pal='t2m',col=NULL,rev=FALSE,n=10,
                                breaks=NULL,type="p",cex=2,h=0.6, v=1,
                                pos=0.1,show=TRUE),
                   new=TRUE,verbose=FALSE,...) {
  if(verbose) print(paste('sphere:',lonR,latR))
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
    #map <- x
    if(length(x)==length(lon(x))) map <- x else map <- rep(1, length(lon(x)))
  }
  
  ## Initialise colbar
  colbar <- colbar.ini(map, colbar=colbar, verbose=verbose)
  breaks <- colbar$breaks
  colb <- colbar$col
  col <- colb[findInterval(map, breaks)]
  bg <- col
  nc <- length(colb)
  
  # Rotation:
  # longitudinal rotation
  #if(!is.null(xlim)) {
  #  lonR <- mean(xlim, na.rm=TRUE)
  #} else if (is.null(lonR)) {
  #  lonR <- mean(lon, na.rm=TRUE)  
  #}
  # latitudinal rotation
  #if(!is.null(ylim)) {
  #  latR <- mean(ylim, na.rm=TRUE) 
  #} else if (is.null(latR)) {
  #  latR <- mean(lat, na.rm=TRUE)
  #}
  
  # Rotate xlim and ylim
  if(is.null(xlim)) xlim <- c(-180, 180)
  if(is.null(ylim)) {
    if(latR==90) ylim <- c(0, 90) else
      if(latR==-90) ylim <- c(-90, 90) else
        ylim <- c(latR - 90, latR)
      
  }
  ## KMP 2025-10-09: Create a lon/lat grid from xlim/ylim, then transform 
  ##  the grid to stereographic coordinates and get Xlim/Ylim based on the 
  ##  range of the rotated grid
  xgrid <- seq(min(xlim), max(xlim), diff(range(xlim))/20)
  ygrid <- seq(min(ylim), max(ylim), diff(range(ylim))/20)
  nx <- length(xgrid)
  ny <- length(ygrid)
  xgrid <- rep(xgrid, ny)
  ygrid <- as.vector(sapply(ygrid, rep, nx))
  xy_grid <- cartesian2sphere(xgrid, ygrid, lonR=lonR, latR=latR)
  Xlim <- range(xy_grid$X[xy_grid$visible])
  Ylim <- range(xy_grid$Y[xy_grid$visible])
  
  ## KMP 2025-10-09: Adding space under the map to fit the color bar
  dy <- 0.3*diff(Ylim)
  if(colbar$show) Ylim <- Ylim + c(-1,0)*dy else Ylim <- Ylim + c(-0.1,0)*dy
  
  # coastline data:
  geoborders <- NULL # KMP 2019-10-11: create dummy to avoid warning during CHECK
  data("geoborders",envir=environment())
  ## If the input data is in dateline format (longitude 0 - 360)
  ## the geoborders longitude also needs to be transformed to avoid problems
  if(max(lon)>180 & min(lon)>=0) {
    gx <- geoborders$x
    gx[!is.na(gx) & gx < 0] <- gx[!is.na(gx) & gx < 0] + 360
    geoborders$x <- gx
  }
  
  gx <- geoborders$x
  gy <- geoborders$y
  ok <- is.finite(gx) & is.finite(gy)
  ## KMP 2024-08-09: Apply xlim and ylim to the geoborders
  if (!is.null(xlim)) ok <- ok & gx>=min(xlim) & gx<=max(xlim)
  if (!is.null(ylim)) ok <- ok & gy>=min(ylim) & gy<=max(ylim)
  
  ## KMP 2025-02-06: Calculating spherical coordinates with a separate function
  xy_geoborders <- cartesian2sphere(gx[ok], gy[ok], lonR=lonR, latR=latR)
  x_geoborders <- xy_geoborders$X[xy_geoborders$visible]
  y_geoborders <- xy_geoborders$Y[xy_geoborders$visible]
  
  ## KMP 2025-02-11: Rotate coordinates using the function cartesian2sphere
  xy_stations <- cartesian2sphere(lon, lat, lonR=lonR, latR=latR)
  X <- xy_stations$X
  Y <- xy_stations$Y
  Visible <- xy_stations$visible
  
  ## Define colour pal:
  if (!is.null(FUN)) {
    if (is.null(col)) col <- colscal(n=n,pal=varid(x0)[1]) else
      if (length(col)==1) {
        col <- colscal(pal=col,n=n)
      }
    nc <- length(col)
    index <- round( nc*( map - min(map) )/
                      ( max(map) - min(map) ) )
  }
  
  # Plot the results:
  if(new) dev.new()
  par(bty="n",xaxt="n",yaxt="n",new=TRUE)
  plot(Xlim, Ylim, pch=".", col="white", xlab="", ylab="",
       Xlim=Xlim, Ylim=Ylim, cex.axis=cex.axis, cex.lab=cex.lab)
  
  points(X[Visible], Y[Visible], cex=cex, pch=pch, 
         col=col[Visible], bg=bg[Visible])
  
  ## Add contour lines
  ## Plot the coast lines
  points(x_geoborders, y_geoborders, pch=".")
  ## Add circle around globe. How does this work in stereographic coordinates?
  #lines(cos(pi/180*1:360),sin(pi/180*1:360),col="black")
  lon_circle <- seq(min(xlim), max(xlim), diff(range(xlim))/100)
  lat_circle <- rep(ylim[which.min(abs(ylim))], 100)
  xy_circle <- cartesian2sphere(lon_circle, lat_circle, lonR=lonR, latR=latR)
  lines(xy_circle$X[xy_circle$visible], xy_circle$Y[xy_circle$visible], col="black")
  
  # Colorbar
  if (!is.null(FUN) & colbar$show) {      
    ## Generate a label (variable name and unit)
    label <- generate_varlabel(x0)
    ## Where to place colorbar
    xlim <- par()$usr[1:2]
    ylim <- par()$usr[3:4]
    dy <- diff(ylim)*0.1
    #below <- c(min(xlim), min(ylim)-dy/2, max(xlim), min(ylim)+dy/2)
    below <- c(min(xlim), min(ylim)+dy, max(xlim), min(ylim)+2*dy)
    dy_below <- below[4]-below[2]
    rect(below[1], below[2], 
         below[3], below[4]-dy_below*0.2, 
         col = "white", border = "white")
    ## Add colorbar and title
    col.bar(below[1],below[2]+dy_below*0.1,below[3],below[4]-dy_below*0.1,
            colbar$breaks,horiz=TRUE,pch=15,v=1,h=1,
            col=colbar$col,cex=2,cex.lab=colbar$cex.lab,
            type=colbar$type,verbose=FALSE,vl=1,border=FALSE)
    title(sub = label, line = -1, cex.sub = cex.lab)
  }
  
  if (inherits(x0,"stationmeta")) {
    result <- data.frame(x=X[Visible], y=Y[Visible])
  } else if (inherits(x0,"station")) {
    result <- data.frame(x=X[Visible], y=Y[Visible], z=map[Visible])
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

#' @export map.stationsummary
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
