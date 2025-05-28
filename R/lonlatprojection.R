# Documentaion in map.R
#' @export
lonlatprojection <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                             xlim=NULL,ylim=NULL,zlim=NULL,lab='default',
                             colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                          pos=0.05,show=TRUE,type="r",cex=2,
                                          srt=0.45,cex.lab=0.7,h=0.6,v=1),
                             type=c("fill","contour"),gridlines=FALSE,
                             col_contour="grey70",
                             verbose=FALSE,geography=TRUE,fancy=TRUE,
                             main=NA,cex.sub=0.8,cex.axis=0.8,
                             fig=NULL,add=FALSE,plot=TRUE,useRaster=TRUE,...) {
  
  if (verbose) {print('lonlatprojection'); str(x)}
  ## REB 2024-04-29
  xargs <- list(...)
  if (!is.null(xargs$showaxis)) showaxis <- xargs$showaxis else showaxis <- TRUE
  if (plot) par0 <- par() 
  attr(x,'source') <- NULL ## REB "2021-12-21: Fed up with problems with silly source information...
  ## Use temperature-palette as default, and check if the variable is precipitation
  ## for precipitation-palette
  if (verbose) print('time')
  ## Set up the attribute 'time'
  if ( (is.null(attr(x,"time"))) & (length(index(x))>1) ) attr(x,"time") <- range(index(x))
  if (!is.null(attr(x,'timescale'))) {
    if (verbose) print(paste('timescale',attr(x,'timescale')))
    timescale <- attr(x,'timescale')
    if (timescale == 'annual') {
      t1 <- year(attr(x,'time'))[1]
      t2 <- year(attr(x,'time'))[2]
    } else if (sum(is.element(c('month','season'),timescale))>0) {
      t1 <- paste0(year(attr(x,'time'))[1],"~",month(attr(x,'time'))[1])
      t2 <- paste0(year(attr(x,'time'))[2],"~",month(attr(x,'time'))[2])
    } else if (grepl("hour|minute",timescale)) {
      t1 <- paste0(strftime(attr(x,"time")[1],format="%Y%m%d"),"~",
                   strftime(attr(x,"time")[1],format="%H:%M"))
      t2 <- paste0(strftime(attr(x,"time")[2],format="%Y%m%d"),"~",
                   strftime(attr(x,"time")[2],format="%H:%M"))
    } else if(inherits(attr(x,"time"),c("Date","POSIXt"))) {
      t1 <- strftime(attr(x,"time")[1],format="%Y%m%d")
      t2 <- strftime(attr(x,"time")[2],format="%Y%m%d")
    } else {
      t1 <- attr(x,"time")[1]
      t2 <- attr(x,"time")[2]
    }
  } else { 
    if (verbose) print('Attribute timescale is not provided')
    if (!is.null(attr(x,'time'))) {t1 <- attr(x,'time')[1]; t2 <- attr(x,'time')[2]} else
      if (is.null(index(x))) {t1 <- min(index(x)); t2 <- max(index(x))} else 
      {t1 <- NA; t2 <- NA}
  }
  if (verbose) print(paste('t1=',t1,'t2=',t2))
  
  colid <- 't2m'; if (is.precip(x)) colid <- 'precip'
  ## If colbar is set to NULL then remember this and do not show the colourbar
  show.colbar <- !is.null(colbar)
  ## Prepare the colourbar nevertheless...
  colbar <- colbar.ini(x,FUN=NULL,colbar=colbar,verbose=FALSE)
  varnm <- varid(x); unitx <- esd::unit(x)
  ## REB 2021-12-21: Sometimes the source information is a bit overwhelming and that too creates a problem
  if (!is.null(src(x)))
    if (nchar(src(x))> 10) {
      first.space <- gregexpr('_',src(x))[[1]]
      if (!is.na(first.space)) attr(x,'source') <- substr(src(x),1,first.space-1)
    }
  
  ## Land contours
  data("geoborders",envir=environment())
  if(!is.null(attr(x,"greenwich"))) if(!attr(x,"greenwich")) {
    gbl <- geoborders$lon
    gbl[gbl < 0] <- gbl[gbl < 0] + 360
    geoborders$lon <- gbl
  }
  if(sum(is.finite(x))==0) stop('No valid data')
  ## To deal with grid-conventions going from north-to-south or east-to-west:
  if(is.null(xlim)) xlim <- range(lon(x))
  if(!any(xlim<0) & any(xlim>180)) {
    greenwich <- TRUE
  } else {
    greenwich <- FALSE
  }
  ## Make sure to use the right arrangement: from dateline or Greenwich 
  if (verbose) {print(dim(x)); print(attr(x,'greenwich')); print(greenwich)}
  if(inherits(x,"matrix") & is.null(attr(x,"dimensions"))) {
    x <- g2dl(x,d=c(length(lon(x)),length(lat(x)),1),
              greenwich=greenwich,verbose=verbose)
  } else {
    x <- g2dl(x,greenwich=greenwich,verbose=verbose)
  }
  if (verbose) print(paste('dimensions of x:',paste(dim(x),collapse=' - '),
                           'tim=',length(index(x)),'lon=',length(lon(x)),
                           'lat=',length(lat(x))))
  dim(x) <- c(length(lon(x)),length(lat(x)))
  ## Make sure the longitudes are ordered correctly
  srtx <- order(lon(x)); lon <- lon(x)[srtx]
  srty <- order(lat(x)); lat <- lat(x)[srty]
  
  if (verbose) print('meta-stuff')
  unit <- attr(x,'unit'); variable <- varid(x); isprecip <- is.precip(x)
  
  ## Not used ? REB
  # if(!is.null(variable)) {
  #   variable <- as.character(variable)
  #  varname <- attr(x,'longname') 
  # }
  if(!is.null(unit)) unit <- as.character(unit)
  if ( (unit=="degC") | (unit=="deg C") | (unit=="degree C") | (unit=="degree Celsius"))
    unit <- "degree*C"
  if (unit=="%") unit <- "'%'"
  if(!is.null(variable)) {
    if ( (tolower(variable)=="t(2m)") | (tolower(variable)=="t2m") |
         (tolower(variable)=="2t") )
      variable <- "T[2*m]"
  }
  if (verbose) print(paste(variable,unit,isprecip,' -> varlabel'))
  if(!is.null(variable)) {
    varlabel <- try(eval(parse(
      text=paste('expression(',gsub(" ","~",variable)," *~(",gsub(" ",
                                                                  "~",unit),"))",sep="")))) 
  } else {
    varlabel <- NULL
  }
  if (inherits(varlabel,'try-error')) varlabel <- NULL
  if (verbose) {print(varlabel); print(src(x))}
  
  #if (is.null(src(x))) attr(x,'source') <- NA
  if (!is.null(src(x))) {
    ## KMP 2019-12-12: Added check of source attribute because basename can't handle NA
    ## But why is source set to NA in as.field as default rather than NULL?
    if (length(src(x))==0) attr(x,'source') <- 'data'
    
    if (is.na(src(x))) {
      sub <- NULL
    } else {
      sub <- paste(basename(src(x)),'*phantom(0)',sep='')
    }
  } else {
    sub <- NULL
  }
  if (sum(is.element(type,'fill'))==0) colbar <- NULL
  
  if (!is.null(t1) & !is.null(t2)) { 
    ##period <- paste('[',t1,', ',t2,']',sep='')  ## REB: square brackets have special role in expressions
    if (is.null(attr(x,'period'))) period <- paste('phantom(0)* (',t1,'-',t2,')',sep='') else
      period <- paste('phantom(0)* (',attr(x,'period'),')',sep='')
    if (verbose) print(period)
  } else {
    if (verbose) print('<No time information provided!>')
    period <- NA
    t1 <- t2 <- NA
  }
  if (verbose) print(paste('period:',period))
  method <- attr(x,'method')
  if (verbose) {
    print(c(dim(x),length(srtx),length(srty)))
    # There is something strange happening with x - in some cases it is filled with NAs (REB)
    # print(srtx); print(srty)
  }
  x <- x[srtx,srty]
  
  if (verbose) {print(xlim); str(x)}
  if (!is.null(xlim)) {
    outside <- (lon < min(xlim)) | (lon > max(xlim))
    if (verbose) print(paste('mask',sum(outside),length(outside)))
    x[outside,] <- NA
  } else xlim <- range(lon)
  
  dy <- 0.2*( max(lat) - min(lat) )
  #print(dy); print(range(lat))
  if (!is.null(ylim)) {
    outside <- (lat < min(ylim)) | (lat > max(ylim))
    if (verbose) print(paste('mask',sum(outside),length(outside)))
    x[,outside] <- NA
    ylim <- ylim + c(-dy,0)
  } else ylim=range(lat) + c(-dy,0)
  
  if (new & plot) {
    if(verbose) print("Create new graphic device")
    dev.new()
    add <- FALSE
  }
  if(!is.null(fig) & plot) par(fig=fig,new=(add & dev.cur()>1))
  if (plot) { 
    par(bty="n",xpd=FALSE)
    
    if (verbose) print('Set up the figure')
    plot(range(lon),range(lat),type="n",xlab="",ylab="", # REB 10.03
         xlim=xlim,ylim=ylim,main=main,
         xaxt="n",yaxt="n") # AM 17.06.2015
    
    if (useRaster) {
      ## ‘useRaster = TRUE’ can only be used with a regular grid
      ## Force the coordinates to be evenly spaced
      if (verbose) print('ensure regular grid')
      lon <- seq(min(lon),max(lon),length=length(lon))
      lat <- seq(min(lat),max(lat),length=length(lat))
    }
    if (sum(is.element(tolower(type),'fill'))>0)   
      image(lon,lat,x,xlab="",ylab="", add=TRUE,useRaster = useRaster,
            col=colbar$col,breaks=colbar$breaks)
    
    if (geography) {
      lines(geoborders$x,geoborders$y,col="darkblue")
      lines(attr(geoborders,'borders')$x,attr(geoborders,'borders')$y,
            col="pink")
      lines(geoborders$x+360,geoborders$y,
            col="darkblue")
    }
    if (sum(is.element(tolower(type),'contour'))>0)
      if(is.null(breaks_contour)) breaks_contour <- colbar$breaks
      contour(lon,lat,x,lwd=1,col=col_contour,levels=breaks_contour,
              add=TRUE,...)
    if (gridlines) grid()
    if (showaxis) { 
      axis(2,at=pretty(lat),col='grey',cex=cex.axis)
      axis(3,at=pretty(lon),col='grey',cex=cex.axis)
    }
    ## REB 2023-01-24
    par(xpd=TRUE)
  }
  dlat <- diff(range(lat))/60
  if (verbose) {print(dlat); print(sub);  print(varlabel)}
  
  if(!is.null(varlabel) & (lab=='default')) label <- paste(varlabel,'*') else label <- ''
  label <- as.expression(parse(text=paste(label,'phantom(0) - phantom(0)')))
  
  if ((!is.null(sub)) & (length(sub)>0)) {
    sub <- paste('pattern derived from',sub)
    label <- try(parse(text=paste(label,'*',as.expression(paste('~ ',
                                                                paste(unlist(strsplit(sub,split=' ')),
                                                                      collapse = ' *~ '), sep='')))))
    if (inherits(label,'try-error')) label <- ''
  } 
  
  if (!is.null(method)) {
    label <- try(parse(text=paste(label,'*',as.expression(method))))
    if (inherits(label,'try-error')) label <- ''
    #title(main = as.expression(method),line = 3, adj =0.5)
    #text(lon[length(lon)],lat[1] - dlat,method,col="grey30",pos=2,cex=0.7)
  }
  if (!is.null(period)) {
    label <- try(parse(text=paste(label,'*',as.expression(period))))
  }
  
  if (inherits(label,'try-error')) label <- ''
  #title(main = as.expression(period),line = 3, adj =1)
  #text(lon[length(lon)],lat[length(lat)] + 0.5*dlat,period,pos=2,cex=0.7,col="grey30")
  #if (lab=='simple')   label <- eval(parse(text=paste('expression(',varnm,')'))) else 
  if (lab=='simple')   label <- eval(parse(text=paste('expression(',varnm,')'))) else 
    if (lab=='unit') label <- eval(parse(text=paste('expression(',varnm,'* phantom0 * (',unitx,')',')'))) else
      if (is.character(lab) & lab!="default") label <- lab
  ## KMP 2023-02-08: removing title here and trying to place it above colorbar instead
  #title(sub = label, line = 0, adj = 0.5, cex.sub = cex.sub)
  
  if (show.colbar & plot) {
    if (verbose) print('Add colourbar')
    ## REB 2023-01-24
    #par(xaxt="s",yaxt="s",las=1,col.axis='grey',col.lab='grey',
    #    cex.lab=0.7,cex.axis=0.7)
    ## REB 2023-01-24
    #par(col.axis='black',col.lab='black',
    #    cex.lab=0.5,cex.axis=0.5)
    if(!is.null(colbar$show)) if (colbar$show) {
      if (verbose) print('Show colourbar')
      if (fancy) {
        if (verbose) print('fancy')
        dy <- diff(ylim)*0.1
        below <- c(min(xlim), min(ylim)-dy/2, max(xlim), min(ylim)+dy/2)
        dy_below <- below[4]-below[2]
        rect(below[1], below[2], 
             below[3], below[4]-dy_below*0.2, 
             col = "white", border = "white")
        col.bar(below[1],below[2]+dy_below*0.1,below[3],below[4]-dy_below*0.1,
                colbar$breaks,horiz=TRUE,pch=15,v=1,h=1,srt=colbar$srt,
                col=colbar$col,cex=2,cex.lab=colbar$cex.lab,
                type=colbar$type,verbose=FALSE,vl=1,border=FALSE)
        title(sub = label, line = 1, cex.sub = cex.sub)
      } else {
        title(sub = label, line = 1, cex.sub = cex.sub)
        image.plot(breaks=colbar$breaks, 
                   lab.breaks=colbar$breaks,horizontal = TRUE,
                   legend.only = TRUE, zlim = range(colbar$breaks),
                   col = colbar$col, legend.width = 1,
                   axis.args = list(cex.axis = 1,hadj = 0.5,mgp = c(0, 0.5, 0)), 
                   border = FALSE, verbose=FALSE)
      }
    }
    if (!is.null(fig)) par(fig=fig)
  }
  if (plot) par(bty=par0$bty,xpd=par0$xpd,col.axis=par0$col.axis,col.lab=par0$col.lab, 
                cex.lab=par0$cex.lab, cex.axis=par0$cex.axis,xaxt=par0$xaxt,yaxt=par0$yaxt,
                new=FALSE)
  if (verbose) print('Add attributes to returned results')
  attr(x,'longitude') <- lon
  attr(x,'latitude') <- lat
  attr(x,'variable') <- variable
  attr(x,'unit') <- unit
  attr(x,'colbar') <- colbar
  if(!is.null(t1) & !is.null(t2)) attr(x,'time') <- c(t1,t2)
  invisible(x)
}

