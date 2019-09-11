# Documentaion in map.R
#' @export
lonlatprojection <- function(x,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                             xlim=NULL,ylim=NULL,zlim=NULL,
                             colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                          pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                             type=c("fill","contour"),gridlines=FALSE,
                             verbose=FALSE,geography=TRUE,fancy=FALSE,
                             main=NA,...) {
  
  if (verbose) {print('lonlatprojection'); str(x)}
  colid <- 't2m'; if (is.precip(x)) colid <- 'precip'
  colorbar <- !is.null(colbar)
  
  colbar <- colbar.ini(x,FUN=NULL,colbar=colbar,verbose=verbose)
  
  fig0 <- c(0,1,0,1)                        # REB 2015-06-25
  data("geoborders",envir=environment())
  if(sum(is.finite(x))==0) stop('No valid data')
  ## To deal with grid-conventions going from north-to-south or east-to-west:
  if(is.null(xlim)) xlim <- range(lon(x))
  if(!any(xlim<0) & any(xlim>180)) {
    greenwich <- TRUE
  } else {
    greenwich <- FALSE
  }
  if(inherits(x,"matrix") & is.null(attr(x,"dimensions"))) {
    x <- g2dl(x,d=c(length(lon(x)),length(lat(x)),1),
              greenwich=greenwich,verbose=verbose)
  } else {
    x <- g2dl(x,greenwich=greenwich,verbose=verbose)
  }
  dim(x) <- c(length(lon(x)),length(lat(x)))
  
  #lon <- lon(x)
  #if(!any(xlim<0) & any(xlim>180)) {
  #  lon[lon<0] <- lon[lon<0]+360
  #} else {
  #  lon[lon>180] <- lon[lon>180]-360
  #}
  srtx <- order(lon(x)); lon <- lon(x)[srtx]
  srty <- order(lat(x)); lat <- lat(x)[srty]
  if (verbose) print('meta-stuff')
  unit <- attr(x,'unit'); variable <- varid(x); varid <- varid(x); isprecip <- is.precip(x)
  
  if(!is.null(variable)) {
    variable <- as.character(variable)
   varname <- attr(x,'longname') 
  }
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
  if (!is.null(attr(x,'source'))) {
    sub <- attr(x,'source')
  } else {
    sub <- NULL
  }
  if (sum(is.element(type,'fill'))==0) colbar <- NULL
  
  if (verbose) print('time')
  if (!is.null(attr(x,'timescale'))) {
    if (verbose) print(attr(x,'timescale'))
    timescale <- attr(x,'timescale')
    if (timescale == 'annual') {
      t1 <- year(attr(x,'time'))[1]
      t2 <- year(attr(x,'time'))[2]
    } else if (sum(is.element(c('month','season'),timescale))>0) {
      t1 <- paste0(year(attr(x,'time'))[1],"~",month(attr(x,'time'))[1])
      t2 <- paste0(year(attr(x,'time'))[2],"~",month(attr(x,'time'))[2])
    } else {
      t1 <- attr(x,'time')[1]  
      t2 <- attr(x,'time')[2]
    }
    ##period <- paste('[',t1,', ',t2,']',sep='')  ## REB: square brackets have special role in expressions
    period <- paste('phantom(0)* (',t1,'-',t2,')',sep='')
  } else period <- NULL
  if (verbose) print(paste('period:',period))
  method <- attr(x,'method')
  if (verbose) {
    print(c(dim(x),length(srtx),length(srty)))
    # There is something strange happening with x - in some cases it is filled with NAs (REB)
    print(srtx); print(srty)
  }
  x <- x[srtx,srty]
  
  if (verbose) {print(xlim); str(x)}
  if (!is.null(xlim)) {
    outside <- (lon < min(xlim)) | (lon > max(xlim))
    if (verbose) print(paste('mask',sum(outside),length(outside)))
    x[outside,] <- NA
  } else xlim <- range(lon)
  
  if (!is.null(ylim)) {
    outside <- (lat < min(ylim)) | (lat > max(ylim))
    if (verbose) print(paste('mask',sum(outside),length(outside)))
    x[,outside] <- NA
  } else ylim=range(lat)
  
  if (new) {
    dev.new()
    par(fig=fig0)
    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE)
  } else {
    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE)
    fig0 <- par()$fig
  }
  
  if (verbose) print('Set up the figure')
  
  plot(range(lon),range(lat),type="n",xlab="",ylab="", # REB 10.03
       xlim=xlim,ylim=ylim,main=main, # to sumerimpose.
       xaxt="n",yaxt="n") # AM 17.06.2015
  ##par0 <- par()
  
  if (sum(is.element(tolower(type),'fill'))>0)   
    image(lon,lat,x,xlab="",ylab="",add=TRUE,
          col=colbar$col,breaks=colbar$breaks,xlim=xlim,ylim=ylim,...)
  
  if (geography) {
    lines(geoborders$x,geoborders$y,col="darkblue")
    lines(attr(geoborders,'borders')$x,attr(geoborders,'borders')$y,col="pink")
    lines(geoborders$x+360,geoborders$y,col="darkblue")
  }
  if (sum(is.element(tolower(type),'contour'))>0)
    contour(lon,lat,x,lwd=1,col="grey70",add=TRUE)
  if (gridlines) grid()
  par(xpd=FALSE)
  dlat <- diff(range(lat))/60
  if (verbose) {print(dlat); print(sub);  print(varlabel)}
  
  if(!is.null(varlabel)) lab <- paste(varlabel,'*') else lab <- ''
  lab <- as.expression(parse(text=paste(lab,'phantom(0) - phantom(0)')))
  ## text(lon[1],lat[length(lat)] - 0.5*dlat,varlabel,pos=4,font=2, cex=0.85)
  
  ## if ((!is.null(sub)) & (length(sub)>0)) text(lon[1],lat[1] - 1.5*dlat,sub,col="grey30",pos=4,cex=0.7)
  if ((!is.null(sub)) & (length(sub)>0)) {
    sub <- paste('pattern derived from',sub)
    lab <- try(parse(text=paste(lab,'*',as.expression(paste('~ ',
                                paste(unlist(strsplit(sub,split=' ')),
				collapse = ' *~ '), sep='')))))
    if (inherits(lab,'try-error')) lab <- ''
  }  #title(main = as.expression(sub),line = 3, adj =0.25)
  
  if (!is.null(method)) {
    lab <- try(parse(text=paste(lab,'*',as.expression(method))))
    if (inherits(lab,'try-error')) lab <- ''
    #title(main = as.expression(method),line = 3, adj =0.5)
    #text(lon[length(lon)],lat[1] - dlat,method,col="grey30",pos=2,cex=0.7)
  }
  if (!is.null(period)) {
    lab <- try(parse(text=paste(lab,'*',as.expression(period))))
  }
  if (inherits(lab,'try-error')) lab <- ''
  #title(main = as.expression(period),line = 3, adj =1)
  #text(lon[length(lon)],lat[length(lat)] + 0.5*dlat,period,pos=2,cex=0.7,col="grey30")

  title(sub = lab,line = 0 , adj = 0.5)
  if (!is.null(colbar)) {
    if (verbose) print('Add colourbar')
    
    par(xaxt="s",yaxt="s",las=1,col.axis='grey',col.lab='grey',
        cex.lab=0.7,cex.axis=0.7)
    axis(2,at=pretty(lat(x)),col='grey')
    axis(3,at=pretty(lon(x)),col='grey')
    if(gridlines) grid()
    
    par(col.axis='black',col.lab='black',
        cex.lab=0.5,cex.axis=0.5)
    
    #if (!is.null(colbar)) {
    if (colbar$show) {
      if (fancy) {
        col.bar(colbar$breaks,horiz=TRUE,pch=21,v=1,h=1,
                col=colbar$col, cex=2,cex.lab=colbar$cex.lab,
                type=type,verbose=FALSE,vl=1,border=FALSE)
        #}
        #  }
      } else {
        #par(fig=par0$fig)
        #op <- par()
        #par(mgp = c(0, 2, 0))
        image.plot(breaks=colbar$breaks,
                   lab.breaks=colbar$breaks,horizontal = TRUE,
                   legend.only = TRUE, zlim = range(colbar$breaks),
                   col = colbar$col, legend.width = 1,
                   axis.args = list(cex.axis = 1,hadj = 0.5,mgp = c(0, 0.5, 0)), border = FALSE)
        #par(op)
      }
    }
  }
  par(fig=fig0)
  
  par(col.axis='black',col.lab='black',cex.lab=1,cex.axis=1,
      xaxt="s",yaxt="s")
  result <- list(x=lon,y=lat,z=x,breaks=colbar$breaks)
  #par(fig=par0$fig)
  invisible(result)
}
