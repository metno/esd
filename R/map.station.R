## Author 	 Rasmus E. Bnestad
## Updated 	 by Abdelkader Mezghani
## Last update   26.07.2013
## Includes	 map.station() ; test.map.station()
## Require 	 geoborders.rda

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
    class(x) <- "stationmeta"
    map.station(x,...)
}

map.station <- function (x = NULL,col = NULL,bg="green",cex=.8, zexpr = "alt",
                         is=list(x=NULL,stid = NULL, param = NULL, lon = NULL,
                             lat = NULL,alt = NULL, cntr = NULL, src = NULL, nmin = NULL),
                         it = NULL,
                         col.subset="darkred",bg.subset="red",cex.subset=1,
                         add.text.subset=FALSE,
                         colbar=list(col=NULL,breaks=NULL,n=10,type="p",cex=2,h=0.6,v=1),
                         showall = FALSE, verbose = FALSE , add.text=FALSE,
                         height=NULL,width=NULL,cex.axis=1,cex.lab=0.6,pch=21,
                         FUN=NULL,from=NULL,to=NULL,showaxis=FALSE,xlim = NULL,
                         ylim = NULL,border=FALSE, full.names=FALSE,
                         full.names.subset=FALSE,new=TRUE,text=FALSE, fancy=FALSE,
                         projection="lonlat",what=NULL,gridlines=FALSE, lonR=NULL,
                         latR=45,axiR=NULL,na.rm=TRUE,...) 
{ 
    
    par0 <- par()
    ##print(par()$fig)
    ## X <- coredata(x)
    ##if (dim(X)[1]==1) X <- coredata(x[1,]) else
    ##if (inherits(X,"matrix")) X <- apply(x,2,FUN=FUN,na.rm=TRUE)
    ##print(length(X))
    ##X <- diag(X,nrow=length(lon(x)),ncol=length(lat(x)))
    ##attr(X,'longitude') <- attr(x,'longitude')
    ##attr(X,'latitude') <- attr(x,'latitude')
    ##attr(X,'variable') <- attr(x,'variable')
    ##  if (attr(x,'unit')=="deg C") attr(X,'unit') <- expression(degree*C) else
    ##attr(X,'unit') <- attr(x,'unit')
    ##attr(X,'source') <- attr(x,'source')
    ##attr(X,'time') <- range(index(x))
    ##attr(X,'method') <- FUN
    ##attr(X,'timescale') <- class(x)[2]
    ##print(length(X)); print(attr(x,'dimensions'))
    ##X <- diag(X,nrow=length(lon(X)),ncol=length(lat(X)))
    ##dim(X) <- c(length(lon(x)),length(lat(x)))
    ##class(X) <- class(x)
    ##str(X)
    
    if (projection=="sphere")
        sphere(x,lonR=lonR,latR=latR,axiR=axiR,
                   gridlines=gridlines,
                   col=col,new=new,FUN=FUN,...)
    else if (projection=="np")
        sphere(x,lonR=lonR,latR=90,axiR=axiR,
                   gridlines=gridlines,
                   col=col,new=new,FUN=FUN,...) else
    if (projection=="sp")
        sphere(x,lonR=lonR,latR=-90,axiR=axiR,
                   ,gridlines=gridlines,
                   col=col,new=new,FUN=FUN,...)
    ## else if (projection=="lonlat")
    ##    lonlatprojection(x=X,xlim=xlim,ylim=ylim, n=colbar$n,col=colbar$col,breaks=colbar$breaks,new=new,
    ##                     what=what,gridlines=gridlines,...)
    else if (projection=="lonlat") {
        
        ## setting default values for the color bar if not specified
                                        #if (is.null(colbar)) {
        #if (length(col) > 1) {
        #    if (is.null(FUN)) {
        #        print("Color vector length is higher than 1, only the first value is used")
        #    }
        #    else if (is.null(colbar$col)) {
        #        colbar$n <- 10
        #        colbar$col <- colscal(n=colbar$n)
        #    }
        #    else {
        #        colbar$col <- col
        #        colbar$n.breaks <- length(col)
        #    }
        #} else col <- col[1]
        
        ##load("~/esd/data/geoborders.rda")
        data("geoborders", envir = environment())
        ##load("station.meta.rda") # data("station.meta", envir = environment())
        ##load("t2m.NORDKLIM.rda") # data("t2m.NORDKLIM", envir = environment())
        ##data("precip.NORDKLIM", envir = environment())
        if (zexpr == "alt") 
            zexpr <- "sqrt( station.meta$alt/max(station.meta$alt,na.rm=TRUE) )"
        ## n <- 100
        
        ## if (!is.null(subset)) {col = "white" ; bg="white"} 
        
        ## if (!is.null(unlist(subset)) & !showall) {col = "white" ; bg="white"} 

        if (!is.null(x)) { 
            if (inherits(x,"stationmeta")) {      
                ss <- x
                ss$variable <-apply(as.matrix(ss$element),1,esd2ele)
            }
            else if (inherits(x,"station")) {
                if (is.null(dim(x))) dim(x) <- c(length(x),1)
                ss <- list(station_id=attr(x,"station_id"),location=attr(x,'location'),country=attr(x,'country'),longitude=attr(x,"longitude"),latitude=attr(x,'latitude'),altitude=attr(x,'altitude'),variable=attr(x,"variable"),longname=attr(x,"longname"),start=rep(start(x),dim(x)[2]),end=rep(end(x),dim(x)[2]),source=attr(x,"source"))
            }
        } else
            ss <- select.station()
        
        if (is.null(attr(ss,"element")))
            ss$element <-apply(as.matrix(ss$variable),1,esd2ele)   

        if (!is.null(unlist(is))) { ## highlight a subset of station
            
            if (is.null(is$x)) {
                highlight <- select.station(x=is$x,loc=is$loc,stid = is$stid, param = is$param, lon = is$lon, lat = is$lat, alt = is$alt, cntr = is$cntr, src = is$src, nmin = is$nmin , it = is$it)
                highlight$variable <-apply(as.matrix(highlight$element),1,esd2ele)
            } else if (inherits(is$x,"station")) {
                if (is.null(dim(is$x)))
                    dim(is$x) <- c(length(is$x),1)
                highlight <- list(station_id=attr(is$x,"station_id"),location=attr(is$x,'location'),country=attr(is$x,'country'),longitude=attr(is$x,"longitude"),latitude=attr(is$x,'latitude'),altitude=attr(is$x,'altitude'),variable=attr(is$x,"variable"),longname=attr(is$x,"longname"), start=rep(start(is$x),dim(is$x)[2]),end=rep(end(is$x),dim(is$x)[2]),source=attr(is$x,"source"))
            }
            highlight$element <- apply(as.matrix(highlight$variable),1,esd2ele)
        }  else highlight <-  NULL
        
        
        ## Set negative altitude to NA
        ss$altitude[ss$altitude < 0] <- NA
        
        tte <- "rwb"
        
        cols <- rgb(seq(0, 0.5, length = 100)^2, seq(0.5, 1, length = 100), 
                    seq(0, 0.5, length = 100)^2)

        ## if (new) dev.new(height=height,width=width)
                                        #par(bty = "n", xaxt = "n", yaxt = "n", xpd = FALSE)

        ## Select a subdomain in the x-axis
        if (is.null(xlim))
            if (is.null(highlight) | showall)
                if (length(is$lon) > 1)
                    xlim <- floor(range(ss$longitude, na.rm = TRUE))
                else
                    xlim <-floor(range(ss$longitude, na.rm = TRUE) + c(-5,5))  # +/- 5 degrees
            else
                if (length(is$lon) > 1)
                    xlim <- floor(range(highlight$longitude, na.rm = TRUE))
                else
                    xlim <-floor(range(highlight$longitude, na.rm = TRUE) + c(-5,5))  # +/- 5 degrees
        ## Select a subdomain in the y-axis
        if (is.null(ylim))
            if (is.null(highlight) | showall)
                if (length(is$lat) > 1)
                    ylim <- floor(range(ss$latitude, na.rm = TRUE))
                else
                    ylim <- floor(range(ss$latitude, na.rm = TRUE) + c(-5,5))
            else
                if (length(is$lat) > 1)
                    ylim <- floor(range(highlight$latitude, na.rm = TRUE))
                else
                    ylim <-floor(range(highlight$latitude, na.rm = TRUE) + c(-5,5))  # +/- 5 degrees

        ## scaling factor to apply on cex ...
        if (!inherits(x,"stationmeta") & !is.null(attr(x,'na')))
            scale <- attr(x,'na')
        else
            scale <- 1
        
##        if ( (par()$mfcol[1]> 1) | (par()$mfcol[2]> 1) ) new <- FALSE
##        if (new) {
##            #dev.new()
##            par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
##            fig=c(0.05,0.95,0.12,0.95),mar=rep(1,4))
##        } else {
##            par(bty="n",xaxt="n",yaxt="n",xpd=TRUE,mar=rep(1,4),new=new)
##        }
        ## browser()
        ##print(par()$fig)
        par(fig=par0$fig)
        if (!is.null(FUN)) col <- "white" 
        if (!is.null(highlight))
            plot(highlight$longitude, highlight$latitude, pch = pch, col = col, bg = bg.all,cex = cex*scale, xlab = "", ylab = "", xlim = xlim, ylim = ylim , axes =FALSE , frame.plot = FALSE)
        else if (!is.null(ss))
            plot(ss$longitude, ss$latitude, pch = pch, col = col, bg = bg[1], cex = cex*scale, xlab = "", ylab = "", xlim = xlim, ylim = ylim , axes = FALSE , frame.plot = FALSE)
        
        ## Add geoborders
        lines(geoborders$x, geoborders$y, col = "black")
        lines(attr(geoborders, "borders")$x, attr(geoborders, "borders")$y, col = "grey90")

        ## par(fig=par0$fig)
        ## print(par()$fig)
        if (showall) {
            ss.all <- select.station(param=is$param)
            points(ss.all$longitude,ss.all$latitude,pch=".",col="grey50",bg="grey",cex=cex/2)
        }
        ## par(fig=par0$fig)
        ## print(par()$fig)
        ## add search info to plot
        
        if (text) {
            if (!is.null(highlight)) {
                title(main=paste("SOURCE(S) : ", paste(levels(factor(highlight$source)),collapse="/" )),line=3,cex.main=.8)
                title(main=paste(length(levels(factor(highlight$location)))),line=2,cex.main=.8)
                title(main=paste(min(highlight$start,na.rm=TRUE),"/",max(highlight$end,na.rm=TRUE)),line=2,cex.main=.8,adj=1)
                if (!is.null(FUN))
                    title(main=paste(paste(toupper(apply(as.matrix(levels(factor(highlight$variable))),1,esd2ele)),collapse="/"),toupper(FUN),sep="/"),line=2,cex.main=.8 , adj = 0)
                else
                    title(main=paste(toupper(apply(as.matrix(levels(factor(highlight$variable))),1,esd2ele)),collapse="/"),line=2,cex.main=.8 , adj = 0)
            }
            else {
                title(main=paste("SOURCE(S) : ", paste(levels(factor(ss$source)),collapse="/" )),line=3,cex.main=.8)
                title(main=paste(length(ss$location)),line=2,cex.main=.8)
                title(main=paste(max(ss$start,na.rm=TRUE),"/",min(ss$end,na.rm=TRUE)),line=2,cex.main=.8,adj=1)
                title(main=paste(paste(toupper(levels(factor(ss$variable))),collapse="/"),toupper(FUN),sep="/"),line=2,cex.main=.8 , adj = 0)
            }
        }
        ## title(main=attr(z,"title"),line=2.2,cex.main=0.7)
        ## add margin text
        if (text) mtext(paste(("ESD package - map.station() - MET Norway 2014"),"(www.met.no)",sep=" "),side=1,line=4,cex=0.6)
        
        ## add grid
        grid()
        ## insert color bar                                    
        if (!is.null(FUN)) {
            if (is.element(FUN,c('lon','lat','alt')))
                eval(parse(text=paste('y <-',FUN,'(x)',sep="")))
            else {
                if (is.element("na.rm",names(formals(FUN))) | is.element("...",names(formals(FUN))) | (is.element(FUN,c("max","min"))))
                    y <- apply(coredata(x),2,FUN=FUN,na.rm=TRUE)
                else if (FUN=="trend")
                    y <- apply(x,2,FUN=FUN,na.omit=FALSE)
                else
                    y <- apply(coredata(x),2,FUN=FUN) ## ,na.rm=TRUE)
            }
            ## y.rng <- floor(range(y,na.rm=TRUE))
            ## browser()           
            if (is.null(colbar$n) & !is.null(colbar$col))
                colbar$n <- length(colbar$col)
            else if (!is.null(colbar$breaks)) {              
                colbar$n <- length(colbar$breaks)
                colbar$col <- colscal(n=colbar$n)
            }
            else { ## set to the default values
                colbar$n <- 10 
            }   

            if (is.null(colbar$breaks)) {
                colbar$breaks <- pretty(y,colbar$n)
                # update colbar$n and colbar$col according to pretty
                colbar$n <- length(colbar$breaks)
                colbar$col <- colscal(n=colbar$n)
            }
            ## reverse the colour for precip
            #print(is.precip(x))
            if (is.precip(x)) colbar$col <- rev(colbar$col)
            
            # find color index in colbar
            icol <- apply(as.matrix(y),2,findInterval,colbar$breaks)
            bg <- colbar$col[icol]
            if (is.null(col)) col <- bg

            bg <- colbar$col[icol]
            if (is.null(col)) col <-bg
         
            if (verbose) print(range(y,na.rm=TRUE))
            
            par(fig=par0$fig,new=TRUE)
            
            ##scale <- apply(y,2,function(x) sum(!is.na(x))/length(x))
            if (!is.null(attr(x,'na'))) ## (!inherits(x,"stationmeta") & 
                scale <- attr(x,'na')
            else
                scale <- 1
            
            points(ss$longitude, ss$latitude, pch = pch, bg=bg , col=col,
                   cex = cex*scale, xlab = "", ylab = "", xlim = xlim, ylim = ylim,...)
            par(fig=par0$fig,new=TRUE)
            ## print(par()$fig)
            
            if (!is.null(highlight)) {
                points(highlight$longitude, highlight$latitude, pch = 21 , col = col.subset,
                       bg=bg.subset, cex = cex.subset,...)

            }

            ## Add geoborders
            lines(geoborders$x, geoborders$y, col = "black")
            lines(attr(geoborders, "borders")$x, attr(geoborders, "borders")$y, col = "grey90")
            
            ## par(fig=par0$fig)
            ## print(par()$fig)
            
            ## add color bar
            if (fancy)
                col.bar(colbar$breaks,horiz=TRUE,pch=21,v=1,h=1,col=colbar$col,
                        cex=2,cex.lab=colbar$cex.lab,type="p",verbose=FALSE,vl=1,border=FALSE)
            else
                image.plot(horizontal = TRUE, legend.only = T, zlim = range(colbar$breaks),
                           col = colbar$col, legend.width = 1, axis.args = list(cex.axis = 0.8),
                           border = FALSE)
        }    
        
        ## add text if TRUE
        if (!is.null(unlist(is)))
            if (add.text.subset)
                if (full.names.subset)
                    text(highlight$longitude, highlight$latitude,highlight$location,pos=3,cex=cex.subset/2)
                else
                    text(highlight$longitude, highlight$latitude,substr(toupper(highlight$location),1,3),pos=3,cex=cex.subset/2)
        
        if (add.text)
            if (full.names)
                text(ss$longitude, ss$latitude,ss$location,pos=3,cex=cex/2)
            else
                text(ss$longitude, ss$latitude,substr(toupper(ss$location),1,3),pos=3,cex=cex/2)
        ##add label text
        
        if (showaxis) title(xlab = "Longitude",line=2.2 , cex.lab = cex.lab) 
        if (showaxis) title(ylab = "Latitude",line=2.2 , cex.lab = cex.lab) 
        ## format axes
        if (showaxis) axis(1,seq(xlim[1],xlim[2],by=10),cex.axis=cex.axis) # 0.7
        if (showaxis) axis(2,seq(ylim[1],ylim[2],by=10),cex.axis=cex.axis)
        if (showaxis) axis(4,seq(ylim[1],ylim[2],by=10),cex.axis=cex.axis)
        if (showaxis) axis(3,seq(xlim[1],xlim[2],by=10),cex.axis=cex.axis)
        
        ## par(fig=par0$fig)
        ## lines(geoborders$x, geoborders$y, col = "black")
        ## lines(attr(geoborders, "borders")$x, attr(geoborders, "borders")$y, col = "grey90")
    }
}

col.bar <- function(breaks,horiz=TRUE,pch=21,v=1,h=1,col=col,cex=2,cex.lab=0.6,type="r",verbose=FALSE,vl=0.5,border=FALSE,...) {
    par0 <- par()
    xleft <- par()$usr[1] 
    xright <- par()$usr[2]
    ybottom <- par()$usr[4] - 1 - h
    ytop <-  par()$usr[4] - 1 
    
    by <- (xright - xleft - v * (length(col)))/(length(breaks))
    steps <-   seq(0, (xright -xleft - v * (length(col))) ,by=by ) # 
    nsteps <- length(steps) 
    
    if (verbose) print(steps)
    if (verbose) print(breaks)
    if (verbose) print(nsteps)
    
    ## if (max(abs(breaks))<=1) breaks <- round(breaks,digits=2)
    
    k <- 1/2
    for (i in 1 :(nsteps-2)) {  
        if (!is.null(v)) 
            if (i == 1) k <- k + v/2 else k <- k + v  
        if (type == "r") { ## "r" for rectangle
            rect(xleft= k  + xleft + steps[i] ,xright= k + xleft + steps[i+1],ybottom=ybottom,ytop=ytop,col=col[i],border=border)
            
            ## text(x = k + xleft + steps[i], y = ybottom - 1,labels=sprintf("%.1f",icn[i]),cex=cex)
        }
        else if (type == "p") { ## "p" points
            points(x= k + xleft + (steps[i]+ steps[i+1])/2, y=(ybottom + ytop)/2,pch=pch, bg=col[i],cex=cex,...)
            
        }
        
        text(x = k + xleft + (steps[i]+ steps[i+1])/2,  y = ybottom - vl, labels=levels(cut(breaks,breaks))[i],col="grey50",cex=cex.lab)
    } 
    par(fig=par0$fig)
}

#trend <- function(x,ns.omit=TRUE,alpha=0.1) {
#    t <- 1:length(x)
#    model <- lm(x ~ t)
#    y <- c(model$coefficients[2]*10)
#    if (ns.omit) if (anova(model)$Pr[1] > alpha) y <- NA
#    names(y) <- c("coefficients")
#    return(y)
#}

colbar2 <- function(x,col) {
    par(mar=c(1,0,0,0),fig=c(0.5,1,0.665,0.695),new=TRUE,cex.axis=0.6)
    nl <- pretty(x)
    n <- length(nl)
    image(cbind(1:n,1:n),col=col) 
    par(xaxt="s",new=new)
    axis(1,at=seq(0,1,length=length(nl)),label=nl)
}


sphere <- function(x,n=30,FUN="mean",lonR=10,latR=45,axiR=0,
                         gridlines=TRUE,col="green",bg="darkgreen",cex=0.2,pch=".",new=TRUE) {
  x0 <- x
  ## Data to be plotted:
  if (inherits(x,"stationmeta")) {
      lon <- x$longitude
      lat <- x$latitude
      param <- param2ele(x$ele)
      unit <- " "
  }
  else if (inherits(x,"station")) {
      lon <- attr(x,'longitude')
      lat <- attr(x,'latitude')
      param <- as.character(levels(factor(attr(x,'parameter'))))
  }
  ## To deal with grid-conventions going from north-to-south or east-to-west:
  ##srtx <- order(attr(x,'longitude')); lon <- lon[srtx]
  ##srty <- order(attr(x,'latitude')); lat <- lat[srty]
 
  if (!is.null(FUN))
      map <- apply(as.matrix(x),2,FUN,na.rm=TRUE) ##map <- x[srtx,srty]
  else
      map <- x
  
  
  # Rotatio:
  if (is.null(lonR)) lonR <- mean(lon)  # logitudinal rotation
  if (is.null(latR)) latR <- mean(lat)  # Latitudinal rotation
  # axiR: rotation of Earth's axis

  # coastline data:
  data("geoborders",envir=environment())
  ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  theta <- pi*geoborders$x[ok]/180; phi <- pi*geoborders$y[ok]/180
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

  ## Define colour palette:
  if (!is.null(FUN)) {
      if (is.null(col)) col <- colscal(n=n) else
      if (length(col)==1) {
          palette <- col
          col <- colscal(palette=palette,n=n)
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
  
# Plot the results:  
  if (new) dev.new()
  par(bty="n",xaxt="n",yaxt="n",new=TRUE)
  plot(x,z,pch=".",col="white",xlab="",ylab="")
 
# plot the grid boxes, but only the gridboxes facing the view point:
  ##Visible <- Y > 0 ##colMeans(Y) > 0
  ##X <- X[,Visible]; Y <- Y[,Visible]; Z <- Z[,Visible]
  ##index <- index[Visible]
  ##apply(rbind(X,Z,index),2,gridbox,cols)
  # c(W,E,S,N, colour)
  # xleft, ybottom, xright, ytop
  
  if (!is.null(FUN)) {
      breaks <- pretty(map,n=n)
      colb <- colscal(n=length(breaks)) 
      col <- colb[findInterval(map,breaks)]
      bg <- col
      nc <- length(colb)
  }
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
      image(seq(breaks[1],breaks[length(breaks)],length=nc),c(1,2),bar,col=col)
      
      par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
          fig=c(0,1,0,1),new=TRUE)
      plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
      text(0.1,0.95,param,cex=1.5,pos=4)
      text(0.72,0.002,unit,pos=4)  
  }
  
  ##result <- data.frame(x=colMeans(Y),y=colMeans(Z),z=c(map))
  if (inherits(x0,"stationmeta"))
      result <- data.frame(x=Y,y=Z)
  else if (inherits(x0,"station"))
      result <- data.frame(x=Y,y=Z,z=map)
  invisible(result)
}
