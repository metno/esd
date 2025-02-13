#' Display a color bar object on an existing plot.
#' 
#' Add a color bar or color points into an exisiting plot or map.
#' For a description of palette choices, see \code{\link{colscal}}.
#' 
#' @aliases col.bar colbar
#' @seealso colbar.ini colscal
#'
#' @importFrom grDevices heat.colors rainbow rgb terrain.colors gray.colors topo.colors col2rgb cm.colors adjustcolor
#' 
#' @param breaks A numeric vector of breakpoints for the colours
#' @param horiz a boolean; if TRUE add horizontal color bar, else add vertical color bar 
#' @param pch see \code{\link{par}} 
#' @param v Vertical space between color bar points
#' @param h horizontal space between color bar points
#' @param col see \code{\link{par}}
#' @param cex A numerical value giving the amount by which plotting text
#'        and symbols should be magnified relative to the default (see \code{\link{par}})
#' @param cex.lab Magnification factor for x and y labels (see \code{\link{par}})
#' @param cex.axis Magnification factor for axis annotations (see \code{\link{par}})
#' @param type r : rectangular shape , p : for points
#' @param verbose a boolean; if TRUE print information about progress
#' @param vl a numerical specifying the relative placement of the vertical lines
#' @param border a boolean; if TRUE show color bar borders
#' @param \dots Additional graphical parameters to be passed on
#' 
#' @export
col.bar <- function(xleft,ybottom,xright,ytop,breaks,horiz=TRUE,
                    pch=15,v=1,h=1,col=col,cex=5,cex.lab=0.6,
                    cex.axis=0.9,type="r",verbose=FALSE,vl=0.5,border="black",...) {
  if (verbose) print('col.bar')
  # par0 <- par(no.readonly=TRUE)
  # xleft <- par0$usr[1] 
  # xright <- par0$usr[2]
  # ybottom <- par0$usr[4] - 1 - h
  # ytop <-  par0$usr[4] - 1 
  
  # by <- (xright - xleft - v * (length(col)))/(length(breaks))
  # steps <-   seq(0, (xright -xleft - v * (length(col))) ,by=by ) # 
  # nsteps <- length(steps) 
  # 
  # if (verbose) print(steps)
  # if (verbose) print(breaks)
  # if (verbose) print(nsteps)
  # 
  # k <- 1/2
  # for (i in 1 :(nsteps-2)) {  
  #   if (!is.null(v)) 
  #     if (i == 1) k <- k + v/2 else k <- k + v  
  #     if (type == "r") { ## "r" for rectangle
  #       rect(xleft= k  + xleft + steps[i] ,xright= k + xleft + steps[i+1],
  #            ybottom=ybottom,ytop=ytop,col=col[i],border=border)
  #     } else if (type == "p") { ## "p" points
  #       points(x= k + xleft + (steps[i]+ steps[i+1])/2, y=(ybottom + ytop)/2,
  #              pch=pch, bg=col[i],cex=cex,...)
  #       
  #     }        
  #     text(x = k + xleft + (steps[i]+ steps[i+1])/2,  y = ybottom - vl,
  #          labels=levels(cut(breaks,breaks))[i],col="grey50",cex=cex.lab)
  # } 
  
  ymid <- 0.5*(ybottom + ytop)
  n <- length(breaks)
  dx <- 0.1*(xright - xleft)
  dy <- 0.1*(ytop - ybottom)
  mids <- seq(xleft+dx,xright-dx,length=length(col)+1)
  
  ## Adjust tick labels
  dm <- diff(mids)[1]*0.5
  db <- dy*0.75
  
  #points(mids,rep(ymid,n-1),col=col,pch=pch,cex=cex)
  #image(0.9*mids,c(ymid,ytop)+c(dy,-dy),cbind(breaks,breaks),col=col,ylim=c(ybottom,ytop),add=TRUE)
  image(mids,c(ymid,ytop)+c(dy,-dy),cbind(breaks,breaks),col=col,ylim=c(ybottom,ytop),add=TRUE)
  #rect(min(mids),ymid,max(mids),ytop,border="black")
  if(n<=11) {
    ii <- rep(TRUE, n)
  } else if(n<=21) {
    ii <- (1:n)%%2 == 1
  } else ii <- (1:n)%%round(n/10) == 1
  #text(mids[ii],rep(ybottom+dy,n)[ii],round(breaks,2)[ii],cex=cex.axis, col='grey30')
  text(mids[ii]+dm, rep(ybottom,n)[ii]+db, round(breaks,2)[ii],cex=cex.axis, col='grey30')
  invisible(list(mids=mids,col=col,breaks=breaks))
}

#' @export
colbar <- function(breaks,col,fig=c(0.15,0.2,0.15,0.3),horiz=FALSE,
                   mar=c(1,0,0,0),new=TRUE,las=1,cex.axis=0.6,...) {
  par0 <- par(no.readonly = TRUE) # save default, for resetting...
  mids <- (breaks[1:(length(breaks)-1)] + breaks[2:length(breaks)])/2
  if (horiz) {
    par(xaxt="s",yaxt="n",fig=fig,mar=mar,new=new,las=las,cex.axis=cex.axis,...)
    #image(breaks,1:2,cbind(breaks,breaks),col=col,cex.axis=cex.axis)
    image(mids,1:2,cbind(mids,mids),col=col,cex.axis=cex.axis,axes=FALSE)
    axis(1, at=breaks, labels=breaks)
  } else {
    par(xaxt="n",yaxt="s",fig=fig,mar=mar,new=new,las=las,cex.axis=cex.axis,...)
    #image(1:2,breaks,rbind(breaks,breaks),col=col,cex.axis=cex.axis)
    image(1:2,mids,rbind(mids,mids),col=col,cex.axis=cex.axis,axes=FALSE)
    axis(4, at=breaks, labels=breaks)
  }
  par(fig=par0$fig, xaxt=par0$xaxt, yaxt=par0$yaxt, 
      mar=par0$mar, las=par0$las, cex.axis=par0$cex.axis)
  invisible(list(mids=mids,col=col,breaks=breaks))
}

#' Display a color bar object on an existing plot.
#' 
#' Generate a color bar list and add information about the breaks of the color scale based on the
#' numerical range of the input data.
#' 
#' @seealso col.bar colbar
#' 
#' @param x an input object, e.g., a 'zoo', 'station' or 'field' object or numerical vector
#' @param FUN a function 
#' @param verbose a boolean; if TRUE print information about progress
#' @param colbar a list: colbar = list(col, breaks, n, type, cex, h, v, pos, show, rev)
#' where 
#' \code{col} is a vector containing the colors corresponding to the values 
#' specified in the numerical vector \code{breaks},
#' \code{n} is the number of breaks (used only if breaks are not specified),
#' \code{show} if TRUE show color bar,
#' \code{rev} if TRUE reverse color scale,
#' \code{cex} see \code{\link{par}},
#' \code{h}, \code{v}, \code{type}: see \code{\link{col.bar}}
#' \code{pos} not in use?
#'
#' @export
colbar.ini <- function(x,FUN=NULL,colbar=NULL,verbose=FALSE) {
  if (verbose) {print('colbar.ini'); print(colbar)}
  par0 <- par(no.readonly = TRUE) # save default, for resetting...
  if (length(x)==0) stop('colbar.ini: x is empty!')
  if (is.null(colbar)) colbar <- list(show=FALSE,n=14,rev=NULL,alpha=NULL)
  if (is.logical(colbar)) colbar <- list(show=colbar)
  if (verbose) print('sort out the colours')
  
  ## Prepare data and calculate range
  if (is.zoo(x)) x <- coredata(x)
  x[!is.finite(x)] <- NA      # REB 2017-09-20: fix to cope with Inf-values
  x.rng <- range(x,na.rm=TRUE)
  if (verbose) {print('Value range:'); print(x.rng)}
  ## If there are bad range values
  if (!is.finite(x.rng[1])) {
    ## If only the first is bad: set to 0 or a value lower than 2nd (negative)
    if (is.finite(x.rng[2])) {
      x.rng[1] <- min(0,x.rng[2]*2) 
    } else {
      x.rng <- c(0,1)
    }
  }
  if (!is.finite(x.rng[2])) {
    ## If only the first is bad: set to 0 or a value higher than 2nd
    if (is.finite(x.rng[1])) {
      x.rng[2] <- max(0,x.rng[1]*2) 
    } else {
      x.rng <- c(0,1)
    }
  }
  if (verbose) print(x.rng)
  nd <- max(0,ndig(x.rng)+2)
  
  ## Set breaks and n
  if (verbose) {print('Set breaks and n:')}
  if (!is.null(colbar$col) & is.null(colbar$breaks)) {
    colbar$n <- length(colbar$col)
    if (is.null(colbar$breaks)) {
      colbar$breaks <- round(seq(x.rng[1],x.rng[2],length.out=length(colbar$col)+1),nd)
    } else if(length(colbar$breaks)!=(colbar$n-1)) {
      colbar$breaks <- pretty(colbar, n=colbar$n-1)
    }
  }
  if (is.null(colbar$breaks)) { 
    if (verbose) print("pretty is used here to set break values ...")
    if (!is.null(colbar$n)) {
      colbar$breaks <- pretty(seq(x.rng[1],x.rng[2],length.out=colbar$n+1),n=colbar$n+1)
    } else {
      colbar$breaks <- pretty(seq(x.rng[1],x.rng[2],length.out=10),n=11)
    }
  } else if(length(colbar$breaks)==2) {
    colbar$breaks <-  pretty(seq(colbar$breaks[1],colbar$breaks[2],length.out=10),n=11)
  }
  colbar$n <- length(colbar$breaks) - 1
  
  ## Activate pallette (pal)
  if (is.null(colbar$pal)) {
    if (is.precip(x)) {
      colbar$pal <- 'precip'
    } else {
      colbar$pal <- 't2m'
    } 
  } 
  
  ## Specify other colbar stuff
  if (is.null(colbar$n)) colbar$n <- length(colbar$breaks) - 1
  if (is.null(colbar$type)) colbar$type <- 'p'
  if (is.null(colbar$cex)) colbar$cex <- 2
  if (is.null(colbar$h)) colbar$h <- 0.6
  if (is.null(colbar$v)) colbar$v <- 1
  if (is.null(colbar$pos)) colbar$pos <- 0.05
  if (is.null(colbar$show)) colbar$show <-TRUE
  if (is.null(colbar$rev)) colbar$rev <- FALSE
  ## Check and define or correct colbar$col
  if (!is.null(colbar$col)) {
    if (is.null(colbar$pal)) colbar$pal <- NA
    if (!is.null(colbar$breaks)) {  
      if (length(colbar$col) != length(colbar$breaks) - 1) {
        # if col and breaks are specified but not consistent, interpolate col to right length:
        col.rgb <- col2rgb(colbar$col)
        col.rgb <- apply(col.rgb,1,function(x) approx(x,
                                                      n=length(colbar$breaks)-1)$y)
        if(is.null(dim(col.rgb))) dim(col.rgb) <- c(1,length(col.rgb))
        colbar$col <- rgb(col.rgb,maxColorValue=255)
      }
    } else if (is.null(colbar$breaks)) {
      colbar$breaks <- round(seq(x.rng[1],x.rng[2],length.out=colbar$n+1),nd)   
    }
  } else {
    if (verbose) print('define col')
    if (verbose) print(paste('colbar$n',colbar$n))
    colbar$col <- colscal(n=colbar$n,pal=colbar$pal,alpha=colbar$alpha,
                          rev=colbar$rev,verbose=verbose)
  }
  #if (verbose) print(colbar)
  if (verbose) print(paste("length(col) =",length(colbar$col),
                           "length(breaks) =",length(colbar$breaks)))
  
  if (length(colbar$col) != length(colbar$breaks)-1) stop('colbar.ini: Error in setting colbar!')
  
  if (verbose) {
    #print(colbar)
    print('exit colbar.ini')
  }
  invisible(colbar)
}

## KMP 2024-02-28: Some of the color scales that span from one color to another are out of balance
## Is this just a feature of the scale, or should the scale be adjusted?
## bwr: the white color is not centered but towards the blue side 
## burd, t2m.kin2100: the blue part of the scale is shorter than the red

#' Generate a color scale
#'
#' Generate a vector of colors
#'
#' Palette options are as follows:
#' 't2m': blue-yellow-red (from seNorge), 
#' 'precip', 'mu' and 'fw': white-blue (from seNorge),
#' 'cold': a cold color scale (from seNorge),
#' 'warm': a warm color scale (from seNorge),
#' 't2m.IPCC': blue-red color scale from the IPCC visual guide for authors from 2019
#' 'precip.IPCC': green-brown colros scale from the IPCC visual guide for authors from 2019
#' 't2m.kin2100': blue-red color scale from Klima i Norge 2100 (KiN2100)
#' 'precip.kin2100': brown-green-blue color scale (from KiN2100)
#' 'precip.trend.kin2100': green-white-blue color scale (from KiN2100)
#' 'cold.kin2100': a cold color scale (from KiN2100),
#' 'warm.kin2100': a warm color scale (from KiN2100),
#' 'bwr' (blue-white-red),
#' 'slp' and 'mslp' (same as 'bwr'),
#' 'rwb' (red-white-blue), 
#' 'faint.rwb': fainter version of 'rwb',
#' 'faint.bwr': fainter version of 'bwr', 
#' 'grmg': green-magenta, 
#' 'brbu': brown-blue,
#' 'budor': blue-orange,
#' 'budrd': blue-red,
#' 'bugr': blue-green,
#' 'bugy': blue-gray,
#' 'buor': blue-orange (brighter colors than 'budor'),
#' 'buorr': blue-green (brighter and more yellow than 'buor' and 'budor'),
#' 'bu': blues,
#' 'rd': reds,
#' 'brgrbu': brown-green-blue, same as precip.kin2100
#' 'grwbu': green-white-blue, same as precip.trend.kin2100
#' 'burd': blue-red, same as t2m.kin2100
#' 'buyrd': blue-yellow-red, same as t2m
#' 'cat': categorical color scale,
#'  color scales from grDevices: 'rainbow', 'gray.colors', 'heat.colors', 'terrain.colors',
#' 'topo.colors', and 'cm.colors'.
#'
#' @seealso colbar col.bar colbar.ini
#'
#' @param n length of color vector
#' @param pal color palette: "bwr",rwb","faint.bwr","faint.rwb","rainbow","gray.colors","heat.colors",
#' "terrain.colors","topo.colors","cm.colors","grmg","brbu","budor","budrd","bugr","bugy","buor",
#' "buorr","bu","rd","brgrbu","grwbu","burd","buyrd","cat","cold", "warm", "t2m", "precip", "fw", "mu", 
#' "precip.kin2100", "precip.trend.kin2100", "t2m.kin2100"
#' @param rev a boolean; if TRUE reverse color scale
#' @param alpha factor defining transparency of color
#' @param test a boolean; if TRUE show a sample of the color scale
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @export
colscal <- function(n=14,pal="t2m",rev=FALSE,alpha=NULL,test=FALSE,verbose=FALSE) {
  
  test.col <- function(r,g,b) {
    #dev.new()
    par(bty="n")
    plot(r,col="red")
    points(b,col="blue")
    points(g,col="green")
  }
  
  if (verbose) print(paste('colscal:',pal,'rev=',rev,'n=',n,'alpha=',alpha))
  par0 <- par(no.readonly = TRUE) # save default, for resetting...
  if ( (is.null(pal)) | (is.na(pal)) ) pal <- 't2m'
  if (is.null(alpha)) alpha <- 1
  # Set up colour-palette
  pal <- tolower(gsub("[[:space:]]", "", pal))#[[:punct:][:space:]]", "", pal))
  x <- 1:n
  r0 <- round(n*0.55)
  g0 <- round(n*0.5)
  b0 <- round(n*0.45)
  s <- -0.1/n
  if (n < 30) sg <- s*2.5 else sg <- s
  n1 <- g0; n2 <- n-n1
  ## Palettes for seNorge
  #R	G	B
  t2m.seNorge <- c(0,   0, 153,
                0,  25, 255,
                0, 153, 255,
                64, 204, 255,
                128, 235, 255,
                179, 255, 255,
                217, 255, 255,
                255, 255, 190,
                255, 255,  64,
                255, 230,  77,
                255, 179,   0,
                255, 102,   0,
                255, 25,    0,	
                204,  0,    0)	
  dim(t2m.seNorge) <- c(3,14)
  
  precip.seNorge <- c(229, 229, 229,
                217, 255, 255,
                179, 255, 255,
                128, 235, 255,
                64, 204, 255,
                0, 153, 255,
                0, 25, 255,
                0, 0, 153)
  dim(precip.seNorge) <- c(3,8)

  ## Palettes for Klima i Norge 2100
  t2m.kin2100 <- c(15,	0,	5,
            66,	1,	20,
            103,	0,	31,
            178,	24,	43,
            214,	96,	77,
            244,	165,	130,
            253,	219,	199,
            247,	247,	247,
            209,	229,	240,
            146,	197,	222,
            67,	147,	195,
            33,	102,	172,
            5,	48,	97)
  dim(t2m.kin2100) <- c(3,13)
  
  precip.kin2100 <- c(69,	80,	138,
            64,	106,	168,
            52,	132,	201,
            26,	160,	237,
            106,	187,	217,
            156,	214,	193,
            181,	222,	159,
            199,	224,	123,
            255,	236,	191,
            230,	201,	145,
            204,	170,	102,
            143,	109,	40,
            110,	79,	18)
  dim(precip.kin2100) <- c(3,13)
  
  dpr.kin2100 <- c(84,	48,	5,
             140,	81,	10,
             191,	129,	45,
             223,	194,	125,
             246,	232,	195,
             245,	245,	245,
             199,	234,	229,
             128,	205,	193,
             53,	151,	143,
             1,	102,	94,
             0,	60,	48)
  dim(dpr.kin2100) <- c(3,11)

  # https://www.ipcc.ch/site/assets/uploads/2019/04/IPCC-visual-style-guide.pdf
  t2m.IPCC <- c(103, 0, 31,
                178, 24, 43,
                214, 96, 77,
                244, 165, 130,
                253, 219, 199,
                247, 247, 247,
                209, 229, 240,
                146, 197, 222,
                67, 147, 195,
                33, 102, 172,
                5, 48, 97)
  dim(t2m.IPCC) <- c(3,11)
  
  precip.IPCC <- c(84, 48, 5,
                   140, 81, 10,
                   191, 129, 45,
                   223, 194, 125,
                   246, 232, 195,
                   245, 245, 245,
                   199, 234, 229,
                   128, 205, 193,
                   53, 151, 143,
                   1, 102, 94,
                   0, 60, 48)
  dim(precip.IPCC) <- c(3,11)
  
  if (!is.null(alpha)) alpha <- rep(alpha[1],n)
  if (tolower(pal[1])=='t2m.ipcc') {
    r <- approx(t2m.IPCC[1,],n=n)$y/255
    g <- approx(t2m.IPCC[2,],n=n)$y/255
    b <- approx(t2m.IPCC[3,],n=n)$y/255
    col <- rev(rgb(r,g,b,alpha))
  } else if (tolower(pal[1])=='precip.ipcc') {
    r <- approx(precip.IPCC[1,],n=n)$y/255
    g <- approx(precip.IPCC[2,],n=n)$y/255
    b <- approx(precip.IPCC[3,],n=n)$y/255
    col <- rgb(r,g,b,alpha)
  } else if ( (pal[1]=="bwr") | (pal[1]=="slp") | (pal[1]=="mslp") |
       (pal[1]=="pressure") ) {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
    g <- exp(sg*(x - g0)^2)^2
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
    col <- rgb(r,g,b,alpha)
  } else if (pal[1]=="rwb") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
    g <- exp(sg*(x - g0)^2)^2
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
    if (is.null(alpha)) col <- rgb(b,g,r)  else
      col <- rgb(r,g,b,alpha)
  } else if (pal[1]=="faint.bwr") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
    g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
    col <- rgb(r,g,b,alpha)
  } else if (pal[1]=="faint.rwb") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
    g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
    col <- rgb(r,g,b,alpha)
  } else if ( (pal[1]=="precip") | (pal[1]=="mu") | (pal[1]=="fw") |
              (pal[1]=="f[w]") | (pal[1]=="tp") | (pal[1]=='rr') | (pal[1]=='prate') ) {
    r <- approx(precip.seNorge[1,],n=n)$y/255
    g <- approx(precip.seNorge[2,],n=n)$y/255
    b <- approx(precip.seNorge[3,],n=n)$y/255
    col <- rgb(r,g,b,alpha)
  } else if (pal[1] %in% c("precip.kin2100","brgrbu")) {
    r <- approx(precip.kin2100[1,],n=n)$y/255
    g <- approx(precip.kin2100[2,],n=n)$y/255
    b <- approx(precip.kin2100[3,],n=n)$y/255
    col <- rev(rgb(r,g,b,alpha))
  } else if (pal[1] %in% c("precip.trend.kin2100","grwbu")) {
    r <- approx(dpr.kin2100[1,],n=n)$y/255
    g <- approx(dpr.kin2100[2,],n=n)$y/255
    b <- approx(dpr.kin2100[3,],n=n)$y/255
    col <- rgb(r,g,b,alpha)
  } else if (pal[1]=="rainbow") {
    col <- rainbow(n,start=0,end=4/6,alpha=alpha[1])
  } else if (pal[1]=="gray.colors") {
    col <- gray.colors(n,alpha=alpha[1])
  } else if (pal[1]=="heat.colors") {
    col <- heat.colors(n,alpha=alpha[1])
  } else if (pal[1]=="terrain.colors") {
    col <- terrain.colors(n,alpha=alpha[1])
  } else if (pal[1]=="topo.colors") {
    col <- topo.colors(n,alpha=alpha[1])
  } else if (pal[1]=="cm.colors") {
    col <- cm.colors(n,alpha=alpha[1])
  } else if (pal[1]==tolower("grmg")) {
    cols <- list(
      r=c(0,0,0,0,0.316,0.526,0.737,1,1,1,1,1,0.947,0.737,0.526,0.316),
      g=c(0.316,0.526,0.737,0.947,1,1,1,1,0.947,0.737,0.526,0.316,0,0,0,0),
      b=c(0,0,0,0,0.316,0.526,0.737,1,1,1,1,1,0.947,0.737,0.526,0.316))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]==tolower("brbu")) {
    cols <- list(
      r=c(0.2,0.4,0.6,0.8,0.85,0.95,0.8,0.6,0.4,0.2,0,0),
      g=c(0.1,0.187,0.379,0.608,0.688,0.855,0.993,0.973,0.94,0.893,0.667,0.48),
      b=c(0,0,0.21,0.480,0.595,0.808,1,1,1,1,0.8,0.6))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]==tolower("budor")) {
    cols <- list(
      r=c(0.12,0.32,0.6,0.7,0.8,0.9,1,1,1,1,0.8,0.6),
      g=c(0.56,0.768,0.98,0.99,0.997,1,0.9,0.793,0.68,0.56,0.347,0.250),
      b=c(0.6,0.8,1,1,1,1,0.8,0.6,0.4,0.2,0,0))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]==tolower("budrd")) {
    cols <- list(
      r=c(0.142,0.097,0.16,0.24,0.34,0.46,0.6,0.74,0.92,1,
          1,1,1,1,1,0.97,0.85,0.65),
      g=c(0,0.112,0.342,0.531,0.692,0.829,0.920,0.978,1,1,
          0.948,0.840,0.676,0.472,0.240,0.155,0.085,0),
      b=c(0.85,0.97,1,1,1,1,1,1,1,0.92,0.74,0.6,0.46,0.34,0.24,0.21,0.187,0.13))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]==tolower("bugr")) {
    cols <- list(
      r=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.9,0.8,0.7,0.6,0.4,0.2,0),
      g=c(0,0.2,0.4,0.6,0.7,0.8,0.9,1,1,1,1,1,1,1),
      b=c(1,1,1,1,1,1,1,0.9,0.8,0.7,0.6,0.4,0.2,0))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]==tolower("bugy")) {
    cols <- list(
      r=c(0,0.4,0.6,0.8,0.9,0.6,0.4,0.2),
      g=c(0.6,0.9,1,1,0.9,0.6,0.4,0.2),
      b=c(0.8,1,1,1,0.9,0.6,0.4,0.2))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]==tolower("buor")) {
    cols <- list(
      r=c(0,0.1,0.2,0.4,0.6,0.8,1,1,1,1,1,1),
      g=c(0.167,0.4,0.6,0.8,0.933,1,1,0.933,0.8,0.6,0.4,0.167),
      b=c(1,1,1,1,1,1,0.8,0.6,0.4,0.2,0.1,0))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]==tolower("buorr")) {
    cols <- list(
      r=c(0.03,0.2,0.35,0.55,0.75,0.9,0.97,1,1,1,1,1,1,1),
      g=c(0.353,0.467,0.567,0.7,0.833,0.933,0.98,1,1,1,0.8,0.6,0.4,0),
      b=c(1,1,1,1,1,1,1,0.8,0.6,0,0,0,0,0))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]=="bu") {
    cols <- list(
      r=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0),
      g=c(1,0.983,0.95,0.9,0.833,0.75,0.65,0.533,0.4,0.250),
      b=c(1,1,1,1,1,1,1,1,1,1))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]=="rd") {
    cols <- list(
      r=c(1,1,1,1,1,1,1,1,1,1),
      g=c(1,0.983,0.95,0.9,0.833,0.75,0.65,0.533,0.4,0.250),
      b=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]=="cat") {
    cols <- list(
      r=c(1,1,1,1,0.7,0.2,0.65,0.1,0.8,0.4,1,0.9),
      g=c(0.75,0.5,1,1,1,1,0.93,0.7,0.75,0.3,0.6,0.10),
      b=c(0.5,0,0.6,0.2,0.55,0,1,1,1,1,0.75,0.2))
    cols <- lapply(cols,function(x) approx(x,n=n)$y)
    col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (pal[1]=="cold") {
    r <- approx(t2m.seNorge[1,1:7],n=n)$y/255
    g <- approx(t2m.seNorge[2,1:7],n=n)$y/255
    b <- approx(t2m.seNorge[3,1:7],n=n)$y/255
    col <- rgb(r,g,b,alpha)    
  }  else if (pal[1]=="warm") {
    r <- approx(t2m.seNorge[1,8:14],n=n)$y/255
    g <- approx(t2m.seNorge[2,8:14],n=n)$y/255
    b <- approx(t2m.seNorge[3,8:14],n=n)$y/255
    col <- rgb(r,g,b,alpha)    
  } else if(pal[1]=="t2m") {
    r <- approx(t2m.seNorge[1,],n=n)$y/255
    g <- approx(t2m.seNorge[2,],n=n)$y/255
    b <- approx(t2m.seNorge[3,],n=n)$y/255
    col <- rgb(r,g,b,alpha)
  } else if(pal[1] %in% c("t2m.kin2100","burd")) {
    r <- approx(t2m.kin2100[1,],n=n)$y/255
    g <- approx(t2m.kin2100[2,],n=n)$y/255
    b <- approx(t2m.kin2100[3,],n=n)$y/255
    col <- rev(rgb(r,g,b,alpha))
  } else if (pal[1]=="warm.kin2100") {
    r <- approx(t2m.kin2100[1,1:7],n=n)$y/255
    g <- approx(t2m.kin2100[2,1:7],n=n)$y/255
    b <- approx(t2m.kin2100[3,1:7],n=n)$y/255
    col <- rgb(r,g,b,alpha)    
  }  else if (pal[1]=="cold.kin2100") {
    r <- approx(t2m.kin2100[1,8:13],n=n)$y/255
    g <- approx(t2m.kin2100[2,8:13],n=n)$y/255
    b <- approx(t2m.kin2100[3,8:13],n=n)$y/255
    col <- rgb(r,g,b,alpha)    
  } else {
    r <- approx(t2m.seNorge[1,],n=n)$y/255
    g <- approx(t2m.seNorge[2,],n=n)$y/255
    b <- approx(t2m.seNorge[3,],n=n)$y/255
    col <- rgb(r,g,b,alpha)
  }
  if (test) { #& !exists("r")) {
    RGB <- col2rgb(col)/255
    r <- RGB[1,]; g <- RGB[2,]; b <- RGB[3,]
  }
  
  if (test) test.col(r,g,b)
  if (rev) col <- rev(col)
  return(col)
}

