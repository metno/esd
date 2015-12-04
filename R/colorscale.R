
##seNorge
#nodata	10000
#nobits	16
#header1	*Temperatur
#header2	*Siste døgn (18-18 UTC)
#legend	*Grader Celsius
#From	To	R	G	B	Forklaring
#2931	10000	204	0	0	*Over 20 
#2881	2931	255	25	0	*15 - 20
#2831	2881	255	102	0	*10 - 15
#2781	2831	255	179	0	*5 - 10
#2761	2781	255	230	77	*3 - 5
#2741	2761	255	255	64	*1 - 3
#2731	2741	255	255	190	*0 - 1
#2721	2731	217	255	255	*÷1 - 0
#2701	2721	179	255	255	*÷3 - ÷1
#2681	2701	128	235	255	*÷5 - ÷3
#2631	2681	64	204	255	*÷10 - ÷5
#2581	2631	0 	153	255	*÷15 - ÷10
#2531	2581	0 	25	255	*÷20 - ÷15
#0	2531	0	0	153	*Under ÷20


#nodata	10000
#nobits	16
#header1	*Nedbør
#header2	*Siste døgn (06-06 UTC)
#legend	*mm
#From	To	R	G	B	Forklaring
#1500	10000	0	0	153	*Over 150 
#750	1500	0 	25	255	*75 - 150 
#500	750	0 	153	255	*50 - 75
#300	500	64	204	255	*30 - 50
#200	300	128	235	255	*20 - 30
#100	200	179	255	255	*10 - 20
#1	100	217	255	255	*Under 10
#0	1	229	229	229	*Ikke nedbør


colbar.ini <- function(x,FUN=NULL,colbar=NULL,verbose=TRUE) {

    ## browser()
    if (is.logical(colbar)) return(NULL)
    ##if (!is.null(colbar)) {
    if (verbose) print('sort out the colours')

    ##if (!is.null(colbar$col)) {
    ##    colbar$n <- length(colbar$col) + 1
    ##    colbar$pal <- NULL
    ##}   
   
    if (is.null(colbar$n))
        if (!is.null(colbar$col))
            colbar$n <- length(colbar$col)
        else
            colbar$n <- 10

    if (is.zoo(x)) x <- coredata(x)
    x.rng <- range(x,na.rm=TRUE)

    ## very easy case if colbar$col and breaks are provided
    if (!is.null(colbar$col)) {
        pal <- NULL ## desactivate pal
        if (!is.null(colbar$breaks)) {  
            if (length(colbar$col) != length(colbar$breaks) - 1)
                stop('Length of breaks must be the lenght of color + 1')   
        } else colbar$breaks <- seq(x.rng[1],x.rng[2],length.out=colbar$n+1)
        ## if only colbar$col is provided, then the breaks are set using seq   
    }

    ## if breaks are null then use pretty
    if (is.null(colbar$breaks)) { 
        if (verbose) print("pretty is used here to set break values ...")
        if (!is.null(colbar$n))
            colbar$breaks <- pretty(seq(x.rng[1],x.rng[2],length.out=colbar$n))
        else
            colbar$breaks <- pretty(seq(x.rng[1],x.rng[2]))
    }
    
    if (is.null(colbar$type)) colbar$type <- 'p'

    if (is.null(colbar$cex)) colbar$cex <- 2

    if (is.null(colbar$h)) colbar$h <- 0.6

    if (is.null(colbar$v)) colbar$v <- 1

    if (is.null(colbar$pos)) colbar$pos <- 0.05

    if (is.null(colbar$show)) colbar$show <-TRUE

    if (is.null(colbar$rev)) colbar$rev <- FALSE
    
    ## activate pallette (pal)
    if (is.null(colbar$pal)) {
      colbar$pal <- varid(x)[1]
      if (!is.precip(x)) {
        colbar$pal <- 't2m'
        if (is.null(colbar$rev)) colbar$rev <- FALSE
      } else  {
        colbar$pal <- 'precip'
        if (is.null(colbar$rev)) colbar$rev <- TRUE
      } 
    }
# REB 2015-12-02: I do not understand these two lines    
#    if (!is.null(FUN)) {
#        if (is.null(colbar$breaks) & !inherits(x,"stationmeta")) {
#            colbar$breaks <- pretty(x,n=colbar$n)
# Replaced by the following line:
    if (is.null(colbar$breaks)) {
      ## If colbar$breaks is unspecified, then set up the levels for colour scale:
      if (verbose) print('define breaks')
      if (is.null(colbar$col)) {
        ## If colbar$col is unspecified, then use pretty to provide pretty numbers
            colbar$breaks <- pretty(x,n=colbar$n)
      } else {
        ## If colbar$col *is* specified, then the numbers are given
        colbar$n <- length(colbar$col)
        colbar$breaks <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=colbar$n+1)
      }
    } else if (length(colbar$breaks)==2) {
      ## If colbar is a vector of two, then it is taken as a range
      if (!is.null(colbar$col)) colbar$n <- length(colbar$col) else
      if (!is.null(colbar$n)) colbar$n <- 11 
      colbar$breaks <- seq(colbar$breaks[1],colbar$breaks[2],
                           length=colbar$n+1)
    }
    colbar$n <- length(colbar$breaks) -1

    
    ## if colbar$col is null
    if (is.null(colbar$col)) {
      if (verbose) print('define col')

        ## colscal is used as default to set the colors
      if (verbose) print(paste('colbar$n',colbar$n))
        colbar$col <- colscal(n=colbar$n,col=colbar$pal,
                              rev=colbar$rev,verbose=verbose)
      }
    if (verbose) print(colbar)
    
    ##    if (verbose) print(paste("length(col) =",length(colbar$col)))
    ##    col <- colscal(n=colbar$n,col=colbar$pal,rev=colbar$rev)       
    
    ## if (!inherits(x,"stationmeta"))
    ##     colbar$col <- colscal(n=colbar$n,col=colbar$pal,rev=colbar$rev,verbose=verbose)
    if (verbose) print(paste("length(col) =",length(colbar$col),
                             "length(breaks) =",length(colbar$breaks)))

    if (length(colbar$col) != length(colbar$breaks)-1)
      stop('colbar.ini: Error in setting colbar!')
    ##}
    invisible(colbar)
}

colscal <- function(n=14,col="t2m",rev=TRUE,alpha=NULL,
                    test=FALSE,verbose=FALSE) {

  test.col <- function(r,g,b) {
    dev.new()
    par(bty="n")
    plot(r,col="red")
    points(b,col="blue")
    points(g,col="green")
  }

  if (verbose) print(paste('colscal:',col))
  if (is.null(col)) col <- 't2m'
  if (is.null(alpha)) alpha <- 1
  # Set up colour-palette
  col <- tolower(col)
  x <- 1:n
  r0 <- round(n*0.55)
  g0 <- round(n*0.5)
  b0 <- round(n*0.45)
  s <- -0.1/n
  if (n < 30) sg <- s*2.5 else sg <- s
  n1 <- g0; n2 <- n-n1
  
#R	G	B
  seNorgeT <- c(204,  0,    0,	
               255, 25,    0,	
               255, 102,   0,	
               255, 179,   0,	
               255, 230,  77,	
               255, 255,  64,	
               255, 255, 190,	
               217, 255, 255,	
               179, 255, 255,	
               128, 235, 255,	
               64, 204, 255,	
               0, 153, 255,	
               0,  25, 255,	
               0,   0, 153)	
  dim(seNorgeT) <- c(3,14)

  seNorgeP <- c(0, 0, 153,
                0, 25, 255,
                0, 153, 255,
                64, 204, 255,
                128, 235, 255,
                179, 255, 255,
                217, 255, 255,
                229, 229, 229)
  dim(seNorgeP) <- c(3,8)

  if (!is.null(alpha)) alpha <- rep(alpha[1],n)
  
  if ( (col[1]=="bwr") | (col[1]=="slp") | (col[1]=="mslp") |
      (col[1]=="pressure") ) {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
    g <- exp(sg*(x - g0)^2)^2
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
    col <- rgb(r,g,b,alpha)
  } else if (col[1]=="rwb") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
    g <- exp(sg*(x - g0)^2)^2
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
    if (is.null(alpha)) col <- rgb(b,g,r)  else
                        col <- rgb(r,g,b,alpha)
  } else if (col[1]=="faint.bwr") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
    g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
    col <- rgb(r,g,b,alpha)
  } else if (col[1]=="faint.rwb") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
    g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
    col <- rgb(r,g,b,alpha)
  } else if ( (col[1]=="precip") | (col[1]=="mu") | (col[1]=="fw") |
              (col[1]=="f[w]") | (col[1]=="tp") | (col[1]=='rr') | (col[1]=='prate') ) {
    r <- approx(seNorgeP[1,],n=n)$y/255
    g <- approx(seNorgeP[2,],n=n)$y/255
    b <- approx(seNorgeP[3,],n=n)$y/255
    col <- rgb(r,g,b,alpha)
    rev <- TRUE
  } else if (col[1]=="rainbow") {
    col <- rainbow(n,start=0,end=4/6,alpha=alpha)
  } else if (col[1]=="gray.colors") {
    col <- gray.colors(n,alpha=alpha)
  } else if (col[1]=="heat.colors") {
    col <- heat.colors(n,alpha=alpha)
  } else if (col[1]=="terrain.colors") {
    col <- terrain.colors(n,alpha=alpha)
  } else if (col[1]=="topo.colors") {
    col <- topo.colors(n,alpha=alpha)
  } else if (col[1]=="cm.colors") {
    col <- cm.colors(n,alpha=alpha)
  } else if (col[1]==tolower("grmg")) {
    cols <- list(
     r=c(0,0,0,0,0.316,0.526,0.737,1,1,1,1,1,0.947,0.737,0.526,0.316),
     g=c(0.316,0.526,0.737,0.947,1,1,1,1,0.947,0.737,0.526,0.316,0,0,0,0),
     b=c(0,0,0,0,0.316,0.526,0.737,1,1,1,1,1,0.947,0.737,0.526,0.316))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]==tolower("brbu")) {
    cols <- list(
     r=c(0.2,0.4,0.6,0.8,0.85,0.95,0.8,0.6,0.4,0.2,0,0),
     g=c(0.1,0.187,0.379,0.608,0.688,0.855,0.993,0.973,0.94,0.893,0.667,0.48),
     b=c(0,0,0.21,0.480,0.595,0.808,1,1,1,1,0.8,0.6))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]==tolower("budor")) {
    cols <- list(
     r=c(0.12,0.32,0.6,0.7,0.8,0.9,1,1,1,1,0.8,0.6),
     g=c(0.56,0.768,0.98,0.99,0.997,1,0.9,0.793,0.68,0.56,0.347,0.250),
     b=c(0.6,0.8,1,1,1,1,0.8,0.6,0.4,0.2,0,0))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]==tolower("budrd")) {
    cols <- list(
     r=c(0.142,0.097,0.16,0.24,0.34,0.46,0.6,0.74,0.92,1,
       1,1,1,1,1,0.97,0.85,0.65),
     g=c(0,0.112,0.342,0.531,0.692,0.829,0.920,0.978,1,1,
       0.948,0.840,0.676,0.472,0.240,0.155,0.085,0),
     b=c(0.85,0.97,1,1,1,1,1,1,1,0.92,0.74,0.6,0.46,0.34,0.24,0.21,0.187,0.13))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]==tolower("bugr")) {
    cols <- list(
     r=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.9,0.8,0.7,0.6,0.4,0.2,0),
     g=c(0,0.2,0.4,0.6,0.7,0.8,0.9,1,1,1,1,1,1,1),
     b=c(1,1,1,1,1,1,1,0.9,0.8,0.7,0.6,0.4,0.2,0))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]==tolower("bugy")) {
    cols <- list(
     r=c(0,0.4,0.6,0.8,0.9,0.6,0.4,0.2),
     g=c(0.6,0.9,1,1,0.9,0.6,0.4,0.2),
     b=c(0.8,1,1,1,0.9,0.6,0.4,0.2))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]==tolower("buor")) {
    cols <- list(
     r=c(0,0.1,0.2,0.4,0.6,0.8,1,1,1,1,1,1),
     g=c(0.167,0.4,0.6,0.8,0.933,1,1,0.933,0.8,0.6,0.4,0.167),
     b=c(1,1,1,1,1,1,0.8,0.6,0.4,0.2,0.1,0))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]==tolower("buorr")) {
    cols <- list(
     r=c(0.03,0.2,0.35,0.55,0.75,0.9,0.97,1,1,1,1,1,1,1),
     g=c(0.353,0.467,0.567,0.7,0.833,0.933,0.98,1,1,1,0.8,0.6,0.4,0),
     b=c(1,1,1,1,1,1,1,0.8,0.6,0,0,0,0,0))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]=="bu") {
    cols <- list(
     r=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0),
     g=c(1,0.983,0.95,0.9,0.833,0.75,0.65,0.533,0.4,0.250),
     b=c(1,1,1,1,1,1,1,1,1,1))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]=="cat") {
    cols <- list(
     r=c(1,1,1,1,0.7,0.2,0.65,0.1,0.8,0.4,1,0.9),
     g=c(0.75,0.5,1,1,1,1,0.93,0.7,0.75,0.3,0.6,0.10),
     b=c(0.5,0,0.6,0.2,0.55,0,1,1,1,1,0.75,0.2))
     cols <- lapply(cols,function(x) approx(x,n=n)$y)
     col <- rgb(cols$r,cols$g,cols$b,alpha)
  } else if (col[1]=="warm") {
    r <- approx(seNorgeT[1,1:7],n=n)$y/255
    g <- approx(seNorgeT[2,1:7],n=n)$y/255
    b <- approx(seNorgeT[3,1:7],n=n)$y/255
    col <- rgb(r,g,b,alpha)    
  }  else if (col[1]=="cold") {
    r <- approx(seNorgeT[1,8:14],n=n)$y/255
    g <- approx(seNorgeT[2,8:14],n=n)$y/255
    b <- approx(seNorgeT[3,8:14],n=n)$y/255
    col <- rgb(r,g,b,alpha)    
  } else {
    r <- approx(seNorgeT[1,],n=n)$y/255
    g <- approx(seNorgeT[2,],n=n)$y/255
    b <- approx(seNorgeT[3,],n=n)$y/255
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

colbar <- function(breaks,col,fig=c(0.15,0.2,0.15,0.3),horiz=FALSE) {
  if (horiz) {
    par(xaxt="s",yaxt="n",fig=fig,mar=c(1,0,0,0),new=TRUE,las=1,cex.axis=0.6)
    image(breaks,1:2,cbind(breaks,breaks),col=col)
  } else {
    par(xaxt="n",yaxt="s",fig=fig,mar=c(0,1,0,0),new=TRUE,las=1,cex.axis=0.6)
    image(1:2,breaks,rbind(breaks,breaks),col=col)
  }
}

col.bar <- function(breaks,horiz=TRUE,pch=21,v=1,h=1,col=col,cex=2,cex.lab=0.6,
                    type="r",verbose=FALSE,vl=0.5,border=FALSE,...) {
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
            rect(xleft= k  + xleft + steps[i] ,xright= k + xleft + steps[i+1],
                 ybottom=ybottom,ytop=ytop,col=col[i],border=border)
            
            ## text(x = k + xleft + steps[i], y = ybottom - 1,labels=sprintf("%.1f",icn[i]),cex=cex)
        }
        else if (type == "p") { ## "p" points
            points(x= k + xleft + (steps[i]+ steps[i+1])/2, y=(ybottom + ytop)/2,
                   pch=pch, bg=col[i],cex=cex,...)
            
        }
        
        text(x = k + xleft + (steps[i]+ steps[i+1])/2,  y = ybottom - vl,
             labels=levels(cut(breaks,breaks))[i],col="grey50",cex=cex.lab)
    } 
    #par(fig=par0$fig)
}

colbar2 <- function(x,col) {
    par(mar=c(1,0,0,0),fig=c(0.5,1,0.665,0.695),new=TRUE,cex.axis=0.6)
    nl <- pretty(x)
    n <- length(nl)
    image(cbind(1:n,1:n),col=col) 
    par(xaxt="s",new=new)
    axis(1,at=seq(0,1,length=length(nl)),label=nl)
}

## colorscheme <- function(pal="BuDOr",n=12,alpha=1) {
##   if (pal=="GrMg") {
##     cols <- list(
##      r=c(0,0,0,0,0.316,0.526,0.737,1,1,1,1,1,0.947,0.737,0.526,0.316),
##      g=c(0.316,0.526,0.737,0.947,1,1,1,1,0.947,0.737,0.526,0.316,0,0,0,0),
##      b=c(0,0,0,0,0.316,0.526,0.737,1,1,1,1,1,0.947,0.737,0.526,0.316))
##   } else if (pal=="BrBu") {
##     cols <- list(
##      r=c(0.2,0.4,0.6,0.8,0.85,0.95,0.8,0.6,0.4,0.2,0,0),
##      g=c(0.1,0.187,0.379,0.608,0.688,0.855,0.993,0.973,0.94,0.893,0.667,0.48),
##      b=c(0,0,0.21,0.480,0.595,0.808,1,1,1,1,0.8,0.6))
##   } else if (pal=="BuDOr") {
##     cols <- list(
##      r=c(0.12,0.32,0.6,0.7,0.8,0.9,1,1,1,1,0.8,0.6),
##      g=c(0.56,0.768,0.98,0.99,0.997,1,0.9,0.793,0.68,0.56,0.347,0.250),
##      b=c(0.6,0.8,1,1,1,1,0.8,0.6,0.4,0.2,0,0))
##   } else if (pal=="BuDRd") {
##     cols <- list(
##      r=c(0.142,0.097,0.16,0.24,0.34,0.46,0.6,0.74,0.92,1,
##        1,1,1,1,1,0.97,0.85,0.65),
##      g=c(0,0.112,0.342,0.531,0.692,0.829,0.920,0.978,1,1,
##        0.948,0.840,0.676,0.472,0.240,0.155,0.085,0),
##      b=c(0.85,0.97,1,1,1,1,1,1,1,0.92,0.74,0.6,0.46,0.34,0.24,0.21,0.187,0.13))
##   } else if (pal=="BuGr") {
##     cols <- list(
##      r=c(0,0.2,0.4,0.6,0.7,0.8,0.9,0.9,0.8,0.7,0.6,0.4,0.2,0),
##      g=c(0,0.2,0.4,0.6,0.7,0.8,0.9,1,1,1,1,1,1,1),
##      b=c(1,1,1,1,1,1,1,0.9,0.8,0.7,0.6,0.4,0.2,0))
##   } else if (pal=="BuGy") {
##     cols <- list(
##      r=c(0,0.4,0.6,0.8,0.9,0.6,0.4,0.2),
##      g=c(0.6,0.9,1,1,0.9,0.6,0.4,0.2),
##      b=c(0.8,1,1,1,0.9,0.6,0.4,0.2))
##   } else if (pal=="BuOr") {
##     cols <- list(
##      r=c(0,0.1,0.2,0.4,0.6,0.8,1,1,1,1,1,1),
##      g=c(0.167,0.4,0.6,0.8,0.933,1,1,0.933,0.8,0.6,0.4,0.167),
##      b=c(1,1,1,1,1,1,0.8,0.6,0.4,0.2,0.1,0))
##   } else if (pal=="BuOrR") {
##     cols <- list(
##      r=c(0.03,0.2,0.35,0.55,0.75,0.9,0.97,1,1,1,1,1,1,1,),
##      g=c(0.353,0.467,0.567,0.7,0.833,0.933,0.98,1,1,1,0.8,0.6,0.4,0),
##      b=c(1,1,1,1,1,1,1,0.8,0.6,0,0,0,0,0))
##   } else if (pal=="Bu") {
##     cols <- list(
##      r=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0),
##      g=c(1,0.983,0.95,0.9,0.833,0.75,0.65,0.533,0.4,0.250),
##      b=c(1,1,1,1,1,1,1,1,1,1))
##   } else if (pal=="Cat") {
##     cols <- list(
##      r=c(1,1,1,1,0.7,0.2,0.65,0.1,0.8,0.4,1,0.9),
##      g=c(0.75,0.5,1,1,1,1,0.93,0.7,0.75,0.3,0.6,0.10),
##      b=c(0.5,0,0.6,0.2,0.55,0,1,1,1,1,0.75,0.2))
##   }
##   cols <- lapply(cols,function(x) approx(x,n=n)$y)
##   return(rgb(cols$r,cols$g,cols$b,alpha))
## }

## from DSE ...
## col.bar <- function(x,horiz=TRUE,v=1,h=1,col=col,cex=0.7,type="r",...) {
    
##     xleft <- par()$usr[1] 
##     xright <- par()$usr[2]
##     ybottom <- par()$usr[4] - 1 - h
##     ytop <-  par()$usr[4] - 1 
    
##     steps <-   seq(0, (xright -xleft - v * (length(col))) , (xright - xleft - v * (length(col)))/(length(col))) # 
##     nsteps <- length(steps) - 1 
##     icn <- seq(0,1,1/nsteps) ; print(icn)
##     k <- 0
##     for (i in 1 :nsteps) {  
##         if (!is.null(v)) 
##             if (i == 1) k <- v/2 else k <- k + v  
##         if (type == "r") { ## "r" for rectangle
##             rect(xleft= k  + xleft + steps[i] ,xright= k + xleft + steps[i+1],ybottom=ybottom,ytop=ytop,col=col[i])
##             text(x = k + xleft + steps[i],  y = ybottom - 1, labels=sprintf("%.1f",icn[i]),cex=cex)
##             text(x = k + xleft + steps[i+1],y = ybottom - 1, labels=sprintf("%.1f",icn[i+1]),cex=cex)
##             ## text(x = k + xleft + steps[i], y = ybottom - 1,labels=sprintf("%.1f",icn[i]),cex=cex)
##         }
##         else if (type == "p") { ## "p" points
##             points(x= k + xleft + (steps[i]+ steps[i+1])/2, y=(ybottom + ytop)/2,pch=21, bg=col[i],cex=v)
##             text(x = k + xleft + steps[i],  y = ybottom - 1, labels=sprintf("%.1f",icn[i]),cex=cex)
##             text(x = k + xleft + steps[i+1],y = ybottom - 1, labels=sprintf("%.1f",icn[i+1]),cex=cex)
##         }
##     }
## }


## copied from infographics.R

#seNorge
#nodata	10000
#nobits	16
#header1	*Temperatur
#header2	*Siste døgn (18-18 UTC)
#legend	*Grader Celsius
#From	To	R	G	B	Forklaring
#2931	10000	204	0	0	*Over 20 
#2881	2931	255	25	0	*15 - 20
#2831	2881	255	102	0	*10 - 15
#2781	2831	255	179	0	*5 - 10
#2761	2781	255	230	77	*3 - 5
#2741	2761	255	255	64	*1 - 3
#2731	2741	255	255	190	*0 - 1
#2721	2731	217	255	255	*÷1 - 0
#2701	2721	179	255	255	*÷3 - ÷1
#2681	2701	128	235	255	*÷5 - ÷3
#2631	2681	64	204	255	*÷10 - ÷5
#2581	2631	0 	153	255	*÷15 - ÷10
#2531	2581	0 	25	255	*÷20 - ÷15
#0	2531	0	0	153	*Under ÷20


#nodata	10000
#nobits	16
#header1	*Nedbør
#header2	*Siste døgn (06-06 UTC)
#legend	*mm
#From	To	R	G	B	Forklaring
#1500	10000	0	0	153	*Over 150 
#750	1500	0 	25	255	*75 - 150 
#500	750	0 	153	255	*50 - 75
#300	500	64	204	255	*30 - 50
#200	300	128	235	255	*20 - 30
#100	200	179	255	255	*10 - 20
#1	100	217	255	255	*Under 10
#0	1	229	229	229	*Ikke nedbør

## colscal <- function(n=14,col="t2m",alpha=NULL,test=FALSE) {

##   test.col <- function(r,g,b) {
##     dev.new()
##     par(bty="n")
##     plot(r,col="red")
##     points(b,col="blue")
##     points(g,col="green")
##   }

##   # Set up colour-palette
##   col <- tolower(col)
##   x <- 1:n
##   r0 <- round(n*0.55)
##   g0 <- round(n*0.5)
##   b0 <- round(n*0.45)
##   s <- -0.1/n
##   if (n < 30) sg <- s*2.5 else sg <- s
##   n1 <- g0; n2 <- n-n1
  
## #R	G	B
##   seNorgeT <- c(204,  0,    0,	
##                255, 25,    0,	
##                255, 102,   0,	
##                255, 179,   0,	
##                255, 230,  77,	
##                255, 255,  64,	
##                255, 255, 190,	
##                217, 255, 255,	
##                179, 255, 255,	
##                128, 235, 255,	
##                64, 204, 255,	
##                0, 153, 255,	
##                0,  25, 255,	
##                0,   0, 153)	
##   dim(seNorgeT) <- c(3,14)

##   seNorgeP <- c(0, 0, 153,
##                 0, 25, 255,
##                 0, 153, 255,
##                 64, 204, 255,
##                 128, 235, 255,
##                 179, 255, 255,
##                 217, 255, 255,
##                 229, 229, 229)
##   dim(seNorgeP) <- c(3,8)

##   ##if (!is.null(col))
##   ##  if ((length(col)==1) & is.character(col) &
##   ##      (sum(is.element(c('t2m','precip','bwr','rwb','mu','fw','tp',
##   ##                        'faint.bwr','faint.rwb','rainbow',
##   ##                        'gray.colors','heat.colors','terrain.colors',
##   ##                        'topo.colors','cm.colors'),col))==0))
##   ##      col <- 'bwr'

##   #if (exists("r")) remove(r)
##   #if (exists("g")) remove(g) 
##   #if (exists("b")) remove(b)

##   if (!is.null(alpha)) alpha <- rep(alpha[1],n)
  
##   if ( (col[1]=="bwr") | (col[1]=="slp") | (col[1]=="mslp") |
##       (col[1]=="pressure") ) {
##     r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
##     g <- exp(sg*(x - g0)^2)^2
##     b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
##     if (is.null(alpha)) col <- rgb(r,g,b) else
##                         col <- rgb(r,g,b,alpha)
##   } else if (col[1]=="rwb") {
##     r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
##     g <- exp(sg*(x - g0)^2)^2
##     b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
##     if (is.null(alpha)) col <- rgb(b,g,r)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if (col[1]=="faint.bwr") {
##     r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
##     g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
##     b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
##     if (is.null(alpha)) col <- rgb(r,g,b)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if (col[1]=="faint.rwb") {
##     r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
##     g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
##     b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
##     if (is.null(alpha)) col <- rgb(b,g,r)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if ( (col[1]=="t2m") | (col[1]=="tmax") | (col[1]=="tmin") |
##              (col[1]=="sst")  | (col[1]=="air") ){
##     r <- approx(seNorgeT[1,],n=n)$y/255
##     g <- approx(seNorgeT[2,],n=n)$y/255
##     b <- approx(seNorgeT[3,],n=n)$y/255
##     if (is.null(alpha)) col <- rgb(b,g,r)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if ((col[1]=="precip") | (col[1]=="mu") | (col[1]=="fw") |
##              (col[1]=="f[w]") | (col[1]=="tp")) {
##     r <- approx(seNorgeP[1,],n=n)$y/255
##     g <- approx(seNorgeP[2,],n=n)$y/255
##     b <- approx(seNorgeP[3,],n=n)$y/255
##     if (is.null(alpha)) col <- rgb(r,g,b)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if (col[1]=="rainbow") {
##     col <- rainbow(n,start=0,end=4/6)
##   } else if (col[1]=="gray.colors") {
##     col <- gray.colors(n)
##   } else if (col[1]=="heat.colors") {
##     col <- heat.colors(n)
##   } else if (col[1]=="terrain.colors") {
##     col <- terrain.colors(n)
##   } else if (col[1]=="topo.colors") {
##     col <- topo.colors(n)
##   } else if (col[1]=="cm.colors") {
##     col <- cm.colors(n)
##   }

##   if (test) { #& !exists("r")) {
##     RGB <- col2rgb(col)/255
##     r <- RGB[1,]; g <- RGB[2,]; b <- RGB[3,]
##   }
  
##   if (test) test.col(r,g,b)
##   return(col)
## }

## colbar <- function(scale,col,fig=c(0.15,0.2,0.15,0.3)) {
##   par(xaxt="n",yaxt="s",fig=fig,mar=c(0,1,0,0),new=TRUE,las=1,cex.axis=0.6)
##   image(1:2,scale,rbind(scale,scale),col=col
##     )

## from vis.R
#seNorge
#nodata	10000
#nobits	16
#header1	*Temperatur
#header2	*Siste døgn (18-18 UTC)
#legend	*Grader Celsius
#From	To	R	G	B	Forklaring
#2931	10000	204	0	0	*Over 20 
#2881	2931	255	25	0	*15 - 20
#2831	2881	255	102	0	*10 - 15
#2781	2831	255	179	0	*5 - 10
#2761	2781	255	230	77	*3 - 5
#2741	2761	255	255	64	*1 - 3
#2731	2741	255	255	190	*0 - 1
#2721	2731	217	255	255	*÷1 - 0
#2701	2721	179	255	255	*÷3 - ÷1
#2681	2701	128	235	255	*÷5 - ÷3
#2631	2681	64	204	255	*÷10 - ÷5
#2581	2631	0 	153	255	*÷15 - ÷10
#2531	2581	0 	25	255	*÷20 - ÷15
#0	2531	0	0	153	*Under ÷20


#nodata	10000
#nobits	16
#header1	*Nedbør
#header2	*Siste døgn (06-06 UTC)
#legend	*mm
#From	To	R	G	B	Forklaring
#1500	10000	0	0	153	*Over 150 
#750	1500	0 	25	255	*75 - 150 
#500	750	0 	153	255	*50 - 75
#300	500	64	204	255	*30 - 50
#200	300	128	235	255	*20 - 30
#100	200	179	255	255	*10 - 20
#1	100	217	255	255	*Under 10
#0	1	229	229	229	*Ikke nedbør

## colscal <- function(n=14,col="t2m",alpha=NULL,test=FALSE) {

##   test.col <- function(r,g,b) {
##     dev.new()
##     par(bty="n")
##     plot(r,col="red")
##     points(b,col="blue")
##     points(g,col="green")
##   }

##   # Set up colour-palette
##   col <- tolower(col)
##   x <- 1:n
##   r0 <- round(n*0.55)
##   g0 <- round(n*0.5)
##   b0 <- round(n*0.45)
##   s <- -0.1/n
##   if (n < 30) sg <- s*2.5 else sg <- s
##   n1 <- g0; n2 <- n-n1
  
## #R	G	B
##   seNorgeT <- c(204,  0,    0,	
##                255, 25,    0,	
##                255, 102,   0,	
##                255, 179,   0,	
##                255, 230,  77,	
##                255, 255,  64,	
##                255, 255, 190,	
##                217, 255, 255,	
##                179, 255, 255,	
##                128, 235, 255,	
##                64, 204, 255,	
##                0, 153, 255,	
##                0,  25, 255,	
##                0,   0, 153)	
##   dim(seNorgeT) <- c(3,14)

##   seNorgeP <- c(0, 0, 153,
##                 0, 25, 255,
##                 0, 153, 255,
##                 64, 204, 255,
##                 128, 235, 255,
##                 179, 255, 255,
##                 217, 255, 255,
##                 229, 229, 229)
##   dim(seNorgeP) <- c(3,8)

##   ##if (!is.null(col))
##   ##  if ((length(col)==1) & is.character(col) &
##   ##      (sum(is.element(c('t2m','precip','bwr','rwb','mu','fw','tp',
##   ##                        'faint.bwr','faint.rwb','rainbow',
##   ##                        'gray.colors','heat.colors','terrain.colors',
##   ##                        'topo.colors','cm.colors'),col))==0))
##   ##      col <- 'bwr'

##   #if (exists("r")) remove(r)
##   #if (exists("g")) remove(g) 
##   #if (exists("b")) remove(b)

##   if (!is.null(alpha)) alpha <- rep(alpha[1],n)
  
##   if ( (col[1]=="bwr") | (col[1]=="slp") | (col[1]=="mslp") |
##       (col[1]=="pressure") ) {
##     r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
##     g <- exp(sg*(x - g0)^2)^2
##     b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
##     if (is.null(alpha)) col <- rgb(r,g,b) else
##                         col <- rgb(r,g,b,alpha)
##   } else if (col[1]=="rwb") {
##     r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
##     g <- exp(sg*(x - g0)^2)^2
##     b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
##     if (is.null(alpha)) col <- rgb(b,g,r)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if (col[1]=="faint.bwr") {
##     r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
##     g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
##     b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
##     if (is.null(alpha)) col <- rgb(r,g,b)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if (col[1]=="faint.rwb") {
##     r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
##     g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
##     b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
##     if (is.null(alpha)) col <- rgb(b,g,r)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if ( (col[1]=="t2m") | (col[1]=="tmax") | (col[1]=="tmin") |
##              (col[1]=="sst")  | (col[1]=="air") ){
##     r <- approx(seNorgeT[1,],n=n)$y/255
##     g <- approx(seNorgeT[2,],n=n)$y/255
##     b <- approx(seNorgeT[3,],n=n)$y/255
##     if (is.null(alpha)) col <- rgb(b,g,r)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if ((col[1]=="precip") | (col[1]=="mu") | (col[1]=="fw") |
##              (col[1]=="f[w]") | (col[1]=="tp")) {
##     r <- approx(seNorgeP[1,],n=n)$y/255
##     g <- approx(seNorgeP[2,],n=n)$y/255
##     b <- approx(seNorgeP[3,],n=n)$y/255
##     if (is.null(alpha)) col <- rgb(r,g,b)  else
##                         col <- rgb(r,g,b,alpha)
##   } else if (col[1]=="rainbow") {
##     col <- rainbow(n,start=0,end=4/6)
##   } else if (col[1]=="gray.colors") {
##     col <- gray.colors(n)
##   } else if (col[1]=="heat.colors") {
##     col <- heat.colors(n)
##   } else if (col[1]=="terrain.colors") {
##     col <- terrain.colors(n)
##   } else if (col[1]=="topo.colors") {
##     col <- topo.colors(n)
##   } else if (col[1]=="cm.colors") {
##     col <- cm.colors(n)
##   }

##   if (test) { #& !exists("r")) {
##     RGB <- col2rgb(col)/255
##     r <- RGB[1,]; g <- RGB[2,]; b <- RGB[3,]
##   }
  
##   if (test) test.col(r,g,b)
##   return(col)
## }

## colbar <- function(scale,col,fig=c(0.15,0.2,0.15,0.3)) {
##   par(xaxt="n",yaxt="s",fig=fig,mar=c(0,1,0,0),new=TRUE,las=1,cex.axis=0.6)
##   image(1:2,scale,rbind(scale,scale),col=col)
## }
