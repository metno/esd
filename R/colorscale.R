
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
    
    if (is.logical(colbar)) colbar <- NULL
    
    ##if (!is.null(colbar)) {
    if (verbose) print('sort out the colours')
    if (is.null(colbar$pal)) colbar$pal <- varid(x)[1]
    if (is.null(colbar$rev)) colbar$rev <- FALSE
    if (is.null(colbar$n)) colbar$n <- 10
    ##if (is.null(colbar$pal)) {
        if (is.null(FUN) | !is.precip(x)) colbar$pal <- 't2m' else
        if ( (is.precip(x)) & ( (FUN=='sum') | (FUN=='trend') |
                               (FUN=='wetmean') | (FUN=='mean')) ) {
            colbar$pal <- 'precip'
            colbar$rev <- TRUE
        } else colbar$pal <- 't2m'
    ##} 
    if (is.zoo(x)) x <- coredata(x)
    if (is.null(colbar$breaks)) {
        colbar$breaks <- pretty(x,n=colbar$n)
    } else if (length(colbar$breaks)==2)
        colbar$breaks <- seq(colbar$breaks[1],colbar$breaks[2],
                             length=colbar$n)
    colbar$n <- length(colbar$breaks)-1
    if (is.null(colbar$type)) colbar$type <- 'p'
    if (is.null(colbar$cex)) colbar$cex <- 2
    if (is.null(colbar$h)) colbar$h <- 0.6
    if (is.null(colbar$v)) colbar$v <- 1
    if (is.null(colbar$pos)) colbar$pos <- 0.05
    if (is.null(colbar$show)) colbar$show <-TRUE 
    if (verbose) print(colbar)
    colbar$col <- colscal(n=colbar$n,col=colbar$pal,
                          rev=colbar$rev,verbose=verbose)
    if (verbose) print(paste("length(col) =",length(colbar$col)))
    col <- colscal(n=colbar$n,col=colbar$pal,rev=colbar$rev)       
    
    if (!is.null(FUN)) {
        if (is.null(colbar$breaks) & !inherits(x,"stationmeta")) {
            colbar$breaks <- pretty(x,n=colbar$n)
        } else if (length(colbar$breaks)==2)
            colbar$breaks <- seq(colbar$breaks[1],colbar$breaks[2],
                                 length=colbar$n)
        colbar$n <- length(colbar$breaks)-1
    }

    if (!inherits(x,"stationmeta"))
        colbar$col <- colscal(n=colbar$n,col=colbar$pal,rev=colbar$rev,verbose=verbose)
    if (verbose) print(paste("length(col) =",length(colbar$col)))

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

  ##if (!is.null(col))
  ##  if ((length(col)==1) & is.character(col) &
  ##      (sum(is.element(c('t2m','precip','bwr','rwb','mu','fw','tp',
  ##                        'faint.bwr','faint.rwb','rainbow',
  ##                        'gray.colors','heat.colors','terrain.colors',
  ##                        'topo.colors','cm.colors'),col))==0))
  ##      col <- 'bwr'

  #if (exists("r")) remove(r)
  #if (exists("g")) remove(g) 
  #if (exists("b")) remove(b)

  if (!is.null(alpha)) alpha <- rep(alpha[1],n)
  
  if ( (col[1]=="bwr") | (col[1]=="slp") | (col[1]=="mslp") |
      (col[1]=="pressure") ) {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0,1,length=n1),rep(1,n2))
    g <- exp(sg*(x - g0)^2)^2
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0,length=n1))
    if (is.null(alpha)) col <- rgb(r,g,b) else
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
    if (is.null(alpha)) col <- rgb(r,g,b)  else
                        col <- rgb(r,g,b,alpha)
  } else if (col[1]=="faint.rwb") {
    r <- exp(s*(x - r0)^2)^0.5 * c(seq(0.5,1,length=n1),rep(1,n2))
    g <- min(exp(sg*(x - g0)^2)^2 + 0.5,1)
    b <- exp(s*(x - b0)^2)^0.5 * c(rep(1,n2),seq(1,0.5,length=n1))
    if (is.null(alpha)) col <- rgb(b,g,r)  else
                        col <- rgb(r,g,b,alpha)
  } else if ( (col[1]=="precip") | (col[1]=="mu") | (col[1]=="fw") |
              (col[1]=="f[w]") | (col[1]=="tp") | (col[1]=='rr') | (col[1]=='prate') ) {
    r <- approx(seNorgeP[1,],n=n)$y/255
    g <- approx(seNorgeP[2,],n=n)$y/255
    b <- approx(seNorgeP[3,],n=n)$y/255
    if (is.null(alpha)) col <- rgb(r,g,b)  else
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
  } else {
    r <- approx(seNorgeT[1,],n=n)$y/255
    g <- approx(seNorgeT[2,],n=n)$y/255
    b <- approx(seNorgeT[3,],n=n)$y/255
    if (is.null(alpha)) col <- rgb(b,g,r)  else
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
