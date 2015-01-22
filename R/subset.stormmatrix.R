
## Author K. Parding
## Based on subset.station.R by Rasmus E. Benestad and A. Mezghani
## Last updated 22.01.2015

subset.stormmatrix <- function(x,it=NULL,is=NULL,verbose=FALSE) {

    x0 <- x
    cls <- class(x)
    if (is.null(it) & is.null(is)) return(x)
 
    l <- dim(x)[1]
    if (is.null(is)) is <- rep(TRUE,l)
    if (is.null(it)) it <- rep(TRUE,l)
    
    ## get time in t
    t <- strptime(x[,colnames(x)=="start"],format="%Y%m%d%H")
    yr <- as.numeric(strftime(t,"%Y"))
    mo <- as.numeric(strftime(t,"%m"))
    dy <- as.numeric(strftime(t,"%d"))
    
    # Generate sequence of days, months or years if range of it value is given
    if (!is.null(it)) {
      is.months <- function(x) (sum(is.element(tolower(substr(x,1,3)),
                                               tolower(month.abb)))>0)
      is.seasons <- function(x) (sum(is.element(tolower(substr(x,1,3)),
                                                names(season.abb())))>0)
      is.dates <- function(x) (!is.months(x) & levels(factor(nchar(x)))==10)
      is.years <- function(x) (!is.months(x) & levels(factor(nchar(x)))==4)

      if (is.months(it)) {
        if (verbose) print('Monthly selected')
        ii <- is.element(mo,(1:12)[is.element(tolower(month.abb),
                                 tolower(substr(it,1,3)))])
      } else if (is.seasons(it)) {
        if (verbose) print("Seasonally selected")
        if (verbose) print(table(mo))
        if (verbose) print(eval(parse(text=paste('season.abb()$',it,sep=''))))
        ii <- is.element(mo,eval(parse(text=paste('season.abb()$',it,sep=''))))
      } else if (is.dates(it)) {
        it <- as.Date(it)
        if ( length(it) == 2 ) {
          if (verbose) print('Between two dates')
          if (verbose) print(it)          
          id <- strftime(seq(it[1],it[2],by='day'),format="%Y%m%d")
          d <- strftime(t,format="%Y%m%d")
          ii <- is.element(d,id)
        } else { 
          if (verbose) print('it is a string of dates')
          ii <- is.element(t,it)
        }
     } else if (is.years(it)) {
        it <- as.integer(it)
        if ( length(it) == 2 ) {
          if (verbose) print('Between two years')
          if (verbose) print(it)          
          ii <- is.element(yr,seq(it[1],it[2],1))
        } else { 
          if (verbose) print('it is a string of years')
          ii <- is.element(yr,it)
        }
      } else {
        ii <- rep(FALSE,length(t))
        warning("subset.station: did not recognise the selection citerion for 'it'")
      }
    } else ii <- rep(TRUE,l)
        
    # is can be a list to select region or according to other criterion
    if (inherits(is,'list')) {
        selx <- rep(TRUE,l); sely <- selx;
        selp <- selx; selF <- selx
        nms <- names(is)
        ix <- grep('lon',tolower(nms))
        iy <- grep('lat',tolower(nms))
        ip <- grep('slp',tolower(nms))
        iF <- grep('FUN',nms)
        if (length(ix)>0) slon <- is[[ix]] else slon <- NULL
        if (length(iy)>0) slat <- is[[iy]] else slat <- NULL
        if (length(ip)>0) sslp <- is[[ip]] else sslp <- NULL        
        if (length(iF)>0) sFUN <- is[[iF]] else sFUN <- NULL
        if (length(slon)==2) {
          fn <- function(x) any(x>=min(slon) & x<=max(slon))
          selx <- apply(x[,colnames(x)=='lon'],1,fn)
        }
        if (length(slat)==2) {
          fn <- function(x) any(x>=min(slat) & x<=max(slat))
          sely <- apply(x[,colnames(x)=='lat'],1,fn)
        }
        if (length(sslp)>0) {
          fn <- function(x) any(x>=min(sslp) & x<=max(sslp))
          selp <- apply(x$slp,1,fn)
          }
        if (length(sFUN)>0) selF <- apply(x,1,sFUN) # Not quite finished...
        ij <- selx & sely & selp & selF
     } else ij <- rep(TRUE,l)

    ist <- (1:l)[(ii & ij)]
    y <- x[ist,]
    
    class(y) <- cls
    y <- attrcp(x,y)
    attr(y,'history') <- history.stamp(x)   
    invisible(y)
}
    
