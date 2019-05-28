#' Visualise trends for multiple overlapping periods
#' 
#' Produce a plot showing trends for multiple periods within a time series. The
#' strength of the trend is represented by the color scale and significant
#' trends are marked with black borders.
#' 
#' 
#' @param x the 'x' argument provides the time series for which the trend
#' analysis is performed. Only zoo objects are accepted.
#' @param minlen minimum time interval to calculate trends for in units of
#' years.
#' @param unitlabel unit of x.
#' @param varlabel name of x.
#' @param vmax upper limit of trend scale.
#' @param show.significance TRUE to mark statistically significant trends.
#' @param pmax maximum p-value of trends marked as significant.
#' @param verbose TRUE or FALSE.
#' @author Kajsa Parding
#' @keywords trend
#' @examples
#' 
#' 
#' t <- seq(as.Date("1955-01-01"),as.Date("2004-12-31"),by=1)
#' x <- zoo(sample(seq(-30,30,1e-1),length(t),rep=TRUE),order.by=t)
#' vis.trends(x,show.significance=FALSE)
#' 
#' data(Oslo)
#' vis.trends(Oslo, unitlabel="oC", varlabel = "Temperature",
#'   pmax = 1e-2, minlen = 40)
#' vis.trends(subset(Oslo,it='jja'), unitlabel="oC",
#'   varlabel = "Temperature JJA",
#'   pmax = 1e-3, vmax=0.5, minlen = 40)
#' vis.trends(subset(Oslo,it='mam'), unitlabel="oC",
#'   varlabel = "Temperature MAM",
#'   pmax = 1e-3, vmax=0.5, minlen = 40)
#' 
#' @export vis.trends
vis.trends <- function(x,...,unitlabel="unit",varlabel="",is=1,
                       pmax=0.01,minlen=15,lwd=NA,vmax=NA,new=TRUE,
                       show.significance=TRUE,verbose=FALSE) {
  if(verbose) print("vis.trends")
  T <- calculate.trends(x,minlen=minlen,is=is,verbose=verbose)
  trends <- T$trends*10
  p <- T$p
  cols <- as.numeric(colnames(trends))
  rows <- as.numeric(rownames(trends))
  significant <- ifelse(p<pmax,trends,NA)
  
  ticks <- seq(1,length(cols),signif(length(cols)/10,1))
  if (is.na(lwd)) lwd <- max(3-0.05*length(cols),0.2)
  
  if (is.na(vmax) | vmax=="q995") vmax <- q995(abs(trends))
  if (vmax=="max") vmax <- max(abs(trends),na.rm=T)
  if (vmax<1) vmax <- signif(vmax,1)
  if (vmax>1) vmax <- signif(vmax,2)
  dv <- signif(vmax/8,1)
  v0 <- 0#signif(dv/2,1)
  vstep <- seq(v0,vmax,dv)
  vstep <- unique(c(-1*vstep,vstep))
  vstep <- vstep[order(vstep)]
  cticks <- vstep[2:length(vstep)]-dv/2
  
  #cstep <- colscal(n=length(vstep)-1,col="t2m")
  cmin <- rgb(239,138,98,maxColorValue=255) # blue
  cmid <- rgb(247,247,247,maxColorValue=255) # white
  cmax <- rgb(103,169,207,maxColorValue=255) # red
  rgb.palette <- colorRampPalette(c(cmax,cmid,cmin),space="rgb")
  cstep <- rgb.palette(n=length(vstep)-1)
  # Plot trend as color
  if (new) dev.new()
  image(cols,rows,t(trends),breaks=vstep,col=cstep,
        xlab='start year',ylab='length of period (years)',
        main=paste(c(varlabel," trend (",unitlabel,"/decade)"),collapse=""))
  
  trends.plus <- t(trends)
  trends.plus[trends.plus<max(vstep)] <- NA
  image(cols,rows,trends.plus,col=cstep[length(cstep)],add=TRUE)
  trends.minus <- t(trends)
  trends.minus[trends.minus>min(vstep)] <- NA
  image(cols,rows,trends.minus,col=cstep[1],add=TRUE)
  
  # Mark significant trends with dark borders
  if(show.significance) {
    if(verbose) print(paste("mark significant trends (p<",pmax,")",sep=""))
    i <- which((is.finite(t(p)) & t(p)<pmax))
    x <- array(sapply(cols,function(x) rep(x,nrow(p))),length(p))[i]
    y <- rep(rows,nrow(p))[i]
    matlines(rbind(x-1/2,x+1/2),rbind(y-1/2,y-1/2),col='black',lwd=lwd,lty=1)
    matlines(rbind(x-1/2,x+1/2),rbind(y+1/2,y+1/2),col='black',lwd=lwd,lty=1)
    matlines(rbind(x-1/2,x-1/2),rbind(y-1/2,y+1/2),col='black',lwd=lwd,lty=1)
    matlines(rbind(x+1/2,x+1/2),rbind(y-1/2,y+1/2),col='black',lwd=lwd,lty=1)
  }
  colbar(cticks,cstep,fig=c(0.85,0.9,0.65,0.85))
}

calculate.trends <- function(x,minlen=15,is=1,verbose=FALSE){
  # Calculate trends of time series x
  if(verbose) print("calculate.trends - calculate trends for all subperiods")
  stopifnot(inherits(x,'zoo'))
  if(!is.null(dim(x))) {
    x <- subset(x,is=is)
    if(!is.null(dim(x))) {
      x <- apply(x,2,mean,na.rm=TRUE)
    }
  }
  if(!inherits(x,c("annual","season"))) {
    xm <- aggregate(x,by=as.yearmon(index(x)),FUN="mean")
    xy <- aggregate(xm,by=strftime(index(xm),"%Y"),FUN="mean")
    ny <- aggregate(xm,by=strftime(index(xm),"%Y"),FUN="nv")
    xy <- xy[ny==max(ny)] # exclude years with missing months
  } else xy <- x
  year <- as.numeric(index(xy))
  firstyear <- min(year):(max(year)-minlen+1)
  trendlen <- minlen:(diff(range(year))+1)
  #lastyear <- firstyear+minlen-1
  n <- length(firstyear)
  trends <- matrix(NA,n,n)
  rownames(trends) <- trendlen#firstyear
  colnames(trends) <- firstyear#lastyear
  p <- trends
  # speed up with apply?
  for (i in firstyear) {
    jvec <- i+trendlen[trendlen<=(max(year)-i+1)]-1#(i+minlen-1):(max(year)+1)
    for (j in jvec) {
      if(!is.na(xy[year==i]) & !is.na(xy[year==j])) {
        #if(verbose) print(paste(i,j))
        ij <- which(year %in% i:j & !is.na(xy))
        ij.model <- lm(xy[ij]~year[ij])
        #ij.kendall <- Kendall(x[ij],year[ij])
        iout <- firstyear==i#as.numeric(colnames(trends))==j
        jout <- trendlen==(j-i+1)#as.numeric(rownames(trends))==i
        trends[jout,iout] <- ij.model$coefficients[2]
        p[jout,iout] <- anova(ij.model)$Pr[1]#ij.kendall$sl[1]
      }
    }
  }  
  return(list("trends"=trends,"p"=p))
}
