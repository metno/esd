#' combined maps with northern and southern stereographic projections respectively.
#' 
#' \code{as.annual} and \code{annual} aggregates time series into annual values (e.g. means).
#'
#' \code{as.monthly} aggregates time series into monthly values (e.g. means).
#'
#' \code{as.daily} aggregates time series into daily values (e.g. means).
#'
#' \code{showmaps} presents a montage of two maps with stereographic projections, one for the northa nd one for the south. 
#' @aliases showmaps
#' @seealso map
#' @param x a 'field' object
#' @param FUN a function, see \code{\link{aggregate.zoo}}
#' @param colbar is a list - see \code{\link{map}}
#' @param nnum how many numbers to be shown along the colour bar
#' 
#' @return list with map output
#'
#' @keywords utilities
#'
#' @examples
#' t2m <- t2m.NCEP()
#' period <- paste(range(year(t2m)),collapse=' - ')
#' attr(t2m,'sub') <- paste0('Average over: ',period,'; source: NCEP')
#' showmaps(t2m,FUN='mean',colbar=list(breaks=seq(-60,20,by=5),show=FALSE))
#'
#' attr(t2m,'sub') <- paste0('Variability (stdv) over: ',period,'; source: NCEP')
#' showmaps(t2m,FUN='sd')
#' 
#' attr(t2m,'sub') <- paste0('Trend over: ',period,'; source: NCEP')
#' showmaps(t2m,FUN='trend',colbar=list(breaks=seq(-1,1,by=0.1),show=FALSE))
#' @export
showmaps <- function(x,FUN='mean',colbar=list(pal='bwr',show=FALSE),nnum=11, ...) { 
  z <- map(x,FUN=FUN,plot=FALSE)
  def.par <- par(no.readonly = TRUE)
  par(mar=c(0.5,0.2,0.5,0.2),xaxt='n',yaxt='n',xpd=TRUE)
  nf <- layout(matrix(c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,0,0,3,3,0,0),6,6,byrow = TRUE), 
               widths=rep(1,36), heights=rep(0.5,36), respect= TRUE)
  
  z1 <- map(z,projection='np',colbar=colbar, ...)
  ## Second plot with same colour scale as previous one
  attr(z,'variable') <- NULL
  z2 <- map(z,projection='sp',colbar=attr(z1,'colbar'), ...)
  if (!is.null(attr(x,'sub'))) text(-1,-1.05,attr(x,'sub'),pos=4)
  ## Use same colour scale as in the plots
  breaks <- attr(z1,'colbar')$breaks
  n <- length(breaks)
  mids <- round(0.5 * (breaks[-1] + breaks[-n]),2)
  par(mar=c(0.25,0.25,0.25,0.25),xaxt='n',yaxt='n')
  plot(mids,rep(0.5,n-1),col=attr(z1,'colbar')$col,pch=19,cex=3,ylim=c(0,1))
  ii <- (1:n)%%(round(n/nnum)) == 1; ii[n] <- TRUE
  text(breaks[ii],rep(0.2,n)[ii],round(breaks,2)[ii],cex=0.7, col='grey30')
  par(def.par)
  invisible(list(z1 = z1, z2 = z2))
}

#t2m <- retrieve('~/data/ERA5/ERA5_t2m_year.nc',it=c(1980,2022))
#t2m <- t2m.NCEP()
# attr(t2m,'sub') <- 'Gjennomsnitt over: 1980-2022; kilde: ERA5'
# showmaps(t2m,FUN='mean',colbar=list(breaks=seq(-60,20,by=5),show=FALSE))

#attr(t2m,'sub') <- 'Trend (C/tiår) over: 1950-2022; kilde: ERA5'
#showmaps(t2m,FUN='trend',colbar=list(breaks=seq(-1,1,by=0.1),show=FALSE))