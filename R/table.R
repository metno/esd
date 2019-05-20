#' Summary showing summary of objects %% ~~function to do ... ~~
#' 
#' Produce a summary table
#' 
#' 
#' @aliases summary.dsensemble summary.station summary.ds summary.eof
#' summary.cca
#' @param x an object of type 'DSensemble', 'station', 'DS', 'EOF' or 'CCA'
#' @param years A set of years for which to produce summary statistics
#' @param im The order of months in the table. Use \code{im=c(10:12,1:9)} for
#' Oct-Sep.
#' @return A matrix containing summary statistics
#' @examples
#' 
#' data("Oslo")
#' summary(Oslo)
#' 
#' @export summary.dsensemble
summary.dsensemble <- function(x,years=seq(1990,2090,by=20)) {
    x0 <- subset(x,it=0)
    djf <- subset(x,it='djf')
    mam <- subset(x,it='mam')
    jja <- subset(x,it='jja')
    son <- subset(x,it='son')

    tab <- rep('',length(years) + 1)
    tab[1] <- paste(loc(x), '  Annual, DFJ, MAM, JJA, SON')
    i <- 1
    ##browser()
    for (yr in years) {
      i <- i + 1
      ##print(i,yr)
      tab[i] <- paste(years[i-1],':    ',
                 round(mean(coredata(subset(x0,it=years[i-1]))),2),
                 ' [',round(quantile(subset(x0,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(x0,it=years[i-1]),0.95),2),'],  ',
                 round(mean(coredata(subset(djf,it=years[i-1]))),2),
                 ' [',round(quantile(subset(djf,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(djf,it=years[i-1]),0.95),2),'],  ',
                 round(mean(coredata(subset(mam,it=years[i-1]))),2),
                 ' [',round(quantile(subset(mam,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(mam,it=years[i-1]),0.95),2),'],  ',
                 round(mean(coredata(subset(jja,it=years[i-1]))),2),
                 ' [',round(quantile(subset(jja,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(jja,it=years[i-1]),0.95),2),'],  ',
                 round(mean(coredata(subset(son,it=years[i-1]))),2),
                 ' [',round(quantile(subset(son,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(son,it=years[i-1]),0.95),2),']',sep='')
      ##print(tab[i])
    }
    tab
}

summary.station <- function(x,im=1:12) {
  tab <- matrix(rep(NA,12*7),12,7)
  for (i in 1:12) {
    y <- subset(x,it=month.abb[i])
    z <- as.numeric(summary(coredata(y)))
    attributes(z) <- NULL
    tab[i,1:length(z)] <-z 
  }
  attn <- attr(summary(coredata(x)),'names')
  if (length(attn)==6) colnames(tab) <- c(attn,"NA's") else colnames(tab) <- attn
  rownames(tab) <- month.abb
  tab[im,]  
}

summary.ds <- function(x) {
}

summary.eof <- function(x) {
}

summary.cca <- function(x) {
}
