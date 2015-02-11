
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

summary.station <- function(x) {
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
  tab  
}

summary.ds <- function(x) {
}

summary.eof <- function(x) {
}

summary.cca <- function(x) {
}
