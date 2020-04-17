## Demo showing how to analyse the number of record-breaking annual rainfall amounts simulated by a random CMIP6 
## global cliamtem odel simulation. @RasmusBenestad. 2020-04-03
library(esd)
meta <- meta.ESGF(param='pr')
pr <- retrieve.ESGF(im=2,meta=meta)
PR <- annual(pr,FUN='sum')
n0 <- 0; for (i in 1:1000) n0 <- n0 + n.records(rnorm(n))$N
n0 <- n0/1000
nrec <- function(x,n=5.095,na.rm=FALSE) 100*as.numeric(n.records(x)$N)/n
map(PR,FUN='nrec',colbar=list(pal='t2m',rev=TRUE))
