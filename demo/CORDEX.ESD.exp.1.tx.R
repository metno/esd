# Script for setting up and running CORDEX ESD experiment 1
# 

# load the predictands: CLARIS precip
load('~/Dropbox/Public/CORDEX-ESDM-data-clumps/claris.Tx.rda')

# Temperature
t2m <- retrieve('data/ERAINT/ERAINT_t2m_mon.nc',lon=c(-80,-40),lat=c(-30,-20))
t2m <- annual(t2m)
eof.t2m <- EOF(t2m)

# Out-going long-wave radiation
olr <- annual(retrieve('data/ERAINT/ERAINT_olr_mon.nc',
                       lon=c(-90,-30),lat=c(-35,-15)),FUN='mean')
attr(olr,'unit') <- "W * m**-2"
eof.olr <- EOF(olr)

# Mean sea level pressure
slp <- annual(retrieve('data/ERAINT/ERAINT_slp_mon.nc',
                       lon=c(-90,-30),lat=c(-35,-15)),FUN='mean')
eof.slp <- EOF(slp)


# Process the precipitation - predictand:
amt <- annual(Tx,FUN='mean')
ast <- annual(Tx,FUN='sd')
pca.amt <- PCA(amt)
pca.ast <- PCA(ast)

# The experiment results are provided by the cross-validation 
z.amt <- DS(pca.amt,eof.t2m,
            m='cordex-esd-exp1',detrend=FALSE,verbose=TRUE)
tx.ds <- pca2station(z.amt)

z.ast <- DS(pca.ast,list(t2m=eof.t2m,slp=eof.slp,olr=eof.olr),
            m='cordex-esd-exp1',detrend=FALSE,verbose=TRUE)
txs.ds <- pca2station(z.ast)

# Extract the predicted values:
d <- dim(attr(tx.ds,'evaluation'))
exp1.tx <- attr(tx.ds,'evaluation')[,seq(2,d[2],by=2)]
exp1.tx <- attrcp(amt,exp1.tx)
class(exp1.tx) <- class(amt)
attr(exp1.tx,'description') <- 'cross-validation'
attr(exp1.tx,'experiment') <- 'CORDEX ESD experiment 1'
attr(exp1.tx,'method') <- 'esd'

i <- 10
plot(subset(amt,is=i),xlim=c(1979,2010))
lines(subset(pca2station(pca.amt),is=i),lwd=2,col="grey")
lines(zoo(subset(exp1.tx,is=i),order.by=year(exp1.tx)),lwd=2,col="black")
