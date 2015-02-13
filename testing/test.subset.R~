require(esd)

it=c(1960,1980)
is <- list(lon=c(-40,30),lat=c(50,70))
T2m.era40 <- t2m.ERA40()
T2m.era40 <- subset(T2m.era40,it=it,is=is,verbose=TRUE)
#T2m.era40 <- subset(T2m.era40,it=it)
#T2m.era40 <- subset(T2m.era40,is=is)

data(ferder)
y <- subset(ferder,it='jul')

NACD <- station(src='nacd')
Y <- subset(NACD,it=c(1930,1990),is=1:10)
