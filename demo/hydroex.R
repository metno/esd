library(esd)
# VIKEDAL (met station no. 46850, LAT 59.5558N, LON 5.9955E
hydroIO <- read.table("~/Dropbox/data/CERAD_to RASMUS METNO.csv",sep=";",header=TRUE)
loc <- 'Vikedal'
lon <- 5.9955
lat <- 59.5558
alt <- 159
stid <- 46850

P <- zoo(x=hydroIO$Precipitation..mm,order.by=as.Date(hydroIO$Date))
t2m <- zoo(x=hydroIO$Temperature.C,order.by=as.Date(hydroIO$Date))
u <- zoo(x=hydroIO$Flow..m3.s,order.by=as.Date(hydroIO$Date))

par(bty="n")
plot(merge(P,u,as.anomaly(t2m)),plot.type='single',col=c("steelblue","black","red"),lwd=c(3,1,2))
dev2bitmap('hydroex-daily.png',res=150)

Pmm <- aggregate(P,as.yearmon,FUN='mean')
Tmm <- aggregate(t2m,as.yearmon,FUN='mean')
umm <- aggregate(u,as.yearmon,FUN='mean')

#plot(coredata(P),coredata(u))

cal <- data.frame(u=coredata(umm),t2m=coredata(Tmm),Pmm=coredata(Pmm))
fit <- lm(u ~ t2m + Pmm,data=cal)
print(summary(fit))

par(bty="n")
plot(coredata(umm),predict(fit),
     main='Monthly river flow predicted from monthly T and P',
     ylim=range(coredata(Pmm)),
     xlab=expression(u[hydro]),ylab=expression(u==f(T[2*m],P)))
lines(range(coredata(Pmm)),range(coredata(Pmm)),col="red",lty=2)
grid()
dev2bitmap(file='hydroex.png',res=150)

temp <- as.station(t2m,lon=lon,lat=lat,alt=alt,unit='deg C',param='t2m',loc=loc,stid=stid)

if (FALSE) {
z.t2m <- try(DSensemble.t2m(temp,biascorrect=TRUE,plot=TRUE))
dev2bitmap('hydroex-dse-z.t2m.png',res=150)
save(file='hydroex-dse-z.t2m.rda',z.t2m)
}

precip <- as.station(P,lon=lon,lat=lat,alt=alt,unit='mm/day',param='precip',loc=loc,stid=stid)
if (FALSE) {
z.mu <- try(DSensemble.precip(precip,biascorrect=TRUE,plot=TRUE))
dev2bitmap('hydroex-dse-z.mu.png',res=150)
save(file='hydroex-dse-z.mu.rda',z.mu)
z.fw <- try(DSensemble.precip(precip,reanalysis=slp,FUN='wetfreq',biascorrect=TRUE,plot=TRUE))
dev2bitmap('hydroex-dse-z.fw.png',res=150)
save(file='hydroex-dse-z.fw.rda',z.fw)

}

flow <- as.station(u,lon=lon,lat=lat,alt=alt,unit='m^3/s',param='u',loc=loc,stid=stid)

wheel(flow)

muam <- annual(precip,FUN='exceedance',threshold=1)
fwam <- annual(precip,FUN='wetfreq',threshold=1)
uamx <- annual(flow,FUN='count',threshold=30)
attr(uamx,'variable') <- 'n(X > m^3/s)'
attr(uamx,'unit') <- ' '

am <- data.frame(u=coredata(uamx),mu=coredata(muam),fw=coredata(fwam))
emu <- glm(u ~ mu + fw, data=am, family='poisson')

plot(uamx,lwd=3,ylim=c(0,70))
lines(zoo(exp(predict(emu)),order.by=index(uamx)))

text(1980,70,pos=4,'Number of annual events: u > 30 m^3/s')
legend(1980,65,c("hydrologic model",expression(f(mu,f[w]))),
       lty=1,col=c("red","black"),lwd=c(3,1),bty="n")
dev2bitmap('hydroex-u.gt.30.png',res=150)
