## Check Global sea-level: combine series from different sources
## in order to get an updated series dating back to the early records

library(esd)
gsl1 <- annual(GSL())
gsl2 <- annual(GSL.aviso())
index(gsl1) <- year(gsl1)
index(gsl2) <- year(gsl2)
gsl12 <- merge(gsl1,gsl2,all=FALSE)
caldat <- data.frame(y=coredata(gsl12$gsl1),x=coredata(gsl12$gsl2))
xdat <- data.frame(x=gsl2)
fit <- lm(y ~ x, data=caldat)
xgsl <- zoo(predict(fit,newdata=xdat),order.by=index(gsl2))
gsl1 <- gsl1[is.finite(gsl1)]
gsl <- c(gsl1,window(xgsl,start=end(gsl1)+1))
plot(gsl,lwd=3,col='grey')
lines(gsl1)
lines(xgsl,col='red',lty=2)
         