# test.retrieve.R

library(esd)
path.data <- "/path/to/file"
file1 <- file.path(path.data,"filename.nc")

file <- file1
x <- retrieve(file, verbose=TRUE)
range(index(x))
class(x)
class(index(x))
map(x)
xsub <- subset(x, it=NULL, is=list(lon=c(-30,30),lat=c(35,75)))
xm <- as.monthly(xsub)
xeof <- EOF(xsub)
plot(xeof)
y <- aggregate.area(xm,FUN='mean')
plot(y)
