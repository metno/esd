## Infographics: Information and Uncertainty
## Rasmus.Benestad@met.no Oslo, Norway, 2017-03-02

ideal <- FALSE
hc <- TRUE
tb <- 0.5    # Transparency of the shape borders
tc <- 0.2   # Transparency of the shape colours
lwd <- 2     # Shape border thickness

## Use a new colour pallette and import it into the rgb to get transparent shapes
library(wesanderson)
wpal <- wes_palette("Zissou")
wr <- round(strtoi(paste('0x',substr(wpal[c(1,3,5)],2,3),sep=''))/255,2)
wg <- round(strtoi(paste('0x',substr(wpal[c(1,3,5)],4,5),sep=''))/255,2)
wb <- round(strtoi(paste('0x',substr(wpal[c(1,3,5)],6,7),sep=''))/255,2)

area <- function(x=0,y=0,col=rgb(0.5,0.5,0.5,0.1),border=rgb(0.5,0.5,0.5,0.3),lwd=1,
                 a=1,w1=0,w2=0,w3=0,w4=0,it=2*pi*(1:180)/180) {
  r <- a*(1+0.25*sin(w1*it) + 0.25*cos(w2*it) + 0.25*sin(w3*it) + 0.25*cos(w4*it))
  polygon(r*cos(it)+x,r*sin(it)+y,col=col,border=border,lwd=lwd)
}

if (ideal) filename <- '~/Desktop/uncertainty-ideal-ensemble.png' else
           filename <- '~/Desktop/uncertainty-adhoc-ensemble.png'
if (hc) {
  png(file=filename,bg='transparent',height = 3000,width = 3000)
  par(cex=4)
}
par(bty='n',xaxt='n',yaxt='n',xpd=TRUE)
if (ideal) main <- 'Ideal ensemble mapping of State Probability Space' else
           main <- 'Ad hoc ensemble mapping of State Probability Space'
plot(c(-5,5),c(-5,5),type='n',pch=19,col=rgb(0.75,0.75,0.75,0.1),cex=1.4,
     main=main,xlab='',ylab='',xlim=c(-5,5),ylim=c(-5,5))

for (a in seq(0.1,2.5,by=0.05))
  area(a=a,col=rgb(wr[1],wg[1],wb[1],0.03),border=rgb(wr[1],wg[1],wb[1],0.03))
x0 <- 0.5*rnorm(1)-0.1; y0 <- 0.5*rnorm(1)+0.2 ## The "truth"
points(x0,y0,cex=1.5,pch=19)

## GCMs
for (i in 1:25) {
  if (!ideal) {x <- 0.75*rnorm(1)+0.05; y <- 0.75*rnorm(1)+0.6; a1=0.65+0.2*rnorm(1); a2=0.4+0.1*rnorm(1)} else
              {x <- 0.75*(i %% 5)-1.5; y <- 0.75*(i %/% 5)-1.5; a1=1.0+0.4*rnorm(1); a2=1.0+0.4*rnorm(1)}
  area(x=x,y=y,col=rgb(wr[2],wg[2],wb[2],tc),border=rgb(wr[2],wg[2],wb[2],tb),lwd=lwd,
       a=a1,w1=0.1+0.01*rnorm(1),w2=2+0.01*rnorm(1),w3=4+0.01*rnorm(1),w4=6+0.01*rnorm(1))
  
  if (!ideal) area(x=x+0.1*rnorm(1),y=y+0.1*rnorm(1),col=rgb(wr[3],wg[3],wb[3],tc),
                   border=rgb(wr[3],wg[3],wb[3],tb),lwd=lwd,
                   a=a2,w1=1+0.01*rnorm(1),w2=2+0.01*rnorm(1),w3=4+0.01*rnorm(1),w4=8+0.01*rnorm(1))
  
}

area(x=-3.5,y=-3.25,col=rgb(wr[2],wg[2],wb[2],tc),border=rgb(wr[2],wg[2],wb[2],tb*1.5),lwd=lwd,
     a=0.65,w1=0.1,w2=2,w3=4,w4=6)
text(x=-4.5,y=-3.25,'GCM',pos=4,cex=1,font=2,col=rgb(0.5*wr[2],0.5*wg[2],0.5*wb[2],0.9))
if (!ideal) {
  area(x=-3.5,y=-4.75,
     col=rgb(wr[3],wg[3],wb[3],tc),border=rgb(wr[3],wg[3],wb[3],tb*1.5),lwd=lwd,
     a=0.45,w1=1,w2=2,w3=4,w4=8)
  text(x=-4.5,y=-4.75,'Downscaled',pos=4,cex=1,font=2,col=rgb(0.5*wr[3],0.5*wg[3],0.5*wb[3],0.9))
}

points(x0,y0,cex=1.5,pch=19,col=rgb(0,0,0,0.3))
if (hc) dev.off()
