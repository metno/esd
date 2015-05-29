# K Parding, 29.05.2015

library(esd)
slp <- slp.ERAINT()

CCI <- function(Z,m=14,nsim=10,it=NULL,is=NULL,
                cyclones=TRUE,accuracy=NULL,verbose=FALSE) {

  stopifnot(inherits(Z,'field'))
  Z <- subset(Z,it=it,is=is)
  resx <- dX(Z,m=m,accuracy=accuracy,verbose=verbose)
  resy <- dY(Z,m=m,accuracy=accuracy,verbose=verbose)

  dslpdx <- resx$dZ
  dslpdx2 <- resx$dZ2
  dslpdy <- resy$dZ
  dslpdy2 <- resy$dZ2
  wind <- dslpdx*0
  
  
#CCI <- function(m=14,lplot=TRUE,nsim=10,fname="data/cyclones.Rdata",
#                fielddata="data/nmc_slp.nc",vname="slp",cyclones=TRUE,
#                force365.25=FALSE,x.rng=c(-180,180),y.rng=c(30,90),tslice=365,
#                rad=NULL,dx=NULL,dy=NULL,times=NULL,label=NULL,rho=1.293,
#                nc.read.method="retrieve.nc",l.scale=TRUE,
#                graph.dir="CCI.graphs/", greenwich = TRUE, 
#                plot.interval=50,EPS=TRUE,verbose=TRUE,accuracy=NULL) {
## if (substr(graph.dir,nchar(graph.dir),nchar(graph.dir))=="/")
##    graph.dir<- substr(graph.dir,1,nchar(graph.dir)-1)
## if ( (graph.dir!="./") & (lplot) &  !file.exists(graph.dir) ) {
##   print(paste("Create new directory:",graph.dir))
##   dir.create( graph.dir )
## }
## graph.dir<- paste(graph.dir,"/",sep="")

## slash <- instring("/",fielddata); slash <- slash[length(slash)]
## back.up.file <- paste("CCI_backup_",substr(fielddata,slash+1,nchar(fielddata)-4),".Rdata",sep="")

##if (verbose) print(paste("Input file:",fielddata," Back-up file:",back.up.file)) 

  ## writeLines("CCI is running",con=".CCI.run")
  ## if (verbose) print(paste("CYCLONES: region ",x.rng[1],"-",x.rng[2],
  ##                        "E/",y.rng[1],"-",y.rng[2],
  ##                         "N. res=", dx,"x",dy,sep=""))
  ## if (verbose) print("run stopCCI() in the same running directory to stop the process")
##if (verbose) {print(times); print(file.exists(fielddata)); print(fielddata)}

## if (file.exists(fielddata)) {
##   if (verbose) print(paste("Checking",fielddata,"for updating..."))
##   ncid <- open.ncdf(fielddata)
##   nv <- ncid$nvars; 
##   cdfvars <- rep("-",nv) 
##   for (i in 1:nv) cdfvars[i] <- ncid$var[[i]]$name
##   if (verbose) print(cdfvars)
##   ipick <- grep(vname,cdfvars)
##   #print(ipick)
##   nd <- ncid$var[[ipick]]$ndims
##   cdfdims <- rep("-",nd)
##   for (i in 1:nd) cdfdims[i] <- ncid$var[[ipick]]$dim[[i]]$name
##   if (verbose) print("Dimensions:");print(cdfdims)  
##   itim <- grep("tim",lower.case(substr(cdfdims,1,3)))
##   tim <- get.var.ncdf(ncid,cdfdims[itim])
##   if (verbose) print(range(tim))
##   arv <- att.get.ncdf(ncid, cdfdims[itim], 'time_unit')
##     if( arv$hasatt ) t.unit <- arv$value else {
##        arv <- att.get.ncdf(ncid, cdfdims[itim], 'unit')
##        if( arv$hasatt ) t.unit <- arv$value else {
##          arv <- att.get.ncdf(ncid, cdfdims[itim], 'units')
##          if( arv$hasatt ) t.unit <- arv$value else t.unit <- NULL
##        } 
##     }
##   if (verbose) print(paste('time unit: ',t.unit))

##   test.tunit <- grep('since',t.unit)
##   if (length(test.tunit)==0) test.tunit <- c(0)
  
##   if (test.tunit[1]==0) {
##     arv <- att.get.ncdf(ncid, cdfdims[itim], 'time_origin')
##       if( arv$hasatt ) t.org <-datestr2num(arv$value) else {
##          arv <- att.get.ncdf(ncid, cdfdims[itim], 'origin')
##          if( arv$hasatt ) t.org <- datestr2num(arv$value) else t.org <- c(1900,01,01)
##        }
##   } else {
##     pos <- instring('since',lower.case(t.unit))
##     t.org <- datestr2num(substr(t.unit,pos+6,pos+15))
##   }
##   if (verbose) print(paste("Time origin:  year=",t.org[1],'  month=', t.org[2],'  day=',t.org[3]))

##   if (substr(lower.case(t.unit),1,4)=="hour") {
##     if (verbose) print("Dividing by 24 to get time units in days")
##     tim <- tim/24
##     if (verbose) print(range(tim))
##   }
##   close.ncdf(ncid)
## # REB 25.01.2006  NT <- length(tim)
## # REB 25.01.2006    times <- seq(1,NT,by=tslice)
##   if (is.null(times)) times <- tim
## # REB 25.01.2006   if (max(times) < NT) times <- c(times,NT)
## # REB 25.01.2006   print("tim:"); print(NT); print(range(tim))
## } else stop(paste('Could not find',fielddata))

# REB 25.01.2006 NT <- max(times)
##   times <- index(Z)
##   NT <- length(times)+1
##   times.physical <- times
##   max.tim <- length(times.physical)
##   if (is.null(label)) label <- paste(fielddata,": ",vname,sep="")
##   if (verbose) print(paste("label=",label))
##   lon <-       matrix(rep(NA,NT*nsim),NT,nsim)
##   lat <-       matrix(rep(NA,NT*nsim),NT,nsim)
##   tim <-       matrix(rep(0,NT*nsim),NT,nsim)
##   yy <-        matrix(rep(NA,NT*nsim),NT,nsim)
##   mm <-        matrix(rep(NA,NT*nsim),NT,nsim)
##   dd <-        matrix(rep(NA,NT*nsim),NT,nsim)
##   psl <-       matrix(rep(NA,NT*nsim),NT,nsim)
##   max.dpsl <-  matrix(rep(NA,NT*nsim),NT,nsim)
##   max.speed <- matrix(rep(NA,NT*nsim),NT,nsim)
##   radius <-    matrix(rep(NA,NT*nsim),NT,nsim)

##   if (cyclones) my.col <- rgb(c(c(seq(0.4,1,length=20),rep(1,21))),
##                             c(c(seq(0.4,1,length=20),rep(1,21))),
##                             c(c(seq(0.4,1,length=20),rep(1,21)))) else
##                 my.col <- rgb(c(c(rep(1,21),seq(1,0.4,length=20))),
##                             c(c(rep(1,21),seq(1,0.4,length=20))),
##                             c(c(rep(1,21),seq(1,0.4,length=20))))

##   i.max <- 0; is0 <- 1
##   plot.now <- TRUE
##   if (file.exists(fname)) {
##     if (verbose) print(paste("Updating",fname))
##     load(fname)
##     ii <- results$i
##     if (length(ii) > NT) {
##       NT <- length(ii)
##       if (verbose) print(paste("NT=",NT))
##       lon <-       matrix(rep(NA,NT*nsim),NT,nsim)
##       lat <-       matrix(rep(NA,NT*nsim),NT,nsim)
##       tim <-       matrix(rep(0,NT*nsim),NT,nsim)
##       yy <-        matrix(rep(NA,NT*nsim),NT,nsim)
##       mm <-        matrix(rep(NA,NT*nsim),NT,nsim)
##       dd <-        matrix(rep(NA,NT*nsim),NT,nsim)
##       psl <-       matrix(rep(NA,NT*nsim),NT,nsim)
##       max.dpsl <-  matrix(rep(NA,NT*nsim),NT,nsim)
##       max.speed <- matrix(rep(NA,NT*nsim),NT,nsim)
##       radius <-    matrix(rep(NA,NT*nsim),NT,nsim)
##     }
  
##   if (verbose) print(paste("length(ii)=",length(ii)," dim(lon)=",dim(lon)[1],"x",
##               dim(lon)[2]," dim(results$lon)=",dim(results$lon)[1],"x",dim(results$lon)[2]," NT=",NT))
##   if (verbose) print(summary(ii))
##   if (verbose) print(c(dim(lon[ii,]),NA,dim(results$lon)))
##   lon[ii,] <- results$lon[ii,]
##   lat[ii,] <- results$lat[ii,]
##   psl[ii,] <- results$psl[ii,]
##   tim[ii] <- results$tim[ii]
##   yy[ii] <- results$yy[ii]
##   mm[ii] <- results$mm[ii]
##   dd[ii] <- results$dd[ii]
##   i.max <- max(results$i,na.rm=TRUE)
##   max.dpsl[ii,] <- results$max.dpsl[ii,]
##   max.speed[ii,] <- results$max.speed[ii,]
##   radius[ii,] <- results$radius[ii,]
  
## # REB 25.01.2006  times <- c(i.max+1,times[times >= i.max])
##   times <- c(max(tim[ii],na.rm=TRUE)+1,times[times >= max(tim[ii],na.rm=TRUE)])
##   if (verbose) print(paste("length(times)=",length(times)))
##   if (verbose) print(summary(times))
##   if (verbose) print(summary(tim[ii]))
   
##   Nx <- sum(times > max(tim[ii],na.rm=TRUE))
##   if (Nx > 0) {
##     if (verbose) print(paste("Detected a need for extending the matrices",
##                              "max(tim[ii])=",max(tim[ii],na.rm=TRUE),
##                              "max(times)=",max(times,na.rm=TRUE)) )
##     lon <- rbind(lon[ii,],matrix(rep(NA,Nx*nsim),Nx,nsim))
##     lat <- rbind(lat[ii,],matrix(rep(NA,Nx*nsim),Nx,nsim))
##     psl <- rbind(psl[ii,],matrix(rep(NA,Nx*nsim),Nx,nsim))
##     radius <- rbind(radius[ii,],matrix(rep(NA,Nx*nsim),Nx,nsim))
##     max.dpsl <- rbind(max.dpsl[ii,],matrix(rep(NA,Nx*nsim),Nx,nsim))
##     max.speed <- rbind(max.speed[ii,],matrix(rep(NA,Nx*nsim),Nx,nsim))
##     tim <- c(tim[ii],rep(0,Nx)) 
##     yy <- c(yy[ii],rep(NA,Nx))
##     mm <- c(mm[ii],rep(NA,Nx))
##     dd <- c(dd[ii],rep(NA,Nx))
##     if (verbose) print(paste("Extending the vectors by ",Nx," elements: dim(psl)= ",dim(psl)[1],", ",dim(psl)[2],sep=""))
##   } 

##   if (verbose) print(paste( "Time in file=",min(tim[ii],na.rm=TRUE),"-",
##                max(tim,na.rm=TRUE)," i.max=",i.max," length(i)=",length(results$i),
##                "length(tim[ii])=",length(tim[ii]),
##                "min(tim[ii])=",min(tim[ii]),"max(tim[ii])=",max(tim[ii]),
##                " extend by ",Nx))
## } else {
##   if (verbose) print(paste(fname,'does not exist...'))
##   tim <- c(0)
##   j.org <- julday(t.org[2],t.org[3],t.org[1])
##   dates <- caldat(tim + j.org)
##   yy <- dates$year
##   mm <- dates$month
##   dd <- dates$day
##   i.max <- 1
##   rm(dates); gc(reset=TRUE)
## }

## ii <- i.max
## if (length(times)>1) time.res <- round(mean(diff(times),na.rm=TRUE),4) else
## if (length(tim[1:ii])>1) time.res <- round(mean(diff(tim[1:ii]),na.rm=TRUE),4) else
##                        time.res <- NA
## if (max.tim==0) times.physical <- tim[1:ii]
## max.tim <- length(times.physical)

## if (!is.finite(1/time.res)) {
##   print("Problem with time information!")
##   if (length(times)>1) {
##     print("times:")
##     print(summary(times))
##     print(summary(diff(times)))
##   } else {
##     print("tim[ii]:")
##     print(summary(tim[1:ii]))
##     print(summary(diff(tim[1:ii])))
##   }     
##   stop("Error in CCI: 1/time.res is not a finite number")
## }
## if (time.res<0) {
##   print("Problem with time information!")
##   if (length(times)>1) {
##     print(summary(times))
##     print(summary(diff(times)))
##   } else {
##     print(summary(tim[1:ii]))
##     print(summary(diff(tim[1:ii])))
##   }  
##   stop("Error in CCI: time.res is negative!")
## }
## if (verbose) print(paste("Time resolution=",1/time.res," times a day: time.res=",time.res,
##                          "max.tim=",max.tim))
## if (verbose) print(summary(tim))

## # Main loop

## #max.tim <- max(times,na.rm=TRUE)
## #times <- times[seq(1,length(times),by=tslice)]
## #if (max(times) < max.tim) times <- c(times,max.tim)

## # Now, times represent the time index rather than the variable 'tim' (unit eg day/hour)
## i.here <- min(seq(1,max.tim,by=1)[times.physical > tim[ii]])

## if (verbose) print(paste("i.max=",i.max,"   i.here=",i.here,"  max.tim=",max.tim,"   tslice=",tslice))
## if (verbose) print(paste("length(times)=",length(times)," min(times)=",min(times)," max(times)=",max(times),"times[1]=",times[1]))

## if ((length(times)> 1) & is.finite(i.here)) {

## times <- seq(i.here,max.tim,by=tslice)
## if (max(times) < max.tim) times <- c(times,max.tim)  # ensure to include the last time slice
## #print(times)
  
## for(is in 1:(length(times)-1)) {
## if (verbose) print(paste("is=",is,"times[is]=",times[is],"times[is+1]-1=",times[is+1]-1))
## slp <- eval(parse(text=paste(nc.read.method,
##           "(fielddata,vname,x.rng=x.rng,y.rng=y.rng,l.scale=",l.scale,
##           ",t.rng=c(times[is],times[is+1]-1),force365.25=force365.25,
##            greenwich=greenwich)",sep="")))
## nx <- length(slp$lon)
## dlon <- slp$lon[2] - slp$lon[1]
## ny <- length(slp$lat)
## nt <- length(slp$tim)
## NT <- max(slp$tim) - min(tim)+1

## if (verbose) print(paste("nt=",nt," range(slp$tim)=",min(slp$tim),"-",max(slp$tim),"; max(tim)=",max(tim),"; NT=",NT))
## if ( (NT > nt) & (file.exists(fname)) ) {
##   Nx <- NT - nt
##   lon <-      rbind(lon,matrix(rep(NA,Nx*nsim),Nx,nsim))
##   lat <-      rbind(lat,matrix(rep(NA,Nx*nsim),Nx,nsim))
##   psl <-      rbind(psl,matrix(rep(NA,Nx*nsim),Nx,nsim))
##   radius <-   rbind(radius,matrix(rep(NA,Nx*nsim),Nx,nsim))
##   max.dpsl <- rbind(max.dpsl,matrix(rep(NA,Nx*nsim),Nx,nsim))
##   max.speed <- rbind(max.speed,matrix(rep(NA,Nx*nsim),Nx,nsim))
##   tim <- c(tim,rep(0,Nx)) 
##   yy <- c(yy,rep(NA,Nx))
##   mm <- c(mm,rep(NA,Nx))
##   dd <- c(dd,rep(NA,Nx))
##   print(paste("Extending the vectors by ",Nx," elements: dim(psl)= ",dim(psl)[1],", ",dim(psl)[2],sep=""))
## } else if ( (NT < nt) & (file.exists(fname))  ){
##   if (verbose) print("### NT <= nt ###!")
##   if (verbose) print("Did not extend vectors")
## }

## #print(c(range(tim),NA,range(slp$tim)))

## jday <- julday(slp$mm,slp$dd,slp$yy)
## i.max <- seq(1,nt,by=1)[slp$tim > max(tim,na.rm=TRUE)][1]

## if (!is.null(dx) & !is.null(dy)) {
##   Interpolate <- TRUE
##   print(paste("Interpolate to grid with",dx,"x",dy," degreee spacing.nt=",nt," i.max=",i.max))
##   lonx <- seq(min(slp$lon),max(slp$lon),by=dx)
##   latx <- seq(min(slp$lat),max(slp$lat),by=dy)
## } else { 
##   Interpolate <- FALSE
##   print(paste("Use original grid with.nt=",nt," i.max=",i.max))
##   lonx <- slp$lon; latx <- slp$lat
## }
## nxx <- length(lonx)
## nyx <- length(latx)
## lon.xy <- rep(slp$lon,ny)
## lat.xy <- sort(rep(slp$lat,nx))

## mslp <- meanField(slp)
## if (Interpolate) {
##   print("interpolate mean values:")
##   mslpmap <- interp(lon.xy,lat.xy,mslp$map,lonx,latx)$z
##   print(dim(mslpmap))
## } else {
##   if (verbose) print("Mean values:")
##   mslpmap <- mslp$map
## }
## mslp$map[mslp$map > 1100] <- -9999
## #postscript(file = paste("cyclones_mean-slp",min(slp$tim),"-",max(slp$tim),".eps",sep=""),
## #            onefile=FALSE,horizontal=FALSE)
## #map(mslp,levels=seq(980,1030,by=5))
## #dev.off()

## #print(paste("TEST! sum(slp$tim > max(tim))>0",(sum(slp$tim > max(tim))>0)))
## #print(sum(is.na(mslp$map)))

## if (verbose) print("---")
## if (verbose) print(paste("i.max=",i.max," nt=",nt," ii=",ii,
##                          "  sum(slp$tim > max(tim,na.rm=TRUE))=",sum(slp$tim > max(tim,na.rm=TRUE))))

## if (sum(slp$tim > max(tim,na.rm=TRUE))>0) {
##  for (it in i.max:nt) {
##   ii <- ii+1
##   if (mod(ii,plot.interval)==0) plot.now <- TRUE
##   P.lowx <- matrix(rep(0,(nxx-1)*(nyx-1)),nxx-1,nyx-1); px <- P.lowx
##   P.lowy <- matrix(rep(0,(nxx-1)*(nyx-1)),nxx-1,nyx-1); py <- P.lowy
##   dpslx <- px; dpsly <- py
##   if (Interpolate) slpmap <- interp(lon.xy,lat.xy,t(slp$dat[it,,])-mslp$map,lonx,latx)$z else
##                    slpmap <- t(slp$dat[it,,])-mslp$map

  resx <- dX(lonx,latx,slpmap,maxhar=maxhar)
  resy <- dY(lonx,latx,slpmap,maxhar=maxhar)

  dslpdx <- resx$dZ
  resx2 <- dX(resx$lon,latx,resx$dZ,maxhar=maxhar)
  dslpdx2 <- resx2$dZ
  dslpdy <- resy$dZ
  resy2 <- dY(lonx,resy$lat,resy$dZ,maxhar=maxhar)
  dslpdy2 <- resy2$dZ
  wind <- dslpdx*0

  NXX <- length(resx$lon); NYY <- length(resy$lat)
  #print(c(nxx,nyx,NXX,NYY))
  lonXY <- rep(0.5*(resx$lon[2:NXX]+resx$lon[1:(NXX-1)]),nyx-1)
  latXY <- sort(rep(0.5*(resx$lat[2:NYY]+resx$lat[1:(NYY-1)]),nxx-1))

  # REB 24.11.2010:
  dim(lonXY) <- c(NXX-1,NYY-1)
  dim(latXY) <- c(NXX-1,NYY-1)
# Masking out the borders to avoid spurious border effects:
#  dslpdx[c(1,length(dslpdx))] <- NA
#  dslpdx2[c(1,length(dslpdx2))] <- NA
#  dslpdy[c(1,length(dslpdy))] <- NA
#  dslpdy[c(1,length(dslpdy2))] <- NA
  
# Search for zero-crossing points, estimate total pressure field interpolated 
# onto finer grid, and gradient:

  #print(c(dim(px),NA,dim(resy$Z.fit),NA,dim(dslpdx),NA,dim(dslpdx2)))
  for (i in 1:(nxx-1)) {
    dy11 <- 0.5*(dslpdy[i,2:nyx]+dslpdy[i+1,2:nyx])
    dy12 <- 0.5*(dslpdy[i,1:(nyx-1)]+dslpdy[i+1,1:(nyx-1)])
    dy21 <- 0.5*(dslpdy2[i,2:nyx]+dslpdy2[i+1,2:nyx])
    dy22 <- 0.5*(dslpdy2[i,1:(nyx-1)]+dslpdy2[i+1,1:(nyx-1)])
    px[i,] <- 0.25*(resy$Z.fit[i,2:nyx]+resy$Z.fit[i+1,2:nyx]+
                    resy$Z.fit[i,1:(nyx-1)]+resy$Z.fit[i+1,1:(nyx-1)]) +
              0.25*(mslpmap[i,2:nyx]+mslpmap[i+1,2:nyx]+
                    mslpmap[i,1:(nyx-1)]+mslpmap[i+1,1:(nyx-1)])
    dpslx[i,] <- 0.5* (dy11 + dy12)/10                       # scale by 100/1000: SLP in hPa -> Pa, but x in km.
    if (cyclones) i.low <- (dy11*dy12 < 0) & (dy21+dy22 > 0) & is.finite(dy11+dy12+dy21+dy22) else
                  i.low <- (dy11*dy12 < 0) & (dy21+dy22 < 0) & is.finite(dy11+dy12+dy21+dy22)
    if (sum(i.low)>0) P.lowy[i,i.low] <- 1
  }


  #x11(); image(dpslx,main="dpslx")
#  resx$Z.fit
  
  
  for (j in 1:(nyx-1)) {
    dx11 <- 0.5*(dslpdx[2:nxx,j]+dslpdx[2:nxx,j+1])
    dx12 <- 0.5*(dslpdx[1:(nxx-1),j]+dslpdx[1:(nxx-1),j+1])
    dx21 <- 0.5*(dslpdx2[2:nxx,j]+dslpdx2[2:nxx,j+1])
    dx22 <- 0.5*(dslpdx2[1:(nxx-1),j]+dslpdx2[1:(nxx-1),j+1])
#    py[,j] <- 0.25*(resx$Z.fit[2:nxx,j]+resx$Z.fit[2:nyx,j+1]+
#                    resx$Z.fit[1:(nxx-1),j]+resx$Z.fit[1:(nyx-1),j+1]) +
#              0.25*(mslpmap[2:nxx,j]+mslpmap[2:nyx,j+1]+
#                    mslpmap[1:(nxx-1),j]+mslpmap[1:(nyx-1),j+1])  
    py[,j] <- 0.25*(resx$Z.fit[2:nxx,j]+resx$Z.fit[2:nxx,j+1]+              # Bug corrected - REB 20.11.2008
                    resx$Z.fit[1:(nxx-1),j]+resx$Z.fit[1:(nxx-1),j+1]) +
              0.25*(mslpmap[2:nxx,j]+mslpmap[2:nxx,j+1]+
                    mslpmap[1:(nxx-1),j]+mslpmap[1:(nxx-1),j+1])  
    dpsly[,j] <- 0.5* (dx11 + dx12)/10                      # scale by  100/1000: SLP in hPa -> Pa, but x in km.
    f <- 0.000147*sin(pi*latx[j]/180)
    wind[,j] <- sqrt(dslpdy[,j]^2 + dslpdx[,j]^2)/(10*f*rho)   # 100/1000: SLP in hPa -> Pa, but x in km.
    if (cyclones) i.low <- (dx11*dx12 < 0) & (dx21+dx22 > 0) & is.finite(dx11+dx12+dx21+dx22) else
                  i.low <- (dx11*dx12 < 0) & (dx21+dx22 < 0) & is.finite(dx11+dx12+dx21+dx22)
    if (sum(i.low)>0) P.lowx[i.low,j] <- 1
  }

  #print("Low-pressure regions")
  #print(dim(P.lowy)); print(dim(P.lowx))
  #x11(); image(P.lowy,main="P.lowy")
  #x11(); image(P.lowx,main="P.lowx")
  lows <- (P.lowy & P.lowx)
  pcent <- 0.5*(px[lows]+py[lows])
  strength <- order(pcent)
  if (!cyclones) strength <- reverse(strength)

  i.sim <- min(c(sum(lows,na.rm=TRUE),nsim),na.rm=TRUE)

  # Tests to provide diagnostics
  test.lon.d <- dim(lon)
  test.lonXY.d <- dim(lonXY)
  test.lows.d <- dim(lows)
  test.strength.l <- length(strength)
  
#  if ( (ii > test.lonXY.d[1]) | (i.sim > test.lonXY.d[2]) |
#       (test.lows.d[1]>test.lonXY.d[1]) | (test.lows.d[2]>test.lonXY.d[2]) |
#       (i.sim > test.strength.l) ) {
  if ( (ii > dim(lon)[1]) | (i.sim > test.lonXY.d[2]) |
       (test.lows.d[1]>test.lonXY.d[1]) | (test.lows.d[2]>test.lonXY.d[2]) |
       (i.sim > test.strength.l) ) {
    print("CCI test indicates subscript out of bounds")
#    print(paste("ii=",ii,test.lonXY.d[1],"i.sim=",i.sim,test.lonXY.d[2],
#                "sum(lows)=",sum(lows),
#                "sum(strength)=",sum(strength)))
    print(paste("ii=",ii,dim(lon)[1],"i.sim=",i.sim,test.lonXY.d[2],
                "sum(lows)=",sum(lows),
                "sum(strength)=",sum(strength)))
    print("Skip this time slice!")
  } else if ( (sum(lows)>0) & (i.sim >=1) ){
    if (length(lon[ii,1:i.sim])==length(lonXY[lows][strength][1:i.sim])) {
      lon[ii,1:i.sim] <- lonXY[lows][strength][1:i.sim]
      lat[ii,1:i.sim] <- latXY[lows][strength][1:i.sim]
      psl[ii,1:i.sim] <- pcent[strength][1:i.sim]
    } else {
      print(paste("CCI warning: inconsistent vector lengths. length(lon[ii,1:i.sim])=",
                  length(lon[ii,1:i.sim]),"length(lonXY[lows][strength][1:i.sim])",
                  length(lonXY[lows][strength][1:i.sim]),"ii=",ii,"i.sim=",i.sim,
                  "sum(lows)=",sum(lows),"sum(strength)=",sum(strength)))
      print("Skip this time slice!")            
    }
    tim[ii] <- slp$tim[it]; yy[ii] <- slp$yy[it]
    mm[ii] <- slp$mm[it]; dd[ii] <- slp$dd[it]
  }
  
# Remove secondary cyclones near a deeper one (same cyclonic system):

  #print(paste("i.sim=",i.sim,"sum(lows)=",sum(lows)))
  del <- rep(FALSE,i.sim)  
  for (i in 1:(i.sim-1)) {
    if ( (length(lon[ii,i])==1 ) & (length(lat[ii,i])==1) ) { 
      d <- distAB(lon[ii,i],lat[ii,i],lon[ii,(i+1):i.sim],lat[ii,(i+1):i.sim])/1000
      del <- del | c(rep(FALSE,i),d < 600)
    }
  }

  psl[ii,del] <- NA
  strength <- order(psl[ii,])
  if (!cyclones) strength <- reverse(strength)
  lon[ii,] <- lon[ii,strength]; lat[ii,] <- lat[ii,strength]; psl[ii,] <- psl[ii,strength]
 
# Gradients, geostrophic windspeed, and radius

  latX <- 0.5*(latx[1:(nyx-1)] + latx[2:nyx]);   latXX <- 0.5*(latX[1:(nyx-2)] + latX[2:(nyx-1)])
  lonX <- 0.5*(lonx[1:(nxx-1)] + lonx[2:nxx]);   lonXX <- 0.5*(lonx[1:(nxx-2)] + lonx[2:(nxx-1)]) 
  dpsl <- sqrt(dpslx^2 + dpsly^2)/10                                     # 100/1000: SLP in hPa -> Pa, but x in km.
  for (i in 1:sum(is.finite(psl[ii,]))) {

# Find points of inflexion (2nd derivative==0) to estimate the storm radius
    #print("Find points of inflexion")
    vec <- dslpdy2[lonX==lon[ii,i],]
    p.infly <- latXX[vec[2:nxx]*vec[1:(nxx-1)]<=0]
    p.infly[abs(p.infly - lat[ii,i]) < 1] <- NA    
    p.low <- p.infly[p.infly < lat[ii,i]]; p.high <- p.infly[p.infly > lat[ii,i]]
    p.low <- reverse(sort(p.low)); p.high <- sort(p.high)
    p.infly <- c(p.low[1],p.high[1])
    if (i==1) {
      #print(c(lon[ii,i],lat[ii,i],p.infly))
      p.infly1 <- p.infly
      y.test <- vec; x.test <- latx
      ilat1 <- (1:nyx)[(latx >= lat[ii,i])][1]
      ilon1 <- (1:nxx)[(lonx >= lon[ii,i])][1]
    }

    vec <- dslpdx2[,latX==lat[ii,i]]
    p.inflx <- lonXX[vec[2:nyx]*vec[1:(nyx-1)]<=0]
    p.inflx[abs(p.inflx - lon[ii,i]) < 1] <- NA
    p.low <- p.inflx[p.inflx < lon[ii,i]]; p.high <- p.inflx[p.inflx > lon[ii,i]]
    p.low <- reverse(sort(p.low)); p.high <- sort(p.high)
    p.inflx <- c(p.low[1],p.high[1])

    #if (i==1) print(c(c(rep(lon[ii,i],2),p.inflx),NA,c(p.infly,rep(lat[ii,i],2))))

    # Bug fixed 20.11.2008
    # - c(rep(lon[ii,i],2),p.inflx) repaced by c(rep(lon[ii,i],length(p.infly)),p.inflx)
    # - c(p.infly,rep(lat[ii,i],2)) replaced by c(p.infly,rep(lat[ii,i],length(p.inflx))))
    #print("HERE");  print(c(i,ii)); print(lon[ii,i]); print(lat[ii,i]);
    #print(p.infly); print(p.inflx)
    if ( (length(lon[ii,i])==1 ) & (length(lat[ii,i])==1) ) {
      vec <- distAB(lon[ii,i],lat[ii,i],c(rep(lon[ii,i],length(p.inflx)),p.inflx),
                    c(p.infly,rep(lat[ii,i],length(p.inflx))))/1000
      vec[(vec < 10) | (vec > 1200)] <- NA
      radius[ii,i] <- min(vec,na.rm=TRUE)
    }
    
# Find speed at points of inflexion:
    #print("Find speed at points of inflexion:")
    #i.near <- ( sqrt((lonXY - lon[ii,i])^2 + (latXY - lat[ii,i])^2) < rad )
    ilat <- c((1:(nyx-1))[latX>=p.infly[1]][1],(1:(nyx-1))[latX>=p.infly[2]][1])
    ilon <- c((1:(nxx-1))[lonX>=p.inflx[1]][1],(1:(nxx-1))[lonX>=p.inflx[2]][1])
    vec <- c(dpsl[lonX==lon[ii,i],ilat], dpsl[ilon,latX==lat[ii,i]])
    max.dpsl[ii,i] <- max(vec[is.finite(vec)])
    vec <- c(wind[lonX==lon[ii,i],ilat], wind[ilon,latX==lat[ii,i]])
    max.speed[ii,i] <- max(vec[is.finite(vec)])
  }

# TEST
#  #print(c(lon[ii,1],lat[ii,1], max.dpsl[ii,1]))
#  #print(c(dim(dpsl),NA,length(lonXY),NA,length(latXY),length(lonx),length(latx)))
#  x11(); image(lonX,latX,dpsl); addland()
#  i.test<- 2; i.near <- ( sqrt((lonXY - lon[ii,i.test])^2 + (latXY - lat[ii,i.test])^2) < rad )
#  vec <- as.vector(dpsl)[i.near]
#  size<-dim(dpsl); dpsl <- as.vector(dpsl); dpsl[!i.near] <- NA; dim(dpsl) <- size
#  #print(c(dim(dpsl),NA,length(lonXY),NA,length(latXY),length(lonX),length(latX)))
#  points(lon[ii,i.test],lat[ii,i.test],pch=20,cex=1.5,col="white")
#  contour(lonX,latX,dpsl,add=TRUE)
#  dpsl <- as.vector(dpsl); dpsl[] <- NA; dpsl[latXY==60] <- 1; dim(dpsl) <- size
#  contour(lonX,latX,dpsl,add=TRUE,col="grey80",lty=2)
#  points(lon[ii,i.test],lat[ii,i.test],pch=20,col="darkblue")
#  stop('break')

  attr(lon,'units') <- 'degrees'
  attr(lat,'units') <- 'degrees'
  attr(psl,'units') <- 'hPa'
  attr(max.dpsl,'units') <- 'Pa/m'
  attr(max.dpsl,'location') <- 'at inflexion points at lon/lat lines through storm center'
  attr(max.speed,'units') <- 'm/s'
  attr(max.speed,'location') <- 'at inflexion points at lon/lat lines through storm center'
  attr(radius,'units') <- 'km'

  results <- list(lon=lon[1:ii,],lat=lat[1:ii,],tim=tim[1:ii],psl=psl[1:ii,],
                  yy=yy[1:ii],mm=mm[1:ii],dd=dd[1:ii],i=1:ii,label=label,
                  max.dpsl=max.dpsl[1:ii,],max.speed=max.speed[1:ii,],
                  radius=radius[1:ii,],rad.max.dpsl=rad,dx=dx,dy=dy,
                  version="cyclones v1.1-6 (after Sep. 8, 2005)",
                  scale.factor.x = resx$scale.factor,scale.factor.y = resy$scale.factor)
  if (mod(ii,30)==0) save(file=back.up.file,results)
  if (mod(ii,10)==0) save(file=fname,results)
  print(paste("ii=",ii," it=",it," yy=",yy[ii]," mm=",mm[ii],
              " dd=",dd[ii]," N.lows=",sum(i.sim)," PSL min=",round(psl[ii,1]),
              " (hPa) dPSL max=",round(max.dpsl[ii,1],2)," max.speed=",
              round(max.speed[ii,1],1),"m/s radius=",round(radius[ii,1],1),"km",
              " lon=",round(lon[ii,1],2)," lat=",round(lat[ii,1],2)),sep="")




#------------------------------------------------------------------------------
# Plotting...

#  print(c((lplot) & (plot.now), lplot,plot.now))

    if ((lplot) & (plot.now)) {

    ii.c <- as.character(ii)
    while (nchar(ii.c) < 5) ii.c <- paste("0",ii.c,sep="")
    
#    newFig()
    if (EPS) postscript(file = paste(graph.dir,"cci",ii.c,".eps",sep=""),onefile=FALSE,horizontal=FALSE) else
            bitmap(file = paste(graph.dir,"cci",ii.c,".png",sep=""),type="png256")
#    bitmap(file = "cyclones.jpg",type="jpeg",width=15, height=15, res=250)
    par(las=1)
    image(lonx,latx,slpmap,zlim=c(-20,20),
          col = my.col,
          main="SLP anomalies",xlab="Longitude (deg E)",ylab="Latitude (deg N)",
          sub=paste(slp$yy[it],"-",slp$mm[it],"-",slp$dd[it]," #",ii,sep=""))
    grid()
    addland()
    contour(slp$lon,slp$lat,t(slp$dat[it,,])-mslp$map,levels=seq(-20,20,by=5),lwd=2,add=TRUE)
    contour(lonx,latx,resy$Z.fit,levels=seq(-20,20,by=5),add=TRUE, col="grey40",lwd=2)
    contour(lonx,latx,resx$Z.fit,levels=seq(-20,20,by=5),add=TRUE, col="grey20",lwd=2,lty=2)
    if (sum(lows)>0) {
      points(lon[ii,],lat[ii,],pch=20,col="black",cex=1.75)
      points(lon[ii,],lat[ii,],pch=20,col="white",cex=1.25)
      points(lon[ii,],lat[ii,],pch=".",col="black",cex=1.25)
    }
    dev.off()
  }


  if ((lplot) & (plot.now)) {
    if (EPS) postscript(file = paste(graph.dir,"cyclones_y_",vname,ii.c,".eps",sep=""),onefile=FALSE,horizontal=FALSE) else
            bitmap(file = paste(graph.dir,"cyclones_y_",vname,ii.c,".png",sep=""),type="png256")
    par(las=1)
    plot(latx,slpmap[ilon1,]+mslpmap[ilon1,],ylim=c(970,1050),type="n",
         main=paste("SLP profile at ",lon[ii,1],"E",sep=""),
         xlab="Latitude (deg N)",ylab="SLP (hPa)",
         sub=paste(slp$yy[it],"-",slp$mm[it],"-",slp$dd[it]," #",ii,sep=""))
    grid()
    polygon(c(rep(p.infly1[1],2),rep(p.infly1[2],2),p.infly1[1]),c(950,1050,1050,950,950),
            col="grey90",border="grey80",lwd=3)
    lines(range(latx),rep(mean(mean(slpmap[ilon1,]+mslpmap[ilon1,])),2),col="grey70")
    lines(rep(latx[ilat1],2),c(950,1050),col="grey70")
    lines(latx,slpmap[ilon1,]+mslpmap[ilon1,],lwd=5,col="grey50")
    points(latx[ilat1],mean(mean(slpmap[ilon1,]+mslpmap[ilon1,])),col="grey40",pch=20,cex=1.5)
    points(latx[ilat1],mean(mean(slpmap[ilon1,]+mslpmap[ilon1,])),pch=21,cex=1.4)
    lines(latx,resx$Z.fit[ilon1,]+mslpmap[ilon1,],col="black",lwd=2,lty=2)
    lines(latx,dslpdy[ilon1,]/sd(c(dslpdy[ilon1,]))*sd(c(slpmap[ilon1,]))
          + mean(slpmap[ilon1,]+mslpmap[ilon1,]),lty=2)
    lines(latx,dslpdy2[ilon1,]/sd(c(dslpdy2[ilon1,]))*sd(c(slpmap[ilon1,]))
          + mean(slpmap[ilon1,]+mslpmap[ilon1,]),lty=3,col="grey20")
    legend(min(latx),1050,c("Original         ","Fitted       ","deepest minima        ",
           "dP(y)/dy      ","d^2P(y)/dy^2      "),
           lty=c(1,2,0,2,3),pch=c(26,26,20,26,26),col=c("grey50","black","grey40","black","grey20"),
           lwd=c(5,2,0,1,1),bg="grey95")

#    lines(x.test,y.test/sd(y.test)*sd(slpmap[ilon1,])+mean(slpmap[ilon1,]+mslpmap[ilon1,]),col="red")
    dev.off()

# X-profile:
    if (EPS) postscript(file = paste(graph.dir,"cyclones_x_",vname,ii.c,".eps",sep=""),
                        onefile=FALSE,horizontal=FALSE) else
            bitmap(file = paste(graph.dir,"cyclones_x_",vname,ii.c,".png",sep=""),type="png256")
    #bitmap(file = "cyclones_x.jpg",type="jpeg",width=15, height=15, res=250)     
    par(las=1)
    plot(lonx,slpmap[,ilat1]+mslpmap[,ilat1],ylim=c(950,1050),
         main=paste("SLP profile at ",latx[ilat1],"N",sep=""),
         xlab="Longitude (deg E)",ylab="SLP (hPa)",
         sub=paste(slp$yy[it],"-",slp$mm[it],"-",slp$dd[it]," #",ii,sep=""))
    grid()
    points(lonx[ilon1],mean(mean(slpmap[,ilat1]+mslpmap[,ilat1])),col="grey60",pch=20)
    points(lonx[ilon1],mean(mean(slpmap[,ilat1]+mslpmap[,ilat1])),,pch=21)
    lines(lonx,resy$Z.fit[,ilat1]+mslpmap[,ilat1],col="grey",lwd=2)
    lines(lonx,resx$Z.fit[,ilat1]+mslpmap[,ilat1],col="darkgrey",lwd=1,lty=2)
    lines(range(lonx),rep(mean(mean(slpmap[,ilat1]+mslpmap[,ilat1])),2),lty=2,col="blue")
    while (dev.cur() > 1) dev.off()

    plot.now <- FALSE
    if (!file.exists(".CCI.run")) stop("Process halted")
  }
 }
 i.max <- it
} else {
  if (verbose) print("### sum(slp$tim > max(tim,na.rm=TRUE)) <= 0 ### !!!")
  print(paste("max(tim,na.rm=TRUE)=",max(tim,na.rm=TRUE),
              " sum(slp$tim <= max(tim,na.rm=TRUE))=",sum(slp$tim <= max(tim,na.rm=TRUE))))
  if (verbose) print(summary(slp$tim))
  if (verbose) print(summary(tim))
} # end of if (sum(slp$tim > max(tim))>0)
} # end of is-loop

file.remove(".CCI.run")

results <- list(lon=lon,lat=lat,tim=tim,psl=psl,
                  yy=yy,mm=mm,dd=dd,i=1:ii,label=label,
                  max.dpsl=max.dpsl,max.speed=max.speed,
                  radius=radius,rad.max.dpsl=rad,dx=dx,dy=dy,
                  version="cyclones v1.1-6 (after Sep. 8, 2005)",
                  scale.factor.x = resx$scale.factor,scale.factor.y = resy$scale.factor)
save(file=fname,results)
file.remove(back.up.file)

} else {
  print("Skipping this netCDF file")
  print(paste('ii=',ii,'  length(times)=',length(times),'  max.tim=',max.tim,
              '  sum(times.physical > tim[ii])=',sum(times.physical > tim[ii]),'tim[ii]=',tim[ii]))
  
  print(i.here)
} # end of if - is.finite(times[is+1])

}

stopCCI <- function() {
  file.remove(".CCI.run")
  print("Please let the present cycle finish")
}

