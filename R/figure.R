
# to load the funtion into R use : source("/home/abdelkaderm/Documents/R_scripts/plot.R")

frame.metno <- function(z=z,select = NULL, col="gray50", col.select = NULL,add.climatology=TRUE,cex.lab=1,...) {

  ## Get data in x
  x  <- as.matrix(coredata(z)) 
  year <- as.numeric(year(index(z)))
  nts <- dim(x)[2]
  
  ## Set Y axis limits
  ylim1 <- round(range(x,na.rm=TRUE)[1]) # round(min(x,na.rm=TRUE))
  ylim2 <- round(range(x,na.rm=TRUE)[2]) # round(max(x,na.rm=TRUE)) 
  
  ## default plot 
  plot(year,x[,1],type="l",col="white",frame.plot=FALSE,axes=FALSE,xlab="YEARS",ylab ="TAS [Deg C]" , ylim=c(ylim1,ylim2),cex.lab=cex.lab,...) #,label.pos=
  ## add line
  ##col <- ifelse(x[,1]>=0,"red","blue")
  lines(year,x[,1],col=col,lw=1.5)
  points(year,x[,1],col=col,type="p",pch=21,bg="grey",cex=0.75) 
  abline(h=0,lty=2)
  ## add title
  title(main=paste("REGION : ", toupper(attr(z,"country")),"(Coord.)"),line=3,cex.main=cex.lab)
  title(main=attr(z,"title"),line=2.2,cex.main=cex.lab)

  ## add margin text
  mtext(paste(toupper("esd package - fig-0.0 - copyright metno 2014"),"(www.met.no)",sep=" "),side=1,line=4,cex=cex.lab)
  ## mtext("Anomaly values relative to 1986-2005",side=2,line=2,cex=cex.lab,las=3)

  ## display axis
  aty <- seq(year[1]-year[1]%%10,year[length(year)]-year[1]%%10+10,10)
  axis(1,at=aty,col ="black",...) # lwd.ticks=2
  axis(2,at=seq(ylim1,ylim2,0.5),col ="black",...) # ,lwd.ticks=2
  axis(3,at=aty,col ="black",...)
  axis(4,at=seq(ylim1,ylim2,0.5),col ="black",...)
  
  ##axis(1,at=seq(year[1],year[length(year)],10),col ="black",col.axis=col.axis[2],padj=0) # lwd.ticks=2
  ##axis(2,at=seq(ylim1,ylim2,0.5),col ="black",col.axis=col.axis[2],padj=0) # ,lwd.ticks=2
  ##axis(3,at=seq(2002,2100,10),col ="black",col.axis="black",cex=0.4) # ,lwd.ticks=2
  ##axis(4,at=seq(1,3,0.5),col ="black",col.axis="black") #,lwd.ticks=2
  ##text(x=2048,y=3,labels="Future Period",col = "black",cex = 1)
  ##text(x=1900,y=0,labels="Hindcast",col = col.axis[1],cex = 1)
  ##text(x=1980,y=-4,labels="Reanalysis Period",col = "Black",cex = 1)
  
}

# Function to generate figures of type fig0 # we should give names later
fig0.zoo <- function(z=z,select = c(1,2), col="gray90", col.select = NULL,add.clim=TRUE,text=FALSE) {

# Get data in x
x  <- as.matrix(coredata(z)) 
year <- as.numeric(year(index(z)))
nts <- dim(x)[2]

# Set Y axis limits
ylim1 <- round(range(x,na.rm=TRUE)[1]) # round(min(x,na.rm=TRUE))
ylim2 <- round(range(x,na.rm=TRUE)[2]) # round(max(x,na.rm=TRUE)) 

# default plot 
plot(year,x[,1],type="l",col="white",frame.plot=FALSE,axes=FALSE,xlab="YEARS",ylab ="TAS [Deg C]" , ylim=c(ylim1,ylim2)) #,label.pos=

if (text & !is.null(attr(z,"country")))
    title(main=paste("REGION : ", toupper(attr(z,"country")),"(Coord.)"),line=3,cex.main=.9)

## add title
txt <- attr(z,"title")
if (text & !is.null(txt)) {
    txtsplit <- unlist(strsplit(txt,split=" "))
    itxt <- grep("for",txtsplit)
    title(main=paste(txtsplit[(itxt+1):length(txtsplit)],collapse=" "),line=2.2,cex.main=0.7)
}
##title(main=attr(z,"title"),line=2.2,cex.main=0.7)

# add margin text
if (text) mtext(paste(toupper("esd package - fig-0.0 - metno 2013"),"(www.met.no)",sep=" "),side=1,line=4,cex=0.6)
## if (attr(anomaly) mtext("Anomaly values relative to 1986-2005",side=2,line=2,cex=0.8,las=3)

# display axis
aty <- seq(year[1]-year[1]%%10,year[length(year)]-year[1]%%10+10,10)
axis(1,at=aty,col ="black",cex=.7) # lwd.ticks=2
axis(2,at=seq(ylim1,ylim2,1),col ="black",cex=.7) # ,lwd.ticks=2
axis(3,at=aty,col ="black",cex=.7)
axis(4,at=seq(ylim1,ylim2,1),col ="black",cex=.7)


# Get climatology if any and add shaded area
## clim <- attr(z,"climatology")
clim <- attr(z,"baseperiod")
if (!is.null(clim)) {
    clim1 <- as.numeric(unlist(strsplit(clim,split="-"))[1])
    clim2 <- as.numeric(unlist(strsplit(clim,split="-"))[2])  
    rect(xleft=clim1,xright=clim2,ybottom=ylim1,ytop=ylim2,col="lightyellow",border=NA)
    text(x=(clim1+clim2)/2-10,y=ylim2-3,labels=paste("Base Period",sep=" "),cex=.6,srt=0,adj=0)
}

for (k in 1:nts) { 
    lines(year,x[,k],col=col,lwd=1)
} 

if (!is.null(attr(z,"aspect"))) if (attr(z,"aspect") == "anomaly") abline(0,0,lty=2)
## abline(2,0,lty=2,col="orange")

# Get colomn names
if (!is.null(attr(z,"model_id")))
    cnames <- attr(z,"model_id")
if (!is.null(colnames(z)))
    cnames <- colnames(z)

# higlight selected lines
if (!is.null(select)) {
   if (is.character(select)) select <- grep(tolower(paste(collapse(select,"|"))),tolower(cnames))
   if (is.null(col.select)) { 
      col.select <- rev(rainbow(length(select)))
   } else if (length(select)!=length(col.select)) { 
     stop("Length of select and col.select vectors must be equal")
   }
for (j in 1:length(select)) lines(year,x[,select[j]],col=col.select[j],lwd=1)
# Add legend 
legend(x="topleft",legend=c(cnames[select],"All models"),col=c(col.select,col),lty=1,bty="n",cex=0.8)
}

} # End of function plot.fig0.0

# Function to generate figures of type fig1
fig1.zoo <- function(z=list(z1=z1,z2=z2),select = NULL, col.select=NULL,tsline = FALSE, col="grey",enveloppe=TRUE,col.enveloppe=NULL,add.clim=FALSE,add.2C=FALSE,error=-1) {

if (is.list(z)) {
   x  <- as.matrix(coredata(z[[1]])) 
   nx <- length(names(z))
   year <- as.numeric(year(index(z[[1]])))
} else if (is.zoo(z)) { 
   x  <- as.matrix(coredata(z)) 
   nx <- 1
   year <- as.integer(year(year(index(z))))
   ny <- length(year)
}
## browser()
# Set Y axis limits
ylim1 <- round(range(x,na.rm=TRUE)[1]) # round(min(x,na.rm=TRUE))
ylim2 <- round(range(x,na.rm=TRUE)[2]) # round(max(x,na.rm=TRUE)) 

# Create the plot and define axis variables
plot(year,x[,1],type="l",col="white",frame.plot=FALSE,axes=FALSE,xlab="YEARS",ylab ="TAS [Deg C]",ylim=c(ylim1,ylim2))

# add region name
title(main=paste("REGION : ", toupper(attr(x,"country")),"(Coord.)"),line=3,cex.main=.9)
# add experiment
if (nx==1) {
   txt <- attr(x,"title")
   if (!is.null(txt)) {
      txtsplit <- unlist(strsplit(txt,split=" "))
      itxt <- grep("for",txtsplit)
      title(main=paste(txtsplit[(itxt+1):length(txtsplit)],collapse=" "),line=2.2,cex.main=0.7)
   }
#title(main=attr(z,"title"),line=2.2,cex.main=0.7)
}

# add margin text
mtext(paste(toupper("esd package - fig-0.0 - metno 2013"),"(www.met.no)",sep=" "),side=1,line=4,cex=0.6)
## mtext("Anomaly values relative to 1986-2005",side=2,line=2,cex=0.8,las=3) 

# display axis
aty <- seq(year[1]-year[1]%%10,year[length(year)]-year[1]%%10+10,10)
axis(1,at=aty,col ="black",cex=.7) # lwd.ticks=2
axis(2,at=seq(ylim1,ylim2,0.5),col ="black",cex=.7) # ,lwd.ticks=2
axis(3,at=aty,col ="black",cex=.7)
axis(4,at=seq(ylim1,ylim2,0.5),col ="black",cex=.7)

col.axis <- c("grey50","black","grey50")

#Get climatology if any
if (add.clim) {
  clim <- attr(x,"climatology")
  if (!is.null(clim) & add.climatology) {
    clim1 <- as.numeric(unlist(strsplit(clim,split="-"))[1])
    clim2 <- as.numeric(unlist(strsplit(clim,split="-"))[2])
    rect(xleft=clim1,xright=clim2,ybottom=ylim1,ytop=ylim2,col="lightyellow",border=NA)
    text(x=c(round(clim1+clim2)/2),y=1,labels=paste("Climatology",sep=" "),cex=.9,srt=90,adj=0)
  }  
  abline(0,0,lty=2)
}

if (add.2C) {
abline(2,0,lty=2,col="red")
text(x=c(round(clim1+clim2)/2),y=2.1,labels = "2°C",col="red")
}
#   axis(1,at=seq(year[1],clim1,10),col ="black",col.axis=col.axis[2]) # lwd.ticks=2
#   axis(2,at=seq(-1.5,0.5,0.5),col ="black",col.axis=col.axis[2]) # ,lwd.ticks=2
#   axis(3,at=seq(clim2,year[length(year)],10),col ="black",col.axis="black",cex=0.4) # ,lwd.ticks=2
#   axis(4,at=seq(1,3,0.5),col ="black",col.axis="black") #,lwd.ticks=2
#   text(x=2048,y=3,labels="Future Period",col = "black",cex = 1)
#   text(x=1900,y=0,labels="Hindcast",col = col.axis[1],cex = 1)
#   text(x=1980,y=-4,labels="Reanalysis Period",col = "Black",cex = 1)
#}

for (k in 1:nx) {
   if (is.list(z)) {
      x  <- as.matrix(coredata(z[[k]])) 
      year <- as.numeric(year(index(z[[k]]))   )
   } else if (is.zoo(z)) {
      x  <- as.matrix(coredata(z))
      year <- as.numeric(year(index(z)))
   } 
   # Compute the number of years
   ny <- length(year)
   if (enveloppe) { # produce enveloppe plot
      if (!is.null(col.enveloppe)) { 
         i1 <- round(col2rgb(col.enveloppe)[1]/255)
         i2 <- round(col2rgb(col.enveloppe)[2]/255)
         i3 <- round(col2rgb(col.enveloppe)[3]/255)
      } else {# Randomly select an rgb colors 
         i1 <- round(runif(1,0,1))
         i2 <- round(runif(1,0,1))
         i3 <- round(runif(1,0,1))
      } 
      probs <- seq(0, 1, 0.25)
      q <- matrix(NA,length(year),length(probs))
      ## browser()
      for (i in 1:length(year)) {
          q1 <-  quantile(x[i,], probs = probs, na.rm = TRUE, names = TRUE, type = 7)
          q[i,] <- as.numeric(rownames(table(q1)))
      }
     
      if (error>0) {
         std <- apply(t(q),2,"sd")
         z3 <- c(q[1:ny,3]-error*std[1:ny],rev(q[1:ny,3]+error*std[1:ny]))
      } else z3 <- c(q[1:ny,1],rev(q[1:ny,length(probs)]))
      t3 <- c(year[1:ny],rev(year[1:ny]))
      polygon(t3,z3, col=rgb(i1,i2,i3,.2),lty=2,border=NA)
      lines(year,q[,3],type="l",col=rgb(i1,i2,i3,1),lwd=2)
   } else # end enveloppe plot
       for (i in 1:dim(x)[2]) lines(year,x[,i],type="l",col=col,lwd=1)

       if (nx>1) {
      if (k == 1) {i11 <- i1 ; i21 <- i2 ; i31 <- i3}
      if (k == 2) {i12 <- i1 ; i22 <- i2 ; i32 <- i3}
      if (k == nx) 
         legend(x="topleft",legend=c("CMIP3","CMIP5"),col=c(rgb(i11,i21,i31,1),rgb(i12,i22,i32,1)),lty=1,bty="n",cex=0.8)
   }
}
if (!is.null(attr(z,"model_id")))
    cnames <- attr(z,"model_id")
if (!is.null(colnames(z)))
    cnames <- colnames(z)

if (!is.null(select)) {
   if (is.character(select)) 
      iselect <- grep(tolower(paste(select,collapse="|")),tolower(cnames))
   else if (is.numeric(select)) 
      iselect <- select
   if (length(iselect)<1) 
      print(paste("Selection",select," has not been found in attributes : ignored",sep=" "))  
   else {
      if (is.null(col.select)) { 
         col.select <- rev(rainbow(length(iselect)))
      } else if (length(select)!=length(col.select)) { 
         stop("Length of select and col.select vectors must be equal")
      }
      for (j in 1:length(iselect)) lines(year,x[,iselect[j]],col=col.select[j],lwd=1)
      # Add legend 
      if (enveloppe)
         legend(x="topleft",legend=c(cnames[iselect],"Median","Enveloppe"),col=c(col.select,rgb(i1,i2,i3,1),rgb(i1,i2,i3,0.2)),lty=1,bty="n",cex=0.8)
      else
         legend(x="topleft",legend=c(cnames[iselect]),col=c(col.select),lty=1,bty="n",cex=0.8) 
   }
}

} # End of fig_1

