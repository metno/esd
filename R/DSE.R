## A. Mezghani - DownScale the CMIP5 ensemble seasonal mean and standard deviation for temperature for several stations
##

DSE <- function(stid="18700",cntr="NORWAY",src="METNOD",param="t2m",lon=c(-10,10),lat=c(-10,10),FUN="mean",path="CMIP5.monthly",rcp="rcp45",biascorrect=TRUE,predictor="ERA40_t2m_mon.nc",email=NULL,save=TRUE,force=FALSE,verbose=FALSE,out.dir="dse.test6",update=FALSE,plot=TRUE,show.map=TRUE,select=1:3,it=NULL,...) {
    source("/home/abdelkaderm/R_scripts/DSensemble.month.R")
    ## email="abdelkaderm@met.no"
    if (!is.null(email))
        if (!library(sendmailR,logical.return=TRUE))
            install.packages("sendmailR")
        else library(sendmailR)
    
    ## Define output directory and create it if does not exist
    out.path <- file.path(path,rcp,out.dir)
    if (!file.exists(out.path)) dir.create(out.path)
    
    ## Search for temporary files ~
    tmpfile <- list.files(pattern="~",out.path,full.names=TRUE)
    file.remove(tmpfile)
    
    ## grep for files
    lfiles <- list.files(pattern="DSE",out.path,full.names=TRUE)
    files <- list.files(pattern="DSE",out.path)
    
    ## check if inventory file exists and is updated
    ## check if log file exists
    logfile <- file.path(out.path,"dse-log.txt") 
    if (!file.exists(logfile) | (file.info(logfile)$size==0)) {
        file.create(logfile)
        if (verbose) print(paste("Creating log file --> ",logfile))
    }
    ## browser()
    ## Checking inv file
    invfile <- file.path(out.path,"dse-inv.txt")
    if ((!file.exists(invfile) | (file.info(invfile)$size==0))) {
        if (verbose) print(paste("Creating/Generating inventory file", invfile))
        meta <- data.frame(NULL)
    }
    else {
        meta <- data.frame(read.csv(invfile,header=TRUE,sep=";"))
        ## open connection to invfile and append data
        ## finv <- file(file.path(out.path,invfile),"at")  
    }
    ## browser()
    
    ## stopifnot(!missing(ss),inherits(ss,"stationmeta"),file.exists(file.path(path,rcp)))
    ## open log and inv file 
    flog <- file(logfile,"at")  
    ## writeLines(paste("#",date(),sep=""),con=flog)
    finv <- file(invfile,"at")
    ## writeLines(paste("#",date(),sep=""),con=finv)
    
    ## browser()
    ## select stations
    ss <- select.station(stid=stid,cntr=cntr,src=src,param=param,it=it,...)
    ## browser()
    if (show.map) {
        col <- rev(terrain.colors(n=10))
        map(ss,cex=0.8,bg="grey",col="grey50",xlim=c(-50,180),ylim=c(0,90),height=10,width=10)
        col.bar(col=col,breaks=seq(0,1,0.1),type="p",v=3,h=3)
        ## rect(xleft=lon[1],xright=lon[2],ybottom=lat[1],ytop=lat[2],border="red",lwd=1)
        ## devid <- dev.cur()
    }
    ## browser()
    if (!is.null(ss)) {
        ## if (is.null(dim(y)) & !is.null(length(y))) d <- 1 else d <- dim(y)[2] 
        d <- length(ss$location)
        k <- 0 ## counter
        for (i in 1:d) { 
            k <- k + 1
            ## open connections to the log and inv files
            flog <- file(logfile,"at")
            finv <- file(invfile,"at")
            ## retrieve station meta data
            loc <- ss$location[i]
            param <- toupper(esd2ele(ss$element[i]))
            src <-  ss$source[i]
            stid <- ss$station_id[i]
            text <- paste(param,loc,stid,src,rcp)
            
            ## retrieve data for selected station
            if (verbose) print(paste("Retrieving",text))      
            eval(parse(text=paste("y <- station.",tolower(src),"(stid=stid,param=tolower(param),...)",sep="")))
            lon.1 <- round( attr(y,'longitude') + lon )
            lat.1 <- round( attr(y,'latitude') + lat )
            ##if (show.map)
            ##  rect(xleft=lon.1[1],xright=lon.1[2],ybottom=lat.1[1],ytop=lat.1[2],border="red",lwd=1)
            
            ## Keep only first stations if many variates
            if (length(attr(y,"station_id")>1)) y <- subset(y,is=1)
            attr(y,"aspect") <- "original"
            
            ## browser()
            if (verbose) str(y)
            
            if (!is.null(y)) {        
                ## look for existing files
                ## generate file name
                rea <-strsplit(predictor,split="_")[[1]][1]
                loc <- gsub(" ","",loc(y),fixed=TRUE)
                loc <- gsub("/","-",loc,fixed=TRUE) ## replace "/" by "-" if any
                loc <- gsub("(","-",loc,fixed=TRUE) ## replace "(" by "-" if any
                loc <- gsub(")","-",loc,fixed=TRUE) ## replace ")" by "-" if any
                type <- "season"
                dom <- list(lon=round((lon(y)+lon)),lat=round(lat(y)+lat))
                reaexp <- paste(rea,rcp,sep="-")
                srcvar <- toupper(paste(attr(y,"source"),substr(type,1,1),FUN,attr(y,"variable"),sep="-"))
                txtdom <- paste("LON",dom$lon[1],"E",dom$lon[2],"E&LAT",dom$lat[1],"N",dom$lat[2],"N",sep="")
                filename = paste(paste("DSE",loc,srcvar,toupper(reaexp),txtdom,sep="_"),".rda",sep="")
                ## browser()
                if (!file.exists(file.path(out.path,filename)) | force) {
                    ## Do the downscaling
                    class(ss) <- "data.frame"
                    print(paste("DOWNSCALING MONTHLY",toupper(FUN),text))
                    ## if (plot) par(new=TRUE)
                    if (tolower(param)=='t2m') {
                        z <- try(DSensemble.t2m.month(y,FUN=FUN,rcp=rcp,predictor=predictor,biascorrect=biascorrect,path=path,lon=lon,lat=lat,select=select,plot=FALSE,...))
                        attr(z,"method") <- FUN
                    } else if (tolower(param)=='precip') {
                        z <- try(DSensemble.precip(y,FUN=FUN,rcp=rcp,predictor=predictor,biascorrect=biascorrect,path=path,lon=lon,lat=lat,select=select,plot=FALSE,...))
                        attr(z,"method") <- FUN
                    } else z <- NULL   
                }
                else {
                    print(paste("Seems like the file",files[i]," exists and that the station has already been downscaled. Please set the argument force to TRUE to overwrite it")) 
                    load(file.path(out.path,filename),envir=environment()) 
                }
                ## close current figure
                ## if (plot) dev.off()
                ## browser()
                ## write errors in flog
                if (inherits(z,"try-error")) {    
                    txt <- as.character(i)
                    txt <- paste(txt,tolower(src(y)),type,FUN,varid(y),stid(y),loc(y),cntr(y),lon(y),lat(y),alt(y),esd2ele(varid(y)),start(y),end(y),sep=";")
                    txt <- paste(text,start(y),end(y),paste(summary(y)[,2],collapse=";"))
                    writeLines(paste(text,start(y),end(y),paste(summary(y)[,2],collapse=";"),length(y),sep=";"),con=flog)
                    writeLines(z[[1]],con=flog)
                    if (show.map) {
                        dev.prev()
                        points(ss$longitude[i],ss$latitude[i],pch=4,col="red",cex=0.8)
                    }
                }
                else if (!is.null(z)) {
                    ## write meta for downscaled station
                    ## writeLines(paste(text,start(y),end(y),paste(summary(y)[,2],collapse=" "),length(y),sep=" "),con=finv)
                    if (save) save(file=file.path(out.path,filename),z)
                    ## write meta data into inv file
                    y <- attr(z,"station")
                    type <- class(y)[2]
                    print(attr(z,"scorestats"))
                    qual <- mean(attr(z,"scorestats"))
                    r2 <- mean(attr(z,"scorestats")[,1])
                    print(r2)
                    method <- FUN
                    txt <- as.character(i)
                    txt <- paste(txt,tolower(src(y)),type,method,varid(y),stid(y),loc(y),cntr(y),lon(y),lat(y),alt(y),esd2ele(varid(y)),start(y),end(y),sprintf("%3.2f", qual),filename,sep=";")
                    writeLines(txt,con=finv)
                    if (!is.na(r2)) {
		       if (show.map) {
                        ## dev.prev()
                        ## browser()
                        k <- 0 ; icn <- seq(0,1,0.1)
                        for (c in 1:10) {
                            id <- (icn[c] <= r2) & (r2 < icn[c+1])
                            if (id) k <- c
                        }
                        points(ss$longitude[i],ss$latitude[i],pch=21,bg=col[k],col="black",cex=1) ## "darkgreen"
                        ## dev.new()
                       }
                    }
		}
                ##} ##else print(paste("Seems like the file",files[i]," exists and that the station has already been downscaled. Please set the argument force to TRUE to overwrite it"))   
                ## close connection with log and inv file
                close(flog)
                close(finv)
            }
            con <- file(logfile,"rt")
            msg <- readLines(con)
            close(con)
        } 
    }
    ## close connection to the log file
    ## capture.output(
    ## browser()
    ## close log and inv files
    ##close(flog)
    ##close(finv)
    
    ## send Notification by email
    if (!is.null(email))
        sendmail(from=paste("<",email,">",sep=''),to=paste("<",email,">",sep=''),subject=paste("DOWNSCALING COMPLETED - THIS IS AN AUTOMATED EMAIL FROM DSE Function",file.path(path,"dse")),msg=msg)
    else return(NULL)
    invisible(z)
    
    ## update the meta file
    if (update) meta.dse(out.path)
}

meta.dse <- function(path="CMIP5.monthly/rcp45/dse",plot=FALSE,verbose=FALSE) {
    ## browser()
    invfile <- file.path(path,"dse-inv.txt")
    options(stringsAsFactors = FALSE)
    if (!file.exists(invfile) | file.info(invfile)$size==0) {
        if (verbose) print(paste("Creating/Generating inventory file", invfile))
        meta <- data.frame(NULL)
    }
    else {
        meta <- data.frame(read.csv(invfile,header=TRUE,sep=";"))
        ## open connection to invfile and append data
        ## finv <- file(file.path(path,invfile),"at")  
    }
    ## browser()
    ## grep for files
    lfiles <- list.files(pattern="DSE",path,full.names=TRUE)
    files <- list.files(pattern="DSE",path)
    ## browser()
    ## Generates inventory file if does not exist or has not been updated
    if ((length(lfiles)>0) & ((dim(meta)[1]==0) | (dim(meta)[1]!=length(lfiles)))) {## create new meta files from dse folder
        ## creates new file
        if (verbose) print("CREATE/UPDATE META FILE")
        file.create(invfile)
        ## open connection to invfile
        finv <- file(invfile,"at")
        ## Insert header to invfile
        if (file.info(invfile)$size==0)
            writeLines(paste(c("id","source","type","fun","variable","station_id","location","country","longitude","latitude","altitude","element","start","end","score","filename"),collapse=";"),con=finv)
        ## loap on dse list of files 
        for (i in 1:length(lfiles)) {
            load(lfiles[i])
            y <- attr(z,"station")
            type <- class(y)[2]
            qual <- mean(attr(z,"scorestats"))
            method <- attr(z,"method") ## deparse(substitute(FUN))
            txt <- as.character(i)
            txt <- paste(txt,tolower(src(y)),type,method,varid(y),stid(y),loc(y),cntr(y),lon(y),lat(y),alt(y),esd2ele(varid(y)),start(y),end(y),sprintf("%3.2f", qual),files[i],sep=";")
            writeLines(txt,con=finv)
            if (verbose) pb <- txtProgressBar(style=3)
            if (verbose) setTxtProgressBar(pb,i/length(lfiles))
        } 
        if (verbose) close(pb)
        close(finv)
        if (verbose) print("DONE !")
    }
    class(meta) <- "stationmeta"
    meta$call <- match.call()

    if (is.null(meta))
        if (plot) map(meta)
    invisible(meta)
}

as.residual.ds <- function (x,detrend=TRUE) {
    yo <- attr(x, "original_data")
    yf <- attr(x, "fitted_values")
    if (detrend)
        y <- trend(yf,result="residual") - trend(yo,result="residual")
    else
        y <- yf - yo
    y <- attrcp(x, y)
    attr(y, "aspect") <- "residual"
    attr(y, "history") <- history.stamp(x)
    class(y) <- class(attr(x, "calibration_data"))
    invisible(y)
}


col.bar <- function(x,horiz=TRUE,v=1,h=1,col=col,cex=0.7,type="r",...) {
    
    xleft <- par()$usr[1] 
    xright <- par()$usr[2]
    ybottom <- par()$usr[4] - 1 - h
    ytop <-  par()$usr[4] - 1 
    
    steps <-   seq(0, (xright -xleft - v * (length(col))) , (xright - xleft - v * (length(col)))/(length(col))) # 
    nsteps <- length(steps) - 1 
    icn <- seq(0,1,1/nsteps) ; print(icn)
    k <- 0
    for (i in 1 :nsteps) {  
        if (!is.null(v)) 
            if (i == 1) k <- v/2 else k <- k + v  
        if (type == "r") { ## "r" for rectangle
            rect(xleft= k  + xleft + steps[i] ,xright= k + xleft + steps[i+1],ybottom=ybottom,ytop=ytop,col=col[i])
            text(x = k + xleft + steps[i],  y = ybottom - 1, labels=sprintf("%.1f",icn[i]),cex=cex)
            text(x = k + xleft + steps[i+1],y = ybottom - 1, labels=sprintf("%.1f",icn[i+1]),cex=cex)
            ## text(x = k + xleft + steps[i], y = ybottom - 1,labels=sprintf("%.1f",icn[i]),cex=cex)
        }
        else if (type == "p") { ## "p" points
            points(x= k + xleft + (steps[i]+ steps[i+1])/2, y=(ybottom + ytop)/2,pch=21, bg=col[i],cex=v)
            text(x = k + xleft + steps[i],  y = ybottom - 1, labels=sprintf("%.1f",icn[i]),cex=cex)
            text(x = k + xleft + steps[i+1],y = ybottom - 1, labels=sprintf("%.1f",icn[i+1]),cex=cex)
        }
    }
}
