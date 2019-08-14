#' Show summary of objects
#' 
#' Produce a summary table 
#'
#' @seealso summary.station summary.ncdf4
#'
#' @param x an object of type 'DSensemble'
#' @param years A set of years for which to produce summary statistics
#'
#' @return A matrix containing summary statistics
#' @examples
#' 
#' data("dse.Oslo")
#' summary(dse.Oslo)
#' 
#' @export summary.dsensemble
summary.dsensemble <- function(object,...,years=seq(1990,2090,by=20),verbose=FALSE) {
    if(verbose) print("summary.dsensemble")
    x <- object
    x0 <- subset(x,it=0)
    djf <- subset(x,it='djf')
    mam <- subset(x,it='mam')
    jja <- subset(x,it='jja')
    son <- subset(x,it='son')

    tab <- rep('',length(years) + 1)
    tab[1] <- paste(loc(x), '  Annual, DFJ, MAM, JJA, SON')
    i <- 1
    ##browser()
    for (yr in years) {
      i <- i + 1
      ##print(i,yr)
      tab[i] <- paste(years[i-1],':    ',
                 round(mean(coredata(subset(x0,it=years[i-1]))),2),
                 ' [',round(quantile(subset(x0,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(x0,it=years[i-1]),0.95),2),'],  ',
                 round(mean(coredata(subset(djf,it=years[i-1]))),2),
                 ' [',round(quantile(subset(djf,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(djf,it=years[i-1]),0.95),2),'],  ',
                 round(mean(coredata(subset(mam,it=years[i-1]))),2),
                 ' [',round(quantile(subset(mam,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(mam,it=years[i-1]),0.95),2),'],  ',
                 round(mean(coredata(subset(jja,it=years[i-1]))),2),
                 ' [',round(quantile(subset(jja,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(jja,it=years[i-1]),0.95),2),'],  ',
                 round(mean(coredata(subset(son,it=years[i-1]))),2),
                 ' [',round(quantile(subset(son,it=years[i-1]),0.05),2),', ',
                 round(quantile(subset(son,it=years[i-1]),0.95),2),']',sep='')
      ##print(tab[i])
    }
    tab
}

#' Show summary of objects
#' 
#' Produce a summary table 
#'
#' @seealso summary.dsensemble summary.ncdf4
#'
#' @param x an object of type 'station'
#' @param im The order of months in the table. Use \code{im=c(10:12,1:9)} for Oct-Sep.
#'
#' @return A matrix containing summary statistics
#' @examples
#' 
#' data("Oslo")
#' summary(Oslo)
#' 
#' @export summary.station
summary.station <- function(object,...,im=1:12,verbose=FALSE) {
  if(verbose) print("summary.station")
  x <- object
  tab <- matrix(rep(NA,12*7),12,7)
  for (i in 1:12) {
    y <- subset(x,it=month.abb[i])
    z <- as.numeric(summary(coredata(y)))
    attributes(z) <- NULL
    tab[i,1:length(z)] <-z 
  }
  attn <- attr(summary(coredata(x)),'names')
  if (length(attn)==6) colnames(tab) <- c(attn,"NA's") else colnames(tab) <- attn
  rownames(tab) <- month.abb
  tab[im,]  
}

summary.ds <- function(x) {
}

summary.eof <- function(x) {
}

summary.cca <- function(x) {
}


#' Summary of netcdf file
#'
#' @param ncfile filename of netcdf file
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @seealso retrieve check.ncdf4
#'
#' @export summary.ncdf4
summary.ncdf4 <- function(object, ..., verbose = TRUE) {
  if(verbose) print("summary.ncdf4")
  ncfile <- object
  if (is.character(ncfile)) {
   if (!file.exists(ncfile)) stop(paste("Sorry, the netcdf file '", ncfile,
      			                "' does not exist or the path has not been set correctly !",sep =""))
   ncid <- nc_open(ncfile)     
  } else if (class(ncfile) == "ncdf4") {
    ncid <- ncfile
  } else {
    stop("ncfile format should be a valid netcdf filename or a netcdf id of class 'ncdf4'")
  }

  model <- ncatt_get(ncid,0)
  if (verbose) print(ncid$filename)
  if (verbose) print(model$title)
  type <- c("Years","Seasons","Months","Days","Hours")
  type.abb <- substr(tolower(type),1,3)
  # Print summary from the netcdf
  for (i in 1:ncid$nvars) {
    vname = ncid$var[[i]]$name
    ndims = ncid$var[[i]]$ndims
    namedims  = names(ncid$dim)
    varstr = paste(vname, " (variable ", as.character(i),") has unit (",ncid$var[[i]]$units,
                   ") and dimension(s)", sep = "")
    for (j in 1:ndims) {
      varstr <- paste(varstr,ncid$var[[i]]$dim[[j]]$name,sep=" ")  
      varstr <- paste(varstr,"(",as.character(ncid$var[[i]]$dim[[j]]$len),")",sep="")
    }
    if (verbose) print(varstr, fill = TRUE)
    for (j in 1:ndims) {
      dimstr = paste(ncid$var[[i]]$dim[[j]]$name, " (dimension ", as.character(j),
                     ") Unit (", ncid$var[[i]]$dim[[j]]$units, ") First point (",
                     as.character(min(ncid$var[[i]]$dim[[j]]$vals)),") Last point (",
	             as.character(max(ncid$var[[i]]$dim[[j]]$vals)),") Step (",
 	             (ncid$var[[i]]$dim[[j]]$vals[2]-ncid$var[[i]]$dim[[j]]$vals[1]),")",sep = "")
      if (verbose) print(dimstr)
    }
    ifreq <- grep("freq",names(model))
    if (length(ifreq)>0) {  
      itype <- grep(tolower(eval(parse(text=paste("model$",names(model)[ifreq],sep="")))),tolower(type))
      if (length(itype)>0) {
        frequency <- model$frequency
        if (verbose) print(paste("Frequency has been found in data attributes: ",
                           model$frequency,"(Values)",sep=" ")) 
      }
    } else if (verbose) {
      print(paste("Warning: Frequency has not been found in data attributes !",sep=" ")) 
      frequency <- NA
    }
    if ((!is.null(grepl("tim",ncid$var[[i]]$dim[[ndims]]$name))) &
        (length(grep("calend",names(ncid$var[[i]]$dim[[ndims]])))>0)) {
      if (verbose) print(paste(ncid$var[[i]]$dim[[ndims]]$calendar,"calendar has been found",sep=" ")) 
    }
  }

  results <- list(title = model$title , name = vname , dims = namedims , frequency = frequency) 
  invisible(results)
} # End of function

