# Name		: test.retrieve.ncdf4
# Description	: Computes the spatial mean surface temperature based on the 92 GCM outputs from the CMIP5 experiment or any other experiment and save the results in an output file if specified.
# Author 	: Abdelkader Mezghani, METNO
# contact 	: abdelkaderm@met.no
# Last Update	: 21-03-2013
# require	: ncdf4,zoo,retrieve.ncdf4,global.tas

#source("esd/R/retrieve.R")
#source("esd/R/spatial.R")
require(ncdf4)
require(zoo)

test.retrieve.ncdf4 <- function(path="CMIP3.Monthly/20c3m-sresa1b/",file="Global_mean_tas_anomaly_cmip3_1986-2005.rda") {
ncfiles <- list.files(path=path,pattern="tas_",full.names=TRUE) # the dash is added to avoid selecting tasmin and tasmax if any in the same folder!

# Initialize
vmodel <- rep("-",length(ncfiles))
for (k in 1:length(ncfiles)) {
    ncid <- nc_open(ncfiles[k])
    print("-------------------------------------------")
    print(ncid$filename)
    gcm <- retrieve.ncdf4(ncfile=ncid,param="Auto", lon=NULL,lat=NULL, lev = NULL,time = NULL, ncdf.check=TRUE,miss2na = TRUE,silent = FALSE)
    model_id <- attr(gcm,"model_id")
    # Spatial averaging along lon and lat
    gcm2 <- spatial.avg.field(gcm)
    gcm.am <- aggregate(x=gcm2,by=format(time(gcm2),"%Y"),FUN=mean) # Annual mean for 1 grid cell
    year <- time(gcm.am)   
    agcm.ave <- gcm.am-mean(gcm.am[is.element(as.numeric(year),c(1986:2005))])
    if (k == 1) {      
       #model_id0 <- model_id
      X <- matrix(NA,length(year),length(ncfiles)) 
      plot(agcm.ave,col="steelblue" ,ylim=c(-5,5)) 
    } else lines(agcm.ave,col="steelblue")
    if (!is.null(model_id)) { 
       vmodel[k+1] <- model_id
    }
    X[,k] <- agcm.ave 
}
# output format
z <- zoo(x=X,order.by=year)
attr(z,"name") <- "tas"
attr(z,"country") <- "Global"
attr(z,"longname") <- "Spatial average of the 2m-Temperature" 
attr(z,"climatology") <- "1986-2005"
attr(z,"models") <- vmodel
attr(z,"experiment") <- "cmip3"
r_function <- list(name = "global.tas.anomaly.cmip3.R" , author = "A. Mezghani, METNO , email : abdelkaderm@met.no" , created = "21-03-2013") 
attr(z,"R_Function_attr") <- r_function
attr(z,"Call") <- match.call()
attr(z,"Created") <- date()
if (!is.null(file)) { 
   attr(z,"file") <- file
   save(z,file=file)
}
}
# e.g.
# 

