# documentation in metaextract.R
#' @export metaextract.cordex
metaextract.cordex <- function(files, file.out=NULL, path.out=NULL, 
                               add=TRUE, force=FALSE, verbose=FALSE) {
  if(verbose) print("metaextract.cordex")
  meta.cordex <- NULL
  if(!is.null(path.out) & !is.null(file.out)) file.out <- file.path(path.out, file.out)
  if(!is.null(file.out)) if(add & file.exists(file.out)) {
    load(file.out)
    if(!force) files <- files[!basename(files) %in% meta.cordex$filename]
  }
  if(length(files)>0) {
    filename <- NA; project_id <- NA; dates <- NA; var <- NA; longname <- NA; vunit <- NA; 
    dim <- NA; res <- NA; res.unit <- NA; lon.rng <- NA; lat.rng <- NA; 
    lon.unit <- NA; lat.unit <- NA; gridtype <- NA; 
    frequency <- NA; creation_date <- NA; varid <- NA;
    rcm <- NA; rcm.v <- NA; rcm.domain <- NA; 
    gcm <- NA; gcm.rip <- NA; experiment_id <- NA
    i <- 0
    for(f in files) {
      ncid <- getatt(f)
      nc <- ncdf4::nc_open(f)
      ncid2 <- check.ncdf4(nc)
      ncdf4::nc_close(nc)
      if(tolower(ncid2$model$project_id)=="cordex") {
        i <- i+1
        if(verbose) print(paste("file",i,"of",length(files)))
        filename[i] <- basename(f)
        project_id[i] <- ncid2$model$project_id
        vnames <- names(ncid$var)
        dnames <- names(ncid$dim)
        lon.id <- vnames[vnames %in% c("lon","longitude")]
        lat.id <- vnames[vnames %in% c("lat","latitude")]
        var.id <- vnames[length(vnames)]
        nc <- ncdf4::nc_open(f)
        lon <- ncdf4::ncvar_get(nc, varid=lon.id)
        lat <- ncdf4::ncvar_get(nc, varid=lat.id)
        vatt <- ncdf4::ncatt_get(nc, varid=var.id)
        ncdf4::nc_close(nc)
        dates[i] <- paste0(range(as.character(ncid2$time$vdate)),collapse=",")
        var[i] <- var.id
        longname[i] <- ncid$var[[var.id]]$longname
        vunit[i] <- ncid$var[[var.id]]$units
        dim[i] <- paste(dnames[!grepl("bnds",dnames)], collapse=",")
        res[i] <- mean(diff(ncid$var[[lat.id]]$dim[[1]]$vals),na.rm=TRUE)
        res.unit[i] <- ncid$var[[lat.id]]$dim[[1]]$units
        gridtype[i] <- vatt$grid_mapping
        lon.rng[i] <- paste(range(lon), collapse=",")
        lat.rng[i] <- paste(range(lat), collapse=",")
        lon.unit[i] <- ncid$var[[lon.id]]$units
        lat.unit[i] <- ncid$var[[lat.id]]$units
        frequency[i] <- ncid2$model$frequency
        creation_date[i] <- ncid2$model$creation_date
        rcm[i] <- ncid2$model$model_id
        rcm.v[i] <- ncid2$model$rcm_version_id
        rcm.domain[i] <- ncid2$model$CORDEX_domain
        gcm[i] <- ncid2$model$driving_model_id
        gcm.rip[i] <- ncid2$model$driving_model_ensemble_member
        exp.i <- ncid2$model$experiment_id
        if(exp.i=="historical" & 
           as.numeric(strftime(max(ncid2$time$vdate),"%Y")) > 
           as.numeric(strftime(Sys.time(),"%Y")) ) {
          history.i <- ncid2$model$history
          i.rcp <- regexpr("rcp[0-9]{2}|ssp[0-9]{2}", history.i)
          if(any(i.rcp)) exp.i <- paste(exp.i, substr(history.i, i.rcp[1], 
                           i.rcp[1]+attr(i.rcp,"match.length")-1), sep="+")
        }
        experiment_id[i] <- exp.i
      }
    }
    m <- data.frame(project_id=project_id, filename=basename(files), 
                    dim=dim, dates=dates, var=var, longname=longname,  
                    unit=vunit, resolution=res, resolution_unit=res.unit, gridtype=gridtype, 
                    lon=lon.rng, lon_unit=lon.unit, lat=lat.rng, lat_unit=lat.unit, 
                    experiment_id=experiment_id, frequency=frequency, creation_date=creation_date, 
                    rcm=rcm, rcm.v=rcm.v, rcm.domain=rcm.domain, 
                    gcm=gcm, gcm_rip=gcm.rip)
    if(is.null(meta.cordex)) {
      meta.cordex <- m
    } else {
      meta.cordex <- rbind(meta.cordex, m)
    }
    meta.cordex <- meta.cordex[order(paste(meta.cordex$var, meta.cordex$experiment_id, 
                                           meta.cordex$rcm, meta.cordex$gcm, meta.cordex$gcm_rip, sep=".")), ]
    if(!is.null(file.out)) save(meta.cordex, file=file.out)
  }
  return(meta.cordex)
}
