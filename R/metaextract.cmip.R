# Documenation in metaextract.R
#' @export metaextract.cmip
metaextract.cmip <- function(x, verbose=FALSE) {
  if(verbose) print("metaextract.cmip")
  ## argument 'x' is input from getCM, getGCMs, getRCMs etc
  
  project_id <- NA; mip_era <- NA; filename <- NA; dim <- NA; dates <- NA
  var <- NA; longname <- NA; vunit <- NA; vid <- NA
  res <- NA; lon.rng <- NA; lon.unit <- NA; lat.rng <- NA; lat.unit <- NA
  experiment <- NA; experiment_id <- NA; 
  frequency <- NA; creation_date <- NA; tracking_id <- NA
  gcm <- NA; gcm.v <- NA; realm <- NA
  gcm.rip <- NA; realization <- NA; initialization <- NA; physics <- NA; forcing <- NA;
  qf <- NULL
  if(!is.null(x$dim)) dim <- paste(names(x$dim),collapse=",")
  if(!is.null(names(x$var))) {
    var <- names(x$var)
    var <- var[!grepl("time|lon|lat|height",var)]
    if(!is.null(x$var[[1]]$longname)) longname <- sapply(var, function(v) x$var[[v]]$longname)
    if(!is.null(x$var[[1]]$units)) vunit <- sapply(var, function(v) x$var[[v]]$units)
    if(!is.null(x$var[[1]]$id$id)) vid <- sapply(var, function(v) x$var[[v]]$id$id)
  }
  if(!is.null(names(x$dim))) {
    if(!is.null(x$dim$lat$vals)) {
      res <- diff(x$dim$lat$vals)[1]
      lat.rng <- paste(range(x$dim$lat$vals),collapse=",")
    }
    if(!is.null(x$dim$lon$vals)) lon.rng <- paste(range(x$dim$lon$vals),collapse=",")
    if(!is.null(x$dim$lat$units)) lat.unit <- x$dim$lat$units
    if(!is.null(x$dim$lon$units)) lon.unit <- x$dim$lon$units
  }

  for(mi in c("filename","dates","frequency",
              "mip_era","project_id","experiment","experiment_id",
              "creation_date","tracking_id")) {
    if(!is.null(x[[mi]])) {
      eval(parse(text=paste(mi," <- x$",mi,sep="")))
    } else if (!is.null(x$model[[mi]])) {
      eval(parse(text=paste(mi," <- x$model$",mi,sep="")))
    }
  }
  
  if(any(grepl("CMIP[0-9]{1}",project_id))) project_id <- project_id[grep("CMIP[0-9]{1}",project_id)]
  for(mi in c("realization","initialization","physics","forcing","realm")) {
    if(any(grepl(mi,names(x$model)))) {
      xi <- x$model[grep(mi,names(x$model))]#[1]]
    } else if(any(grepl(mi,names(x)))) {
      xi <- x[grep(mi,names(x$model))]#[1]]
    } else {
      xi <- NA
    }
    if(length(xi)>1) {
      if(any(grepl("[0,9]{1,2}",xi))) {
        xi <- xi[grep("[0,9]{1,2}",xi)]
      } else {
        xi <- NA
      }
    } 
    eval(parse(text=paste(mi," <- xi",sep="")))
  }
    
  if(frequency=="month") frequency <- "mon"
  
  ## KMP 2022-05-02: If defined, using parent_source_id instead of model_id because 
  ## the latter is wrong in some CMIP6 files. For CMIP5, parent_source_id is not defined.
  if(!is.null(x$model$parent_source_id)) {
    gcm <- x$model$parent_source_id
  } else if(!is.null(x$model$model_id)) {
    gcm <- x$model$model_id
  }
  if(any(grep("parent_experiment_rip",names(x$model)))) {
    nm <- names(x$model)[grep("parent_experiment_rip",names(x$model))[1]]
    gcm.rip <- x$model[[nm]]
  } else if(any(grep("variant",names(x$model)))) {
    nm <- names(x$model)[grep("variant",names(x$model))[1]]
    gcm.rip <- x$model[[nm]]
  }
  
  ## Version number:
  if(any(grepl("version",names(x$model)))) {
    iv <- which(grepl("version",names(x$model)) & !grepl("physics_version",names(x$model)))
    gcm.v <- paste(paste(names(x$model)[iv],paste(x$model[iv])),collapse=", ")
  }

  ## If information is missing in netCDF header, use the model history to obtain it. 
  if(!is.null(x$model$history)) {
    h <- x$model$history
    creation_date <- substr(h, 1, regexpr("[0-9]{4}",h)+3)
    h <- strsplit(h," ")
    fname <- unlist(h)[grep("Amon.*.nc",unlist(h))[1]]
    if(is.null(x$model$model_id)) {
      qf <- c(qf, paste("Information about model missing from netCDF header.",
                        "Metadata recovered from history attribute."))
      frequency <- "mon"
      longname <- switch(var,"tas"="Near-Surface Air Temperature",
                         "pr"="Precipitation")
    }
    if(is.na(gcm)) gcm <- gsub("_.*","",gsub(".*_Amon_","",fname))
    if(is.na(gcm.rip)) gcm.rip <- unlist(strsplit(substr(fname,regexpr("r[0-9]{1,2}i[0-9]{1,2}p[0-9]{1,2}",fname)[1],
                      nchar(fname)),split="_"))[1]
    if(is.na(experiment)) {
      if(project_id=="CMIP6" & grepl("ssp",fname)) {
        experiment <- unlist(strsplit(substr(fname,regexpr("ssp",fname)[1],nchar(fname)),split="_"))[1]
        experiment <- toupper(experiment)
      } else if(project_id=="CMIP5" & grepl("rcp",fname)) {
        experiment <- unlist(strsplit(substr(fname,regexpr("rcp",fname)[1],nchar(fname)),split="_"))[1]
        N <- nchar(experiment)
        experiment <- paste(substr(experiment,1,N-1),substr(experiment,N,N),sep=".")
        experiment <- toupper(experiment)
      } else if(min(x$dates)<as.Date("2010-01-01")) {
        experiment <- "historical"
      }
    }
    if(is.na(experiment_id)) experiment_id <- experiment
  } 
  experiment_id <- tolower(sub("[.]","",experiment_id))
  if(min(x$dates)<as.Date("2010-01-01") & !grepl("historical",experiment_id)) {
    experiment_id <- paste("historical",experiment_id,sep="+")
  }

  ## Check and correct rip - some simulations have the wrong rip attached.
  #if (is.na(gcm.rip)) qf <- c(qf,"Missing experiment_rip in netCDF header.")
  if (is.na(gcm.rip)) qf <- c(qf,"Missing experiment_rip in netCDF header.")
  if(!is.na(realization) & !is.na(initialization) & !is.na(physics)) {
    gcm.rip2 <- paste0("r", realization,
                       "i", initialization,
                       "p", physics)
    if(project_id=="CMIP6" & !is.na(forcing)) gcm.rip2 <- paste0(gcm.rip2,"f",forcing)
    if (is.na(gcm.rip) | gcm.rip!=gcm.rip2) {
      gcm.rip <- gcm.rip2
      qf <- c(qf,paste("Discrepancy in experiment_rip in netCDF header.",
                       "Replaced experiment_rip with realization, intialization, physics",
                       "(and forcing for CMIP6 files)."))
    }
  }
  
  ## An extra check of the experiment_rip:
  h <- strsplit(x$model$history," ")
  fname <- unlist(h)[grep(".*.nc",unlist(h))[1]]
  if(grepl("r[0-9]{1,2}i[0-9]{1,2}p[0-9]{1,2}",fname)) {
    gcm.rip3 <- unlist(strsplit(substr(fname,regexpr("r[0-9]{1,2}i[0-9]{1,2}p[0-9]{1,2}",fname)[1],
                       nchar(fname)),split="_"))[1]
    if(is.na(gcm.rip)) {
      gcm.rip <- gcm.rip3
      qf <- c(qf,"Missing experiment_rip in netCDF header. Using information from model history.")
    } else if(gcm.rip!=gcm.rip3) {
      gcm.rip <- gcm.rip3
      qf <- c(qf,paste("Discrepancy in experiment_rip in netCDF header.",
                       "Replaced experiment_rip with information from model history."))
    }
  }
  if(!is.na(filename)) filename <- gsub(".*/","",filename)
  mx <- data.frame(project_id=paste(project_id,collapse=","), 
                   filename=filename, 
                   dim=paste(dim,collapse=","), dates=dates, var=paste(var,collapse=","),
                   longname=paste(longname,collapse=","), unit=paste(vunit,collapse=","),
                   resolution=res, lon=lon.rng, lon_unit=lon.unit, 
                   lat=lat.rng, lat_unit=lat.unit,
                   experiment=experiment, 
                   experiment_id=experiment_id, 
                   frequency=frequency, 
                   creation_date=creation_date, 
                   gcm=gcm, gcm_rip=gcm.rip, qf=paste(qf,collapse="; "))
  X <- matrix(sapply(mx,as.character), nrow=1, ncol=length(mx))
  colnames(X) <- colnames(mx)
  return(X)
}
