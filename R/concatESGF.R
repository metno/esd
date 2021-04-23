#' concat.ESGF
#' Concatenate netcdf files from ESGF into a combined netcdf file of historical + projected 
#' climate simulations using CDO command line
#' 
#' This function assumes that the data has been downloaded locally using the functions \code{meta.ESGF} and \code{retrieve.ESGF}. 
#' 
#' @import ncdf4
#'
#' @seealso meta.ESGF retrieve.ESGF
#'
#' @param meta Output data.frame from meta.ESGF() function. 
#' The full path to the netCDF file can either be included in 'ncfile' or entered as a separate input ('path').
#' @param path Path to netcdf files
#' @param param Climate parameter as defined in ESGF.
#' 
#' @return NetCDF file output containing both historical and future climate simulations.
#'
#' @export
concat.ESGF <- function(meta, param='tas', path='/lustre/storeB/project/CMIP/CMIP6.monthly/from_synda/CMIP6',
            expid = 'ssp585', use = 'CDO',verbose = FALSE) {
  # Group by model and member.id 
  splitlist <- split(meta,f=list(meta$member.id,meta$model),drop=T)
  # Print out the number of groups/experiments
  if (verbose) print(paste('Number of experiments found is', length(splitlist)))
  # Loop over all model-model.id groups
  for (i in 1:length(splitlist)) {
    # Get the list of files for a specific experiment                     
    lfiles <- splitlist[[i]]$title
    nfiles <- length(lfiles)
    urls <- splitlist[[i]]$OpenDap
    urls.fix <- sapply(1:nfiles,function(x) return(paste(unlist(strsplit(paste(urls[x]),split = '/'))[7:16],collapse='/')))
    urls.upd <- file.path(path,urls.fix)
                       
    browser()
    if (!file.exists(urls.upd)) retrieve.ESGF(meta = meta, path = urls.upd,)
    
    if (verbose) 
      print(paste('Processing',nfiles,'files for', unique(splitlist[[i]]$model),unique(splitlist[[i]]$member.id)))
                       
    # Get the periods
    it <- splitlist[[i]]$period
    it.len <- length(it)
    # rank the periods
    srt <- order(it)
    it.srt <- it[srt]
    it1 <- substr(it.srt[1],1,6)
    it2 <- substr(it.srt[it.len],8,nchar(as.character(it.srt[it.len])))
    # order the files according to the periods
    lfiles.srt <- lfiles[srt]
    # create output filename
    ofile <- file.path(path, paste(gsub('historical',paste('historical',expid,sep='-'),
                                        substr(lfiles.srt[1],1,nchar(as.character(lfiles.srt[1]))-16)),
                                   paste(it1,it2,sep='-'),'.nc',sep=''))
    if (verbose) print(paste('Saving to file', ofile))
    # Excecute CDO job
    #eval(parse(text = paste("system('cdo", paste(lfiles.srt,collapse = ' '), ofile,"')")
  }
}
