#' Combine or append netCDF files for stations into one. 
#' Rasmus.Benestad@met.no & chatGPT, , 2025-10-09
#' 
#' @seealso write2ncdf4.station as.station retrieve.station 
#' 
#' @param input_files a string of names of netCDF station files to combine
#' @param output_file netCDF file with the combined station data
#' @param verbose used for debugging and diagnostics of executing the code
#'
#' @keywords utilities
#'
#' @examples
#' \dontrun{
#' concat_nc_along_stid(ncfiles, paste0('~/Downloads/',variable,'.ecad_blend.nc'),verbose=TRUE)
#' }
#'
#' @export
concat_nc_along_stid <- function(input_files, output_file, verbose = TRUE) {
  # library(ncdf4)
  require(abind)
  # library(zoo)
  
  if (verbose) cat('concat_nc_along_stid \n')
  if (length(input_files) < 2)
    stop("Please provide at least two NetCDF files.")
  
  ncs <- lapply(input_files, nc_open)
  #on.exit(lapply(ncs, nc_close), add = TRUE)
  
  template <- ncs[[1]]
  dims <- template$dim
  vars <- template$var
  
  # Collect all stid and time values
  all_stids <- unlist(lapply(ncs, function(nc) ncvar_get(nc, "stid")))
  all_stids <- 1:length(all_stids)
  if (verbose) print(all_stids)
  all_times <- unique(unlist(lapply(ncs, function(nc) {
    ncvar_get(nc, "time") 
  })))
  all_times <- sort(na.omit(all_times))
  if (verbose) {cat('Times: '); cat(range(all_times)); cat('\n')}
  # --- Dynamically define ALL dimensions from template ---
  if (verbose) cat('Define all dimensions: \n')
  dim_defs <- list()
  for (dname in names(dims)) {
    if (verbose) cat(dname,'... \n')
    d <- dims[[dname]]
    if (dname == "stid") {
      dim_defs[[dname]] <- ncdim_def("stid", d$units, vals = all_stids, create_dimvar = TRUE)
    } else if (dname == "time") {
      # Convert all_times to numeric
      all_times <- sort(as.numeric(na.omit(all_times)))
      dim_defs[[dname]] <- ncdim_def("time", d$units, vals = all_times, create_dimvar = TRUE)
    } else {
      dim_defs[[dname]] <- ncdim_def(dname, d$units, vals = d$vals, create_dimvar = TRUE)
    }
  }
  
  # # --- Dynamically define ALL dimensions from template ---
  # dim_defs <- list()
  # for (dname in names(dims)) {
  #   d <- dims[[dname]]
  #   if (dname == "stid") {
  #     dim_defs[[dname]] <- ncdim_def("stid", d$units, vals = all_stids, create_dimvar = TRUE)
  #   } else if (dname == "time") {
  #     dim_defs[[dname]] <- ncdim_def("time", d$units, vals = all_times, create_dimvar = TRUE)
  #   } else {
  #     dim_defs[[dname]] <- ncdim_def(dname, d$units, vals = d$vals, create_dimvar = TRUE)
  #   }
  # }
  
  # --- Define variables ---
  if (verbose) cat('Define all variables: \n')
  var_defs <- list()
  for (vname in names(vars)) {
    if (verbose) cat(vname,'... \n')
    v <- vars[[vname]]
    var_dims <- lapply(v$dim, function(d) {
      if (!d$name %in% names(dim_defs))
        stop(sprintf("Missing dimension '%s' for variable '%s'", d$name, vname))
      dim_defs[[d$name]]
    })
    
    missval <- if (!is.null(v$missval)) v$missval else NA
    if (v$prec=='int') v$prec <- 'integer'
    var_defs[[vname]] <- ncvar_def(v$name, v$units, dim = var_dims,
                                   longname = v$longname, prec = v$prec,
                                   missval = missval)
  }
  
  # --- Create output file ---
  if (verbose) cat("Creating output file:", output_file, "\n")
  nc_out <- nc_create(output_file, var_defs, force_v4 = TRUE)
  #on.exit(try(nc_close(nc_out), silent = TRUE), add = TRUE)
  
  # --- Merge and write variables ---
  for (vname in names(vars)) {
    if (verbose) cat("Merging variable:", vname, "\n")
    data_list <- list()
    for (i in seq_along(ncs)) {
      nc <- ncs[[i]]
      if (!(vname %in% names(nc$var))) next
      vals <- ncvar_get(nc, vname)
      time <- ncvar_get(nc, 'time')
      torg <- sub('days since ','',ncatt_get(nc,'time','units')$value)
      time <- as.Date(time,origin = torg)
      attr(vals,'time') <- time
      data_list[[i]] <- vals
    }
    
    # Concatenate only if variable has stid dimension
    if (verbose) cat("Synchronise and pad the data from different nc files \n")
    dims <- sapply(vars[[vname]]$dim, function(d) d$name)
    ## Make sure that all data chuncks cover the whole period
    # 1. Get the union of all time indices
    all_times <- as.Date(sort(unique(unlist(lapply(data_list, function(m) attr(m, "time"))))))
    if (verbose) {cat('Combined time span: \n'); print(range(all_times))}
    # 2. Align each matrix to the full time vector
    if ('time' %in% dims) { 
      data_list_synced <- lapply(data_list, function(m) {
        time_m <- attr(m, "time")
        ny <- nrow(m)
        
        # create output matrix filled with NAs
        out <- matrix(NA, ncol = length(all_times), nrow = ny)
        colnames(out) <- all_times
        
        # find matching time indices
        match_idx <- match(time_m, all_times)
        if (verbose) cat('length(out[,match_idx])',length(out[,match_idx]),' length(m)',length(m),' \n')
        out[,match_idx] <- m
        
        # keep the aligned time attribute
        attr(out, "time") <- all_times
        out
      })
      data_list <- data_list_synced
    }
    if (verbose) cat("Check for dimension 'time' \n")
    if ("time" %in% dims) {
      stid_idx <- which(dims == "stid")
      combined <- abind(data_list, along = stid_idx)
    } else {
      combined <- unlist(data_list)
    }
    if (verbose) cat("ncvar_put ", vname, " length=",length(combined),"=",nc_out$dim$stid$len,"\n")
    ncvar_put(nc_out, vname, combined)
  }
  
  nc_close(nc_out)
  cat("Combined NetCDF written to", output_file, "\n")
  invisible(output_file)
}
