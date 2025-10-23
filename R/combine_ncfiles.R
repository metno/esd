#' Combine or append netCDF files for stations into one. 
#' Rasmus.Benestad@met.no & chatGPT, , 2025-10-09
#' 
#' @seealso write2ncdf4.station as.station retrieve.station 
#' 
#' @param input_files a string of names of netCDF station files to combine
#' @param output_file netCDF file with the combined station data
#' @param verbose used for debugging and diagnostics of executing the code
#' @param stid_dim identifier for station IDs
#' @param time_dimidentifier for time dimension
#'
#' @keywords utilities
#'
#' @examples
#' \dontrun{
#' Files <- list.files("~/Downloads/ecad2nc4/", pattern = "\\.nc$", full.names = TRUE)
#' 
#' params <- c('t2m','tmax','tmin','precip','sd')
#' for (param in params) { 
#'   files <- Files[grep(param,Files)]
#'   print(files)
#'   output <- paste0(param,".ecad_blend.nc")
#'   if (file.exists(output)) file.remove(output)
#'   combine_nc_stid(files, output)
#' }
#' }
#'
#' @export
combine_nc_stid <- function(input_files,
                            output_file = "combined.nc",
                            stid_dim = "stid",
                            time_dim = "time",
                            verbose = TRUE) {
  require(ncdf4)
  
  if (length(input_files) < 1) stop("No input files provided.")
  if (verbose) message("Combining ", length(input_files), " netCDF files...")
  
  ## ---- 1. Read structure from first file ----
  nc0 <- nc_open(input_files[1])
  on.exit(nc_close(nc0), add = TRUE)
  
  # Dimensions
  dims <- nc0$dim
  vars <- nc0$var
  global_atts <- nc0$gatts
  
  if (!(stid_dim %in% names(dims)))
    stop("stid dimension '", stid_dim, "' not found in file.")
  
  if (!(time_dim %in% names(dims)))
    warning("No time dimension found, assuming only stid dimension variables.")
  
  # Store time values from first file (used for synchronization check)
  time_vals_ref <- if (time_dim %in% names(vars)) ncvar_get(nc0, time_dim) else NULL
  time_len <- length(time_vals_ref)
  
  # Count total stid length across files
  stid_lengths <- numeric(length(input_files))
  time <- rep(NA,2); torg <- c()
  for (i in seq_along(input_files)) {
    nci <- nc_open(input_files[i])
    stid_lengths[i] <- nci$dim[[stid_dim]]$len
    time[1] <- min(time,nci$dim[[time_dim]]$val,na.rm=TRUE)
    time[2] <- max(time,nci$dim[[time_dim]]$val,na.rm=TRUE)
    torg <- c(torg,nci$dim[[time_dim]]$units)
    nc_close(nci)
  }
  stid_total <- sum(stid_lengths)
  if (verbose) message("Total stations (stid): ", stid_total)
  if (verbose) {print(time); print(table(torg))}
  ## Define a time dimension that spans all files
  dims[['time']]$vals <- seq(time[1],time[2],by=1)
  
  ## ---- 2. Define dimensions for output ----
  dims_out <- list()
  for (d in names(dims)) {
    cat('dimension',d,dims[[d]]$units,'time:',range(dims[[d]]$vals),'\n')
    if (d == stid_dim) {
      dims_out[[d]] <- ncdim_def(name = d,
                                 units = dims[[d]]$units,
                                 vals = seq_len(stid_total))
    } else {
      dims_out[[d]] <- ncdim_def(name = d,
                                 units = dims[[d]]$units,
                                 vals = dims[[d]]$vals)
    }
  }
  
  ## ---- 3. Define variables dynamically ----
  vars_out <- list()
  for (vname in names(vars)) {
    v <- vars[[vname]]
    v_dims <- lapply(v$dim, function(d) dims_out[[d$name]])
    if (v$prec=='int') v$prec <- 'integer'
    vars_out[[vname]] <- ncvar_def(name = vname,
                                   units = v$units,
                                   dim = v_dims,
                                   missval = v$missval,
                                   longname = v$longname,
                                   prec = v$prec,
                                   compression = 4)
  }
  
  ## ---- 4. Create output file ----
  nc_out <- nc_create(output_file, vars = vars_out)
  
  # Copy global attributes
  for (att_name in names(global_atts)) {
    ncatt_put(nc_out, 0, att_name, global_atts[[att_name]])
  }
  ncatt_put(nc_out, 0, "history",
            paste("Combined by combine_nc_stid on", Sys.time()))
  
  ## ---- 5. Loop through files, appending along stid ----
  stid_start <- 1
  
  for (i in seq_along(input_files)) {
    f <- input_files[i]
    if (verbose) message("[", i, "/", length(input_files), "] ", f)
    nc_in <- nc_open(f)
    
    # --- Synchronize time ---
    if (!is.null(time_vals_ref) && time_dim %in% names(nc_in$var)) {
      time_vals <- ncvar_get(nc_in, time_dim)
      if (length(time_vals) != length(time_vals_ref) ||
          any(time_vals != time_vals_ref)) {
        stop("Time coordinate mismatch between ", input_files[1],
             " and ", f, ". Files must have identical time values.")
      }
    }
    
    # --- Determine stid slice ---
    nst <- nc_in$dim[[stid_dim]]$len
    stid_end <- stid_start + nst - 1
    
    # --- Copy every variable ---
    for (vname in names(vars_out)) {
      vdef <- vars_out[[vname]]
      vinfo <- nc_in$var[[vname]]
      if (is.null(vinfo)) next  # skip if var not in current file
      
      dims_v <- sapply(vdef$dim, `[[`, "name")
      start <- rep(1, length(dims_v))
      count <- sapply(vdef$dim, `[[`, "len")
      
      # Adjust for stid dimension
      if (stid_dim %in% dims_v) {
        stid_index <- which(dims_v == stid_dim)
        start[stid_index] <- stid_start
        count[stid_index] <- nst
      }
      
      # Read from input (only this slice)
      data_in <- ncvar_get(nc_in, vname)
      
      # Write into output
      if (length(dim(data_in))==2) if (count[2] != dim(data_in)[2]) count[2] <- dim(data_in)[2]
      cat('ncvar_put:',vname,'dim:',dim(data_in),' start: ',start,'count:',count,'\n')
      ncvar_put(nc_out, vname, data_in,
                start = start,
                count = count)
    }
    
    nc_close(nc_in)
    stid_start <- stid_end + 1
  }
  
  ## ---- 6. Close output ----
  nc_close(nc_out)
  if (verbose) message("Combined file written to ", output_file)
  
  invisible(output_file)
}

