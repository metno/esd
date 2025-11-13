#' Combine or append netCDF files for stations into one. 
#' 
#' This revised version correctly concatenates and writes the station coordinate
#' variables (stid, lat, lon) to preserve geographical context.
#' 
#' @param input_files a string of names of netCDF station files to combine
#' @param output_file netCDF file with the combined station data
#' @param stid_dim identifier for station IDs
#' @param time_dim identifier for time dimension
#' @param verbose used for debugging and diagnostics of executing the code
#'
#' @export
combine_nc_stid <- function(input_files,
                            output_file = "combined.nc",
                            stid_dim = "stid",
                            time_dim = "time",
                            verbose = FALSE) {
  require(ncdf4)
  
  if (length(input_files) < 1) stop("No input files provided.")
  if (verbose) message("Combining ", length(input_files), " netCDF files...")
  
  # --- METADATA STORAGE ---
  stid_vals_all <- c() # Stores actual station IDs (e.g., "LA0004")
  lat_vals_all <- c()  # Stores latitude coordinates
  lon_vals_all <- c()  # Stores longitude coordinates
  stid_lengths <- numeric(length(input_files))
  
  ## ---- 1. Read structure from first file and collect ALL metadata ----
  nc0 <- nc_open(input_files[1])
  on.exit(nc_close(nc0), add = TRUE)
  
  dims <- nc0$dim
  vars <- nc0$var
  global_atts <- nc0$gatts
  
  if (!(stid_dim %in% names(dims)))
    stop("stid dimension '", stid_dim, "' not found in file.")
  
  time <- rep(NA,2); torg <- c()
  
  # Loop to collect total stid length, time range, AND coordinate values
  for (i in seq_along(input_files)) {
    nci <- nc_open(input_files[i])
    
    # Check for required coordinates
    if (!("lat" %in% names(nci$var) && "lon" %in% names(nci$var))) {
      warning("File ", input_files[i], " is missing 'lat' or 'lon' variables. Skipping.")
      #nc_close(nci)
      next
    }
    
    # Collect metadata
    nst <- nci$dim[[stid_dim]]$len
    stid_lengths[i] <- nst
    
    stid_vals_all <- c(stid_vals_all, ncvar_get(nci, stid_dim))
    lat_vals_all <- c(lat_vals_all, ncvar_get(nci, "lat"))
    lon_vals_all <- c(lon_vals_all, ncvar_get(nci, "lon"))
    
    # Time collection (for range check)
    if (time_dim %in% names(nci$dim)) {
      time[1] <- min(time[1], nci$dim[[time_dim]]$val, na.rm=TRUE)
      time[2] <- max(time[2], nci$dim[[time_dim]]$val, na.rm=TRUE)
      torg <- c(torg, nci$dim[[time_dim]]$units)
    }
    nc_close(nci)
  }
  
  stid_total <- sum(stid_lengths)
  if (verbose) message("Total stations (stid): ", stid_total)
  
  # Define a time dimension that spans all files
  if (time_dim %in% names(dims)) {
    # Use the collected time range
    dims[[time_dim]]$vals <- seq(time[1], time[2], by = 1) 
  }
  
  ## ---- 2. Define dimensions for output ----
  dims_out <- list()
  for (d in names(dims)) {
    if (d == stid_dim) {
      # --- FIX: Use the actual concatenated station IDs (stid_vals_all) ---
      # This is critical to maintain the geographical context.
      dims_out[[d]] <- ncdim_def(name = d,
                                 units = dims[[d]]$units,
                                 vals = stid_vals_all, # Use actual station ID values
                                 create_dimvar = TRUE) # Ensure stid variable is created
    } else {
      dims_out[[d]] <- ncdim_def(name = d,
                                 units = dims[[d]]$units,
                                 vals = dims[[d]]$vals)
    }
  }
  
  ## ---- 3. Define variables dynamically ----
  vars_out <- list()
  
  # Define coordinate variables (lat and lon) which are NOT included in nc0$var automatically
  stid_dim_out <- dims_out[[stid_dim]]
  
  # vars_out[["lat"]] <- ncvar_def(name = "lat",
  #                                units = "degrees_north",
  #                                dim = list(stid_dim_out),
  #                                longname = "latitude of station",
  #                                prec = "float")
  # 
  # vars_out[["lon"]] <- ncvar_def(name = "lon",
  #                                units = "degrees_east",
  #                                dim = list(stid_dim_out),
  #                                longname = "longitude of station",
  #                                prec = "float")
  
  # Define all main data variables
  if (verbose) cat('Define ')
  for (vname in names(vars)) {
    if (verbose) cat(vname)
    v <- vars[[vname]]
    v_dims <- lapply(v$dim, function(d) dims_out[[d$name]])
    if (v$prec=='int') v$prec <- 'integer'
    if (v$prec=='short') v$prec <- 'float' # Use float/real for output to store scaled data
    
    vars_out[[vname]] <- ncvar_def(name = vname,
                                   units = v$units,
                                   dim = v_dims,
                                   missval = v$missval,
                                   longname = v$longname,
                                   prec = v$prec,
                                   compression = 4)
  }
  if (verbose) cat('... al variables defined \n')
  
  ## ---- 4. Create output file ----
  nc_out <- nc_create(output_file, vars = vars_out)
  
  # --- FIX: Write the collected coordinate data immediately ---
  # ncvar_put(nc_out, "lat", lat_vals_all, start = 1, count = stid_total)
  # ncvar_put(nc_out, "lon", lon_vals_all, start = 1, count = stid_total)
  
  # Copy global attributes
  if (verbose) cat('Add global attributes... ')
  for (att_name in names(global_atts)) {
    ncatt_put(nc_out, 0, att_name, global_atts[[att_name]])
  }
  ncatt_put(nc_out, 0, "history",
            paste("Combined by combine_nc_stid on", Sys.time()))
  if (verbose) cat('\n')
  
  ## ---- 5. Loop through files, appending along stid ----
  stid_start <- 1
  
  if (verbose) cat('Add from ncfile ')
  for (i in seq_along(input_files)) {
    f <- input_files[i]
    if (verbose) message("[", i, "/", length(input_files), "] ", f)
    nc_in <- nc_open(f)
    #on.exit(nc_close(nc_in), add = TRUE) # Ensure closure if loop fails
    
    # --- Determine stid slice ---
    nst <- nc_in$dim[[stid_dim]]$len
    stid_end <- stid_start + nst - 1
    
    # --- Copy every data variable (excluding lat/lon/stid, which were written above) ---
    for (vname in names(vars)) { # Only iterate over original data variables
      if (verbose) cat(vname,': ')
      vdef <- vars_out[[vname]]
      vinfo <- nc_in$var[[vname]]
      if (is.null(vinfo)) next
      
      dims_v <- sapply(vdef$dim, `[[`, "name")
      start <- rep(1, length(dims_v))
      count <- sapply(vdef$dim, `[[`, "len")
      
      # Adjust for stid dimension
      if (stid_dim %in% dims_v) {
        stid_index <- which(dims_v == stid_dim)
        start[stid_index] <- stid_start
        count[stid_index] <- nst
      }
      
      # Read data
      data_in <- ncvar_get(nc_in, vname, raw_datavals=TRUE)
      
      # --- Data Scaling/Offset Handling (Retained from original script) ---
      # This block handles reading raw scaled integers and applying offset/scaling
      if (vname==names(vars[1])) {
        offset <- try(ncatt_get(nc_in, vname,"add_offset")$value)
        scaling <- try(ncatt_get(nc_in, vname,"scale_factor")$value)
        missing <- try(ncatt_get(nc_in, vname,"missing_value")$value)
        if (!inherits(scaling,'missing')) data_in[round(data_in)==round(missing)] <- NA
        if (!inherits(scaling,'try-error')) data_in <- data_in *scaling
        if (!inherits(offset,'try-error')) data_in <- data_in + offset
        if (!inherits(scaling,'missing')) data_in[round(data_in)==round(missing)] <- NA
      }
      if (verbose) cat('Data range: ',range(c(data_in),na.rm=TRUE),'\n')
      
      # Write into output
      # Ensure count matches data dimensions, especially if time dimension varies slightly
      if (time_dim %in% dims_v) {
        time_index <- which(dims_v == time_dim)
        count[time_index] <- dim(data_in)[time_index]
      }
      #if (is.numeric(data_in)) plot(data_in,main=vname,type='b')
      
      if (verbose) cat('Start: ',start,' Count: ',count,'\n')
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