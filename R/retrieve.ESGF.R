#' Retrieve CMIP data directly from the Earth System Grid Federation (ESGF)
#' https://earthsystemcog.org/projects/cog/esgf_search_restful_api
#' meta.ESGF returns a data.frame with the model metadata and the OpenDAP URL that 
#' can be used with retrieve.ESGF. The function \code{retrieve.ESGF} is a wraparound
#' for \code{retrieve} that reads several files belonging to the same model and run.
#'
#' @aliases retrieve.ESGF meta.ESGF
#'
#' @import ncdf4
#'
#' @seealso retrieve 
#'
#' @param param Name of parameter
#' @param verbose Logical value defaulting to FALSE. If FALSE, do not display
#' comments (silent mode). If TRUE, displays extra information on progress.
#' @param url The base URL for ESGF
#' @param n Number of models to read - NULL reads everything
#' @param expid Name of the MIP experiment ('historical' for historical run)
#' @param mip CMIP5 or CMIP6
#' @param freq Frequency of data 
#' @return A data.frame for \code{meta.ESGF} and a "zoo" "field" object with additional attributes used for further
#' processing for +code{retrieve.ESGF}.
#'
#' @examples 
#' \dontrun{
#' meta <- meta.ESGF(n=3)
#' X <- retrieve.ESGF(im=3,lon=c(-30,40),lat=c(50,70),meta=meta)
#' map(X)
#'}
#' @export retrieve.ESGF
retrieve.ESGF <- function(im=1,meta=NULL,verbose=FALSE,...) { 
  if (is.null(meta)) meta <- meta.ESGF(verbose=verbose,...)
  if (inherits(im,c('numeric','integer'))) model <- as.character(meta$model[im]) else
    if (is.character(im)) { 
      if (verbose) print(im)
      ## Assume the shape '<model>_<expid>_<ensid>'
      im <- strsplit(gsub('_',' ',im),' ')
      model <- im[[1]]; expid <- im[[2]]; ensid <- im[[3]]
      im <- intersect( grep(model,meta$model), grep(ensid,meta$member.id) )
      model <- meta$model[im]
      if (verbose) print(c(im,model,ensid))
      if (is.na(im)) return(NULL)
    }
  mem <- as.character(meta$member.id[im])
  jm <- (1:length(meta$model))[is.element(as.character(meta$model),model) &
                               is.element(as.character(meta$member.id),mem)]
  for (j in jm) { 
    print(as.character(meta$period[j]))
    opendap <- as.character(meta$OpenDap[j])
    if (verbose) print(opendap)
    x <- try(retrieve(file=opendap,param=attr(meta,'variable'),verbose=verbose,...))
    if ( (!inherits(x,"try-error")) & (length(jm)>1) ) { 
      if (j == jm[1]) X <- x else X <- c(zoo(X),zoo(x))
      X <- attrcp(x,X); class(X) <- class(x)
    } else if ( (length(jm)==1) | (inherits(x,"try-error")) ) X <- x 
  }
  invisible(X)
}

#' Read metadata from ESGF
#'
#' @seealso retrieve.ESGF
#'
#' @param url URL
#' @param mip CMIP5 or CMIP6
#' @param param meteorological variable name, e.g., "tas"
#' @param freq temporal frequency, e.g., "mon"=monthly
#' @param expid scenario (ssp for CMIP6, rcp for CMIP5)
#' @param verbose if TRUE print progress
#' @param n number of datasets
#'
#' @export meta.ESGF
meta.ESGF <- function(url="https://esgf-data.dkrz.de/esg-search/search/",mip="CMIP6",param="tas",
                      freq="mon",expid="ssp585",verbose=FALSE,n=NULL) {
  if (verbose) print('meta.ESFG - this function uses the jsonlite package to read metadata from ESGF')

  ## Check if the JSON library is installed
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' needed to use 'meta.ESGF'. Please install it.")
  } else {
    
    ## Get number of available datasets
    ### Facet switches: CMIP5 | CMIP6 | else
    facet_cmip <- switch(mip,"CMIP5" = "project", "CMIP6" = "mip_era", "project") 
    facet_var  <- switch(mip,"CMIP5" = "variable", "CMIP6" = "variable_id")
    facet_freq <- switch(mip,"CMIP5" = "time_frequency", "CMIP6" = "frequency")
    facet_exp  <- switch(mip,"CMIP5" = "experiment", "CMIP6" = "experiment_id")
    facet_mod  <- switch(mip,"CMIP5" = "model", "CMIP6" = "source_id")
    facet_ens  <- switch(mip,"CMIP5" = "ensemble", "CMIP6" = "member_id")
    facet_ins  <- switch(mip,"CMIP5" = "institute", "CMIP6" = "institution_id")

    search_string <- paste("type=Dataset&replica=false&latest=true&",facet_cmip,"=",mip,"&",facet_var,"=",param,
                         "&",facet_freq,"=",freq,"&",facet_exp,"=",expid,"&format=application%2Fsolr%2Bjson",sep="")
  
    if (verbose) print(search_string)

    ## Get number of available datasets
    nof_datasets <- jsonlite::fromJSON(paste(url,"?",search_string,sep=""))$response$numFound

    if (!is.null(n)) nof_datasets <- n
    print(paste('Found',nof_datasets,'datasets'))
    if (nof_datasets==0) {
      print('Please check the search URL (CMIP, experiment, variable, facets, etc.):')
      print(paste(url,search_string,sep=""))      
      print('If all looks OK, maybe the connection is down - try again later!')
      return(NULL)
    }
    
    ## Query ESGF server
    ESGF_query <- jsonlite::fromJSON(paste(url,"?limit=",nof_datasets,"&",search_string,sep=""))
    if (verbose) str(ESGF_query)
    
    ## List files in each dataset and link to files from the corresponding historical dataset
    #nof_datasets <- 7 ## test
    results <- list()
    for (i in 1:nof_datasets) {
      nof_files <- jsonlite::fromJSON(paste(url,"?limit=0&type=File&replica=false&latest=true&",facet_cmip,"=",mip,"&",facet_var,"=",param,"&",facet_freq,"=",freq,"&",facet_exp,"=",expid,"&",facet_mod,"=",ESGF_query$response$docs[i,facet_mod],"&",facet_ens,"=",ESGF_query$response$docs[i,facet_ens],"&",facet_ins,"=",ESGF_query$response$docs[i,facet_ins],"&format=application%2Fsolr%2Bjson",sep=""))$response$numFound
      
      ESGF_file_query <- jsonlite::fromJSON(paste(url,"?limit=",nof_files,"&type=File&replica=false&latest=true&",facet_cmip,"=",mip,"&",facet_var,"=",param,"&",facet_freq,"=",freq,"&",facet_exp,"=",expid,"&",facet_mod,"=",ESGF_query$response$docs[i,facet_mod],"&",facet_ens,"=",ESGF_query$response$docs[i,facet_ens],"&",facet_ins,"=",ESGF_query$response$docs[i,facet_ins],"&format=application%2Fsolr%2Bjson",sep=""))
      
      nof_hist_files <- jsonlite::fromJSON(paste(url,"?limit=0&type=File&replica=false&latest=true&",facet_cmip,"=",mip,"&",facet_var,"=",param,"&",facet_freq,"=",freq,"&",facet_exp,"=historical&",facet_mod,"=",ESGF_query$response$docs[i,facet_mod],"&",facet_ens,"=",ESGF_query$response$docs[i,facet_ens],"&",facet_ins,"=",ESGF_query$response$docs[i,facet_ins],"&format=application%2Fsolr%2Bjson",sep=""))$response$numFound
      
      ESGF_hist_file_query <- jsonlite::fromJSON(paste(url,"?limit=",nof_hist_files,"&type=File&replica=false&latest=true&",facet_cmip,"=",mip,"&",facet_var,"=",param,"&",facet_freq,"=",freq,"&",facet_exp,"=historical&",facet_mod,"=",ESGF_query$response$docs[i,facet_mod],"&",facet_ens,"=",ESGF_query$response$docs[i,facet_ens],"&",facet_ins,"=",ESGF_query$response$docs[i,facet_ins],"&format=application%2Fsolr%2Bjson",sep=""))
      
      if (verbose) {
        print(paste(rep("=",nchar(ESGF_query$response$docs$id[i])+9),collapse=""))
        print(paste("Dataset:",ESGF_query$response$docs$id[i]))
        print("Files:")
        print(ESGF_file_query$response$docs$title)
        print(paste("Historical dataset:",ESGF_hist_file_query$response$docs$dataset_id[1]))
        print("Files:")
        print(ESGF_hist_file_query$response$docs$title)
        print(paste(rep("=",nchar(ESGF_query$response$docs$id[i])+9),collapse=""))
      } else cat('.')
      
      ic <- as.character(i); if (i < 100) ic <- paste('0',ic,sep=''); if (i < 10) ic <- paste('0',ic,sep='')
      results[[paste('dataset.query:',ic,sep='_')]] <- ESGF_query$response$docs[i,]
      
      for (j in 1:nof_files) {
        jc <- as.character(j); if (j < 100) jc <- paste('0',jc,sep=''); if (j < 10) jc <- paste('0',jc,sep='')
        
        results[[paste('file.query:',ic,jc,sep='_')]] <- ESGF_file_query$response$docs[j,]
        
        opendap_idx <- grep("OPENDAP",ESGF_file_query$response$docs$url[[j]])
        http_idx <- grep("HTTPServer",ESGF_file_query$response$docs$url[[j]])
        http_url <- unlist(strsplit(ESGF_file_query$response$docs$url[[j]][http_idx],"|",fixed=TRUE))[1]
        #results[[paste('http',ic,jc,sep='_')]] <- http_url
        results[[paste('https',ic,jc,sep='_')]] <- http_url  ## REB 2021-04-12 - modification
        # if (verbose) print(http_url)
        opendap_url <- unlist(strsplit(ESGF_file_query$response$docs$url[[j]][opendap_idx],"|",fixed=TRUE))[1]
        # if (verbose) print(opendap_url)
        results[[paste('OpenDAP',ic,jc,sep='_')]] <- gsub("http://","https://",gsub(".nc.html",".nc",opendap_url))
        results[[paste('grid',ic,jc,sep='_')]] <- ESGF_file_query$response$docs$grid[[j]]
        results[[paste('member_id',ic,jc,sep='_')]] <- ESGF_file_query$response$docs[[j,facet_ens]]
        results[[paste('source_id',ic,jc,sep='_')]] <- ESGF_file_query$response$docs[[j,facet_mod]]
        results[[paste('type',ic,jc,sep='_')]] <- ESGF_file_query$response$docs$source_type[[j]]
        results[[paste('title',ic,jc,sep='_')]] <- ESGF_file_query$response$docs$title[[j]]
        results[[paste('timestamp',ic,jc,sep='_')]] <- ESGF_file_query$response$docs$timestamp[[j]]
        
        #Set facets that are not available in CMIP5 to NA
        if (is.null(results[[paste('grid',ic,jc,sep='_')]])) results[[paste('grid',ic,jc,sep='_')]] <- "NA"
        if (is.null(results[[paste('type',ic,jc,sep='_')]])) results[[paste('type',ic,jc,sep='_')]] <- "NA"
      }
      
      if (nof_hist_files == 0)
      {
        if (verbose) cat('\n','No historical dataset found for',ESGF_query$response$docs$id[i])
      } else  {
        for (j in 1:nof_hist_files) {
          jc <- as.character(j+nof_files); if (j+nof_files < 100) jc <- paste('0',jc,sep=''); if (j+nof_files < 10) jc <- paste('0',jc,sep='')
          
          results[[paste('file.query:',ic,jc,sep='_')]] <- ESGF_hist_file_query$response$docs[j,]
          
          opendap_idx <- grep("OPENDAP",ESGF_hist_file_query$response$docs$url[[j]])
          http_idx <- grep("HTTPServer",ESGF_hist_file_query$response$docs$url[[j]])
          http_url <- unlist(strsplit(ESGF_hist_file_query$response$docs$url[[j]][http_idx],"|",fixed=TRUE))[1]
          #results[[paste('http',ic,jc,sep='_')]] <- http_url
          results[[paste('https',ic,jc,sep='_')]] <- http_url  ## REB 2021-04-12 - modification
          # if (verbose) print(http_url)
          opendap_url <- unlist(strsplit(ESGF_hist_file_query$response$docs$url[[j]][opendap_idx],"|",fixed=TRUE))[1]
          # if (verbose) print(opendap_url)
          results[[paste('OpenDAP',ic,jc,sep='_')]] <- gsub("http://","https://",gsub(".nc.html",".nc",opendap_url))
          results[[paste('grid',ic,jc,sep='_')]] <- ESGF_hist_file_query$response$docs$grid[[j]]
          results[[paste('member_id',ic,jc,sep='_')]] <- ESGF_hist_file_query$response$docs[[j,facet_ens]]
          results[[paste('source_id',ic,jc,sep='_')]] <- ESGF_hist_file_query$response$docs[[j,facet_mod]]
          results[[paste('type',ic,jc,sep='_')]] <- ESGF_hist_file_query$response$docs$source_type[[j]]
          results[[paste('title',ic,jc,sep='_')]] <- ESGF_hist_file_query$response$docs$title[[j]]
          results[[paste('timestamp',ic,jc,sep='_')]] <- ESGF_hist_file_query$response$docs$timestamp[[j]]
          
          #Set facets that are not available in CMIP5 to NA
          if (is.null(results[[paste('grid',ic,jc,sep='_')]])) results[[paste('grid',ic,jc,sep='_')]] <- "NA"
          if (is.null(results[[paste('type',ic,jc,sep='_')]])) results[[paste('type',ic,jc,sep='_')]] <- "NA"
        }
      }
    }
    
    elements <- names(results)
    opendap <- grep('OpenDAP',names(results))
    #http <- grep('http',names(results))
    http <- grep('https',names(results))  ## REB 2021-04-12: modification
    mem <- grep('member_id',names(results))
    grid <- grep('grid',names(results))
    model <- grep('source_id',names(results))
    type <- grep('type',names(results))
    title <- grep('title',names(results))
    timestamp <- grep('timestamp',names(results))
    period <- substr(as.character(results[title]),nchar(as.character(results[title]))-15,
                     nchar(as.character(results[title]))-3)
    
    meta <- data.frame(OpenDap=as.character(results[opendap]),http=as.character(results[http]),
                       member.id=as.character(results[mem]),grid=as.character(results[grid]),
                       model=as.character(results[model]),type=as.character(results[type]),
                       title=as.character(results[title]), period=as.character(period),
                       timestamp=as.character(results[timestamp]))
    
    attr(meta,'variable') <- param
    attr(meta,'expid') <- expid
    attr(meta,'file.query.data') <- results[grep('file.query',names(results))]
    attr(meta,'dataset.query.data') <- results[grep('dataset.query',names(results))]
    attr(meta,'history') <- history.stamp(meta)
    
    return(meta)
  }
}
