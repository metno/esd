#' Retrieve CMIP data directly from the Earth System Grid Federation (ESGF)
#' https://earthsystemcog.org/projects/cog/esgf_search_restful_api
#' meta.ESGF returns a data.frame with the model metadata and the OpenDAP URL that 
#' can be used with retrieve.ESGF.

#' @aliases retrieve.ESGF meta.ESGF

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
  if (is.null(meta)) meta <- meta.ESGF()
  opendap <- as.character(meta$OpenDap[im])
  if (verbose) print(opendap)
  x <- try(retrieve(ncfile=opendap,param=attr(meta,'variable'),verbose=verbose,...))
  invisible(x)
}

#' @export meta.ESGF
meta.ESGF <- function(url="https://esgf-data.dkrz.de/esg-search/search/",mip="CMIP6",param="tas",
                      freq="mon",expid="ssp585",verbose=FALSE,n=NULL) {
  if (verbose) print('meta.ESFG - this function uses the rjson-package to read metadata from ESGF')
  ## Check if the JSON library is installed
  ## Check if you need to get the devtools-package:
  not.installed.rjson <- ("rjson" %in% rownames(installed.packages()) == FALSE)
  
  if (not.installed.rjson) {
    print('Need to install the rjson package')
    ## You need online access.
    print("e.g. install.packages('rjson')")
  }
  
  require(rjson)
  ## Get number of available datasets
  URL <- paste(url,"?type=Dataset&replica=false&latest=true&mip_era=",mip,"&variable_id=",param,
               "&frequency=",freq,"&experiment_id=",expid,"&format=application%2Fsolr%2Bjson",sep="")
  if (verbose) print(URL)

  nof_datasets <- fromJSON(file=URL)$response$numFound
  if (!is.null(n)) nof_datasets <- n
  print(paste('Found',nof_datasets,'datasets'))
  if (nof_datasets==0) {
    print('The connection seems to be down - try again later')
    return(NULL)
  }
  
  ## Query ESGF server
  ESGF_query <- fromJSON(file=paste(url,"?limit=",nof_datasets,"&type=Dataset&replica=false&latest=true&mip_era=",mip,"&variable_id=",param,
                                    "&frequency=",freq,"&experiment_id=,",expid,"&format=application%2Fsolr%2Bjson&facets=source_id",sep=""))
  if (verbose) str(ESGF_query)
  
  ## List files in each dataset
  #nof_datasets <- 7 ## test
  results <- list()
  for (i in 1:nof_datasets) {
    nof_files <- ESGF_query$response$docs[[i]]$number_of_files
    ESGF_file_query <- fromJSON(file=paste("https://esgf-data.dkrz.de/esg-search/search/?limit=",nof_files,"&type=File&format=application%2Fsolr%2Bjson&dataset_id=",ESGF_query$response$docs[[i]]$id,sep=""))
    
    if (verbose) { 
      print(paste(rep("=",nchar(ESGF_query$response$docs[[i]]$id)+9),collapse=""))
      print(paste("Dataset:",ESGF_query$response$docs[[i]]$id))
      print(paste(rep("=",nchar(ESGF_query$response$docs[[i]]$id)+9),collapse=""))
    } else cat('.')
    
    for (j in 1:nof_files) {
      ic <- as.character(i); if (i < 100) ic <- paste('0',ic,sep=''); if (i < 10) ic <- paste('0',ic,sep='')
      jc <- as.character(j); if (j < 100) jc <- paste('0',jc,sep=''); if (j < 10) jc <- paste('0',jc,sep='')
      results[[paste('file.query:',ic,jc,sep='_')]] <- ESGF_query$response$docs[[i]]
      url <- strsplit(ESGF_file_query$response$docs[[j]]$url[1],"|",fixed=TRUE)[[1]][1]
      results[[paste('http',ic,jc,sep='_')]] <- url
      if (verbose) print(url)
      url <- strsplit(ESGF_file_query$response$docs[[j]]$url[1],"|",fixed=TRUE)[[1]][1]
      results[[paste('OpenDAP',ic,jc,sep='_')]] <- gsub("fileServer","dodsC",url)
      #print(paste('OpenDAP',ic,jc,sep='_'))
      results[[paste('member_id',ic,jc,sep='_')]] <- ESGF_file_query$response$docs[[j]]$member_id
      results[[paste('grid',ic,jc,sep='_')]] <- ESGF_file_query$response$docs[[j]]$grid
      results[[paste('source_id',ic,jc,sep='_')]] <- ESGF_file_query$response$docs[[j]]$source_id
      results[[paste('type',ic,jc,sep='_')]] <- ESGF_file_query$response$docs[[j]]$source_type
      results[[paste('title',ic,jc,sep='_')]] <- ESGF_file_query$response$docs[[j]]$title
    }
  }
  elements <- names(results)
  opendap <- grep('OpenDAP',names(results))
  http <- grep('http',names(results))
  mem <- grep('member_id',names(results))
  grid <- grep('grid',names(results))
  model <- grep('source_id',names(results))
  type <- grep('type',names(results))
  title <- grep('title',names(results))
  period <- substr(as.character(results[title]),nchar(as.character(results[title]))-15,
                   nchar(as.character(results[title]))-3)
  meta <- data.frame(OpenDap=as.character(results[opendap]),http=as.character(results[http]),
                     member.id=as.character(results[mem]),grid=as.character(results[grid]),
                     model=as.character(results[model]),type=as.character(results[type]),
                     title=as.character(results[title]),period=as.character(period))
  attr(meta,'variable') <- param
  attr(meta,'all.query.data') <- results[grep('file.query',names(results))]
  attr(meta,'history') <- history.stamp(meta)
  return(meta)
  
}
