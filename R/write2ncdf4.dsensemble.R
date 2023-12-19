## ================================================
## Plan for the function write2ncdf4.dsensemble
##
## The input argument 'type' is used to redirect to different functions:
##
## type = c('dsensemble') -> ncdf4_dsensemble -> write PCA based dsensemble to a netCDF file, including an extraction method as an attribute
## type = c('ensemblestatistics', 'station') -> ncdf4_dsensemblestatistics -> convert PCA based dsensemble to station objects, calculate ensemble statistics, write to a netCDF file
## type = c('ensemblestatistics', 'field') -> ncdf4_dsensemblestatistics -> convert PCA based dsensemble to field objects, calculate ensemble statistics, write to a netCDF file
## type = c('model', 'field') -> ncdf4_model -> write downscaled results to netCDF files for individual models as gridded data
## type = c('model', 'station') -> ncdf4_model -> write downscaled results to netCDF files for individual models as station data
#
## the functions above return nc objects
##
## add_ncattributes is then used to add attributes
## Finally, the nc object is closed.
##
## Progress: 
##
## ncdf4_dsensemble and ncdf4_model exist and work (mostly)
##
## ncdf4_dsensemblestatistics is empty
##
## ncdf_model works for writing data as gridded fields, but not for station data.
##
## ================================================

#' Saves climate data as netCDF.
#' 
#' Method to save station data as netCDF, making sure to include the data
#' structure and meta-data (attributes). The code tries to follow the netCDF
#' 'CF' convention. The method is built on the \code{ncdf4} package.
#'
#' To save space, the values are saved as short (16-bit signed integer that
#' can hold values between -32768 and 32767).
#' (see NC_SHORT in \url{https://www.unidata.ucar.edu/software/netcdf/docs/data_type.html}).
#'
#' @seealso write2ncdf4
#' 
#' @param x a 'dsensemble' object containing empirically-statistically downscaled results (output from \code{'DSensemble'})
#' @param file name of output file 
#' @param path path to output file
#' @param force If TRUE (default), overwrite existing files. If FALSE, do not write data to a file if there is already one with the same file name.  
#' @param prec Precision: see \code{\link[ncdf4]{ncvar_def}}
#' @param scale Sets the attribute 'scale_factor' which is used to scale
#' (multiply) the values stored (to save space may be represented as 'short').
#' @param im a model index, see \code{\link{subset}}
#' @param im a numerical or numerical vector with indices of the ensemble members to be included. If NULL include all.
#' @param type a vector defining the type of data to write to the netCDF file. 
#' The default type is "dsensemble" which writes dsensemble data in its original form.
#' Other options are "ensemblestatistics" to calculate and write ensemble statistics,
#' or "model" which writes downscaled results for a single model simulation   
#' (the ensemble member is then selected with the input argument \code{'im'}).
#' The options 'ensemblestatistics' and 'model' can be used in conjunction with the additional argument, 
#' e.g. c('model','field') or c('model','station') to save the downscaled results as a field (gridded data) 
#' or a 'station' object (the stations that were used as predictand in the downscaling).
#' @param offset Sets the attribute 'add_offset' which is added to the values
#' stored (to save space may be represented as 'short').
#' @param torg Time origin
#' @param missval Missing value: see \code{\link[ncdf4]{ncvar_def}}
#' @param region Name of the region of the downscaled data, if not NULL (which is the default) used in the filename that is generated if \code{file} is not specified.
#' @param ensemblename Name of the ensemble - goes into the filename that is generated if \code{file} is not specified. The default is NULL and the ensemble name does not have to be set unless a subset of models are selected (e.g., with the input \code{im}).
#' @param id An identifier for the data set, provided by and unique within its naming authority. The combination of the "naming authority" and the "id" should be globally unique, but the id can be globally unique by itself also. IDs can be URLs, URNs, DOIs, meaningful text strings, a local key, or any other unique string of characters. The id should not include white space characters.
#' @param naming_authority	The organization that provides the initial id (see above) for the dataset. The naming authority should be uniquely specified by this attribute. We recommend using reverse-DNS naming for the naming authority; URIs are also acceptable. Example: 'edu.ucar.unidata'.
#' @param history	Provides an audit trail for modifications to the original data. This attribute is also in the NetCDF Users Guide: 'This is a character array with a line for each invocation of a program that has modified the dataset. Well-behaved generic netCDF applications should append a line containing: date, time of day, user name, program name and command arguments.' To include a more complete description you can append a reference to an ISO Lineage entity; see NOAA EDM ISO Lineage guidance.
#' @param source	The method of production of the original data. If it was model-generated, source should name the model and its version. If it is observational, source should characterize it. This attribute is defined in the CF Conventions. Examples: 'temperature from CTD #1234'; 'world model v.0.1'.
#' @param processing_level	A textual description of the processing (or quality control) level of the data.
#' @param comment	Miscellaneous information about the data, not captured elsewhere. This attribute is defined in the CF Conventions.
#' @param acknowledgement	A place to acknowledge various types of support for the project that produced this data.
#' @param license	Provide the URL to a standard or specific license, enter "Freely Distributed" or "None", or describe any restrictions to data access and distribution in free text.
#' @param standard_name_vocabulary	The name and version of the controlled vocabulary from which variable standard names are taken. (Values for any standard_name attribute must come from the CF Standard Names vocabulary for the data file or product to comply with CF.) Example: 'CF Standard Name Table v27'.
#' @param date_created	The date on which this version of the data was created. (Modification of values implies a new version, hence this would be assigned the date of the most recent values modification.) Metadata changes are not considered when assigning the date_created. The ISO 8601:2004 extended date format is recommended, as described in the Attribute Content Guidance section.
#' @param creator_name	The name of the person (or other creator type specified by the creator_type attribute) principally responsible for creating this data.
#' @param creator_email	The email address of the person (or other creator type specified by the creator_type attribute) principally responsible for creating this data.
#' @param creator_url	The URL of the person (or other creator type specified by the creator_type attribute) principally responsible for creating this data.
#' @param institution	The name of the institution principally responsible for originating this data. This attribute is recommended by the CF convention.
#' @param project	The name of the project(s) principally responsible for originating this data. Multiple projects can be separated by commas, as described under Attribute Content Guidelines. Examples: 'PATMOS-X', 'Extended Continental Shelf Project'.
#' @param publisher_name	The name of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.
#' @param publisher_email	The email address of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.
#' @param publisher_url	The URL of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.
#' @param creator_type	Specifies type of creator with one of the following: 'person', 'group', 'institution', or 'position'. If this attribute is not specified, the creator is assumed to be a person.
#' @param creator_institution	The institution of the creator; should uniquely identify the creator's institution. This attribute's value should be specified even if it matches the value of publisher_institution, or if creator_type is institution.
#' @param publisher_type	Specifies type of publisher with one of the following: 'person', 'group', 'institution', or 'position'. If this attribute is not specified, the publisher is assumed to be a person.
#' @param publisher_institution	The institution that presented the data file or equivalent product to users; should uniquely identify the institution. If publisher_type is institution, this should have the same value as publisher_name.
#' @param product_version	Version identifier of the data file or product as assigned by the data creator. For example, a new algorithm or methodology could result in a new product_version.
#' @param summary Summary
#' @param summary_no Summary in Norwegian
#' @param keywords Keywords
#' @param keywords_vocabulary Vocabulary of keywords. Default: GCMDSK:GCMD Science Keywords:https://gcmd.earthdata.nasa.gov/kms/concepts/concept_scheme/sciencekeywords, GCMDPROV:GCMD
#' @param references	Published or web-based references that describe the data or methods used to produce it. Recommend URIs (such as a URL or DOI) for papers or other references. This attribute is defined in the CF conventions.
#' @param title Title of data set
#' @param title_no Title of data set in Norwegian
#' @param additional_attributes List of attributes to add to the netCDF file. The default is NULL 
#' but it could be something like list("distribution_statement"="Free", "title"="Temperature projections for Norway")
#' @param verbose TRUE If TRUE, print information on progress.
#' @param \dots additional arguments
#' 
#' @return None
#' 
#' @keywords netcdf ncdf4 save
#' 
#' @exportS3Method
#' @export write2ncdf4.dsensemble
write2ncdf4.dsensemble <- function(x,...,file=NULL,path=NULL,force=TRUE,
                                   prec='short',offset=0,scale=0.1 ,
                                   type='dsensemble',im=NULL,is=NULL,it=NULL,
                                   torg="1970-01-01",missval=-999,
                                   method="metnoESD",region=NULL,ensemblename=NULL,
                                   conventions="ACDD-1.3",
                                   id=NA,naming_authority=NA,
                                   source=NA,processing_level="Scientific",
                                   comment=NA,acknowledgement=NA,
                                   license="Freely distributed",
                                   standard_name_vocabulary="CF Standard Name Table v27",
                                   creator_name=NA,creator_email=NA,creator_url=NA,
                                   creator_institution=NA,creator_type=NA,
                                   institution=NA,institution_short_name=NA,
                                   project=NA,project_short_name=NA,
                                   publisher_name=NA,publisher_email=NA,publisher_url=NA,
                                   publisher_type=NA,publisher_institution=NA,
                                   product_version=NA,summary=NA,keywords=NA,
                                   keywords_vocabulary="GCMDSK:GCMD Science Keywords:https://gcmd.earthdata.nasa.gov/kms/concepts/concept_scheme/sciencekeywords, GCMDPROV:GCMD",
                                   references=NA,title=NA,title_no=NA,
                                   additional_attributes=NA,verbose=FALSE) {
  
  if(verbose) print('write2ncdf4.dsensemble')
  if("model" %in% type & length(im)>1) {
    if(verbose) print(paste("Writing",length(im),"netCDF files."))
    if(length(im)>1) file <- NULL
    for(i in im) {
      write2ncdf4(x,...,file=file,path=path,
                  im=i,is=is,it=it,
                  prec=prec,scale=scale,offset=offset,
                  type=type,method=method,region=region,ensemblename=ensemblename,
                  torg=torg,missval=missval,
                  conventions=conventions,
                  id=id,naming_authority=naming_authority,
                  source=source,processing_level=processing_level,
                  comment=comment,acknowledgement=acknowledgement,license=license,
                  standard_name_vocabulary=standard_name_vocabulary,
                  creator_name=creator_name,creator_email=creator_email,creator_url=creator_url,
                  creator_institution=creator_institution,creator_type=creator_type,
                  institution=institution,institution_short_name=institution_short_name,
                  project=project,project_short_name=project_short_name,
                  publisher_name=publisher_name,publisher_email=publisher_email,publisher_url=publisher_url,
                  publisher_type=publisher_type,publisher_institution=publisher_institution,
                  product_version=product_version,summary=summary,keywords=keywords,
                  references=references,title=title,title_no=title_no,
                  additional_attributes=additional_attributes,
                  force=force,verbose=verbose)
    } 
  } else {
    if("dsensemble" %in% type) {
      file <- ncdf4_dsensemble(x,...,file=file,path=path,prec=prec,
                             im=im,is=is,it=it,
                             scale=scale,type=type,
                             offset=offset,torg=torg,missval=missval,
                             method=method,region=region,ensemblename=ensemblename,
                             force=force,verbose=verbose)
    } else if ("ensemblestatistics" %in% type) {
      file <- ncdf4_dsensemblestatistics(x,...,file=file,path=path,prec=prec,
                                       im=im,is=is,it=it,
                                       scale=scale,type=type,
                                       offset=offset,torg=torg,missval=missval,
                                       method=method,region=region,ensemblename=ensemblename,
                                       force=force,verbose=verbose)
    } else if("model" %in% type) {
      file <- ncdf4_dsmodel(x,...,file=file,path=path,
                          im=im,is=is,it=it,
                          prec=prec,scale=scale,type=type,
                          offset=offset,torg=torg,missval=missval,
                          method=method,region=region,
                                       force=force,verbose=verbose)
    } else file <- NULL
   
    if(!is.null(file)) if(file.exists(file)) {
      nc <- try(nc_open(file, write=TRUE))
      if(inherits(nc, "ncdf4")) {
        nc <- add_ncattributes(nc, conventions=conventions,id=id,naming_authority=naming_authority,
                         source=source,processing_level=processing_level,
                         comment=comment,acknowledgement=acknowledgement,
                         license=license,standard_name_vocabulary=standard_name_vocabulary,
                         creator_name=creator_name,creator_email=creator_email,creator_url=creator_url,
                         creator_institution=creator_institution,creator_type=creator_type,
                         institution=institution,institution_short_name=institution_short_name,
                         project=project,project_short_name=project_short_name,
                         publisher_name=publisher_name,publisher_email=publisher_email,publisher_url=publisher_url,
                         publisher_type=publisher_type,publisher_institution=publisher_institution,
                         product_version=product_version,summary=summary,keywords=keywords,
                         references=references,title=title,title_no=title_no,
                         additional_attributes=additional_attributes)
        nc_close(nc)
      }
      if (verbose) print(paste('Successfully finished writing dsensemble data and attributes to file', file))
    }
  }
  if (verbose) print(paste('exit write2ncdf4.dsensemble'))
  gc(reset=TRUE)
}

add_ncattributes <- function(nc, conventions="ACDD-1.3",
                             id=NA,naming_authority=NA,
                             source=NA,processing_level="Scientific",
                             comment=NA,acknowledgement=NA,
                             license="Freely distributed",
                             standard_name_vocabulary="CF Standard Name Table v27",
                             creator_name=NA,creator_email=NA,creator_url=NA,
                             creator_institution=NA,creator_type=NA,
                             institution=NA,institution_short_name=NA,
                             project=NA,project_short_name=NA,
                             publisher_name=NA,publisher_email=NA,publisher_url=NA,
                             publisher_type=NA,publisher_institution=NA,
                             product_version=NA,summary=NA,keywords=NA,
                             references=NA,title=NA,title_no=NA,
                             additional_attributes=NA,verbose=FALSE) {
  if(verbose) print("add_ncattributes")
  stopifnot("The input nc must be a 'ncdf4' object"=inherits(nc, "ncdf4"))
  if(verbose) print("Adding attributes from input arguments.")
  ncatt_put( nc, 0, "title", as.character(title), prec="text")
  ncatt_put( nc, 0, "title_no", as.character(title_no), prec="text")
  ncatt_put( nc, 0, "summary", as.character(summary), prec="text")
  ncatt_put( nc, 0, "keywords", as.character(keywords), prec="text")
  ncatt_put( nc, 0, "conventions", as.character(conventions), prec="text")
  ncatt_put( nc, 0, "id", as.character(id), prec="text")
  ncatt_put( nc, 0, "naming_authority", as.character(naming_authority), prec="text")
  ncatt_put( nc, 0, "comment", as.character(comment), prec="text")
  ncatt_put( nc, 0, "acknowledgement", as.character(acknowledgement), prec="text")
  ncatt_put( nc, 0, "license", as.character(license), prec="text")
  ncatt_put( nc, 0, "standard_name_vocabulary", as.character(standard_name_vocabulary), prec="text")
  ncatt_put( nc, 0, "creator_name", as.character(creator_name), prec="text")
  ncatt_put( nc, 0, "creator_email", as.character(creator_email), prec="text")
  ncatt_put( nc, 0, "creator_institution", as.character(creator_institution), prec="text")
  ncatt_put( nc, 0, "creator_type", as.character(creator_type), prec="text")
  ncatt_put( nc, 0, "institution", as.character(institution), prec="text")
  ncatt_put( nc, 0, "institution_short_name", as.character(institution_short_name), prec="text")
  ncatt_put( nc, 0, "project", as.character(project), prec="text")
  ncatt_put( nc, 0, "project_short_name", as.character(project_short_name), prec="text")
  ncatt_put( nc, 0, "publisher_name", as.character(publisher_name), prec="text")
  ncatt_put( nc, 0, "publisher_email", as.character(publisher_email), prec="text")
  ncatt_put( nc, 0, "publisher_url", as.character(publisher_url), prec="text")
  ncatt_put( nc, 0, "publisher_institution", as.character(publisher_institution), prec="text")
  ncatt_put( nc, 0, "publisher_type", as.character(publisher_type), prec="text")
  ncatt_put( nc, 0, "product_version", as.character(product_version), prec="text")
  ncatt_put( nc, 0, "references", as.character(references), prec="text")
  
  ## Add additional attributes from the input list 'additional_attributes'
  if(is.list(additional_attributes) & length(additional_attributes)>0) {
    for(a in names(additional_attributes)) {
      ncatt_put( nc, 0, a, as.character(additional_attributes[[a]]), prec="text")
    }
  }
  
  ## Return nc object with added attributes
  if (verbose) print(paste('Finished writing attributes to file', file))
  if (verbose) print(paste('exit add_ncattributes'))
  return(nc)
}

ncdf4_dsensemble <- function(x,...,file=NULL,path=NULL,force=TRUE,
                             prec='short',offset=0,scale=0.1,
                             torg="1970-01-01",missval=-99,
                             im=NULL, is=NULL, it=NULL,
                             method="metnoESD",region=NULL,ensemblename=NULL,
                             verbose=TRUE) {
  if (verbose) print('ncdf4_dsensemble')
  if(!is.null(im) | !is.null(is) | !is.null(it)) x <- subset(x, im=im, is=is, it=it)
  
  ## Put the PCA and EOF patterns into a list
  x0 <- x
  info <- x$info
  pattern <- list()
  if(!is.null(x$pca)) pattern$pca <- x$pca
  if(!is.null(x$eof)) {
    pattern$eof <- x$eof
    res <- paste0(round(diff(lon(x$eof))[1], digits=2),"x",
                  round(diff(lat(x$eof))[1], digits=2),"deg")
  } else res <- "stations"
  x$pca <- NULL
  x$info <- NULL
  x$eof <- NULL
  
  if(is.null(file)) {
    file <- "dsensemble" 
    file <- paste0(file,"_",attr(x,"scenario"))
    if(!is.null(method)) file <- paste0(file,"_",method)
    if(!is.null(region)) file <- paste0(file,"_",region)
    file <- paste0(file,"_",res)
    file <- paste0(file,"_",attr(x,"variable")[1])
    if(inherits(x0[[2]],"season")) {
      file <- paste0(file,"_sem")
      file <- paste0(file,"_",toupper(season(x0[[2]])[1]))
    } else if(inherits(x0[[2]],"annual")) {
      file <- paste0(file,"_annual-mean")
    } else if(inherits(x0[[2]],"month")) {
      file <- paste0(file,"_monthly-mean")
    }
    file <- paste0(file,"_",paste(range(year(x[[1]])),collapse="-"))
    file <- paste0(file,"_ngcm",length(attr(x, "model_id")))
    if(!is.null(ensemblename)) file <- paste0(file, "_", ensemblename)
    file <- paste0(file, ".nc")
  }
  if(!is.null(path)) file <- file.path(path, basename(file))
  
  if(file.exists(file) & !force) {
    nc <- NULL
    if (verbose) print(paste('File',file,'already exists. No new file created.'))
  } else {
    ## Get the two first elements of the list which contain information common to the rest
    if(inherits(x,"pca")) {
      ## Organise the time index
      tx <- index(x[[1]])
      if(inherits(tx,"numeric")) {
        yr <- seq(min(tx), max(tx))
        mn <- month(attr(x[[1]], "evaluation"))[1]
        tx <- as.Date(paste(yr,mn,"01",sep="-"), torg=torg)
      }
      yr <- year(tx)
      
      ## Put the eigenvectors into an array Y with dimensions [PC mode, model simulation, time step]
      Y <- array(NA, dim=c(max(sapply(x,ncol)), length(x), length(tx)))
      for(j in seq(1,length(x))) {
        i <- 1:ncol(x[[j]])
        k <- which(yr %in% year(x[[j]]))
        l <- which(year(x[[j]]) %in% yr)
        Y[i,j,k] <- x[[j]][l]
      }
      d <- dim(Y)
      offset <- mean(Y,na.rm=TRUE)
      scale <- round2magnitude((max(abs(c(Y)),na.rm=TRUE) - offset)/10000)
      Y <- round((Y-offset)/scale)

      ## Create dimensions and define a variable for the eigenvectors
      dimpc <- ncdim_def("mode", "index", seq(dim(Y)[1]),
                         longname="Mode of variability")
      dimim <- ncdim_def("ensemble_member", "index", seq(dim(Y)[2]),
                          longname="Ensemble member index")
      dimtim <- ncdim_def("time", paste("days since",torg),
                          as.numeric(as.Date(tx,origin=torg)) )
      x4nc <- ncvar_def("eigenvectors", "unitless", list(dimpc,dimim,dimtim), -1, 
                        longname=paste("Eigenvectors of the",
                                       gsub("_"," ",attr(x0,'longname'))), prec=prec)
      
      ## Add eigenvalues as a variable
      nc <- nc_create( file, x4nc )
      ncvar_put( nc, "eigenvectors", Y)
      ncatt_put( nc, "eigenvectors", "add_offset", offset, prec="float" )
      ncatt_put( nc, "eigenvectors", "scale_factor", scale, prec="float" ) 
      ncatt_put( nc, "eigenvectors", "_FillValue", missval, prec="float" ) 
      ncatt_put( nc, "eigenvectors", "missing_value", missval, prec="float" )
      ncatt_put( nc, "eigenvectors", "time_coverage_start", as.character(min(tx)), prec="text")
      ncatt_put( nc, "eigenvectors", "time_coverage_end", as.character(max(tx)), prec="text")
      if(!is.null(attr(x,"model_id"))) ncatt_put( nc, "eigenvectors", "ensemble_members", 
                                                  attr(x,"model_id"), prec="text" ) 
      nc_close(nc)
      
      ## Add eigenvalues as a separate variable
      if(!is.null(pattern$pca)) {
        eig <- attr(pattern$pca, "eigenvalues") 
      } else if (!is.null(pattern$eof)) {
        eig <- attr(pattern$eof, "eigenvalues") 
      } else eig <- NULL
      if(!is.null(eig)) {
        offset <- 0
        scale <- round2magnitude((max(abs(c(eig)),na.rm=TRUE) - offset)/10000)
        eig <- round((eig-offset)/scale)
        dimpc <- ncdim_def("mode", "index", seq(length(eig)),
                           longname="Mode of variability")
        var_new <- ncvar_def("eigenvalues", "unitless", list(dimpc), 
                             missval, prec=prec,
                             longname=paste("Eigenvalues of the",
                                            gsub("_"," ",attr(x,'longname')[1])))
        nc <- nc_open(file, write=TRUE)
        nc <- ncvar_add( nc, var_new )
        nc_close(nc)
        nc <- nc_open(file, write=TRUE)
        ncvar_put( nc, "eigenvalues", eig)
        ncatt_put( nc, "eigenvalues", "add_offset", offset, prec="float" )
        ncatt_put( nc, "eigenvalues", "scale_factor", scale, prec="float" ) 
        ncatt_put( nc, "eigenvalues", "_FillValue", missval, prec="float" ) 
        ncatt_put( nc, "eigenvalues", "missing_value", missval, prec="float" )
      }
      
      ## Add PCA pattern as a separate variable
      if(!is.null(pattern$pca)) {
        pca <- attr(pattern$pca, "pattern")
        offset <- 0
        scale <- round2magnitude((max(abs(c(pca)),na.rm=TRUE) - offset)/10000)
        pca <- round((pca-offset)/scale)
        lons <- attr(pattern$pca, "longitude")
        lats <- attr(pattern$pca, "latitude")
        alts <- attr(pattern$pca, "altitude")
        stid <- as.numeric(attr(pattern$pca, "station_id"))
        locs <- attr(pattern$pca, "location")
        dimst <- ncdim_def("stid", "id", stid, longname="Station id")
        dimpc <- ncdim_def("mode", "index", seq(ncol(pca)),
                           longname="Mode of variability")
        var_new <- ncvar_def("PCA_pattern", "unitless", list(dimst,dimpc), missval, 
                             longname=paste("PCA pattern of the",
                                            gsub("_"," ",attr(x,'longname')[1])),
                             prec=prec)
        nc <- nc_open(file, write=TRUE)
        nc <- ncvar_add( nc, var_new )
        nc_close(nc)
        nc <- nc_open(file, write=TRUE)
        ncvar_put( nc, "PCA_pattern", pca)
        ncatt_put( nc, "PCA_pattern", "add_offset", offset, prec="float" )
        ncatt_put( nc, "PCA_pattern", "scale_factor", scale, prec="float" ) 
        ncatt_put( nc, "PCA_pattern", "_FillValue", missval, prec="float" ) 
        ncatt_put( nc, "PCA_pattern", "missing_value", missval, prec="float" )
        ncatt_put( nc, "PCA_pattern", "station_id", stid, prec="float" )
        ncatt_put( nc, "PCA_pattern", "location", locs, prec="text" )
        ncatt_put( nc, "PCA_pattern", "longitude", lons, prec="float" )
        ncatt_put( nc, "PCA_pattern", "longitude_unit", "degrees E", prec="text" )
        ncatt_put( nc, "PCA_pattern", "latitude", lats, prec="float" )
        ncatt_put( nc, "PCA_pattern", "latitude_unit", "degrees N", prec="text" )
        ncatt_put( nc, "PCA_pattern", "altitude", alts, prec="float" )
        ncatt_put( nc, "PCA_pattern", "altitude_unit", "masl", prec="text" )
        ncatt_put( nc, "PCA_pattern", "spatial_representation", "vector", prec="text")
        nc_close(nc)
      }
      
      ## Add PCA pattern as a separate variable
      if(!is.null(pattern$eof)) {
        eof <- attr(pattern$eof, "pattern")
        offset <- 0
        scale <- round2magnitude((max(abs(c(eof)),na.rm=TRUE) - offset)/10000)
        eof <- round((eof-offset)/scale)
        lons <- attr(pattern$eof, "longitude")
        lats <- attr(pattern$eof, "latitude")
        dimlon <- ncdim_def( "longitude", "degree_east", lons )
        dimlat <- ncdim_def( "latitude", "degree_north", lats )
        dimpc <- ncdim_def("mode", "index", seq(dim(eof)[3]),
                           longname="Mode of variability")
        var_new <- ncvar_def("EOF_pattern", "unitless", list(dimlon, dimlat, dimpc), 
                             missval, prec=prec,
                             longname=paste("EOF pattern of the",
                                            gsub("_"," ",attr(x,'longname')[1])))
        nc <- nc_open(file, write=TRUE)
        nc <- ncvar_add( nc, var_new )
        nc_close(nc)
        nc <- nc_open(file, write=TRUE)
        ncvar_put( nc, "EOF_pattern", eof)
        ncatt_put( nc, "EOF_pattern", "add_offset", offset, prec="float" )
        ncatt_put( nc, "EOF_pattern", "scale_factor", scale, prec="float" ) 
        ncatt_put( nc, "EOF_pattern", "_FillValue", missval, prec="float" ) 
        ncatt_put( nc, "EOF_pattern", "missing_value", missval, prec="float" )
        ncatt_put( nc, "EOF_pattern", "geospatial_lon_min", as.character(round(min(lons))), prec="text")
        ncatt_put( nc, "EOF_pattern", "geospatial_lon_max", as.character(round(max(lons))), prec="text")
        ncatt_put( nc, "EOF_pattern", "geospatial_lon_units", "degrees E", prec="text")
        ncatt_put( nc, "EOF_pattern", "geospatial_lon_resolution", as.character(diff(lons)[[1]]), prec="text")
        ncatt_put( nc, "EOF_pattern", "geospatial_lat_min", as.character(round(min(lats))), prec="text")
        ncatt_put( nc, "EOF_pattern", "geospatial_lat_max", as.character(round(max(lats))), prec="text")
        ncatt_put( nc, "EOF_pattern", "geospatial_lat_units", "degrees N", prec="text")
        ncatt_put( nc, "EOF_pattern", "geospatial_lat_resolution", as.character(diff(lats)[[1]]), prec="text")
        ncatt_put( nc, "EOF_pattern", "spatial_representation", "grid", prec="text")
        nc_close(nc)
      }
      
      if(!is.null(attr(x, "r.xval"))) {
        rx <- attr(x,"r.xval")
        offset <- 0
	scale <- 1
        rx <- round((rx-offset)/scale)
        dimpc <- ncdim_def("mode", "index", seq(ncol(rx)),
                           longname="Mode of variability")
        dimim <- ncdim_def("ensemble_member", "index", seq(nrow(rx)),
                            longname="Ensemble member index")
        
        var_new <- ncvar_def("r_xval", "unitless", list(dimim,dimpc), missval, 
                             longname="Cross-validation correlation coefficient",
                             prec=prec)
        nc <- nc_open(file, write=TRUE)
        nc <- ncvar_add( nc, var_new )
        nc_close(nc)
        nc <- nc_open(file, write=TRUE)
        ncvar_put( nc, "r_xval", rx)
        ncatt_put( nc, "r_xval", "add_offset", offset, prec="float" )
        ncatt_put( nc, "r_xval", "scale_factor", scale, prec="float" ) 
        ncatt_put( nc, "r_xval", "_FillValue", missval, prec="float" ) 
        ncatt_put( nc, "r_xval", "missing_value", missval, prec="float" )
        nc_close(nc)
      }
    } else {
      ## What should we do when applying this function to a non PCA-based DSensemble object?
      print(paste("write2netcdf4.dsensemble cannot yet be applied to this type of input object.",
                  "The function has not been adapted to write objects of class",
                  paste(class(x), sep=", ")))
      browser()
    }
    
    ## Add global attributes
    if(file.exists(file)) {
      nc <- nc_open(file, write=TRUE)
      ncatt_put( nc, 0, "description", 
                 paste("Saved from esd on",date(),"using write2ncdf4.dsensemble to",
                       "write DSensemble data for all ensemble members in its original format."))
      ncatt_put( nc, 0, "esd-version", attr(x,'history')$session$esd.version)
      ncatt_put( nc, 0, "date_created", as.character(as.Date(Sys.time())), prec="text")
      ncatt_put( nc, 0, 'class', class(x))
      ncatt_put( nc, 0, "variable", attr(x,"variable")[1], prec="text")
      ncatt_put( nc, 0, "unit", attr(x,"unit")[1], prec="text")
      if(!is.null(attr(x,"model_id"))) ncatt_put( nc, 0, "ensemble_members", attr(x, "model_id"), prec="text" )
      if(!is.null(attr(x,"domain"))) ncatt_put( nc, 0, "predictor_domain", 
                                                paste0(paste(attr(x,"domain")$lon, collapse="-"), "E/",
                                                       paste(attr(x,"domain")$lat, collapse="-"), "N"), prec="text")
      if(!is.null(ensemblename)) ncatt_put( nc, 0, "ensemble_name", ensemblename, prec="text" )
      
      # Add the attributes of object x
      attnames <- names(attributes(x))
      if (verbose) print(attnames)
      exclude <- c('history','units','variable','dim','index','longitude',
                   'latitude','greenwich','domain','call','original_data',
                   'r.xval')
      for(a in exclude) attnames <- attnames[!grepl(a,attnames)]
      for (ia in seq_along(attnames)) {
        if(!inherits(attr(x,attnames[ia]),c("matrix","array","data.frame","list"))) {
          if (verbose) print(paste(attnames[ia], paste(attr(x,attnames[ia]), collapse=", ")))
          ncatt_put( nc, 0, attnames[ia], as.character(attr(x,attnames[ia])), prec="text")
        }
      }
      nc_close(nc)
    }
    if (verbose) print(paste('Finished writing dsensemble to file', file))
  }
  if (verbose) print(paste('exit ncdf4_dsensemble'))
  return(file)
}

ncdf4_dsensemblestatistics <- function(x,...,file=NULL,path=NULL,force=TRUE,
                                  prec='short',offset=0,scale=0.01,qc=TRUE,
                                  method="metnoESD",region=NULL,ensemblename=NULL,
                                  FUNX=c("mean","median","max","min","q95","q5","sd"),
                                  im=NULL, is=NULL, it=NULL, eof=TRUE,
                                  torg="1970-01-01",missval=-99,verbose=TRUE) {
  if (verbose) print('ncdf4_dsensemblestatistics')
  if(verbose) print("Function for writing ensemble statistics of downscaled results (output from DSensemble) to a netCDF file.")
  stopifnot(inherits(x,"dsensemble") & inherits(x,"list"))
  if(!is.null(im) | !is.null(is) | !is.null(it)) x <- subset(x, is=is, it=it)
  
  if(is.null(file)) {
    file <- "dsensemblestatistics"
    file <- paste0(file,"_",attr(x,"scenario"))
    if(!is.null(method)) file <- paste0(file,"_",method)
    if(!is.null(region)) file <- paste0(file,"_",region)
    if(!is.null(x$eof)) {
      res <- paste0(round(diff(lon(x$eof))[1], digits=2),"x",
                    round(diff(lat(x$eof))[1], digits=2),"deg")
    } else res <- "stations"
    file <- paste0(file,"_",res)
    file <- paste0(file,"_",attr(x,"variable")[1])
    if(inherits(x[[2]],"season")) {
      file <- paste0(file,"_sem")
      file <- paste0(file,"_",toupper(season(x[[2]])[1]))
    } else if(inherits(x[[2]],"annual")) {
      file <- paste0(file,"_annual-mean")
    } else if(inherits(x[[2]],"month")) {
      file <- paste0(file,"_monthly-mean")
    }
    file <- paste0(file,"_",paste(range(year(x[[3]])),collapse="-"))
    if(is.null(im)) ngcm <- length(attr(x,"model_id")) else if(is.logical(im)) ngcm <- sum(im) else ngcm <- length(im)
    file <- paste0(file,"_ngcm",ngcm)
    if(!is.null(ensemblename)) file <- paste0(file, "_", ensemblename)
    file <- paste0(file, ".nc")
  }
  if(!is.null(path)) file <- file.path(path, basename(file))
  
  if(file.exists(file) & !force) {
    nc <- NULL
    if (verbose) print(paste('File',file,'already exists. No new file created.'))
  } else {
    if(verbose) print("Calculate ensemble statistics.")
    ## KMP 2023-05-23: Moved the selection of models ('im') from subset (on line 536) to aggregate (below) 
    ## so that same ensemble member can be subset multiple times calculating weighted ensemble statistics
    xstats <- aggregate(x, FUNX=FUNX, eof=eof, im=im, qc=qc, verbose=verbose)
    if(length(FUNX)==1 & !is.list(xstats)) {
      xf <- xstats
      xstats <- list()
      xstats[[paste0(FUNX)]] <- xf
    }
    ## Generate time index
    if(is.list(xstats)) tx <- index(xstats[[1]]) else tx <- index(xstats)
    if (inherits(tx,c('numeric','integer'))) {
      yr <- tx
      if(!is.null(x$pca)) {
        mn <- switch(season(x$pca)[[1]], "djf"="01", "mam"="03", 
                     "jja"="06", "son"="10", "01") 
      } else {
        mn <- switch(season(x)[[1]], "djf"="01", "mam"="03", 
                     "jja"="06", "son"="10", "01")
      }
      tx <- as.Date(paste(yr,'-',mn,'-01',sep=''))
    }

    ## Replace missing values, organize dimensions of data, and subtract offset and divide by scale 
    if (is.null(offset)) offset <- mean(xstats$mean,na.rm=TRUE)
    if (is.null(scale)) scale <- (max(abs(c(xstats$mean)),na.rm=TRUE) - offset)/10000
    for(fn in names(xstats)) {
      xf <- xstats[[fn]]
      xf <- t(xf)
      if(eof) {
        d <- c(dim(attr(x$eof,"pattern"))[1:2], dim(xf)[2])
      } else {
        d <- dim(xf)
      }
      xf <- as.vector(xf)
      xf <- round((xf-offset)/scale)
      xf[!is.finite(xf)] <- missval
      dim(xf) <- d
      xstats[[fn]] <- xf
    }
    
    if (verbose) {
      print(attr(xstats[[1]],'dimensions')); print(c(scale,offset))
      print(range(c(xstats[[1]]))); print(range(c(xstats[[1]]),na.rm=TRUE))
    }

    ## Set dimensions for the netCDF file
    dimlon <- ncdim_def( "longitude", "degree_east", lon(xstats) )
    dimlat <- ncdim_def( "latitude", "degree_north", lat(xstats) )
    dimtim <- ncdim_def( "time", paste("days since",torg),
                         as.numeric(as.Date(tx,origin=torg)) )
    
    for(i in seq_along(names(xstats))) {
      fn <- names(xstats)[i]
      varname <- paste0("ens",fn)
      x4nc <- ncvar_def(varname, attr(xstats,"unit")[1], list(dimlon,dimlat,dimtim), -1, 
                        longname=paste("ensemble",fn,"of the",
                                     gsub("_"," ",attr(xstats,'longname'))), 
                                     prec=prec)
      # Create a netCDF file with this variable or put it in an existing file
      if(i>1) {
        nc <- nc_open(file, write=TRUE)
        nc <- ncvar_add(nc, x4nc)
        nc_close(nc)
        nc <- nc_open(file, write=TRUE)
      } else {
        nc <- nc_create( file, x4nc )
        ncatt_put( nc, 0, "description", 
                   paste("Saved from esd on",date(),"using write2ncdf4.dsensemble to",
                         "write DSensemble data as ensemble statistics."))
        #if (verbose) print(attr(x,'history'))
        ncatt_put( nc, 0, "esd-version", attr(x,'history')$session$esd.version)
        # Add standard attributes
        ncatt_put( nc, 0, "date_created", as.character(as.Date(Sys.time())), prec="text")
        ncatt_put( nc, 0, "units", attr(x,"unit")[[1]], prec="text")
        ncatt_put( nc, 0, "geospatial_lon_min", as.character(round(min(lon(xstats)))), prec="text")
        ncatt_put( nc, 0, "geospatial_lon_max", as.character(round(max(lon(xstats)))), prec="text")
        ncatt_put( nc, 0, "geospatial_lon_units", "degrees E", prec="text")
        ncatt_put( nc, 0, "geospatial_lon_resolution", as.character(diff(lon(xstats))[[1]]), prec="text")
        ncatt_put( nc, 0, "geospatial_lat_min", as.character(round(min(lat(xstats)))), prec="text")
        ncatt_put( nc, 0, "geospatial_lat_max", as.character(round(max(lat(xstats)))), prec="text")
        ncatt_put( nc, 0, "geospatial_lat_units", "degrees N", prec="text")
        ncatt_put( nc, 0, "geospatial_lat_resolution", as.character(diff(lat(xstats))[[1]]), prec="text")
        ncatt_put( nc, 0, "time_coverage_start", as.character(min(tx)), prec="text")
        ncatt_put( nc, 0, "time_coverage_end", as.character(max(tx)), prec="text")
        if(!is.null(ensemblename)) ncatt_put( nc, 0, "ensemble_name", ensemblename, prec="text" )
        # Add the attributes of object x
        attnames <- names(attributes(x))
        if (verbose) print(attnames)
        attnames <- attnames[!grepl('history',attnames)]
        attnames <- attnames[!grepl('units',attnames)]
        attnames <- attnames[!grepl('variable',attnames)]
        attnames <- attnames[!grepl('dim',attnames)]
        attnames <- attnames[!grepl('index',attnames)]
        attnames <- attnames[!grepl('longitude',attnames)]
        attnames <- attnames[!grepl('latitude',attnames)]
        attnames <- attnames[!grepl('greenwich',attnames)]
        attnames <- attnames[!grepl('domain',attnames)]
        attnames <- attnames[!grepl('call',attnames)]
        attnames <- attnames[!grepl('original_data',attnames)]
        attnames <- attnames[!grepl('scorestats',attnames)]
        attnames <- attnames[!grepl('r.xval',attnames)]
        for (ia in seq_along(attnames)) {
          if(!inherits(attr(x,attnames[ia]),c("matrix","array","data.frame","list"))) {
            if (verbose) print(paste(attnames[ia], paste(attr(x,attnames[ia]), collapse=", ")))
            ncatt_put( nc, 0, attnames[ia], as.character(attr(x,attnames[ia])), prec="text")
          }
        }
        if(!is.null(attr(x, "domain"))) ncatt_put( nc, 0, "predictor_domain", 
                                                   paste0(paste(attr(x,"domain")$lon, collapse="to"), "E/",
                                                          paste(attr(x,"domain")$lat, collapse="-"), "N"), prec="text")
        if(!is.null(attr(x, "r.xval"))) ncatt_put( nc, 0, "r.xval", attr(x,"r.xval")[im,], prec="double")
        
        if(inherits(xstats[[fn]],c("array","field"))) {
          ncatt_put( nc, 0, "spatial_representation", "grid", prec="text")
        } else if(inherits(xstats[[fn]],c("vector","station"))) {
          ncatt_put( nc, 0, "spatial_representation", "station", prec="text")
        }
      }
    
      # Write some values to this variable on disk.
      ncvar_put( nc, varname, xstats[[fn]])
      ncatt_put( nc, varname, "add_offset", offset, prec="float" )
      ncatt_put( nc, varname, "scale_factor", scale, prec="float" ) 
      ncatt_put( nc, varname, "_FillValue", missval, prec="float" ) 
      ncatt_put( nc, varname, "missing_value", missval, prec="float" ) 

      history <- toString(attr(x,'history')$call)
      ncatt_put( nc, varname, "history", history, prec="text" )

      if(inherits(xstats[[fn]],c("array","field"))) {
        ncatt_put( nc, varname, "spatial_representation", "grid", prec="text")
      } else if(inherits(xstats[[fn]],c("vector","station"))) {
        ncatt_put( nc, varname, "spatial_representation", "station", prec="text")
      }

      nc_close(nc)
      if (verbose) print(paste('Finished writing ensemble',fn,'to file', file))
    }
    if (verbose) print(paste('Finished writing downscaled ensemble statistics to file', file))
  }
  if (verbose) print(paste('exit ncdf4_dsensemblestatistics'))
  return(file)
}

ncdf4_dsmodel <- function(x,...,file=NULL,path=NULL,force=TRUE,
                          method="metnoESD",region=NULL,
                          im=NULL, is=NULL, it=NULL, type="field",
                          prec='short',scale=NULL,offset=NULL,
                          torg="1970-01-01",missval=-999,verbose=FALSE) {
  if (verbose) {print('write2ncdf4_dsmodel'); print(names(attributes(x)))}
  if(verbose) print("Function for writing downscaled results for a single model simulation (output from DSensemble) to a netCDF file.")
  stopifnot(inherits(x,"dsensemble") & inherits(x,"list"))
  if(is.null(im)) im <- 1
  if(inherits(x, "pca")|inherits(x, "eof")) {
    x0 <- x
    if("station" %in% type) {
      x <- as.station(x, im=im, is=is, it=it)
    } else if("field" %in% type) {
      x <- as.field(x, im=im, is=is, it=it)
    } else {
      x <- as.field(x, im=im, is=is, it=it) # Transform to field if the format is not specified
    }
  } else {
    x <- subset(x, im=im, is=is, it=it)
  }
  if(inherits(x0[[2]],"season")) {
    class(x) <- c(class(x), "season")
    attr(x, "season") <- season(x0[[2]])[1]
  } else if(inherits(x0[[2]],"annual")) {
    class(x) <- c(class(x), "annual")
  }
  
  if(is.null(file)) {
    #forcing_scenario_bias-baseline_impactmodel_region_1km_ variablename_datetime_years.nc4
    #Example: ecearth-r12i1p1-cclm_rcp26_metnoESD_nordic_0.08x0.08deg_tas_seasavg_djf_1950-2100.nc4
    file <- NULL
    if(length(attr(x, "model_id"))==1) {
      file <- attr(x, "model_id")
    } else file <- attr(x, "model_id")[im]
    file <- paste0(file,"_",attr(x,"scenario"))
    if(!is.null(method)) file <- paste0(file,"_",method)
    if(!is.null(region)) file <- paste0(file,"_",region)
    file <- paste0(file,"_",round(diff(lon(x[[1]]))[1], digits=2),"x",
                   round(diff(lat(x[[1]]))[1], digits=2),"deg")
    file <- paste0(file,"_",attr(x,"variable")[1])
    if(inherits(x,"season")) {
      file <- paste0(file,"_sem")
      if(!is.null(attr(x,"season"))) {
        file <- paste0(file,"_",toupper(attr(x, "season")))
      } else {
        file <- paste0(file,"_",season(x0[[2]]))
      }
    } else if(inherits(x[[1]],"annual")) {
      file <- paste0(file,"_annual-mean")
    } else if(inherits(x[[1]],"month")) {
      file <- paste0(file,"_monthly-mean")
    }
    file <- paste0(file,"_",paste(range(year(x[[1]])),collapse="-"))
    if("station" %in% type) file <- paste0(file,"_station")
    file <- paste0(file, ".nc4")
  }
  if(!is.null(path)) file <- file.path(path, basename(file))
  if(file.exists(file) & !force) {
    nc <- NULL
    if (verbose) print(paste('File',file,'already exists. No new file created.'))
  } else {
    x0 <- x[[1]]
    y <- coredata(x0)
    if (is.null(offset)) offset <- mean(y,na.rm=TRUE)
    if (is.null(scale)) scale <- (max(abs(c(y)),na.rm=TRUE) - offset)/10000
    y <- t(y)
    y <- round((y-offset)/scale)
    y[!is.finite(y)] <- missval
    if (verbose) {
      print(attr(y,'dimensions')); print(c(scale,offset))
      print(range(c(y))); print(range(c(x0),na.rm=TRUE))
    }
    dim(y) <- attr(x0,'dimensions')
    
    dimlon <- ncdim_def( "longitude", "degree_east", lon(x0) )
    dimlat <- ncdim_def( "latitude", "degree_north", lat(x0) )
    if (inherits(index(x),c('numeric','integer'))) {
      mn <- switch(season(x0)[[1]], "djf"="01", "mam"="03", 
                   "jja"="06", "son"="10", "01")
      tx <- as.Date(paste(index(x0),'-',mn,'-01',sep=''))
    }
    
    dimtim <- ncdim_def( "time", paste("days since",torg),
                         as.numeric(as.Date(tx,origin=torg)) )
    x4nc <- ncvar_def(varid(x0)[1], attr(x0,"unit")[1], list(dimlon,dimlat,dimtim), -1, 
                      longname=attr(x0,'longname'), prec=prec)
    
    # Create a netCDF file with this variable
    nc <- nc_create( file, x4nc )
    
    # Write some values to this variable on disk.
    ncvar_put( nc, x4nc, round(y) )
    ncatt_put( nc, x4nc, "add_offset", offset, prec="float" )
    ncatt_put( nc, x4nc, "scale_factor", scale, prec="float" )
    ncatt_put( nc, x4nc, "_FillValue", missval, prec="float" )
    ncatt_put( nc, x4nc, "missing_value", missval, prec="float" )
    history <- toString(attr(x,'history')$call)
    ncatt_put( nc, x4nc, "history", history, prec="text" ) 
    #ncatt_put( nc, 0, 'class', class(x))
    ncatt_put( nc, 0, "description", 
               paste("Saved from esd on",date(),"using write2ncdf4.dsensemble to",
                     "write DSensemble data for a single ensemble member."))
    #if (verbose) print(attr(x,'history'))
    ncatt_put( nc, 0, "esd-version", attr(x,'history')$session$esd.version)
    
    # Add standard attributes
    ncatt_put( nc, 0, "date_created", as.character(as.Date(Sys.time())), prec="text")
    ncatt_put( nc, 0, "units", attr(x0,"unit")[[1]], prec="text")
    ncatt_put( nc, 0, "geospatial_lon_min", as.character(round(min(lon(x0)))), prec="text")
    ncatt_put( nc, 0, "geospatial_lon_max", as.character(round(max(lon(x0)))), prec="text")
    ncatt_put( nc, 0, "geospatial_lon_units", "degrees E", prec="text")
    ncatt_put( nc, 0, "geospatial_lon_resolution", as.character(diff(lon(x0))[[1]]), prec="text")
    ncatt_put( nc, 0, "geospatial_lat_min", as.character(round(min(lat(x0)))), prec="text")
    ncatt_put( nc, 0, "geospatial_lat_max", as.character(round(max(lat(x0)))), prec="text")
    ncatt_put( nc, 0, "geospatial_lat_units", "degrees N", prec="text")
    ncatt_put( nc, 0, "geospatial_lat_resolution", as.character(diff(lat(x0))[[1]]), prec="text")
    ncatt_put( nc, 0, "time_coverage_start", as.character(min(index(x0))), prec="text")
    ncatt_put( nc, 0, "time_coverage_end", as.character(max(index(x0))), prec="text")
    #ncatt_put( nc, 0, "lon_bnds", as.character(round(min(lon(x0)))), prec="text")
    #ncatt_put( nc, 0, "lat_bnds", as.character(round(min(lon(x0)))), prec="text")
    
    # Add the attributes of object x
    attnames <- names(attributes(x))
    if (verbose) print(attnames)
    attnames <- attnames[!grepl('history',attnames)]
    attnames <- attnames[!grepl('units',attnames)]
    attnames <- attnames[!grepl('variable',attnames)]
    attnames <- attnames[!grepl('dim',attnames)]
    attnames <- attnames[!grepl('index',attnames)]
    attnames <- attnames[!grepl('longitude',attnames)]
    attnames <- attnames[!grepl('latitude',attnames)]
    attnames <- attnames[!grepl('greenwich',attnames)]
    attnames <- attnames[!grepl('domain',attnames)]
    attnames <- attnames[!grepl('call',attnames)]
    attnames <- attnames[!grepl('original_data',attnames)]
    attnames <- attnames[!grepl('scorestats',attnames)]
    attnames <- attnames[!grepl('r.xval',attnames)]
    for (ia in seq_along(attnames)) {
      if(!inherits(attr(x,attnames[ia]),c("matrix","array","data.frame","list"))) {
        if (verbose) print(paste(attnames[ia], paste(attr(x,attnames[ia]), collapse=", ")))
        ncatt_put( nc, 0, attnames[ia], as.character(attr(x,attnames[ia])), prec="text")
      }
    }
    if(!is.null(attr(x, "domain"))) ncatt_put( nc, 0, "predictor_domain", 
                                               paste0(paste(attr(x,"domain")$lon, collapse="to"), "E/",
                                                      paste(attr(x,"domain")$lat, collapse="-"), "N"), prec="text")
    if(!is.null(attr(x, "r.xval"))) ncatt_put( nc, 0, "r.xval", attr(x,"r.xval")[im,], prec="double")
    
    if(inherits(x,"field")) {
      ncatt_put( nc, 0, "spatial_representation", "grid", prec="text")
    } else if(inherits(x,"station")) {
      ncatt_put( nc, 0, "spatial_representation", "station", prec="text")
    }

    if (verbose) print(paste('Finished writing dsensemble for a single model to file', file))
    nc_close(nc)
  }
  if (verbose) print(paste('exit ncdf4_dsmodel'))
  return(file)
}

