## ================================================
## Plan for the function write2ncdf4.dsensemble
##
## The input argument 'type' is used to redirect to different functions:
##
## type = c('dsensemble') -> ncdf4_dsensemble -> write PCA based dsensemble to a netCDF file, including an extraction method as an attribute
## type = c('dsensemblestats', 'station') -> ncdf4_dsensemblestats -> convert PCA based dsensemble to station objects, calculate ensemble statistics, write to a netCDF file
## type = c('dsensemblestats', 'field') -> ncdf4_dsensemblestats -> convert PCA based dsensemble to field objects, calculate ensemble statistics, write to a netCDF file
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
## ncdf4_dsensemblestats is empty
##
## ncdf_model works for writing data as gridded fields, but not for station data.
##
## ================================================

#' Saves climate data as netCDF.
#' 
#' Method to save station data as netCDF, making sure to include the data
#' structure and meta-data (attributes). The code tries to follow the netCDf
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
#' @param prec Precision: see \code{\link[ncdf4]{ncvar_def}}
#' @param scale Sets the attribute 'scale_factor' which is used to scale
#' (multiply) the values stored (to save space may be represented as 'short').
#' @param im a model index, see \code{\link{subset}}
#' @param im a numerical or numerical vector with indices of the ensemble members to be included. If NULL include all.
#' @param type a vector defining the type of data to write to the netCDF file. 
#' The default type is "dsensemble" which writes dsensemble data in its original form.
#' Other options are "ensemblestats" to calculate and write ensemble statistics,
#' or "model" which writes downscaled results for a single model simulation   
#' (the ensemble member is then selected with the input argument \code{'im'}).
#' The options 'ensemblestats' and 'model' can be used in conjunction with the additional argument, 
#' e.g. c('model','field') or c('model','station') to save the downscaled results as a field (gridded data) 
#' or a 'station' object (the stations that were used as predictand in the downscaling).
#' @param offset Sets the attribute 'add_offset' which is added to the values
#' stored (to save space may be represented as 'short').
#' @param torg Time origin
#' @param missval Missing value: see \code{\link[ncdf4]{ncvar_def}}
#' @param id	An identifier for the data set, provided by and unique within its naming authority. The combination of the "naming authority" and the "id" should be globally unique, but the id can be globally unique by itself also. IDs can be URLs, URNs, DOIs, meaningful text strings, a local key, or any other unique string of characters. The id should not include white space characters.
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
#' @param verbose TRUE - clutter the screen.
#' @param \dots additional arguments
#' 
#' @return None
#' 
#' @keywords netcdf ncdf4 save
#' 
#' @exportS3Method
#' @export write2ncdf4.dsensemble
write2ncdf4.dsensemble <- function(x,...,file=NULL,path=NULL,
                                   prec='short',scale=NULL,offset=NULL,
                                   type='dsensemble',im=1,
                                   torg="1970-01-01",missval=-999,
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
    for(i in im) {
      write2ncdf4(x,...,im=i,file=file,path=path,
                  prec=prec,scale=scale,offset=offset,
                  type=type,
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
                  additional_attributes=additional_attributes,verbose=verbose)
    }
  } else {
    if("dsensemble" %in% type) {
      nc <- ncdf4_dsensemble(x,...,file=file,path=path,prec=prec,scale=scale,type=type,
                             offset=offset,torg=torg,missval=missval,verbose=verbose)
    } else if ("model" %in% type) {
      nc <- ncdf4_dsmodel(x,...,file=file,path=path,prec=prec,scale=scale,type=type,
                          offset=offset,torg=torg,missval=missval,im=im,verbose=verbose)
    } else if ("ensemblestats" %in% type) {
      nc <- ncdf4_dsensemblestats(x,...,file=file,path=path,prec=prec,scale=scale,type=type,
                                  offset=offset,torg=torg,missval=missval,im=im,verbose=verbose)
    }
    
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
      if (verbose) print(paste('Successfully finished writing dsensemble data and attributes to file', file))
    }
  }
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
  return(nc)
}

ncdf4_dsensemble <- function(x,...,file='esd.dsensemble.nc',prec='short',offset=0,scale=0.1,
                             torg="1970-01-01",missval=-99,verbose=TRUE) {
  ## prec - see http://james.hiebert.name/blog/work/2015/04/18/NetCDF-Scale-Factors/
  if (verbose) print('ncdf4_dsensemble')
  class.x <- class(x)
  ngcms <- length(x) - 2
  ## Get the two first elements of the list which contain information common to the rest
  info <- x$info
  if (!is.null(x$pca)) pca <- x$pca
  if (!is.null(x$eof)) pca <- x$eof
  lons <- lon(pca)
  lats <- lat(pca)
  
  ## Clear these so that the rest are zoo objects describing the model runs
  x$info <- NULL
  x$pca <- NULL
  x$eof <- NULL
  names.x <- names(x)
  if (verbose) print(names.x)
  
  ## Define the variables and dimensions
  if (verbose) {print('Check index type - set to year'); print(index(x[[1]]))}
  if (is.list(x)) {
    if (class(index(x))=='Date') tim <- year(x[[1]]) + (month(x[[1]])-1)/12 else
      tim <- year(x[[1]])
    dimgcm <- dim(x[[1]])
    nloc <- length(x)
  } else {
    if (class(index(x))=='Date') tim <- year(x) + (month(x)-1)/12 else
      tim <- year(x)
    dimgcm <- dim(x)
    nloc <- 1
  }
  if (!is.null(pca)) {
    pcaatts <- names(attributes(pca))
    pattern <- attr(pca,'pattern')
    dpat <- dim(pattern)
  } else dpat <- 1
  if (verbose) {print(pcaatts); print(dim(pattern))}
  
  ## Set dimensions  
  dimtim <- ncdim_def( "time", 'year', as.integer(tim) )
  dimens <- ncdim_def( "ensemble_member", 'index', 1:nloc )
  if (!is.null(pca)) {
    if (class(index(pca))=='Date') index(pca) <- year(pca) + (month(pca)-1)/12 else
      index(pca) <- year(pca)
    dimpca <- ncdim_def( "i_pca", "index", 1:dimgcm[2] )
    dimtimpca <- ncdim_def( "time_pca", "year", index(pca) )
    if (length(dpat)==2) {
      dimxy <- ncdim_def( "space_pca", "index", 1:dim(pattern)[1] )
    } else {
      dimx <- ncdim_def( "longitude", "degrees_east", lons )
      dimy <- ncdim_def( "latitude", "degrees_north", lats )
    }
  }
  dimsea <- ncdim_def( "season", "index", 1:4 )
  if (verbose) print(tim)
  varlist <- list() # For single-station dsensemble objects
  varlist$gcm <- ncvar_def("gcm", "weights", list(dimtim,dimpca,dimens), missval, 
                           longname='principal components', prec=prec)
  
  ## set up the dimensions of the PCA
  if (!is.null(pca)) {
    varlist$pca <- ncvar_def("PC", "weights", list(dimtimpca,dimpca), missval, 
                             longname='principal components', prec=prec)
    ## EOFs and PCA have different number of dimensions
    if (length(dpat)==2) {
      if (verbose) print('--- PCA ---')
      varlist$pat <- ncvar_def("pattern", "weights", list(dimxy,dimpca), missval, 
                               longname='principal component analysis patterns', prec=prec)
      varlist$lon <- ncvar_def("longitude", "degree_east", dimxy,missval, 
                               longname='longitude', prec='float')
      varlist$lat <- ncvar_def("latitude", "degree_north", dimxy,missval, 
                               longname='latitude', prec='float')
      varlist$alt <- ncvar_def("altitude", "m", dimxy,missval, 
                               longname='altitude', prec='float')
      varlist$stid <- ncvar_def("station_id", "number", dimxy,"NA", 
                                longname='station ID', prec="char")
      varlist$loc <- ncvar_def("location", "name", dimxy,"NA", 
                               longname='location name', prec="char")
      varlist$cntr <- ncvar_def("country", "name", dimxy,"NA", 
                                longname='country name', prec="char")
      varlist$src <- ncvar_def("src", "name", dimxy,"NA", 
                               longname='source', prec="char")
    } else
      if (length(dpat)==3) {
        if (verbose) print('--- EOF ---')
        varlist$pat <- ncvar_def("pattern", "weights", list(dimx,dimy,dimpca), missval, 
                                 longname='principal component analysis patterns', prec=prec)
      }
    varlist$lambda <- ncvar_def("lambda", "number", dimpca,missval, 
                                longname='eigenvalues', prec="float")
  }
  
  if (verbose) print(names(varlist))
  
  ## Create a netCDF file with this variable
  if (verbose) print('Create netCDF-file')
  
  nc <- nc_create( file, varlist,verbose=verbose)
  
  if (verbose) print('write pca/eof data')                   
  ## Add the information stored in the list elements as 2D zoo objects
  X <- unlist(lapply(x,function(x) x[1:239,]))
  X <- round((X - offset)/scale)
  if (verbose) print(c(length(X),dim(x[[1]]),length(x)))
  if (verbose) print(table(unlist(lapply(x,dim))))
  dim(X) <- c(dim(x[[1]]),length(x))
  if (verbose) print('write zoo data')     
  ## Write some values to this variable on disk: GCM results
  ncvar_put( nc, varlist$gcm, X )
  ncatt_put( nc, varlist$gcm, "add_offset", offset, prec="float" )
  ncatt_put( nc, varlist$gcm, "scale_factor", scale, prec="float" ) 
  ncatt_put( nc, varlist$gcm, "_FillValue", missval, prec="float" ) 
  ncatt_put( nc, varlist$gcm, "missing_value", missval, prec="float" )
  ## GCM names
  if (verbose) print('write GCM names')  
  ncatt_put( nc, varlist$gcm, "GCM runs", paste(names.x,collapse=','),prec="char" )
  ## PCA/EOF variables:
  if (verbose) print('EOF/PCA variables')
  ncvar_put( nc, varlist$pca, round((pca - offset)/scale) )
  ncatt_put( nc, varlist$pca, "add_offset", offset, prec="float" )
  ncatt_put( nc, varlist$pca, "scale_factor", scale, prec="float" ) 
  ncatt_put( nc, varlist$pca, "_FillValue", missval, prec="float" ) 
  ncatt_put( nc, varlist$pca, "missing_value", missval, prec="float" ) 
  ncvar_put( nc, varlist$pat, round((pattern - offset)/scale) )
  ncatt_put( nc, varlist$pat, "add_offset", offset, prec="float" )
  ncatt_put( nc, varlist$pat, "scale_factor", scale, prec="float" ) 
  ncatt_put( nc, varlist$pat, "_FillValue", missval, prec="float" ) 
  ncatt_put( nc, varlist$pat, "missing_value", missval, prec="float" ) 
  ncatt_put( nc, varlist$pat, "dimensions_pca", paste(attr(pattern,'dimensions'),collapse=', '), prec="char" )
  ## If the object contains PCAs
  if (length(dpat)==2) {
    if (verbose) print('PCA only variables')
    ncvar_put( nc, varlist$lon, lon(x) )
    ncvar_put( nc, varlist$lat, lat(x) )
    ncvar_put( nc, varlist$alt, alt(x) )
    ncvar_put( nc, varlist$loc, loc(x) )
    ncvar_put( nc, varlist$stid, as.character(stid(x)) )
    ncvar_put( nc, varlist$cntr, cntr(x) )
    ncvar_put( nc, varlist$src, src(x) )
    ncvar_put( nc, varlist$lambda, attr(pca,'eigenvalues') )
  } else if (length(dpat)==3)
    ncvar_put( nc, varlist$lambda, attr(pca,'eigenvalues') )
  if (verbose) print('history')
  ncatt_put( nc, varlist$pca, "history", paste(attr(x,'history'),collapse=';'), prec="char" )  
  ## Global attributes:
  ncatt_put( nc, 0, 'class', class(x))
  ncatt_put( nc, 0, "description", "Saved from esd using write2ncdf4.dsensemble")
  #ncatt_put( nc, 0, "class", paste(class.x,collapse='-'))
  ncatt_put( nc, 0, "esd-version", attr(x,'history')$session$esd.version)
  if (verbose) print(paste('Finished writing dsensemble to file', file))
  return(nc)
}

ncdf4_dsensemblestats <- function(x,...,file='esd.dsensemble.nc',prec='short',offset=0,scale=0.1,
                                  type="field",torg="1970-01-01",missval=-99,verbose=TRUE) {
  ## prec - see http://james.hiebert.name/blog/work/2015/04/18/NetCDF-Scale-Factors/
  if (verbose) print('ncdf4_dsensemblestats')
  print("Function ncdf4_dsensemblestats is not finished - returning nothing!")
  nc <- NULL
  return(nc)
}

ncdf4_dsmodel <- function(x,...,file=NULL,path=NULL,
                          method="metnoESD",region="nordic",
                          im=1, is=NULL, it=NULL, type="field",
                          prec='short',scale=NULL,offset=NULL,
                          torg="1970-01-01",missval=-999,verbose=FALSE) {
  if (verbose) {print('write2ncdf4_dsmodel'); print(names(attributes(x)))}
  if(verbose) print("Function for writing downscaled results for a single model simulation (output from DSensemble) to a netCDF file.")
  stopifnot(inherits(x,"dsensemble") & inherits(x,"list"))
  
  if(inherits(x, "pca")|inherits(x, "eof")) {
    if("station" %in% type) {
      x <- as.station(x, im=im)
    } else if("field" %in% type) {
      x <- as.field(x, im=im)
    } else {
      x <- as.field(x, im=im) # Transform to field if the format is not specified
    }
  }
  
  if(is.null(file)) {
    #forcing_scenario_bias-baseline_impactmodel_region_1km_ variablename_datetime_years.nc4
    #Example: ecearth-r12i1p1-cclm_rcp26_eqm-sn2018v2005_rawbc_norway_1km_tas_daily_2021.nc4
    #ESD example: ecearth-r12i1p1-cclm_rcp26_metnoESD_nordic_0.08x0.08deg_tas_seasavg_djf_1950-2100.nc4
    file <- NULL
    if(length(attr(x, "model_id"))==1) {
      file <- attr(x, "model_id")
    } else file <- attr(x, "model_id")[im]
    file <- paste0(file,"_",attr(x,"scenario"))
    file <- paste0(file,"_",method)
    file <- paste0(file,"_",region)
    file <- paste0(file,"_",round(diff(lon(x[[1]]))[1], digits=2),"x",
                   round(diff(lat(x[[1]]))[1], digits=2),"deg")
    file <- paste0(file,"_",attr(x,"variable")[1])
    if(inherits(x[[1]],"season")) {
      file <- paste0(file,"_sem")
      file <- paste0(file,"_",toupper(season(x[[1]])[1]))
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
             paste("Saved from esd using write2ncdf4",date()))
  #if (verbose) print(attr(x,'history'))
  ncatt_put( nc, 0, "esd-version", attr(x,'history')$session$esd.version)
   
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
    if(inherits(attr(x,attnames[ia]),"matrix")) {
      browser()
    } else {
      if (verbose) print(paste(attnames[ia], paste(attr(x,attnames[ia]), collapse=", ")))
      ncatt_put( nc, 0, attnames[ia], as.character(attr(x,attnames[ia])), prec="text")
    }
  }
  if(!is.null(attr(x, "domain"))) ncatt_put( nc, 0, "predictor_domain", 
                                             paste0(paste(attr(x,"domain")$lon, collapse="-"), "E/",
                                                    paste(attr(x,"domain")$lat, collapse="-"), "N"), prec="text")
  if(!is.null(attr(x, "r.xval"))) ncatt_put( nc, 0, "r.xval", attr(x,"r.xval"), prec="double")

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
  
  if(inherits(x,"field")) {
    ncatt_put( nc, 0, "spatial_representation", "grid", prec="text")
  } else if(inherits(x,"station")) {
    ncatt_put( nc, 0, "spatial_representation", "station", prec="text")
  }
  
  if (verbose) print(paste('Finished writing dsensemble for a single model to file', file))
  return(nc)
}



## OLD VERSION OF write2ncdf4.dsensemble - like ncdf4_dsensemble
# Saves climate data as netCDF.
#
# Method to save 'dsensemble' data as netCDF, making sure to include the data
# structure and meta-data (attributes). The code tries to follow the netCDf
# 'CF' convention. The method is built on the \code{ncdf4} package.
#
# To save space, the values are saved as short (16-bit signed integer that
# can hold values between -32768 and 32767).
# (see NC_SHORT in \url{https://www.unidata.ucar.edu/software/netcdf/docs/data_type.html}).
#
# @seealso write2ncdf4
#
# @param x data object
# @param file filename
# @param prec Precision: see \code{\link[ncdf4]{ncvar_def}}
# @param scale Sets the atttribute 'scale_factor' which is used to scale
# (multiply) the values stored (to save space may be represented as 'short').
# @param offset Sets the attribute 'add_offset' which is added to the values
# stored (to save space may be represented as 'short').
# @param torg Time origin
# @param missval Missing value: see \code{\link[ncdf4]{ncvar_def}}
# @param verbose If TRUE print progress
# @param \dots additional arguments
# 
# @return None
# 
# @keywords netcdf ncdf4 save
# 
# @exportS3Method
# @export write2ncdf4.dsensemble 
#write2ncdf4.dsensemble <- function(x,...,file='esd.dsensemble.nc',prec='short',offset=0,scale=0.1,
write2ncdf4_dsensemble_old <- function(x,...,file='esd.dsensemble.nc',prec='short',offset=0,scale=0.1,
                                       torg="1970-01-01",missval=-99,verbose=TRUE) {
  ## prec - see http://james.hiebert.name/blog/work/2015/04/18/NetCDF-Scale-Factors/
  if (verbose) print('write2ncdf4.field')
  class.x <- class(x)
  ngcms <- length(x) - 2
  ## Get the two first elements of the list which contain information common to the rest
  info <- x$info
  if (!is.null(x$pca)) pca <- x$pca
  if (!is.null(x$eof)) pca <- x$eof
  lons <- lon(pca)
  lats <- lat(pca)
  
  ## Clear these so that the rest are zoo objects describing the model runs
  x$info <- NULL
  x$pca <- NULL
  x$eof <- NULL
  names.x <- names(x)
  if (verbose) print(names.x)
  
  ## Define the variables and dimensions
  if (verbose) {print('Check index type - set to year'); print(index(x[[1]]))}
  if (is.list(x)) {
    if (class(index(x))=='Date') tim <- year(x[[1]]) + (month(x[[1]])-1)/12 else
      tim <- year(x[[1]])
    dimgcm <- dim(x[[1]])
    nloc <- length(x)
  } else {
    if (class(index(x))=='Date') tim <- year(x) + (month(x)-1)/12 else
      tim <- year(x)
    dimgcm <- dim(x)
    nloc <- 1
  }
  if (!is.null(pca)) {
    pcaatts <- names(attributes(pca))
    pattern <- attr(pca,'pattern')
    dpat <- dim(pattern)
  } else dpat <- 1
  if (verbose) {print(pcaatts); print(dim(pattern))}
  
  ## Set dimensions  
  dimtim <- ncdim_def( "time", 'year', as.integer(tim) )
  dimens <- ncdim_def( "ensemble_member", 'index', 1:nloc )
  if (!is.null(pca)) {
    if (class(index(pca))=='Date') index(pca) <- year(pca) + (month(pca)-1)/12 else
      index(pca) <- year(pca)
    dimpca <- ncdim_def( "i_pca", "index", 1:dimgcm[2] )
    dimtimpca <- ncdim_def( "time_pca", "year", index(pca) )
    if (length(dpat)==2) {
      dimxy <- ncdim_def( "space_pca", "index", 1:dim(pattern)[1] )
    } else {
      dimx <- ncdim_def( "longitude", "degrees_east", lons )
      dimy <- ncdim_def( "latitude", "degrees_north", lats )
    }
  }
  dimsea <- ncdim_def( "season", "index", 1:4 )
  if (verbose) print(tim)
  varlist <- list() # For single-station dsensemble objects
  varlist$gcm <- ncvar_def("gcm", "weights", list(dimtim,dimpca,dimens), missval, 
                           longname='principal components', prec=prec)
  
  ## set up the dimensions of the PCA
  if (!is.null(pca)) {
    varlist$pca <- ncvar_def("PC", "weights", list(dimtimpca,dimpca), missval, 
                             longname='principal components', prec=prec)
    ## EOFs and PCA have different number of dimensions
    if (length(dpat)==2) {
      if (verbose) print('--- PCA ---')
      varlist$pat <- ncvar_def("pattern", "weights", list(dimxy,dimpca), missval, 
                               longname='principal component analysis patterns', prec=prec)
      varlist$lon <- ncvar_def("longitude", "degree_east", dimxy,missval, 
                               longname='longitude', prec='float')
      varlist$lat <- ncvar_def("latitude", "degree_north", dimxy,missval, 
                               longname='latitude', prec='float')
      varlist$alt <- ncvar_def("altitude", "m", dimxy,missval, 
                               longname='altitude', prec='float')
      varlist$stid <- ncvar_def("station_id", "number", dimxy,"NA", 
                                longname='station ID', prec="char")
      varlist$loc <- ncvar_def("location", "name", dimxy,"NA", 
                               longname='location name', prec="char")
      varlist$cntr <- ncvar_def("country", "name", dimxy,"NA", 
                                longname='country name', prec="char")
      varlist$src <- ncvar_def("src", "name", dimxy,"NA", 
                               longname='source', prec="char")
    } else
      if (length(dpat)==3) {
        if (verbose) print('--- EOF ---')
        varlist$pat <- ncvar_def("pattern", "weights", list(dimx,dimy,dimpca), missval, 
                                 longname='principal component analysis patterns', prec=prec)
      }
    varlist$lambda <- ncvar_def("lambda", "number", dimpca,missval, 
                                longname='eigenvalues', prec="float")
  }
  
  if (verbose) print(names(varlist))
  
  ## Create a netCDF file with this variable
  if (verbose) print('Create netCDF-file')
  
  ncnew <- nc_create( file, varlist,verbose=verbose)
  
  if (verbose) print('write pca/eof data')                   
  ## Add the information stored in the list elements as 2D zoo objects
  X <- unlist(lapply(x,function(x) x[1:239,]))
  X <- round((X - offset)/scale)
  if (verbose) print(c(length(X),dim(x[[1]]),length(x)))
  if (verbose) print(table(unlist(lapply(x,dim))))
  dim(X) <- c(dim(x[[1]]),length(x))
  if (verbose) print('write zoo data')     
  ## Write some values to this variable on disk: GCM results
  ncvar_put( ncnew, varlist$gcm, X )
  ncatt_put( ncnew, varlist$gcm, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, varlist$gcm, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, varlist$gcm, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, varlist$gcm, "missing_value", missval, prec="float" )
  ## GCM names
  if (verbose) print('write GCM names')  
  ncatt_put( ncnew, varlist$gcm, "GCM runs", paste(names.x,collapse=','),prec="char" )
  ## PCA/EOF variables:
  if (verbose) print('EOF/PCA variables')
  ncvar_put( ncnew, varlist$pca, round((pca - offset)/scale) )
  ncatt_put( ncnew, varlist$pca, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, varlist$pca, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, varlist$pca, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, varlist$pca, "missing_value", missval, prec="float" ) 
  ncvar_put( ncnew, varlist$pat, round((pattern - offset)/scale) )
  ncatt_put( ncnew, varlist$pat, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, varlist$pat, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, varlist$pat, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, varlist$pat, "missing_value", missval, prec="float" ) 
  ncatt_put( ncnew, varlist$pat, "dimensions_pca", paste(attr(pattern,'dimensions'),collapse=', '), prec="char" )
  ## If the object contains PCAs
  if (length(dpat)==2) {
    if (verbose) print('PCA only variables')
    ncvar_put( ncnew, varlist$lon, lon(x) )
    ncvar_put( ncnew, varlist$lat, lat(x) )
    ncvar_put( ncnew, varlist$alt, alt(x) )
    ncvar_put( ncnew, varlist$loc, loc(x) )
    ncvar_put( ncnew, varlist$stid, as.character(stid(x)) )
    ncvar_put( ncnew, varlist$cntr, cntr(x) )
    ncvar_put( ncnew, varlist$src, src(x) )
    ncvar_put( ncnew, varlist$lambda, attr(pca,'eigenvalues') )
  } else if (length(dpat)==3)
    ncvar_put( ncnew, varlist$lambda, attr(pca,'eigenvalues') )
  if (verbose) print('history')
  ncatt_put( ncnew, varlist$pca, "history", paste(attr(x,'history'),collapse=';'), prec="char" )  
  ## Global attributes:
  ncatt_put( ncnew, 0, 'class', class(x))
  ncatt_put( ncnew, 0, "description", "Saved from esd using write2ncdf4.dsensemble")
  #ncatt_put( ncnew, 0, "class", paste(class.x,collapse='-'))
  ncatt_put( ncnew, 0, "esd-version", attr(x,'history')$session$esd.version)
  nc_close(ncnew)
  if (verbose) print(paste('Finished successfully - file', file))  
}
