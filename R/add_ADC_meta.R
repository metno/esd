#https://adc.met.no/submit-data-as-netcdf-cf
add_ADC_meta <- function(ncid,x,conventions=NA,title=NA,summary=NA,project=NA,license=NA,
                         featuretype="timeSeries",keywords=NA,keywordvoc=NA,
                         signature_file='~/.esd_add_ADC_meta2nctCDF',verbose=FALSE) {
  if (verbose) cat('add_ADC_meta \n')
  if (!file.exists(signature_file)) {
    cat("You haven't specified your details for netCDF metadata neededfor write2ncdf4 \n")
    cat("See https://adc.met.no/submit-data-as-netcdf-cf for guidance \n")
    name <- readline("Name: ")
    email <- readline("e-mail: ")
    institution <- readline("institution: ")
    type <- readline("creator_type ('person', 'group', 'institution', or 'position'): ")
    signature <- data.frame(V1=c("creator_name","creator_email","creator_institution","creator_type"),
                            V2=c(name,email,institution,type))
    write.table(signature,file=signature_file,sep=':')
  } 
  sign <- read.table(signature_file,sep=':')
  if (!is.na(conventions)) ncatt_put( ncid, 0, "Conventions",conventions)
  ncatt_put( ncid, 0, "creator_name", sign$V2[grep('name',sign$V1)])
  ncatt_put( ncid, 0, "creator_email", sign$V2[grep('mail',sign$V1)])
  ncatt_put( ncid, 0, "creator_institution",sign$V2[grep('institution',sign$V1)])
  ncatt_put( ncid, 0, "creator_type",sign$V2[grep('type',sign$V1)])
  ncatt_put( ncid, 0, "date_created",as.Date(Sys.time()))
  ncatt_put( ncid, 0, "featureType",featuretype)
  ncatt_put( ncid, 0, "geospatial_lat_max",max(lat(x)))
  ncatt_put( ncid, 0, "geospatial_lat_min",min(lat(x)))
  ncatt_put( ncid, 0, "geospatial_lon_max",max(lon(x)))
  ncatt_put( ncid, 0, "geospatial_lon_min",min(lon(x)))
  ncatt_put( ncid, 0, "history",paste(attr(x,'history')$call,
                                      collapse='; '))
  if (!is.na(keywords)) ncatt_put( ncid, 0, "keywords",keywords)
  if (!is.na(keywordvoc)) ncatt_put( ncid, 0, "keywords_vocabulary",keywordvoc)
  if (!is.na(license)) ncatt_put( ncid, 0, "license",license)
  if (!is.na(project)) ncatt_put( ncid, 0, "project",project)
  if (!is.na(title)) ncatt_put( ncid, 0, "summary",summary)
  if (inherits(x,'zoo')) { 
    ncatt_put( ncid, 0, "time_coverage_end",max(index(x)))
    ncatt_put( ncid, 0, "time_coverage_start",min(index(x)))
  } else { 
    ncatt_put( ncid, 0, "time_coverage_end",'NA')
    ncatt_put( ncid, 0, "time_coverage_start",'NA')
  }
  if (is.na(title)) title <- paste(attr(x,'longname'),'from',attr(x,'source'))
  if (is.na(title))ncatt_put( ncid, 0, "title",title)
}
