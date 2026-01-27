
add_ADC_meta <- function(ncid,x,conventions="CF-hybrid",email='NA',institution='met.no',
                         creatorname='NA',creatortype='NA',featuretype='NA',keywords=NULL,
                         keywordvoc='NA',license='NA',project='NA',summary='NA',title='NA',verbose=FALSE) {
  if (verbose) cat('add_ADC_meta \n')
  ncatt_put( ncid, 0, "Conventions",conventions)
  ncatt_put( ncid, 0, "creator_email", email)
  ncatt_put( ncid, 0, "creator_institution",institution)
  ncatt_put( ncid, 0, "creator_name", creatorname)
  ncatt_put( ncid, 0, "creator_type",creatortype)
  ncatt_put( ncid, 0, "date_created",as.Date(Sys.time()))
  ncatt_put( ncid, 0, "featureType",featuretype)
  ncatt_put( ncid, 0, "geospatial_lat_max",max(lat(x)))
  ncatt_put( ncid, 0, "geospatial_lat_min",min(lat(x)))
  ncatt_put( ncid, 0, "geospatial_lon_max",max(lon(x)))
  ncatt_put( ncid, 0, "geospatial_lon_min",min(lon(x)))
  ncatt_put( ncid, 0, "history",paste(as.character(attr(x,'history')),collapse=' - '))
  ncatt_put( ncid, 0, "keywords",keywords)
  ncatt_put( ncid, 0, "keywords_vocabulary",keywordvoc)
  ncatt_put( ncid, 0, "license",license)
  ncatt_put( ncid, 0, "project",project)
  ncatt_put( ncid, 0, "summary",summary)
  if (inherits(x,'zoo')) { 
    ncatt_put( ncid, 0, "time_coverage_end",max(index(x)))
    ncatt_put( ncid, 0, "time_coverage_start",min(index(x)))
  } else { 
    ncatt_put( ncid, 0, "time_coverage_end",'NA')
    ncatt_put( ncid, 0, "time_coverage_start",'NA')
  }
  ncatt_put( ncid, 0, "title",title)
}
