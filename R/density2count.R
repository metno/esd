density2count <- function(y,it=NULL,is=NULL,verbose=TRUE) {
  if(verbose) print("density2number - estimate the storm count based on cyclone density")
  if(inherits(y,"trajectory")) y <- as.events(y)
  if(inherits(y,"events")) y <- as.field(y)
  stopifnot(inherits(y,"field"))
  
  if(!is.null(attr(y,"unitarea"))) { unitarea <- attr(y,"unitarea")
  } else unitarea <- pi*7E5**2
  if(verbose) print(paste("Unit area of cyclone density:",round(unitarea*1E-6),"km2"))
  
  y <- subset(y,is=is,it=it)
  y.mean <- aggregate.area(y,FUN="mean")
  A <- area.lonlat(c(min(longitude(y)),min(longitude(y)),max(longitude(y)),max(longitude(y))),
                   c(max(latitude(y)),min(latitude(y)),min(latitude(y)),max(latitude(y))))
  if(verbose) print(paste("Area of predictand domain:",round(A*1E-6),"km2"))
  N <- y.mean*A/unitarea
  attr(N,"variable") <- "storm~count"
  attr(N,"longname") <- "monthly mean number of cyclones"
  attr(N,"unit") <- "events/month"
  attr(N,"aspect") <- "density2count"
  invisible(N)
}
