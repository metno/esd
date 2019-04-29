density.trajectory <- function(x,it=NULL,is=NULL,dx=2,dy=2,radius=5E5,verbose=FALSE) {
  if(verbose) print("density.trajectory")
  y <- subset(x,it=it,is=is)
  if(!is.null(dim(y))) {
    A <- apply(y,1,function(x) trackdensity(x[colnames(y)=='lon'],
                                          x[colnames(y)=='lat'],
                                          dx=dx,dy=dy,radius=radius))
    lon <- unlist(lapply(A,function(x) factor2numeric(x$lon)))
    lat <- unlist(lapply(A,function(x) factor2numeric(x$lat)))
    hits <- as.data.frame(table(lon,lat))
    names(hits)[names(hits)=="Freq"] <- "density"
  } else {
    A <- trackdensity(y[names(y)=='lon'],y[names(y)=='lat'],dx=dx,dy=dy,radius=radius)
    lon <- factor2numeric(A$lon)
    lat <- factor2numeric(A$lat)
    hits <- as.data.frame(table(lon,lat))
    names(hits)[names(hits)=="Freq"] <- "density"
  }
  invisible(hits)
}
