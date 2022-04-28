#' 
#' A function that lists files on thredds catalogue.
#' @aliases retrieve
#'  
#' @param caturl url of the thredds ctalogue containing the files
#' @param extension the extension ot the netCDF files (default is ".nc")
#' @param pattern text pattern to filter file names
#' @param verbose if true promt test lines
#'
#' @examples
#' 
#' list_thredds()
#' list_thredds(pattern=c("rcp45","r1i1p1","RACMO"))
#' list_thredds(caturl="http://esgf3.dkrz.de/thredds/catalog/esgcet/12/CMIP6.ScenarioMIP.DKRZ.MPI-ESM1-2-HR.ssp585.r1i1p1f1.day.pr.gn.v20190710.html")
#' caturl <- "https://thredds.met.no/thredds/catalog/KSS/Klima_i_Norge_2100/bias_corrected/3DBC/cross-validation/noresm-r1i1p1-remo/tasmin/catalog.html"
#' extension <- '.nc4'
#' list_thredds(caturl,extension)
#' @export

list_thredds <- function(caturl="https://thredds.met.no/thredds/catalog/KSS/Klima_i_Norge_2100/seasonal_RCM/catalog.html", 
                         extension=".nc", pattern = NULL, verbose=FALSE) {
  
  if (verbose) print(paste('list_thredds',caturl))
  
  URL <- url(caturl)
  readLines(URL) -> xmlcode
  close(URL)
  files<- xmlcode[grep('.nc',xmlcode,fixed=TRUE)]
  files <- gsub('.*<tt>','',files)
  files <- gsub('.*<code>','',files)
  files <- gsub('</tt></a></td>','',files)
  files <- gsub('</code></a>','',files)
  if (!is.null(pattern)) for (i in 1:length(pattern)) files <- files[grep(pattern[i],files,ignore.case = TRUE)]
  attr(files,'caturl') <- caturl
  attr(files,'pattern') <- pattern
  attr(files,'history') <- history.stamp()
  return(files)
}

