#' 
#' A function that lists files on thredds catalogue.
#' @aliases retrieve
#'  
#' @param caturl
#' @param extension
#' @param pattern
#' @param verbose
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



#     ## Need to extract information on the file names from the XML mess:
#     data_heads <- sub('.*catalog/','catalog.html?dataset=',sub( 'catalog.html','',caturl))
#     t1 <- gregexpr(data_heads,xmlcode,fixed=TRUE)
#     t2 <- gregexpr(paste0('.',extension,'</tt>'),xmlcode,fixed=TRUE)
#     i1 <- t1[[1]] + nchar(data_heads)
#     i2 <- t2[[1]] + nchar(extension)
#     n1 <- length(i1); n2 <- length(i2); n <- min(c(n1,n2))
#     print(c(n1,n2,n))
#     files <- rep('',n)
#     for (i in 1:n){
#       files[i] <- substr(xmlcode,i1[i],i2[i])
#       t3 <- regexpr("'><tt>",files[i],fixed=TRUE)
#       files[i] <- substr(files[i],1,t3-1)
#     }
#     return(files)
#   }
# }
# # 
# library(esd)
# require("RCurl")
# #result <- getURL("https://thredds.met.no/thredds/catalog/KSS/Klima_i_Norge_2100/seasonal_RCM/catalog.html",
# #                 verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE)
# xmlcode <- getURLContent("https://thredds.met.no/thredds/catalog/KSS/Klima_i_Norge_2100/seasonal_RCM/catalog.html")
# 
# ## Need to extract information on the file names from the XML mess:
# t1 <- gregexpr('catalog.html?dataset=KSS/Klima_i_Norge_2100/',xmlcode,fixed=TRUE)
# t2 <- gregexpr('.nc</tt>',xmlcode,fixed=TRUE)
# i1 <- t1[[1]] + nchar('catalog.html?dataset=KSS/Klima_i_Norge_2100/')
# i2 <- t2[[1]] + 2
# n1 <- length(i1); n2 <- length(i2); n <- min(c(n1,n2))
# print(c(n1,n2,n))
# files <- rep('',n)
# for (i in 1:n) files[i] <- substr(xmlcode,i1[i+1],i2[i])
# files <- sub('seasonal_RCM/','',files)
# for (i in 1:n) {t3 <- regexpr("'><tt>",files[i],fixed=TRUE); files[i] <- substr(files[i],1,t3-1)}
# ## Check the result
# print(files)



