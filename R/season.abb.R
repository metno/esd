#' Season abbreviation
#'
#' @return a list with season abbreviations and their corresponding months: list("annual"=1:12, "djf"=c(12,1,2),...)
#' 
#' @export season.abb
season.abb <- function() {
  season.abb <- c('annual','djf','jfm','fma','mam','amj',
                  'mjj','jja','jas','aso',
                  'son','ond','ndj','ondjfm','amjjas',
                  'ndjf','jjas','mjjas',
                  'djfm','djfma','ndjfma',
                  'ndjfmam','ondjfma')
  season<-list(1:12,c(12,1,2),1:3,2:4,3:5,4:6,5:7,6:8,7:9,8:10,9:11,10:12,
               c(11,12,1),c(10:12,1:3),4:9,c(11,12,1,2),6:9,5:9,
               c(12,1,2,3),c(12,1,2,3,4),c(11,12,1,2,3,4),c(11,12,1,2,3,4,5),
               c(10,11,12,1,2,3,4))
  season.abb <- c(season.abb,toupper(season.abb))
  season <- rep(season,2)
  names(season) <- season.abb
  return(season)
}