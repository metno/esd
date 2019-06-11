#' Help and assistance
#'
#' @aliases ABC4ESD esd.tips esd.issues esd.tips esd.issues downscaling.about element.kode
#'
#' @export
manual <- function(url='https://ndownloader.figshare.com/files/2126237',browser='firefox')
  system(paste(browser,url))

#' @export
ABC4ESD <- function(url='https://ndownloader.figshare.com/files/9703462',browser='firefox')
  system(paste(browser,url))

#' @export
esd.tips <- function(url='https://github.com/metno/esd/wiki',browser='firefox')
  system(paste(browser,url))

#' @export
esd.issues <- function(url='https://github.com/metno/esd/issues',browser='firefox')
  system(paste(browser,url))

#' @export
downscaling.about <- function(url='http://climatescience.oxfordre.com/view/10.1093/acrefore/9780190228620.001.0001/acrefore-9780190228620-e-27',browser='firefox')
  system(paste(browser,url))

#' @export
element.kode <- function(url='http://klapp/metnopub/production/metno?re=26&ct=text/plain&del=semicolon&tab=T_ELEM_DIURNAL',browser='firefox')
  system(paste(browser,url)) 

