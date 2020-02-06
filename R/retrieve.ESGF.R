#' Retrieve CMIP data directly from the Earth System Grid Federation (ESGF)
#' https://earthsystemcog.org/projects/cog/esgf_search_restful_api

#' @aliases retrieve.ESGF

#' @import ncdf4
#'
#' @seealso retrieve
#'
#' @param param Name of parameter
#' @param verbose Logical value defaulting to FALSE. If FALSE, do not display
#' comments (silent mode). If TRUE, displays extra information on progress.
#' @param url The base URL for ESGF
#' @param CMIP Which generation of CMIP
#' @param MIP Name of activity
#' @param inst Institution ID
#' @param model Name of climate model
#' @param scen Name of scenario
#' @param ID The name of the run
#' @param type The table ID
#' @param grid Grid ID
#' @param version Version
#' @param period The time segment of the file
#' @return A "zoo" "field" object with additional attributes used for further
#' processing.
#'
#' @examples 
#' retrieve.ESGF(lon=c(-30,40),lat=c(50,70),period='205001-210012') -> X
#' map(X)

#' @export retrieve.ESGF
retrieve.ESGF <- function(param='tas',verbose=FALSE,...,
                          url='https://esgf-data3.ceda.ac.uk/thredds/dodsC',CMIP='esg_cmip6/CMIP6',
                          MIP='ScenarioMIP',inst='MOHC',model='UKESM1-0-LL',scen='ssp585',
                          ID='r2i1p1f2',type='Amon',grid='gn',version='v20190507',period='2015-2100') {
  fname <- paste0(paste(param,type,model,scen,ID,grid,period,sep='_'),'.nc')
  opendap <- paste(url,CMIP,MIP,inst,model,scen,ID,type,param,grid,version,fname,sep='/')
  if (verbose) print(opendap)
  x <- try(retrieve(ncfile=opendap,verbose=verbose,...))
  invisible(x)
  
}
