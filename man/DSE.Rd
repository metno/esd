%\name{DSE}
%\alias{DSE}
%\alias{DSE.metno}
%\alias{DSE.ecad}
%\alias{DSE.ghcnd}
%
%\title{Downscale monthly climate variables and parameters for several stations.}
%\description{Performs a complete downscaling job based on \link{DSensemble}
  for several stations and save the downscaled results locally. DSE can
  be used to downscale climate variables, parameters or statistics for
  various stations. Each station is downscaled separately for the whole rcp/gcm %ensemble and stored the results in a
%  seperate rda file. DSE() will also generate two files: an inventory file
%  containing the meta data and some statisitcs for the successfully
%  downscaled stations and log file containing encountered errors and a
%  list of stations that have not been downscaled.}
%\usage{

%DSE.default(cntr, src, param, FUN, lon, lat, path, rcp, biascorrect,reanalysis, 
%    email, save, out.dir, force, update, ...) 
%DSE.metno(x,...) # Downscale METNO stations
%DSE.ecad(x,...)
%DSE.ghcnd(x,...)
%meta.dse(path) # update the inventory file with additional downscaled results.
%}
%\arguments{
%  \item{stid}{A string of characters as an identifier of the weather/climate st%ation.}
%  \item{param}{Parameter or element type or variable identifier. There are seve%ral core parameters or elements as well as a number of additional parameters. T%he parameters or elements are:
%    prcp, pr     = Precipitation (mm)
%    tas, tavg    = 2m-surface temperature (in degrees Celcius)
%    tmax, tasmax = Maximum temperature (in degrees Celcius)
%    tmin, tasmin = Minimum temperature (in degrees Celcius)
%  }
%  \item{FUN}{A mathematical transformation of the parameter e.g. "mean",
%    "sd", "r1", etc ...} 
%  \item{rcp}{A character string corresponding to the name of the RCP
%    experiment. Possible values for CMIP5 are rcp20, rcp45, rcp65,
%    rcp85. For the CMIP3, they are: sresa1b, sresa2.}
%  \item{src}{Source: limit the downscaling to a specific data set ("NARP",
%    "NACD", "NORDKLIMA", "GHCNM", "METNOM", "ECAD", "GHCND" and "METNOD")}
%  \item{stid}{A string of characters as an identifier of the weather/climate st%ation.}
%  \item{lon}{A numeric vector containing the range of longitude values
%    in the form of c(lon.min,lon.max) which corresponds to the predictor
%    longitude extension}
%  \item{lat}{A numeric vector containing the range of latitude values in
%    the form of c(lat.min,lat.max) which corresponds to the predictor
%    latitude domain extension}
%  
%  \item{biascorrect}{Logical value. If TRUE, a bias correction method is
%    applied} 
%  \item{save}{Logical value. If TRUE, the results are save in rda files
%    and stored locally in directory specified by path.}
%  \item{out.dir}{A character string of the full path specifying the name of the output
%    directory. If the directory does not exist, il will be created
%    automaitcally.}
%  \item{cntr}{A string or a vector of strings of the full name of the
%    country. Select the stations from a specified country or a
%    set of countries.}
% \item{force}{Logical value. If TRUE, the DSE is performed from scratch
%   and all existing DSE files will be re-written. If FALSE, DSE will skip
%   already downscaled stations and do the downscaling only for non
%   existing file stations.} 
% \item{path}{A character string corresponding to the full path name of
%   the CMIP3/5 results.}
% \item{reanalysis}{A character string with the full path file name of
%   the reanlysis to be used. Only netcdf files are accepted.}        update=
% \item{email}{A character string with an email address. If specified, a
%   notification email is sent when the downscaling is complete. The
%   email text message contains a list of non-downscaled stations
%   because of errors. If any error, please forward the email to abdelkadem@met.no}
% \item{update}{A logical value. If TRUE, \link{meta.dse} is called and
%   the inventory file is updated if additional DSE files are copied to
%   output directory.}
% \item{...}{Additional arguments to be passed in the function}
%}
%\value{A "zoo" "DSensemble" object}
%
%\seealso{
%  \code{\link{DSensemble}}
%}
%\author{A. Mezghani, MET Norway}
%\examples{
%\dontrun{
%DSE.metno(cntr = "NORWAY", src = "METNO", param = "t2m",nmin=100, 
%    lon = c(-10, 10), lat = c(-10, 10), FUN = "mean", path = "CMIP5.monthly", 
%    rcp = "rcp45", biascorrect = TRUE, reanalysis = "ERA40_t2m_mon.nc", 
%    email = NULL, save = TRUE, force = FALSE, verbose = FALSE, 
%    out.dir = "dse.100", update = FALSE) 
%}
%}
%\keyword{parameter , element}