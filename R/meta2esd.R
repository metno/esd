# author Rasmus E. Benestad
# Last-Update 2013-07-23 by Abdelkader Mezghani
 
decryptcn <- function(codes,src="ECAD") {
  n <- length(codes)
  country <- rep("",n)
  for (i in 1:n) {
    country[i]= switch(codes[i],'SE'='Sweden','AT'=,'Austria','BE'='Belgium',
            'HR'='Croatia','CY'='Cypros','CZ'='Chzeck Republic',
            'FI'='Finland','FR'='France','DE'='Germany',
            'IS'='Iceland','RU'='Russia','DK'='Denmark',
            'IE'='Ireland','NL'='the Netherlands','IT'='Italy',
            'NO'='Norway','LV'='Latvia','LT'='Lituania',
            'PT'='Portugal','RO'='Romania','SK'='Slovakia',
            'SI'='Slovenia','ES'='Spain','CH'='Switzerland',
            'RS'='Serbia','EE'='Estonia','MK'='Makedonia',
            'GB'='Great Britain','BA'='Bosnia-Hertsogovina',
            'AL'='Albania','DZ'='Algeria','LU'='Luxemburg',
            'AM'='Armenia','GL'='Greenland','AZ'='Azerbaijan',
            'EG'='Egypt','GR'='Greece','PL'='Poland','IL'='Israel',
            'BY'='Belarus','GE'='Georgia','HU'='Hungary',
            'IQ'='Iraq','KZ'='Khazakstan','LY'='Libya',
            'MD'='Moldova','MO'='Morocco','MT'='Malta',
            'SA'='Saudi Arabia','SY'='Syria','TJ'='Tajikistan',
            'TR'="Turkey",'UA'='Ukraina','UZ'='Uzbekistan',
            'B'='Belgia','FIN'='Finland','FRI'='Faroe Islands',
            'G'='Greenland','IRL'='Ireland','IS'='Iceland',
            'N'='Norway','S'='Sweden')
  }
  invisible(country)
}

meta2esd <- function(silent=TRUE) {
    require(clim.pact)
    ## These data sets will be re-saved in esd and consolidated in terms of parameter names
    data("nacd.meta",envir=environment())
    data("narp.meta",envir=environment())
    ##data("narp.meta",envir=environment())
    ##narp.meta <- NARP
    ##data("nordklim.meta",envir=environment())
    load("nordklim.meta.rda")
    ##load("eca.meta.rda")
    load("ecad.meta.rda") ; ecsn.meta <- ecad.meta
    load("ghcnd.meta.rda")
    load("ghcnm.meta.rda")
    
    n.narp <- length(narp.meta$stnr)
    n.nacd <- length(nacd.meta$station.number)
    n.nork <- length(nordklim.meta$station_id)
  
    n.nacd <- length(nacd.meta$station.number)
    n.ecad <- length(ecad.meta$station_id)
    n.ghcnm <- length(ghcnm.meta$station_id)
    n.ghcnd <- length(ghcnd.meta$station_id)
    ##nacd.lon <- nacd.meta$degE + nacd.meta$minE/60
    nacd.lon <- (nacd.meta$degE + nacd.meta$minE/60)*ifelse(as.character(nacd.meta$E.W)==" E",1,-1)
    nacd.lat <- nacd.meta$degN + nacd.meta$minN/60

    nacd.cn <- strip(as.character(nacd.meta$country))
    nacd.cn[is.element(nacd.cn,'FR')] <- 'FRI'
    ##r.script <- readLines("meta2esd.R")
    
    ## rearrange narp meta
    ele <- as.numeric(ele2param(src='narp')$element)
    narp.ele.len <- length(ele)
    narp.ele <- rep(ele,each=length(narp.meta$stnr))
    narp.stid <- rep(narp.meta$stnr,narp.ele.len)
    narp.loc <- rep(narp.meta$names,narp.ele.len)
    narp.lon <- rep(narp.meta$lons,narp.ele.len)
    narp.lat <- rep(narp.meta$lats,narp.ele.len)
    narp.cntr <- rep(narp.meta$countries,narp.ele.len)
    narp.wmo <-  rep(narp.meta$WMO.number,narp.ele.len)
    n.narp <- length(narp.stid)
    ##  data("observation.meta",envir=environment()
    station.meta <- list(
    	station_id=as.character(c(nacd.meta$station.number,narp.stid,nordklim.meta$station_id,ecad.meta$station_id,ghcnm.meta$station_id,ghcnd.meta$station_id)),
    	location=c(nacd.meta$location,narp.loc,nordklim.meta$location,ecad.meta$location,ghcnm.meta$location,ghcnd.meta$location),
	country=c(toupper(decryptcn(nacd.cn)),toupper(narp.cntr),nordklim.meta$country,ecad.meta$country,ghcnm.meta$country,ghcnd.meta$country),
    	longitude=c(nacd.lon,narp.lon,nordklim.meta$longitude,ecad.meta$lon,ghcnm.meta$longitude,ghcnd.meta$longitude),
    	latitude=c(nacd.lat,narp.lat,nordklim.meta$latitude,ecad.meta$lat,ghcnm.meta$latitude,ghcnd.meta$latitude),
    	altitude=c(nacd.meta$alt,rep(NA,n.narp),nordklim.meta$altitude,ecad.meta$alt,ghcnm.meta$altitude,ghcnd.meta$altitude),
	element=c(nacd.meta$element,narp.ele,nordklim.meta$element,ecad.meta$element,rep(NA,n.ghcnm),rep(NA,n.ghcnd)),
	start=c(nacd.meta$start,rep(1870,n.narp),nordklim.meta$start,ecad.meta$start,ghcnm.meta$start,ghcnd.meta$start),
	end=c(rep(NA,n.nacd),rep(2000,n.narp),nordklim.meta$end,ecad.meta$end,ghcnm.meta$end,ghcnd.meta$end),    	
	source=c(rep("NACD",n.nacd),rep("NARP",n.narp),rep("NORDKLIM",n.nork),rep("ECAD",n.ecad),rep("GHCNM",n.ghcnm),rep("GHCND",n.ghcnd)),
	wmo=c(nacd.meta$wmo.number,narp.wmo,rep(NA,n.nork),rep(NA,n.ecad),rep(NA,n.ghcnm),rep(NA,n.ghcnd)),  
    	quality=c(nacd.meta$quality,rep(NA,n.narp),rep(NA,n.nork),rep(NA,n.ecad),rep(NA,n.ghcnm),rep(NA,n.ghcnd))
)

  if (!silent) print(str(station.meta))
  if (!silent) print(table(station.meta$source))
  attr(station.meta,'history') <- c('meta2esd.R - data taken from the clim.pact package and consilidated for NACD and NARP',"nordklima.meta.rda","ecad.meta.rda","ghcnd.meta.rda","ghcnm.meta.rda")
  #attr(station.meta,'R-script') <- r.script
  attr(station.meta,'date') <- date
  attr(station.meta,'call') <- match.call()
  attr(station.meta,'author') <- 'R.E. Benestad & A. Mezghani'
  attr(station.meta,'URLs') = c(
         "www.dmi.dk/dmi/sr96-1.pdf",
         "http://www.norden.org/en/publications/publikationer/2005-450",
         "http://www.smhi.se/hfa_coord/nordklim/",
         "http://eca.knmi.nl/")
  save(file="station.meta.rda",station.meta)
return(station.meta)
}
