## Frost API to read METNO data

test.API <- function(keyfile='~/.FrostAPI.key') {
  ## This test function works
  if (file.exists(keyfile)) {
    if (verbose) print(paste('Read client ID from',keyfile))
    frostID <- readLines(keyfile) 
  } else { 
    if (verbose) print('Generate new client ID')  
    system(paste(browser,'https://frost.met.no/auth/newclientid.html'))
    frostID <- rep("",2)
    frostID[1] <- readline('Please give me the first key:')
    frostID[2] <- readline('Please give me the second key:')
    writeLines(frostID,con=keyfile)
  }
  url <- paste0('https://',frostID[1],':',frostID[2],
                '@frost.met.no/observations/v0.csv?sources=sn18700',
                '&referencetime=2019-07-10%2F2019-07-11&elements=air_temperature')
  ## The output is a table with data
  View(readLines(url))
}

stripblanks <- function(x) {
  while(substr(x,1,1)==' ') x <- substr(x,2,nchar(x))
  while(substr(x,nchar(x),nchar(x))==' ') x <- substr(x,1,nchar(x)-1)
  return(x)
}

metafrostAPI <- function(keyfile='~/.FrostAPI.key',verbose=FALSE,
                         fields='id,name,masl,county,municipality,wmoid,geometry,type',
                         url='frost.met.no/sources/v0.jsonld?country=NO') {
  
  require(jsonlite)
  t0 <- Sys.time()
  
  if (file.exists(keyfile)) {
    if (verbose) print(paste('Read client ID from',keyfile))
    frostID <- readLines(keyfile) 
  } else { 
    if (verbose) print('Generate new client ID')  
    system(paste(browser,'https://frost.met.no/auth/newclientid.html'))
    frostID <- rep("",2)
    frostID[1] <- readline('Please give me the first key:')
    frostID[2] <- readline('Please give me the second key:')
    writeLines(frostID,con=keyfile)
  }
  url <- paste0('https://',frostID[1],':',frostID[2],'@',url,'&fields=',fields)
  if (verbose) print(url)
  
  # ## Read the station metadata - it's in JSON format
  # json <- readLines(url)
  # ## Distill useful information and discard all the garbage
  # il <- grep('"name"',json)
  # ia <- grep('"masl"',json)
  # ic <- grep('"county"',json)
  # im <- grep('"municipality"',json)
  # #is <- grep('""',json)
  # id <- grep('"id"',json)
  # ik <- grep('"coordinates"',json)
  # if (verbose) print(paste('IL:',length(il),' IA:',length(ia),' IC:',length(ic),
  #                    ' IM: ',length(im),' ID:',length(id),' IK:',length(ik)))
  # browser()
  # locs <- stripblanks(sub(':','',gsub('"','',gsub('"name"','',json[il]))))
  # alt <- stripblanks(sub(':','',gsub('"','',gsub('"masl"','',json[ia]))))
  # fylke <- stripblanks(sub(':','',gsub('"','',gsub('"county"','',json[ic]))))
  # kommune <- stripblanks(sub(':','',gsub('"','',gsub('"municipality"','',
  #                                                    json[im]))))
  # # start <- stripblanks(sub(':','',gsub('"','',gsub('"validFrom"','',
  # #                                                      json[is]))))
  # stid <- stripblanks(sub('SN','',sub(':','',gsub('"','',
  #                                                 gsub('"id"','',json[id])))))
  # 
  # koords <- stripblanks(sub(':','',gsub('"','',gsub('"coordinates"','',
  #                                                   json[ik]))))
  # koords <- sub('],','',sub('[ ','',koords,fixed=TRUE),fixed=TRUE)
  # n <- length(koords)
  # end <- rep(NA,n)
  # ele <- rep(NA,n)
  # varid <- rep(NA,n)
  # koords <- as.numeric(scan(text=koords,sep=','))
  # print(summary(koords))
  # koords <- matrix(data=c(koords),nrow=2,ncol=n)
  # if (verbose) print(table(fylke))
  # if (verbose) print(koords)
  # if (verbose) {par(bty='n');plot(t(koords),xlab='',ylab=''); grid()}
  # if (verbose) View(cbind(locs,stid,fylke,kommune,alt,t(koords)))
  
  xs <- try(fromJSON(URLencode(url),flatten=TRUE))
  if (class(xs) != 'try-error') {
    print("Data retrieved from frost.met.no!")
    data <- xs$data
    if (verbose) View(data)
    n <- length(data$id)
    start <- rep(NA,n)
    end <- rep(NA,n)
    ele <- rep(NA,n)
    varid <- rep(NA,n)
    lons <- varid; lats <- varid
    for (i in 1:n) {
      try <- lapply(data$geometry.coordinates,function(x) x[2])[[i]]
      if (length(try)==1) lons[i] <- try 
      try <- lapply(data$geometry.coordinates,function(x) x[1])[[i]]
      if (length(try)==1) lats[i] <- try
    }
    alts <- as.numeric(data$masl)
    X <- data.frame(station_id=data$id,location=data$name,country=rep('Norway',n),
                    longitude=lons,latitude=lats,altitude=alts,
                    element=ele,start=start,end=end,source=rep('METNO',n),
                    wmo=data$wmoId,quality=rep(NA,n),variable=varid,
                    fylke=data$county,kommune=data$municipality)
    class(X) <- c("stationmeta", "data.frame")
  } else {
    print("Error: the data retrieval was not successful!")
    X <- NULL
  }
  print(paste('time of request was',round(Sys.time()-t0),'seconds'))
  invisible(X)
}


frostAPI <- function(param='mean(air_temperature P1D)',stid=18700,
                     type='observations',keyfile='~/.FrostAPI.key',it = NULL,
                     browser='firefox',verbose=FALSE,level='2',
                     fields='referenceTime%2Cvalue') {
  ## Use saved client ID
  if (file.exists(keyfile)) {
    if (verbose) print(paste('Read client ID from',keyfile))
    frostID <- readLines(keyfile) 
  } else { 
    if (verbose) print('Generate new client ID')  
    system(paste(browser,'https://frost.met.no/auth/newclientid.html'))
    frostID <- rep("",2)
    frostID[1] <- readline('Please give me the first key:')
    frostID[2] <- readline('Please give me the second key:')
    writeLines(frostID,con=keyfile)
  }
  ## Build the URL
  element <- switch(param,'t2m'='air_temperature')
  stid <- paste0('SN',stid,collapse=',')
  print(param)
  param <- gsub('(','%28',param,fixed=TRUE) 
  param <- gsub(')','%29',param,fixed=TRUE) 
  param <- gsub(' ','%20',param)
  #URLencode()
  
  if (is.null(it)) it <- paste0('1960-01-01/',format(Sys.time()-24*3600,"%Y-%m-%d"))
  it <- sub('/','%2F',it)
  url <- paste0('https://',frostID[1],':',frostID[2],'@frost.met.no/',type,'/v0.csv?sources=',
                stid,'&referencetime=',it,'&elements=',param,
                '&fields=',fields,collapse=NULL)
  if (param=='mean(air_temperature P1D)') url <- paste0(url,'&level=',level)
  if (verbose) print(url)
  varid <- switch(param,
                  'mean(air_temperature P1D)'='t2m')
  y <- read.table(url,sep=',',header=TRUE)
  t <- as.character(y[,1])
  print('fix the string')
  t <- sub('T',' ',t); t <- sub('.000Z','',t)
  y <- zoo(as.numeric(as.character(y[,2])),order.by=as.POSIXlt(t))
  y <- as.station(y,stid=stid,param=varid)
  #print(y)
  plot(y)
  invisible(y)
}