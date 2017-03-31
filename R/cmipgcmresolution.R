## GCMresolution
## https://climatedataguide.ucar.edu/climate-model-evaluation/common-spectral-model-grid-resolutions
## Truncation 	lat x lon 	km at Eq 	deg at Eq
## T21 	32x64 	625 	5.625
## T42 	64x128 	310 	2.8125
## T62 	94x192 	210 	1.875
## T63 	96x192 	210 	1.875
## T85 	128x256 	155 1.4
## T106 	160x320 	125 	1.125
## T255 	256x512 	60 	0.703125
## T382 	576x1152 	38 	0.313
## T799 	800x1600 	25 	0.225

##http://www.ecmwf.int/en/forecasts/documentation-and-support/gaussian_n640

cmipgcmresolution <- function(what='deg') {
  data("IPCC.AR5.Table.9.A.1")
  tX <- c('T21',   'T31',  'T42',     'T62',    'T63',   'T85',      'T106',   'T255',     'T382',       'T799',      'T127',
          'N96',   'T159', 'C48',     'N48',    'T126',  NA,         'C360',   'C180',     'R42',        'T959',      'M45', 
           'unknown')
  dX <- c(5.625,    3.75,   2.8,      1.875,    1.875,    1.4,        1.125,    0.703,      0.313,       0.225,        0.95,
          1.875,    1.125,  1.8,      3.75,      0.9375,  1.9,        0.225,    0.449,      2.81,        0.1875,       2.5,
          1.25)   
  dn <- c('32x64', '96x48', '128x64', '94x192', '96x192', '128x256', '160x320', '256x512',  '576x1152',  '800x1600',   NA,
          '192x145','320x160', NA,    '3.75x2.5', NA,     '96x95',    NA,       NA,         '128x108',  '1920x960',    NA,
          '144x143') 
  dr <- c(625,      NA,      310,      210,      210,      155,       125,      60,         38,          25,          106,
          209,      125.23,  200,      417,      104,      212,       25,       50,         313,         20,          278,
          139)
  
  
  x <- as.character(IPCC.AR5.Table.9.A.1$Horizontal.Grid)
  n <- length(x)
  x <- gsub('TL','T',x,fixed=TRUE)
  x <- gsub(' deg latitude','', x,fixed=TRUE)
  x <- gsub(' deg longitude','', x,fixed=TRUE)
  x <- gsub(' deg in latitude','', x,fixed=TRUE)
  x <- gsub(' deg in longitude','', x,fixed=TRUE)
  x <- gsub('Nominally ','', x,fixed=TRUE)
  x <- gsub('Finite Volume ','', x,fixed=TRUE)
  y <- rep(NA,n)
  
  ## Search on T/N specs.
  for (i in 1:length(tX)) {
    im <- grep(tX[i],x)
    if (what=='deg') y[im] <- dX[i] else if (what=='km') y[im] <- dr[i] else y[im] <- dn[i]
    x[im] <- NA
  }
  for (i in 1:length(tX)) {
    im <- grep(dn[i],x)
    if (what=='deg') y[im] <- dX[i] else if (what=='km') y[im] <- dr[i] else y[im] <- dn[i]
    x[im] <- NA
  }
  xs <- grep('x',x)
  for (i in 1:length(xs)) {
    y[xs][i] <- as.numeric(substr(x[xs][i],1,regexpr('x',x[xs][i])-1))
    x[xs][i] <- NA
  }
  xs <- grep(' deg',x)
  for (i in 1:length(xs)) {
    y[xs][i] <- as.numeric(substr(x[xs][i],1,regexpr(' deg',x[xs][i])-1))
    x[xs][i] <- NA
  }
  xs <- grep(',',x)
  for (i in 1:length(xs)) {
    y[xs][i] <- as.numeric(substr(x[xs][i],1,regexpr(',',x[xs][i])-1))
    x[xs][i] <- NA
  }
  #print(x)
  return(y)
}
