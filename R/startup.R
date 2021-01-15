# ## Start up - to fix the problem with S3 methods for R-versions after 4.0
# ## This script is run when the package loads.
# print('You have made an excellent choise to use the esd-package')
# print('--------------------------------------------------------')
# #if ( (file.exists('~/git/esd/R')) & (if(getRversion() >= "3.6.0")) ) { 
# system('grep UseMethod ~/git/esd/R/*.R',intern=TRUE) -> esds3methods
# ns3 <- length(esds3methods)
# for (i in 1:ns3) { 
#   method <- substr(esds3methods[i],regexec('UseMethod',esds3methods[i])[[1]],nchar(esds3methods[i]))
#   method <- sub('UseMethod','',method)
#   method <- sub('(','',method,fixed=TRUE)
#   method <- sub(')','',method,fixed=TRUE)
#   method <- gsub('\\','',method,fixed=TRUE)
#   method <- gsub('"','',method,fixed=TRUE)
#   print(method)
#   #eval(parse(text=paste('.S3method(',method,',',)))
#   }
# #}