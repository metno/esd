
# This function adds a stamp in the history of x
# with 'sys.call', 'date()', and 'src (source)

history.stamp <-function(x=NULL,y=NULL,...) UseMethod("history.stamp")

history.stamp.default <- function(x=NULL,y=NULL) {
#  if (!is.null(attr(x,'source')))
#    src <- attr(x,'source') else
#    src <- "unknown source"
  si <- list(R.version=sessionInfo()$R.version$version.string,
             esd.version=paste(sessionInfo()$otherPkgs$esd$Package,
                               sessionInfo()$otherPkgs$esd$Version,sep="_"),
             platform=sessionInfo()$platform)
  if (is.null(x)) x <- 0
  if (!is.null(attr(x,'history'))) {
    history <- attr(x,'history')
    call <- c(sys.call(sys.parent(n = 1)),history$call)
    sessioninfo <- c(si,history$sessioninfo)
    timestamp <- c(date(),history$timestamp)
    newhistory <- list(call=call,timestamp=timestamp,
                       session=sessioninfo)
    
  } else {
    newhistory <- list(call=sys.call(sys.parent(n = 1)),
                    timestamp=date(),
                    sessioninfo=si)
  }
  #print(sys.call(sys.parent(n = 1)))
  return(newhistory)
}



