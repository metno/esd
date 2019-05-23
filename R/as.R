
# KMP 2019-05-14: seems obsolete and does not return a ds object
#as.ds <- function(x,verbose=FALSE) {
#  if(verbose) print(as.ds)
#  y <- zoo(x,order.by=index(x))
#  attr(y,'location') <- attr(x,'location')
#  attr(y,'variable') <- attr(x,'variable')
#  attr(y,'unit') <- attr(x,'unit')
#  attr(y,'longitude') <- attr(x,'longitude')
#  attr(y,'latitude') <- attr(x,'latitude')
#  attr(y,'altitude') <- attr(x,'altitude')
#  attr(y,'country') <- attr(x,'country')
#  attr(y,'longname') <- attr(x,'longname')
#  attr(y,'station_id') <- attr(x,'station_id')
#  attr(y,'quality') <- attr(x,'quality') 
#  attr(y,'calendar') <- attr(x,'calendar')
#  attr(y,'source') <- attr(x,'source')
#  attr(y,'URL') <- attr(x,'URL')
#  #attr(y,'history') <- attr(x,'history')
#  #attr(y,'date-stamp') <- attr(x,'date-stamp')
#  attr(y,'type') <- attr(x,'type')
#  attr(y,'aspect') <- attr(x,'aspect')
#  attr(y,'reference') <- attr(x,'reference')
#  attr(y,'info') <- attr(x,'info')
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- c("station","month","zoo")
#  return(y)
#}


