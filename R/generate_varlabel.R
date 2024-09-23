#' Generate label for maps and plots
#' 
#' Create a label containing the variable name and unit and optional additional information to be used in maps and plots.
#' 
#' @param x The object to be plotted. If x is provided and the variable and unit names are not explicitly stated, this information is extracted from the attributes of x.
#' @param variabel Variable name. If provided, variable is used rather than the attr(x, 'variable')
#' @param unit Unit name. If provided, unit is used rather than the attr(x, 'unit')
#' @param method Method. If provided, this information is added to the label
#' @param period Time period. If provided, this information is added to the label
#' @param verbose If TRUE, print information about the progress
#'
#' @return An expression to be used as a label in maps and plots
#' 
#' @keywords label
#' @examples
#' 
#' data(bjornholt)
#' label <- generate_varlabel(bjornholt)
#' print(label)
#' 
#' @export
generate_varlabel <- function(x=NULL, variable=NULL, unit=NULL, method=NULL, 
                           period=NULL, verbose=FALSE) {
  if(verbose) print("generate_varlabel")
  if(is.null(unit)) unit <- attr(x,'unit')
  if(is.null(variable)) variable <- varid(x)
  
  if(!is.null(unit)) unit <- as.character(unit)
  
  if ( (unit=="degC") | (unit=="deg C") | (unit=="degree C") | (unit=="degree Celsius"))
    unit <- "degree*C"
  if (unit=="%") unit <- "'%'"
  if(!is.null(variable)) {
    if ( (tolower(variable)=="t(2m)") | (tolower(variable)=="t2m") |
         (tolower(variable)=="2t") )
      variable <- "T[2*m]"
  }
  if (verbose) print(paste(variable, unit, ' -> varlabel'))

  if(!is.null(variable)) {
    varlabel <- try(eval(parse(
      text=paste('expression(', gsub(" ", "~", variable), 
                 " *~(", gsub(" ", "~", unit), "))", sep="")))) 
  } else {
    varlabel <- NULL
  }

  if (inherits(varlabel,'try-error')) varlabel <- NULL
  if (verbose) {print(varlabel); print(src(x))}

  if(!is.null(varlabel)) label <- paste(varlabel,'*')
  label <- as.expression(parse(text=paste(label,'phantom(0) - phantom(0)')))
  if (inherits(label,'try-error')) label <- ''
  
  if (!is.null(method)) {
    label2 <- try(parse(text=paste(label,'*',as.expression(method))))
    if (!inherits(label2,'try-error')) label <- label2
  }
  
  if (!is.null(period)) {
    label2 <- try(parse(text=paste(label,'*',as.expression(period))))
    if (!inherits(label2,'try-error')) label <- label2
  }

  return(label)
}
  