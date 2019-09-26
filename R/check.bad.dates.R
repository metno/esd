#' Check dates for impossibilities
#'
#' @param yyyy year(s)
#' @param mm month(s)
#' @param dd day(s)
#'
#' @export
check.bad.dates <- function(yyyy,mm,dd) {
  yyyy <- as.numeric(yyyy); mm <- as.numeric(mm); dd <- as.numeric(dd)
  ndpy <- as.numeric(table(yyyy))  ## Number of days per year
  bad.dates <- (max(ndpy) > 366) | (min(mm) < 1) | (max(mm) > 12) |
                (min(dd) < 0) | (max(dd)>31)
  if (is.na(bad.dates)) bad.dates <- TRUE
  return(bad.dates)
}