# Estimate the Kling-Gupta model efficiency
#
# This function returns the 2012 Kling-Gupta model efficiency using CV for the variability of x in predicting y.
# An efficiency of one corresponds to a perfect match, while the lowest score is -infinity.
# x: is a 1-D zoo object corresponding to the prediction
# y: is a 1-D zoo object corresponding to the reference
# scc, sv, sm: are weights within (0,1) allowing an informed change from the default weight of 1.
# author: Helene B. Erlandsen
# details: This functions return the dimensionless, multi-objective Kling-Gupta efficiency metric,
# KGE, which is based on the correlation, variability ratio, and mean ratio
# between the prediction and reference objects and is a modified verision of the
# Nash-Sutcliffe efficiency (which again corresponds to the unbiased \eqn{R^2}).
# The weight of the sub-metrics (scc, sv, sm) may be adjusted from 1 in the arguments
# kge= 1 - sqrt( (scc(cc-1)^2 + sv(CV_m/CV_o - 1)^2 + sm(mean(mod)/mean(obs)-1)^2 )) 
# It is calculated for pairwise complete observations
# references: Gupta HV, Kling H, Yilmaz KK and Martinez GF (2009).,
# "Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling.", Journal of Hydrology, 377(1), pp. 80–91.
#
# Warning:
#  The mean ratio may not be approriate for Celcius (can blow up around 0)
#  or other variables that may have a mean within (-1,1). Though the metric is dimensionless, keep in mind that
#  variables with highere avarages (e.g. temperature in Kelvin, longwave radiation) will be less penalized
#  for mean devaitons than variables with lower averages (temp. in Celsius, winter SW radiation)



#' Model efficiency evaluation
#' 
#' The KGE function returns the 2012 Kling-Gupta model efficiency using CV for
#' the variability of x in predicting y.  The NSE function returns the
#' dimensionless Nash Sutcliffe model efficiency, evalutaing x as a prediction
#' of y. An efficiency of one corresponds to a perfect match, while the lowest
#' score is -infinity.
#' 
#' The NSE function returns the multi-objective Nash-Sutcliffe efficiency
#' metric, which corresponds to the unbiased \eqn{R^2}. The metric is scaled by
#' the observed variance.  The metric is given as nse = 1 -
#' mean_square_error/observed_variance.  It is calculated for pairwise complete
#' observations
#' 
#' Warning: Keep in mind that in regions with higher variance (e.g.
#' seasonality) equal absolute deviations will be penalized less than for
#' regions with a lower variance.
#' 
#' The KGE function returns the dimensionless, multi-objective Kling-Gupta
#' efficiency metric, which is based on the correlation, variability ratio, and
#' mean ratio between the prediction and reference objects and is a modified
#' version of the Nash-Sutcliffe efficiency.  The weight of the sub-metrics
#' (scc, sv, sm) may be adjusted from 1 in the arguments kge = 1 - sqrt(
#' (scc(cc-1)^2 + sv(CV_m/CV_o - 1)^2 + sm(mean(mod)/mean(obs)-1)^2 )).
#' 
#' Warning: The mean ratio may not be approriate for Celcius (can blow up
#' around 0) or other variables that may have a mean within (-1,1). Though the
#' KGE etric is dimensionless, keep in mind that variables with higher averages
#' (e.g. temperature in Kelvin, longwave radiation) will be less penalized for
#' mean deviations than variables with lower averages (temp. in Celsius, winter
#' SW radiation).
#' 
#' @aliases KGE Kling-Gupta nse NSE Nash-Sutcliffe
#' @param x a 1-D zoo object corresponding to the prediction
#' @param y a 1-D zoo object corresponding to the reference
#' @param return_all Directional resolution in degrees
#' @param scc weight of (0-1)
#' @param sv weight of (0-1)
#' @param sm weight of (0-1)
#' @author H.B. Erlandsen
#' @references Gupta HV, Kling H, Yilmaz KK and Martinez GF (2009).,
#' "Decomposition of the mean squared error and NSE performance criteria:
#' Implications for improving hydrological modelling.", Journal of Hydrology,
#' 377(1), 80-91. Nash JE and Sutcliffe JV (1970), River flow forecasting
#' through conceptual models part I-A discussion of principles.???, Journal of
#' hydrology, 10(3), pp. 282-290. Gupta HV, Kling H, Yilmaz KK and Martinez GF
#' (2009)., "Decomposition of the mean squared error and NSE performance
#' criteria: Implications for improving hydrological modelling.", Journal of
#' Hydrology, 377(1), pp. 80-91.
#' @keywords "validation" "model-efficiency"
#' @examples
#' 
#'   \dontrun{
#'     data('Oslo')
#'     x <- subset(Oslo,it=!is.na(Oslo))
#'     mydata <- data.frame(x=x,t=index(x))
#'     fit <- lm(x ~ t, data=mydata)
#'     y <- zoo(fitted(fit,index(x)),order.by=index(x))
#'     KGE(x,y)
#'   }
#' 
#' @export kge
kge <- function (x,y,return_all=FALSE, scc=1., sv=1.,sm=1.) {
  d <- merge.zoo(x,y, all=F, fill=0)
  x <- d[,1]
  y <- d[,2]
  cc     <- cor(x, y)
  sd_ratio <- sd(x)/ sd(y)
  mean_ratio  <- mean(x)/mean(y)
  cv_ratio <- ((sd(x)/mean(x))/(sd(y)/mean(y)))
  if ((abs(mean(y)) < 1) || (abs(mean(x))) < 1) {
    print("Mean within (-1,1) may lead to blowing up the score, please rescale.")}
  eds <- sqrt( scc*(cc-1)^2 + sv*(cv_ratio-1)^2 + sm*(mean_ratio-1)^2 )
  z <- 1-eds
  kgel <- list(kge = z, r=cc, cv_ratio=cv_ratio, mean_ratio=mean_ratio)
  if(return_all) return(kgel) else return(z)
}
KGE <- function(x,y) return(kge(x,y))


# Estimate the Nash-Sutcliffe model efficiency
#
# This function returns the dimensionless Nash Sutcliffe model efficiency, evalutaing x as a prediction of y.
# An efficiency of one corresponds to a perfect match, while the lowest score is -infinity.
# x: is a 1-D zoo object corresponding to the prediction
# y: is a 1-D zoo object corresponding to the reference
# author: Helene B. Erlandsen
# details: This functions return the multi-objective Nash-Sutcliffe efficiency metric
# NSE, which corresponds to the unbiased \eqn{R^2}. The metric is scaled by the observed variance. The metric is given as
# nse= 1 - mean_square_error/observed_variance
# and varies thus form 1 to minus infinity. It is calculated for pairwise complete observations
# references: Nash JE and Sutcliffe JV (1970),
# River flow forecasting through conceptual models part I-A discussion of principles.”,
#  Journal of hydrology, 10(3), pp. 282–290.,
#
#  Gupta HV, Kling H, Yilmaz KK and Martinez GF (2009).,
# "Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling.", Journal of Hydrology, 377(1), pp. 80–91.
# Warning: Keep in mind that in regions with higher variance (e.g. seasonality) equal absolute deviations will be penalized less than for regions with a lower variance

## Estimate the Nash-Sutcliffe efficiency
## nse=1 -mse/var_O
nse <- function(x,y) {
  MSE <- mean((x - y)^2, na.rm=TRUE)
  var_O <- var(x, y = NULL, na.rm = TRUE, use='pairwise.complete.obs')
  z <- 1-MSE/var_O
  return(z)
}
NSE <- function(x,y) return(nse(x,y))

