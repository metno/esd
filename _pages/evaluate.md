---
layout: home
toc: true
permalink: /evaluate/
---

# Evaluation, assessment & validation

<details markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

## Central limit theorem
The central limit theorem states that the distribution of the sum or average of a large number of independent, identically distributed variables (iid) will be approximately normal, regardless
of the underlying distribution.
(src. http://www.math.uah.edu/stat/sample/CLT.html)

## Huth’s dilemma
Huth’s dilemma is the situation where ESD is calibrated on short-term (fast) fluctuation but is not capable of predicting long-term (slow) changes. This is a problem if there are different
processes responsible for variations at different time scales. One test to provide a diagnose on whether this is the case is to calibrate the model on de-trended data and then use the original
field as input to see if it is able to predict the long-term trend over the calibration period.
In ‘esd’, the default process de-trends the data before calibration. The downscaled results are derived with the original field, however, and trends are compared against the original ob-
servations. Furthermore, the question of stationary signal should be examined in the context of climate model results, and the past changes seen in the observations should be compared with
corresponding downscaled results from GCM historic runs. Ensembles of GCMs are particularly suitable for assessing the ESD model’s ability to capture the long-term changes.

## Non-stationarity check
The non-stationarity check in ‘esd’ treats the GCM results as pseudo-reality, where results interpolated to the coordinates as a given station is used as the predictand. The same time
interval as the reanalysis data used for calibrating the model is extracted from the GCM for representing the calibration predictors. The model trained on the subset of GCM results is then
used to predict future values for the predictand, and the predictions from the downscaling of the pseudo-reality are compared with corresponding results derived by the interpolation.
The downscaled pseudo-reality can also be compared with the results downscaled for the actual station data. This is done in ‘DSensemble’ when the argument ‘non.stationarity.check’ is
set to ‘TRUE’.

### iid.test
The iid.test included in the esd package is used to test whether a variable is independent and identically distributed (iid) mainly for daily station records (Benestad , 2003, 2004). It has to be
noted that this test is sensitive to missing data (NA) and can produce an under-count. A non i.i.d. behaviour appears when the ’forward’ (resp. ’backward’) analysis indicates higher (resp.
lower) number of record-events than the confidence interval. 

## Diagnose
The ‘diagnose’ method provides a range of different approaches for the evaluation of ’esd’ results depending on the type of objects. For instance, it returns diagnostics of common EOFs
or cross-validation of DS results. The diagnostic also includes combined fields, MVR, and CCA results by applying appropriate internal tests to check for consistency. For instance, when applied
to a combined ‘eof’ object, ‘diagnose’ investigates the difference in the mean between the PCs of the calibration data and the PCs of the GCMs over a common period in addition to the
ratio of standard deviations and lag-one auto-correlation. A bias correction method can then be applied if the difference is significant as described in (Imbert and Benestad , 2005).
Figure 26: The results of ‘diagnose’ applied to four different data objects: (a) common field, (b) to
common EOFs, (c) downscaled results for an ensemble, and (d) a station.

### Station
The ‘diagnose’ for station object offers another way of plotting the information contents and is designed for a group of stations to indicate the availability of valid data (Example 6.1,
Figure 26d).
### Combined fields
The application of ‘diagnose’ to a combined pair of field objects will produce two comparable curves that show their spatial mean values (Example 6.2, Figure 26a). This way, the spatial
statistics of the two fields can be compared.
### Common EOFs
For common EOFs, ‘diagnose’ returns a set of different diagnostics, such as mean differences, the standard deviation ratio, and different auto-correlation for the principal components over
an overlapping interval (6.3).

Example 6.1.
```R
# Retrieve NACD data sets
nacd <- station(src=’nacd’)
# Diagnose for data availability
diagnose(nacd)
Example 6.2.
# Sample 2m air temperature ERAINT field from EOF’s data
t2m <- t2m.ERAINT(lon=c(-30,40),lat=c(50,70))
# Retrieve GCM air temperature from file
# Sample 2m air temperature NorESM field from EOF’s data
gcm <- t2m.NorESM.M()
# Combine annual fields
X <- combine(annual(t2m),annual(gcm))
# Diagnose the result
diagnose(X)
Example 6.3.
# Compute EOF, X is taken from previous example
eof <- EOF(X)
# Diagnose eof
a <- diagnose(eof)
# Plot the results
plot(a)
```
## Downscaled results
The quality of the downscaled results can be investigated and checked with the function `diagnose.ds` (Example 6.4). In addition cross-validation and residuals can be used to evaluate
the downscaling procedure.

### Cross-validation
A default option is to perform a five-fold cross-validation test (’crossval’)), but it can also be set to perform the CORDEX-ESD experiment 1 or VALUE experiment 1 set-up. The results of the cross-validation consist of a zoo-object that is stored in the attribute ‘attr(*,’evaluation’).
### Residuals
The residuals of the downscaling should ideally be a white noise that has no auto-correlation, no trend, and is normally distributed. These conditions can be tested using auto-correlation
functions and time series and qqnorm plots.

### Downscaled ensemble results
For ensembles ‘diagnose’ shows a comparison between DS results derived from GCMs and the past variations. It compares the number of observed values outside the predicted 90%
confidence interval based on fitting a normal distribution to the ensemble member values for each year, in addition to carrying out an evaluation of trends based on the rank the observed
trend has compared to corresponding trends simulated by the ensemble for the past. The diagnostics produces an info-graphic in the form of a target (Example 6.5, Fig 26d).
The count of the cases where observations fall outside the 90% confidence interval is expected to follow a binomial distribution if the model results reproduce the correct statistics (H0). Hence,
p = 0.1 for falling outside, and the test is to see if the actual count of cases outside the confidence interval is consistent with the H0 binomial distribution (Benestad and Mezghani , 2015).

Example 6.4. # Load temperature data for Oslo
```R
data(Oslo)
# Get the 2m air temperature from NCEP reanalysis
t2m <- t2m.NCEP(lon=c(-5,20),lat=c(50,65))
# Compute the EOF on annual temperature values
eof <- EOF(annual(t2m))
# Downscale the annual values
z <- DS(annual(Oslo),eof)
# Do the diagnostic
diagnose(z,plot=TRUE)
Example 6.5. # Download the file from the figshare link as
download.file(url="http://files.figshare.com/2081407/dse.oslo.rda",destfile="dseOslo.rda")
# Load the file into R session
load(’dseOslo.rda’)
# Do the diagnostic
diagnose(dse.oslo)
```
