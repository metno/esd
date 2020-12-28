---
layout: home
title : Downscaling
nav_order: 6
has_children: false
toc: true
permalink: /downscale/
---

# Downscaling
{: .no_toc }
Downscaling tries to identify (statistical) existing relationships between large-scale spatial climate patterns (called predictors) and local climate variations (predictands) in climate variables, parameters, and/or indicators (indices).

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

## Defining predictands
In downscaling based on ‘esd’, the predictand can be a station data (‘station’ class), a groups of stations (‘pca’), gridded data (‘eof’), or trajectory statistics. The strategy differs from that of ‘clim.pact’ by downscaling seasonal or annual mean statistics, rather than monthly and daily. The reason for this is the idea that there is a stronger connection between local and large-scale climate than between local and large-scale weather. Here ‘climate’ is defined as weather statistics and involves parameters derived from a sample of some variable, e.g. the mean or frequency over a season or a year. It follows that statistics such as the mean is estimated more accurately with larger sample sizes according to the central limit theorem, reducing the effect of both sampling fluctuations and random errors. Furthermore, according to the central limit theorem, the mean converges to being normally distributed with increasing sample size.
The seasonal and annual scales also imply more data points, especially for precipitation - it doesn’t rain every day, and a monthly estimate may typically involve about 10 rain events - or 10 data points. One motivation for using PCA to represent the predictand is that it brings forth the coherent variations and the signal in the data that is expected to have the strongest connection with the large scales (Benestad et al., 2015). Furthermore, the use of PCA takes care of spatial cross correlation among the different stations.
The downscaling model only considers variability - the mean level is prescribed based on the observations and a given reference period. Hence, it is not affected by biases in the mean (a form for systematic errors). The climatology of the end-result is set by the climatology of the observations (base period) and the downscaled results describe the deviation from this climatology.
The downscaling uses common EOFs (Barnett, 1999) as a basis to ensure that the covariates used in the regression represent the exact same spatial structures (modes) when calibrating the regression model against reanalysis and when using the model to make predictions based on GCM results. The use of common EOFs is described in Benestad (2001). 

## Defining predictors
The large-scale predictors are represented in terms of EOFs and common EOFs in ‘esd’, where a regression analysis is used to estimate best-fit weights (regression coefficients βi) whose products with the PCs (ˆy(t) = Pn β i=1 iVi(t)) give the closest representation of the original time series y.
n
X
ˆ
y(t) = β0 +
βiVi(t)
(3)
i=1
The term ‘covariate’ will henceforth be used to refer to the time-dependent variables on the right hand side of the multiple regression equation (Vi(t)). 
The default is that the regression is based on ordinary linear model, but there is an option in the argument ‘method’ to use generalised linear models or any other method that produces results with similar structure. The downscaling will also involve a forward-backward stepping, removing the PCs that do not increase the skill.
The downscaling can combine several predictors of different types by using a list object as the predictor argument, containing the EOFs describing the different sets of predictors. In Example 5.1, a combination of the T(2m) and SLP are used as predictors for downscaling.
The EOFs contained in the list are combined to a single ‘eof’ class, and then subject to the default downscaling as if the predictor represented one single variable. The combination of the predictors involves a weighting of each principal component according to its eigenvalue, and then a singular value decomposition (Strang, 1988; Press et al., 1989b) is used to distill the most pronounced and common variability embedded in the different data sets, in a similar way ‘mixed’ EOFs were used in Benestad et al. (2002) (or ‘mixed’ EOFs in ‘clim.pact’), described
as ‘CPCA’ in Bretherton et al. (1992). Sets of several different types of predictors are combined in terms of their principal components (PC) according to the following strategy: the PCs for each set are weighted by a set of normalised weights to ensure that each type of predictor carry similar weight in addition to giving more weight to the leading EOFs that explain more of the variance (For each respective EOF, P w i i = 1). The weighted PCs are then combined and a singular value decomposition (SVD) is used to extract the common variability within the chosen
predictors, and then the 20 leading modes are used to synthesise a data object that resembles the EOF object. This synthesised matrix is then used as a predictor that describes the most pronounced common/combined variability of the chosen predictor types. 

### Tuning
The downscaling can be tuned after the set of predictors is chosen. The predictor domain is by default set to a fixed spatial window (width to be specified by the user), but different choices for the domain will affect the results (Benestad , 2001). Also the number of EOFs included will affect the outcome, and the default is set to include the 7 leading modes. The number of EOFs can be determined manually as the first components with the maximum of variability (This can be achieved by looking at ‘plot(eof)’).

## Options for downscaling
There are several downscaling approaches depending on the temporal resolution and the final objective of the impact study i.e. one can perform a DS for each year, season or month then combine the 12 months. For instance, annual downscaling is based on a new strategy which means extract annual statistics from daily or monthly data and do the downscaling.

### Downscaling for single station and a single model
The basic downscaling is done for a single predictand (station) and a single predictor (EOF), and there are various ‘wrappers’ that use this elementary part as pieces of a more complicated process. The predictor can be a common EOF, constructed using ‘combine’ on two fields and then an estimation of the EOFs. The downscaling can also be applied to a group of stations or gridded fields.

### Downscaling a group of stations - PCA
One benefit of downscaling ’PCA’ objects (using ‘DS.pca’ built-in function) rather than downscaling each single station separately is that the PCA-based analysis places more emphasis on the common signal found in adjacent stations and takes care about any strong spatial pattern that could be found in nearby stations. In other words, it may enhance the signal-to-noise ratio by separating out the signal (leading PCAs) from the noise (high-level PCAs). Another aspect is that the downscaling takes into account the spatial correlation structure, embedded in the
PCA where eachmode is orthogonal (Benestad et al., 2015).

Example 5.2 and Figure 23 give a demonstration and a test of the downscaling of PCA.

### Downscaling an ensemble of climate models
The method ‘DSensemble’ is designed to deal with downscaling an ensemble of global or regional climate models such as CMIP3 and CMIP5 runs. The function ‘DSensemble’ reads all data specified in a common folder (specified in the input arguments) to do the downscaling, and each result is saved separately for each climate model. Each file in the climate models’ folder must be unique and must correspond to one climate model output stored in a NetCDF format.
The GCM or RCM files should not be stored in several files spanning, for instance, different time slices. In this case, the files have to be concatenated before applying ‘DSensemble’. 
‘DSensemble’ has also been extended to deal with a group of stations’ object (‘pca’), with a speed-up and benefits associated with PCA.
The design of the ‘DSensemble’ methods tries to meet some of the criticism that Estrada et al. (2013) presented against ESD that are summarised in the following issues.

* a. The underlying probabilistic model assumed by an ESD method does not reflect the distribution of the underlying data.
* b. The method requires the variables to be stationary.
* c. ESD focuses on the short-term variability (Huth’s dilemma).
* d. Estimated relations of trend and/or auto-correlated series are at risk of being spurious.
* e. Statistical downscaling is an empirical method and the statistical adequacy should be evaluated using appropriate tools.
* f. Downscaling tool boxes seldom give any information regarding the significance of the estimated coefficients.

The downscaling is applied to samples of data aiming at predicting parameters describing the sample distribution. For temperature, the sample size is N ≈90 data points (seasonal aggregates) whereas for precipitation, the downscaling is applied to annual aggregates because it doesn’t rain every day (typical size N ≈100). The idea is to downscale climate (defined as weather statistics) rather than weather, hence calibrating the models on aggregated parameters rather than day-to-day situations (issue ’c’). Furthermore, the primary parameters are the seasonal
means for temperature and annual wet-day means and frequency for precipitation. According to the central limit theorem, the distribution of mean estimates is expected to converge to being normal with sample size (Wheelan, 2013). Hence, the default method (‘lm’) is appropriate for such parameters (issue ’a’), although other methods such as GLM may be used for other types of distributions.
The usual requirement for ESD is that the transfer coefficients describing the link between the large and small scales are stationary, and not the variables themselves. Climate change by definition implies non-stationary variables, such a trend in temperature. In both esd and clim.pact the default is to de-trend the data before model calibration, in order to reduce the risk of spurious results (issues ‘b’ & ‘d’). The auto-correlation between successive winters (or any other season) or years is low and forecasts for local temperature and precipitation aggregated for a year ahead is notoriously difficult as a result, which implies that auto-correlation has little effect on the outcome of the downscaled results in the default setting (issue ‘d’). When it comes to model evaluation tools (issue ‘e’), the method ’DS’, on which ‘DSensemble’ is based, uses cross-validation by default for assessing its skill. The R-environment also offers a set of tools for testing whether the data used for model calibration are normally distributed: e.g. ’qqnorm’. A Shapiro-Wilk test of normality is also included in ’diagnose.ds’ . Furthermore, ’DS’ invokes a step-wise screening that removes the predictors which do not make a statistical significant contribution to the regression (issue ‘f’), and esd includes a function that carries out an iid-test which can be applied to de-trended series: ‘iid.test’ (Benestad , 2008).

## Downscaling of trajectory objects
The function ‘DS’ can be applied to trajectory objects, but it is not the paths themselves that are downscaled but rather some statistical measure of each trajectory (Example 5.4, Figure 25). By default, the downscaling is applied to the monthly trajectory count. Trajectory objects may also be downscaled with regards to other aspects by adding a parameter name (‘param’) and a function to apply to the parameter for each path (‘FUN’). For example, we can calculate and downscale the genesis latitude (param=’lat’, FUN=’first’) or the mean longitude (param=’lon’,FUN=’mean’).

Example 5.1.
```R
# Load eofs of NCEP 2m air temperature
data(eof.t2m.NCEP)
# Load eofs of NCEP sea level pressure
data(eof.slp.NCEP)
# Load temperature data for Oslo
data(Oslo)
# Do the downscaling
z <- DS(Oslo,list(t2m=eof.t2m.NCEP,slp=eof.slp.NCEP),mon=1)
# Plot the results
plot(z)
```
Figure 22: An example of downscaled annual mean temperature in Oslo, based on the ERAINT reanalysis and using PCA for downscaling a group of stations simultaneously.

Example 5.2.
```R
# Load ERAINT 2m air temperature from file
T2m <- t2m.ERAINT(lon=c(-10,30),lat=c(50,70))
# Compute the EOFs on annual values
eof <- EOF(annual(T2m))
# Retrieve Nordklim data sets
y <- station(src=’nordklim’)
# Get rid of stations with little valid data
Y <- allgood(annual(y))
# Do the PCA
pca <- PCA(Y)
# DO the downscaling
z <- DS(pca,eof)
# Plot the result
plot(z)
# Get the observations corresponding to the first station in the list
obs <- subset(Y,is=1)
# Get the predictions
pre <- subset(as.station(z),is=1)
# Match dates between obs and pre
obs <- matchdate(obs,pre)
# Update the attributes
attr(obs,’location’) <- ’Observed’
attr(pre,’location’) <- ’Predicted’
# Combine the two
obspre <- combine.stations(obs,pre)
# Plot the result in a single plot
plot(obspre,plot.type="single")
```
Figure 23: An example of downscaled annual mean temperature in Oslo, based on the ERAINT reanalysis and using PCA for downscaling a group of stations simultaneously.

Example 5.3.
```R
# Get the predictand: Maximum temperature for some Norwegian sttions
download.file(url="http://files.figshare.com/2073466/Norway.Tx.rda",destfile="NorwayTx.rda")
load("NorwayTx.rda")
# Process the predictand (annual mean)
ma <- annual(Tx)
pca <- PCA(ma)
# Get the ERAINT T(2m):
t2m <- t2m.ERAINT(lon=c(-30,30),lat=c(50,75))
# Process the predictor: Compute EOFs of annual means
eof <- EOF(annual(X))
# Test: CCA:
cca <- CCA(pca,eof)
# Do the downscaling
z <- DS(pca,eof,rmtrend=FALSE)
# Plot the results for one station:
plot(subset(ma,is=1))
lines(subset(as.station(pca),is=is),lty=2,col="darkred")
lines(subset(as.station(z),is=is),lwd=2)
```
Figure 24: An example of downscaled CMIP5 RCP8.5 DJF seasonal mean temperature eastern Norway. The ERA40 reanalysis is used here for calibration. The inner plot in the right bottom of the figure shows a comparison between the estimated linear trend from both simulations (in terms of distribution) and observations (black point). While, the inner plot on the top right of the figure shows the distribution of the number of observations lying outside the 90% of the confidence interval based on simulations.

Example 5.4.
```R
# Load storm trajectories
data(imilast.M03)
# Load NCEP sea level pressure data for the northern hemisphere
slp <- slp.NCEP(lat=c(0,90))
# Compute eof
eof.slp <- EOF(slp)
# Do downscaling for the storm count
(default trajectory downscaling)
ds <- DS(imilast.M03,eof.slp)
# Do downscaling from slp
ds.slp <- DS(imilast.M03,eof.slp,param=’slp’,FUN=’min’,unit=’hPa’)
# Plot downscaled results
plot(ds)
plot(ds.slp)
```
Figure 25: Two examples of downscaling the storm tracks of the North Atlantic region with regards to the storm count (a) and minimum slp (b).

## Downscaling probabilities and number of events
The original time series can be summarised as the probability distribution of a climate variable or parameter, the latter can then be downscaled instead of the original time series. For instance, the 24-hour wet-day precipitation amount is expected to follows an exponential or gamma distribution, whereas the temperature is more likely to follow a normal distribution, and the wind-speed a Weibull distribution. The count of events is expected to follow a binomial or Poisson distribution and spell-length statistics tends to have a geometric distribution. The predictands here are the parameters of the probability distribution function rather than the original climate variable. These are mainly the first order and second order moments such as the mean (for exponential distribution), the mean and the standard deviation for the normal distribution, and the mean and the coefficient of variation for the gamma distribution. 
Downscaling of µ, fw, and ncwd (wet spell length) require different choices in the type of regression model because they involve different types of data (Estrada et al., 2013). The set-up may involve continuous values, fractions (f ∈ [0, 1]), and discrete counts. sometimes it may be difficult to decide which category to use. For example, the annual wet-day frequency is a mean estimated over a set of zeros and ones and exhibits a character closer to normal than logistic distribution, but while each day would be subject to a logistic distribution, its annual mean is subject to the central limit theorem.

## Weather generators - synthesising daily time series
Weather generators are usually used to simulate realistic sequences of daily weather and climate variables such as precipitation, maximum and minimum temperature for present (Semenov and Barrow , 1997; Wilks and Wilby, 1999) and future climate conditions (Wilks, 1992). Weather generators can be seen as complementary tools and are used in conjunction with downscaling techniques. They use (downscaled) parameters of the several distribution functions of climate variables. The generation process differs between climate parameters, variables and statistics. For precipitation, most weather generators will generate first dry and wet sequences, and then for the wet days, the precipitation amount is generated following different probability and using different random numbers form different distributions (Mezghani and Hingray, 2009). The simulation process takes into account the various specific aspects of climate variability and change (Soltani and Hoogenboom, 2003). The simulated sequences must share the main statistical properties of the observed or original time series used for calibration (Semenov and Brooks, 1999).
The general form of a WG can be written as follows:
$$ Y = f (X, ǫ)$$, (4)
where $ǫ$ is an ensemble of random variables which contains the remaining information that has not been taken into account by X and the random processes that is involved.

There are various ways a weather generator can be designed, and ’esd’ may include more types in the future. Presently, there are two different types, one for precipitation and one for temperature. These are designed so that they are conditioned upon downscaled results for seasonal of annual scales, using vital parameters to disaggregate the results to daily values. The results will span a similar time interval as the provided observations x but may be shifted in time if a set of projections is provided (taking an interval that ends on the last projected year).

### Stochastic precipitation
The stochastic or generation process for precipitation uses information about wet-day mean, frequency, and consecutive wet days to generate stochastic precipitation series. The call is ‘WG.fw.day.precip’ and its features are (roughly):
- Estimate wet-day mean and frequency as well as the wet- and dry-spell statistics for a station.
- Uses phase-scrambled Fourier transform of annual observed/downscaled wet-day mean and frequency.
- Estimate the number of wet days per year.
- Wet days: amounts prescribed according to the exponential distribution and µ.
- Try to distribute wet days according to both the number of wet days per year and with a wet-spell statistics according to annual mean number of consecutive wet days assuming a geometric distribution.

### Stochastic temperature
The weather generator for daily temperature is ‘WG.FT.day.t2m’, currently designed for single stations. The key features of this weather generator are (roughly):
- Uses phase-scrambled Fourier transform of observed or downscaled annual mean temperature anomaly to produce a random order of annual data.
- Uses phase-scrambled Fourier transform of observed or downscaled annual standard deviation of daily temperature anomalies to produce a random order of annual data.
- Uses quantile-quantile mapping to generate normally distributed stochastic variables conditioned on prescribed mean and standard deviation.
