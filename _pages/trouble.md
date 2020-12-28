---
layout: home
title: Trouble
toc: true
permalink: /trouble/
---

# Trouble shooting
The ‘esd’ comes as it is and there is no guarantee that it’s free from bugs or incorrect coding. The user therefore must make sure that the functions work as expected through thorough testing.
Any tool like this should never be used blindly. Neither is there any guarantee that the code will work flawlessly. The open-source nature, however, means that the problems can be diagnosed
or fixed by anyone with insight and programming skills. Here we try to make potential trouble shooting easier.

## General
In case of error messages, it may be wise to check all the attributes of the data object. The way to proceed is to activate the ‘verbose’ argument if it is included in the function options,
and repeat the analysis with the sample data provided by ‘esd’. The objects can be assigned new attributes through the ‘attr’ command. Alternatively, the source code of the function can
be modified locally, by typing the name of the function without ‘()’ and copied into a local file (Example 7.1). The code in the local file can then be modified and when it is called override
the original ‘esd’ version.

Example 7.1.
```R
# get the source code in ’myeof.R’ R script file
dput(EOF.default,file="myeof.R")
# Modify this code, but keep the same function name by adding ’EOF.default’ at the top of
# the file.
# Source the new code as
source("myeof.R",local=environment())
# now the call "EOF()" will use the modified version.
```
Debugging can be done by inserting ‘browser()’ inside the function code.

## Functions ‘annual’ and ‘aggregate’
Missing data (‘NA’) can be a problem in the analysis and the downscaling of local temperature
and precipitation. One example of an error is:
```R
y <- annual(x)
Error in aggregate.data.frame(as.data.frame(x), ...) :
no rows to aggregate
```

However, if the number of missing data is small compared to the sample size, they may not have a large effect on the calibration of downscaling models or the aggregated statistics. For the
temperature, a fix can be to use interpolation to fill in the missing values as in Example 7.2. For daily precipitation, a better option can be to replace missing data with zero.

Example 7.2.
```R
# Retrieve 2m air temperature at Ferder Fyr weather station
data(ferder)
# Extract the coredata
z <- coredata(ferder)
# identify the missing days in the record
md <- is.na(ferder)
# Compute the climatology
clim <- as.climatology(ferder)
# Replace the missing values with the climatology as approximations
ferder[md] <- coredata(clim)[as.POSIXlt(index(ferder)[md])$yday+1]
# Note that as.POSIXlt()$yday function returns the day in the year
```

## DSensemble
The method ‘DSensemble’ uses ‘try’ to avoid that the process stops due to errors when looping through the ensemble of GCM results. However, it cannot recover from errors in For-
tran or C-based external modules, such as those based on external libraries. For instance, the following errors have been encountered:
Error in La.svd(x, nu, nv) : error code 1 from Lapack routine ’dgesdd’ and/or Error in x[rep(1:NROW(x), length.out = length(index)), , drop = FALSE] : subscript out of bounds
The problem is that the data matrix to which an SVD is applied is close to singular. One work-around fix is to use a different domain selection (e.g. set ‘lon’ and ‘lat’ arguments to
new values). Another way to get around this problem is to use the argument ‘select’ to skip a particular model with the problem.
### Poor fit
Sometimes a sub-selection of the predictand can improve the calibration of the downscaling model, especially if part of the data (older) has a lower quality than more recent measurements.
One way to examine the data for quality can be to plot the time series of the annual wet-day mean and frequency for precipitation.
```R
plot(annual(y,FUN=’wetfreq’))
plot(annual(y,FUN=’wetmean’))
```
### Other error messages
There may occasionally be problems associated with different types of time stamps, e.g. when combining daily, monthly, seasonal, or annual data. The time index for annual data is a numeric
- the year, and some of the methods handling the date may fail for annual data. E.g: Error in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L,: invalid ’trim’ argument

A couple of solutions are: 
```R
index(y) <- as.Date(paste(year(y),’01-01’,sep=’-’))
## or
index(y) <- year(y)
```
## Validate
The ‘esd’ package also has a method ‘validate’ which provides a simple test statistics describing how well the results match the expectations. The null hypothesis depends on the
type of object.
