## Author A. Mezghani
## Description Contains some rtools ...
## Created 14.11.2014

as.decimal <- function(x=NULL) {
    ## converts from degree min sec format to degrees ...
    ##x is in the form "49°17´38´´"
    if (!is.null(x)) {
        deg <-as.numeric(substr(x,1,2)) 
        min <- as.numeric(substr(x,4,5))
        sec <- as.numeric(substr(x,7,8))     
        x <- deg + min/60 + sec/3600
    }
    return(x)
}


## compute the percentage of missing data in x
missval <- function(x) sum(is.na(coredata(x)))/length(coredata(x))

## compute the quantile 95% of x
q95 <- function(x) quantile(x,probs=.95,na.rm=TRUE)

## compute the quantile 5% of x
q5 <- function(x) quantile(x,probs=.05,na.rm=TRUE)

## compute the quantile 5% of x
q995 <- function(x) quantile(x,probs=.995,na.rm=TRUE)

## compute the quantile 5% of x
q975 <- function(x) quantile(x,probs=.975,na.rm=TRUE)

## count the number of valid data points
nv <- function(x) sum(is.finite(x))

## Compute the coefficient of variation of x
cv <- function(x,na.rm=TRUE) {sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm)}
