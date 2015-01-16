clean.station <- function(x,miss=.1,verbose=TRUE) {

    stopifnot(inherits(x,"station"))
    loc <- loc(x)
    if (verbose) par(mfrow=c(2,1))
    if (verbose) image(coredata(annual(x)))
    if (verbose) title('Original data')
    # compute number of finite values
    n <- apply(coredata(x),2,FUN=nv)

    if (verbose) summary(n)
    keep <- n > (dim(x)[1]*miss)
    x <- subset(x,is=keep)
    if (verbose) print("Removed stations ...")
    if (verbose) print(loc[!keep])
    if (verbose) image(coredata(annual(x)))
    if (verbose) title('Filtered data')
    return(x)

}
