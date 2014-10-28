image.plot <- function (..., add = FALSE, nlevel = 64, horizontal = FALSE, 
    legend.shrink = 0.9, legend.width = 1.2, legend.mar = ifelse(horizontal, 
        3.1, 5.1), legend.lab = NULL, legend.line = 2, graphics.reset = FALSE, 
    bigplot = NULL, smallplot = NULL, legend.only = FALSE, col = tim.colors(nlevel), 
    lab.breaks = NULL, axis.args = NULL, legend.args = NULL, 
    midpoint = FALSE, border = NA, lwd = 1) 
{
    
    old <- par()
    ## print("old") ; print(old$fig)
    old.par <- par(no.readonly = TRUE)
    info <- imageplot.info(...)
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (!info$poly.grid) {
            image(..., add = add, col = col)
        }
        else {
            poly.image(..., add = add, col = col, midpoint = midpoint, 
                border = border, lwd.poly = lwd)
        }
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!is.null(breaks) & !is.null(lab.breaks)) {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
            at = breaks, labels = lab.breaks), axis.args)
    }
    else {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
            axis.args)
    }
    if (!horizontal) {
        if (is.null(breaks)) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col)
        }
        else {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks = breaks)
        }
    }
    else {
        if (is.null(breaks)) {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col)
        }
        else {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks = breaks)
        }
    }
    do.call("axis", axis.args)
    box()
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
            1, 4), line = legend.line)
    }
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    ## print(old$fig)
    par(fig = old$fig)
    par(mfcol=old$mfcol)
}



imageplot.info <- function (...) 
{
    temp <- list(...)
    xlim <- NA
    ylim <- NA
    zlim <- NA
    poly.grid <- FALSE
    if (is.list(temp[[1]])) {
        xlim <- range(temp[[1]]$x, na.rm = TRUE)
        ylim <- range(temp[[1]]$y, na.rm = TRUE)
        zlim <- range(temp[[1]]$z, na.rm = TRUE)
        if (is.matrix(temp[[1]]$x) & is.matrix(temp[[1]]$y) & 
            is.matrix(temp[[1]]$z)) {
            poly.grid <- TRUE
        }
    }
    if (length(temp) >= 3) {
        if (is.matrix(temp[[1]]) & is.matrix(temp[[2]]) & is.matrix(temp[[3]])) {
            poly.grid <- TRUE
        }
    }
    if (is.matrix(temp[[1]]) & !poly.grid) {
        xlim <- c(0, 1)
        ylim <- c(0, 1)
        zlim <- range(temp[[1]], na.rm = TRUE)
    }
    if (length(temp) >= 3) {
        if (is.matrix(temp[[3]])) {
            xlim <- range(temp[[1]], na.rm = TRUE)
            ylim <- range(temp[[2]], na.rm = TRUE)
            zlim <- range(temp[[3]], na.rm = TRUE)
        }
    }
    if (is.matrix(temp$x) & is.matrix(temp$y) & is.matrix(temp$z)) {
        poly.grid <- TRUE
    }
    xthere <- match("x", names(temp))
    ythere <- match("y", names(temp))
    zthere <- match("z", names(temp))
    if (!is.na(zthere)) 
        zlim <- range(temp$z, na.rm = TRUE)
    if (!is.na(xthere)) 
        xlim <- range(temp$x, na.rm = TRUE)
    if (!is.na(ythere)) 
        ylim <- range(temp$y, na.rm = TRUE)
    if (!is.null(temp$zlim)) 
        zlim <- temp$zlim
    if (!is.null(temp$xlim)) 
        xlim <- temp$xlim
    if (!is.null(temp$ylim)) 
        ylim <- temp$ylim
    list(xlim = xlim, ylim = ylim, zlim = zlim, poly.grid = poly.grid)
}

poly.image <- function (x, y, z, col = tim.colors(64), breaks, transparent.color = "white", 
    midpoint = FALSE, zlim = range(z, na.rm = TRUE), xlim = range(x), 
    ylim = range(y), add = FALSE, border = NA, lwd.poly = 1, 
    ...) 
{
    Dx <- dim(x)
    Dy <- dim(y)
    if (any((Dx - Dy) != 0)) {
        stop(" x and y matrices should have same dimensions")
    }
    Dz <- dim(z)
    if (all((Dx - Dz) == 0) & !midpoint) {
        x <- poly.image.regrid(x)
        y <- poly.image.regrid(y)
    }
    if (missing(breaks)) {
        breaks <- NA
    }
    zcol <- drape.color(z, col = col, midpoint = midpoint, zlim = zlim, 
        transparent.color = transparent.color, breaks = breaks)$color.index
    if (!add) {
        plot(xlim, ylim, type = "n", ...)
    }
    N <- ncol(x)
    Nm1 <- N - 1
    M <- nrow(x)
    Mm1 <- M - 1
    for (i in (1:Mm1)) {
        xp <- cbind(x[i, 1:Nm1], x[i + 1, 1:Nm1], x[i + 1, 2:N], 
            x[i, 2:N], rep(NA, Nm1))
        yp <- cbind(y[i, 1:Nm1], y[i + 1, 1:Nm1], y[i + 1, 2:N], 
            y[i, 2:N], rep(NA, Nm1))
        xp <- c(t(xp))
        yp <- c(t(yp))
        pcol <- c(zcol[i, 1:Nm1])
        polygon(xp, yp, border = pcol, col = pcol, lwd = lwd.poly)
        if (!is.na(border)) {
            polygon(xp, yp, border = border, col = NA, lwd = lwd.poly)
        }
    }
}
imageplot.setup <- function (x, add = FALSE, legend.shrink = 0.9, legend.width = 1, 
    horizontal = FALSE, legend.mar = NULL, bigplot = NULL, smallplot = NULL, 
    ...) 
{
    old.par <- par(no.readonly = TRUE)
    if (is.null(smallplot)) 
        stick <- TRUE
    else stick <- FALSE
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    char.size <- ifelse(horizontal, par()$cin[2]/par()$din[2], 
        par()$cin[1]/par()$din[1])
    offset <- char.size * ifelse(horizontal, par()$mar[1], par()$mar[4])
    legend.width <- char.size * legend.width
    legend.mar <- legend.mar * char.size
    if (is.null(smallplot)) {
        smallplot <- old.par$plt
        if (horizontal) {
            smallplot[3] <- legend.mar
            smallplot[4] <- legend.width + smallplot[3]
            pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink)/2)
            smallplot[1] <- smallplot[1] + pr
            smallplot[2] <- smallplot[2] - pr
        }
        else {
            smallplot[2] <- 1 - legend.mar
            smallplot[1] <- smallplot[2] - legend.width
            pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
            smallplot[4] <- smallplot[4] - pr
            smallplot[3] <- smallplot[3] + pr
        }
    }
    if (is.null(bigplot)) {
        bigplot <- old.par$plt
        if (!horizontal) {
            bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
        }
        else {
            bottom.space <- old.par$mar[1] * char.size
            bigplot[3] <- smallplot[4] + offset
        }
    }
    if (stick & (!horizontal)) {
        dp <- smallplot[2] - smallplot[1]
        smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
        smallplot[2] <- smallplot[1] + dp
    }
    return(list(smallplot = smallplot, bigplot = bigplot))
}
