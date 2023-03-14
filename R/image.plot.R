#' Plot image
#'
#' @param breaks A numeric vector of breakpoints for the colours
#' @param add if TRUE add to current plot
#' @param nlevel the number of breaks (used only if breaks are not specified),
#' @param horizontal if TRUE legend is horizontal
#' @param legend.shrink shrinkage factor for legend
#' @param legend.mar margins of legend (see \code{link{par}})
#' @param legend.lab legend label
#' @param legend.line margin line of the legend
#' @param graphics.reset a boolean
#' @param bigplot A vector of the form ‘c(x1, x2, y1, y2)’ giving the
#'        coordinates of the plot region as fractions of the current
#'        figure region.
#' @param smallplot Same as bigplot but for a second plot
#' @param legend.only if TRUE show only legend
#' @param col Specification of the plotting color (a single color or a vector)
#' @param pal color palette, used only if col is NULL (see \code{\link{colscal}})
#' @param lab.breaks labels of breakpoints for the colours (argument \code{"breaks"})
#' @param axis.args list containing arguments for axis
#' @param legend.args list containing arguments for legend
#' @param midpoint If TRUE color scale is formed for midpoints by averaging 4 corners.
#' @param border the color to draw the border. Use NA to omit the border. 
#' @param lwd width of line
#' @param rev if TRUE reverse palette
#' @param verbose if TRUE print progress
#' @param \dots additional arguments
#'
#' @importFrom graphics box
#'
#' @export image.plot
image.plot <- function (..., breaks=NULL, add = FALSE, nlevel = 64, horizontal = FALSE, 
                        legend.shrink = 0.9, legend.width = 1.2, 
                        legend.mar = ifelse(horizontal, 3.1, 5.1), 
                        legend.lab = NULL, legend.line = 2, graphics.reset = FALSE, 
                        bigplot = NULL, smallplot = NULL, legend.only = FALSE, col = NULL, pal="heat", 
                        lab.breaks = NULL, axis.args = NULL, legend.args = NULL, 
                        midpoint = FALSE, border = NA, lwd = 1, rev=FALSE,verbose=FALSE) {
  
  if(verbose) print("image.plot")
  par0 <- par(no.readonly = TRUE) # save default, for resetting...
  args <- names(list(...))
  if (verbose) print(args)
  if(!is.null(breaks)) nlevel <- length(breaks)-1
  if(is.null(col)) col <- colscal(n=nlevel, pal=pal)
  if (rev) col <- rev(col)
  if (length(breaks) != (length(col)+1)) {
    print(breaks)
    print(col)
    browser()
  }
  if (verbose) print(c(nlevel,pal))
  info <- imageplot.info(verbose=verbose,...)
  if (verbose) print(info)
  if (add) {
    big.plot <- par0$plt
  }
  if (legend.only) {
    graphics.reset <- TRUE
  }
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  if (verbose) print('set up plot')
  temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
                          legend.width = legend.width, legend.mar = legend.mar, 
                          horizontal = horizontal, bigplot = bigplot, smallplot = smallplot, verbose=verbose)
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot
  
  if (!legend.only) {
    if (!add) {
      par(plt = bigplot)
    }
    if (!info$poly.grid) {
      image(..., add = add, col = col)
    } else {
      poly.image(..., add = add, col = col, midpoint = midpoint, 
                 border = border, lwd.poly = lwd, verbose=verbose)
    }
    #big.par <- par(no.readonly = TRUE)  # This has been commented out...
  }
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(par0)
    stop("plot region too small to add legend\n")
  }
  
  ix <- 1
  minz <- info$zlim[1]
  maxz <- info$zlim[2]
  ## KMP 2015-09-23: for unevenly spaced breaks
  if(is.null(breaks)) {
    if (verbose) print('unevenly spaced breaks')
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
  } else {
    if (verbose) print('evenly spaced breaks')
    z <- unique(c(minz,breaks,maxz))
    binwidth <- diff(z)
    midpoints <- z[1:(length(z)-1)]+binwidth/2
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
  }
  par(new = FALSE, pty = "m", plt = smallplot, err = -1)
  if (!is.null(breaks) & !is.null(lab.breaks)) {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(2, 0.5, 0), las = ifelse(horizontal, 0, 2), 
                        at = breaks, labels = lab.breaks), axis.args)
  } else {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(2, 0.5, 0), las = ifelse(horizontal, 0, 2)), 
                   axis.args)
  }
  if (verbose) str(axis.args)
  if (!horizontal) {
    if (verbose) print("not horizontal")
    if (is.null(breaks)) {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
    } else {
      if (verbose) print("horizontal")
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col, breaks = breaks)
    }
  } else {
    if (is.null(breaks)) {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
    } else {
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
  if (verbose) {print('image.plot - par0()'); print(par0$mfcol); print(par0$mfrow); print(par0$mfg)}
  par(pty=par0$pty, plt=par0$plt, err=par0$err) # reset to default
}


# internal function - no need to export
imageplot.info <- function (verbose=FALSE,...) {
  if(verbose) print("imageplot.info")
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

# internal function - no need to export
poly.image <- function (x, y, z, col = colscal(n=64,pal="heat"), breaks, transparent.color = "white", 
                        midpoint = FALSE, zlim = range(z, na.rm = TRUE), xlim = range(x), 
                        ylim = range(y), add = FALSE, border = NA, lwd.poly = 1, verbose=FALSE, ...) {
  if(verbose) print("poly.image")
  if(!requireNamespace("fields",quietly=TRUE)) {
    stop("Package \"fields\" needed to regrid image. Please install it.")
  } else {
    Dx <- dim(x)
    Dy <- dim(y)
    if (any((Dx - Dy) != 0)) {
      stop(" x and y matrices should have same dimensions")
    }
    Dz <- dim(z)
    if (all((Dx - Dz) == 0) & !midpoint) {
      x <- fields::poly.image.regrid(x)
      y <- fields::poly.image.regrid(y)
    }
    if (missing(breaks)) {
      breaks <- NA
    }
    zcol <- fields::drape.color(z, col = col, midpoint = midpoint, zlim = zlim, 
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
}

# internal function - no need to export
imageplot.setup <- function (x, add = FALSE, legend.shrink = 0.9, legend.width = 1, 
                             horizontal = FALSE, legend.mar = NULL, bigplot = NULL, smallplot = NULL, 
                             verbose=FALSE, ...) {
  if(verbose) print("imageplot.setup")
  par0 <- par() 
  if (is.null(smallplot)) {
    stick <- TRUE
  } else {
    stick <- FALSE
  }
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  char.size <- ifelse(horizontal, par0$cin[2]/par0$din[2], 
                      par0$cin[1]/par0$din[1])
  offset <- char.size * ifelse(horizontal, par0$mar[1], par0$mar[4])
  legend.width <- char.size * legend.width
  legend.mar <- legend.mar * char.size
  if (verbose) print(paste('legend.mar=',legend.mar,'char.size=',char.size,'offset=',offset,'legend.width=',legend.width))
  if (is.null(smallplot)) {
    smallplot <- par0$plt
    if (horizontal) {
      smallplot[3] <- legend.mar
      smallplot[4] <- legend.width + smallplot[3]
      pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink)/2)
      smallplot[1] <- smallplot[1] + pr
      smallplot[2] <- smallplot[2] - pr
    } else {
      smallplot[2] <- 1 - legend.mar
      smallplot[1] <- smallplot[2] - legend.width
      pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
      smallplot[4] <- smallplot[4] - pr
      smallplot[3] <- smallplot[3] + pr
    }
  }
  if (is.null(bigplot)) {
    bigplot <- par0$plt
    if (!horizontal) {
      bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
    } else {
      bottom.space <- par0$mar[1] * char.size
      bigplot[3] <- smallplot[4] + offset
    }
  }
  if (stick & (!horizontal)) {
    dp <- smallplot[2] - smallplot[1]
    smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
    smallplot[2] <- smallplot[1] + dp
  }
  if (verbose) {print(list(smallplot = smallplot, bigplot = bigplot)); print('exit imageplot.setup')}
  return(list(smallplot = smallplot, bigplot = bigplot))
}
